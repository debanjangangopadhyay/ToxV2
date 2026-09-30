import json
import os
import hashlib
from datetime import datetime, timezone
from io import BytesIO
import numpy as np
from typing import Dict, List, Any, Optional

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw, Descriptors, rdMolDescriptors, rdFingerprintGenerator
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import HRFlowable, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle, Image as PlatypusImage
from reportlab.pdfgen import canvas

# ============================================================================
# DEMPSTER-SHAFER INTEGRATION ENGINE
# ============================================================================
class EvidenceMass:
    def __init__(self, safe: float, toxic: float, uncertain: float):
        total = safe + toxic + uncertain
        self.safe = safe / total
        self.toxic = toxic / total
        self.uncertain = uncertain / total

def combine_evidence(m1: EvidenceMass, m2: EvidenceMass) -> EvidenceMass:
    safe_int = (m1.safe * m2.safe) + (m1.safe * m2.uncertain) + (m1.uncertain * m2.safe)
    toxic_int = (m1.toxic * m2.toxic) + (m1.toxic * m2.uncertain) + (m1.uncertain * m2.toxic)
    unc_int = m1.uncertain * m2.uncertain
    
    conflict_k = (m1.safe * m2.toxic) + (m1.toxic * m2.safe)
    if conflict_k >= 0.999: return EvidenceMass(0.0, 0.0, 1.0)
    
    norm = 1.0 - conflict_k
    return EvidenceMass(safe_int / norm, toxic_int / norm, unc_int / norm)

# ============================================================================
# CORE ANALYTICAL MODULES
# ============================================================================
def execute_provenance_audit(mol: Chem.Mol, schema_path: str = "regulatory_rules.json") -> Dict[str, Any]:
    res = {"valid": False, "violations": [], "warnings": []}
    if not os.path.exists(schema_path): return res
    try:
        with open(schema_path, "r") as f: rules = json.load(f).get("regulatory_rules", [])
        for rule in rules:
            pat = Chem.MolFromSmarts(rule["smarts"])
            if pat and mol.HasSubstructMatch(pat):
                flag = {k: rule[k] for k in ["rule_id", "target_class", "status", "legal_instrument"]}
                if "Prohibited" in rule["status"]: res["violations"].append(flag)
                else: res["warnings"].append(flag)
        res["valid"] = True
    except Exception: pass
    return res

def run_structural_filter(mol: Chem.Mol) -> List[str]:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    return [match.GetDescription() for match in catalog.GetMatches(mol)]

def simulate_skin_metabolism(parent_mol: Chem.Mol) -> List[Dict[str, Any]]:
    rxns = {
        "Ester Hydrolysis": "[CX3:1](=[OX1:2])[OX2][C:4]>>[CX3:1](=[OX1:2])[OH].[C:4][OH]",
        "Catechol Oxidation": "[c:1]1[c:2]([OH:7])[c:3]([OH:8])[c:4][c:5][c:6]1>>[C:1]1=[C:6][C:5]=[C:4][C:3](=[O:8])[C:2]1=[O:7]",
        "Allylic Oxidation": "[C:1]=[C:2]-[CH2:3][OH]>>[C:1]=[C:2]-[CH:3]=O"
    }
    records = []
    for name, smarts in rxns.items():
        try:
            products = AllChem.ReactionFromSmarts(smarts).RunReactants((parent_mol,))
            for p_tuple in products:
                for met in p_tuple:
                    Chem.SanitizeMol(met)
                    alerts = run_structural_filter(met)
                    records.append({
                        "pathway": name,
                        "smiles": Chem.MolToSmiles(met),
                        "alerts": alerts,
                        "is_prohapten": len(alerts) > 0
                    })
        except Exception: continue
    return records

def calculate_probabilistic_mos(pod: float, conc_pct: float, abs_pct: float, n: int = 10000) -> Dict[str, Any]:
    np.random.seed(42)
    af_target = 10.0 * 10.0 * 3.0  # Base inter/intra/duration AF
    
    bw = np.random.normal(60.0, 10.2, n).clip(40.0, 120.0)
    da = np.random.normal(abs_pct, abs_pct * 0.3, n).clip(0.1, 100.0) / 100.0
    amt = np.random.normal(1.54, 0.23, n).clip(0.5, 3.0) * 1000.0
    pod_dist = np.random.normal(pod, pod * 0.1, n).clip(pod * 0.5, None)
    
    sed = (amt * (conc_pct / 100.0) * da) / bw
    mos = pod_dist / sed
    fail_prob = float(np.sum(mos < af_target) / n)
    
    return {
        "median_mos": round(float(np.median(mos)), 1),
        "ci_05_mos": round(float(np.percentile(mos, 5)), 1),
        "failure_probability": fail_prob,
        "target_af": af_target
    }

def execute_multidimensional_read_across(target_mol: Chem.Mol) -> List[Dict[str, Any]]:
    refs = [
        {"name": "Niacinamide", "smiles": "NC(=O)c1cccnc1"},
        {"name": "Squalane", "smiles": "CCCCCCCCCCCC(C)CCCC(C)CCCC(C)CCCC(C)C"},
        {"name": "Allantoin", "smiles": "O=C1NC(=O)NC1NC(=O)N"}
    ]
    try:
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        t_fp = gen.GetFingerprint(target_mol)
        t_logp, t_mw = Descriptors.MolLogP(target_mol), Descriptors.MolWt(target_mol)
        
        cands = []
        for r in refs:
            r_mol = Chem.MolFromSmiles(r["smiles"])
            if not r_mol: continue
            dist_top = 1.0 - DataStructs.TanimotoSimilarity(t_fp, gen.GetFingerprint(r_mol))
            dist_logp = min(abs(t_logp - Descriptors.MolLogP(r_mol)) / 5.0, 1.0)
            dist_mw = min(abs(t_mw - Descriptors.MolWt(r_mol)) / 300.0, 1.0)
            
            score = max((1.0 - ((dist_top * 0.6) + (dist_logp * 0.3) + (dist_mw * 0.1))) * 100.0, 0.0)
            cands.append({"name": r["name"], "homology_score": round(score, 1), "delta_logp": round(abs(t_logp - Descriptors.MolLogP(r_mol)), 2)})
        return sorted(cands, key=lambda x: x["homology_score"], reverse=True)
    except Exception: return []

# ============================================================================
# MASTER ORCHESTRATION
# ============================================================================
def execute_full_compound_audit(smiles: str, pod: float, conc: float, da: float) -> Dict[str, Any]:
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return {"valid": False, "error": "Invalid SMILES"}
    
    reg = execute_provenance_audit(mol)
    alerts = run_structural_filter(mol)
    met = simulate_skin_metabolism(mol)
    mos = calculate_probabilistic_mos(pod, conc, da)
    analogs = execute_multidimensional_read_across(mol)
    
    # Dempster-Shafer Fusion
    m_mos = EvidenceMass((1.0 - mos["failure_probability"]) * 0.9, mos["failure_probability"], 0.1 - (mos["failure_probability"] * 0.1))
    m_reg = EvidenceMass(0.0, 0.98, 0.02) if reg["violations"] else EvidenceMass(0.0, 0.60, 0.40) if alerts else EvidenceMass(0.40, 0.0, 0.60)
    top_score = analogs[0]["homology_score"] / 100.0 if analogs else 0.0
    m_ana = EvidenceMass(top_score * 0.8, 0.0, 1.0 - (top_score * 0.8))
    
    final_mass = combine_evidence(combine_evidence(m_mos, m_reg), m_ana)
    
    return {
        "valid": True, "mol": mol, "smiles": smiles,
        "dst_metrics": {
            "belief_safe": round(final_mass.safe * 100, 1),
            "plausibility_safe": round((final_mass.safe + final_mass.uncertain) * 100, 1),
            "epistemic_uncertainty": round(final_mass.uncertain * 100, 1)
        },
        "reg": reg, "alerts": alerts, "met": met, "mos": mos, "analogs": analogs,
        "assessment_id": f"TX-{hashlib.sha256(smiles.encode()).hexdigest()[:8].upper()}"
    }

# ============================================================================
# PAGINATED PDF GENERATION
# ============================================================================
class NumberedCanvas(canvas.Canvas):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._saved_page_states = []
    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()
    def save(self):
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.saveState()
            self.setFont("Helvetica", 8)
            self.setFillColor(colors.HexColor("#718096"))
            self.drawString(36, 20, "Confidential Pre-Clinical Dossier | Lead Computational Biologist: Debanjan Gangopadhyay")
            self.drawRightString(576, 20, f"Page {self._pageNumber} of {len(self._saved_page_states)}")
            self.restoreState()
            super().showPage()
        super().save()

def generate_enterprise_pdf(data: Dict[str, Any]) -> BytesIO:
    buf = BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=letter, rightMargin=36, leftMargin=36, topMargin=36, bottomMargin=36)
    sty = getSampleStyleSheet()
    h1 = ParagraphStyle("H1", parent=sty["Normal"], fontName="Helvetica-Bold", fontSize=12, textColor=colors.HexColor("#1B365D"), spaceBefore=10)
    c_b = ParagraphStyle("CB", parent=sty["Normal"], fontSize=8, fontName="Helvetica-Bold")
    c_n = ParagraphStyle("CN", parent=sty["Normal"], fontSize=8)
    
    img_buf = BytesIO(); Draw.MolToImage(data["mol"], size=(180, 180)).save(img_buf, format="PNG"); img_buf.seek(0)
    
    story = [
        Paragraph("IN-SILICO TOXICOLOGY & BIOACTIVATION DOSSIER", ParagraphStyle("T", fontName="Helvetica-Bold", fontSize=18)),
        HRFlowable(width="100%", thickness=1, color=colors.HexColor("#1B365D")), Spacer(1, 10),
        Table([[PlatypusImage(img_buf, width=120, height=120), Paragraph(f"<b>ID:</b> {data['assessment_id']}<br/><b>Belief (Safe):</b> {data['dst_metrics']['belief_safe']}%<br/><b>Uncertainty:</b> {data['dst_metrics']['epistemic_uncertainty']}%", sty["Normal"])]], colWidths=[130, 400], style=[('BACKGROUND',(0,0),(-1,-1),colors.HexColor('#F7FAFC')),('BOX',(0,0),(-1,-1),0.5,colors.HexColor('#CBD5E0'))]),
        Spacer(1, 10), Paragraph("1. Monte Carlo Margin of Safety (10,000 Iterations)", h1),
        Table([[Paragraph("Median MoS", c_b), Paragraph("5th Percentile MoS", c_b), Paragraph("Failure Probability", c_b)], [Paragraph(str(data["mos"]["median_mos"]), c_n), Paragraph(str(data["mos"]["ci_05_mos"]), c_n), Paragraph(f"{data['mos']['failure_probability']*100}%", c_n)]], colWidths=[170, 170, 170], style=[('GRID',(0,0),(-1,-1),0.5,colors.HexColor('#CBD5E0'))])
    ]
    doc.build(story, canvasmaker=NumberedCanvas)
    buf.seek(0)
    return buf
    
