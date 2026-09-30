"""Enterprise in-silico toxicology screening engine.

Research-use computational screening only. The engine deliberately separates
structural alerts, transformation hypotheses, exposure calculations,
regulatory-rule matches, read-across similarity, provenance, and uncertainty.
"""
from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from io import BytesIO
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import requests

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import HRFlowable, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle, Image as PlatypusImage
from reportlab.pdfgen import canvas

ENGINE_VERSION = "2.0.0"
REGULATORY_RULESET_VERSION = "research-ruleset-2026.09"
DEFAULT_BODY_WEIGHT_KG = 60.0
MOS_REFERENCE_THRESHOLD = 100.0


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe_float(value: Any, default: float = np.nan) -> float:
    try:
        x = float(value)
        return x if np.isfinite(x) else default
    except (TypeError, ValueError):
        return default


def _html(value: Any) -> str:
    text = str(value if value is not None else "")
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def _validate_molecule(smiles: str) -> Tuple[Optional[Chem.Mol], Optional[str]]:
    if not isinstance(smiles, str) or not smiles.strip():
        return None, "SMILES input is empty."
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
    except Exception as exc:
        return None, f"RDKit could not parse the SMILES: {exc}"
    if mol is None:
        return None, "Invalid SMILES string passed to the engine."
    try:
        Chem.SanitizeMol(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception as exc:
        return None, f"Molecular sanitization failed: {exc}"
    return mol, None


# ---------------------------------------------------------------------------
# 1. CUTANEOUS TRANSFORMATION HYPOTHESES
# ---------------------------------------------------------------------------
CUTANEOUS_REACTION_SMARTS = {
    "Cutaneous Carboxylesterase (Hydrolysis)": "[CX3:1](=[OX1:2])[OX2][C:4]>>[CX3:1](=[OX1:2])[OH].[C:4][OH]",
    "Cutaneous Phenol/Catechol Oxidation Hypothesis": "[c:1]1[c:2]([OH:7])[c:3]([OH:8])[c:4][c:5][c:6]1>>[C:1]1=[C:6][C:5]=[C:4][C:3](=[O:8])[C:2]1=[O:7]",
    "Cutaneous Alcohol Dehydrogenase (Allylic Alcohol Oxidation Hypothesis)": "[C:1]=[C:2]-[CH2:3][OH]>>[C:1]=[C:2]-[CH:3]=O",
    "Cutaneous Epoxidation (Alkene Activation Hypothesis)": "[C:1]=[C:2]>>[C:1]1[O][C:2]1",
    "Cutaneous Aromatic Amine N-Hydroxylation Hypothesis": "[c:1][NH2:2]>>[c:1][NH:2][OH]",
}

CUTANEOUS_REACTION_RULES: Dict[str, Any] = {}
for _name, _smarts in CUTANEOUS_REACTION_SMARTS.items():
    try:
        _rxn = AllChem.ReactionFromSmarts(_smarts)
        if _rxn: CUTANEOUS_REACTION_RULES[_name] = _rxn
    except Exception:
        continue


def simulate_skin_metabolism(parent_mol: Optional[Chem.Mol]) -> List[Dict[str, Any]]:
    if parent_mol is None: return []
    try: parent_smiles = Chem.MolToSmiles(parent_mol)
    except Exception: return []

    records: List[Dict[str, Any]] = []
    seen_smiles = {parent_smiles}
    
    for pathway_name, rxn in CUTANEOUS_REACTION_RULES.items():
        try:
            product_sets = rxn.RunReactants((parent_mol,))
        except Exception as exc:
            continue

        pathway_products = 0
        for product_tuple in product_sets:
            for metabolite in product_tuple:
                try:
                    Chem.SanitizeMol(metabolite)
                    met_smiles = Chem.MolToSmiles(metabolite)
                    if met_smiles in seen_smiles or len(met_smiles) <= 2: continue
                    seen_smiles.add(met_smiles)
                    pathway_products += 1
                    alerts = run_structural_filter(metabolite)
                    has_alert = bool(alerts) and not all("Passed" in a for a in alerts)
                    records.append({
                        "pathway": pathway_name,
                        "metabolite_smiles": met_smiles,
                        "is_prohapten": False,
                        "prohapten_status": "Not established from structure alone",
                        "secondary_alerts": alerts,
                        "structural_alert_detected": has_alert,
                        "status": "Transformation hypothesis generated",
                        "evidence_level": "Computational hypothesis",
                    })
                except Exception:
                    continue
                    
    return records


# ---------------------------------------------------------------------------
# 2. STRUCTURAL AND COMMERCIAL RULES
# ---------------------------------------------------------------------------
RETAILER_RESTRICTIONS = {
    "Paraben Ester": {"smarts": "O=C(Oc1ccc(O)cc1)[#6]", "sephora": "Banned", "credo": "Banned", "eu_status": "Restricted"},
    "Ortho-Phthalate Ester": {"smarts": "O=C(c1ccccc1C(=O)O[#6])O[#6]", "sephora": "Banned", "credo": "Banned", "eu_status": "Restricted"},
    "Formaldehyde Releaser (Bronopol Core)": {"smarts": "OCC(Br)([$([N+](=O)[O-]),$(N(=O)=O)])CO", "sephora": "Banned", "credo": "Banned", "eu_status": "Restricted"},
    "Butylated Hydroxytoluene (BHT) Substructure": {"smarts": "Oc1c(C(C)(C)C)cc(C(C)(C)C)cc1", "sephora": "Restricted (<0.1%)", "credo": "Banned", "eu_status": "Restricted"},
    "Cyclic Siloxane (D4/D5 core)": {"smarts": "[Si]1O[Si]O[Si]O[Si]O1", "sephora": "Banned", "credo": "Banned", "eu_status": "Restricted"},
    "Triclosan Phenoxy Core": {"smarts": "c1cc(Cl)c(Oc2ccc(Cl)cc2Cl)cc1O", "sephora": "Banned", "credo": "Banned", "eu_status": "Restricted"},
    "Alkyl Sulfate Surfactant (SLS/SLES)": {"smarts": "[#6]OS(=O)(=O)[O-,OH]", "sephora": "Banned", "credo": "Restricted", "eu_status": "Monitored"},
}

EU_EXPANDED_ALLERGENS = {
    "Limonene/Terpene Hydrocarbon": "CC1=CCC(CC1)C(=C)C",
    "Linalool Core": "CC(=CCCC(C)(C=C)O)C",
    "Cinnamaldehyde Derivative": "c1ccccc1C=CC=O",
    "Cinnamyl Alcohol Derivative": "c1ccccc1C=CCO",
    "Eugenol/Isoeugenol Core": "c1cc(c(cc1CC=C)OC)O",
    "Coumarin Core": "O=C1OC2=CC=CC=C2C=C1",
    "Citral (Geranial/Neral)": "CC(=CCCC(=CC=O)C)C",
    "Benzyl Salicylate": "c1ccccc1C(=O)OCc2ccccc2O",
    "Hydroxycitronellal": "CC(CCCC(C)(C)O)CC=O",
    "Farnesol Core": "CC(=CCCC(=CCCC(=CCO)C)C)C",
}


def run_structural_filter(mol: Optional[Chem.Mol]) -> List[str]:
    if mol is None: return ["Invalid molecular entity"]
    try:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
        catalog = FilterCatalog(params)
        alerts = [match.GetDescription() for match in catalog.GetMatches(mol)]
        return alerts or ["Passed: no BRENK/PAINS/NIH catalog alert detected"]
    except Exception as exc:
        return [f"Structural catalog evaluation error: {exc}"]


def _compile_smarts(smarts: str) -> Optional[Chem.Mol]:
    try: return Chem.MolFromSmarts(smarts)
    except Exception: return None


def audit_retailer_standards(mol: Optional[Chem.Mol]) -> List[Dict[str, Any]]:
    findings: List[Dict[str, Any]] = []
    if mol is None: return findings
    for rule_name, data in RETAILER_RESTRICTIONS.items():
        pattern = _compile_smarts(data["smarts"])
        if pattern is not None and mol.HasSubstructMatch(pattern):
            findings.append({"liability": rule_name, "sephora": data["sephora"], "credo": data["credo"], "eu_status": data["eu_status"]})
    return findings


def audit_eu_allergens(mol: Optional[Chem.Mol]) -> List[str]:
    detected: List[str] = []
    if mol is None: return detected
    for name, smarts in EU_EXPANDED_ALLERGENS.items():
        pattern = _compile_smarts(smarts)
        if pattern is not None and mol.HasSubstructMatch(pattern):
            detected.append(name)
    return detected


# ---------------------------------------------------------------------------
# 3. EXPOSURE / MoS
# ---------------------------------------------------------------------------
SCCS_EXPOSURE_DEFAULTS = {
    "Face Cream (Leave-on)": {"daily_amount_g": 1.54, "retention_factor": 1.0},
    "Body Lotion (Leave-on)": {"daily_amount_g": 7.82, "retention_factor": 1.0},
    "Facial Cleanser (Rinse-off)": {"daily_amount_g": 20.0, "retention_factor": 0.01},
    "Shampoo (Rinse-off)": {"daily_amount_g": 10.46, "retention_factor": 0.01},
    "Eye Cream (Leave-on)": {"daily_amount_g": 0.50, "retention_factor": 1.0},
    "Lip Balm (Leave-on)": {"daily_amount_g": 0.06, "retention_factor": 1.0},
}

def calculate_margin_of_safety(product_type: str, concentration_pct: float, noael_mg_kg_day: float, dermal_absorption_pct: float = 50.0, body_weight_kg: float = DEFAULT_BODY_WEIGHT_KG) -> Dict[str, Any]:
    defaults = SCCS_EXPOSURE_DEFAULTS.get(product_type, {"daily_amount_g": 1.54, "retention_factor": 1.0})
    concentration = _safe_float(concentration_pct)
    noael = _safe_float(noael_mg_kg_day)
    absorption = _safe_float(dermal_absorption_pct)
    bw = _safe_float(body_weight_kg)
    
    if not all(np.isfinite(x) for x in (concentration, noael, absorption, bw)) or bw <= 0 or concentration < 0 or noael <= 0 or not 0 <= absorption <= 100:
        return {"sed_mg_kg_day": np.nan, "mos": np.nan, "verdict": "Indeterminate", "is_safe": False, "valid": False}

    applied_mg_day = defaults["daily_amount_g"] * 1000.0
    amount_in_product_mg = applied_mg_day * (concentration / 100.0) * defaults["retention_factor"]
    systemic_absorbed_mg = amount_in_product_mg * (absorption / 100.0)
    sed = systemic_absorbed_mg / bw
    
    if sed <= 0: return {"sed_mg_kg_day": 0.0, "mos": float("inf"), "verdict": "Indeterminate", "is_safe": False, "valid": True}
    
    mos = noael / sed
    return {
        "sed_mg_kg_day": round(sed, 6),
        "mos": round(mos, 1),
        "verdict": "Threshold met" if mos >= MOS_REFERENCE_THRESHOLD else "Threshold failed",
        "is_safe": bool(mos >= MOS_REFERENCE_THRESHOLD),
        "valid": True,
        "threshold": MOS_REFERENCE_THRESHOLD,
    }


# ---------------------------------------------------------------------------
# 4. READ-ACROSS
# ---------------------------------------------------------------------------
REFERENCE_APPROVED_ACTIVES = [
    {"name": "Bakuchiol", "inci": "Bakuchiol", "category": "Antioxidant / Retinoid Alternative", "smiles": "CC(=CCC/C(=C/CC1=CC=C(C=C1)O)/C)C"},
    {"name": "Niacinamide", "inci": "Niacinamide", "category": "Skin Barrier / Restorative", "smiles": "NC(=O)c1cccnc1"},
    {"name": "Squalane", "inci": "Squalane", "category": "Emollient", "smiles": "CCCCCCCCCCCC(C)CCCC(C)CCCC(C)CCCC(C)C"},
    {"name": "Allantoin", "inci": "Allantoin", "category": "Soothing / Keratolytic", "smiles": "O=C1NC(=O)NC1NC(=O)N"},
    {"name": "Panthenol", "inci": "Panthenol", "category": "Humectant", "smiles": "OCC(C)(C)C(O)C(=O)NCCCCO"},
]

def find_read_across_analogs(target_mol: Optional[Chem.Mol], top_k: int = 3) -> List[Dict[str, Any]]:
    if target_mol is None: return []
    try:
        morgan_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        target_fp = morgan_generator.GetFingerprint(target_mol)
    except Exception: return []
    
    candidates = []
    for item in REFERENCE_APPROVED_ACTIVES:
        ref_mol = Chem.MolFromSmiles(item["smiles"])
        if ref_mol is None: continue
        ref_fp = morgan_generator.GetFingerprint(ref_mol)
        similarity = float(DataStructs.TanimotoSimilarity(target_fp, ref_fp))
        candidates.append({
            "name": item["name"], "inci": item["inci"], "category": item["category"],
            "similarity": round(similarity * 100, 1), "interpretation": "Structural similarity only."
        })
        
    candidates.sort(key=lambda x: x["similarity"], reverse=True)
    return candidates[:max(1, int(top_k))]


# ---------------------------------------------------------------------------
# 5. PUBCHEM & ORCHESTRATION
# ---------------------------------------------------------------------------
def get_pubchem_regulatory_summary(smiles: str) -> Dict[str, Any]:
    result = {"status": "Unavailable", "cid": None}
    try:
        url_cid = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/JSON"
        response = requests.post(url_cid, data={"smiles": smiles}, timeout=6.0)
        if response.status_code == 200:
            cid = response.json()["IdentifierList"]["CID"][0]
            result["cid"] = cid
            safety_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON?heading=Safety%20and%20Hazards"
            safety_response = requests.get(safety_url, timeout=6.0)
            result["status"] = "Dossier located" if safety_response.ok else "Compound located"
    except Exception:
        pass
    return result

def _source_hash(payload: Any) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode("utf-8")).hexdigest()[:16]

def execute_full_compound_audit(smiles: str, product_type: str, concentration: float, noael: float, dermal_absorption: float, body_weight_kg: float = DEFAULT_BODY_WEIGHT_KG) -> Dict[str, Any]:
    mol, error = _validate_molecule(smiles)
    if mol is None: return {"valid": False, "error": error}

    clean_smiles = Chem.MolToSmiles(mol)
    alerts = run_structural_filter(mol)
    retailer = audit_retailer_standards(mol)
    allergens = audit_eu_allergens(mol)
    metabolic = simulate_skin_metabolism(mol)
    mos = calculate_margin_of_safety(product_type, concentration, noael, dermal_absorption, body_weight_kg)
    analogs = find_read_across_analogs(mol)
    pubchem = get_pubchem_regulatory_summary(clean_smiles)

    structural_flag = any("Passed" not in a for a in alerts)
    metabolic_flag = any(r.get("structural_alert_detected") for r in metabolic)
    mos_failed = bool(mos.get("valid") and not mos.get("is_safe"))

    if retailer or mos_failed or (structural_flag and metabolic_flag): verdict = "CRITICAL RISK: REFORMULATION ADVISED"
    elif structural_flag or metabolic_flag or allergens: verdict = "MODERATE LIABILITY: MONITORING REQUIRED"
    else: verdict = "LOW RISK: ACCEPTABLE FOR PRE-CLINICAL BATCHING"

    return {
        "valid": True, "mol": mol, "clean_smiles": clean_smiles, "alerts": alerts,
        "retailer_violations": retailer, "eu_allergens": allergens, "metabolic_records": metabolic,
        "mos_metrics": mos, "analogs": analogs, "pubchem_status": pubchem, "composite_verdict": verdict,
        "assessment_id": f"TX-{_source_hash({'smiles': clean_smiles, 'time': _utc_now()})}",
        "exposure_params": {"product_type": product_type, "concentration": concentration, "noael": noael, "dermal_absorption": dermal_absorption}
    }


# ---------------------------------------------------------------------------
# 6. ENTERPRISE PDF RENDERER (MULTI-PAGE)
# ---------------------------------------------------------------------------
class NumberedCanvas(canvas.Canvas):
    def __init__(self, *args, **kwargs):
        canvas.Canvas.__init__(self, *args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.draw_page_number(num_pages)
            canvas.Canvas.showPage(self)
        canvas.Canvas.save(self)

    def draw_page_number(self, page_count):
        self.saveState()
        self.setFont("Helvetica", 8)
        self.setFillColor(colors.HexColor("#718096"))
        self.drawString(36, 20, "Confidential Pre-Clinical Dossier | Lead Computational Biologist: Debanjan Gangopadhyay")
        self.drawRightString(576, 20, f"Page {self._pageNumber} of {page_count}")
        self.restoreState()

def generate_enterprise_pdf(audit_result: Dict[str, Any]) -> BytesIO:
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter, rightMargin=36, leftMargin=36, topMargin=36, bottomMargin=36)
    styles = getSampleStyleSheet()
    title = ParagraphStyle("DocTitle", parent=styles["Normal"], fontName="Helvetica-Bold", fontSize=18, leading=22, textColor=colors.HexColor("#1B365D"))
    h1 = ParagraphStyle("H1", parent=styles["Normal"], fontName="Helvetica-Bold", fontSize=12, leading=15, textColor=colors.HexColor("#1B365D"), spaceBefore=10, spaceAfter=5)
    body = ParagraphStyle("Body", parent=styles["Normal"], fontSize=8.5, leading=11, textColor=colors.HexColor("#2D3748"))
    cell = ParagraphStyle("Cell", parent=body, fontSize=8)
    bold_cell = ParagraphStyle("BoldCell", parent=cell, fontName="Helvetica-Bold")

    story = [Paragraph("IN-SILICO TOXICOLOGY & BIOACTIVATION DOSSIER", title), Spacer(1, 8), HRFlowable(width="100%", thickness=1, color=colors.HexColor("#1B365D")), Spacer(1, 10)]
    
    img = Draw.MolToImage(audit_result["mol"], size=(220, 220))
    img_buf = BytesIO(); img.save(img_buf, format="PNG"); img_buf.seek(0)
    structure = PlatypusImage(img_buf, width=120, height=120)
    
    summary = Paragraph(f"<b>Assessment ID:</b> {_html(audit_result['assessment_id'])}<br/><b>Target SMILES:</b> {_html(audit_result['clean_smiles'])}<br/><br/><b>Composite Verdict:</b> {_html(audit_result['composite_verdict'])}", body)
    
    top = Table([[structure, summary]], colWidths=[135, 405])
    top.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'MIDDLE'),('BACKGROUND',(0,0),(-1,-1),colors.HexColor('#F7FAFC')),('BOX',(0,0),(-1,-1),0.8,colors.HexColor('#CBD5E0')),('PADDING',(0,0),(-1,-1),8)]))
    story += [top, Spacer(1, 10)]

    story.append(Paragraph("1. Cutaneous Bioactivation & Pro-Hapten Screening", h1))
    rows = [[Paragraph("Enzymatic Pathway", bold_cell), Paragraph("Metabolite Generated", bold_cell), Paragraph("Secondary Alert Flag", bold_cell)]]
    for rec in audit_result["metabolic_records"][:10]:
        rows.append([Paragraph(_html(rec.get("pathway")), cell), Paragraph(_html(rec.get("metabolite_smiles") or "-"), cell), Paragraph(_html(", ".join(rec.get("secondary_alerts", [])) or "None"), cell)])
    table = Table(rows, colWidths=[180, 180, 180])
    table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#EDF2F7')),('GRID',(0,0),(-1,-1),0.4,colors.HexColor('#CBD5E0')),('PADDING',(0,0),(-1,-1),4)]))
    story += [table, Spacer(1, 10)]

    story.append(Paragraph("2. Safety & Exposure", h1))
    exp = audit_result["exposure_params"]
    mos = audit_result["mos_metrics"]
    rows = [[Paragraph("Concentration", cell), Paragraph(f"{_html(exp['concentration'])} %", cell)], [Paragraph("Margin of Safety (MoS)", cell), Paragraph(str(mos.get('mos')), cell)]]
    table = Table(rows, colWidths=[270, 270])
    table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#EDF2F7')),('GRID',(0,0),(-1,-1),0.4,colors.HexColor('#CBD5E0')),('PADDING',(0,0),(-1,-1),4)]))
    story += [table, Spacer(1, 10)]

    story.append(Paragraph("3. Commercial & Regulatory Screening", h1))
    if audit_result["retailer_violations"]:
        rows = [[Paragraph("Liability", bold_cell), Paragraph("Sephora", bold_cell), Paragraph("Credo", bold_cell), Paragraph("EU Status", bold_cell)]]
        for v in audit_result["retailer_violations"]:
            rows.append([Paragraph(_html(v["liability"]), cell), Paragraph(_html(v["sephora"]), cell), Paragraph(_html(v["credo"]), cell), Paragraph(_html(v["eu_status"]), cell)])
        table = Table(rows, colWidths=[180, 120, 120, 120])
        table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#EDF2F7')),('GRID',(0,0),(-1,-1),0.4,colors.HexColor('#CBD5E0')),('PADDING',(0,0),(-1,-1),4)]))
        story.append(table)
    else:
        story.append(Paragraph("Passed: No local retailer structural-rule match detected.", body))
    
    story.append(Spacer(1, 5))
    story.append(Paragraph("EU Allergen Matches: " + _html(", ".join(audit_result["eu_allergens"]) or "None detected"), body))

    doc.build(story, canvasmaker=NumberedCanvas)
    buffer.seek(0)
    return buffer
