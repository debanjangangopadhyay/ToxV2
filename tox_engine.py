"""
Enterprise Probabilistic In-Silico Toxicology & IATA Engine (Engine 4 - Refactored)
Features scenario-aware Monte Carlo Margin of Safety simulations, Shafer-discounted
Dempster-Shafer evidence fusion, structural similarity, and paginated PDF generation.
"""

import os
import json
import html
import math
import hashlib
from datetime import datetime, timezone
from io import BytesIO
from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw, Descriptors, rdFingerprintGenerator
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import HRFlowable, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle, Image as PlatypusImage
from reportlab.pdfgen import canvas


class ValidationError(ValueError):
    """Raised when numerical inputs or physical boundaries are violated."""
    pass


# ============================================================================
# 1. INPUT CONTRACT & SCENARIO EXPOSURE DEFAULTS
# ============================================================================

SCCS_PRODUCT_EXPOSURE = {
    "Face Cream (Leave-on)": {"daily_amount_g": 1.54, "sigma_g": 0.23, "retention_factor": 1.0},
    "Body Lotion (Leave-on)": {"daily_amount_g": 7.82, "sigma_g": 1.10, "retention_factor": 1.0},
    "Facial Cleanser (Rinse-off)": {"daily_amount_g": 20.0, "sigma_g": 2.50, "retention_factor": 0.01},
    "Shampoo (Rinse-off)": {"daily_amount_g": 10.46, "sigma_g": 1.80, "retention_factor": 0.01},
    "Eye Cream (Leave-on)": {"daily_amount_g": 0.50, "sigma_g": 0.08, "retention_factor": 1.0},
    "Lip Balm (Leave-on)": {"daily_amount_g": 0.06, "sigma_g": 0.01, "retention_factor": 1.0}
}


@dataclass(frozen=True)
class AnalysisInputContract:
    smiles: str
    product_type: str
    concentration_pct: float
    pod_noael_mg_kg_day: float
    dermal_absorption_pct: float = 50.0
    body_weight_kg: float = 60.0
    mc_samples: int = 10000
    random_seed: int = 42

    def __post_init__(self):
        if not isinstance(self.smiles, str) or not self.smiles.strip():
            raise ValidationError("Input SMILES must be a non-empty string.")

        if self.product_type not in SCCS_PRODUCT_EXPOSURE:
            raise ValidationError(
                f"Invalid product_type '{self.product_type}'. "
                f"Allowed scenarios: {list(SCCS_PRODUCT_EXPOSURE.keys())}"
            )

        floats = {
            "concentration_pct": self.concentration_pct,
            "pod_noael_mg_kg_day": self.pod_noael_mg_kg_day,
            "dermal_absorption_pct": self.dermal_absorption_pct,
            "body_weight_kg": self.body_weight_kg
        }
        for name, val in floats.items():
            if not isinstance(val, (int, float)):
                raise ValidationError(f"Parameter '{name}' must be numeric. Got: {type(val)}")
            if math.isnan(val) or math.isinf(val):
                raise ValidationError(f"Analytical Integrity Failure: '{name}' contains NaN or Inf.")

        if not (0.0 < self.concentration_pct <= 100.0):
            raise ValidationError(f"concentration_pct must be in (0, 100]. Got: {self.concentration_pct}")
        if not (0.0 < self.dermal_absorption_pct <= 100.0):
            raise ValidationError(f"dermal_absorption_pct must be in (0, 100]. Got: {self.dermal_absorption_pct}")
        if self.pod_noael_mg_kg_day <= 0.0:
            raise ValidationError(f"pod_noael_mg_kg_day must be > 0. Got: {self.pod_noael_mg_kg_day}")
        if self.body_weight_kg <= 0.0:
            raise ValidationError(f"body_weight_kg must be > 0. Got: {self.body_weight_kg}")
        if not (1000 <= self.mc_samples <= 1000000):
            raise ValidationError(f"mc_samples must be between 1,000 and 1,000,000. Got: {self.mc_samples}")


# ============================================================================
# 2. EVIDENCE FUSION ENGINE (SHAFER-DISCOUNTED DEMPSTER-SHAFER)
# ============================================================================

@dataclass
class EvidenceMass:
    safe: float
    toxic: float
    uncertain: float

    def __post_init__(self):
        total = self.safe + self.toxic + self.uncertain
        if not math.isclose(total, 1.0, abs_tol=1e-5):
            self.safe /= total
            self.toxic /= total
            self.uncertain /= total

    def apply_discounting(self, alpha: float) -> 'EvidenceMass':
        """Applies Shafer Discounting Operator: alpha in [0, 1] represents evidence reliability/independence."""
        alpha = max(0.0, min(1.0, alpha))
        return EvidenceMass(
            safe=alpha * self.safe,
            toxic=alpha * self.toxic,
            uncertain=1.0 - alpha * (self.safe + self.toxic)
        )


def combine_evidence(m1: EvidenceMass, m2: EvidenceMass) -> EvidenceMass:
    """Fuses two independent evidence masses via Dempster's Rule of Combination."""
    k_conflict = (m1.safe * m2.toxic) + (m1.toxic * m2.safe)
    if k_conflict >= 0.999:
        return EvidenceMass(0.0, 0.0, 1.0)  # Total conflict fallback

    norm = 1.0 - k_conflict
    safe_comb = ((m1.safe * m2.safe) + (m1.safe * m2.uncertain) + (m1.uncertain * m2.safe)) / norm
    toxic_comb = ((m1.toxic * m2.toxic) + (m1.toxic * m2.uncertain) + (m1.uncertain * m2.toxic)) / norm
    unc_comb = (m1.uncertain * m2.uncertain) / norm

    return EvidenceMass(safe_comb, toxic_comb, unc_comb)


# ============================================================================
# 3. CORE ANALYTICAL MODULES
# ============================================================================

def execute_provenance_audit(mol: Chem.Mol, schema_path: str = "regulatory_rules.json") -> Dict[str, Any]:
    """Scans compound graph against JSON regulatory schema with explicit file status tracking."""
    if not os.path.exists(schema_path):
        return {
            "valid": False,
            "status": "REGULATORY_SCHEMA_UNAVAILABLE",
            "violations": [],
            "warnings": [],
            "details": f"Schema file '{schema_path}' not found."
        }

    try:
        with open(schema_path, "r") as f:
            rules = json.load(f).get("regulatory_rules", [])

        violations, warnings = [], []
        for rule in rules:
            pat = Chem.MolFromSmarts(rule["smarts"])
            if pat and mol.HasSubstructMatch(pat):
                flag = {k: rule[k] for k in ["rule_id", "target_class", "status", "legal_instrument"] if k in rule}
                if "Prohibited" in rule.get("status", ""):
                    violations.append(flag)
                else:
                    warnings.append(flag)

        return {"valid": True, "status": "ACTIVE", "violations": violations, "warnings": warnings}
    except Exception as e:
        return {"valid": False, "status": "SCHEMA_READ_ERROR", "violations": [], "warnings": [], "details": str(e)}


def run_structural_filter(mol: Chem.Mol) -> List[str]:
    """Evaluates molecule against RDKit native BRENK and PAINS catalogs."""
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    return [match.GetDescription() for match in catalog.GetMatches(mol)]


def simulate_skin_metabolism(parent_mol: Chem.Mol) -> List[Dict[str, Any]]:
    """Generates putative bioactivation hypotheses using reaction SMARTS."""
    rxns = {
        "Cutaneous Ester Hydrolysis": "[CX3:1](=[OX1:2])[OX2][C:4]>>[CX3:1](=[OX1:2])[OH].[C:4][OH]",
        "Catechol Oxidation": "[c:1]1[c:2]([OH:7])[c:3]([OH:8])[c:4][c:5][c:6]1>>[C:1]1=[C:6][C:5]=[C:4][C:3](=[O:8])[C:2]1=[O:7]",
        "Allylic Oxidation": "[C:1]=[C:2]-[CH2:3][OH]>>[C:1]=[C:2]-[CH:3]=O"
    }
    records = []
    for name, smarts in rxns.items():
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((parent_mol,))
            for p_tuple in products:
                for met in p_tuple:
                    try:
                        Chem.SanitizeMol(met)
                        alerts = run_structural_filter(met)
                        records.append({
                            "pathway": name,
                            "smiles": Chem.MolToSmiles(met),
                            "secondary_alerts": alerts,
                            "hypothesis": "Putative Bioactivation Alert" if len(alerts) > 0 else "Stable Metabolite"
                        })
                    except Exception:
                        continue
        except Exception:
            continue
    return records


def calculate_probabilistic_mos(contract: AnalysisInputContract, logp: float) -> Dict[str, Any]:
    """Executes scenario-aware Monte Carlo Margin of Safety vectorization with finiteness checks."""
    rng = np.random.default_rng(contract.random_seed)
    scenario = SCCS_PRODUCT_EXPOSURE[contract.product_type]
    n = contract.mc_samples

    target_af = 100.0  # Standard SCCS uncertainty target (10x inter * 10x intra)

    # Lipophilicity-driven variance
    cv_da = max(0.15, 0.40 - (0.05 * abs(logp - 2.5)))
    da_mean = contract.dermal_absorption_pct / 100.0

    bw = rng.normal(contract.body_weight_kg, 10.2, n).clip(40.0, 120.0)
    da = rng.normal(da_mean, da_mean * cv_da, n).clip(0.001, 1.0)
    applied_g = rng.normal(scenario["daily_amount_g"], scenario["sigma_g"], n).clip(0.01, None)
    applied_mg = applied_g * 1000.0 * scenario["retention_factor"]
    pod_dist = rng.normal(contract.pod_noael_mg_kg_day, contract.pod_noael_mg_kg_day * 0.1, n).clip(0.001, None)

    sed = (applied_mg * (contract.concentration_pct / 100.0) * da) / bw
    mos = pod_dist / sed

    # Assert Array Finiteness (Fail Closed on NaN/Inf)
    if np.isnan(mos).any() or np.isinf(mos).any():
        raise ArithmeticError("Monte Carlo execution produced NaN or Inf values in MoS distribution.")

    fail_prob = float(np.sum(mos < target_af) / n)

    return {
        "median_sed_mg_kg_day": round(float(np.median(sed)), 6),
        "median_mos": round(float(np.median(mos)), 1),
        "ci_05_mos": round(float(np.percentile(mos, 5)), 1),
        "failure_probability": fail_prob,
        "target_af": target_af,
        "dynamic_cv_da": round(cv_da, 3)
    }


def execute_similarity_read_across(target_mol: Chem.Mol) -> List[Dict[str, Any]]:
    """Calculates structural and physicochemical similarity scores against reference compounds."""
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
            if not r_mol:
                continue

            tanimoto = DataStructs.TanimotoSimilarity(t_fp, gen.GetFingerprint(r_mol))
            dist_logp = min(abs(t_logp - Descriptors.MolLogP(r_mol)) / 5.0, 1.0)
            dist_mw = min(abs(t_mw - Descriptors.MolWt(r_mol)) / 300.0, 1.0)

            similarity_score = max((1.0 - (((1.0 - tanimoto) * 0.6) + (dist_logp * 0.3) + (dist_mw * 0.1))) * 100.0, 0.0)
            cands.append({
                "name": r["name"],
                "similarity_score": round(similarity_score, 1),
                "tanimoto": round(tanimoto * 100, 1)
            })

        return sorted(cands, key=lambda x: x["similarity_score"], reverse=True)
    except Exception:
        return []


# ============================================================================
# 4. MASTER ORCHESTRATION PIPELINE
# ============================================================================

def execute_full_compound_audit(
    smiles: str,
    product_type: str,
    concentration_pct: float,
    pod_noael_mg_kg_day: float,
    dermal_absorption_pct: float = 50.0,
    body_weight_kg: float = 60.0,
    mc_samples: int = 10000,
    random_seed: int = 42
) -> Dict[str, Any]:
    """Executes strict contract validation, stochastic exposure, and discounted Dempster-Shafer fusion."""
    
    # 1. Enforce Input Contract
    try:
        contract = AnalysisInputContract(
            smiles=smiles,
            product_type=product_type,
            concentration_pct=concentration_pct,
            pod_noael_mg_kg_day=pod_noael_mg_kg_day,
            dermal_absorption_pct=dermal_absorption_pct,
            body_weight_kg=body_weight_kg,
            mc_samples=mc_samples,
            random_seed=random_seed
        )
    except ValidationError as e:
        return {"valid": False, "error": str(e)}

    # 2. Canonicalization
    mol = Chem.MolFromSmiles(contract.smiles)
    if not mol:
        return {"valid": False, "error": f"RDKit SMILES parsing failure: {contract.smiles!r}"}

    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    logp = Descriptors.MolLogP(mol)

    # 3. Compute Analytics
    reg = execute_provenance_audit(mol)
    alerts = run_structural_filter(mol)
    met = simulate_skin_metabolism(mol)
    mos = calculate_probabilistic_mos(contract, logp)
    analogs = execute_similarity_read_across(mol)

    # 4. Construct Evidence Masses with Shafer Discounting
    # Primary Exposure Mass (Independent)
    exp_safe_prob = (1.0 - mos["failure_probability"]) * 0.90
    m_mos = EvidenceMass(
        safe=exp_safe_prob,
        toxic=mos["failure_probability"],
        uncertain=1.0 - exp_safe_prob - mos["failure_probability"]
    )

    # Structural Rules Mass (Discounted alpha = 0.50 due to graph correlation)
    if reg.get("violations"):
        m_struct_raw = EvidenceMass(safe=0.0, toxic=0.85, uncertain=0.15)
    elif alerts:
        m_struct_raw = EvidenceMass(safe=0.0, toxic=0.65, uncertain=0.35)
    else:
        m_struct_raw = EvidenceMass(safe=0.50, toxic=0.0, uncertain=0.50)

    m_struct = m_struct_raw.apply_discounting(alpha=0.50)

    # Similarity Read-Across Mass (Discounted alpha = 0.30 due to graph overlap)
    top_sim = (analogs[0]["similarity_score"] / 100.0) if analogs else 0.0
    m_ana_raw = EvidenceMass(safe=top_sim * 0.80, toxic=0.0, uncertain=1.0 - (top_sim * 0.80))
    m_ana = m_ana_raw.apply_discounting(alpha=0.30)

    # Fused Evidence
    fused_mass = combine_evidence(combine_evidence(m_mos, m_struct), m_ana)

    # 5. Deterministic Assessment Hash ID
    hash_payload = (
        f"{canonical_smiles}|{contract.product_type}|{contract.concentration_pct}|"
        f"{contract.pod_noael_mg_kg_day}|{contract.dermal_absorption_pct}|{contract.random_seed}"
    )
    assessment_id = f"IATA-V4-{hashlib.sha256(hash_payload.encode()).hexdigest()[:10].upper()}"

    return {
        "valid": True,
        "assessment_id": assessment_id,
        "mol": mol,
        "canonical_smiles": canonical_smiles,
        "contract": contract,
        "dst_metrics": {
            "belief_safe": round(fused_mass.safe * 100, 1),
            "belief_toxic": round(fused_mass.toxic * 100, 1),
            "epistemic_uncertainty": round(fused_mass.uncertain * 100, 1),
            "plausibility_safe": round((fused_mass.safe + fused_mass.uncertain) * 100, 1)
        },
        "reg": reg,
        "alerts": alerts,
        "met": met,
        "mos": mos,
        "analogs": analogs
    }


# ============================================================================
# 5. PAGINATED REPORTLAB PDF GENERATOR
# ============================================================================

class NumberedCanvas(canvas.Canvas):
    """Custom canvas handling dynamic pagination and legal disclaimers."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.saveState()
            self.setFont("Helvetica", 8)
            self.setFillColor(colors.HexColor("#718096"))
            self.drawString(36, 20, "Research-Use Computational Assessment | Non-Standalone Dossier")
            self.drawRightString(576, 20, f"Page {self._pageNumber} of {num_pages}")
            self.restoreState()
            super().showPage()
        super().save()


def generate_enterprise_pdf(audit_data: Dict[str, Any]) -> BytesIO:
    """Renders a PDF ledger with custom pagination and XML escaping."""
    if not audit_data.get("valid", False):
        raise ValueError(f"Cannot render PDF for invalid audit: {audit_data.get('error')}")

    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter, rightMargin=36, leftMargin=36, topMargin=36, bottomMargin=36)
    styles = getSampleStyleSheet()

    h1 = ParagraphStyle("H1", parent=styles["Normal"], fontName="Helvetica-Bold", fontSize=11, textColor=colors.HexColor("#1B365D"), spaceBefore=8, spaceAfter=4)
    cell_bold = ParagraphStyle("CB", parent=styles["Normal"], fontSize=8, fontName="Helvetica-Bold", leading=10)
    cell_norm = ParagraphStyle("CN", parent=styles["Normal"], fontSize=8, leading=10)

    def _esc(txt: Any) -> str:
        return html.escape(str(txt))

    img_buf = BytesIO()
    Draw.MolToImage(audit_data["mol"], size=(200, 200)).save(img_buf, format="PNG")
    img_buf.seek(0)

    dst = audit_data["dst_metrics"]
    summary_html = f"""
    <b>Assessment ID:</b> {_esc(audit_data['assessment_id'])}<br/>
    <b>Canonical SMILES:</b> {_esc(audit_data['canonical_smiles'])}<br/>
    <b>Belief (Safe):</b> {dst['belief_safe']}% | <b>Belief (Toxic):</b> {dst['belief_toxic']}%<br/>
    <b>Epistemic Uncertainty:</b> {dst['epistemic_uncertainty']}% | <b>Plausibility (Safe):</b> {dst['plausibility_safe']}%
    """

    story = [
        Paragraph("COMPUTATIONAL TOXICOLOGY & ELEVATED EVIDENCE DOSSIER", ParagraphStyle("T", fontName="Helvetica-Bold", fontSize=16, leading=20, textColor=colors.HexColor("#1B365D"))),
        Paragraph("IATA-Aligned Probabilistic Methodological Framework", styles["Normal"]),
        HRFlowable(width="100%", thickness=1.5, color=colors.HexColor("#1B365D"), spaceAfter=8),
        Table([[PlatypusImage(img_buf, width=100, height=100), Paragraph(summary_html, styles["Normal"])]], colWidths=[110, 430], style=[
            ('BACKGROUND', (0,0), (-1,-1), colors.HexColor('#F7FAFC')),
            ('BOX', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E0')),
            ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
            ('PADDING', (0,0), (-1,-1), 6)
        ]),
        Spacer(1, 8),
        Paragraph("1. Stochastic Scenario Margin of Safety (Monte Carlo Analysis)", h1),
        Table([
            [Paragraph("Median SED (mg/kg/d)", cell_bold), Paragraph("Median MoS", cell_bold), Paragraph("5th Percentile MoS", cell_bold), Paragraph("Failure Probability", cell_bold)],
            [Paragraph(str(audit_data["mos"]["median_sed_mg_kg_day"]), cell_norm), Paragraph(str(audit_data["mos"]["median_mos"]), cell_norm), Paragraph(str(audit_data["mos"]["ci_05_mos"]), cell_norm), Paragraph(f"{audit_data['mos']['failure_probability']*100:.2f}%", cell_norm)]
        ], colWidths=[135, 135, 135, 135], style=[('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E0')), ('PADDING', (0,0), (-1,-1), 4)])
    ]

    doc.build(story, canvasmaker=NumberedCanvas)
    buffer.seek(0)
    return buffer
            
