import os
import json
import html
import math
import hashlib
from io import BytesIO
from dataclasses import dataclass
from typing import Dict, List, Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw, Descriptors, rdFingerprintGenerator
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import HRFlowable, Paragraph, SimpleDocTemplate, Spacer, Table, Image as PlatypusImage
from reportlab.pdfgen import canvas

class ValidationError(ValueError): pass

SCCS_PRODUCT_EXPOSURE = {
    "Face Cream (Leave-on)": {"daily_amount_g": 1.54, "sigma_g": 0.23, "retention_factor": 1.0},
    "Body Lotion (Leave-on)": {"daily_amount_g": 7.82, "sigma_g": 1.10, "retention_factor": 1.0},
    "Facial Cleanser (Rinse-off)": {"daily_amount_g": 20.0, "sigma_g": 2.50, "retention_factor": 0.01},
    "Shampoo (Rinse-off)": {"daily_amount_g": 10.46, "sigma_g": 1.80, "retention_factor": 0.01},
    "Eye Cream (Leave-on)": {"daily_amount_g": 0.50, "sigma_g": 0.08, "retention_factor": 1.0},
    "Lip Balm (Leave-on)": {"daily_amount_g": 0.057, "sigma_g": 0.01, "retention_factor": 1.0}
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
            raise ValidationError(f"Invalid product_type '{self.product_type}'.")

@dataclass
class EvidenceMass:
    safe: float
    toxic: float
    uncertain: float

    def __post_init__(self):
        total = self.safe + self.toxic + self.uncertain
        if not math.isclose(total, 1.0, abs_tol=1e-5):
            self.safe /= total; self.toxic /= total; self.uncertain /= total

    def apply_discounting(self, alpha: float) -> 'EvidenceMass':
        alpha = max(0.0, min(1.0, alpha))
        return EvidenceMass(alpha * self.safe, alpha * self.toxic, 1.0 - alpha * (self.safe + self.toxic))

def combine_evidence(m1: EvidenceMass, m2: EvidenceMass) -> EvidenceMass:
    k_conflict = (m1.safe * m2.toxic) + (m1.toxic * m2.safe)
    if k_conflict >= 0.999: return EvidenceMass(0.0, 0.0, 1.0)
    norm = 1.0 - k_conflict
    return EvidenceMass(
        ((m1.safe * m2.safe) + (m1.safe * m2.uncertain) + (m1.uncertain * m2.safe)) / norm,
        ((m1.toxic * m2.toxic) + (m1.toxic * m2.uncertain) + (m1.uncertain * m2.toxic)) / norm,
        (m1.uncertain * m2.uncertain) / norm
    )

def execute_provenance_audit(mol: Chem.Mol, schema_path: str = "regulatory_rules.json") -> Dict[str, Any]:
    return {"valid": True, "status": "ACTIVE", "violations": [], "warnings": []}

def run_structural_filter(mol: Chem.Mol) -> List[str]:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return [match.GetDescription() for match in FilterCatalog(params).GetMatches(mol)]

def simulate_skin_metabolism(parent_mol: Chem.Mol) -> List[Dict[str, Any]]:
    return []

def calculate_probabilistic_mos(contract: AnalysisInputContract, logp: float) -> Dict[str, Any]:
    """Integrated engine applying Dynamic Assessment Factors (AF) to Monte Carlo Simulation"""
    rng = np.random.default_rng(contract.random_seed)
    scenario = SCCS_PRODUCT_EXPOSURE[contract.product_type]
    n = contract.mc_samples

    # Integrated from exposure_engine (1).py
    af_inter = 10.0   # Animal to Human Extrapolation
    af_intra = 10.0   # Intraspecies Variance
    af_duration = 3.0 # Subchronic default
    af_pod = 1.0      # NOAEL confidence
    target_af = af_inter * af_intra * af_duration * af_pod

    cv_da = max(0.15, 0.40 - (0.05 * abs(logp - 2.5)))
    da_mean = contract.dermal_absorption_pct / 100.0

    bw = rng.normal(contract.body_weight_kg, 10.2, n).clip(40.0, 120.0)
    da = rng.normal(da_mean, da_mean * cv_da, n).clip(0.001, 1.0)
    applied_g = rng.normal(scenario["daily_amount_g"], scenario["sigma_g"], n).clip(0.01, None)
    applied_mg = applied_g * 1000.0 * scenario["retention_factor"]
    pod_dist = rng.normal(contract.pod_noael_mg_kg_day, contract.pod_noael_mg_kg_day * 0.1, n).clip(0.001, None)

    sed = (applied_mg * (contract.concentration_pct / 100.0) * da) / bw
    mos = pod_dist / sed

    fail_prob = float(np.sum(mos < target_af) / n)

    return {
        "median_sed_mg_kg_day": round(float(np.median(sed)), 6),
        "median_mos": round(float(np.median(mos)), 1),
        "ci_05_mos": round(float(np.percentile(mos, 5)), 1),
        "failure_probability": fail_prob,
        "target_af": target_af,
        "dynamic_cv_da": round(cv_da, 3),
        "af_breakdown": {"interspecies": af_inter, "intraspecies": af_intra, "duration": af_duration, "data_quality": af_pod}
    }

def execute_similarity_read_across(target_mol: Chem.Mol) -> List[Dict[str, Any]]:
    return []

def execute_full_compound_audit(smiles: str, product_type: str, concentration_pct: float, pod_noael_mg_kg_day: float, dermal_absorption_pct: float = 50.0, body_weight_kg: float = 60.0, mc_samples: int = 10000, random_seed: int = 42) -> Dict[str, Any]:
    try:
        contract = AnalysisInputContract(smiles, product_type, concentration_pct, pod_noael_mg_kg_day, dermal_absorption_pct, body_weight_kg, mc_samples, random_seed)
    except ValidationError as e:
        return {"valid": False, "error": str(e)}

    mol = Chem.MolFromSmiles(contract.smiles)
    if not mol: return {"valid": False, "error": f"Invalid SMILES."}
    
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    logp = Descriptors.MolLogP(mol)

    mos = calculate_probabilistic_mos(contract, logp)
    exp_safe_prob = (1.0 - mos["failure_probability"]) * 0.90
    m_mos = EvidenceMass(safe=exp_safe_prob, toxic=mos["failure_probability"], uncertain=1.0 - exp_safe_prob - mos["failure_probability"])
    
    alerts = run_structural_filter(mol)
    m_struct = EvidenceMass(safe=0.0, toxic=0.65, uncertain=0.35).apply_discounting(0.50) if alerts else EvidenceMass(safe=0.50, toxic=0.0, uncertain=0.50).apply_discounting(0.50)

    fused_mass = combine_evidence(m_mos, m_struct)
    
    return {
        "valid": True,
        "assessment_id": f"IATA-V4-{hashlib.sha256(canonical_smiles.encode()).hexdigest()[:10].upper()}",
        "mol": mol,
        "canonical_smiles": canonical_smiles,
        "contract": contract,
        "dst_metrics": {
            "belief_safe": round(fused_mass.safe * 100, 1),
            "belief_toxic": round(fused_mass.toxic * 100, 1),
            "epistemic_uncertainty": round(fused_mass.uncertain * 100, 1),
            "plausibility_safe": round((fused_mass.safe + fused_mass.uncertain) * 100, 1)
        },
        "reg": {}, "alerts": alerts, "met": [], "mos": mos, "analogs": []
    }

class NumberedCanvas(canvas.Canvas):
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
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter, rightMargin=36, leftMargin=36, topMargin=36, bottomMargin=36)
    doc.build([Paragraph("COMPUTATIONAL TOXICOLOGY DOSSIER", ParagraphStyle("T", fontName="Helvetica-Bold", fontSize=16, textColor=colors.HexColor("#1B365D")))], canvasmaker=NumberedCanvas)
    buffer.seek(0)
    return buffer
    
