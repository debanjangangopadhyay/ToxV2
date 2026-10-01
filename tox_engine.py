import os
import json
import math
import hashlib
from io import BytesIO
from dataclasses import dataclass
from typing import Dict, List, Any, Optional

import numpy as np
from PIL import Image as PILImage
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import (
    HRFlowable, Paragraph, SimpleDocTemplate, Spacer, Table, Image, PageBreak, KeepTogether
)
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
    pod_type: str = "NOAEL"
    species: str = "Rat"
    duration: str = "Subchronic"
    dermal_absorption_pct: float = 50.0  
    body_weight_kg: float = 60.0
    mc_samples: int = 10000
    random_seed: int = 42

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
    if k_conflict >= 0.999: 
        return EvidenceMass(0.0, 0.0, 1.0)
    norm = 1.0 - k_conflict
    return EvidenceMass(
        ((m1.safe * m2.safe) + (m1.safe * m2.uncertain) + (m1.uncertain * m2.safe)) / norm,
        ((m1.toxic * m2.toxic) + (m1.toxic * m2.uncertain) + (m1.uncertain * m2.toxic)) / norm,
        (m1.uncertain * m2.uncertain) / norm
    )

def render_2d_molecule_png(mol: Chem.Mol, width: int = 350, height: int = 350) -> BytesIO:
    """Renders high-resolution 2D structure image buffer with Cairo/Pillow fallback."""
    img_buffer = BytesIO()
    try:
        img = Draw.MolToImage(mol, size=(width, height))
        img.save(img_buffer, format="PNG")
    except Exception:
        # Drawing Fallback Mode
        pil_img = PILImage.new('RGB', (width, height), color=(255, 255, 255))
        pil_img.save(img_buffer, format="PNG")
    img_buffer.seek(0)
    return img_buffer

def run_structural_filter(mol: Chem.Mol) -> List[str]:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return [match.GetDescription() for match in FilterCatalog(params).GetMatches(mol)]

def calculate_probabilistic_mos(contract: AnalysisInputContract, logp: float) -> Dict[str, Any]:
    rng = np.random.default_rng(contract.random_seed)
    scenario = SCCS_PRODUCT_EXPOSURE[contract.product_type]
    n = contract.mc_samples

    af_inter = 1.0 if contract.species.lower() == "human" else 10.0
    af_intra = 10.0 
    duration_lower = contract.duration.lower()
    if duration_lower == "chronic": af_duration = 1.0
    elif duration_lower == "subchronic": af_duration = 3.0
    else: af_duration = 10.0
        
    pod_lower = contract.pod_type.lower()
    if pod_lower == "loael": af_pod = 3.0
    elif "bmdl" in pod_lower: af_pod = 1.0
    else: af_pod = 1.0
        
    target_af = float(af_inter * af_intra * af_duration * af_pod)

    da_max = max(contract.dermal_absorption_pct / 100.0, 0.0001)
    target_mu = 0.70  
    cv_da = 0.30 
    target_var = (target_mu * cv_da) ** 2

    phi = max((target_mu * (1.0 - target_mu) / target_var) - 1.0, 0.001)
    alpha_beta = target_mu * phi
    beta_beta = (1.0 - target_mu) * phi

    da_raw = rng.beta(alpha_beta, beta_beta, n)
    da = da_raw * da_max

    bw = rng.normal(contract.body_weight_kg, 10.2, n).clip(40.0, 120.0)
    applied_g = rng.normal(scenario["daily_amount_g"], scenario["sigma_g"], n).clip(0.01, None)
    applied_mg = applied_g * 1000.0 * scenario["retention_factor"]
    pod_dist = rng.normal(contract.pod_noael_mg_kg_day, contract.pod_noael_mg_kg_day * 0.1, n).clip(0.001, None)

    sed = np.clip((applied_mg * (contract.concentration_pct / 100.0) * da) / bw, 1e-9, None)
    mos = pod_dist / sed

    fail_prob = float(np.sum(mos < target_af) / n)

    return {
        "median_sed_mg_kg_day": float(round(float(np.median(sed)), 6)),
        "median_mos": float(round(float(np.median(mos)), 1)),
        "ci_05_mos": float(round(float(np.percentile(mos, 5)), 1)),
        "failure_probability": fail_prob,
        "target_af": target_af,
        "af_breakdown": {
            "interspecies": float(af_inter), 
            "intraspecies": float(af_intra), 
            "duration": float(af_duration), 
            "data_quality": float(af_pod)
        }
    }

def execute_full_compound_audit(
    smiles: str, 
    product_type: str, 
    concentration_pct: float, 
    pod_noael_mg_kg_day: float, 
    pod_type: str = "NOAEL", 
    species: str = "Rat", 
    duration: str = "Subchronic", 
    dermal_absorption_pct: float = 50.0, 
    body_weight_kg: float = 60.0, 
    mc_samples: int = 10000, 
    random_seed: int = 42,
    baseline_safe_belief: float = 0.70,
    structural_discount_rate: float = 0.50
) -> Dict[str, Any]:
    try:
        contract = AnalysisInputContract(
            smiles, product_type, concentration_pct, pod_noael_mg_kg_day, 
            pod_type, species, duration, dermal_absorption_pct, body_weight_kg, 
            mc_samples, random_seed
        )
    except ValidationError as e:
        return {"valid": False, "error": str(e)}

    raw_mol = Chem.MolFromSmiles(contract.smiles)
    if not raw_mol: 
        return {"valid": False, "error": "Invalid SMILES structure provided."}
    
    mol = rdMolStandardize.LargestFragmentChooser().choose(raw_mol)
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    logp = float(Descriptors.MolLogP(mol))

    mos = calculate_probabilistic_mos(contract, logp)
    exp_safe_prob = (1.0 - mos["failure_probability"]) * 0.90
    m_mos = EvidenceMass(
        safe=exp_safe_prob, 
        toxic=mos["failure_probability"], 
        uncertain=1.0 - exp_safe_prob - mos["failure_probability"]
    )
    
    alerts = run_structural_filter(mol)
    alert_count = len(alerts)
    if alert_count > 0:
        toxic_mass = min(0.85, 0.40 + (0.15 * alert_count))
        m_struct = EvidenceMass(safe=0.0, toxic=toxic_mass, uncertain=1.0-toxic_mass).apply_discounting(structural_discount_rate)
    else:
        m_struct = EvidenceMass(safe=baseline_safe_belief, toxic=0.0, uncertain=1.0-baseline_safe_belief).apply_discounting(structural_discount_rate)

    fused_mass = combine_evidence(m_mos, m_struct)
    
    hash_payload = f"{canonical_smiles}_{contract.product_type}_{contract.concentration_pct}_{contract.dermal_absorption_pct}"
    
    return {
        "valid": True,
        "assessment_id": f"IATA-V5-{hashlib.sha256(hash_payload.encode()).hexdigest()[:10].upper()}",
        "mol": mol,
        "canonical_smiles": canonical_smiles,
        "contract": contract,
        "dst_metrics": {
            "belief_safe": float(round(fused_mass.safe * 100.0, 1)),
            "belief_toxic": float(round(fused_mass.toxic * 100.0, 1)),
            "epistemic_uncertainty": float(round(fused_mass.uncertain * 100.0, 1)),
            "plausibility_safe": float(round((fused_mass.safe + fused_mass.uncertain) * 100.0, 1))
        },
        "alerts": alerts, 
        "mos": mos
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
            
            # Running Header
            self.setFont("Helvetica", 8)
            self.setFillColor(colors.HexColor("#718096"))
            self.drawString(36, 760, "IATA Computational Toxicology Dossier")
            self.setStrokeColor(colors.HexColor("#E2E8F0"))
            self.setLineWidth(0.5)
            self.line(36, 752, 576, 752)
            
            # Running Footer
            self.line(36, 35, 576, 35)
            self.drawString(36, 22, "CONFIDENTIAL & PROPRIETARY — RESEARCH & PRE-VALIDATION USE ONLY")
            self.drawRightString(576, 22, f"Page {self._pageNumber} of {num_pages}")
            self.restoreState()
            super().showPage()
        super().save()


def generate_enterprise_pdf(audit_data: Dict[str, Any], diep_data: Optional[Dict[str, Any]] = None) -> BytesIO:
    buffer = BytesIO()
    doc = SimpleDocTemplate(
        buffer, 
        pagesize=letter, 
        rightMargin=36, 
        leftMargin=36, 
        topMargin=54, 
        bottomMargin=54
    )
    styles = getSampleStyleSheet()
    
    title_style = ParagraphStyle("TitleStyle", fontName="Helvetica-Bold", fontSize=16, leading=20, textColor=colors.HexColor("#1B365D"))
    h1 = ParagraphStyle("H1", fontName="Helvetica-Bold", fontSize=11, leading=14, textColor=colors.HexColor("#1B365D"), spaceBefore=12, spaceAfter=6)
    h2 = ParagraphStyle("H2", fontName="Helvetica-Bold", fontSize=9, leading=12, textColor=colors.HexColor("#2B6CB0"), spaceBefore=8, spaceAfter=4)
    cell_bold = ParagraphStyle("CB", fontName="Helvetica-Bold", fontSize=8, leading=10, textColor=colors.HexColor("#2D3748"))
    cell_norm = ParagraphStyle("CN", fontName="Helvetica", fontSize=8, leading=10, textColor=colors.HexColor("#4A5568"))
    code_style = ParagraphStyle("Code", fontName="Courier", fontSize=7, leading=9, textColor=colors.HexColor("#2C5282"))

    story = []

    # Section 1: Compound Identification
    story.append(Paragraph("COMPUTATIONAL TOXICOLOGY & BIOAVAILABILITY DOSSIER", title_style))
    story.append(HRFlowable(width="100%", thickness=2, color=colors.HexColor("#1B365D"), spaceAfter=10))
    
    contract = audit_data["contract"]
    
    summary_table_data = [
        [Paragraph("Assessment ID", cell_bold), Paragraph(audit_data["assessment_id"], cell_norm), Paragraph("Date Generated", cell_bold), Paragraph("October 2026", cell_norm)],
        [Paragraph("Target SMILES", cell_bold), Paragraph(audit_data["canonical_smiles"], code_style), Paragraph("Product Scenario", cell_bold), Paragraph(contract.product_type, cell_norm)],
        [Paragraph("Concentration (%)", cell_bold), Paragraph(f"{contract.concentration_pct:.2f}%", cell_norm), Paragraph("PoD Benchmark", cell_bold), Paragraph(f"{contract.pod_noael_mg_kg_day} mg/kg/day ({contract.pod_type})", cell_norm)]
    ]
    story.append(Table(summary_table_data, colWidths=[100, 170, 100, 170], style=[
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E0')),
        ('BACKGROUND', (0,0), (0,-1), colors.HexColor('#F7FAFC')),
        ('BACKGROUND', (2,0), (2,-1), colors.HexColor('#F7FAFC')),
        ('PADDING', (0,0), (-1,-1), 4)
    ]))
    story.append(Spacer(1, 10))

    mol_png_buffer = render_2d_molecule_png(audit_data["mol"], 350, 350)
    mol_image = Image(mol_png_buffer, width=180, height=180)
    
    structure_table = Table([
        [mol_image, Paragraph(
            f"<b>2D Chemical Structure Graph</b><br/><br/>"
            f"<b>Molecular Weight:</b> {Descriptors.MolWt(audit_data['mol']):.2f} g/mol<br/>"
            f"<b>Calculated LogP:</b> {Descriptors.MolLogP(audit_data['mol']):.2f}<br/>"
            f"<b>Heavy Atom Count:</b> {audit_data['mol'].GetNumHeavyAtoms()}<br/>"
            f"<b>Rotatable Bonds:</b> {Descriptors.NumRotatableBonds(audit_data['mol'])}<br/>"
            f"<b>Aromatic Rings:</b> {Descriptors.NumAromaticRings(audit_data['mol'])}", 
            cell_norm
        )]
    ], colWidths=[200, 340], style=[
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#E2E8F0')),
        ('PADDING', (0,0), (-1,-1), 6)
    ])
    story.append(structure_table)
    story.append(Spacer(1, 14))

    # Section 2: Deterministic Biophysics & Cramer Decision Tree
    story.append(Paragraph("1. Deterministic Biophysical & Transdermal Flux Bounds (DIEP-MoS)", h1))
    if diep_data:
        diep_table_data = [
            [Paragraph("Unionized Fraction ($f_{ui}$)", cell_bold), Paragraph(f"{diep_data['f_ui']:.4f}", cell_norm), Paragraph("Transdermal Absorption", cell_bold), Paragraph(f"{diep_data['da_pct_applied']:.2f}%", cell_norm)],
            [Paragraph("Systemic Dose (SED)", cell_bold), Paragraph(f"{diep_data['sed_ug_day']:.2f} µg/day", cell_norm), Paragraph("EFSA TTC Limit", cell_bold), Paragraph(f"{diep_data['ttc_limit_ug']:.2f} µg/day", cell_norm)],
            [Paragraph("Cramer Topological Class", cell_bold), Paragraph(diep_data['cramer_class'], cell_norm), Paragraph("Systemic Gate Status", cell_bold), Paragraph(diep_data['status'], cell_bold)]
        ]
        story.append(Table(diep_table_data, colWidths=[135, 135, 135, 135], style=[
            ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E0')),
            ('PADDING', (0,0), (-1,-1), 4)
        ]))
        
        story.append(Paragraph("Computational Cramer Decision Tree Traversal Log:", h2))
        log_paragraphs = [Paragraph(f"• {log_entry}", cell_norm) for log_entry in diep_data.get("cramer_tree_log", [])]
        story.append(KeepTogether(log_paragraphs))
    else:
        story.append(Paragraph("Deterministic flux module disabled during audit run.", cell_norm))

    story.append(Spacer(1, 14))

    # Section 3: Monte Carlo Probabilistic Exposure
    story.append(Paragraph("2. Stochastic Probabilistic Risk Engine (Monte Carlo Simulation)", h1))
    mos_data = audit_data["mos"]
    mc_table_data = [
        [Paragraph("Simulation Iterations", cell_bold), Paragraph(f"{contract.mc_samples:,}", cell_norm), Paragraph("Target Assessment Factor", cell_bold), Paragraph(str(mos_data["target_af"]), cell_norm)],
        [Paragraph("Median SED (mg/kg/day)", cell_bold), Paragraph(str(mos_data["median_sed_mg_kg_day"]), cell_norm), Paragraph("Median Margin of Safety", cell_bold), Paragraph(str(mos_data["median_mos"]), cell_norm)],
        [Paragraph("5th Percentile MoS", cell_bold), Paragraph(str(mos_data["ci_05_mos"]), cell_norm), Paragraph("Formulation Failure Prob.", cell_bold), Paragraph(f"{mos_data['failure_probability']*100:.2f}%", cell_norm)]
    ]
    story.append(Table(mc_table_data, colWidths=[135, 135, 135, 135], style=[
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E0')),
        ('PADDING', (0,0), (-1,-1), 4)
    ]))
    story.append(Spacer(1, 14))

    # Section 4: Dempster-Shafer Evidence Fusion
    story.append(Paragraph("3. Dempster-Shafer Evidentiary Mass Fusion", h1))
    dst = audit_data["dst_metrics"]
    dst_table_data = [
        [Paragraph("Belief (Safe)", cell_bold), Paragraph("Belief (Toxic)", cell_bold), Paragraph("Epistemic Uncertainty", cell_bold), Paragraph("Plausibility (Safe)", cell_bold)],
        [Paragraph(f"{dst['belief_safe']}%", cell_norm), Paragraph(f"{dst['belief_toxic']}%", cell_norm), Paragraph(f"{dst['epistemic_uncertainty']}%", cell_norm), Paragraph(f"{dst['plausibility_safe']}%", cell_norm)]
    ]
    story.append(Table(dst_table_data, colWidths=[135, 135, 135, 135], style=[
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E0')),
        ('PADDING', (0,0), (-1,-1), 4)
    ]))
    story.append(Spacer(1, 14))

    # Section 5: TRACE-Onco Integration
    story.append(Paragraph("4. TRACE-Onco Hepatic Stress & Risk Vector", h1))
    if diep_data:
        trace_table_data = [
            [Paragraph("Hepatic Burden Ratio (HBR)", cell_bold), Paragraph(f"{diep_data['hepatic_burden_ratio']:.4f}", cell_norm)],
            [Paragraph("Oncogenic Risk Index", cell_bold), Paragraph(f"{diep_data['oncogenic_risk_index']:.4f}", cell_norm)],
            [Paragraph("Structural Toxicophore Alerts", cell_bold), Paragraph(", ".join(audit_data["alerts"]) if audit_data["alerts"] else "None Detected", cell_norm)]
        ]
        story.append(Table(trace_table_data, colWidths=[180, 290], style=[
            ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E0')),
            ('PADDING', (0,0), (-1,-1), 4)
        ]))

    doc.build(story, canvasmaker=NumberedCanvas)
    buffer.seek(0)
    return buffer
    
