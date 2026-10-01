import math
from typing import Dict, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

class ToxicophoreMatchException(Exception): pass
class PhysicochemicalConstraintError(ValueError): pass

# --- 1. SCCS Mechanistic Product Registry ---
SCCS_MECHANISTIC_REGISTRY = {
    "Lip Balm (Leave-on)": (57.0, 4.8, 1.0),
    "Body Lotion (Leave-on)": (7820.0, 15670.0, 1.0),
    "Face Cream (Leave-on)": (1540.0, 565.0, 1.0),
}

TOXICOPHORES = {
    "Tropane_Alkaloid_Core": "[C]1[C][C]2[C][C]1[N]2", 
    "Potent_Anesthetic_Ester": "O=C(OCCN(CC)CC)c1ccc(N)cc1" 
}

def standardize_api(smiles: str) -> Tuple[float, float, Chem.Mol]:
    """Isolates the Active Pharmaceutical Ingredient (API) by stripping salts/counterions."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError("Invalid SMILES.")
    clean_mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
    return Descriptors.MolWt(clean_mol), Descriptors.MolLogP(clean_mol), clean_mol

def calculate_dual_pathway_flux(mw: float, logp: float, ph: float, pka: float, is_base: bool, conc_pct: float, product_type: str) -> float:
    """
    Computes flux (mg/cm^2/hr) accounting for both unionized (lipophilic) and ionized (appendageal/pore) pathways.
    """
    # 1. Ionization State (Henderson-Hasselbalch)
    f_ui = 1 / (1 + 10**(pka - ph)) if is_base else 1 / (1 + 10**(ph - pka))
    f_i = 1.0 - f_ui
    
    # 2. Base Permeability (Potts & Guy for Stratum Corneum)
    kp_ui_cm_hr = 10 ** ((0.71 * logp) - (0.0061 * mw) - 2.72)
    
    # Ionized species permeate approx 100x slower through intact skin (appendageal route)
    kp_i_cm_hr = kp_ui_cm_hr * 0.01 

    # 3. Mucosal Membrane Permeability Dynamics
    if "lip" in product_type.lower() or "oral" in product_type.lower():
        # Mucosa lacks stratum corneum; ionized/aqueous pore transport increases significantly
        kp_ui_cm_hr *= 10.0   # Lipophilic pathway resistance drops
        kp_i_cm_hr *= 100.0   # Aqueous pore pathway opens dramatically
        
    # 4. Volumetric Vehicle Concentration (assuming aqueous/lipid emulsion density ~ 1.0 g/cm^3)
    c_vehicle_mg_cm3 = (conc_pct / 100.0) * 1000.0
    
    # 5. Dual-Pathway Fickian Flux (J = Kp * C)
    flux_ui = kp_ui_cm_hr * (c_vehicle_mg_cm3 * f_ui)
    flux_i = kp_i_cm_hr * (c_vehicle_mg_cm3 * f_i)
    
    return flux_ui + flux_i

def execute_vmtb_safety_gate(smiles: str, conc_pct: float, ph: float, pka: float, is_base: bool, product_type: str, noael: float, body_weight: float) -> Dict:
    mw, logp, clean_mol = standardize_api(smiles)
    
    for alert, smarts in TOXICOPHORES.items():
        if clean_mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            raise ToxicophoreMatchException(f"Critical Systemic Hazard: {alert} detected.")
            
    amount_mg, area_cm2, retention = SCCS_MECHANISTIC_REGISTRY[product_type]
    total_flux = calculate_dual_pathway_flux(mw, logp, ph, pka, is_base, conc_pct, product_type)
    
    # Total Systemic Exposure Dose (SED)
    absorbed_mass_mg = total_flux * area_cm2 * 24.0 * retention
    sed_mg_kg_day = absorbed_mass_mg / body_weight
    
    # Deterministic Margin of Safety
    mos = noael / sed_mg_kg_day if sed_mg_kg_day > 0 else float('inf')
    
    # Systemic Bioavailability Scalar (For TRACE-Onco integration)
    # Scales 0 to >1.0. A score >1.0 means the patient's exposure exceeds the 100x safety factor.
    s_bio_scalar = 100.0 / mos if mos > 0 else float('inf')
    
    return {
        "api_metrics": {"MW": round(mw, 2), "LogP": round(logp, 2), "f_ui": round(1 / (1 + 10**(pka - ph)) if is_base else 1 / (1 + 10**(ph - pka)), 4)},
        "sed_mg_kg_day": round(sed_mg_kg_day, 6),
        "deterministic_mos": round(mos, 2),
        "trace_onco_payload": {
            "systemic_bioavailability_scalar": round(s_bio_scalar, 4),
            "clinical_verdict": "SUITABLE" if mos >= 100 else "REJECTED: SYSTEMIC TOXICITY"
        }
    }
    
