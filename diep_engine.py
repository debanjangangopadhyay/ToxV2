import math
from typing import Dict, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

class ToxicophoreMatchException(Exception): pass
class PhysicochemicalConstraintError(ValueError): pass

# --- 1. SCCS Mechanistic Product Registry ---
# Maps product to (Amount applied mg/day, Surface Area cm2, Retention Factor)
SCCS_MECHANISTIC_REGISTRY = {
    "Lip Balm (Leave-on)": (57.0, 4.8, 1.0),
    "Body Lotion (Leave-on)": (7820.0, 15670.0, 1.0),
    "Face Cream (Leave-on)": (1540.0, 565.0, 1.0),
    "Deodorant (Leave-on)": (1500.0, 200.0, 1.0),
    "Shampoo (Rinse-off)": (10460.0, 1440.0, 0.01)
}

EFSA_TTC_SYSTEMIC_UG_DAY = {
    "Cramer_I": 1800.0 * 0.8,
    "Cramer_II": 540.0 * 0.6,
    "Cramer_III": 90.0 * 0.5
}

TOXICOPHORES = {
    "Tropane_Alkaloid_Core": "[C]1[C][C]2[C][C]1[N]2", 
    "Potent_Anesthetic_Ester": "O=C(OCCN(CC)CC)c1ccc(N)cc1" 
}

def standardize_api(smiles: str) -> Tuple[float, float, Chem.Mol]:
    """Strips salts and returns canonical MW and LogP for the active API."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    
    # Strip counterions (e.g., HCl, Na+) to analyze only the active parent
    fragment_chooser = rdMolStandardize.LargestFragmentChooser()
    clean_mol = fragment_chooser.choose(mol)
    
    mw = Descriptors.MolWt(clean_mol)
    logp = Descriptors.MolLogP(clean_mol)
    return mw, logp, clean_mol

def evaluate_structural_alerts(mol: Chem.Mol) -> None:
    for alert_name, smarts in TOXICOPHORES.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            raise ToxicophoreMatchException(f"Critical Alert: Structural match for {alert_name}. Systemic risk too high for topical batching.")

def calculate_fickian_absorption(mw: float, logp: float, f_ui: float, conc_pct: float, product_type: str) -> Tuple[float, float]:
    """
    Calculates exact Dermal Absorption % using Fick's First Law and SCCS Surface Area.
    """
    if product_type not in SCCS_MECHANISTIC_REGISTRY:
        raise PhysicochemicalConstraintError(f"Product '{product_type}' not in mechanistic registry.")
        
    amount_mg, area_cm2, retention = SCCS_MECHANISTIC_REGISTRY[product_type]
    
    # Base Permeability (Potts & Guy)
    log_kp = (0.71 * logp) - (0.0061 * mw) - 2.72
    
    # Mucosal Stratum Corneum Bypass mechanism
    if "lip" in product_type.lower() or "oral" in product_type.lower():
        # Mucosa is 10-100x more permeable, highly dependent on MW, less on LogP
        log_kp = -0.0061 * mw - 1.5 
    
    kp_cm_hr = 10 ** log_kp
    
    # Volumetric concentration (assuming vehicle density ~ 1000 mg/cm3)
    c_vehicle_mg_cm3 = (conc_pct / 100.0) * 1000.0 
    
    # Fick's Law: J (Flux) = Kp * C_vehicle
    flux_mg_cm2_hr = kp_cm_hr * c_vehicle_mg_cm3 * f_ui
    
    # Total mass absorbed over 24 hours
    total_absorbed_mg = flux_mg_cm2_hr * area_cm2 * 24.0 * retention
    
    # Theoretical maximum applied mass of the API
    total_api_applied_mg = amount_mg * (conc_pct / 100.0) * retention
    
    # Dermal Absorption % (capped at 100%)
    da_pct = min((total_absorbed_mg / total_api_applied_mg) * 100.0, 100.0) if total_api_applied_mg > 0 else 0.0
    
    return da_pct, (total_absorbed_mg * 1000.0) # return DA% and SED in micrograms

def run_diep_gatekeeper(smiles: str, conc_pct: float, ph: float, pka: float, is_base: bool, product_type: str) -> Dict:
    mw, logp, clean_mol = standardize_api(smiles)
    evaluate_structural_alerts(clean_mol)
    
    if not (0 <= ph <= 14):
        raise PhysicochemicalConstraintError("Formulation pH must be 0-14.")
    f_ui = 1 / (1 + 10**(pka - ph)) if is_base else 1 / (1 + 10**(ph - pka))
    
    da_pct, sed_ug_day = calculate_fickian_absorption(mw, logp, f_ui, conc_pct, product_type)
    
    ttc_limit = EFSA_TTC_SYSTEMIC_UG_DAY["Cramer_III"]
    hepatic_burden = sed_ug_day / ttc_limit
    
    return {
        "api_mw": mw,
        "api_logp": logp,
        "f_ui": f_ui,
        "da_pct_applied": da_pct,
        "sed_ug_day": sed_ug_day,
        "ttc_limit_ug": ttc_limit,
        "hepatic_burden_ratio": hepatic_burden,
        "status": "PASS" if hepatic_burden <= 1.0 else "FAIL"
    }
