import math
from typing import Dict, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

class ToxicophoreMatchException(Exception): pass
class PhysicochemicalConstraintError(ValueError): pass

SCCS_MECHANISTIC_REGISTRY = {
    "Face Cream (Leave-on)": (1540.0, 565.0, 1.0),
    "Body Lotion (Leave-on)": (7820.0, 15670.0, 1.0),
    "Facial Cleanser (Rinse-off)": (20000.0, 565.0, 0.01),
    "Shampoo (Rinse-off)": (10460.0, 1440.0, 0.01),
    "Eye Cream (Leave-on)": (500.0, 16.0, 1.0),
    "Lip Balm (Leave-on)": (57.0, 4.8, 1.0)
}

EFSA_TTC_SYSTEMIC_UG_DAY = {"Cramer_I": 1440.0, "Cramer_II": 324.0, "Cramer_III": 45.0}

TOXICOPHORES = {
    "Tropane_Alkaloid_Core": "[C]1[C][C]2[C][C]1[N]2", 
    "Potent_Anesthetic_Ester": "O=C(OCCN(CC)CC)c1ccc(N)cc1" 
}

def standardize_api(smiles: str) -> Tuple[float, float, Chem.Mol]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError("Invalid SMILES.")
    clean_mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
    return Descriptors.MolWt(clean_mol), Descriptors.MolLogP(clean_mol), clean_mol

def evaluate_cramer_class(mol: Chem.Mol) -> str:
    """Dynamically assigns Cramer Class based on topological complexity."""
    mw = Descriptors.MolWt(mol)
    rings = Descriptors.NumAromaticRings(mol)
    heteroatoms = Descriptors.NumHeteroatoms(mol)
    
    if mw < 250 and rings == 0 and heteroatoms <= 2:
        return "Cramer_I"
    elif mw > 400 or rings > 1 or heteroatoms > 4:
        return "Cramer_III"
    return "Cramer_II"

def calculate_dual_pathway_flux(mw: float, logp: float, ph: float, pka: float, is_base: bool, conc_pct: float, product_type: str) -> float:
    f_ui = 1 / (1 + 10**(pka - ph)) if is_base else 1 / (1 + 10**(ph - pka))
    f_i = 1.0 - f_ui
    
    kp_ui_cm_hr = 10 ** ((0.71 * logp) - (0.0061 * mw) - 2.72)
    kp_i_cm_hr = kp_ui_cm_hr * 0.01 

    if "lip" in product_type.lower() or "oral" in product_type.lower() or "eye" in product_type.lower():
        pore_multiplier = max(1.0, 150.0 * math.exp(-0.01 * mw)) 
        lipo_multiplier = max(1.0, 10.0 * math.exp(-0.005 * mw))
        kp_ui_cm_hr *= lipo_multiplier   
        kp_i_cm_hr *= pore_multiplier   
        
    k_vsc = 1.0 if logp < 3.0 else 0.5
    vehicle_density_g_cm3 = 1.0
    c_vehicle_mg_cm3 = (conc_pct / 100.0) * (vehicle_density_g_cm3 * 1000.0) * k_vsc
    
    return (kp_ui_cm_hr * c_vehicle_mg_cm3 * f_ui) + (kp_i_cm_hr * c_vehicle_mg_cm3 * f_i)

def run_diep_gatekeeper(smiles: str, conc_pct: float, ph: float, pka: float, is_base: bool, product_type: str) -> Dict:
    mw, logp, clean_mol = standardize_api(smiles)
    
    for alert, smarts in TOXICOPHORES.items():
        if clean_mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            raise ToxicophoreMatchException(f"Critical Systemic Hazard: {alert} detected.")
            
    if product_type not in SCCS_MECHANISTIC_REGISTRY:
        raise PhysicochemicalConstraintError(f"Product '{product_type}' not found in registry.")
        
    amount_mg, area_cm2, retention = SCCS_MECHANISTIC_REGISTRY[product_type]
    total_flux = calculate_dual_pathway_flux(mw, logp, ph, pka, is_base, conc_pct, product_type)
    
    total_api_applied_mg = amount_mg * (conc_pct / 100.0) * retention
    theoretical_absorbed_mg = total_flux * area_cm2 * 24.0 * retention
    absorbed_mass_mg = min(theoretical_absorbed_mg, total_api_applied_mg)
    
    da_pct = (absorbed_mass_mg / total_api_applied_mg * 100.0) if total_api_applied_mg > 0 else 0.0
    sed_ug_day = absorbed_mass_mg * 1000.0
    
    cramer_class = evaluate_cramer_class(clean_mol)
    ttc_limit = EFSA_TTC_SYSTEMIC_UG_DAY[cramer_class]
    hepatic_burden = sed_ug_day / ttc_limit if ttc_limit > 0 else float('inf')
    f_ui = 1 / (1 + 10**(pka - ph)) if is_base else 1 / (1 + 10**(ph - pka))
    
    return {
        "api_mw": mw,
        "api_logp": logp,
        "f_ui": f_ui,
        "da_pct_applied": da_pct,
        "cramer_class": cramer_class,
        "sed_ug_day": sed_ug_day,
        "ttc_limit_ug": ttc_limit,
        "hepatic_burden_ratio": hepatic_burden,
        "status": "PASS" if hepatic_burden <= 1.0 else "FAIL"
    }
    
