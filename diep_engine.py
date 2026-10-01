import math
from typing import Dict, Tuple, List, Any
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

class ToxicophoreMatchException(Exception): pass
class PhysicochemicalConstraintError(ValueError): pass

# SCCS Mechanistic Product Exposure Registry (Daily Amount mg, Surface Area cm², Retention Factor)
SCCS_MECHANISTIC_REGISTRY = {
    "Face Cream (Leave-on)": (1540.0, 565.0, 1.0),
    "Body Lotion (Leave-on)": (7820.0, 15670.0, 1.0),
    "Facial Cleanser (Rinse-off)": (20000.0, 565.0, 0.01),
    "Shampoo (Rinse-off)": (10460.0, 1440.0, 0.01),
    "Eye Cream (Leave-on)": (500.0, 16.0, 1.0),
    "Lip Balm (Leave-on)": (57.0, 4.8, 1.0)
}

# EFSA Threshold of Toxicological Concern (TTC) Systemic Human Limits (µg/person/day)
EFSA_TTC_SYSTEMIC_UG_DAY = {
    "Cramer_I": 1440.0, 
    "Cramer_II": 324.0, 
    "Cramer_III": 45.0
}

# Critical Mechanistic Toxicophores (SMARTS)
TOXICOPHORES = {
    "Tropane_Alkaloid_Core": "[C]1[C][C]2[C][C]1[N]2", 
    "Potent_Anesthetic_Ester": "O=C(OCCN(CC)CC)c1ccc(N)cc1",
    "Alkylating_Sulfonate": "COS(=O)(=O)C",
    "Aromatic_Nitrosamine": "c1ccccc1N(C)N=O"
}

# Decision Tree Structural SMARTS Definitions
CRAMER_TREE_SMARTS = {
    "heavy_metals": "[#12,#13,#20,#26,#27,#28,#29,#30,#33,#48,#50,#80,#82]",
    "quaternary_nitrogen": "[N+]",
    "organophosphorus": "P(=O)(O)(O)",
    "reactive_halo_alkyl": "[CX4][Cl,Br,I]",
    "aromatic_ring": "a",
    "heterocyclic_ring": "[!#6;R]",
    "ester_or_lactone": "[CX3](=O)[OX2H0][#6]",
    "aliphatic_chain_long": "[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]",
    "reactive_electrophile": "[#6]=O;![#6]O"
}

class ComputationalCramerDecisionTree:
    """
    Formal 8-Stage Computational Cramer Decision Tree based on explicit topological queries,
    SMARTS pattern matching, heteroatom ratios, and molecular complexity thresholds.
    """
    def __init__(self, mol: Chem.Mol):
        self.mol = mol
        self.log: List[str] = []

    def evaluate(self, has_critical_alerts: bool) -> Tuple[str, List[str]]:
        self.log.clear()
        
        # Stage 0: Structural Alert Override
        if has_critical_alerts:
            self.log.append("Stage 0: Critical toxicophore alert detected -> Escalated to Cramer Class III")
            return "Cramer_III", self.log

        # Stage 1: Organometallic & Quaternary Nitrogen Screening
        if self._matches("heavy_metals") or self._matches("quaternary_nitrogen") or self._matches("organophosphorus"):
            self.log.append("Stage 1: Organometallic center, organophosphate, or Quaternary Nitrogen present -> Cramer Class III")
            return "Cramer_III", self.log
        self.log.append("Stage 1: Passed organometallic and quaternary nitrogen screening.")

        # Compute Topological Descriptors
        mw = float(Descriptors.MolWt(self.mol))
        logp = float(Descriptors.MolLogP(self.mol))
        n_rings = Descriptors.RingCount(self.mol)
        n_hetero = Descriptors.NumHeteroatoms(self.mol)
        rot_bonds = Descriptors.NumRotatableBonds(self.mol)
        is_aromatic = self._matches("aromatic_ring")
        is_hetero_ring = self._matches("heterocyclic_ring")

        # Stage 2: Highly Reactive Electrophiles & Alkylating Agents
        if self._matches("reactive_halo_alkyl") and self._matches("reactive_electrophile"):
            self.log.append("Stage 2: Co-occurrence of electrophilic carbonyl and alkyl halide -> Cramer Class III")
            return "Cramer_III", self.log
        self.log.append("Stage 2: Passed reactive electrophile screening.")

        # Stage 3: High Polycyclic Complexity or Extreme Substituted Heterocycles
        if is_aromatic and (n_rings >= 3 or n_hetero > 5 or mw > 400.0):
            self.log.append(f"Stage 3: Polycyclic aromatic system (Rings: {n_rings}, Heteroatoms: {n_hetero}, MW: {mw:.1f}) -> Cramer Class III")
            return "Cramer_III", self.log
        self.log.append("Stage 3: Passed polycyclic complexity screening.")

        # Stage 4: Simple Acyclic / Carbocyclic Aliphatics (Class I Pathway)
        if not is_aromatic and not is_hetero_ring and mw <= 250.0 and n_hetero <= 3:
            self.log.append(f"Stage 4: Simple non-aromatic carbocyclic or acyclic structure (MW: {mw:.1f}, Heteroatoms: {n_hetero}) -> Cramer Class I")
            return "Cramer_I", self.log
        self.log.append("Stage 4: Structure exceeds Class I simple aliphatic boundaries.")

        # Stage 5: Common Functionalized Monocyclic Systems & Esters (Class II Pathway)
        if (self._matches("ester_or_lactone") or n_rings == 1) and mw <= 300.0 and n_hetero <= 4 and rot_bonds <= 6:
            self.log.append(f"Stage 5: Monocyclic system with standard functional groups (MW: {mw:.1f}, RotBonds: {rot_bonds}) -> Cramer Class II")
            return "Cramer_II", self.log
        self.log.append("Stage 5: Compound structural features exceed Class II monocyclic criteria.")

        # Stage 6: High Flexibility & Heteroatom Density Escalation
        if rot_bonds > 8 or n_hetero > 4 or logp > 4.5:
            self.log.append(f"Stage 6: High rotatable bond count ({rot_bonds}) or heteroatom density ({n_hetero}) -> Cramer Class III")
            return "Cramer_III", self.log

        # Stage 7: Fallback Default Classification
        if mw > 300.0:
            self.log.append(f"Stage 7: High molecular mass threshold exceeded (MW: {mw:.1f} > 300.0) -> Cramer Class III")
            return "Cramer_III", self.log

        self.log.append("Stage 7: Intermediate structural complexity -> Cramer Class II")
        return "Cramer_II", self.log

    def _matches(self, key: str) -> bool:
        smarts = CRAMER_TREE_SMARTS[key]
        query = Chem.MolFromSmarts(smarts)
        return self.mol.HasSubstructMatch(query) if query else False


def standardize_api(smiles: str) -> Tuple[float, float, Chem.Mol]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: 
        raise ValueError("Invalid SMILES structure provided.")
    clean_mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
    return float(Descriptors.MolWt(clean_mol)), float(Descriptors.MolLogP(clean_mol)), clean_mol


def calculate_multiprotic_fui(ph: float, acid_pkas: List[float], base_pkas: List[float]) -> float:
    """Calculates total unionized fraction (f_ui) across multiple ionization centers."""
    f_ui_total = 1.0
    for pka in acid_pkas:
        f_ui_total *= 1.0 / (1.0 + 10.0 ** (ph - pka))
    for pka in base_pkas:
        f_ui_total *= 1.0 / (1.0 + 10.0 ** (pka - ph))
    return max(float(f_ui_total), 1e-9)


def calculate_dual_pathway_flux(
    mw: float, 
    logp: float, 
    ph: float, 
    acid_pkas: List[float], 
    base_pkas: List[float], 
    conc_pct: float, 
    product_type: str
) -> float:
    """Calculates biophysical transdermal flux considering unionized & ionized species pathways."""
    f_ui = calculate_multiprotic_fui(ph, acid_pkas, base_pkas)
    f_i = max(0.0, 1.0 - f_ui)
    
    # Potts & Guy permeability coefficient model (cm/hr)
    kp_ui_cm_hr = 10.0 ** ((0.71 * logp) - (0.0061 * mw) - 2.72)
    kp_i_cm_hr = kp_ui_cm_hr * 0.01

    if any(k in product_type.lower() for k in ["lip", "oral", "eye"]):
        pore_multiplier = max(1.0, 150.0 * math.exp(-0.01 * mw)) 
        lipo_multiplier = max(1.0, 10.0 * math.exp(-0.005 * mw))
        kp_ui_cm_hr *= lipo_multiplier   
        kp_i_cm_hr *= pore_multiplier   
        
    k_vsc = 1.0 if logp < 3.0 else 0.5
    c_vehicle_mg_cm3 = (max(0.001, conc_pct) / 100.0) * 1000.0 * k_vsc
    
    return float((kp_ui_cm_hr * c_vehicle_mg_cm3 * f_ui) + (kp_i_cm_hr * c_vehicle_mg_cm3 * f_i))


def run_diep_gatekeeper(
    smiles: str, 
    conc_pct: float, 
    ph: float, 
    acid_pkas: List[float], 
    base_pkas: List[float], 
    product_type: str
) -> Dict[str, Any]:
    """Executes biophysical skin absorption & TTC hazard gating."""
    mw, logp, clean_mol = standardize_api(smiles)
    
    has_alert = False
    alert_name = None
    for alert, smarts in TOXICOPHORES.items():
        if clean_mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            has_alert = True
            alert_name = alert
            break
            
    if product_type not in SCCS_MECHANISTIC_REGISTRY:
        raise PhysicochemicalConstraintError(f"Product '{product_type}' missing from SCCS registry.")
        
    amount_mg, area_cm2, retention = SCCS_MECHANISTIC_REGISTRY[product_type]
    total_flux = calculate_dual_pathway_flux(mw, logp, ph, acid_pkas, base_pkas, conc_pct, product_type)
    
    total_api_applied_mg = amount_mg * (conc_pct / 100.0) * retention
    theoretical_absorbed_mg = total_flux * area_cm2 * 24.0 * retention
    absorbed_mass_mg = min(theoretical_absorbed_mg, total_api_applied_mg)
    
    da_pct = (absorbed_mass_mg / total_api_applied_mg * 100.0) if total_api_applied_mg > 0 else 0.0
    sed_ug_day = absorbed_mass_mg * 1000.0
    
    # Computational Cramer Decision Tree Execution
    tree_evaluator = ComputationalCramerDecisionTree(clean_mol)
    cramer_class, tree_log = tree_evaluator.evaluate(has_alert)
    
    ttc_limit = EFSA_TTC_SYSTEMIC_UG_DAY[cramer_class]
    hepatic_burden = sed_ug_day / ttc_limit if ttc_limit > 0 else float('inf')
    f_ui = calculate_multiprotic_fui(ph, acid_pkas, base_pkas)
    
    # Mathematical Oncogenic Hazard Index
    oncogenic_risk_index = (hepatic_burden * 0.40) + ((1.0 - f_ui) * 0.20) + (1.5 if has_alert else 0.0)

    if has_alert and alert_name in ["Tropane_Alkaloid_Core", "Alkylating_Sulfonate"]:
        raise ToxicophoreMatchException(f"Critical Systemic Hazard: High reactivity toxicophore '{alert_name}' detected.")

    return {
        "api_mw": mw,
        "api_logp": logp,
        "f_ui": f_ui,
        "da_pct_applied": float(da_pct),
        "cramer_class": cramer_class,
        "cramer_tree_log": tree_log,
        "sed_ug_day": float(sed_ug_day),
        "ttc_limit_ug": float(ttc_limit),
        "hepatic_burden_ratio": float(hepatic_burden),
        "oncogenic_risk_index": float(round(oncogenic_risk_index, 4)),
        "status": "PASS" if hepatic_burden <= 1.0 else "FAIL"
    }
    
