import numpy as np
from typing import Dict, List, Any, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, rdFingerprintGenerator
from rdkit import DataStructs

# Extended reference database. In a production environment, this is loaded from a JSON/SQL schema.
REFERENCE_DATABASE = [
    {"name": "Bakuchiol", "inci": "Bakuchiol", "category": "Antioxidant", "smiles": "CC(=CCC/C(=C/CC1=CC=C(C=C1)O)/C)C"},
    {"name": "Niacinamide", "inci": "Niacinamide", "category": "Skin Barrier", "smiles": "NC(=O)c1cccnc1"},
    {"name": "Squalane", "inci": "Squalane", "category": "Emollient", "smiles": "CCCCCCCCCCCC(C)CCCC(C)CCCC(C)CCCC(C)C"},
    {"name": "Allantoin", "inci": "Allantoin", "category": "Keratolytic", "smiles": "O=C1NC(=O)NC1NC(=O)N"},
    {"name": "Phenoxyethanol", "inci": "Phenoxyethanol", "category": "Preservative", "smiles": "OCCc1ccccc1"},
    {"name": "Bisabolol", "inci": "Bisabolol", "category": "Anti-inflammatory", "smiles": "CC(=CCCC(C)(C1CCC(=CC1)C)O)C"}
]

def _calculate_physicochemical_profile(mol: Chem.Mol) -> Dict[str, float]:
    """Extracts critical parameters governing dermal absorption and systemic transport."""
    return {
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol)
    }

def execute_multidimensional_read_across(
    target_mol: Optional[Chem.Mol], 
    top_k: int = 3,
    weight_topology: float = 0.50,
    weight_logp: float = 0.25,
    weight_tpsa: float = 0.15,
    weight_mw: float = 0.10
) -> List[Dict[str, Any]]:
    """
    Computes a Bio-Structural Homology Score integrating 2D Morgan fingerprints 
    with 3D physicochemical continuous variables to satisfy OECD RAAF requirements.
    """
    if target_mol is None:
        return []

    try:
        # Generate 2D Topology Fingerprint
        morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        target_fp = morgan_gen.GetFingerprint(target_mol)
        
        # Calculate 3D/Physicochemical Profile
        target_physio = _calculate_physicochemical_profile(target_mol)
    except Exception:
        return []

    candidates = []

    for ref in REFERENCE_DATABASE:
        ref_mol = Chem.MolFromSmiles(ref["smiles"])
        if ref_mol is None:
            continue

        try:
            # 1. Topological Tanimoto Similarity
            ref_fp = morgan_gen.GetFingerprint(ref_mol)
            tanimoto_sim = DataStructs.TanimotoSimilarity(target_fp, ref_fp)
            
            # Convert similarity to a distance penalty (0 = identical topology)
            dist_topology = 1.0 - tanimoto_sim

            # 2. Physicochemical Delta Calculations
            ref_physio = _calculate_physicochemical_profile(ref_mol)
            
            # Normalization Denominators (Max expected biological variance in cosmetics)
            norm_logp = 5.0   # Lipophilicity variance constraint
            norm_tpsa = 100.0 # Polarity variance constraint
            norm_mw = 300.0   # Mass variance constraint

            dist_logp = min(abs(target_physio["logp"] - ref_physio["logp"]) / norm_logp, 1.0)
            dist_tpsa = min(abs(target_physio["tpsa"] - ref_physio["tpsa"]) / norm_tpsa, 1.0)
            dist_mw = min(abs(target_physio["mw"] - ref_physio["mw"]) / norm_mw, 1.0)

            # 3. Weighted Composite Distance
            composite_distance = (
                (dist_topology * weight_topology) +
                (dist_logp * weight_logp) +
                (dist_tpsa * weight_tpsa) +
                (dist_mw * weight_mw)
            )

            # Convert final distance back to a 0-100% Homology Score
            homology_score = max((1.0 - composite_distance) * 100.0, 0.0)

            candidates.append({
                "name": ref["name"],
                "inci": ref["inci"],
                "category": ref["category"],
                "raw_tanimoto_pct": round(tanimoto_sim * 100, 1),
                "bio_homology_score_pct": round(homology_score, 1),
                "physio_deltas": {
                    "delta_logp": round(abs(target_physio["logp"] - ref_physio["logp"]), 2),
                    "delta_tpsa": round(abs(target_physio["tpsa"] - ref_physio["tpsa"]), 1),
                    "delta_mw": round(abs(target_physio["mw"] - ref_physio["mw"]), 1)
                },
                "interpretation": "Mechanistically plausible analog" if homology_score >= 75.0 else "Weak mechanistic homology"
            })
            
        except Exception:
            continue

    # Rank by the integrated Bio-Structural Homology Score, not raw Tanimoto
    candidates.sort(key=lambda x: x["bio_homology_score_pct"], reverse=True)
    return candidates[:max(1, int(top_k))]
  
