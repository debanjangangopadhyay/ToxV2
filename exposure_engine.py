import numpy as np
from typing import Dict, Any

def calculate_probabilistic_mos(
    pod_value: float, 
    pod_type: str = "NOAEL", 
    species: str = "Rat", 
    duration: str = "Subchronic",
    product_amount_g: float = 1.54, 
    concentration_pct: float = 2.0, 
    retention_factor: float = 1.0,
    dermal_abs_mean: float = 50.0, 
    n_iterations: int = 10000,
    random_seed: int = 42
) -> Dict[str, Any]:
    """
    Executes a Monte Carlo probabilistic exposure & Margin of Safety (MoS) simulation.
    All scalar outputs are strictly cast to Python native float types.
    """
    # Toxicological Assessment Factors (AF)
    af_inter = 1.0 if species.lower() == "human" else 10.0
    af_intra = 10.0 
    
    duration_lower = duration.lower()
    if duration_lower == "chronic": af_duration = 1.0
    elif duration_lower == "subchronic": af_duration = 3.0
    else: af_duration = 10.0
        
    pod_lower = pod_type.lower()
    if pod_lower == "loael": af_pod = 3.0
    elif "bmdl" in pod_lower: af_pod = 1.0
    else: af_pod = 1.0
        
    target_af = float(af_inter * af_intra * af_duration * af_pod)

    # Vectorized Monte Carlo Distribution Sampling
    np.random.seed(random_seed)
    
    bw_dist = np.random.normal(60.0, 10.2, n_iterations).clip(40.0, 120.0)
    da_std = dermal_abs_mean * 0.30
    da_dist = np.random.normal(dermal_abs_mean, da_std, n_iterations).clip(0.1, 100.0) / 100.0
    amount_dist = np.random.normal(product_amount_g, product_amount_g * 0.15, n_iterations).clip(product_amount_g * 0.5, product_amount_g * 2.0)
    pod_dist = np.random.normal(pod_value, pod_value * 0.10, n_iterations).clip(pod_value * 0.5, None)

    # Vectorized SED & MoS Evaluation
    conc_fraction = concentration_pct / 100.0
    applied_mg_dist = amount_dist * 1000.0
    
    sed_dist = (applied_mg_dist * conc_fraction * retention_factor * da_dist) / bw_dist
    mos_dist = pod_dist / sed_dist

    # Statistical Extraction & Risk Verdict
    median_mos = float(np.median(mos_dist))
    ci_05 = float(np.percentile(mos_dist, 5))
    ci_95 = float(np.percentile(mos_dist, 95))
    
    failure_probability = float(np.sum(mos_dist < target_af) / n_iterations)
    
    if failure_probability < 0.05:
        verdict = "High Confidence of Safety"
        is_safe = True
    elif failure_probability < 0.15:
        verdict = "Marginal Safety: Formulation Refinement Suggested"
        is_safe = False
    else:
        verdict = "High Exposure Risk: Threshold Failure"
        is_safe = False

    return {
        "valid": True,
        "target_assessment_factor": target_af,
        "af_breakdown": {
            "interspecies": float(af_inter),
            "intraspecies": float(af_intra),
            "duration": float(af_duration),
            "data_quality": float(af_pod)
        },
        "simulation_results": {
            "iterations": int(n_iterations),
            "median_sed_mg_kg_day": float(round(float(np.median(sed_dist)), 6)),
            "median_mos": float(round(median_mos, 1)),
            "mos_5th_percentile": float(round(ci_05, 1)),
            "mos_95th_percentile": float(round(ci_95, 1)),
            "failure_probability_pct": float(round(failure_probability * 100.0, 2))
        },
        "verdict": verdict,
        "is_safe": is_safe
    }
    
