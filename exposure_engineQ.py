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
    n_iterations: int = 10000
) -> Dict[str, Any]:
    """
    Executes a Monte Carlo simulation for Systemic Exposure Dose (SED) and Margin of Safety (MoS),
    incorporating toxicological Assessment Factors (AF) for dynamic thresholding.
    """
    
    # ---------------------------------------------------------
    # 1. DYNAMIC ASSESSMENT FACTOR (AF) MATRIX
    # ---------------------------------------------------------
    # Interspecies extrapolation (Animal to Human)
    af_inter = 1.0 if species.lower() == "human" else 10.0
    
    # Intraspecies variance (Human to Human variability)
    af_intra = 10.0 
    
    # Exposure Duration Extrapolation to Chronic
    duration_lower = duration.lower()
    if duration_lower == "chronic":
        af_duration = 1.0
    elif duration_lower == "subchronic":
        af_duration = 3.0
    else:  # Subacute or acute
        af_duration = 10.0
        
    # Point of Departure (PoD) Quality Penalty
    pod_lower = pod_type.lower()
    if pod_lower == "loael":
        af_pod = 3.0  # Penalty for using Lowest Observed Adverse Effect Level
    elif "bmdl" in pod_lower:
        af_pod = 1.0  # Benchmark Dose Lower Confidence Limit is highly reliable
    else:
        af_pod = 1.0  # Standard NOAEL
        
    # The new target threshold (Replaces the hardcoded 100)
    target_af = af_inter * af_intra * af_duration * af_pod

    # ---------------------------------------------------------
    # 2. MONTE CARLO PROBABILITY DISTRIBUTIONS
    # ---------------------------------------------------------
    np.random.seed(42)  # Set seed for reproducibility in legal/audit contexts
    
    # Body Weight (Normal Dist: Mean 60kg, StdDev 10.2kg, clipped to biological limits)
    bw_dist = np.random.normal(60.0, 10.2, n_iterations).clip(40.0, 120.0)
    
    # Dermal Absorption (Normal Dist: 30% Coefficient of Variation)
    da_std = dermal_abs_mean * 0.30
    da_dist = np.random.normal(dermal_abs_mean, da_std, n_iterations).clip(0.1, 100.0) / 100.0
    
    # Consumer Usage Amount (Normal Dist: 15% Variation around SCCS default)
    amount_dist = np.random.normal(product_amount_g, product_amount_g * 0.15, n_iterations).clip(product_amount_g * 0.5, product_amount_g * 2.0)
    
    # Point of Departure Variance (Simulating experimental bench error, CV = 10%)
    pod_dist = np.random.normal(pod_value, pod_value * 0.10, n_iterations).clip(pod_value * 0.5, None)

    # ---------------------------------------------------------
    # 3. VECTORIZED EXPOSURE CALCULATION
    # ---------------------------------------------------------
    conc_fraction = concentration_pct / 100.0
    applied_mg_dist = amount_dist * 1000.0
    
    # SED Array: (Amount * Concentration * Retention * Absorption) / Body Weight
    sed_dist = (applied_mg_dist * conc_fraction * retention_factor * da_dist) / bw_dist
    
    # MoS Array: PoD / SED
    mos_dist = pod_dist / sed_dist

    # ---------------------------------------------------------
    # 4. STATISTICAL EXTRACTION & VERDICT LOGIC
    # ---------------------------------------------------------
    median_mos = np.median(mos_dist)
    ci_05 = np.percentile(mos_dist, 5)  # Conservative 5th percentile lower bound
    ci_95 = np.percentile(mos_dist, 95)
    
    # Probability of Failure: What percentage of the 10,000 simulated consumers fall below the safe threshold?
    failure_probability = np.sum(mos_dist < target_af) / n_iterations
    
    if failure_probability < 0.05:
        verdict = "High Confidence of Safety"
        is_safe = True
    elif failure_probability < 0.15:
        verdict = "Marginal Safety: Formulation Refinement Suggested"
        is_safe = False
    else:
        verdict = "High Regulatory Liability: Severe Exposure Risk"
        is_safe = False

    return {
        "valid": True,
        "target_assessment_factor": target_af,
        "af_breakdown": {
            "interspecies": af_inter,
            "intraspecies": af_intra,
            "duration": af_duration,
            "data_quality": af_pod
        },
        "simulation_results": {
            "iterations": n_iterations,
            "median_sed_mg_kg_day": round(np.median(sed_dist), 4),
            "median_mos": round(median_mos, 1),
            "mos_5th_percentile": round(ci_05, 1),
            "mos_95th_percentile": round(ci_95, 1),
            "failure_probability_pct": round(failure_probability * 100, 2)
        },
        "verdict": verdict,
        "is_safe": is_safe
    }
