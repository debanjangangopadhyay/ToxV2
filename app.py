"""
Streamlit Interface for IATA Computational Screener (v3.0)
Integrates TRACE-Onco compatible Deterministic Biophysics (DIEP-MoS).
"""
import streamlit as st
import pandas as pd

# External modules (ensure these are in the same directory)
from tox_engine import execute_full_compound_audit, generate_enterprise_pdf
from diep_engine import (
    run_diep_gatekeeper,
    SCCS_MECHANISTIC_REGISTRY,
    ToxicophoreMatchException,
    PhysicochemicalConstraintError
)

st.set_page_config(page_title="IATA Computational Screener", layout="wide", page_icon="🧬")
st.title("In-Silico Cutaneous Bioactivation & Toxicological Screener")

col1, col2 = st.columns([1, 3])

with col1:
    st.subheader("Formulation Inputs")
    smiles = st.text_input("Target SMILES", "CC(=O)Oc1ccccc1C(=O)O")
    
    # Maps dynamically to the DIEP registry, eliminating the previous hardcoded bug
    product_type = st.selectbox("Product Scenario", list(SCCS_MECHANISTIC_REGISTRY.keys()), index=0)
    conc = st.number_input("Concentration (%)", min_value=0.01, max_value=100.0, value=2.0)
    pod = st.number_input("NOAEL / POD (mg/kg/day)", min_value=0.1, max_value=10000.0, value=250.0)
    
    st.markdown("---")
    st.markdown("### Deterministic Biophysics (DIEP-MoS)")
    use_deterministic_flux = st.checkbox("Enable Fickian Deterministic Flux", value=True)
    formulation_ph = st.number_input("Formulation pH", min_value=0.0, max_value=14.0, value=5.5)
    pka = st.number_input("Compound pKa", min_value=0.0, max_value=14.0, value=8.6)
    is_base = st.checkbox("Molecule is a Weak Base?", value=True)
    
    da = st.slider("Mean Dermal Absorption (%)", 0.1, 100.0, 50.0, disabled=use_deterministic_flux)
    bw = st.number_input("Body Weight (kg)", 10.0, 150.0, 60.0)
    mc_samples = st.number_input("Monte Carlo Iterations", 1000, 100000, 10000)
    
    run = st.button("Run Research Audit", type="primary")

with col2:
    if run and smiles:
        with st.spinner("Executing biophysical gates and Monte Carlo simulation..."):
            
            diep_results = None
            if use_deterministic_flux:
                try:
                    diep_results = run_diep_gatekeeper(smiles, conc, formulation_ph, pka, is_base, product_type)
                    da = diep_results["da_pct_applied"] # Overwrite slider with Fickian math
                except ToxicophoreMatchException as e:
                    st.error(f"❌ **CRITICAL FATAL ALERT:** {str(e)}")
                    st.stop()
                except Exception as e:
                    st.error(f"❌ **Biophysical Calculation Error:** {str(e)}")
                    st.stop()

            # Execute probabilistic Monte Carlo engine with the newly bounded DA parameter
            res = execute_full_compound_audit(
                smiles=smiles,
                product_type=product_type,
                concentration_pct=conc,
                pod_noael_mg_kg_day=pod,
                dermal_absorption_pct=da, 
                body_weight_kg=bw,
                mc_samples=mc_samples
            )
            
            if not res.get("valid", False): 
                st.error(f"Audit Execution Failed: {res.get('error')}")
            else: 
                st.session_state["audit"] = res
                st.session_state["diep"] = diep_results

if "audit" in st.session_state:
    data = st.session_state["audit"]
    diep = st.session_state.get("diep")
    
    t1, t2 = st.tabs(["Monte Carlo & Biophysics", "TRACE-Onco Synergy Payload"])
    
    with t1:
        st.subheader("Systemic Bioavailability & Simulation Limits")
        
        if diep:
            st.markdown("#### 1. Fickian Deterministic Bounds (DIEP-MoS)")
            st.write(f"**Canonical API Properties:** MW: `{diep['api_mw']:.2f}` | LogP: `{diep['api_logp']:.2f}` (Salts Stripped)")
            st.write(f"**Calculated Unionized Fraction ($f_{{ui}}$):** `{diep['f_ui']:.4f}`")
            st.write(f"**Fickian Dermal Absorption (DA%):** `{diep['da_pct_applied']:.2f}%`")
            st.write(f"**Max Systemic Exposure Dose (SED):** `{diep['sed_ug_day']:.2f} µg/day` vs TTC Limit `{diep['ttc_limit_ug']:.2f} µg/day`")
            
            if diep['status'] == "FAIL":
                st.error("❌ **DETERMINISTIC FAILURE:** Absolute systemic exposure exceeds safe EFSA thresholds.")
            else:
                st.success("✓ **DETERMINISTIC PASS:** Systemic exposure is within safe limits.")
        
        st.markdown("#### 2. Probabilistic Exposure (Monte Carlo)")
        st.write(f"**Median MoS:** `{data['mos']['median_mos']}` | **Failure Prob:** `{data['mos']['failure_probability']*100:.2f}%`")

    with t2:
        st.subheader("TRACE-Onco / VMTB Output Vector")
        if diep:
            st.json({
                "patient_hepatic_burden_ratio": diep["hepatic_burden_ratio"],
                "bioavailability_status": diep["status"],
                "recommendation": "Integrate ratio into Lifelines Cox-PH model to adjust overall survival curves based on systemic tolerability."
            })
        else:
            st.info("Enable Deterministic Biophysics to generate the TRACE-Onco payload.")
            
