"""
Streamlit Interface for IATA Computational Screener (v5.3)
Integrates TRACE-Onco compatible Deterministic Biophysics (DIEP-MoS) & Dynamic AF Matrix.
"""
import streamlit as st
import pandas as pd

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
    product_type = st.selectbox("Product Scenario", list(SCCS_MECHANISTIC_REGISTRY.keys()), index=0)
    conc = st.number_input("Concentration (%)", min_value=0.01, max_value=100.0, value=2.0)
    
    st.markdown("#### Toxicological Point of Departure (PoD)")
    pod = st.number_input("PoD Value (mg/kg/day)", min_value=0.1, max_value=10000.0, value=250.0)
    
    pod_col1, pod_col2, pod_col3 = st.columns(3)
    with pod_col1:
        pod_type = st.selectbox("PoD Type", ["NOAEL", "LOAEL", "BMDL"], index=0)
    with pod_col2:
        species = st.selectbox("Test Species", ["Rat", "Mouse", "Dog", "Rabbit", "Human"], index=0)
    with pod_col3:
        duration = st.selectbox("Study Duration", ["Chronic", "Subchronic", "Subacute"], index=1)
    
    st.markdown("---")
    st.markdown("### Deterministic Biophysics (DIEP-MoS)")
    use_deterministic_flux = st.checkbox("Enable Fickian Deterministic Flux", value=True)
    formulation_ph = st.number_input("Formulation pH", min_value=0.0, max_value=14.0, value=5.5)
    pka = st.number_input("Compound pKa", min_value=0.0, max_value=14.0, value=8.6)
    is_base = st.checkbox("Molecule is a Weak Base?", value=True)
    
    da = st.slider("Mean Dermal Absorption (%)", 0.1, 100.0, 50.0, disabled=use_deterministic_flux)
    bw = st.number_input("Body Weight (kg)", 10.0, 150.0, 60.0)
    
    # Cast to integer and set safe step intervals for the Monte Carlo sampler
    mc_samples = int(st.number_input("Monte Carlo Iterations", min_value=1000, max_value=100000, value=10000, step=1000))
    
    run = st.button("Run Research Audit", type="primary")

with col2:
    if run and smiles:
        with st.spinner("Executing biophysical gates and Monte Carlo simulation..."):
            
            diep_results = None
            if use_deterministic_flux:
                try:
                    diep_results = run_diep_gatekeeper(smiles, conc, formulation_ph, pka, is_base, product_type)
                    da = diep_results["da_pct_applied"] 
                except ToxicophoreMatchException as e:
                    st.error(f"❌ **CRITICAL FATAL ALERT:** {str(e)}")
                    st.stop()
                except Exception as e:
                    st.error(f"❌ **Biophysical Calculation Error:** {str(e)}")
                    st.stop()

            res = execute_full_compound_audit(
                smiles=smiles,
                product_type=product_type,
                concentration_pct=conc,
                pod_noael_mg_kg_day=pod,
                pod_type=pod_type,
                species=species,
                duration=duration,
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
    
    st.info(f"**Assessment ID:** `{data['assessment_id']}` | **Canonical SMILES:** `{data['canonical_smiles']}`")
    
    t1, t2, t3 = st.tabs(["Monte Carlo & Biophysics", "TRACE-Onco Synergy Payload", "Export Audit Ledger"])
    
    with t1:
        st.subheader("Systemic Bioavailability & Simulation Limits")
        
        if diep:
            st.markdown("#### 1. Fickian Deterministic Bounds (DIEP-MoS)")
            st.write(f"**Canonical API Properties:** MW: `{diep['api_mw']:.2f}` | LogP: `{diep['api_logp']:.2f}`")
            st.write(f"**Topological Classification:** `{diep['cramer_class']}`")
            st.write(f"**Calculated Unionized Fraction ($f_{{ui}}$):** `{diep['f_ui']:.4f}`")
            st.write(f"**Fickian Dermal Absorption (DA%):** `{diep['da_pct_applied']:.2f}%`")
            st.write(f"**Max Systemic Exposure Dose (SED):** `{diep['sed_ug_day']:.2f} µg/day` vs TTC Limit `{diep['ttc_limit_ug']:.2f} µg/day`")
            
            if diep['status'] == "FAIL":
                st.error("❌ **DETERMINISTIC FAILURE:** Absolute systemic exposure exceeds safe EFSA thresholds.")
            else:
                st.success("✓ **DETERMINISTIC PASS:** Systemic exposure is within safe limits.")
        
        st.markdown("#### 2. Probabilistic Exposure (Monte Carlo)")
        st.write(f"**Median MoS:** `{data['mos']['median_mos']}` | **Failure Prob:** `{data['mos']['failure_probability']*100:.2f}%` against AF target of `{data['mos']['target_af']}`")
        st.json(data["mos"]["af_breakdown"])

    with t2:
        st.subheader("TRACE-Onco / VMTB Output Vector")
        if diep:
            st.json({
                "patient_hepatic_burden_ratio": diep["hepatic_burden_ratio"],
                "bioavailability_status": diep["status"],
                "structural_alerts": [alert for alert in data.get("alerts", [])],
                "recommendation": "Integrate ratio directly into decentralized Lifelines Cox-PH model, dynamically weighted against patient baseline De Ritis ratio to account for hepatic stress."
            })
        else:
            st.info("Enable Deterministic Biophysics to generate the downstream clinical integration payload.")
            
    with t3:
        st.subheader("Generate & Download PDF Ledger")
        try:
            pdf_bytes = generate_enterprise_pdf(data).getvalue()
            st.download_button(
                "📥 Download IATA-Aligned Computational Assessment (PDF)", 
                data=pdf_bytes, 
                file_name=f"{data['assessment_id']}.pdf", 
                mime="application/pdf",
                type="primary"
            )
        except Exception as e:
            st.error(f"PDF Generation Error: {str(e)}")
            
