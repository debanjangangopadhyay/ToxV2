"""
Streamlit Interface for IATA Computational Screener.
Provides reactive execution parameters, product scenario selection, ablation testing controls, and PDF retrieval.
"""

import streamlit as st
import pandas as pd
from tox_engine_4 import (
    execute_full_compound_audit, 
    generate_enterprise_pdf, 
    SCCS_PRODUCT_EXPOSURE
)

st.set_page_config(page_title="IATA Computational Screener", layout="wide", page_icon="🧬")

st.title("In-Silico Cutaneous Bioactivation & Toxicological Screener")
st.caption("Research-Use Computational Assessment | Not for Regulatory Submission")
st.markdown("Integrates Monte Carlo MoS vectorization, bio-structural homology, and empirically-calibrated Dempster-Shafer evidentiary fusion.")

col1, col2 = st.columns([1, 3])

with col1:
    st.subheader("Formulation Inputs")
    
    # 1. Chemical Structure Input
    smiles = st.text_input("Target SMILES", "CC(=O)Oc1ccccc1C(=O)O")
    
    # 2. Exposure Scenario Selection (Fixes invalid product_type error)
    product_type = st.selectbox(
        "Product Exposure Scenario",
        options=list(SCCS_PRODUCT_EXPOSURE.keys()),
        index=0,
        help="Select SCCS product category to load standardized daily application rate and retention factor."
    )
    
    # 3. Quantitative Formulation Parameters
    conc = st.number_input("Concentration (%)", min_value=0.01, max_value=100.0, value=2.0, step=0.1)
    pod = st.number_input("NOAEL / POD (mg/kg/day)", min_value=0.1, max_value=10000.0, value=250.0, step=10.0)
    da = st.slider("Mean Dermal Absorption (%)", min_value=0.1, max_value=100.0, value=50.0)
    bw = st.number_input("Body Weight (kg)", min_value=10.0, max_value=150.0, value=60.0)
    
    with st.expander("Advanced Validation Mode (Ablation)"):
        use_ds = st.checkbox("Enable Dempster-Shafer Fusion", value=True)
        use_met = st.checkbox("Enable Cutaneous Metabolism", value=True)
        use_ana = st.checkbox("Enable Read-Across", value=True)
        mc_samples = st.number_input("Monte Carlo Iterations", min_value=1000, max_value=100000, value=10000, step=1000)
        
    run = st.button("Run Research Audit", type="primary")

with col2:
    if run and smiles:
        with st.spinner("Executing scenario-aware Monte Carlo & evidence fusion..."):
            # Call engine using EXPLICIT KEYWORD ARGUMENTS to eliminate positional mismatches
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
            st.error(f"**Audit Execution Failed:** {res.get('error', 'Unknown validation error')}")
        else: 
            st.session_state["audit"] = res

if "audit" in st.session_state:
    data = st.session_state["audit"]
    
    st.info(f"**Assessment ID:** `{data['assessment_id']}` | **Canonical SMILES:** `{data['canonical_smiles']}`")
    
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Integrated Belief (Safe)", f"{data['dst_metrics']['belief_safe']}%")
    m2.metric("Epistemic Uncertainty", f"{data['dst_metrics']['epistemic_uncertainty']}%")
    m3.metric("5th Percentile MoS", data["mos"]["ci_05_mos"])
    m4.metric("Failure Probability", f"{data['mos']['failure_probability']*100:.2f}%")
    
    t1, t2, t3, t4, t5 = st.tabs([
        "Evidence Calibration & DS Fusion", 
        "Monte Carlo MoS", 
        "Bio-Homology & Read-Across", 
        "Regulatory Provenance", 
        "Export Assessment"
    ])
    
    with t1:
        st.subheader("Evidentiary Mass Distribution")
        
        col_ds1, col_ds2 = st.columns(2)
        with col_ds1:
            st.json(data["dst_metrics"])
        with col_ds2:
            st.markdown("""
            ### Evidence Fusion Interpretation
            * **Belief (Safe):** Discounted lower bound of certainty supporting chemical safety.
            * **Belief (Toxic):** Lower bound of certainty indicating toxic hazard / exposure failure.
            * **Epistemic Uncertainty:** Unassigned probability mass due to incomplete data or model conflict.
            * **Plausibility (Safe):** Maximum possible belief in safety ($1.0 - \\text{Belief(Toxic)}$).
            """)
        
    with t2:
        st.subheader("Probabilistic Exposure Simulation")
        st.write(f"**Product Scenario:** `{data['contract'].product_type}`")
        st.write(f"**Median Systemic Exposure Dose (SED):** `{data['mos']['median_sed_mg_kg_day']} mg/kg/day`")
        st.write(f"**Median Margin of Safety (MoS):** `{data['mos']['median_mos']}`")
        st.write(f"**Failure Probability:** `{data['mos']['failure_probability']*100:.2f}%` against target Assessment Factor of `{data['mos']['target_af']}`.")
        st.write(f"**Dynamic Variance (CV_DA):** `{data['mos']['dynamic_cv_da']}` (Scaled to molecular LogP).")
        
    with t3:
        st.subheader("Multidimensional Read-Across")
        if data.get("analogs"): 
            st.dataframe(pd.DataFrame(data["analogs"]), use_container_width=True)
        else:
            st.warning("No read-across analogs calculated.")

        if data.get("met"):
            st.subheader("Cutaneous Bioactivation Hypotheses")
            st.dataframe(pd.DataFrame(data["met"]), use_container_width=True)
        
    with t4:
        st.subheader("Decoupled Legal Provenance")
        reg_status = data.get("reg", {})
        if reg_status.get("violations"): 
            st.error("Prohibited Structural Features Detected:")
            st.dataframe(pd.DataFrame(reg_status["violations"]), use_container_width=True)
        else: 
            st.success("✓ No prohibited regulatory structures identified in active schema.")
            
        if reg_status.get("warnings"):
            st.warning("Regulatory Warnings / Restrictions:")
            st.dataframe(pd.DataFrame(reg_status["warnings"]), use_container_width=True)

        if data.get("alerts"):
            st.subheader("Structural Alerts (PAINS / BRENK)")
            for alert in data["alerts"]:
                st.write(f"- ⚠️ `{alert}`")
        
    with t5:
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
            
