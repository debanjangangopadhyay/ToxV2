import streamlit as st
import pandas as pd
from tox_engine import execute_full_compound_audit, generate_enterprise_pdf

st.set_page_config(page_title="In-Silico Safety Architecture", layout="wide", page_icon="🧬")

st.title("In-Silico Cutaneous Bioactivation & Toxicological Screener")
st.markdown("Integrates Monte Carlo MoS vectorization, bio-structural homology, and Dempster-Shafer evidentiary fusion.")

col1, col2 = st.columns([1, 3])
with col1:
    st.subheader("Formulation Inputs")
    smiles = st.text_input("Target SMILES", "CC(=O)Oc1ccccc1C(=O)O")
    conc = st.number_input("Concentration (%)", 0.0, 100.0, 2.0)
    pod = st.number_input("NOAEL (mg/kg/day)", 0.1, 5000.0, 250.0)
    da = st.slider("Dermal Absorption (%)", 0.1, 100.0, 50.0)
    run = st.button("Run Enterprise Audit", type="primary")

with col2:
    if run and smiles:
        with st.spinner("Executing mathematical integrations..."):
            res = execute_full_compound_audit(smiles, pod, conc, da)
            
        if not res["valid"]: st.error(res["error"])
        else:
            st.session_state["audit"] = res

if "audit" in st.session_state:
    data = st.session_state["audit"]
    
    m1, m2, m3 = st.columns(3)
    m1.metric("Integrated Belief (Safe)", f"{data['dst_metrics']['belief_safe']}%")
    m2.metric("Epistemic Uncertainty", f"{data['dst_metrics']['epistemic_uncertainty']}%")
    m3.metric("5th Percentile MoS", data["mos"]["ci_05_mos"])
    
    t1, t2, t3, t4, t5 = st.tabs(["Dempster-Shafer Fusion", "Monte Carlo MoS", "Bio-Homology", "Regulatory Provenance", "Export"])
    
    with t1:
        st.subheader("Evidentiary Mass Distribution")
        st.json(data["dst_metrics"])
        
    with t2:
        st.subheader("Probabilistic Exposure Simulation")
        st.write(f"**Failure Probability:** {data['mos']['failure_probability']*100}% against dynamic Assessment Factor target of {data['mos']['target_af']}.")
        
    with t3:
        st.subheader("Multidimensional Read-Across")
        if data["analogs"]: st.dataframe(pd.DataFrame(data["analogs"]), use_container_width=True)
        
    with t4:
        st.subheader("Decoupled Legal Provenance")
        if data["reg"]["violations"]: st.dataframe(pd.DataFrame(data["reg"]["violations"]), use_container_width=True)
        else: st.success("No prohibited regulatory structures identified.")
        
    with t5:
        st.download_button("📥 Download Official IATA Dossier", data=generate_enterprise_pdf(data).getvalue(), file_name=f"Toxicology_{data['assessment_id']}.pdf", mime="application/pdf")
