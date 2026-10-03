"""
Streamlit Interface for Multi-Domain In-Silico Computational Screener
Integrates Cutaneous Bioactivation (IATA/DIEP-MoS/TRACE-Onco) &
FSANZ Standard 2.9.4 Food Science Engines via Strategy Registry Architecture.
"""

import streamlit as st
import pandas as pd

# Strategy Registry Architecture & Domain Engine Imports
from engines import REGISTRY
from diep_engine import ToxicophoreMatchException
from tox_engine import generate_enterprise_pdf
from engines.fsanz_engine import build_factory_spec_pdf

# =====================================================================
# 1. STREAMLIT UI CONFIGURATION & DOMAIN REGISTRY ROUTING
# =====================================================================

st.set_page_config(
    page_title="Multi-Domain In-Silico Screener", 
    layout="wide", 
    page_icon="🧬"
)

st.sidebar.title("Domain Engine Registry")
st.sidebar.markdown("Select computational engine for complete domain isolation:")

registered_engines = REGISTRY.list_engines()
engine_options = {f"{e['domain']} → {e['name']}": e['id'] for e in registered_engines}

selected_label = st.sidebar.selectbox("Active Domain Engine:", list(engine_options.keys()))
active_engine_id = engine_options[selected_label]
active_engine = REGISTRY.get(active_engine_id)

st.sidebar.markdown("---")
st.sidebar.info(
    f"**Active Engine ID:** `{active_engine.engine_id}`\n\n"
    f"**Domain Category:** `{active_engine.domain_category}`"
)

st.title(active_engine.engine_name)

# =====================================================================
# 2. DYNAMIC INPUTS & EXECUTION WORKFLOW
# =====================================================================

col1, col2 = st.columns([1, 2.5])

with col1:
    inputs = active_engine.render_inputs(st)
    run = st.button("Run Research Audit", type="primary", use_container_width=True)

with col2:
    if run:
        with st.spinner("Executing active mathematical domain engine..."):
            try:
                results = active_engine.execute(inputs)
                st.session_state["results"] = results
                st.session_state["engine_id"] = active_engine.engine_id
            except ToxicophoreMatchException as e:
                st.error(f"❌ **CRITICAL FATAL ALERT:** {str(e)}")
            except Exception as e:
                st.error(f"❌ **Execution Error:** {str(e)}")

    if "results" in st.session_state and st.session_state.get("engine_id") == active_engine.engine_id:
        res = st.session_state["results"]

        # -------------------------------------------------------------
        # DOMAIN 1: CUTANEOUS BIOACTIVATION & TOXICOLOGICAL SCREENER
        # -------------------------------------------------------------
        if active_engine.engine_id == "cutaneous_iata_diep":
            data = res["audit"]
            diep = res["diep"]

            st.info(f"**Assessment ID:** `{data['assessment_id']}` | **Canonical SMILES:** `{data['canonical_smiles']}`")

            t1, t2, t3 = st.tabs(["Monte Carlo & Decision Tree", "TRACE-Onco Synergy Payload", "Export Audit Ledger"])

            with t1:
                st.subheader("Systemic Bioavailability & Decision Tree Log")

                if diep:
                    st.markdown("#### 1. Deterministic Bounds & Cramer Decision Tree (DIEP-MoS)")
                    st.write(f"**Canonical API Properties:** MW: `{diep['api_mw']:.2f}` g/mol | LogP: `{diep['api_logp']:.2f}`")
                    st.write(f"**Topological Classification:** `{diep['cramer_class']}`")
                    st.write(f"**Calculated Unionized Fraction ($f_{{ui}}$):** `{diep['f_ui']:.4f}`")
                    st.write(f"**Fickian Dermal Absorption (DA%):** `{diep['da_pct_applied']:.2f}%`")
                    st.write(f"**Max Systemic Exposure Dose (SED):** `{diep['sed_ug_day']:.2f} µg/day` vs TTC Limit `{diep['ttc_limit_ug']:.2f} µg/day`")

                    st.markdown("**Computational Cramer Tree Execution Trail:**")
                    for step in diep["cramer_tree_log"]:
                        st.text(f"  └── {step}")

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
                        "oncogenic_risk_index": diep["oncogenic_risk_index"],
                        "bioavailability_status": diep["status"],
                        "structural_alerts": [alert for alert in data.get("alerts", [])],
                        "recommendation": "Integrate ratio directly into decentralized Lifelines Cox-PH model, dynamically weighted against patient baseline De Ritis ratio to account for hepatic stress."
                    })
                else:
                    st.info("Enable Deterministic Biophysics to generate the downstream clinical integration payload.")

            with t3:
                st.subheader("Generate & Download PDF Ledger")
                try:
                    pdf_bytes = generate_enterprise_pdf(data, diep).getvalue()
                    st.download_button(
                        "📥 Download Multi-Section Computational Assessment Dossier (PDF)",
                        data=pdf_bytes,
                        file_name=f"{data['assessment_id']}.pdf",
                        mime="application/pdf",
                        type="primary"
                    )
                except Exception as e:
                    st.error(f"PDF Generation Error: {str(e)}")

        # -------------------------------------------------------------
        # DOMAIN 2: FSANZ STANDARD 2.9.4 SPORTS DRINK ENGINE
        # -------------------------------------------------------------
        elif active_engine.engine_id == "fsanz_294_sports_drink":
            if "PASS" in res["overall_status"]:
                st.success(f"✓ **FSANZ 2.9.4 COMPLIANCE VERDICT:** {res['overall_status']}")
            elif "WARNING" in res["overall_status"]:
                st.warning(f"⚠️ **FSANZ 2.9.4 COMPLIANCE VERDICT:** {res['overall_status']}")
            else:
                st.error(f"❌ **FSANZ 2.9.4 COMPLIANCE VERDICT:** {res['overall_status']}")

            m1, m2, m3 = st.columns(3)
            m1.metric("Calculated Osmolality", f"{res['osmolality_mOsm_kg']} mOsm/kg")
            m2.metric("Prepared Sodium Concentration", f"{res['sodium_mmol_l']:.2f} mmol/L")
            m3.metric("Hydration Profile", res["osmo_classification"])

            st.markdown("#### 1. Schedule 29 Active Yield & Dosage Audit Matrix")
            df_audit = pd.DataFrame(res["audit_table"])
            st.dataframe(df_audit, use_container_width=True)

            st.markdown("#### 2. Mandatory Package Label Warning Statements")
            for warning in res["mandatory_warnings"]:
                st.warning(f"⚠️ **REQUIRED STATEMENT:** {warning}")

            st.markdown("#### 3. Factory Specification PDF Generation")
            try:
                pdf_buffer = build_factory_spec_pdf(res)
                st.download_button(
                    "📥 Download Factory Specification PDF Dossier",
                    data=pdf_buffer.getvalue(),
                    file_name=f"{res['product_name'].replace(' ', '_')}_FSANZ_Spec.pdf",
                    mime="application/pdf",
                    type="primary"
                )
            except Exception as e:
                st.error(f"Failed to generate FSANZ PDF specification: {str(e)}")
                
