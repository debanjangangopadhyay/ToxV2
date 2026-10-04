"""
===============================================================================
Production-Grade Streamlit Application (app.py)
Multi-Engine Regulatory, Physical Chemistry & Toxicology Framework
===============================================================================
"""

import streamlit as st
import pandas as pd
import io
import zipfile
from typing import Dict, Any, List

# Engine Registry and Factory Imports
from engines import REGISTRY
from engines.deliverables_factory import build_individual_deliverable_pdfs, build_consolidated_zip

# IMPORT THE MISSING PDF GENERATOR FOR THE TOX ENGINE
from tox_engine import generate_enterprise_pdf

# -----------------------------------------------------------------------------
# 1. PAGE CONFIGURATION & CUSTOM CSS STYLING
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="Computational Science Engine",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# (आपका बायोलॉजिकल DOM/CSS कोड यहाँ मौजूद रहेगा जो पिछले स्टेप में दिया गया था)
st.markdown("""
<style>
    .status-pass { color: #1E8449; background-color: #EAFAF1; padding: 6px 12px; border-radius: 4px; font-weight: 600; border: 1px solid #A9DFBF; }
    .status-fail { color: #922B21; background-color: #FDEDEC; padding: 6px 12px; border-radius: 4px; font-weight: 600; border: 1px solid #F9E79F; }
    .export-box { background-color: rgba(0, 255, 204, 0.1); border: 1px solid #00ffcc; padding: 16px; border-radius: 8px; margin-top: 15px; margin-bottom: 15px; }
</style>
""", unsafe_allow_html=True)

# -----------------------------------------------------------------------------
# 2. SIDEBAR ENGINE ROUTING & REGISTRY CONTROL
# -----------------------------------------------------------------------------
st.sidebar.markdown("## **Engine Control Hub**")

raw_engines = REGISTRY.list_engines() if hasattr(REGISTRY, "list_engines") else REGISTRY.get_all()
registered_engines = {}
if isinstance(raw_engines, dict):
    registered_engines = raw_engines
elif isinstance(raw_engines, list):
    for item in raw_engines:
        if hasattr(item, "get_metadata"):
            meta = item.get_metadata()
            registered_engines[meta.get("id", str(item))] = meta

engine_options = {}
for eid, meta in registered_engines.items():
    category = meta.get("category", "General")
    name = meta.get("name", eid)
    label = f"{category} → {name}"
    engine_options[label] = eid

selected_label = st.sidebar.selectbox("Active Computational Engine:", list(engine_options.keys()))
active_engine_id = engine_options[selected_label]
active_engine = REGISTRY.get(active_engine_id) if hasattr(REGISTRY, "get") else REGISTRY.get_engine(active_engine_id)

# -----------------------------------------------------------------------------
# 3. MAIN DASHBOARD DISPLAY & INPUT FORM
# -----------------------------------------------------------------------------
st.markdown(f"<h1>{active_engine.engine_name}</h1>", unsafe_allow_html=True)
st.markdown("### Formulation Parameters & Batch Input Control")

inputs = active_engine.render_inputs(st)

st.markdown("---")
if st.button("🚀 Run Mathematical Audit & Build Deliverables", type="primary", use_container_width=True):
    with st.spinner("Executing structural chemistry and statutory audit..."):
        try:
            results = active_engine.execute(inputs)
            st.session_state["results"] = results
            st.session_state["engine_id"] = active_engine.engine_id
            st.success("✅ Computational audit completed successfully.")
        except Exception as e:
            st.error(f"❌ Calculation Execution Error: {str(e)}")
            st.exception(e)

# -----------------------------------------------------------------------------
# 4. CONDITIONAL RESULTS DASHBOARD & DELIVERABLES EXPORT
# -----------------------------------------------------------------------------
if "results" in st.session_state and st.session_state.get("engine_id") == active_engine.engine_id:
    res = st.session_state["results"]
    
    st.markdown("---")
    
    # =========================================================================
    # CONDITIONAL BRANCH A: FSANZ ENGINE
    # =========================================================================
    if active_engine.engine_id == "fsanz_294_sports_drink":
        overall_status = res.get("overall_status", "PASS")
        
        st.markdown("### 📊 **FSANZ Audit Summary Dashboard**")
        if "PASS" in overall_status:
            st.markdown(f"<div class='status-pass'>✓ OVERALL VERDICT: {overall_status}</div>", unsafe_allow_html=True)
        else:
            st.markdown(f"<div class='status-fail'>❌ OVERALL VERDICT: {overall_status}</div>", unsafe_allow_html=True)

        m1, m2, m3 = st.columns(3)
        m1.metric("Calculated Osmolality", f"{res.get('osmolality_mOsm_kg', 0):.1f} mOsm/kg")
        m2.metric("Sodium Concentration", f"{res.get('sodium_mmol_l', 0):.2f} mmol/L")
        m3.metric("Hydration Class", res.get("osmo_classification", "N/A"))

        st.markdown("### 📦 **FSANZ Deliverable Dossiers Hub (D1 – D10)**")
        with st.spinner("Compiling FSANZ PDF deliverables..."):
            pdf_deliverables = build_individual_deliverable_pdfs(res)
            zip_bytes = build_consolidated_zip(pdf_deliverables)
            
        st.download_button(
            label="📥 DOWNLOAD FSANZ DELIVERABLES PACKAGE (ZIP)",
            data=zip_bytes,
            file_name="FSANZ_Complete_Deliverables.zip",
            mime="application/zip",
            type="primary",
            use_container_width=True
        )

    # =========================================================================
    # CONDITIONAL BRANCH B: CUTANEOUS TOXICOLOGY ENGINE (THE FIX)
    # =========================================================================
    elif active_engine.engine_id == "cutaneous_iata_diep":
        audit_res = res.get("audit", {})
        diep_res = res.get("diep")
        
        st.markdown("### 🧬 **Cutaneous Biophysics & Toxicology Dashboard**")
        
        if not audit_res.get("valid"):
            st.error(f"Audit Failed: {audit_res.get('error', 'Unknown Error')}")
        else:
            is_safe = audit_res["mos"]["failure_probability"] < 0.05
            verdict = "High Confidence of Safety" if is_safe else "Marginal Safety / High Exposure Risk"
            
            if is_safe:
                st.markdown(f"<div class='status-pass'>✓ VERDICT: {verdict}</div>", unsafe_allow_html=True)
            else:
                st.markdown(f"<div class='status-fail'>❌ VERDICT: {verdict}</div>", unsafe_allow_html=True)
            
            # Metrics Display
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Margin of Safety (MoS)", audit_res["mos"]["median_mos"])
            m2.metric("Failure Probability", f"{audit_res['mos']['failure_probability']*100:.2f}%")
            if diep_res:
                m3.metric("Cramer Class", diep_res.get("cramer_class", "N/A"))
                m4.metric("Dermal Flux Absorption", f"{diep_res.get('da_pct_applied', 0):.2f}%")
            
            st.markdown("---")
            st.markdown("### 📄 **Toxicology & Bioavailability Dossier Export**")
            st.info("The missing PDF generator feature has been successfully wired. Download your comprehensive ReportLab dossier below.")
            
            # THE MISSING PDF GENERATOR INTEGRATION
            with st.spinner("Rendering 2D Molecules and Generating PDF Report..."):
                try:
                    pdf_buffer = generate_enterprise_pdf(audit_data=audit_res, diep_data=diep_res)
                    
                    st.markdown("<div class='export-box'>", unsafe_allow_html=True)
                    st.download_button(
                        label="📥 DOWNLOAD COMPUTATIONAL TOXICOLOGY DOSSIER (PDF)",
                        data=pdf_buffer,
                        file_name=f"Toxicology_Dossier_{audit_res.get('assessment_id', 'Report')}.pdf",
                        mime="application/pdf",
                        type="primary",
                        use_container_width=True
                    )
                    st.markdown("</div>", unsafe_allow_html=True)
                except Exception as pdf_err:
                    st.error(f"Failed to generate PDF Document: {pdf_err}")
                    
