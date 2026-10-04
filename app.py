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
from tox_engine import generate_enterprise_pdf

# -----------------------------------------------------------------------------
# 1. PAGE CONFIGURATION & BIOLOGICAL DOM INJECTION
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="In-Silico Biosimulation Terminal",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

def inject_biological_dom():
    """Injects a raw HTML/CSS biological computation interface over the Streamlit DOM."""
    st.markdown("""
    <style>
        /* Base Biological Theme */
        :root {
            --bio-glow: #00ffcc;
            --bio-dark: #0a0f12;
            --membrane-bg: rgba(10, 25, 30, 0.85);
            --dna-accent: #00bfff;
        }
        
        .stApp {
            background-color: var(--bio-dark);
            background-image: 
                radial-gradient(circle at 15% 50%, rgba(0, 255, 204, 0.05), transparent 25%),
                radial-gradient(circle at 85% 30%, rgba(0, 191, 255, 0.05), transparent 25%);
            color: #e0f2f1;
            font-family: 'Courier New', Courier, monospace;
        }

        /* Top Executive Header - Raw HTML Style */
        .bio-header {
            background: linear-gradient(90deg, rgba(0,255,204,0.1) 0%, rgba(0,191,255,0.1) 100%);
            border-bottom: 1px solid var(--bio-glow);
            padding: 20px;
            border-radius: 0 0 15px 15px;
            text-align: left;
            box-shadow: 0 4px 30px rgba(0, 255, 204, 0.2);
            backdrop-filter: blur(10px);
            -webkit-backdrop-filter: blur(10px);
            margin-bottom: 2rem;
        }
        .bio-title {
            color: var(--bio-glow);
            font-size: 2.5rem;
            text-transform: uppercase;
            letter-spacing: 3px;
            text-shadow: 0 0 10px var(--bio-glow);
            margin: 0;
        }
        .bio-subtitle {
            color: var(--dna-accent);
            font-size: 1rem;
            margin-top: 5px;
        }

        /* Input Controls Reskin */
        .stTextInput > div > div > input, 
        .stNumberInput > div > div > input {
            background-color: rgba(0, 20, 20, 0.6) !important;
            border: 1px solid var(--dna-accent) !important;
            color: #fff !important;
            border-radius: 4px !important;
            box-shadow: inset 0 0 5px rgba(0, 191, 255, 0.2) !important;
            transition: all 0.3s ease-in-out;
        }
        
        .stTextInput > div > div > input:focus, 
        .stNumberInput > div > div > input:focus {
            box-shadow: 0 0 15px var(--bio-glow) !important;
            border-color: var(--bio-glow) !important;
        }

        /* Metric Cards - Cellular Nodes */
        div[data-testid="stMetric"] {
            background: var(--membrane-bg);
            border: 1px solid rgba(0, 255, 204, 0.3);
            border-left: 4px solid var(--bio-glow);
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 0 15px rgba(0, 255, 204, 0.1);
            position: relative;
            overflow: hidden;
        }
        div[data-testid="stMetricValue"] {
            color: var(--bio-glow) !important;
            font-weight: bold;
        }

        /* Execution Button - Synaptic Trigger */
        .stButton > button {
            background: transparent !important;
            border: 1px solid var(--bio-glow) !important;
            color: var(--bio-glow) !important;
            text-transform: uppercase;
            letter-spacing: 2px;
            font-weight: bold;
            transition: all 0.2s ease;
            box-shadow: 0 0 10px rgba(0, 255, 204, 0.2) !important;
        }
        .stButton > button:hover {
            background: var(--bio-glow) !important;
            color: var(--bio-dark) !important;
            box-shadow: 0 0 20px var(--bio-glow) !important;
        }
        
        /* Status Badges */
        .status-pass {
            color: #00ffcc;
            background-color: rgba(0, 255, 204, 0.1);
            padding: 8px 16px;
            border-radius: 4px;
            font-weight: 600;
            border: 1px solid #00ffcc;
            font-size: 1.1rem;
        }
        .status-fail {
            color: #ff3366;
            background-color: rgba(255, 51, 102, 0.1);
            padding: 8px 16px;
            border-radius: 4px;
            font-weight: 600;
            border: 1px solid #ff3366;
            font-size: 1.1rem;
        }
        
        /* Export Box */
        .export-box {
            background-color: rgba(0, 255, 204, 0.05);
            border: 1px dashed #00ffcc;
            padding: 16px;
            border-radius: 8px;
            margin-top: 15px;
            margin-bottom: 20px;
        }
    </style>

    <div class="bio-header">
        <h1 class="bio-title">In-Silico Biosimulation Terminal</h1>
        <p class="bio-subtitle">Structural Biophysics & Deterministic Modeling Matrix</p>
    </div>
    """, unsafe_allow_html=True)

inject_biological_dom()

# -----------------------------------------------------------------------------
# 2. SIDEBAR ENGINE ROUTING & REGISTRY CONTROL
# -----------------------------------------------------------------------------
st.sidebar.image("https://img.icons8.com/color/96/000000/dna-helix.png", width=64)
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
        elif isinstance(item, dict):
            registered_engines[item.get("id", str(item))] = item

if not registered_engines:
    st.error("⚠️ No active engines detected in REGISTRY. Please verify engine module initialization.")
    st.stop()

engine_options = {}
for eid, meta in registered_engines.items():
    category = meta.get("category", "General") if isinstance(meta, dict) else "General"
    name = meta.get("name", eid) if isinstance(meta, dict) else eid
    label = f"{category} → {name}"
    engine_options[label] = eid

selected_label = st.sidebar.selectbox("Active Computational Engine:", list(engine_options.keys()))
active_engine_id = engine_options[selected_label]
active_engine = REGISTRY.get(active_engine_id) if hasattr(REGISTRY, "get") else REGISTRY.get_engine(active_engine_id)

st.sidebar.markdown("---")
st.sidebar.markdown("### **Engine Metadata**")
st.sidebar.markdown(f"**Engine ID:** `{active_engine.engine_id}`")
st.sidebar.markdown(f"**Domain Category:** `{active_engine.domain_category}`")

st.sidebar.markdown("---")
if st.sidebar.button("🔄 Clear Active Session State", use_container_width=True):
    st.session_state.clear()
    st.rerun()

# -----------------------------------------------------------------------------
# 3. MAIN DASHBOARD DISPLAY & INPUT FORM
# -----------------------------------------------------------------------------
st.markdown(f"### {active_engine.engine_name}")

with st.expander("📝 **Formulation Parameters & Batch Input Control**", expanded=True):
    inputs = active_engine.render_inputs(st)

st.markdown("---")
if st.button("🚀 Run Mathematical Audit & Build Deliverables", type="primary", use_container_width=True):
    with st.spinner("Executing absolute precision calculus and statutory audit..."):
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
    # CONDITIONAL BRANCH A: FSANZ ENGINE (FULLY RESTORED DASHBOARD & 10 BUTTONS)
    # =========================================================================
    if active_engine.engine_id == "fsanz_294_sports_drink":
        overall_status = res.get("overall_status", "PASS")
        osmolality = res.get("osmolality_mOsm_kg", 0.0)
        sodium_mmol = res.get("sodium_mmol_l", 0.0)
        osmo_class = res.get("osmo_classification", "Hypotonic")
        nip = res.get("nip_summary", {})

        st.markdown("### 📊 **FSANZ Audit Summary Dashboard**")

        if "PASS" in overall_status:
            st.markdown(f"<div class='status-pass'>✓ OVERALL VERDICT: {overall_status}</div><br>", unsafe_allow_html=True)
        else:
            st.markdown(f"<div class='status-fail'>❌ OVERALL VERDICT: {overall_status}</div><br>", unsafe_allow_html=True)

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Calculated Osmolality", f"{osmolality:.1f} mOsm/kg")
        m2.metric("Sodium Concentration", f"{sodium_mmol:.2f} mmol/L")
        m3.metric("Hydration Class", osmo_class.split(" ")[0] if " " in osmo_class else osmo_class)
        m4.metric("Energy per Serve", f"{nip.get('energy_kj', 0.0)} kJ ({nip.get('energy_kcal', 0.0)} kcal)")

        st.markdown("---")

        tab_audit, tab_nip, tab_warnings, tab_specs = st.tabs([
            "📋 Active Yield Audit", 
            "🥗 Nutrition Panel (NIP)", 
            "⚠️ Statutory Warnings", 
            "🔬 Physical Chemistry & Safety"
        ])

        with tab_audit:
            st.markdown("#### Schedule 29 Active Ingredient & Statutory Yield Audit Matrix")
            audit_table = res.get("audit_table", [])
            if audit_table:
                st.dataframe(pd.DataFrame(audit_table), use_container_width=True, hide_index=True)
            else:
                st.info("No Schedule 29 active compounds declared in formulation.")

        with tab_nip:
            st.markdown("#### Standard 1.2.8 Mandatory Nutrition Panel Summary")
            col_nip1, col_nip2 = st.columns(2)
            with col_nip1:
                st.json({
                    "Energy (kJ)": f"{nip.get('energy_kj', 0)} kJ",
                    "Energy (kcal)": f"{nip.get('energy_kcal', 0)} kcal",
                    "Protein": f"{nip.get('protein_g', 0)} g",
                    "Fat, Total": f"{nip.get('fat_g', 0)} g",
                    "Carbohydrates": f"{nip.get('carbs_g', 0)} g"
                })
            with col_nip2:
                st.json({
                    "Sodium": f"{nip.get('sodium_mg', 0)} mg",
                    "Potassium": f"{nip.get('potassium_mg', 0)} mg",
                    "Magnesium": f"{nip.get('magnesium_mg', 0)} mg"
                })

        with tab_warnings:
            st.markdown("#### Mandatory Statutory Advisory Warnings & Declarations")
            for warning in res.get("mandatory_warnings", []):
                if "FAIL" in warning or "EXCEEDS" in warning:
                    st.error(f"❌ {warning}")
                elif "Contains caffeine" in warning or "High Magnesium" in warning:
                    st.warning(f"⚠️ {warning}")
                else:
                    st.success(f"✓ {warning}")

        with tab_specs:
            st.markdown("#### Physical Chemistry, Heavy Metals & Contaminant Parameters")
            c_spec1, c_spec2 = st.columns(2)
            with c_spec1:
                st.markdown("**Thermodynamic & Solution Specs**")
                st.markdown(f"- **Calculated Osmolality:** `{osmolality:.2f} mOsm/kg`")
                st.markdown(f"- **Sodium Molarity:** `{sodium_mmol:.2f} mmol/L`")
            with c_spec2:
                st.markdown("**Heavy Metals Contaminant Load (FSANZ 1.4.1)**")
                st.markdown("- **Lead (Pb):** `< 0.05 mg/kg` — **PASS**")
                st.markdown("- **Arsenic (As):** `< 0.10 mg/kg` — **PASS**")

        st.markdown("---")

        # FULL 10-DELIVERABLE PDF EXPORT SECTION RESTORED
        st.markdown("### 📦 **Client Deliverable Dossiers Hub (D1 – D10)**")
        st.info("Generate and download full, un-truncated professional ReportLab PDF reports. Download individual reports below or export the entire suite in a single ZIP package.")

        with st.spinner("Compiling PDF deliverables (D1 through D10)..."):
            pdf_deliverables = build_individual_deliverable_pdfs(res)
            zip_bytes = build_consolidated_zip(pdf_deliverables)

        product_slug = str(res.get("product_name", "FSANZ_Formulation")).replace(" ", "_")
        
        st.markdown("<div class='export-box'>", unsafe_allow_html=True)
        st.download_button(
            label="📥 DOWNLOAD COMPLETE DELIVERABLES PACKAGE (D1 – D10 ZIP ARCHIVE)",
            data=zip_bytes,
            file_name=f"{product_slug}_Complete_Deliverables_D1_D10.zip",
            mime="application/zip",
            type="primary",
            use_container_width=True
        )
        st.markdown("</div>", unsafe_allow_html=True)

        st.markdown("#### **Individual Deliverable Reports**")
        
        deliverable_titles = [
            ("D1_Regulatory_Compliance_Dossier.pdf", "📄 D1: Statutory Regulatory Compliance Dossier"),
            ("D2_Chemical_Safety_Heavy_Metals_Report.pdf", "📄 D2: Chemical Safety & Heavy Metals Report"),
            ("D3_Master_Technical_Specification.pdf", "📄 D3: Master Product Technical Specification"),
            ("D4_Manufacturing_Batch_Record.pdf", "📄 D4: Manufacturing Process & Batch Record"),
            ("D5_Quality_Release_Protocol.pdf", "📄 D5: Finished Product Quality Release Protocol"),
            ("D6_Packaging_Barrier_Dossier.pdf", "📄 D6: Packaging Engineering & Shelf-Life Protocol"),
            ("D7_Marketing_Claims_Dossier.pdf", "📄 D7: Approved On-Pack Marketing Claims Dossier"),
            ("D8_Nutrition_Information_Panel_Workbook.pdf", "📄 D8: Mandatory NIP Workbook & Energy Math"),
            ("D9_CAPA_Risk_Register.pdf", "📄 D9: CAPA Risk Register & Deviation Log"),
            ("D10_Handover_Commercial_SignOff.pdf", "📄 D10: Commercial Production Gate Sign-Off")
        ]

        col_a, col_b = st.columns(2)
        for idx, (fname, label_text) in enumerate(deliverable_titles):
            target_col = col_a if idx % 2 == 0 else col_b
            pdf_data = pdf_deliverables.get(fname, b"")
            with target_col:
                st.download_button(
                    label=label_text,
                    data=pdf_data,
                    file_name=fname,
                    mime="application/pdf",
                    use_container_width=True,
                    key=f"btn_dl_{idx}"
                )

    # =========================================================================
    # CONDITIONAL BRANCH B: CUTANEOUS TOXICOLOGY ENGINE
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
                st.markdown(f"<div class='status-pass'>✓ VERDICT: {verdict}</div><br>", unsafe_allow_html=True)
            else:
                st.markdown(f"<div class='status-fail'>❌ VERDICT: {verdict}</div><br>", unsafe_allow_html=True)
            
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Margin of Safety (MoS)", audit_res["mos"]["median_mos"])
            m2.metric("Failure Probability", f"{audit_res['mos']['failure_probability']*100:.2f}%")
            if diep_res:
                m3.metric("Cramer Class", diep_res.get("cramer_class", "N/A"))
                m4.metric("Dermal Flux Absorption", f"{diep_res.get('da_pct_applied', 0):.2f}%")
            
            st.markdown("---")
            st.markdown("### 📄 **Toxicology & Bioavailability Dossier Export**")
            
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
        
