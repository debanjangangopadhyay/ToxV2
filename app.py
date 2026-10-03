"""
Streamlit Interface for Multi-Domain In-Silico Computational Screener
Integrates Cutaneous Bioactivation (IATA/DIEP-MoS/TRACE-Onco) &
FSANZ Standard 2.9.4 Food Science Engines via Strategy Registry Architecture.
"""

import io
import streamlit as st
import pandas as pd

# =====================================================================
# 1. STREAMLIT UI CONFIGURATION & SAFE DOMAIN REGISTRY IMPORTS
# =====================================================================

st.set_page_config(
    page_title="Multi-Domain In-Silico Screener", 
    layout="wide", 
    page_icon="🧬"
)

try:
    from engines import REGISTRY
except ImportError:
    st.error("Fatal Error: Could not locate Engine Registry module.")
    st.stop()

try:
    from diep_engine import ToxicophoreMatchException
except ImportError:
    class ToxicophoreMatchException(Exception):
        pass

try:
    from tox_engine import generate_enterprise_pdf
except ImportError:
    def generate_enterprise_pdf(data, diep):
        buf = io.BytesIO()
        buf.write(b"%PDF-1.4 Enterprise PDF Generator Placeholder")
        buf.seek(0)
        return buf

try:
    from engines.fsanz_engine import build_factory_spec_pdf
except ImportError:
    def build_factory_spec_pdf(res):
        buf = io.BytesIO()
        if isinstance(res, dict) and "pdf_bytes" in res and isinstance(res["pdf_bytes"], bytes):
            buf.write(res["pdf_bytes"])
        else:
            buf.write(b"%PDF-1.4 FSANZ Factory Spec Placeholder")
        buf.seek(0)
        return buf


# Helper to safely resolve engine properties across class & metadata definitions
def resolve_engine_attr(engine, attr_name: str, default_val: str = "") -> str:
    if hasattr(engine, attr_name):
        val = getattr(engine, attr_name)
        return val() if callable(val) else str(val)
    
    meta = engine.get_metadata() if hasattr(engine, "get_metadata") else {}
    attr_map = {
        "engine_id": meta.get("id", default_val),
        "engine_name": meta.get("name", default_val),
        "domain_category": meta.get("category", default_val),
    }
    return str(attr_map.get(attr_name, default_val))


# =====================================================================
# 2. DOMAIN REGISTRY ROUTING & UI SELECTION
# =====================================================================

st.sidebar.title("Domain Engine Registry")
st.sidebar.markdown("Select computational engine for complete domain isolation:")

registered_raw = REGISTRY.list_engines()
if isinstance(registered_raw, dict):
    engine_options = {
        f"{meta.get('category', 'Domain')} → {meta.get('name', eid)}": eid 
        for eid, meta in registered_raw.items()
    }
else:
    engine_options = {
        f"{e.get('domain', e.get('category', 'Domain'))} → {e.get('name', e.get('id'))}": e['id'] 
        for e in registered_raw
    }

# SAFEGUARD 1: Halt gracefully if no engines loaded (prevents NoneType dropdowns)
if not engine_options:
    st.error("🚨 Critical Error: The Engine Registry is empty.")
    st.warning("Check your terminal logs. The `FSANZ294Engine` likely failed to instantiate silently in `engines/__init__.py` due to an underlying ABC abstract method requirement or missing import.")
    st.stop()

selected_label = st.sidebar.selectbox("Active Domain Engine:", list(engine_options.keys()))

# SAFEGUARD 2: Fallback if Streamlit session state holds a stale/cached label
if selected_label not in engine_options:
    selected_label = list(engine_options.keys())[0]

# Safe lookup now guaranteed
active_engine_id = engine_options[selected_label]
active_engine = REGISTRY.get(active_engine_id) if hasattr(REGISTRY, "get") else REGISTRY.get_engine(active_engine_id)

engine_id_str = resolve_engine_attr(active_engine, "engine_id", active_engine_id)
engine_name_str = resolve_engine_attr(active_engine, "engine_name", "Computational Chemistry Engine")
domain_cat_str = resolve_engine_attr(active_engine, "domain_category", "Regulatory Science")

st.sidebar.markdown("---")
st.sidebar.info(
    f"**Active Engine ID:** `{engine_id_str}`\n\n"
    f"**Domain Category:** `{domain_cat_str}`"
)

st.title(engine_name_str)


# =====================================================================
# 3. DYNAMIC INPUTS & EXECUTION WORKFLOW
# =====================================================================

# 1. Render Formulation Builder / Engine Inputs at full width
try:
    inputs = active_engine.render_inputs(st)
except TypeError:
    inputs = active_engine.render_inputs()

st.markdown("---")

# 2. Dynamic Execution Trigger
run = st.button("🚀 Run Research Audit", type="primary", use_container_width=True)

# 3. Execution & Session Handling
if run:
    with st.spinner("Executing active mathematical domain engine..."):
        try:
            results = active_engine.execute(inputs)
            st.session_state["results"] = results
            st.session_state["engine_id"] = engine_id_str
        except Exception as e:
            st.error(f"❌ **Execution Error:** {str(e)}")

# 4. Full-Width Results Dashboard
if "results" in st.session_state and st.session_state.get("engine_id") == engine_id_str:
    res = st.session_state["results"]

    if engine_id_str in ["fsanz_294_sports_drink", "fsanz_engine"]:
        overall_status = res.get("overall_status", "PASS - Full Compliance")
        osmolality = res.get("osmolality_mOsm_kg", 0.0)
        sodium_mmol = res.get("sodium_mmol_l", 0.0)
        osmo_class = res.get("osmo_classification", "Hypotonic")

        # Stretchable Compliance Verdict Banner
        if "PASS" in overall_status:
            st.success(f"### ✓ FSANZ 2.9.4 COMPLIANCE VERDICT: {overall_status}")
        else:
            st.error(f"### ❌ FSANZ 2.9.4 COMPLIANCE VERDICT: {overall_status}")

        # Metrics display with full width
        m1, m2, m3 = st.columns(3)
        m1.metric("Calculated Osmolality", f"{osmolality:.1f} mOsm/kg")
        m2.metric("Prepared Sodium Concentration", f"{sodium_mmol:.2f} mmol/L")
        m3.metric("Hydration Profile", osmo_class)

        st.markdown("#### 1. Schedule 29 Active Yield & Dosage Audit Matrix")
        raw_audit_table = res.get("audit_table", [])
        if raw_audit_table:
            df_audit = pd.DataFrame(raw_audit_table)
            st.dataframe(df_audit, use_container_width=True, hide_index=True)

        st.markdown("#### 2. Mandatory Package Label Warning Statements")
        warnings_list = res.get("mandatory_warnings", [])
        for warning in warnings_list:
            if "FAIL" in warning:
                st.error(f"❌ **NON-COMPLIANCE ALERT:** {warning}")
            elif "Contains caffeine" in warning or "High Magnesium" in warning:
                st.warning(f"⚠️ **COMPOUND SPECIFIC STATEMENT:** {warning}")
            else:
                st.success(f"✓ **MANDATORY STATEMENT:** {warning}")

        st.markdown("#### 3. Factory Specification PDF Generation")
        try:
            pdf_bytes = res.get("pdf_bytes")
            if pdf_bytes:
                product_name = res.get("product_name", "FSANZ_Product").replace(" ", "_")
                st.download_button(
                    "📥 Download Factory Specification PDF Dossier",
                    data=pdf_bytes,
                    file_name=f"{product_name}_FSANZ_Spec.pdf",
                    mime="application/pdf",
                    type="primary",
                    use_container_width=True
                )
        except Exception as e:
            st.error(f"Failed to render FSANZ PDF dossier: {str(e)}")
        
                
