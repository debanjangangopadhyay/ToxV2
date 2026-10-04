"""
Streamlit UI Application for Multi-Domain In-Silico Computational Screener.
Direct integration with FSANZ294Engine for real-time calculation, NIP rendering,
deliverables execution (D1-D10), and client PDF document delivery.
"""

import io
import streamlit as st
import pandas as pd

# Page Configuration
st.set_page_config(
    page_title="Multi-Domain In-Silico Screener Engine", 
    layout="wide", 
    page_icon="🧬"
)

# Registry Imports
try:
    from engines import REGISTRY
except ImportError:
    st.error("Fatal Error: Engine Registry module could not be imported.")
    st.stop()

# Helper attribute resolver
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


# Sidebar Domain Routing
st.sidebar.title("Domain Engine Registry")
st.sidebar.markdown("Select computational engine for domain-isolated evaluation:")

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

if not engine_options:
    st.error("🚨 Critical Error: The Engine Registry is empty. Ensure FSANZ294Engine is registered.")
    st.stop()

selected_label = st.sidebar.selectbox("Active Domain Engine:", list(engine_options.keys()))
active_engine_id = engine_options[selected_label]
active_engine = REGISTRY.get(active_engine_id) if hasattr(REGISTRY, "get") else REGISTRY.get_engine(active_engine_id)

engine_id_str = resolve_engine_attr(active_engine, "engine_id", active_engine_id)
engine_name_str = resolve_engine_attr(active_engine, "engine_name", "FSANZ Regulatory Engine")
domain_cat_str = resolve_engine_attr(active_engine, "domain_category", "Food Science & Regulatory Chemistry")

st.sidebar.markdown("---")
st.sidebar.info(
    f"**Active Engine ID:** `{engine_id_str}`\n\n"
    f"**Domain Category:** `{domain_cat_str}`"
)

st.title(engine_name_str)

# Render Input Components
try:
    inputs = active_engine.render_inputs(st)
except TypeError:
    inputs = active_engine.render_inputs()

st.markdown("---")

# Execution Action Trigger
run = st.button("🚀 Run Mathematical Audit & Deliverables Pipeline", type="primary", use_container_width=True)

if run:
    with st.spinner("Executing mathematical physical chemistry & regulatory audit..."):
        try:
            results = active_engine.execute(inputs)
            st.session_state["results"] = results
            st.session_state["engine_id"] = engine_id_str
        except Exception as e:
            st.error(f"❌ Execution Error: {str(e)}")

# Dashboard Results Display
if "results" in st.session_state and st.session_state.get("engine_id") == engine_id_str:
    res = st.session_state["results"]

    if engine_id_str in ["fsanz_294_sports_drink", "fsanz_engine"]:
        overall_status = res.get("overall_status", "PASS - Full Compliance")
        osmolality = res.get("osmolality_mOsm_kg", 0.0)
        sodium_mmol = res.get("sodium_mmol_l", 0.0)
        osmo_class = res.get("osmo_classification", "Hypotonic")
        nip = res.get("nip_summary", {})

        if "PASS" in overall_status:
            st.success(f"### ✓ FSANZ AUDIT VERDICT: {overall_status}")
        else:
            st.error(f"### ❌ FSANZ AUDIT VERDICT: {overall_status}")

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Calculated Osmolality", f"{osmolality:.1f} mOsm/kg")
        m2.metric("Prepared Sodium Concentration", f"{sodium_mmol:.2f} mmol/L")
        m3.metric("Hydration Classification", osmo_class.split(" ")[0])
        m4.metric("Energy per Serve", f"{nip.get('energy_kj', 0.0)} kJ ({nip.get('energy_kcal', 0.0)} kcal)")

        st.markdown("#### 1. Schedule 29 Active Ingredient & Yield Audit Matrix")
        raw_audit_table = res.get("audit_table", [])
        if raw_audit_table:
            st.dataframe(pd.DataFrame(raw_audit_table), use_container_width=True, hide_index=True)

        st.markdown("#### 2. Standard 1.2.8 Nutrition Information Panel (NIP) Summary")
        c1, c2 = st.columns(2)
        with c1:
            st.json({
                "Energy (kJ / kcal)": f"{nip.get('energy_kj', 0)} kJ / {nip.get('energy_kcal', 0)} kcal",
                "Protein": f"{nip.get('protein_g', 0)} g",
                "Fat, Total": f"{nip.get('fat_g', 0)} g",
                "Carbohydrates": f"{nip.get('carbs_g', 0)} g",
                "Dietary Fibre": f"{nip.get('fibre_g', 0)} g"
            })
        with c2:
            st.json({
                "Polyols": f"{nip.get('polyols_g', 0)} g",
                "Organic Acids": f"{nip.get('org_acids_g', 0)} g",
                "Sodium": f"{nip.get('sodium_mg', 0)} mg",
                "Potassium": f"{nip.get('potassium_mg', 0)} mg",
                "Magnesium": f"{nip.get('magnesium_mg', 0)} mg"
            })

        st.markdown("#### 3. Mandatory Advisory & Warning Statements")
        warnings_list = res.get("mandatory_warnings", [])
        for warning in warnings_list:
            if "FAIL" in warning:
                st.error(f"❌ **NON-COMPLIANCE ALERT:** {warning}")
            elif "Contains caffeine" in warning or "High Magnesium" in warning:
                st.warning(f"⚠️️ **COMPOUND ADVISORY:** {warning}")
            else:
                st.success(f"✓ **MANDATORY STATEMENT:** {warning}")

        st.markdown("#### 4. Complete Client Deliverables Matrix (D1 - D10, Step 12, Step 13)")
        deliverables = res.get("deliverables_summary", {})
        if deliverables:
            st.dataframe(pd.DataFrame([deliverables]).T.rename(columns={0: "Deliverable Client Data Payload"}), use_container_width=True)

        st.markdown("#### 5. Instant Dossier Delivery to Client")
        try:
            pdf_bytes = res.get("pdf_bytes")
            if pdf_bytes:
                product_name = res.get("product_name", "FSANZ_Product").replace(" ", "_")
                st.download_button(
                    "📥 Direct Download Factory Specification PDF Dossier",
                    data=pdf_bytes,
                    file_name=f"{product_name}_FSANZ_Regulatory_Dossier.pdf",
                    mime="application/pdf",
                    type="primary",
                    use_container_width=True
                )
        except Exception as e:
            st.error(f"Failed to export client PDF dossier: {str(e)}")
            
