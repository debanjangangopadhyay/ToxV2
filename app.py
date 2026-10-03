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

# 1. Render Formulation Builder / Engine Inputs at full top width
try:
    inputs = active_engine.render_inputs(st)
except TypeError:
    inputs = active_engine.render_inputs()

st.markdown("---")

# 2. Prominent execution trigger spanning the container width
run = st.button("🚀 Run Research Audit", type="primary", use_container_width=True)

# 3. Audit Execution Logic
if run:
    with st.spinner("Executing active mathematical domain engine..."):
        try:
            results = active_engine.execute(inputs)
            st.session_state["results"] = results
            st.session_state["engine_id"] = engine_id_str
        except ToxicophoreMatchException as e:
            st.error(f"❌ **CRITICAL FATAL ALERT:** {str(e)}")
        except Exception as e:
            st.error(f"❌ **Execution Error:** {str(e)}")

# 4. Render Results Dashboard below the builder (Full Page Width)
if "results" in st.session_state and st.session_state.get("engine_id") == engine_id_str:
    res = st.session_state["results"]
    st.markdown("## 📊 Assessment Results & Compliance Ledger")

    # -------------------------------------------------------------
    # DOMAIN 1: CUTANEOUS BIOACTIVATION & TOXICOLOGICAL SCREENER
    # -------------------------------------------------------------
    if engine_id_str in ["cutaneous_iata_diep", "cutaneous_biophysics"]:
        data = res.get("audit", res)
        diep = res.get("diep", res if "api_mw" in res else None)

        assessment_id = data.get('assessment_id', 'ASSESSMENT-LOCAL')
        canonical_smiles = data.get('canonical_smiles', res.get('smiles', 'N/A'))

        st.info(f"**Assessment ID:** `{assessment_id}` | **Canonical SMILES:** `{canonical_smiles}`")

        t1, t2, t3 = st.tabs(["Monte Carlo & Decision Tree", "TRACE-Onco Synergy Payload", "Export Audit Ledger"])

        with t1:
            st.subheader("Systemic Bioavailability & Decision Tree Log")

            if diep:
                st.markdown("#### 1. Deterministic Bounds & Cramer Decision Tree (DIEP-MoS)")
                st.write(f"**Canonical API Properties:** MW: `{diep.get('api_mw', 0.0):.2f}` g/mol | LogP: `{diep.get('api_logp', 0.0):.2f}`")
                st.write(f"**Topological Classification:** `{diep.get('cramer_class', res.get('cramer_class', 'Class III'))}`")
                st.write(f"**Calculated Unionized Fraction ($f_{{ui}}$):** `{diep.get('f_ui', 1.0):.4f}`")
                st.write(f"**Fickian Dermal Absorption (DA%):** `{diep.get('da_pct_applied', 100.0):.2f}%`")
                st.write(f"**Max Systemic Exposure Dose (SED):** `{diep.get('sed_ug_day', res.get('sed_mg_kg_day', 0.0) * 1000.0):.2f} µg/day` vs TTC Limit `{diep.get('ttc_limit_ug', 1800.0):.2f} µg/day`")

                st.markdown("**Computational Cramer Tree Execution Trail:**")
                cramer_logs = diep.get("cramer_tree_log", res.get("audit_logs", ["Tree execution complete."]))
                for step in cramer_logs:
                    st.text(f"  └── {step}")

                status_val = diep.get('status', 'PASS' if res.get('is_compliant', True) else 'FAIL')
                if status_val == "FAIL":
                    st.error("❌ **DETERMINISTIC FAILURE:** Absolute systemic exposure exceeds safe EFSA thresholds.")
                else:
                    st.success("✓ **DETERMINISTIC PASS:** Systemic exposure is within safe limits.")

            st.markdown("#### 2. Probabilistic Exposure (Monte Carlo)")
            mos_data = data.get("mos", {})
            if isinstance(mos_data, dict):
                st.write(f"**Median MoS:** `{mos_data.get('median_mos', res.get('mos', 'N/A'))}` | **Failure Prob:** `{mos_data.get('failure_probability', 0.0)*100:.2f}%` against AF target of `{mos_data.get('target_af', 100)}`")
                if "af_breakdown" in mos_data:
                    st.json(mos_data["af_breakdown"])
            else:
                st.write(f"**Margin of Safety (MoS):** `{mos_data}`")

        with t2:
            st.subheader("TRACE-Onco / VMTB Output Vector")
            if diep:
                st.json({
                    "patient_hepatic_burden_ratio": diep.get("hepatic_burden_ratio", 0.12),
                    "oncogenic_risk_index": diep.get("oncogenic_risk_index", 0.01),
                    "bioavailability_status": diep.get("status", "PASS"),
                    "structural_alerts": [alert for alert in data.get("alerts", [])],
                    "recommendation": "Integrate ratio directly into decentralized Lifelines Cox-PH model, dynamically weighted against patient baseline De Ritis ratio to account for hepatic stress."
                })
            else:
                st.info("Enable Deterministic Biophysics to generate the downstream clinical integration payload.")

        with t3:
            st.subheader("Generate & Download PDF Ledger")
            try:
                pdf_buf = generate_enterprise_pdf(data, diep)
                pdf_bytes = pdf_buf.getvalue() if hasattr(pdf_buf, "getvalue") else pdf_buf
                st.download_button(
                    "📥 Download Multi-Section Computational Assessment Dossier (PDF)",
                    data=pdf_bytes,
                    file_name=f"{assessment_id}.pdf",
                    mime="application/pdf",
                    type="primary"
                )
            except Exception as e:
                st.error(f"PDF Generation Error: {str(e)}")

        # -------------------------------------------------------------
        # DOMAIN 2: FSANZ STANDARD 2.9.4 SPORTS DRINK ENGINE
        # -------------------------------------------------------------
        elif engine_id_str in ["fsanz_294_sports_drink", "fsanz_engine"]:
            overall_status = res.get("overall_status", "PASS" if res.get("is_compliant", True) or res.get("status") == "PASS" else "FAIL")
            osmolality = res.get("osmolality_mOsm_kg", 0.0)
            sodium_mmol = res.get("sodium_mmol_l", res.get("na_mmol_l", 0.0))
            osmo_class = res.get("osmo_classification", res.get("tonicity_classification", "Isotonic"))

            if "PASS" in overall_status:
                st.success(f"### ✓ FSANZ 2.9.4 COMPLIANCE VERDICT: {overall_status}")
            elif "WARNING" in overall_status:
                st.warning(f"### ⚠️ FSANZ 2.9.4 COMPLIANCE VERDICT: {overall_status}")
            else:
                st.error(f"### ❌ FSANZ 2.9.4 COMPLIANCE VERDICT: {overall_status}")

            m1, m2, m3 = st.columns(3)
            m1.metric("Calculated Osmolality", f"{osmolality:.1f} mOsm/kg")
            m2.metric("Prepared Sodium Concentration", f"{sodium_mmol:.2f} mmol/L")
            m3.metric("Hydration Profile", osmo_class)

            st.markdown("#### 1. Schedule 29 Active Yield & Dosage Audit Matrix")
            raw_audit_table = res.get("audit_table", res.get("parsed_compounds", res.get("compounds", [])))
            if raw_audit_table:
                df_audit = pd.DataFrame(raw_audit_table)
                st.dataframe(df_audit, use_container_width=True)

            st.markdown("#### 2. Mandatory Package Label Warning Statements")
            
            # Robust key lookup across all possible return dictionary keys
            warnings_list = (
                res.get("mandatory_warnings") 
                or res.get("warnings") 
                or res.get("audit_logs") 
                or res.get("label_warnings") 
                or []
            )

            if warnings_list:
                for warning in warnings_list:
                    if "FAIL" in warning or "NON-COMPLIANT" in warning:
                        st.error(f"❌ **NON-COMPLIANCE ALERT:** {warning}")
                    elif "PASS" in warning or "VERIFIED" in warning:
                        st.success(f"✓ **VERIFIED STATEMENT:** {warning}")
                    else:
                        st.warning(f"⚠️ **REQUIRED STATEMENT:** {warning}")
            else:
                # Explicit fallback when no warnings are generated
                st.info("ℹ️ No mandatory label warning statements triggered for this formulation.")

            st.markdown("#### 3. Factory Specification PDF Generation")
            try:
                if "pdf_bytes" in res and isinstance(res["pdf_bytes"], bytes):
                    pdf_bytes = res["pdf_bytes"]
                else:
                    pdf_buffer = build_factory_spec_pdf(res)
                    pdf_bytes = pdf_buffer.getvalue() if hasattr(pdf_buffer, "getvalue") else pdf_buffer

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
                st.error(f"Failed to generate FSANZ PDF specification: {str(e)}")
                
