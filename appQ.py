"""
Streamlit Multi-Domain In-Silico Computational Platform
Integrates Cutaneous Bioactivation (IATA/DIEP-MoS/TRACE-Onco) and 
FSANZ Standard 2.9.4 Food Science Engines via Engine Registry Architecture.
"""
import streamlit as st
import pandas as pd
from io import BytesIO

# Import existing Cutaneous Core Modules (UNTOUCHED)
from tox_engine import execute_full_compound_audit, generate_enterprise_pdf
from diep_engine import (
    run_diep_gatekeeper,
    SCCS_MECHANISTIC_REGISTRY,
    ToxicophoreMatchException,
    PhysicochemicalConstraintError
)

# Import FSANZ 2.9.4 Core Modules
from fsanz_294_engine import (
    StickPackFormula,
    FormulaIngredient,
    FSANZ294Auditor,
    build_factory_spec_pdf
)

# Import Strategy Registry Architecture
from engine_registry import BaseComputationalEngine, REGISTRY


# =====================================================================
# 1. ENGINE ADAPTER WRAPPERS (ISOLATED MATHEMATICAL DOMAINS)
# =====================================================================

class CutaneousBioactivationEngine(BaseComputationalEngine):
    @property
    def engine_id(self) -> str:
        return "cutaneous_iata_diep"

    @property
    def engine_name(self) -> str:
        return "In-Silico Cutaneous Bioactivation & Toxicological Screener"

    @property
    def domain_category(self) -> str:
        return "Cutaneous Biophysics & Dermal Toxicology"


class FSANZSportsFoodEngine(BaseComputationalEngine):
    @property
    def engine_id(self) -> str:
        return "fsanz_294_sports_drink"

    @property
    def engine_name(self) -> str:
        return "FSANZ Standard 2.9.4 Sports Powder Auditor"

    @property
    def domain_category(self) -> str:
        return "Food Science & Oral Nutraceuticals"


# Register Engines into Strategy Architecture
REGISTRY.register(CutaneousBioactivationEngine())
REGISTRY.register(FSANZSportsFoodEngine())


# =====================================================================
# 2. STREAMLIT UI CONFIGURATION & DOMAIN ROUTING
# =====================================================================

st.set_page_config(page_title="Multi-Domain In-Silico Screener", layout="wide", page_icon="🧬")

st.sidebar.title("Domain Engine Registry")
st.sidebar.markdown("Select computational engine for complete domain isolation:")

registered_engines = REGISTRY.list_engines()
engine_options = {f"{e['domain']} → {e['name']}": e['id'] for e in registered_engines}

selected_label = st.sidebar.selectbox("Active Domain Engine:", list(engine_options.keys()))
active_engine_id = engine_options[selected_label]
active_engine = REGISTRY.get(active_engine_id)

st.sidebar.markdown("---")
st.sidebar.info(f"**Active Engine ID:** `{active_engine.engine_id}`\n\n**Domain Category:** `{active_engine.domain_category}`")


# =====================================================================
# DOMAIN 1: IN-SILICO CUTANEOUS BIOACTIVATION & TOXICOLOGICAL SCREENER
# =====================================================================
if active_engine.engine_id == "cutaneous_iata_diep":
    st.title("In-Silico Cutaneous Bioactivation & Toxicological Screener")
    st.caption("Integrates TRACE-Onco compatible Deterministic Biophysics (DIEP-MoS) & Multiprotic Ionization.")

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
        
        st.markdown("#### Multiprotic Ionization Centers")
        acid_pkas_str = st.text_input("Acidic pKa values (comma-separated)", "4.2")
        base_pkas_str = st.text_input("Basic pKa values (comma-separated)", "")
        
        da = st.slider("Mean Dermal Absorption (%)", 0.1, 100.0, 50.0, disabled=use_deterministic_flux)
        bw = st.number_input("Body Weight (kg)", 10.0, 150.0, 60.0)
        mc_samples = int(st.number_input("Monte Carlo Iterations", min_value=1000, max_value=100000, value=10000, step=1000))
        
        with st.expander("Advanced Dempster-Shafer Prior Calibration"):
            baseline_safe = st.slider("Baseline Safe Prior Mass", 0.1, 0.9, 0.70, 0.05)
            discount_rate = st.slider("Source Discount Rate", 0.1, 0.9, 0.50, 0.05)

        run = st.button("Run Research Audit", type="primary")

    with col2:
        if run and smiles:
            with st.spinner("Executing decision tree and Monte Carlo simulation..."):
                acid_pkas = [float(x.strip()) for x in acid_pkas_str.split(",") if x.strip()]
                base_pkas = [float(x.strip()) for x in base_pkas_str.split(",") if x.strip()]
                
                diep_results = None
                if use_deterministic_flux:
                    try:
                        diep_results = run_diep_gatekeeper(smiles, conc, formulation_ph, acid_pkas, base_pkas, product_type)
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
                    mc_samples=mc_samples,
                    baseline_safe_belief=baseline_safe,
                    structural_discount_rate=discount_rate
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


# =====================================================================
# DOMAIN 2: FSANZ STANDARD 2.9.4 SPORTS POWDER FOOD SCIENCE ENGINE
# =====================================================================
elif active_engine.engine_id == "fsanz_294_sports_drink":
    st.title("FSANZ Standard 2.9.4 Powdered Sports Drink Auditor")
    st.caption("Evaluates powdered stick-pack formulations against Australian Schedule 29 limits, active yields, and mandatory labeling requirements.")

    c1, c2 = st.columns([1, 2])

    with c1:
        st.subheader("Product & Formulation Configuration")
        product_name = st.text_input("Product Name", "Electrolyte Hydration Stick Pack")
        flavor = st.text_input("Flavor Variant", "Lemon Lime")
        stick_weight_g = st.number_input("Stick Pack Net Weight (g)", min_value=1.0, max_value=50.0, value=7.5)
        serves_per_day = st.number_input("Max Daily Servings", min_value=1, max_value=6, value=2)

        st.markdown("---")
        st.markdown("### Active Ingredients Configuration")
        
        # Salt raw materials with yield conversions
        sod_citrate = st.number_input("Sodium Citrate Dihydrate (mg/serve)", value=1000.0)
        pot_chloride = st.number_input("Potassium Chloride (mg/serve)", value=300.0)
        mag_glycinate = st.number_input("Magnesium Glycinate (mg/serve)", value=400.0)
        vit_c = st.number_input("Ascorbic Acid / Vitamin C (mg/serve)", value=45.0)
        caffeine = st.number_input("Natural Caffeine (mg/serve)", value=35.0)

        run_fsanz = st.button("Run FSANZ 2.9.4 Compliance Audit", type="primary")

    with c2:
        if run_fsanz:
            ingredients_list = [
                FormulaIngredient("Sodium Citrate Dihydrate", sod_citrate, "Sodium", 23.5),
                FormulaIngredient("Potassium Chloride", pot_chloride, "Potassium", 52.4),
                FormulaIngredient("Magnesium Glycinate", mag_glycinate, "Magnesium", 14.1),
                FormulaIngredient("Ascorbic Acid", vit_c, "Vitamin C", 100.0),
                FormulaIngredient("Natural Caffeine", caffeine, "Caffeine", 100.0),
            ]

            formula = StickPackFormula(
                product_name=product_name,
                stick_pack_weight_g=stick_weight_g,
                serves_per_day=serves_per_day,
                flavor_variant=flavor,
                ingredients=ingredients_list
            )

            auditor = FSANZ294Auditor(formula)
            audit_output = auditor.execute_audit()
            st.session_state["fsanz_results"] = audit_output

        if "fsanz_results" in st.session_state:
            res = st.session_state["fsanz_results"]

            if "PASS" in res["overall_status"]:
                st.success(f"✓ **FSANZ 2.9.4 COMPLIANCE VERDICT:** {res['overall_status']}")
            else:
                st.error(f"❌ **FSANZ 2.9.4 COMPLIANCE VERDICT:** {res['overall_status']}")

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
