import streamlit as st
from typing import Dict, Any
from engines.base_engine import BaseComputationalEngine

# Direct imports from root-level scripts
from tox_engine import execute_full_compound_audit, generate_enterprise_pdf
from diep_engine import (
    run_diep_gatekeeper,
    SCCS_MECHANISTIC_REGISTRY,
    ToxicophoreMatchException
)

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

    def render_inputs(self, st_ctx: Any) -> Dict[str, Any]:
        st_ctx.subheader("Formulation Inputs")
        smiles = st_ctx.text_input("Target SMILES", "CC(=O)Oc1ccccc1C(=O)O")
        product_type = st_ctx.selectbox("Product Scenario", list(SCCS_MECHANISTIC_REGISTRY.keys()), index=0)
        conc = st_ctx.number_input("Concentration (%)", min_value=0.01, max_value=100.0, value=2.0)
        
        st_ctx.markdown("#### Toxicological Point of Departure (PoD)")
        pod = st_ctx.number_input("PoD Value (mg/kg/day)", min_value=0.1, max_value=10000.0, value=250.0)
        
        pod_col1, pod_col2, pod_col3 = st_ctx.columns(3)
        with pod_col1:
            pod_type = st_ctx.selectbox("PoD Type", ["NOAEL", "LOAEL", "BMDL"], index=0)
        with pod_col2:
            species = st_ctx.selectbox("Test Species", ["Rat", "Mouse", "Dog", "Rabbit", "Human"], index=0)
        with pod_col3:
            duration = st_ctx.selectbox("Study Duration", ["Chronic", "Subchronic", "Subacute"], index=1)
        
        st_ctx.markdown("---")
        st_ctx.markdown("### Deterministic Biophysics (DIEP-MoS)")
        use_deterministic_flux = st_ctx.checkbox("Enable Fickian Deterministic Flux", value=True)
        formulation_ph = st_ctx.number_input("Formulation pH", min_value=0.0, max_value=14.0, value=5.5)
        
        st_ctx.markdown("#### Multiprotic Ionization Centers")
        acid_pkas_str = st_ctx.text_input("Acidic pKa values (comma-separated)", "4.2")
        base_pkas_str = st_ctx.text_input("Basic pKa values (comma-separated)", "")
        
        da = st_ctx.slider("Mean Dermal Absorption (%)", 0.1, 100.0, 50.0, disabled=use_deterministic_flux)
        bw = st_ctx.number_input("Body Weight (kg)", 10.0, 150.0, 60.0)
        mc_samples = int(st_ctx.number_input("Monte Carlo Iterations", min_value=1000, max_value=100000, value=10000, step=1000))
        
        with st_ctx.expander("Advanced Dempster-Shafer Prior Calibration"):
            baseline_safe = st_ctx.slider("Baseline Safe Prior Mass", 0.1, 0.9, 0.70, 0.05)
            discount_rate = st_ctx.slider("Source Discount Rate", 0.1, 0.9, 0.50, 0.05)

        return {
            "smiles": smiles, "product_type": product_type, "conc": conc,
            "pod": pod, "pod_type": pod_type, "species": species, "duration": duration,
            "use_deterministic_flux": use_deterministic_flux, "formulation_ph": formulation_ph,
            "acid_pkas_str": acid_pkas_str, "base_pkas_str": base_pkas_str,
            "da": da, "bw": bw, "mc_samples": mc_samples,
            "baseline_safe": baseline_safe, "discount_rate": discount_rate
        }

    def execute(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        acid_pkas = [float(x.strip()) for x in inputs["acid_pkas_str"].split(",") if x.strip()]
        base_pkas = [float(x.strip()) for x in inputs["base_pkas_str"].split(",") if x.strip()]
        
        diep_results = None
        da = inputs["da"]

        if inputs["use_deterministic_flux"]:
            diep_results = run_diep_gatekeeper(
                inputs["smiles"], inputs["conc"], inputs["formulation_ph"], 
                acid_pkas, base_pkas, inputs["product_type"]
            )
            da = diep_results["da_pct_applied"]

        audit_res = execute_full_compound_audit(
            smiles=inputs["smiles"], product_type=inputs["product_type"],
            concentration_pct=inputs["conc"], pod_noael_mg_kg_day=inputs["pod"],
            pod_type=inputs["pod_type"], species=inputs["species"],
            duration=inputs["duration"], dermal_absorption_pct=da,
            body_weight_kg=inputs["bw"], mc_samples=inputs["mc_samples"],
            baseline_safe_belief=inputs["baseline_safe"],
            structural_discount_rate=inputs["discount_rate"]
        )

        return {"audit": audit_res, "diep": diep_results}
      
