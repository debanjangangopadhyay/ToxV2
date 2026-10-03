"""
Cutaneous Bioactivation & Dermal Biophysics Engine Adapter
Full feature-parity wrapper for diep_engine and tox_engine.
"""

import streamlit as st
import pandas as pd
from typing import Dict, Any
from engines.base_engine import BaseComputationalEngine

from tox_engine import execute_full_compound_audit
from diep_engine import (
    run_diep_gatekeeper,
    SCCS_MECHANISTIC_REGISTRY,
    ToxicophoreMatchException
)


class CutaneousEngineAdapter(BaseComputationalEngine):

    @property
    def engine_id(self) -> str:
        return "cutaneous_iata_diep"

    @property
    def engine_name(self) -> str:
        return "Cutaneous Biophysics & Toxicological Screener"

    @property
    def domain_category(self) -> str:
        return "Dermal Safety & Bioactivation"

    def get_metadata(self) -> Dict[str, str]:
        return {
            "id": self.engine_id,
            "name": self.engine_name,
            "category": self.domain_category,
            "version": "2.1.0",
            "description": "IATA Computational Screener integrating TRACE-Onco Deterministic Biophysics (DIEP-MoS), Multiprotic Ionization, and Monte Carlo Probabilistic Risk."
        }

    def render_inputs(self, st_ctx: Any = st) -> Dict[str, Any]:
        """Renders all original IATA input fields in the sidebar."""
        st_ctx.sidebar.subheader("Formulation Inputs")
        smiles = st_ctx.sidebar.text_input("Target SMILES", "CC(=O)Oc1ccccc1C(=O)O")
        product_type = st_ctx.sidebar.selectbox("Product Scenario", list(SCCS_MECHANISTIC_REGISTRY.keys()), index=0)
        conc = st_ctx.sidebar.number_input("Concentration (%)", min_value=0.01, max_value=100.0, value=2.0)

        st_ctx.sidebar.markdown("#### Toxicological Point of Departure (PoD)")
        pod = st_ctx.sidebar.number_input("PoD Value (mg/kg/day)", min_value=0.1, max_value=10000.0, value=250.0)
        
        pod_type = st_ctx.sidebar.selectbox("PoD Type", ["NOAEL", "LOAEL", "BMDL"], index=0)
        species = st_ctx.sidebar.selectbox("Test Species", ["Rat", "Mouse", "Dog", "Rabbit", "Human"], index=0)
        duration = st_ctx.sidebar.selectbox("Study Duration", ["Chronic", "Subchronic", "Subacute"], index=1)

        st_ctx.sidebar.markdown("---")
        st_ctx.sidebar.markdown("### Deterministic Biophysics (DIEP-MoS)")
        use_deterministic_flux = st_ctx.sidebar.checkbox("Enable Fickian Deterministic Flux", value=True)
        formulation_ph = st_ctx.sidebar.number_input("Formulation pH", min_value=0.0, max_value=14.0, value=5.5)

        st_ctx.sidebar.markdown("#### Multiprotic Ionization Centers")
        acid_pkas_str = st_ctx.sidebar.text_input("Acidic pKa values (comma-separated)", "4.2")
        base_pkas_str = st_ctx.sidebar.text_input("Basic pKa values (comma-separated)", "")

        da = st_ctx.sidebar.slider("Mean Dermal Absorption (%)", 0.1, 100.0, 50.0, disabled=use_deterministic_flux)
        bw = st_ctx.sidebar.number_input("Body Weight (kg)", 10.0, 150.0, 60.0)
        mc_samples = int(st_ctx.sidebar.number_input("Monte Carlo Iterations", min_value=1000, max_value=100000, value=10000, step=1000))

        with st_ctx.sidebar.expander("Advanced Dempster-Shafer Prior Calibration"):
            baseline_safe = st_ctx.sidebar.slider("Baseline Safe Prior Mass", 0.1, 0.9, 0.70, 0.05)
            discount_rate = st_ctx.sidebar.slider("Source Discount Rate", 0.1, 0.9, 0.50, 0.05)

        return {
            "smiles": smiles,
            "product_type": product_type,
            "conc": conc,
            "pod": pod,
            "pod_type": pod_type,
            "species": species,
            "duration": duration,
            "use_deterministic_flux": use_deterministic_flux,
            "ph": formulation_ph,
            "acid_pkas_str": acid_pkas_str,
            "base_pkas_str": base_pkas_str,
            "da": da,
            "bw": bw,
            "mc_samples": mc_samples,
            "baseline_safe": baseline_safe,
            "discount_rate": discount_rate
        }

    def execute(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Executes the dual DIEP biophysics gatekeeper and full compound audit."""
        smiles = inputs["smiles"]
        product_type = inputs["product_type"]
        conc = inputs["conc"]
        pod = inputs["pod"]
        pod_type = inputs["pod_type"]
        species = inputs["species"]
        duration = inputs["duration"]
        use_deterministic_flux = inputs["use_deterministic_flux"]
        ph = inputs["ph"]
        da = inputs["da"]
        bw = inputs["bw"]
        mc_samples = inputs["mc_samples"]
        baseline_safe = inputs["baseline_safe"]
        discount_rate = inputs["discount_rate"]

        acid_pkas = [float(x.strip()) for x in inputs["acid_pkas_str"].split(",") if x.strip()]
        base_pkas = [float(x.strip()) for x in inputs["base_pkas_str"].split(",") if x.strip()]

        diep_results = None
        if use_deterministic_flux:
            # Propagates ToxicophoreMatchException directly to app.py handler
            diep_results = run_diep_gatekeeper(
                smiles=smiles,
                conc_pct=conc,
                ph=ph,
                acid_pkas=acid_pkas,
                base_pkas=base_pkas,
                product_type=product_type
            )
            da = diep_results["da_pct_applied"]

        audit_res = execute_full_compound_audit(
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

        if not audit_res.get("valid", False):
            raise ValueError(f"Audit Execution Failed: {audit_res.get('error')}")

        return {
            "audit": audit_res,
            "diep": diep_results
        }
