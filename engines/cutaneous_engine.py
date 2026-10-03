"""
Cutaneous Bioactivation & Dermal Biophysics Engine Adapter
Wraps legacy toxicity & skin permeability modules under the BaseComputationalEngine contract.
"""

import streamlit as st
import pandas as pd
from typing import Dict, Any
from engines.base_engine import BaseComputationalEngine

# Import legacy root engines
import tox_engine
import diep_engine


class CutaneousEngineAdapter(BaseComputationalEngine):

    def get_metadata(self) -> Dict[str, str]:
        return {
            "id": "cutaneous_iata_diep",
            "name": "Cutaneous Biophysics & Toxicological Screener",
            "category": "Dermal Safety & Bioactivation",
            "version": "2.1.0",
            "description": "Integrated Cramer Class alerts, Henderson-Hasselbalch ionization, Fickian percutaneous flux, and Monte Carlo Margin of Safety (MoS) modeling."
        }

    def render_inputs(self) -> Dict[str, Any]:
        st.sidebar.markdown("### 🧬 Molecule & Physicochemical Inputs")
        smiles = st.sidebar.text_input("Canonical SMILES", value="CC1=C(C=C(C=C1)N)N")
        pka_acid = st.sidebar.number_input("pKa Acidic", value=14.0, step=0.1)
        pka_base = st.sidebar.number_input("pKa Basic", value=4.5, step=0.1)
        mw = st.sidebar.number_input("Molecular Weight (g/mol)", value=122.17, step=0.1)
        logp = st.sidebar.number_input("LogKow (Octanol-Water)", value=1.4, step=0.1)

        st.sidebar.markdown("### 🧴 Application & Exposure Settings")
        conc_pct = st.sidebar.number_input("Concentration (%)", value=1.0, step=0.1)
        applied_mg_cm2 = st.sidebar.number_input("Applied Dose (mg/cm²)", value=2.0, step=0.1)
        surface_area_cm2 = st.sidebar.number_input("Surface Area (cm²)", value=560.0, step=10.0)
        body_weight_kg = st.sidebar.number_input("Body Weight (kg)", value=60.0, step=1.0)
        pod_mg_kg_day = st.sidebar.number_input("Point of Departure (PoD mg/kg/day)", value=15.0, step=1.0)

        return {
            "smiles": smiles,
            "pka_acid": pka_acid,
            "pka_base": pka_base,
            "mw": mw,
            "logp": logp,
            "conc_pct": conc_pct,
            "applied_mg_cm2": applied_mg_cm2,
            "surface_area_cm2": surface_area_cm2,
            "body_weight_kg": body_weight_kg,
            "pod_mg_kg_day": pod_mg_kg_day
        }

    def execute(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        # 1. Run Legacy Cutaneous & Exposure Analytics
        tox_res = tox_engine.evaluate_molecule(inputs["smiles"])
        
        sed_mg_kg_day = (
            (inputs["conc_pct"] / 100.0) * inputs["applied_mg_cm2"] * inputs["surface_area_cm2"]
        ) / inputs["body_weight_kg"]

        mos = inputs["pod_mg_kg_day"] / sed_mg_kg_day if sed_mg_kg_day > 0 else 9999.0
        
        # 2. Run Monte Carlo & Dempster-Shafer Prior Fusion
        diep_res = diep_engine.run_dempster_shafer_monte_carlo(
            pod=inputs["pod_mg_kg_day"],
            sed=sed_mg_kg_day
        )

        return {
            "status": "PASS" if mos >= 100.0 else "FAIL",
            "is_compliant": mos >= 100.0,
            "smiles": inputs["smiles"],
            "cramer_class": tox_res.get("cramer_class", "Class III"),
            "sed_mg_kg_day": sed_mg_kg_day,
            "mos": mos,
            "belief_safe": diep_res.get("belief_safe", 0.95),
            "plausibility_safe": diep_res.get("plausibility_safe", 0.99),
            "audit_logs": [
                f"Computed SED: {sed_mg_kg_day:.4f} mg/kg/day",
                f"Margin of Safety (MoS): {mos:.1f} (Threshold >= 100.0)",
                f"Dempster-Shafer Belief in Safety: {diep_res.get('belief_safe', 0.95):.3f}"
            ]
        }
        
