import streamlit as st
import pandas as pd
from io import BytesIO
from dataclasses import dataclass
from typing import List, Dict, Any
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, HRFlowable, KeepTogether
from reportlab.pdfgen import canvas
from engines.base_engine import BaseComputationalEngine

FSANZ_SCHEDULE_29_LIMITS = {
    "Sodium": {"max_daily": 520.0, "unit": "mg"},
    "Potassium": {"max_daily": 1300.0, "unit": "mg"},
    "Magnesium": {"max_daily": 320.0, "unit": "mg"},
    "Vitamin C": {"max_daily": 100.0, "unit": "mg"},
    "Caffeine": {"max_daily": 100.0, "unit": "mg"}
}

@dataclass
class FormulaIngredient:
    name: str
    amount_per_serve_mg: float
    active_compound: str
    active_yield_pct: float

class FSANZ294Engine(BaseComputationalEngine):
    @property
    def engine_id(self) -> str:
        return "fsanz_294_sports_drink"

    @property
    def engine_name(self) -> str:
        return "FSANZ Standard 2.9.4 Sports Drink Auditor"

    @property
    def domain_category(self) -> str:
        return "Food Science & Oral Nutraceuticals"

    def render_inputs(self, st_ctx: Any) -> Dict[str, Any]:
        st_ctx.subheader("Product Parameters")
        product_name = st_ctx.text_input("Product Name", "Hydration Electrolyte Powder")
        stick_weight = st_ctx.number_input("Stick Pack Weight (g)", value=7.5)
        serves_day = st_ctx.number_input("Servings / Day", min_value=1, max_value=6, value=2)

        st_ctx.markdown("#### Formula Ingredients (mg/serve)")
        sod_citrate = st_ctx.number_input("Sodium Citrate Dihydrate (mg)", value=1000.0)
        pot_chloride = st_ctx.number_input("Potassium Chloride (mg)", value=300.0)
        mag_glycinate = st_ctx.number_input("Magnesium Glycinate (mg)", value=400.0)
        vit_c = st_ctx.number_input("Ascorbic Acid (mg)", value=45.0)
        caffeine = st_ctx.number_input("Natural Caffeine (mg)", value=35.0)

        return {
            "product_name": product_name, "stick_weight": stick_weight, "serves_day": serves_day,
            "ingredients": [
                FormulaIngredient("Sodium Citrate Dihydrate", sod_citrate, "Sodium", 23.5),
                FormulaIngredient("Potassium Chloride", pot_chloride, "Potassium", 52.4),
                FormulaIngredient("Magnesium Glycinate", mag_glycinate, "Magnesium", 14.1),
                FormulaIngredient("Ascorbic Acid", vit_c, "Vitamin C", 100.0),
                FormulaIngredient("Natural Caffeine", caffeine, "Caffeine", 100.0),
            ]
        }

    def execute(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        serves = inputs["serves_day"]
        audit_table = []
        overall_status = "PASS"
        warnings = set()

        for ing in inputs["ingredients"]:
            active_mg_serve = ing.amount_per_serve_mg * (ing.active_yield_pct / 100.0)
            active_mg_daily = active_mg_serve * serves
            limits = FSANZ_SCHEDULE_29_LIMITS.get(ing.active_compound, {})
            max_daily = limits.get("max_daily", 99999.0)

            status = "PASS"
            if active_mg_daily > max_daily:
                status = "FAIL"
                overall_status = "FAIL"
            
            if ing.active_compound == "Caffeine":
                warnings.add("Contains caffeine. Not recommended for children or pregnant women.")
                if active_mg_serve > 80.0:
                    status = "FAIL (Exceeds 80mg single serve limit)"
                    overall_status = "FAIL"

            audit_table.append({
                "Compound": ing.active_compound,
                "Per Serve": f"{active_mg_serve:.1f} mg",
                "Daily Intake": f"{active_mg_daily:.1f} mg",
                "Schedule 29 Cap": f"{max_daily} mg",
                "Status": status
            })

        return {
            "product_name": inputs["product_name"],
            "overall_status": overall_status,
            "audit_table": audit_table,
            "warnings": list(warnings)
        }
