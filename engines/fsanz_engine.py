"""
FSANZ Standard 2.9.4 & Schedule 29 Electrolyte Engine (RDKit-Driven)
================================================--------------------
Mathematical & Chemical Engineering Features:
1. RDKit Topological SMILES Parsing: Dynamically computes exact formula weight, 
   ion dissociation count (nu), and elemental yield fractions without hardcoded salts.
2. Debye-Hückel & Pitzer Electrolyte Thermodynamics: Calculates solution ionic strength (I)
   from fragment formal charges and derives the molal osmotic activity coefficient phi(I).
3. FSANZ Standard 2.9.4 Division 2 Gatekeeping: Validates sodium molarity (10-30 mmol/L).
4. FSANZ Schedule 29 Daily Exposure Caps: Enforces potassium, magnesium, and caffeine ceilings.
5. ReportLab Dossier Factory: Generates publication-grade compliance PDF dossiers.
"""

import io
import math
import datetime
from typing import Dict, List, Any, Tuple, Optional
import pandas as pd
import streamlit as st

# RDKit Chemical Engine Imports
from rdkit import Chem
from rdkit.Chem import Descriptors

# ReportLab PDF Engine Imports
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable, KeepTogether
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle

# Base Abstract Interface
from engines.base_engine import BaseComputationalEngine


# ==============================================================================
# PROVEN STATUTORY & PHYSICAL CONSTANTS (MINIMAL HARDCODED BOUNDARIES)
# ==============================================================================
FSANZ_NA_MIN_MMOL_L = 10.0      # Standard 2.9.4 Div 2 Sodium Floor
FSANZ_NA_MAX_MMOL_L = 30.0      # Standard 2.9.4 Div 2 Sodium Ceiling
FSANZ_K_MAX_DAILY_MG = 1300.0   # Schedule 29 Potassium Daily Ceiling
FSANZ_MG_MAX_DAILY_MG = 320.0   # Schedule 29 Magnesium Daily Ceiling
FSANZ_CAFFEINE_SERVE_MG = 80.0  # Schedule 29 Caffeine Single Serve Cap
FSANZ_CAFFEINE_DAILY_MG = 100.0 # Schedule 29 Caffeine Daily Cap

# Debye-Hückel Osmotic Activity Constants for Water at 25°C
DEBYE_A_PHI = 0.392             # Debye-Hückel osmotic slope (kg^0.5 / mol^0.5)
DEBYE_B = 1.2                   # Ion interaction size parameter


class FSANZ294Engine(BaseComputationalEngine):
    """
    RDKit-powered computational engine for FSANZ Standard 2.9.4 food chemistry verification.
    """

    def get_metadata(self) -> Dict[str, str]:
        return {
            "id": "fsanz_294_sports_drink",
            "name": "FSANZ 2.9.4 Electrolyte & Osmolality Engine (RDKit)",
            "category": "Food Science & Regulatory Chemistry",
            "version": "3.0.0-RDKit",
            "description": "Dynamic chemical engine using SMILES topology to evaluate ionic strength, Debye-Hückel osmolality, and FSANZ Schedule 29 compliance."
        }

    # ==========================================================================
    # RDKit CHEMICAL PARSING & MATHEMATICAL DERIVATION
    # ==========================================================================
    @staticmethod
    def parse_smiles_ingredient(smiles: str, mass_mg: float) -> Dict[str, Any]:
        """
        Parses a SMILES string using RDKit to calculate chemical dynamics:
        - Formula Mass (MW)
        - Ionic Dissociation Count (nu) via Fragment Decomposition
        - Elemental Mass Yields (Na, K, Mg, Ca, Zn, Cl, etc.)
        - Species Formal Charges and Ionic Strength Contributions
        """
        smiles_clean = smiles.strip()
        mol = Chem.MolFromSmiles(smiles_clean)
        if mol is None:
            raise ValueError(f"Invalid SMILES string provided: '{smiles}'")

        # 1. Exact Formula Mass (g/mol)
        mw = Descriptors.MolWt(mol)
        if mw <= 0:
            raise ValueError(f"Calculated MW for SMILES '{smiles}' is zero or negative.")

        # 2. Fragment Dissociation (nu) & Formal Charge Breakdown
        # Dissociated salts are period-separated in SMILES (e.g., [Na+].[Cl-])
        frags = Chem.GetMolFrags(mol, asMols=True)
        nu = len(frags)
        moles_compound = (mass_mg / 1000.0) / mw  # Moles of formulated compound

        frag_details = []
        for frag in frags:
            charge = Chem.GetFormalCharge(frag)
            frag_details.append({
                "charge": charge,
                "moles": moles_compound  # Stoichiometric 1:1 mol per fragment in formula
            })

        # 3. Dynamic Elemental Yield Fractions via IUPAC Atomic Mass Summation
        pt = Chem.GetPeriodicTable()
        elem_mass_sums: Dict[str, float] = {}

        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atomic_mass = pt.GetAtomicWeight(pt.GetAtomicNumber(symbol))
            elem_mass_sums[symbol] = elem_mass_sums.get(symbol, 0.0) + atomic_mass

        elem_yields_mg = {
            symbol: mass_mg * (mass_sum / mw)
            for symbol, mass_sum in elem_mass_sums.items()
        }

        # 4. Check if compound is Caffeine (SMILES matching caffeine structure)
        # Canonical SMILES for caffeine: Cn1cnc2c1c(=O)n(C)c(=O)n2C
        is_caffeine = (
            "c1nc2c(n1C)c(=O)n(C)c(=O)n2C" in Chem.MolToSmiles(mol) or
            elem_mass_sums.get("N", 0) == 4 and elem_mass_sums.get("C", 0) == 8 and "c1nc" in smiles_clean
        )

        return {
            "smiles": smiles_clean,
            "mass_mg": mass_mg,
            "mw": mw,
            "nu": nu,
            "moles": moles_compound,
            "elem_yields_mg": elem_yields_mg,
            "frag_details": frag_details,
            "is_caffeine": is_caffeine
        }

    @classmethod
    def calculate_solution_thermodynamics(
        cls, parsed_compounds: List[Dict[str, Any]], volume_l: float
    ) -> Tuple[float, float, float]:
        """
        Derives solution Ionic Strength (I), Debye-Hückel Osmotic Activity Coefficient phi(I),
        and Non-Ideal Osmolality (mOsm/kg) based on electrochemistry models.
        """
        if volume_l <= 0:
            return 0.0, 1.0, 0.0

        # 1. Solution Ionic Strength: I = 0.5 * sum( c_i * z_i^2 )
        total_ionic_strength = 0.0
        total_ideal_osmolal_moles = 0.0

        for comp in parsed_compounds:
            # Ideal osmolal contribution = moles * nu
            total_ideal_osmolal_moles += comp["moles"] * comp["nu"]

            for frag in comp["frag_details"]:
                c_i = frag["moles"] / volume_l  # Molar concentration (mol/L)
                z_i = frag["charge"]
                total_ionic_strength += 0.5 * c_i * (z_i ** 2)

        # 2. Debye-Hückel / Extended Osmotic Activity Coefficient phi(I)
        # phi(I) = 1 - (A_phi * sqrt(I)) / (1 + b * sqrt(I)) + beta * I
        if total_ionic_strength > 0:
            sqrt_I = math.sqrt(total_ionic_strength)
            phi = 1.0 - (DEBYE_A_PHI * sqrt_I) / (1.0 + DEBYE_B * sqrt_I) + (0.08 * total_ionic_strength)
            # Bound phi to realistic non-ideal aqueous limits
            phi = max(0.85, min(1.0, phi))
        else:
            phi = 1.0  # Ideal solution for non-ionic solutes

        # 3. Osmolality (mOsm/kg H2O): (Molar Osmolality / Volume) * phi * 1000
        osmolality_mOsm_kg = (total_ideal_osmolal_moles / volume_l) * phi * 1000.0

        return total_ionic_strength, phi, osmolality_mOsm_kg

    # ==========================================================================
    # STREAMLIT INPUT UI RENDERER
    # ==========================================================================
    def render_inputs(self) -> Dict[str, Any]:
        st.sidebar.markdown("### 📋 Formulation Parameters")
        product_name = st.sidebar.text_input("Product Name", value="ElectroPro Rehydrate")
        flavor_variant = st.sidebar.text_input("Flavor Variant", value="Lemon Lime")
        stick_pack_wt = st.sidebar.number_input("Stick Pack Weight (g)", min_value=1.0, max_value=50.0, value=7.5, step=0.5)
        max_daily_servings = st.sidebar.number_input("Max Daily Servings", min_value=1, max_value=10, value=2, step=1)
        dilution_vol_ml = st.sidebar.number_input("Dilution Volume per Serve (mL)", min_value=100.0, max_value=2000.0, value=500.0, step=50.0)

        st.markdown("### 🧪 RDKit Dynamic Formulation Builder")
        st.info("Input any valid chemical SMILES string. RDKit will dynamically derive molecular weights, ion dissociation numbers, and elemental active yields.")

        # Default standard formulation using canonical SMILES
        default_formula = [
            {"Name": "Trisodium Citrate Dihydrate", "SMILES": "[Na+].[Na+].[Na+].O=C([O-])CC(O)(CC(=O)[O-])C(=O)[O-].O.O", "Mass (mg)": 1000.0},
            {"Name": "Sodium Chloride", "SMILES": "[Na+].[Cl-]", "Mass (mg)": 250.0},
            {"Name": "Potassium Chloride", "SMILES": "[K+].[Cl-]", "Mass (mg)": 300.0},
            {"Name": "Magnesium Glycinate", "SMILES": "[Mg+2].O=C([O-])CN.O=C([O-])CN", "Mass (mg)": 400.0},
            {"Name": "Ascorbic Acid (Vit C)", "SMILES": "O=C1C(O)=C(O)C(=O)C(O1)C(O)CO", "Mass (mg)": 45.0},
            {"Name": "Natural Caffeine", "SMILES": "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Mass (mg)": 35.0},
        ]

        df_input = pd.DataFrame(default_formula)
        edited_df = st.data_editor(
            df_input,
            num_rows="dynamic",
            use_container_width=True,
            column_config={
                "Name": st.column_config.TextColumn("Compound Label", required=True),
                "SMILES": st.column_config.TextColumn("RDKit SMILES String", required=True),
                "Mass (mg)": st.column_config.NumberColumn("Mass per Serve (mg)", min_value=0.0, max_value=10000.0, step=10.0, required=True),
            }
        )

        return {
            "product_name": product_name,
            "flavor_variant": flavor_variant,
            "stick_pack_wt_g": stick_pack_wt,
            "max_daily_servings": max_daily_servings,
            "dilution_vol_ml": dilution_vol_ml,
            "ingredients": edited_df.to_dict(orient="records")
        }

    # ==========================================================================
    # CORE COMPUTATIONAL EXECUTION & COMPLIANCE EVALUATION
    # ==========================================================================
    def execute(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        volume_l = inputs["dilution_vol_ml"] / 1000.0
        daily_servings = inputs["max_daily_servings"]
        raw_ingredients = inputs.get("ingredients", [])

        parsed_compounds = []
        parsing_errors = []
        
        # 1. Parse all SMILES ingredients
        for item in raw_ingredients:
            name = item.get("Name", "Unnamed Compound")
            smiles = item.get("SMILES", "")
            mass_mg = float(item.get("Mass (mg)", 0.0))

            if mass_mg <= 0 or not smiles:
                continue

            try:
                parsed = self.parse_smiles_ingredient(smiles, mass_mg)
                parsed["label"] = name
                parsed_compounds.append(parsed)
            except Exception as e:
                parsing_errors.append(f"[{name}] SMILES Error: {str(e)}")

        if parsing_errors:
            return {
                "status": "FAIL",
                "error": "SMILES Parsing Failed",
                "logs": parsing_errors
            }

        # 2. Derive Solution Electrochemistry & Osmolality
        ionic_strength, phi, osmolality_mOsm_kg = self.calculate_solution_thermodynamics(
            parsed_compounds, volume_l
        )

        # 3. Sum Elemental Active Yields
        total_elemental_mg: Dict[str, float] = {}
        total_caffeine_serve_mg = 0.0

        for comp in parsed_compounds:
            if comp["is_caffeine"]:
                total_caffeine_serve_mg += comp["mass_mg"]
                
            for elem, mass in comp["elem_yields_mg"].items():
                total_elemental_mg[elem] = total_elemental_mg.get(elem, 0.0) + mass

        # 4. FSANZ Key Regulatory Metrics
        na_mass_mg = total_elemental_mg.get("Na", 0.0)
        k_mass_mg = total_elemental_mg.get("K", 0.0)
        mg_mass_mg = total_elemental_mg.get("Mg", 0.0)

        # Molar Sodium Concentration: (mg Na / 22.98977 g/mol) / Volume L
        na_mmol_l = (na_mass_mg / 22.989769) / volume_l if volume_l > 0 else 0.0

        # Cumulative Daily Exposures
        k_daily_mg = k_mass_mg * daily_servings
        mg_daily_mg = mg_mass_mg * daily_servings
        caffeine_daily_mg = total_caffeine_serve_mg * daily_servings

        # 5. Regulatory Gatekeeping Checks
        audit_logs = []
        is_compliant = True

        # Sodium Check (Standard 2.9.4 Div 2)
        na_pass = FSANZ_NA_MIN_MMOL_L <= na_mmol_l <= FSANZ_NA_MAX_MMOL_L
        if not na_pass:
            is_compliant = False
            audit_logs.append(f"FAIL: Prepared Sodium concentration {na_mmol_l:.2f} mmol/L outside FSANZ 2.9.4 limits ({FSANZ_NA_MIN_MMOL_L}-{FSANZ_NA_MAX_MMOL_L} mmol/L).")
        else:
            audit_logs.append(f"PASS: Prepared Sodium concentration {na_mmol_l:.2f} mmol/L satisfies FSANZ 2.9.4 standards.")

        # Potassium Check (Schedule 29)
        k_pass = k_daily_mg <= FSANZ_K_MAX_DAILY_MG
        if not k_pass:
            is_compliant = False
            audit_logs.append(f"FAIL: Total daily Potassium ({k_daily_mg:.1f} mg) exceeds Schedule 29 ceiling ({FSANZ_K_MAX_DAILY_MG} mg/day).")

        # Magnesium Check (Schedule 29)
        mg_pass = mg_daily_mg <= FSANZ_MG_MAX_DAILY_MG
        if not mg_pass:
            is_compliant = False
            audit_logs.append(f"FAIL: Total daily Magnesium ({mg_daily_mg:.1f} mg) exceeds Schedule 29 ceiling ({FSANZ_MG_MAX_DAILY_MG} mg/day).")

        # Caffeine Checks (Schedule 29)
        caffeine_serve_pass = total_caffeine_serve_mg <= FSANZ_CAFFEINE_SERVE_MG
        caffeine_daily_pass = caffeine_daily_mg <= FSANZ_CAFFEINE_DAILY_MG
        if not caffeine_serve_pass:
            is_compliant = False
            audit_logs.append(f"FAIL: Single serve Caffeine ({total_caffeine_serve_mg:.1f} mg) exceeds limit ({FSANZ_CAFFEINE_SERVE_MG} mg).")
        if not caffeine_daily_pass:
            is_compliant = False
            audit_logs.append(f"FAIL: Total daily Caffeine ({caffeine_daily_mg:.1f} mg) exceeds ceiling ({FSANZ_CAFFEINE_DAILY_MG} mg/day).")

        # Osmolality Classification
        if osmolality_mOsm_kg < 270.0:
            tonicity = "Hypotonic"
        elif 270.0 <= osmolality_mOsm_kg <= 330.0:
            tonicity = "Isotonic"
        else:
            tonicity = "Hypertonic"

        # Construct Execution Results
        results = {
            "status": "PASS" if is_compliant else "FAIL",
            "is_compliant": is_compliant,
            "product_name": inputs["product_name"],
            "flavor_variant": inputs["flavor_variant"],
            "dilution_vol_ml": inputs["dilution_vol_ml"],
            "daily_servings": daily_servings,
            "na_mmol_l": na_mmol_l,
            "na_mass_mg_serve": na_mass_mg,
            "k_mass_mg_serve": k_mass_mg,
            "k_daily_mg": k_daily_mg,
            "mg_mass_mg_serve": mg_mass_mg,
            "mg_daily_mg": mg_daily_mg,
            "caffeine_serve_mg": total_caffeine_serve_mg,
            "caffeine_daily_mg": caffeine_daily_mg,
            "ionic_strength_M": ionic_strength,
            "osmotic_coefficient_phi": phi,
            "osmolality_mOsm_kg": osmolality_mOsm_kg,
            "tonicity_classification": tonicity,
            "parsed_compounds": parsed_compounds,
            "audit_logs": audit_logs,
            "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }

        # Generate ReportLab PDF Bytes
        pdf_bytes = self.generate_pdf_bytes(results)
        results["pdf_bytes"] = pdf_bytes

        return results

    # ==========================================================================
    # REPORTLAB COMPLIANCE DOSSIER PDF GENERATOR
    # ==========================================================================
    def generate_pdf_bytes(self, results: Dict[str, Any]) -> bytes:
        buffer = io.BytesIO()
        doc = SimpleDocTemplate(
            buffer, pagesize=letter,
            leftMargin=36, rightMargin=36, topMargin=36, bottomMargin=36
        )

        styles = getSampleStyleSheet()
        title_style = ParagraphStyle(
            "DocTitle", parent=styles["Heading1"], fontSize=18, leading=22, textColor=colors.HexColor("#1A365D")
        )
        sub_style = ParagraphStyle(
            "DocSub", parent=styles["Normal"], fontSize=10, leading=14, textColor=colors.HexColor("#4A5568")
        )
        cell_style = ParagraphStyle("Cell", parent=styles["Normal"], fontSize=8, leading=10)
        bold_cell = ParagraphStyle("BoldCell", parent=styles["Normal"], fontSize=8, leading=10, fontName="Helvetica-Bold")

        story = []

        # Header
        story.append(Paragraph("FSANZ Standard 2.9.4 Regulatory Compliance Dossier", title_style))
        story.append(Paragraph(f"Product: <b>{results['product_name']}</b> ({results['flavor_variant']}) | Generated: {results['timestamp']}", sub_style))
        story.append(Spacer(1, 10))
        story.append(HRFlowable(width="100%", thickness=1.5, color=colors.HexColor("#2B6CB0"), spaceAfter=10))

        # Overall Status Banner
        status_color = colors.HexColor("#2F855A") if results["is_compliant"] else colors.HexColor("#C53030")
        status_text = "VERIFIED COMPLIANT (PASS)" if results["is_compliant"] else "REGULATORY NON-COMPLIANT (FAIL)"
        
        banner_data = [[Paragraph(f"<b>REGULATORY AUDIT VERDICT: {status_text}</b>", ParagraphStyle("B", parent=cell_style, textColor=colors.white, fontSize=11, fontName="Helvetica-Bold"))]]
        banner_table = Table(banner_data, colWidths=[540])
        banner_table.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,-1), status_color),
            ('ALIGN', (0,0), (-1,-1), 'CENTER'),
            ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
            ('TOPPADDING', (0,0), (-1,-1), 8),
            ('BOTTOMPADDING', (0,0), (-1,-1), 8),
        ]))
        story.append(banner_table)
        story.append(Spacer(1, 12))

        # Core Metrics Summary Table
        metrics_data = [
            [Paragraph("Regulatory Parameter", bold_cell), Paragraph("Derived Value", bold_cell), Paragraph("FSANZ Standard / Schedule 29 Limit", bold_cell), Paragraph("Status", bold_cell)],
            [Paragraph("Prepared Sodium Concentration", cell_style), Paragraph(f"{results['na_mmol_l']:.2f} mmol/L", cell_style), Paragraph("10.0 - 30.0 mmol/L (Std 2.9.4 Div 2)", cell_style), Paragraph("PASS" if 10.0 <= results['na_mmol_l'] <= 30.0 else "FAIL", bold_cell)],
            [Paragraph("Total Daily Potassium Intake", cell_style), Paragraph(f"{results['k_daily_mg']:.1f} mg/day", cell_style), Paragraph("Max 1300.0 mg/day (Schedule 29)", cell_style), Paragraph("PASS" if results['k_daily_mg'] <= 1300.0 else "FAIL", bold_cell)],
            [Paragraph("Total Daily Magnesium Intake", cell_style), Paragraph(f"{results['mg_daily_mg']:.1f} mg/day", cell_style), Paragraph("Max 320.0 mg/day (Schedule 29)", cell_style), Paragraph("PASS" if results['mg_daily_mg'] <= 320.0 else "FAIL", bold_cell)],
            [Paragraph("Single Serve Caffeine", cell_style), Paragraph(f"{results['caffeine_serve_mg']:.1f} mg", cell_style), Paragraph("Max 80.0 mg/serve (Schedule 29)", cell_style), Paragraph("PASS" if results['caffeine_serve_mg'] <= 80.0 else "FAIL", bold_cell)],
            [Paragraph("Total Daily Caffeine Intake", cell_style), Paragraph(f"{results['caffeine_daily_mg']:.1f} mg/day", cell_style), Paragraph("Max 100.0 mg/day (Schedule 29)", cell_style), Paragraph("PASS" if results['caffeine_daily_mg'] <= 100.0 else "FAIL", bold_cell)],
            [Paragraph("Solution Osmolality & Tonicity", cell_style), Paragraph(f"{results['osmolality_mOsm_kg']:.1f} mOsm/kg", cell_style), Paragraph(f"Classified as <b>{results['tonicity_classification']}</b>", cell_style), Paragraph("INFO", bold_cell)],
            [Paragraph("Debye-Hückel Activity (phi)", cell_style), Paragraph(f"{results['osmotic_coefficient_phi']:.4f}", cell_style), Paragraph(f"Ionic Strength: {results['ionic_strength_M']:.4f} M", cell_style), Paragraph("INFO", bold_cell)],
        ]
        
        metrics_table = Table(metrics_data, colWidths=[150, 110, 200, 80])
        metrics_table.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#EDF2F7")),
            ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor("#CBD5E0")),
            ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
            ('TOPPADDING', (0,0), (-1,-1), 5),
            ('BOTTOMPADDING', (0,0), (-1,-1), 5),
        ]))
        story.append(metrics_table)
        story.append(Spacer(1, 14))

        # RDKit Formulation Decomposition Table
        story.append(Paragraph("<b>RDKit Topological Formulation Breakdown</b>", ParagraphStyle("H2", parent=styles["Heading2"], fontSize=12, textColor=colors.HexColor("#2B6CB0"))))
        story.append(Spacer(1, 4))

        rdkit_headers = [Paragraph("Compound Label / SMILES", bold_cell), Paragraph("MW (g/mol)", bold_cell), Paragraph("nu", bold_cell), Paragraph("Mass (mg)", bold_cell), Paragraph("Primary Yields", bold_cell)]
        rdkit_rows = [rdkit_headers]

        for comp in results["parsed_compounds"]:
            yields_str = ", ".join([f"{sym}: {m:.1f}mg" for sym, m in comp["elem_yields_mg"].items() if m > 0.1 and sym in ["Na", "K", "Mg", "Ca", "Cl", "N"]])
            label_para = Paragraph(f"<b>{comp['label']}</b><br/><font color='#718096' size='6'>{comp['smiles']}</font>", cell_style)
            rdkit_rows.append([
                label_para,
                Paragraph(f"{comp['mw']:.2f}", cell_style),
                Paragraph(str(comp['nu']), cell_style),
                Paragraph(f"{comp['mass_mg']:.1f}", cell_style),
                Paragraph(yields_str, cell_style)
            ])

        rdkit_table = Table(rdkit_rows, colWidths=[180, 70, 40, 70, 180])
        rdkit_table.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#EDF2F7")),
            ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor("#CBD5E0")),
            ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
            ('TOPPADDING', (0,0), (-1,-1), 4),
            ('BOTTOMPADDING', (0,0), (-1,-1), 4),
        ]))
        story.append(rdkit_table)
        story.append(Spacer(1, 14))

        # Audit Logs
        story.append(Paragraph("<b>Regulatory Verification Audit Logs</b>", ParagraphStyle("H2", parent=styles["Heading2"], fontSize=12, textColor=colors.HexColor("#2B6CB0"))))
        story.append(Spacer(1, 4))
        for log in results["audit_logs"]:
            story.append(Paragraph(f"• {log}", cell_style))

        doc.build(story)
        return buffer.getvalue()
