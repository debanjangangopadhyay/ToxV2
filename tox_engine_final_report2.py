"""
Integrated Probabilistic Toxicology & Evidence Engine
Version 2.0 — research-use computational screening

Scope:
- Cosmetic / chemical pre-screening
- Probabilistic MoS and uncertainty propagation
- Dependency-aware evidence fusion
- Mechanism-aware read-across
- Cutaneous transformation hypotheses
- Applicability-domain diagnostics
- Experimental test-priority / value-of-information framework

Important:
This software is NOT a regulatory decision engine, does NOT establish human
safety, and does NOT replace validated toxicological testing.
"""

from __future__ import annotations
import hashlib
import json
import os
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from io import BytesIO
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, Draw, rdFingerprintGenerator
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import HRFlowable, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle, Image as PlatypusImage

RDLogger.DisableLog("rdApp.error")
RDLogger.DisableLog("rdApp.warning")

ENGINE_VERSION = "2.0.0"
DEFAULT_SEED = 42
DEFAULT_MC_ITERATIONS = 10000
MOS_SCREENING_THRESHOLD = 100.0

# These are configuration priors, not empirical estimates. They are deliberately
# conservative placeholders until fitted on a locked reference dataset.
DEFAULT_PRIORS = {
    "mos_epistemic": 0.10,
    "structural_alert": 0.50,
    "read_across": 0.50,
    "regulatory_flag": 0.95,
}

@dataclass
class EvidenceMass:
    support: float
    concern: float
    uncertainty: float
    source: str = ""

    def __post_init__(self):
        vals = np.asarray([self.support, self.concern, self.uncertainty], dtype=float)
        if not np.all(np.isfinite(vals)) or np.any(vals < 0):
            raise ValueError("Evidence masses must be finite and non-negative.")
        total = float(vals.sum())
        if total <= 0:
            raise ValueError("Evidence masses must have positive total mass.")
        self.support, self.concern, self.uncertainty = (vals / total).tolist()

def fuse_independent_mass(m1: EvidenceMass, m2: EvidenceMass) -> EvidenceMass:
    """Dempster combination for genuinely independent evidence sources."""
    conflict = m1.support*m2.concern + m1.concern*m2.support
    if conflict >= 0.999999:
        return EvidenceMass(0.0, 0.0, 1.0, "conflict")
    denom = 1.0 - conflict
    support = (m1.support*m2.support + m1.support*m2.uncertainty +
               m1.uncertainty*m2.support) / denom
    concern = (m1.concern*m2.concern + m1.concern*m2.uncertainty +
               m1.uncertainty*m2.concern) / denom
    uncertainty = (m1.uncertainty*m2.uncertainty) / denom
    return EvidenceMass(support, concern, uncertainty, f"{m1.source}+{m2.source}")

def dependency_adjusted_fusion(masses: Sequence[EvidenceMass],
                               dependency_penalty: float = 0.35) -> EvidenceMass:
    """
    Conservative sequential fusion. Correlated streams receive reduced weight
    before fusion; this prevents structural alert/read-across double counting.
    """
    if not masses:
        return EvidenceMass(1/3, 1/3, 1/3, "neutral")
    penalty = float(np.clip(dependency_penalty, 0.0, 1.0))
    first = masses[0]
    for m in masses[1:]:
        # Shrink the second stream toward ignorance according to estimated dependence.
        shrink = penalty
        adjusted = EvidenceMass(
            support=m.support*(1-shrink) + shrink/3,
            concern=m.concern*(1-shrink) + shrink/3,
            uncertainty=m.uncertainty*(1-shrink) + shrink/3,
            source=m.source
        )
        first = fuse_independent_mass(first, adjusted)
    return first

def _validate_positive(name: str, value: float) -> float:
    x = float(value)
    if not np.isfinite(x) or x <= 0:
        raise ValueError(f"{name} must be finite and > 0.")
    return x

def _validate_fraction(name: str, value: float, upper: float = 100.0) -> float:
    x = float(value)
    if not np.isfinite(x) or x < 0 or x > upper:
        raise ValueError(f"{name} must be between 0 and {upper}.")
    return x

def _seeded_rng(seed: int = DEFAULT_SEED) -> np.random.Generator:
    return np.random.default_rng(int(seed))

def run_structural_filter(mol: Chem.Mol) -> List[str]:
    if mol is None:
        return ["Invalid molecular entity"]
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    return [m.GetDescription() for m in catalog.GetMatches(mol)]

CUTANEOUS_REACTION_RULES = {
    "Ester hydrolysis": "[CX3:1](=[OX1:2])[OX2][C:4]>>[CX3:1](=[OX1:2])[OH].[C:4][OH]",
    "Catechol oxidation hypothesis": "[c:1]1[c:2]([OH])[c:3]([OH])[c:4][c:5][c:6]1>>[c:1]1[c:2](=O)[c:3](=O)[c:4][c:5][c:6]1",
    "Allylic alcohol oxidation hypothesis": "[C:1]=[C:2]-[CH2:3][OH]>>[C:1]=[C:2]-[CH:3]=O",
    "Alkene epoxidation hypothesis": "[C:1]=[C:2]>>[C:1]1[O][C:2]1",
    "Aromatic amine N-hydroxylation hypothesis": "[c:1][NH2:2]>>[c:1][NH:2][OH]",
}

def simulate_skin_metabolism(parent_mol: Chem.Mol) -> List[Dict[str, Any]]:
    if parent_mol is None:
        return []
    seen = {Chem.MolToSmiles(parent_mol)}
    records = []
    for pathway, smarts in CUTANEOUS_REACTION_RULES.items():
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((parent_mol,))
            for tup in products:
                for met in tup:
                    try:
                        Chem.SanitizeMol(met)
                        smi = Chem.MolToSmiles(met)
                        if smi in seen or len(smi) < 3:
                            continue
                        seen.add(smi)
                        alerts = run_structural_filter(met)
                        records.append({
                            "pathway": pathway,
                            "metabolite_smiles": smi,
                            "secondary_alerts": alerts,
                            "reactive_alert_hypothesis": bool(alerts),
                            "interpretation": (
                                "Transformation generated a secondary structural alert; "
                                "hypothesis requires experimental confirmation."
                                if alerts else
                                "No monitored secondary alert generated."
                            )
                        })
                    except Exception:
                        continue
        except Exception:
            continue
    return records

def calculate_deterministic_sed(pod: float, concentration_pct: float,
                                daily_amount_g: float, retention: float,
                                dermal_absorption_pct: float, body_weight_kg: float) -> Dict[str, float]:
    pod = _validate_positive("POD", pod)
    concentration_pct = _validate_fraction("concentration_pct", concentration_pct)
    daily_amount_g = _validate_positive("daily_amount_g", daily_amount_g)
    retention = _validate_fraction("retention", retention, 1.0)
    dermal_absorption_pct = _validate_fraction("dermal_absorption_pct", dermal_absorption_pct)
    body_weight_kg = _validate_positive("body_weight_kg", body_weight_kg)
    applied_mg = daily_amount_g * 1000.0
    amount_mg = applied_mg * concentration_pct/100.0 * retention
    sed = amount_mg * dermal_absorption_pct/100.0 / body_weight_kg
    mos = float("inf") if sed == 0 else pod/sed
    return {"sed_mg_kg_day": sed, "mos": mos}

def calculate_probabilistic_mos(
    pod: float,
    concentration_pct: float,
    dermal_absorption_pct: float,
    daily_amount_g: float = 1.54,
    retention: float = 1.0,
    body_weight_kg: float = 60.0,
    n: int = DEFAULT_MC_ITERATIONS,
    seed: int = DEFAULT_SEED,
    pod_cv: float = 0.10,
    absorption_cv: float = 0.25,
) -> Dict[str, Any]:
    """
    Propagates uncertainty using lognormal positive distributions.
    The absorption CV is an explicit model assumption, not a fabricated
    empirical calibration. It must be replaced by fitted data when available.
    """
    pod = _validate_positive("POD", pod)
    concentration_pct = _validate_fraction("concentration_pct", concentration_pct)
    dermal_absorption_pct = _validate_fraction("dermal_absorption_pct", dermal_absorption_pct)
    daily_amount_g = _validate_positive("daily_amount_g", daily_amount_g)
    retention = _validate_fraction("retention", retention, 1.0)
    body_weight_kg = _validate_positive("body_weight_kg", body_weight_kg)
    n = int(n)
    if n < 1000:
        raise ValueError("n must be >= 1000 for stable percentile estimation.")
    if pod_cv < 0 or absorption_cv < 0:
        raise ValueError("CV values must be non-negative.")

    rng = _seeded_rng(seed)

    def lognormal_from_mean_cv(mean, cv, size):
        if cv == 0:
            return np.full(size, mean, dtype=float)
        sigma2 = np.log1p(cv**2)
        sigma = np.sqrt(sigma2)
        mu = np.log(mean) - 0.5*sigma2
        return rng.lognormal(mu, sigma, size)

    bw = np.clip(lognormal_from_mean_cv(body_weight_kg, 0.17, n), 35.0, 150.0)
    da = np.clip(lognormal_from_mean_cv(dermal_absorption_pct, absorption_cv, n), 0.01, 100.0)
    amount = np.clip(lognormal_from_mean_cv(daily_amount_g, 0.15, n), 0.05, 100.0)
    pod_dist = lognormal_from_mean_cv(pod, pod_cv, n)

    sed = (amount*1000.0*(concentration_pct/100.0)*retention*(da/100.0))/bw
    mos = np.divide(pod_dist, sed, out=np.full_like(sed, np.inf), where=sed > 0)

    finite = mos[np.isfinite(mos)]
    if finite.size == 0:
        raise ValueError("Monte Carlo produced no finite MoS values.")

    return {
        "iterations": n,
        "seed": int(seed),
        "median_mos": float(np.median(finite)),
        "p05_mos": float(np.percentile(finite, 5)),
        "p95_mos": float(np.percentile(finite, 95)),
        "p_mos_below_100": float(np.mean(finite < MOS_SCREENING_THRESHOLD)),
        "mean_mos": float(np.mean(finite)),
        "absorption_cv": float(absorption_cv),
        "pod_cv": float(pod_cv),
        "interpretation": (
            "P(MoS < 100) is the simulated probability that the modeled MoS "
            "falls below the selected screening threshold under the stated "
            "uncertainty assumptions. It is NOT a probability of toxicity or harm."
        )
    }

REFERENCE_ANALOGS = [
    {"name":"Niacinamide","smiles":"NC(=O)c1cccnc1","category":"skin barrier"},
    {"name":"Allantoin","smiles":"O=C1NC(=O)NC1NC(=O)N","category":"soothing"},
    {"name":"Phenethyl alcohol","smiles":"OCCc1ccccc1","category":"preservative/booster"},
    {"name":"Ferulic acid","smiles":"COc1cc(/C=C/C(=O)O)ccc1O","category":"antioxidant"},
    {"name":"Panthenol","smiles":"OCC(C)(C)C(O)C(=O)NCCCCO","category":"humectant"},
]

def _descriptor_vector(mol):
    return np.array([Descriptors.MolLogP(mol), Descriptors.MolWt(mol),
                     Descriptors.TPSA(mol), Descriptors.NumHDonors(mol),
                     Descriptors.NumHAcceptors(mol)], dtype=float)

def execute_multidimensional_read_across(target_mol: Chem.Mol, top_k: int = 3) -> List[Dict[str, Any]]:
    if target_mol is None:
        return []
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    tfp = gen.GetFingerprint(target_mol)
    tv = _descriptor_vector(target_mol)
    out = []
    for ref in REFERENCE_ANALOGS:
        rm = Chem.MolFromSmiles(ref["smiles"])
        if rm is None:
            continue
        sim = DataStructs.TanimotoSimilarity(tfp, gen.GetFingerprint(rm))
        rv = _descriptor_vector(rm)
        scale = np.array([5.0, 300.0, 100.0, 5.0, 10.0])
        desc_dist = float(np.linalg.norm((tv-rv)/scale)/np.sqrt(len(scale)))
        mech = max(0.0, 1.0-desc_dist)
        composite = 0.65*sim + 0.35*mech
        out.append({
            "name": ref["name"], "category": ref["category"],
            "tanimoto": float(sim), "descriptor_similarity": float(mech),
            "mechanistic_similarity_proxy": float(mech),
            "read_across_confidence": float(composite),
            "warning": "Analogue similarity is not toxicological equivalence; endpoint data must be verified."
        })
    return sorted(out, key=lambda x: x["read_across_confidence"], reverse=True)[:top_k]

def applicability_domain(target_mol: Chem.Mol, analogs: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    if target_mol is None:
        return {"status":"invalid"}
    best = max((a["tanimoto"] for a in analogs), default=0.0)
    if best >= 0.70:
        status = "inside"
    elif best >= 0.50:
        status = "borderline"
    else:
        status = "outside"
    return {"status": status, "best_tanimoto": best,
            "rule": "heuristic structural domain; requires dataset-specific validation"}

def value_of_information_prioritization(result: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Transparent experimental prioritization heuristic.
    This is a prioritization aid, not an experimentally validated VOI estimate.
    """
    p = result["mos"]["p_mos_below_100"]
    u = result["fusion"]["epistemic_uncertainty"]
    ad = result["applicability_domain"]["status"]
    ad_penalty = {"inside":0.2, "borderline":0.6, "outside":1.0}[ad]
    candidates = [
        ("dermal_permeation", 0.45, 0.50),
        ("repeated_exposure_or_POD_refinement", 0.60, 0.60),
        ("bioactivation_metabolite_confirmation", 0.55, 0.45),
        ("endpoint_specific_in_vitro_assay", 0.75, 0.80),
    ]
    rows = []
    for test, information_gain, cost_proxy in candidates:
        priority = (0.45*p + 0.35*u + 0.20*ad_penalty) * information_gain / cost_proxy
        rows.append({"test":test, "priority_index":float(priority),
                     "interpretation":"Transparent research-prioritization index; not a validated clinical/regulatory score."})
    return sorted(rows, key=lambda x: x["priority_index"], reverse=True)

def execute_provenance_audit(mol: Chem.Mol, schema_path: str = "regulatory_rules.json") -> Dict[str, Any]:
    result = {"available":False, "violations":[], "warnings":[], "source":schema_path}
    if not os.path.exists(schema_path):
        return result
    try:
        with open(schema_path, "r", encoding="utf-8") as fh:
            rules = json.load(fh).get("regulatory_rules", [])
        for rule in rules:
            smarts = rule.get("smarts")
            if not smarts:
                continue
            pat = Chem.MolFromSmarts(smarts)
            if pat and mol.HasSubstructMatch(pat):
                flag = {k:rule.get(k) for k in ("rule_id","target_class","status","legal_instrument")}
                if "prohibited" in str(rule.get("status","")).lower():
                    result["violations"].append(flag)
                else:
                    result["warnings"].append(flag)
        result["available"] = True
    except (OSError, ValueError, json.JSONDecodeError):
        result["available"] = False
    return result

def _assessment_id(smiles: str) -> str:
    return "TOX2-" + hashlib.sha256(smiles.encode("utf-8")).hexdigest()[:12].upper()

def execute_full_compound_audit(
    smiles: str,
    pod: float,
    concentration_pct: float,
    dermal_absorption_pct: float,
    daily_amount_g: float = 1.54,
    retention: float = 1.0,
    n: int = DEFAULT_MC_ITERATIONS,
    seed: int = DEFAULT_SEED,
    ablation_config: Optional[Dict[str,bool]] = None,
) -> Dict[str,Any]:
    ablation_config = ablation_config or {"use_ds":True,"use_metabolism":True,"use_read_across":True}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"valid":False, "error":"Invalid SMILES string."}
    Chem.SanitizeMol(mol)
    clean = Chem.MolToSmiles(mol)
    alerts = run_structural_filter(mol)
    reg = execute_provenance_audit(mol)
    met = simulate_skin_metabolism(mol) if ablation_config.get("use_metabolism",True) else []
    mos = calculate_probabilistic_mos(
        pod=pod,
        concentration_pct=concentration_pct,
        dermal_absorption_pct=dermal_absorption_pct,
        daily_amount_g=daily_amount_g,
        retention=retention,
        n=n,
        seed=seed,
    )
    analogs = execute_multidimensional_read_across(mol) if ablation_config.get("use_read_across",True) else []
    ad = applicability_domain(mol, analogs)

    masses = []
    p = mos["p_mos_below_100"]
    masses.append(EvidenceMass(1-p, p, DEFAULT_PRIORS["mos_epistemic"], "probabilistic_mos"))
    if alerts:
        masses.append(EvidenceMass(0, DEFAULT_PRIORS["structural_alert"], 1-DEFAULT_PRIORS["structural_alert"], "structural_alert"))
    if reg["violations"]:
        masses.append(EvidenceMass(0, DEFAULT_PRIORS["regulatory_flag"], 1-DEFAULT_PRIORS["regulatory_flag"], "regulatory_screen"))
    if analogs:
        a = analogs[0]["read_across_confidence"]
        masses.append(EvidenceMass(1-a, 0, a, "read_across"))

    if ablation_config.get("use_ds",True):
        fused = dependency_adjusted_fusion(masses, dependency_penalty=0.35)
    else:
        fused = EvidenceMass(1/3,1/3,1/3,"ablation-neutral")

    fusion = {"belief_support":fused.support,
              "belief_concern":fused.concern,
              "epistemic_uncertainty":fused.uncertainty,
              "caveat":"Evidence streams are dependency-adjusted; this is not a calibrated probability of safety or toxicity."}

    result = {
        "valid":True, "engine_version":ENGINE_VERSION, "assessment_id":_assessment_id(clean),
        "mol": mol,
        "timestamp_utc":datetime.now(timezone.utc).isoformat(), "smiles":clean,
        "logp":float(Descriptors.MolLogP(mol)), "alerts":alerts, "reg":reg,
        "metabolism":met, "mos":mos, "analogs":analogs,
        "applicability_domain":ad, "fusion":fusion,
    }
    result["test_priorities"] = value_of_information_prioritization(result)
    result["limitations"] = [
        "Structural alerts are hazard hypotheses, not endpoint predictions.",
        "Read-across scores are similarity measures, not proof of toxicological equivalence.",
        "Transformation rules are mechanistic hypotheses requiring experimental confirmation.",
        "Monte Carlo uncertainty distributions are assumptions unless fitted to empirical data.",
        "Regulatory screening is limited to the supplied ruleset and is not a legal compliance determination.",
        "The integrated fusion is not externally calibrated and must not be interpreted as a probability of safety or toxicity.",
        "Applicability-domain thresholds are heuristic until validated on a reference dataset."
    ]
    return result

def generate_enterprise_pdf(data: Dict[str,Any]) -> BytesIO:
    buf = BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=letter, rightMargin=36,leftMargin=36,topMargin=36,bottomMargin=36)
    styles = getSampleStyleSheet()
    h1 = ParagraphStyle("H1", parent=styles["Normal"], fontName="Helvetica-Bold",fontSize=12,spaceBefore=10,spaceAfter=6)
    body = ParagraphStyle("Body", parent=styles["Normal"], fontSize=8.5, leading=11)
    story = [Paragraph("PROBABILISTIC COMPUTATIONAL TOXICOLOGY ASSESSMENT", ParagraphStyle("T",parent=styles["Normal"],fontName="Helvetica-Bold",fontSize=18)),
             Paragraph("Research-use pre-screening | Engine 2.0", body),
             HRFlowable(width="100%", thickness=1), Spacer(1,10)]
    if data.get("valid"):
        img = BytesIO(); Draw.MolToImage(data["mol"],size=(180,180)).save(img,format="PNG"); img.seek(0)
        story += [PlatypusImage(img,width=120,height=120),
                  Paragraph(f"<b>Assessment:</b> {data['assessment_id']}<br/><b>MoS median:</b> {data['mos']['median_mos']:.2f}<br/><b>MoS 5th percentile:</b> {data['mos']['p05_mos']:.2f}<br/><b>P(MoS &lt; 100):</b> {data['mos']['p_mos_below_100']:.3f}<br/><b>AD:</b> {data['applicability_domain']['status']}",body),
                  Spacer(1,8), Paragraph("Evidence and uncertainty",h1),
                  Paragraph(f"Belief support={data['fusion']['belief_support']:.3f}; concern={data['fusion']['belief_concern']:.3f}; epistemic uncertainty={data['fusion']['epistemic_uncertainty']:.3f}. These are evidence masses, not probabilities of safety/toxicity.",body),
                  Paragraph("Cutaneous transformation hypotheses",h1)]
        for m in data["metabolism"][:6]:
            story.append(Paragraph(f"{m['pathway']}: {m['metabolite_smiles']} — {m['interpretation']}",body))
        story += [Paragraph("Testing-priority candidates",h1)]
        for t in data["test_priorities"]:
            story.append(Paragraph(f"{t['test']}: index {t['priority_index']:.3f}",body))
        story += [Spacer(1,8), Paragraph("Limitations",h1)]
        for x in data["limitations"]: story.append(Paragraph("• "+x,body))
    doc.build(story); buf.seek(0); return buf
