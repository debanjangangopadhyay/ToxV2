"""
Engine Module Initializer
Auto-registers computational engines into the central registry singleton
and exposes registered_engines as a dictionary map for app.py UI routing.
"""

from engines.registry import REGISTRY

# 1. Register Cutaneous Engine Adapter
try:
    from engines.cutaneous_engine import CutaneousEngineAdapter
    REGISTRY.register(CutaneousEngineAdapter())
except Exception as e:
    print(f"[Warning] Cutaneous Engine registration error: {e}")

# 2. Register FSANZ 2.9.4 Food Science Engine
try:
    from engines.fsanz_engine import FSANZ294Engine
    REGISTRY.register(FSANZ294Engine())
except Exception as e:
    print(f"[Warning] FSANZ Engine registration error: {e}")


def _build_registered_engines_map():
    """
    Safely converts registered engine instances from REGISTRY into a dictionary 
    mapping engine_id -> metadata dictionary expected by app.py.
    """
    engine_map = {}
    
    # Check if REGISTRY has get_all() or get_all_metadata()
    if hasattr(REGISTRY, "get_all"):
        raw_engines = REGISTRY.get_all()
    elif hasattr(REGISTRY, "engines"):
        raw_engines = REGISTRY.engines
    else:
        raw_engines = []

    # Handle dictionary vs list storage within REGISTRY
    if isinstance(raw_engines, dict):
        for eid, item in raw_engines.items():
            if hasattr(item, "get_metadata"):
                engine_map[eid] = item.get_metadata()
            elif isinstance(item, dict):
                engine_map[eid] = item
            else:
                engine_map[eid] = {"id": eid, "instance": item}
    elif isinstance(raw_engines, list):
        for item in raw_engines:
            if hasattr(item, "get_metadata"):
                meta = item.get_metadata()
                eid = meta.get("id", str(item))
                engine_map[eid] = meta
            elif isinstance(item, dict):
                eid = item.get("id", str(item))
                engine_map[eid] = item

    return engine_map


# Expose dictionary mapping for app.py (Line 117 loop)
registered_engines = _build_registered_engines_map()
