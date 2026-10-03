"""
Engine Module Initializer
Auto-registers computational engines into the central registry singleton.
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
    
