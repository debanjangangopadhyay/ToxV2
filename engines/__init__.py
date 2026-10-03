"""
Engine Module Initializer
Auto-registers computational engines into the central registry.
"""

from engines.registry import REGISTRY

# Register FSANZ 2.9.4 Engine
try:
    from engines.fsanz_engine import FSANZ294Engine
    REGISTRY.register(FSANZ294Engine())
except Exception as e:
    print(f"[Warning] FSANZ Engine registration error: {e}")

# Register Cutaneous Bioactivation Adapter
try:
    from engines.cutaneous_engine import CutaneousEngineAdapter
    REGISTRY.register(CutaneousEngineAdapter())
except Exception as e:
    print(f"[Warning] Cutaneous Engine registration error: {e}")
    
