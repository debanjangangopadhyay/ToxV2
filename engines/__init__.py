"""
Engine Module Initializer
Auto-registers computational engines into the central registry.
"""

from engines.registry import EngineRegistry
from engines.cutaneous_engine import CutaneousEngineAdapter
from engines.fsanz_engine import FSANZ294Engine

# Central Registry Singleton
registry = EngineRegistry()

# Register FSANZ 2.9.4 Engine
try:
    fsanz_instance = FSANZ294Engine()
    registry.register("fsanz_294_sports_drink", fsanz_instance)
except Exception as e:
    print(f"[Warning] FSANZ Engine registration error: {e}")

# Register Cutaneous Bioactivation Adapter
try:
    cutaneous_instance = CutaneousEngineAdapter()
    registry.register("cutaneous_iata_diep", cutaneous_instance)
except Exception as e:
    print(f"[Warning] Cutaneous Engine registration error: {e}")

# Expose both uppercase and lowercase aliases to eliminate import mismatch errors
REGISTRY = registry
