"""
Engine Module Initializer
Auto-registers computational engines into the central registry.
"""

from engines.registry import EngineRegistry
from engines.cutaneous_engine import CutaneousEngineAdapter
from engines.fsanz_engine import FSANZ294Engine

# Singleton instance
registry = EngineRegistry()

# Register FSANZ 2.9.4 Engine
try:
    fsanz_instance = FSANZ294Engine()
    registry.register("fsanz_294_sports_drink", fsanz_instance)
except TypeError as e:
    # Fallback if __init__ expects parameters or needs default args
    print(f"[Warning] Instantiation failed with default args: {e}")

# Register Cutaneous Adapter
try:
    cutaneous_instance = CutaneousEngineAdapter()
    registry.register("cutaneous_iata_diep", cutaneous_instance)
except Exception as e:
    print(f"[Warning] Cutaneous engine registration skipped: {e}")

# Expose both uppercase and lowercase aliases for module compatibility
REGISTRY = registry
