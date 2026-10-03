"""
Engine Module Initializer
Auto-registers computational engines into the central registry.
"""

from engines.registry import EngineRegistry
from engines.cutaneous_engine import CutaneousEngineAdapter
from engines.fsanz_engine import FSANZ294Engine

# Initialize singleton instance
registry = EngineRegistry()

# Register computational strategy handlers
registry.register("fsanz_294_sports_drink", FSANZ294Engine())
registry.register("cutaneous_iata_diep", CutaneousEngineAdapter())
