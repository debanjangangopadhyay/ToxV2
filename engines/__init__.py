from engines.registry import REGISTRY
from engines.cutaneous_engine import CutaneousBioactivationEngine
from engines.fsanz_engine import FSANZ294Engine

REGISTRY.register(CutaneousBioactivationEngine())
REGISTRY.register(FSANZ294Engine())
