from typing import Dict, List
from engines.base_engine import BaseComputationalEngine

class EngineRegistry:
    def __init__(self):
        self._engines: Dict[str, BaseComputationalEngine] = {}

    def register(self, engine: BaseComputationalEngine) -> None:
        self._engines[engine.engine_id] = engine

    def get(self, engine_id: str) -> BaseComputationalEngine:
        if engine_id not in self._engines:
            raise KeyError(f"Engine '{engine_id}' is not registered in Strategy Registry.")
        return self._engines[engine_id]

    def list_engines(self) -> List[Dict[str, str]]:
        return [
            {
                "id": eng.engine_id,
                "name": eng.engine_name,
                "domain": eng.domain_category
            }
            for eng in self._engines.values()
        ]

REGISTRY = EngineRegistry()
