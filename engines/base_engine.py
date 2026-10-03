from abc import ABC, abstractmethod
from typing import Dict, Any

class BaseComputationalEngine(ABC):
    @property
    @abstractmethod
    def engine_id(self) -> str:
        """Unique key identifying the engine."""
        pass

    @property
    @abstractmethod
    def engine_name(self) -> str:
        """Human-readable display title for UI rendering."""
        pass

    @property
    @abstractmethod
    def domain_category(self) -> str:
        """Category domain string."""
        pass

    @abstractmethod
    def render_inputs(self, st_ctx: Any) -> Dict[str, Any]:
        """Renders dynamic, domain-specific Streamlit input widgets."""
        pass

    @abstractmethod
    def execute(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Executes domain-specific mathematical algorithms."""
        pass
        
