from dataclasses import dataclass
import numpy as np
from math import sin, cos, radians


@dataclass
class Molecule:
    """
    Base class of a molecule
    """

    sites: np.ndarray
    labels: list[str]
    name: str
    is_water: bool = False

    def __repr__(self) -> str:
        return (
            f"Molecule(name={self.name!r}, "
            f"n_sites={len(self.sites)}, "
            f"labels={self.labels}, "
            f"is_water={self.is_water})"
        )
