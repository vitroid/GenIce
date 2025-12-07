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
