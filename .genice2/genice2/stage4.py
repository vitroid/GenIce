from dataclasses import dataclass
import numpy as np
import networkx as nx
from genice2.decorators import banner, timeit
import genice_core
from genice2 import ConfigurationError


@dataclass
class Stage4Output:
    """Stage3の出力データ"""

    digraph: nx.DiGraph  # 有向グラフ


class Stage4:
    """有向氷グラフの生成を行うステージ"""

    def __init__(
        self, graph: nx.Graph, reppositions: np.ndarray, fixedEdges: nx.DiGraph
    ):
        self.graph = graph
        self.reppositions = reppositions
        self.fixedEdges = fixedEdges

    @timeit
    @banner
    def execute(self, depol: str = "none") -> Stage4Output:
        """Makes a directed graph."""
        iter = 0 if depol == "none" else 1000
        digraph = genice_core.ice_graph(
            self.graph.to_undirected(),
            vertexPositions=self.reppositions,
            isPeriodicBoundary=True,
            dipoleOptimizationCycles=iter,
            fixedEdges=self.fixedEdges,
        )
        if not digraph:
            raise ConfigurationError("Failed to generate a directed graph.")
        return Stage4Output(digraph=digraph)
