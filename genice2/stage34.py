from dataclasses import dataclass
import numpy as np
import networkx as nx
from genice2.decorators import banner, timeit
import genice_core


@dataclass
class Stage34EOutput:
    """Stage34Eの出力データ"""

    digraph: nx.DiGraph  # 有向グラフ


class Stage34E:
    """有向氷グラフの生成を行うステージ"""

    def __init__(
        self, graph: nx.Graph, reppositions: np.ndarray, fixedEdges: nx.DiGraph
    ):
        self.graph = graph
        self.reppositions = reppositions
        self.fixedEdges = fixedEdges

    @timeit
    @banner
    def execute(self, depol: str = "none") -> Stage34EOutput:
        """Makes a directed graph."""
        iter = 0 if depol == "none" else 1000
        digraph = genice_core.ice_graph(
            self.graph.to_undirected(),
            vertexPositions=self.reppositions,
            isPeriodicBoundary=True,
            dipoleOptimizationCycles=iter,
            fixedEdges=self.fixedEdges,
        )
        return Stage34EOutput(digraph=digraph)
