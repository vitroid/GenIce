# coding: utf-8

import itertools as it
import sys
from logging import getLogger

import numpy as np
import networkx as nx
import plotly.graph_objects as go
from genice3.genice import GenIce3


desc = {
    "ref": {"Codes": "https://github.com/vitroid/Yaplot"},
    "brief": "Plot the topology with Matplotlib.",
    "usage": """
Usage: genice2 icename -f matplotlib[options]

options:
    type=x  Set the plot type
            full or digraph: fully arranged.
            fixed: fixed edges and undecided edges.
            frame or graph: undirected graph
    ax      matplotlib axes object
""",
}


def draw_graph(g: nx.Graph, pos: dict, fixed: nx.DiGraph = None, dopant: list = []):
    # draw the graph with edges and labels in 3D
    pos = np.array([pos[i] for i in sorted(g.nodes())])

    # 通常の辺のリスト
    normal_edges = [
        (i, j)
        for i, j in g.edges()
        if fixed is None or not fixed.has_edge(i, j) or not fixed.has_edge(j, i)
    ]
    normal_edge_vectors = []
    for i, j in normal_edges:
        d = pos[j] - pos[i]
        d -= np.floor(d + 0.5)
        normal_edge_vectors.append([pos[i], pos[i] + d])

    # 固定された辺のリスト
    fixed_edges = [edge for edge in fixed.edges()] if fixed is not None else []
    fixed_edge_vectors = []
    for i, j in fixed_edges:
        d = pos[j] - pos[i]
        d -= np.floor(d + 0.5)
        fixed_edge_vectors.append([pos[i], pos[i] + d])

    normal_edge_vectors = np.array(normal_edge_vectors)
    fixed_edge_vectors = np.array(fixed_edge_vectors)

    normal_nodes = [i for i in g.nodes() if i not in dopant]
    dopant_nodes = [i for i in dopant]

    normal_pos = pos[normal_nodes]
    dopant_pos = pos[dopant_nodes]

    return go.Figure(
        data=[
            # ノードの表示
            go.Scatter3d(
                x=normal_pos[:, 0],
                y=normal_pos[:, 1],
                z=normal_pos[:, 2],
                mode="markers+text",
                marker=dict(size=2, color="blue"),
                text=[str(i) for i in normal_nodes],
                textposition="top center",
            ),
            go.Scatter3d(
                x=dopant_pos[:, 0],
                y=dopant_pos[:, 1],
                z=dopant_pos[:, 2],
                mode="markers+text",
                marker=dict(size=4, color="red"),
                text=[str(i) for i in dopant_nodes],
                textposition="top center",
            ),
            # 通常の辺の表示
            *(
                [
                    go.Scatter3d(
                        x=edge_vector[:, 0],
                        y=edge_vector[:, 1],
                        z=edge_vector[:, 2],
                        mode="lines",
                        line=dict(color="gray", width=2),
                        hoverinfo="none",
                    )
                    for edge_vector in normal_edge_vectors
                ]
                if len(normal_edge_vectors) > 0
                else []
            ),
            # 固定された辺（矢印）の表示
            *(
                [
                    go.Scatter3d(
                        x=edge_vector[:, 0],
                        y=edge_vector[:, 1],
                        z=edge_vector[:, 2],
                        mode="lines",
                        line=dict(color="green", width=3),
                        hoverinfo="none",
                    )
                    for edge_vector in fixed_edge_vectors
                ]
                if len(fixed_edge_vectors) > 0
                else []
            ),
            # 矢印の先端（コーン）の表示
            *(
                [
                    go.Cone(
                        x=fixed_edge_vectors[:, 1, 0],
                        y=fixed_edge_vectors[:, 1, 1],
                        z=fixed_edge_vectors[:, 1, 2],
                        u=fixed_edge_vectors[:, 1, 0] - fixed_edge_vectors[:, 0, 0],
                        v=fixed_edge_vectors[:, 1, 1] - fixed_edge_vectors[:, 0, 1],
                        w=fixed_edge_vectors[:, 1, 2] - fixed_edge_vectors[:, 0, 2],
                        colorscale=[[0, "green"], [1, "green"]],
                        showscale=False,
                        anchor="tip",
                    )
                ]
                if len(fixed_edge_vectors) > 0
                else []
            ),
        ]
    )


def figure(genice: GenIce3, **options):
    "Draw the topology in the cell with Plotly."
    logger = getLogger()

    graph_type = options.get("type", "full")
    if graph_type == "full" or graph_type == "digraph":
        fixed_edges = genice.digraph
    elif graph_type == "fixed":
        fixed_edges = genice.fixedEdges
    elif graph_type == "frame" or graph_type == "graph":
        fixed_edges = nx.DiGraph()
    else:
        raise ValueError(f"Invalid graph type: {graph_type}")

    # plotのしかたはgenice-coreにあったぞ。
    return draw_graph(
        genice.graph,
        genice.lattice_sites,
        fixed_edges,
        genice.substitutional_ions().keys(),
    )
