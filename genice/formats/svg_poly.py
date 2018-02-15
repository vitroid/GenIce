# coding: utf-8
"""
"""

import numpy as np
from countrings import countrings_nx as cr
import networkx as nx
import re
import svgwrite as sw
import colorsys
from math import sin,cos

offset = np.zeros(3)

sun = np.array([1., -10., 5.])  # right, down, front
sun /= np.linalg.norm(sun)

proj = np.array([[-3**0.5/2, 1./2.],
                 [+3**0.5/2, 1./2.],
                 [0.,       -1.]])


proj = np.array([[1., -1., 0.], [1., 1., -2.], [1., 1., 1.]])
theta = 0.1
smallrot = np.array([[cos(theta),-sin(theta),0.],
                     [+sin(theta),cos(theta),0.],
                     [0.,0.,1.0]])

for i in range(3):
    proj[i] /= np.linalg.norm(proj[i])
proj = np.dot(proj,smallrot)
proj = np.linalg.inv(proj)

hue_sat = {3:(60., 1.0), 4:(180, 0.8), 5:(0, 0.5), 6:(0, 0.0), 7:(240, 0.5), 8:(300, 0.5)}


def face(center, rpos, svg):
    group = svg.g()
    n = rpos.shape[0]
    
    # normalize relative vectors
    normalized = np.zeros_like(rpos)
    for i in range(n):
        normalized[i] = rpos[i] / np.linalg.norm(rpos[i])
    #normal for each triangle
    normals = np.zeros_like(rpos)
    for i in range(n):
        normals[i] = np.cross(normalized[i-1], normalized[i])
    # central normal
    c_normal = np.sum(normals, axis=0)
    c_normal /= np.linalg.norm(c_normal)
    cosine = abs(np.dot(c_normal, sun))
    
    hue, sat = hue_sat[n]
    bri = cosine*0.5+0.5
    if sat < 0.2:
        bri *= 0.9
    if cosine > 0.8:
        sat *= (1 - (cosine-0.8)*3)

    r,g,b = colorsys.hsv_to_rgb(hue/360., sat, bri)
    rgb = "#{0:x}{1:x}{2:x}".format(int(r*15.9), int(g*15.9), int(b*15.9))

    pos = (rpos + center)[:,:2]*100 + 400
    group.add(sw.shapes.Polygon(pos, fill=rgb, stroke_linejoin="round", fill_opacity="1.0", stroke="#444", stroke_width=3))
    return group


def hook4(lattice):
    global offset
    lattice.logger.info("Hook4: SVG (polyhedral expressions).")
    graph = nx.Graph(lattice.spacegraph) #undirected
    cellmat = lattice.repcell.mat
    projected = np.dot(cellmat, proj)
    queue = []
    for ring in cr.CountRings(graph).rings_iter(8):
        deltas = np.zeros((len(ring),3))
        d2 = np.zeros(3)
        for k,i in enumerate(ring):
            d = lattice.reppositions[i] - lattice.reppositions[ring[0]]
            d -= np.floor(d+0.5)
            deltas[k] = d
            dd = lattice.reppositions[ring[k]] - lattice.reppositions[ring[k-1]]
            dd -= np.floor(dd+0.5)
            d2 += dd
        # d2 must be zeros
        if np.all(np.absolute(d2) < 1e-5):
            comofs = np.sum(deltas, axis=0) / len(ring)
            deltas -= comofs
            com = lattice.reppositions[ring[0]] + comofs + offset
            com -= np.floor(com)
            # rel to abs
            com    = np.dot(com,    projected)
            deltas = np.dot(deltas, projected)
            queue.append((com, deltas))
    # translate = -(cellmat[0,:]+cellmat[1,:]+cellmat[2,:])/2
    svg = sw.Drawing()
    for com, deltas in sorted(queue, key=lambda x: x[0][2]): #key is z of com
        svg.add(face(com, deltas, svg))
    print(svg.tostring())
    lattice.logger.info("Hook4: end.")


def argparser(arg):
    global offset
    assert re.match("^[-+0-9,.]+$", arg) is not None, "Argument must be three floating points separated by commas."
    offset = np.array([float(x) for x in arg.split(",")]).reshape(3)
        

hooks = {4:hook4}
