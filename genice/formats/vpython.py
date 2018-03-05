# coding: utf-8
"""
Direct visualization with VPython.
"""

from collections import defaultdict
import colorsys
import numpy as np
import networkx as nx
import vpython as vp
from countrings import countrings_nx as cr

hue_sat = {3:(60., 1.0), 4:(180, 0.8), 5:(0, 0.5), 6:(0, 0.0), 7:(240, 0.5), 8:(300, 0.5)}

sun = np.array([1., -10., 5.])  # right, down, front
sun /= np.linalg.norm(sun)

def face(center, rpos):
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
    if np.dot(c_normal, sun) < 0.0:
        c_normal = - c_normal
        normals  = - normals
    
    hue, sat = hue_sat[n]
    bri = 1
    r,g,b = colorsys.hsv_to_rgb(hue/360., sat, bri)
    pos = rpos + center
    v_center = vp.vertex( pos=vp.vector(*center), normal=vp.vector(*c_normal), color=vp.vector(r,g,b))
    vertices = [vp.vertex( pos=vp.vector(*p), normal=vp.vector(*(normals[i])), color=vp.vector(r,g,b) ) for i,p in enumerate(pos)]
    faces = set()
    for i in range(n):
        faces.add(vp.triangle(v0=vertices[i-1], v1=vertices[i], v2=v_center,  ))
#    group.add(sw.shapes.Polygon(pos, fill=rgb, stroke_linejoin="round", fill_opacity="1.0", stroke="#444", stroke_width=3))
    return faces


def hook4(lattice):
    lattice.logger.info("Hook4: VPython (polyhedral expressions).")
    graph = nx.Graph(lattice.spacegraph) #undirected
    cellmat = lattice.repcell.mat
    lattice.vpobjects = defaultdict(set)
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
            com = lattice.reppositions[ring[0]] + comofs
            com -= np.floor(com)
            # rel to abs
            com    = np.dot(com-0.5,    cellmat)
            deltas = np.dot(deltas, cellmat)
            ringsize = '{0}'.format(len(ring))
            lattice.vpobjects[ringsize] |= face(com, deltas)
    lattice.logger.info("Hook4: end.")








def draw(lattice):
    lattice.logger.info("Hook6: Display water molecules with VPython.")
    lattice.logger.info("  Total number of atoms: {0}".format(len(lattice.atoms)))
    cellmat = lattice.repcell.mat
    offset = (cellmat[0] + cellmat[1] + cellmat[2]) / 2
    # prepare the reverse dict
    waters = defaultdict(dict)
    for atom in lattice.atoms:
        resno, resname, atomname, position, order = atom
        if "O" in atomname:
            waters[order]["O"] = position - offset
        elif "H" in atomname:
            if "H0" not in waters[order]:
                waters[order]["H0"] = position - offset
            else:
                waters[order]["H1"] = position - offset
    for order, water in waters.items():
        O = water["O"]        
        H0 = water["H0"]        
        H1 = water["H1"]        
        lattice.vpobjects['w'].add(vp.simple_sphere(radius=0.03, pos=vp.vector(*O), color=vp.vector(1,0,0)))
        lattice.vpobjects['w'].add(vp.simple_sphere(radius=0.02, pos=vp.vector(*H0), color=vp.vector(0,1,1)))
        lattice.vpobjects['w'].add(vp.simple_sphere(radius=0.02, pos=vp.vector(*H1), color=vp.vector(0,1,1)))
        lattice.vpobjects['w'].add(vp.cylinder(radius=0.015, pos=vp.vector(*O), axis=vp.vector(*(H0-O))))
        lattice.vpobjects['w'].add(vp.cylinder(radius=0.015, pos=vp.vector(*O), axis=vp.vector(*(H1-O))))
        lattice.vpobjects['l'].add(vp.label(pos=vp.vector(*O), xoffset=30, text="{0}".format(order), visible=False))
    for i,j in lattice.spacegraph.edges(data=False):
        if i in waters and j in waters:  # edge may connect to the dopant
            O = waters[j]["O"]
            H0 = waters[i]["H0"]
            H1 = waters[i]["H1"]
            d0 = H0 - O
            d1 = H1 - O
            rr0 = np.dot(d0,d0)
            rr1 = np.dot(d1,d1)
            if rr0 < rr1 and rr0 < 0.27**2:
                lattice.vpobjects['a'].add(vp.arrow(shaftwidth=0.015, pos=vp.vector(*H0), axis=vp.vector(*(O-H0)), color=vp.vector(1,1,0)))
            if rr1 < rr0 and rr1 < 0.245**2:
                lattice.vpobjects['a'].add(vp.arrow(shaftwidth=0.015, pos=vp.vector(*H1), axis=vp.vector(*(O-H1)), color=vp.vector(1,1,0)))
    lattice.logger.info("  Tips: use keys to draw/hide layers. [3 4 5 6 7 8 a w l]")
    lattice.logger.info("  Tips: Type ctrl-C twice at the terminal to stop.")
    lattice.logger.info("Hook6: end.")


def hook6(lattice):
    draw(lattice)
    # Toggle visibility
    visible = dict()
    for L in lattice.vpobjects:
        visible[L] = False
    visible['l'] = True
    while True:
        ev = vp.scene.waitfor("keydown")
        if ev.key in lattice.vpobjects:
            for L in lattice.vpobjects[ev.key]:
                L.visible = visible[ev.key]
            visible[ev.key] = not visible[ev.key]



hooks = {6:hook6, 4:hook4}
