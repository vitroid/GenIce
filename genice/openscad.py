import numpy as np


def use(f):
    return "use <{0}>;\n".format(f)

def defvar(name, value):
    return "{0}={1};\n".format(name,value)

def scale(sc, objs):
    s = "scale([{0},{1},{2}]){{\n".format(*sc)
    for o in objs:
        s += o
    s += "} //scale\n"
    return s

def translate(sc, objs):
    s = "translate([{0},{1},{2}]){{\n".format(*sc)
    for o in objs:
        s += o
    s += "} //translate\n"
    return s

def intersection(objs):
    s = "intersection(){\n"
    for o in objs:
        s += o
    s += "} //intersection\n"
    return s

def union(objs):
    s = "union(){"
    for o in objs:
        s += o
    s += "} //union\n"
    return s

def rhomb(cell):
    origin = np.zeros(3)
    points = [x+y+z for x in (origin,cell[0]) for y in (origin,cell[1]) for z in (origin,cell[2])]
    s = "polyhedron(["
    for point in points:
        s += "[{0},{1},{2}],\n".format(*point)
    s += "],\n"
    faces = [[0,1,3,2],[0,4,5,1],[0,2,6,4],[5,4,6,7],[6,2,3,7],[3,1,5,7]]
    s += "["
    for face in faces:
        s += "[{3},{2},{1},{0}],\n".format(*face)
    s += "]);\n"
    return s;

def sphere(r=1.0):
    return "sphere(r={0});\n".format(r)

def bond(s1,s2,r=1.0):
    return "bond([{0},{1},{2}],[{3},{4},{5}],r={6});\n".format(*s1,*s2,r)



