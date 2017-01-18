#!/usr/bin/env python2

import read_cif
import itertools as it

def shortest_distance(atoms):
    dmin = 1e99
    for a1,a2 in it.combinations(atoms,2):
        name1,x1,y1,z1 = a1
        name2,x2,y2,z2 = a2
        d = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
        if d < dmin:
            dmin = d
    return dmin**0.5


# python format for GenIce.
def write_py(atoms, box, f, matchfunc=None):
    filtered = []
    if matchfunc is not None:
        for a in atoms:
            if matchfunc(a[0]):
                filtered.append(a)
    else:
        filtered = atoms
    dmin = shortest_distance(filtered)
    scale = 2.76 / dmin

    s = ""
    if (len(box) == 6):
        logger.error("Non-orthogonal cell is not supported.")
        sys.exit(1)
    else:
        s += 'celltype = "rect"\n'
        s += 'cell = """\n{0} {1} {2}\n"""\n'.format(box[0]*scale,
                                                     box[1]*scale,
                                                     box[2]*scale)
    s += 'waters = """\n'
    for name,x,y,z in filtered:
        s += "{0} {1} {2}\n".format(x*scale,y*scale,z*scale)
    s += '"""\n'
    s += 'coord = "absolute"\n'
    s += 'bondlen = 3\n'
    density = len(filtered)*18.0/(box[0]*box[1]*box[2]*scale**3*1e-24*6.022e23)
    s += 'density = {0}\n'.format(density)
    f.write(s)
    
if __name__ == "__main__":
    fNameIn, fNameOut, Nbox, make_rect_box = read_cif.parse_commandline(filetypes = ['.py',])
    atoms, box = read_cif.read_and_process(fNameIn, fNameOut, Nbox, make_rect_box)
    fOut = open(fNameOut, "w")
    write_py(atoms, box, fOut, matchfunc=lambda x: x[0] == "T")
    
