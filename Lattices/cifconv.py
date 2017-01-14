#!/usr/bin/env python2

import read_cif

if __name__ == "__main__":
    fNameIn, fNameOut, Nbox, make_rect_box = read_cif.parse_commandline()
    atoms, box = read_cif.read_and_process(fNameIn, fNameOut, Nbox, make_rect_box)
    read_cif.write(fNameOut, atoms, box, fNameIn=fNameIn)
    
