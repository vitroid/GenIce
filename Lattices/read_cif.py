#!/usr/bin/python
#
# Python script that converts a CIF file (Crystallographic Information File)
# into a configuration file for Gromacs or LAMMPS.
#
# Copyright (C) 2016  Erik Lascaris
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#                              Written by Erik Lascaris (erikl-AT-bu-DOT-edu)
#                              Version 24-July-2014
#
# NOTE: this script uses the PyCifRW code from  https://bitbucket.org/jamesrhester/pycifrw/
#
#
# TODO:
# - Double-check that it produces a good Gromacs file, and indicate on
#   website how to use it in a Gromacs simulation.
# - Double-check that it produces a good LAMMPS file, and indicate on
#   website how to use it in a LAMMPS simulation.
#
# =============================================================================

import sys

from math import *

# Import python files from the "pycifrw-4.1.1-min" directory.
sys.path.insert(0, 'pycifrw-4.1.1-min')
from CifFile import CifFile, CifBlock


# =============================================================================
# =============================================================================

# Tell the user how to use this script, and exits the script.
def print_usage():

    print('*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *')
    print('This script reads a Crystallographic Information File (CIF) that describes a crystal,')
    print('and creates a configuration file that can be used to start a Gromacs or LAMMPS simulation.')
    print('')
    print('Copyright (C) 2014  Erik Lascaris')
    print('')
    print('This program comes with ABSOLUTELY NO WARRANTY; see details inside this python script.')
    print('This is free software, and you are welcome to redistribute it under certain conditions.')
    print('')
    print('                             --- HOW TO USE  ---')
    print('Examples:')
    print('  %s  -i crystal.cif  -o unitcell.xyz  -r' % sys.argv[0])
    print('  %s  -i crystal.cif  -o conf.gro  -b 5 5 5  -r' % sys.argv[0])
    print('  %s  -i crystal.cif  -o unitcell.cif' % sys.argv[0])
    print('')
    print('List of optional arguments:')
    print('  -i  filename     CIF file with description of crystal')
    print('  -o  filename     Output configuration file.  Extension must be one of the following:')
    print('                     .xyz         XYZ file')
    print('                     .lammpstrj   LAMMPS trajectory file')
    print('                     .gro         Gromacs file')
    print('                     .cif         copy of original CIF file with additional atoms')
    print('                                  (looks nice with Jmol)')
    print('  -b int int int   Box size in terms of the unit cell. For example, if the unit cell')
    print('                   consists of 2 atoms and has size 1x2x3 A, then -b 4 5 6 creates a box')
    print('                   of size 4x10x18 A consisting of 2x4x5x6=240 atoms. Default is -b 1 1 1')
    print('  -r               Make the box rectangular (default is shape of the unit cell)')
    print('*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *')

    exit(1)


# =============================================================================

# Shows an error message and the usage.
def print_error(msg):

    print('')
    print('    ERROR:  %s' % msg)
    print('')
    print_usage()


# =============================================================================

# Converts an "_atom_type_label" into an element name.
def extract_element(label):

    elem2 = ['He','Li','Be','Ne','Na','Mg','Al','Si','Cl','Ar','Ca','Sc','Ti',
             'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
             'Rb','Sr','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
             'Sb','Te','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',
             'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','Re','Os','Ir','Pt',
             'Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',
             'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']

    if (label[0:2] in elem2):
        return label[0:2]

    elem1 = ['H','B','C','N','O','F','P','S','K','V','Y','I','W','U']

    if (label[0] in elem1):
        return label[0]

    print('WARNING: could not convert "%s" into element name!' % label)
    return label

# =============================================================================

# Basic XYZ format.
def write_xyz(atoms, box, f):

    # Write the number of atoms.
    N = len(atoms)
    f.write('%d\n' % N)

    # Write a comment line with the box size.
    if (len(box) == 6):
        # box = (ax,bx,by,cx,cy,cz)
        ax = box[0]
        bx = box[1]
        by = box[2]
        cx = box[3]
        cy = box[4]
        cz = box[5]
        f.write('Crystal created from CIF file. Box vectors:')
        f.write(' a= %10.5f %10.5f %10.5f' % (ax, 0.0, 0.0))
        f.write(' b= %10.5f %10.5f %10.5f' % (bx, by, 0.0))
        f.write(' c= %10.5f %10.5f %10.5f\n' % (cx, cy, cz))
    else:
        # box = (ax,by,cz)
        f.write('Crystal created from CIF file. Box size:') 
        f.write(' %10.5f %10.5f %10.5f\n' % box)

    # Write atom data (units are Angstroms).
    # The argument "atoms" has format ('Si', x, y, z) for example
    for i in range(N):
        f.write('%-10s %10.6f %10.6f %10.6f\n' % atoms[i])


# =============================================================================


# LAMMPS *.lammpstrj format.  Can be displayed in VMD, or used to start a new
# simulation in LAMMPS using read_dump().
def write_lammpstrj(atoms, box, f):


    # Triclinic crystal structures are often defined using three lattice constants a,
    # b, and c, and three angles alpha, beta and gamma.  Note that in this
    # nomenclature, the a, b, and c lattice constants are the scalar lengths of the
    # edge vectors a, b, and c defined above.  The relationship between these 6
    # quantities (a,b,c,alpha,beta,gamma) and the LAMMPS box sizes (lx,ly,lz) =
    # (xhi-xlo,yhi-ylo,zhi-zlo) and tilt factors (xy,xz,yz) is as follows: 
    #
    #   a = lx
    #   b^2 = ly^2 + xy^2
    #   c^2 = lz^2 + xz^2 + yz^2
    #   cos(alpha) = (xy*xz +ly*yz) / (b*c)
    #   cos(beta) = xz / c
    #   cos(gamma) = xy / b
    #
    # The inverse relationship can be written as follows: 
    #
    #   lx = a
    #   xy = b*cos(gamma)
    #   xz = c*cos(beta)
    #   ly^2 = b^2 - xy^2
    #   yz = (b*c*cos(alpha) - xy*xz) / ly
    #   lz^2 = c^2 - xz^2 - yz^2
    #
    #   ****  Info for lammpstrj file:
    #
    # ITEM: BOX BOUNDS xy xz yz
    # xlo_bound xhi_bound xy
    # ylo_bound yhi_bound xz
    # zlo_bound zhi_bound yz 
    #
    # This bounding box is convenient for many visualization programs and is
    # calculated from the 9 triclinic box parameters
    # (xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz) as follows:
    #
    #   xlo_bound = xlo + MIN(0.0,xy,xz,xy+xz)
    #   xhi_bound = xhi + MAX(0.0,xy,xz,xy+xz)
    #   ylo_bound = ylo + MIN(0.0,yz)
    #   yhi_bound = yhi + MAX(0.0,yz)
    #   zlo_bound = zlo
    #   zhi_bound = zhi 
    

    # Write timestep (always zero).
    f.write('ITEM: TIMESTEP\n')
    f.write('0\n')

    # Write the number of atoms.
    N = len(atoms)
    f.write('ITEM: NUMBER OF ATOMS\n')
    f.write('%d\n' % N)


    # Write the box size.  See
    #    lammps.sandia.gov/doc/Section_howto.html#howto_12
    # for more information about the box size in *.lammpstrj files.
    if (len(box) == 6):
        # box = (ax,bx,by,cx,cy,cz)
        # For triclinic boxes LAMMPS uses the format
        #   ITEM: BOX BOUNDS xy xz yz
        #   xlo_bound xhi_bound xy
        #   ylo_bound yhi_bound xz
        #   zlo_bound zhi_bound yz 
        # This bounding box is convenient for many visualization programs and
        # is calculated from the 9 triclinic box parameters
        # (xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz) as follows:
        #   xlo_bound = xlo + MIN(0.0,xy,xz,xy+xz)
        #   xhi_bound = xhi + MAX(0.0,xy,xz,xy+xz)
        #   ylo_bound = ylo + MIN(0.0,yz)
        #   yhi_bound = yhi + MAX(0.0,yz)
        #   zlo_bound = zlo
        #   zhi_bound = zhi 
        # where the "tilt factors" (xy,xz,yz) are given by
        #   lx = a                 = ax
        #   xy = b*cos(gamma)      = bx
        #   xz = c*cos(beta)       = cx
        #   ly^2 = b^2 - xy^2      = by^2
        #   yz = (b*c*cos(alpha) - xy*xz) / ly   = cy
        #   lz^2 = c^2 - xz^2 - yz^2             = cz^2
        # and
        #   xlo = 0,  xhi = lx,  etc.
        lx = box[0]  # ax
        xy = box[1]  # bx
        ly = box[2]  # by
        xz = box[3]  # cx
        yz = box[4]  # cy
        lz = box[5]  # cz
        f.write('ITEM: BOX BOUNDS xy xz yz\n')
        f.write('%f %f %f\n' % (min(0.0,xy,xz,xy+xz), lx + max(0.0,xy,xz,xy+xz), xy))
        f.write('%f %f %f\n' % (min(0.0,yz), max(0.0,yz), xz))
        f.write('%f %f %f\n' % (0.0, lz, yz))

    else:
        # box = (ax,by,cz)
        f.write('ITEM: BOX BOUNDS pp pp pp\n')
        f.write('0.0 %f\n' % box[0])  # Format:  "xlo xhi"
        f.write('0.0 %f\n' % box[1])  # Format:  "ylo yhi"
        f.write('0.0 %f\n' % box[2])  # Format:  "zlo zhi"

    # Write atom data.
    f.write('ITEM: ATOMS element id x y z\n')
    for i in range(N):
        (name,x,y,z) = atoms[i]  # = ('Si', x, y, z) for example
        f.write('%-5s %5d %10.6f %10.6f %10.6f\n' % (name, i+1, x, y, z))

    f.close()


# =============================================================================


# Gromacs *.gro format.
def write_gro(atoms, box, f):

    # Write the header.
    f.write('Crystal created from CIF file, t= 0\n')

    # Write the number of atoms.
    N = len(atoms)
    f.write('%d\n' % N)

    # Write atom data.  Note that Gromacs uses nm, while CIF uses Angstrom.
    for i in range(N):
        (name,x,y,z) = atoms[i]  # = ('Si', x, y, z) for example
        f.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n'
               % (i+1, name, name, i+1, 0.1*x, 0.1*y, 0.1*z))

    # Write the box size.  In a *.gro file it has the (free) format
    #    ax  by  cz  0  0  bx  0  cx  cy
    # and the last 6 columns are optional.  Here ax actually represents the
    # x-component of the unit cell a-vector times the number of unit cells in
    # the x-direction, i.e. ax = a[x] * Nx.
    if (len(box) == 6):
        # box = (ax,bx,by,cx,cy,cz)
        ax = 0.1*box[0]
        bx = 0.1*box[1]
        by = 0.1*box[2]
        cx = 0.1*box[3]
        cy = 0.1*box[4]
        cz = 0.1*box[5]
        f.write('%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n'
                                    % (ax,by,cz, 0.0,0.0,bx, 0.0,cx,cy))
    else:
        # box = (ax,by,cz)
        f.write('%10.5f%10.5f%10.5f\n' % (0.1*box[0], 0.1*box[1], 0.1*box[2]))


# =============================================================================


# Writes an extended CIF file.  It basically makes a copy of the original file,
# but adds additional atoms to complete the unit cell.
# Note that the "atoms" are in fact the fractional coordinates of the atoms,
# not the actual positions.
def write_cif(fNameIn, atoms, fNameOut):

    # Simply copy each line of the input file "fIn" to the output file "fOut",
    # except for the few lines we want to replace.  Those lines is the info
    # such ashere:
    #
    #       'x-y,-y,-z'
    # loop_
    # _atom_site_label
    # _atom_site_fract_x
    # _atom_site_fract_y
    # _atom_site_fract_z
    # Si   0.46970   0.00000   0.00000     <-- REPLACE THIS
    # O   0.41350   0.26690   0.11910      <-- REPLACE THIS
    # loop_
    # _atom_site_aniso_label
    #    etc.
    try:

        fIn = open(fNameIn, 'r')
        fOut = open(fNameOut, 'w')

        # Keep track of where we are: inside the loop with "_atom_site_" and
        # its data, or outside of this loop.
        inside_atom_site_loop = False

        # Check each line of the input file.
        for line in fIn:

            # Split line into columns.
            cols = line.split()

            # Search for the keyword "_atom_site_label".  We expect this to be
            # the start of the atom_site data.  Instead of writing this key,
            # write the 4 keys and all the data right here and now.
            if (len(cols)>0  and  cols[0] == '_atom_site_label'):

                fOut.write('_atom_site_label\n')
                fOut.write('_atom_site_fract_x\n')
                fOut.write('_atom_site_fract_y\n')
                fOut.write('_atom_site_fract_z\n')

                # Remember that atom[i] = (label_i, x_i, y_i, z_i).
                for atom in atoms:
                    fOut.write('%s   %.5f   %.5f   %.5f\n' % atom)

                # Clearly, we are inside the loop with the atom_site data.
                inside_atom_site_loop = True


            # Else, if we are inside the atom_site loop...
            elif (inside_atom_site_loop):

                # ... then ignore all lines with 4 or more columns, and ignore
                # all lines that start with "_atom_site_".  But any other line
                # is an indication that we've reached the end of the loop.
                if (len(cols) < 4  and  (len(line) < 11  or  line[:11] != '_atom_site_')):

                    inside_atom_site_loop = False
                    fOut.write(line)


            # In all other cases, just write the line to the output file.
            else:
                fOut.write(line)


        fIn.close()
        fOut.close()

    except:
        print_error('Failed to write to output file')


# =============================================================================


# Read CIF file, and extract the necessary info in the form of a dictionary.
# E.g., the value of "_cell_volume" can be found with data['_cell_volume'].
def read_cif(fNameIn):

    data = {}


    # Open the CIF file and read all the lines into a list of strings.
    try:
        f = open(fNameIn, 'r')
        lines = []
        for line in f:
            stripped = line.strip()
            if (len(stripped) > 0):  lines.append(stripped)
    except:
        print "Failed to open CIF file '{0}'".format(fNameIn)
        sys.exit()

    # Use the CifFile parser to extract the data.  Although there might be
    # multiple data blocks, we'll only use the first one.
    cif_file = CifFile(fNameIn)
    
    for db in cif_file:
        data_block = db
        break


    try:

        # Extract some parameters, and convert them to floats.
        data['_cell_length_a']    = float(data_block['_cell_length_a'])
        data['_cell_length_b']    = float(data_block['_cell_length_b'])
        data['_cell_length_c']    = float(data_block['_cell_length_c'])
        data['_cell_angle_alpha'] = float(data_block['_cell_angle_alpha'])
        data['_cell_angle_beta']  = float(data_block['_cell_angle_beta'])
        data['_cell_angle_gamma'] = float(data_block['_cell_angle_gamma'])
        data['_cell_volume']      = float(data_block['_cell_volume'])


        # Get the symbolic operations that define the space group.  In a CIF file
        # that's the part that looks like:
        #
        # loop_
        # _symmetry_equiv_pos_as_xyz
        #   'x,y,z'
        #   'y,x,2/3-z'
        #   '-y,x-y,2/3+z'
        #   '-x,-x+y,1/3-z'
        #   '-x+y,-x,1/3+z'
        #   'x-y,-y,-z'
        #
        # In some cases it's called "_space_group_symop_operation_xyz" apparently?!?!
        data['_symmetry_equiv_pos_as_xyz'] = []

        try:
            xyz = data_block["_symmetry_equiv_pos_as_xyz"]

        except KeyError:
            try:
                xyz = data_block["_space_group_symop_operation_xyz"]
            except KeyError:
                print "Missing item in CIF file: need either '_symmetry_equiv_pos_as_xyz' or '_space_group_symop_operation_xyz'."
                sys.exit()


        # Copy the x,y,z symmetry group operations.  Remove the quotes if there
        # are any.
        for op_xyz in xyz:

            if (op_xyz[0] == '\''):
                data['_symmetry_equiv_pos_as_xyz'].append(op_xyz[1:-1])
            else:
                data['_symmetry_equiv_pos_as_xyz'].append(op_xyz)


        # Add x,y,z of the atoms to "data", but make sure to convert
        # e.g. "0.1549(8)" to "0.1549".
        data['_atom_site_label'] = data_block['_atom_site_label']

        data['_atom_site_fract_x'] = []
        for str_x in data_block['_atom_site_fract_x']:
            data['_atom_site_fract_x'].append( float(str_x.split('(')[0]) )

        data['_atom_site_fract_y'] = []
        for str_y in data_block['_atom_site_fract_y']:
            data['_atom_site_fract_y'].append( float(str_y.split('(')[0]) )

        data['_atom_site_fract_z'] = []
        for str_z in data_block['_atom_site_fract_z']:
            data['_atom_site_fract_z'].append( float(str_z.split('(')[0]) )

    
    except KeyError as e:
        print "Error!  Missing item in file."
        print e
        sys.exit()


    #print "ALL DATA:"
    #print data
    #print

    # Return the extracted data.
    return data


# =============================================================================
# =============================================================================

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Read arguments given.


# Default settings.
fNameIn = ''
fNameOut = ''
Nx = 1
Ny = 1
Nz = 1
make_rect_box = False


# Read the arguments.  We expect at least 4.
if (len(sys.argv) <= 4):
    print_usage()

i = 1
while (i < len(sys.argv)):


    # Check if the name of the input file was given.
    if (sys.argv[i] == '-i'):
        
        # Make sure a file name is given.
        if (i+1 == len(sys.argv)):
            print_error('no input file name given')
        
        fNameIn = sys.argv[i+1]
        i = i + 2


    # Check if the name of the output file was given.
    elif (sys.argv[i] == '-o'):

        # Make sure a file name is given.
        if (i+1 == len(sys.argv)):
            print_error('no output file name given')

        # Check we have a valid file extension.
        fNameOut = sys.argv[i+1]
        unknown = True

        for ext in ['.xyz', '.lammpstrj', '.gro', '.cif']:
            if (fNameOut.endswith(ext)):
                unknown = False

        if (unknown):
            print_error('unknown file extension of output file')

        i = i + 2


    # Check if the box size was given.
    elif (sys.argv[i] == '-b'):

        # Make sure 3 integers are given.
        if (i+3 >= len(sys.argv)):
            print_error('need 3 integers to indicate box size')

        Nx = int(sys.argv[i+1])
        Ny = int(sys.argv[i+2])
        Nz = int(sys.argv[i+3])

        if (Nx == 0  or  Ny == 0  or  Nz == 0):
            print_error('box size integers need to be larger than zero')

        i = i + 4


    # Check if the final configuration should be in a rectangular shape, or in
    # the same shape as the unit cell.
    elif (sys.argv[i] == '-r'):

        make_rect_box = True
        i = i + 1


    # Anything else is wrong.
    else:
        print_error('invalid argument "%s"' % sys.argv[i])


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Read input file.


# Make sure an input file was given.
if (fNameIn == ''):
    print_error('no input file given.  Use:  -i filename')


# Open the CIF file and read the data.
data = read_cif(fNameIn)


# Extract lengths and angles from the CIF file.
La = float(data['_cell_length_a'])
Lb = float(data['_cell_length_b'])
Lc = float(data['_cell_length_c'])
alpha = radians(float(data['_cell_angle_alpha']))
beta = radians(float(data['_cell_angle_beta']))
gamma = radians(float(data['_cell_angle_gamma']))
volume = float(data['_cell_volume'])


# Extract the symmetry operations.  This will be a list of strings such as:
#    ['x,y,z', 'y,x,2/3-z', '-y,x-y,2/3+z', '-x,-x+y,1/3-z', ... ]
ops = data['_symmetry_equiv_pos_as_xyz']

# For proper evaluation, we need to convert "2/3" into "2./3", etc. to prevent
# integer division which would turn e.g. 2/3 into 0.
for i in range(len(ops)):
    ops[i] = ops[i].replace("0/", "0./") # also for e.g. 10/9
    ops[i] = ops[i].replace("1/", "1./")
    ops[i] = ops[i].replace("2/", "2./")
    ops[i] = ops[i].replace("3/", "3./")
    ops[i] = ops[i].replace("4/", "4./")
    ops[i] = ops[i].replace("5/", "5./")
    ops[i] = ops[i].replace("6/", "6./")
    ops[i] = ops[i].replace("7/", "7./")
    ops[i] = ops[i].replace("8/", "8./")
    ops[i] = ops[i].replace("9/", "9./")
#    ops[i] = ops[i].replace("/", "./")


# Get the atom labels and coordinates.
labels = data['_atom_site_label']
fX = [ float(s) for s in data['_atom_site_fract_x'] ]
fY = [ float(s) for s in data['_atom_site_fract_y'] ]
fZ = [ float(s) for s in data['_atom_site_fract_z'] ]

# Create a list of 4-tuples, where each tuple is an atom:
#   [ ('Si', 0.4697, 0.0, 0.0),  ('O', 0.4135, 0.2669, 0.1191),  ... ]
atoms = [ (labels[i], fX[i], fY[i], fZ[i]) for i in range(len(labels)) ]

# Make sure that all atoms lie within the unit cell.  Also convert names such
# as 'Oa1' into 'O'.
for i in range(len(atoms)):
    (name,xn,yn,zn) = atoms[i]
    xn = (xn + 10.0) % 1.0
    yn = (yn + 10.0) % 1.0
    zn = (zn + 10.0) % 1.0
    name = extract_element(name)
    atoms[i] = (name,xn,yn,zn)


# Update the user.
print 'Loaded a CIF file with %d atom coordinates and %d symmetry operations.' % (len(atoms), len(ops))


# Just for reference, here is a typical example of a CIF file:
"""
_cell_length_a 4.916
_cell_length_b 4.916
_cell_length_c 5.4054
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 120
_cell_volume 113.131
_exptl_crystal_density_diffrn      2.646
_symmetry_space_group_name_H-M 'P 32 2 1'
loop_
_space_group_symop_operation_xyz
  'x,y,z'
  'y,x,2/3-z'
  '-y,x-y,2/3+z'
  '-x,-x+y,1/3-z'
  '-x+y,-x,1/3+z'
  'x-y,-y,-z'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si   0.46970   0.00000   0.00000
O   0.41350   0.26690   0.11910
"""


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Use symmetry operations to create the unit cell.


# The CIF file consists of a few atom positions plus several "symmetry
# operations" that indicate the other atom positions within the unit cell.  So
# using these operations, create copies of the atoms until no new copies can be
# made.


# Two atoms are on top of each other if they are less than "eps" away.
eps = 0.01  # in Angstrom


# For each atom, apply each symmetry operation to create a new atom.
imax = len(atoms)
i=0
while (i < imax):

    label,x,y,z = atoms[i]

    for op in ops:

        # Python is awesome: calling e.g. eval('x,y,1./2+z') will convert the
        # string into a 3-tuple using the current values for x,y,z!
        xn,yn,zn = eval(op)

        # Make sure that the new atom lies within the unit cell.
        xn = (xn + 10.0) % 1.0
        yn = (yn + 10.0) % 1.0
        zn = (zn + 10.0) % 1.0

        # Check if the new position is actually new, or the same as a previous
        # atom.
        new_atom = True
        for at in atoms:
            if (abs(at[1]-xn) < eps  and  abs(at[2]-yn) < eps  and  abs(at[3]-zn) < eps):
                new_atom = False

                # Check that this is the same atom type.
                if (at[0] != label):
                    print_error('invalid CIF file: atom of type %s overlaps with atom of type %s' % (at[0],label))

        # If the atom is new, add it to the list!
        if (new_atom):
            atoms.append( (label,xn,yn,zn) )  # add a 4-tuple


    # Update the loop iterator.
    i = i + 1
    imax = len(atoms)


# Sort the atoms according to type alphabetically.
atoms = sorted(atoms, key=lambda at: at[0])
atoms.reverse()


# Done with creating the unit cell.  Update the user.
print('Created a unit cell consisting of %d atoms.' % len(atoms))

print('Fractional coordinates:')
for atom in atoms:
    print('%10s  %.3f  %.3f  %.3f' % atom)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create a larger box made of several unit cells: the super cell.


atomlist = []

for atom in atoms:

    # Get label and fractional coordinates.
    label,xf,yf,zf = atom

    for i in range(Nx):
            x = i+xf

            for j in range(Ny):
                y = j+yf

                for k in range(Nz):
                    z = k+zf
                    atomlist.append( (label,x,y,z) ) # add 4-tuple

atoms = atomlist


# If the user wants us to create a copy of the current CIF file, with
# additional atoms, then do that.  Note that the atoms here have *fractional*
# coordinates!
if (fNameOut.endswith('.cif')):
    
    write_cif(fNameIn, atoms, fNameOut)

    print('Done writing extended CIF file (%d atoms in total).' % len(atoms))
    exit(0)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Convert the fractional coordinates into real coordinates.


# The primitive vectors a,b,c are such that 
#
#   cos(alpha) = b.c / |b||c|
#   cos(beta)  = a.c / |a||c|
#   cos(gamma) = a.b / |a||b|
#
# with the convention
#
#   a = La*xhat
#   b = bx*xhat + by*yhat
#   c = cx*xhat + cy*yhat + cz*zhat
#
cosa = cos(alpha)
sina = sin(alpha)
cosb = cos(beta)
sinb = sin(beta)
cosg = cos(gamma)
sing = sin(gamma)

cosa2 = cosa * cosa
cosb2 = cosb * cosb
sing2 = sing * sing

ax = La

bx = Lb * cosg
by = Lb * sing

cx = Lc * cosb
cy = Lc * (cosa - cosg*cosb) / sing
cz = Lc * sqrt( 1 - (cosa2 + cosb2 - 2*cosg*cosb*cosa) / sing2 )


# Use the volume to check if we did the vectors right.
V = ax*by*cz
if ( abs(V - volume) > 0.1):
    print_error('volume does not match that calculated from primitive vectors')


# Check if we have a rectangular box.
if (bx < eps  and  cx < eps  and cy < eps):
    make_rect_box = True


# Update the user.
print('The primitive unit cell vectors are:')
print('   a = [%6.4f, %6.4f, %6.4f]' % (ax,0,0))
print('   b = [%6.4f, %6.4f, %6.4f]' % (bx,by,0))
print('   c = [%6.4f, %6.4f, %6.4f]' % (cx,cy,cz))
print('This gives a volume of %f A^3 (CIF file indicates it is %f A^3)' % (V,volume))


# Determine the box size.
Lx = Nx * La
Ly = Ny * Lb
Lz = Nz * Lc


for i in range(len(atoms)):

    # Get label and fractional coordinates.
    label,xf,yf,zf = atoms[i]

    xa = xf * ax  # contribution of a-vector to the x-coordinate of this atom
    #ya = 0       # a-vector has no y-component, so does not affect y of atom
    #za = 0       # a-vector has no z-component, so does not affect z of atom
    
    xb = yf * bx  # contribution of b-vector to the x-coordinate of this atom
    yb = yf * by  # contribution of b-vector to the y-coordinate of this atom
    #zb = 0       # b-vector has no z-component, so does not affect z of atom

    xc = zf * cx  # contribution of c-vector to the x-coordinate of this atom
    yc = zf * cy  # contribution of c-vector to the y-coordinate of this atom
    zc = zf * cz  # contribution of c-vector to the z-coordinate of this atom

    # Add all contributions.
    xn = xa + xb + xc
    yn = yb + yc
    zn = zc

    if (make_rect_box):
        xn = (xn + Lx) % Lx
        yn = (yn + Ly) % Ly
        zn = (zn + Lz) % Lz

    atoms[i] = (label, xn, yn, zn)


# Determine the box-vector.
if (make_rect_box):
    box = (Lx, Ly, Lz)
else:
    box = (Lx, Ly, Lz, Nx*cx, Ny*cy, Nz*cz)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create the output file.

try:
    fOut = open(fNameOut, 'w')

    if (fNameOut.endswith('.xyz')):
        write_xyz(atoms, box, fOut)

    elif (fNameOut.endswith('.lammpstrj')):
        write_lammpstrj(atoms, box, fOut)

    elif (fNameOut.endswith('.gro')):
        write_gro(atoms, box, fOut)

except:
    print_error('Failed to write to output file')


print('Created output file %s (%d atoms in total).' % (fNameOut, len(atoms)))
fOut.close()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


