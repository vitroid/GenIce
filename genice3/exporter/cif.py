# coding: utf-8
"""
Crude Cif file format
"""

from logging import getLogger
from io import TextIOWrapper
import sys

import numpy as np

from genice3.util import cellshape
from genice3.genice import GenIce3
from genice3.exporter import (
    _parse_guest_option,
    spot_guest_processor,
    water_model_processor,
)


def _format_cell_shape(a, b, c, A, B, C):
    aL, bL, cL, alpha, beta, gamma = a, b, c, A, B, C

    s = ""
    s += "data_genice\n"
    s += "#" + "\n#Command line: " + " ".join(sys.argv) + "\n"
    s += "_cell_length_a                {0}\n".format(aL * 10)
    s += "_cell_length_b                {0}\n".format(bL * 10)
    s += "_cell_length_c                {0}\n".format(cL * 10)
    s += "_cell_angle_alpha             {0}\n".format(alpha)
    s += "_cell_angle_beta              {0}\n".format(beta)
    s += "_cell_angle_gamma             {0}\n".format(gamma)
    s += "\n"
    rights = np.array([alpha, beta, gamma]) - 90
    if np.allclose(rights, 0):
        s += "_symmetry_cell_setting        'orthorhombic'\n"
        s += "_symmetry_space_group_name_H-M   'P 1 '\n"
    else:
        # for now it is always triclinic
        s += "_symmetry_cell_setting        'triclinic'\n"
        s += "_symmetry_space_group_name_H-M   'P 1 '\n"
    return s


def dump(
    genice: GenIce3,
    file: TextIOWrapper = sys.stdout,
    guest: dict = {},
    spot_guest: dict = {},
    water_model: str = "3site",
    name: str = "",
):
    "Output in CIF format."
    logger = getLogger()
    assert name in ["cif", ""]

    a, b, c, A, B, C = cellshape(genice.cell)
    s = _format_cell_shape(a, b, c, A, B, C)

    s += """
loop_
_symmetry_equiv_pos_as_xyz
X,Y,Z

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
"""
    # オプションの処理
    guest_info = parse_guest_option(guest)
    spot_guest_info = parse_spot_guest_option(spot_guest)
    water_model = parse_water_model_option(water_model)

    # 分子の取得
    waters = genice.water_molecules(water_model=water_model)
    guests = genice.guest_molecules(guests=guest_info, spot_guests=spot_guest_info)
    ions = genice.substitutional_ions()

    # セル行列の逆行列を計算（絶対座標から相対座標への変換用）
    cell_matrix = genice.cell
    cell_inv = np.linalg.inv(cell_matrix)

    # 原子を収集
    atoms = []
    for water in waters.values():
        for name, position in zip(water.labels, water.sites):
            atoms.append([water.name, name, position])
    for guest in guests:
        for name, position in zip(guest.labels, guest.sites):
            atoms.append([guest.name, name, position])
    for ion in ions.values():
        for name, position in zip(ion.labels, ion.sites):
            atoms.append([ion.name, name, position])

    # CIF形式で出力
    for i, atom in enumerate(atoms):
        resname, atomname, position = atom
        # 絶対座標から相対座標に変換
        pos = position @ cell_inv
        pos = pos - np.floor(pos)  # 0-1の範囲に正規化
        label = f"{atomname}{i}"
        s += f"{label:>6}{atomname:>6}{pos[0]:10.4f}{pos[1]:10.4f}{pos[2]:10.4f}\n"
    s += "\n"
    file.write(s)
