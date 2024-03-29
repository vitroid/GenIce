#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
a sample for general alkyl group with methyls
tree | branch: [root,branch1,branch2,branch3]
or             [root,branch1,branch2]
or             [root,branch1]
or             [leaf]
or             leaf
v1 and v2 must be given as a unit vector.


一応、normal alkylを生成する関数を切りわけてみたものの、あまり汎用性がある気がしない。
むしろ、Use caseとして、ケージ内の残基の形を最適化する手順を例として示す?
それはそれでけっこう難しいかもしれない。(Gromacsなどを動員することになる)


"""

from logging import getLogger
from math import cos, radians, sin
from types import SimpleNamespace

import numpy as np

from genice2.cell import rel_wrap


def find_bondlen(bondlen, a, b):
    for ia in range(len(a), 0, -1):
        aa = a[:ia]
        for ib in range(len(b), 0, -1):
            bb = b[:ib]
            if (aa, bb) in bondlen:
                return bondlen[aa, bb]
            elif (bb, aa) in bondlen:
                return bondlen[bb, aa]
    assert False, f"No corresponding bond in bondlen{{}}: {a}-{b}."


def Alkyl(cage_center, root_position, cell, molname, tree, bondlen, origin_atom="N"):
    """
    put a normal-alkyl group rooted at root_position and growing toward cage_center.
    """

    def alkyl_(e1, aim, tree, origin=np.zeros(3), origin_atom="N"):
        """
        put a normal-alkyl group rooted at origin toward the aim.

        origin: 残基がくっつく基点となる原子の位置
        e1: 最初の結合の向きの単位ベクトル
        aim: 残基の向かうべき位置。通常はケージの中心
        tree: 残基のトポロジー
        bondlen: 2つの原子のtupleをキーとし、結合長を値とする辞書。
        """
        logger = getLogger()

        if not isinstance(tree, list):
            tree = [tree]

        # v2 is a vector from pivot to the aim
        v2 = cell.abs_wrap(aim - (origin + e1))
        e2 = v2 / np.linalg.norm(v2)
        e2d = e1 @ e2
        while abs(e2d) > 0.999:
            # They are inline. It is not safe to determine the orientation.
            v2 = np.random.random(3)
            e2 = v2 / np.linalg.norm(v2)
            e2d = e1 @ e2
        v2 = e2 - e2d * e1
        e2 = v2 / np.linalg.norm(v2)

        e3 = np.cross(e1, e2)  # the thild unit vector

        c = cos(radians(120))
        s = sin(radians(120))
        e4 = e2 * c + e3 * s  # a branch vector
        e5 = e2 * c - e3 * s  # another branch vector

        c = cos(radians(109.5))
        s = sin(radians(109.5))
        e2 = -e1 * c + e2 * s
        e4 = -e1 * c + e4 * s
        e5 = -e1 * c + e5 * s

        atomname = tree[0]
        length = find_bondlen(bondlen, origin_atom, atomname)
        displace = length * e1
        atoms = [(atomname, origin + displace)]
        for vec, topo in zip([e2, e4, e5], tree[1:]):
            atoms += alkyl_(
                vec, aim, topo, origin=origin + displace, origin_atom=atomname
            )
        return atoms

    logger = getLogger()
    # logger.info("  Put butyl at {0}".format(molname))
    v1abs = cell.rel2abs(rel_wrap(cage_center - root_position))
    e1 = v1abs / np.linalg.norm(v1abs)

    origin = cell.rel2abs(root_position)
    rawatoms = alkyl_(
        e1, cell.rel2abs(cage_center), tree, origin=origin, origin_atom=origin_atom
    )

    atomnames = []
    atompos = []
    order = []
    # atoms = []
    for i, atom in enumerate(rawatoms):
        atomname, pos = atom
        atomnames.append(atomname)
        atompos.append(cell.abs_wrapf(pos))
        order.append(i)
    mols = SimpleNamespace(
        resname=molname,
        atomnames=atomnames,
        positions=[atompos],  # atomic positions
        orig_order=order,  #
    )
    return mols
