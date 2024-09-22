"""
Iterative file loader for analice.
"""

from logging import getLogger
import re
import os
from pathlib import Path
from genice2.plugin import safe_import
from genice2.cell import rel_wrap
import numpy as np
from types import SimpleNamespace
import pairlist as pl


def str2rangevalues(s):
    values = s.split(":")
    assert len(values) > 0
    if len(values) == 1:
        return [0, int(values[0]), 1]
    elif len(values) == 2:
        return [int(values[0]), int(values[1]), 1]
    else:
        return [int(v) for v in values]


def str2range(s):
    return range(*str2rangevalues(s))


def iterate(filename, oname, hname, filerange, framerange, suffix=None):
    logger = getLogger()
    rfile = str2range(filerange)
    rframe = str2rangevalues(framerange)
    logger.info("  file number range: {0}:{1}:{2}".format(*str2rangevalues(filerange)))
    logger.info("  frame number range: {0}:{1}:{2}".format(*rframe))
    # test whether filename has a regexp for enumeration
    logger.info(filename)
    m = re.search("%[0-9]*d", filename)
    # prepare file list
    if m is None:
        filelist = [
            filename,
        ]
    else:
        filelist = []
        for num in rfile:
            fname = filename % num
            if os.path.exists(fname):
                filelist.append(fname)
    logger.debug("File list: {0}".format(filelist))
    frame = 0
    for fname in filelist:
        logger.info("  File name: {0}".format(fname))
        # single file may contain multiple frames
        if suffix is None:
            suffix = Path(fname).suffix[1:]
        loader = safe_import("loader", suffix)
        file = open(fname)
        for oatoms, hatoms, cellmat in loader.load_iter(file, oname=oname, hname=hname):
            if frame == rframe[0]:
                logger.info("Frame: {0}".format(frame))
                yield oatoms, hatoms, cellmat
                rframe[0] += rframe[2]
                if rframe[1] <= rframe[0]:
                    return
            else:
                logger.info("Skip frame: {0}".format(frame))
            frame += 1


def average(load_iter, span=0):
    logger = getLogger()
    if span <= 1:
        # pass through.
        yield from load_iter()
        return
    ohist = []  # center-of-mass position (relative)
    for oatoms, hatoms, cellmat in load_iter():
        # center of mass
        if len(ohist) == 0:
            # first ohist; just store
            ohist.append(oatoms)
        elif span > 1:
            # displacements
            d = oatoms - ohist[-1]
            d -= np.floor(d + 0.5)
            # new positions
            ohist.append(ohist[-1] + d)
            # if too many storage
            if len(ohist) > span:
                # drop the oldest one.
                ohist.pop(0)
            # overwrite the water positions with averaged ones.
            oatoms = np.average(np.array(ohist), axis=0)
        else:
            ohist[0] = oatoms

        yield oatoms, hatoms, cellmat


# def history(load_iter, span=0):


def make_lattice_info(oatoms, hatoms, cellmat):
    logger = getLogger()

    assert oatoms.shape[0] > 0
    assert hatoms is None or oatoms.shape[0] * 2 == hatoms.shape[0]
    coord = "relative"
    density = oatoms.shape[0] / (np.linalg.det(cellmat) * 1e-21) * 18 / 6.022e23

    if hatoms is None:
        return SimpleNamespace(
            waters=oatoms, coord=coord, density=density, bondlen=0.3, cell=cellmat
        )

    rotmat = np.zeros((oatoms.shape[0], 3, 3))
    for i in range(oatoms.shape[0]):
        ro = oatoms[i]
        rdh0 = rel_wrap(hatoms[i * 2] - ro)
        rdh1 = rel_wrap(hatoms[i * 2 + 1] - ro)
        dh0 = rdh0 @ cellmat
        dh1 = rdh1 @ cellmat
        y = dh0 - dh1
        y /= np.linalg.norm(y)
        z = dh0 + dh1
        z /= np.linalg.norm(z)
        x = np.cross(y, z)
        rotmat[i] = np.vstack([x, y, z])
        # 重心位置を補正。
        oatoms[i] += (rdh0 + rdh1) * 1.0 / 18.0
    # remove intramolecular OHs
    # 水素結合は原子の平均位置で定義している。
    pairs = []
    logger.debug("  Make pair list.")
    for o, h in pl.pairs_iter(
        oatoms, maxdist=0.245, cell=cellmat, pos2=hatoms, distance=False
    ):
        if not (h == o * 2 or h == o * 2 + 1):
            # hとoは別の分子の上にあって近い。
            # register a new intermolecular pair
            pairs.append((h // 2, o))
    logger.debug("  # of pairs: {0} {1}".format(len(pairs), oatoms.shape[0]))

    return SimpleNamespace(
        waters=oatoms,
        coord=coord,
        density=density,
        pairs=pairs,
        rotmat=rotmat,
        cell=cellmat,
        __doc__=None,
    )
