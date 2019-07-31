from logging import getLogger
import re
from pathlib import Path
from genice.importer import safe_import
import numpy as np


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
        filelist = [filename, ]
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
        loader = safe_import("loader", suffix).Loader(fname, oname, hname)
        for lat in loader.load_iter():
            if frame == rframe[0]:
                logger.info("Frame: {0}".format(frame))
                lat.lattice_type = "analice"
                yield lat
                rframe[0] += rframe[2]
                if rframe[1] <= rframe[0]:
                    return
            else:
                logger.info("Skip frame: {0}".format(frame))
            frame += 1


def average(load_iter, span=0):
    logger = getLogger()
    if span < 1:
        span = 1
    com = [] # center-of-mass position (relative)
    rot = [] # rotation matrices
    for lat in load_iter():

        # center of mass
        if len(com) == 0:
            # first com; just store
            com.append(lat.waters)
        else:
            # displacements
            d = lat.waters - com[-1]
            d -= np.floor(d+0.5)
            # new positions
            com.append(com[-1]+d)
            # if too many storage
            if len(com) > span:
                # drop the oldest one.
                com.pop(0)
            # overwrite the water positions with averaged ones.
            lat.waters = np.average(np.array(com), axis=0)

        yield lat
