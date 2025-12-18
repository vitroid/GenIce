# coding: utf-8

from logging import getLogger
import numpy as np
from genice3.genice import GenIce3
from io import TextIOWrapper

desc = {
    "ref": {"K1939": "J. G. Kirkwood, J. Chem. Phys. 7, 911 (1939)."},
    "brief": "Kirkworrd G factor.",
    "usage": "Kirkwood G is a convenient index for testing the long-range order.",
}


# It should be expressed as a function of distance.


def calculate(genice: GenIce3) -> np.ndarray:
    "Kirkwood G."
    logger = getLogger("genice3.exporter._KG")
    zvec = np.array([0.0, 0.0, 1.0])
    dipoles = zvec @ genice.orientations
    n = dipoles.shape[0]

    logger.info(f"  Inner products of the dipoles. {n=}")
    # inner products of dipole pairs: ip[i, j] = dipoles[i] @ dipoles[j]
    ip = (dipoles @ dipoles.T).reshape(n,n)

    denominator = np.mean(np.diag(ip))

    logger.info("  Pair distances.")
    # pair distances: delta[i, j] = coord[j] - coord[i] (with periodic wrapping)
    coord = genice.lattice_sites
    # Calculate differences for all pairs: delta[i, j, :] = coord[j, :] - coord[i, :]
    delta = coord[:, np.newaxis, :] - coord[np.newaxis, :, :]
    # Apply periodic wrapping
    delta -= np.floor(delta + 0.5)
    # Convert to Cartesian coordinates and calculate distances
    delta = delta @ genice.cell
    delta = np.sqrt(np.sum(delta * delta, axis=2)).reshape(n, n)

    logger.info("  G(r)")
    # Distance bin parameters (in units of lattice spacing)
    # These values are arbitrary but commonly used for radial distribution analysis
    BIN_START = 0.04  # Starting distance for binning
    BIN_WIDTH = 0.08  # Width of each distance bin

    mx = np.max(delta)
    # Generate distance bins
    x_bins = np.arange(BIN_START, mx, BIN_WIDTH)

    # Calculate G(r) for each distance bin
    y_values = []
    for x in x_bins:
        mask = delta < x
        selected_ip = ip[mask]
        n = np.sum(mask)
        y_values.append(np.mean(selected_ip) / denominator)

    return np.array([x_bins, y_values]).T


def dumps(genice: GenIce3) -> str:
    result = calculate(genice)
    x_bins = result[:, 0]
    y_values = result[:, 1]
    # Format output as string
    return "\n".join(["{0:.5f} {1:.5f}".format(x, y) for x, y in zip(x_bins, y_values)])


def dump(genice: GenIce3, file: TextIOWrapper, **kwargs):
    file.write(dumps(genice))
