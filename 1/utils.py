import math
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pygad as pg
from tqdm import tqdm

figsize = np.array([8.268, 8.268 * 2**0.5])
tick_labelsize = mpl.rcParamsOrig["xtick.labelsize"]
axes_labelsize = mpl.rcParamsOrig["axes.labelsize"]


def prepare_step(y):
    """
    Prepare data for a nice step plot.

    Args:
        y (list):   List containing bin values.

    Returns:
        y (list):   List containing corrected bin values.
    """
    y = np.insert(y, 0, y[0])

    return y


"""Generate combinations of halos with all tracing types."""
halos = (
    "M0408",
    "M0501",
    "M0616",
    "M0664",
    "M0858",
    "M0959",
    "M0977",
    "M1192",
    "M1196",
    "M1646",
    "M1859",
    "M2283",
)
# types =  ('disc-Uebler', 'disc', 'ball', 'ism')
types = ("ism",)
combinations = []
for halo in halos:
    for type in types:
        combinations.append((halo, type))
# shuffle(combinations)


def save(path, filename):
    """
    Create directory if it does not exist.

    Args:
        path (str): Path for the new directory.

    Returns:
        None
    """
    if not os.path.isdir(path):
        os.mkdir(path)
    plt.savefig("{}/{}".format(path, filename))


def finite(data):
    return data[np.isfinite(data)]


def norm_hist(data, maximum):
    return np.array(data) * float(maximum) / max(data)


def round_to_n(x, n):
    "Round x to n significant figures"
    return round(x, -int(math.floor(np.sign(x) * np.log10(abs(x)))) + n)


def str_fmt(x, n=1):
    "Format x into nice Latex rounding to n"
    if x == 0:
        return 0
    power = int(np.log10(round_to_n(x, 0)))
    f_SF = round_to_n(x, n) * pow(10, -power)
    return r"{}\cdot 10^{{{}}}".format(f_SF, power)


def flatten_list(array):
    return [item for sublist in array for item in sublist]


def getsize(obj):
    """sum size of object & members in MB."""
    import sys
    from gc import get_referents
    from types import FunctionType, ModuleType

    # Custom objects know their class.
    # Function objects seem to know way too much, including modules.
    # Exclude modules as well.
    BLACKLIST = type, ModuleType, FunctionType

    if isinstance(obj, BLACKLIST):
        raise TypeError("getsize() does not take argument of type: " + str(type(obj)))
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size / (1024**2)


def radial_shell(r0, dr):
    """Returns a mask including particles within r0 + dr"""
    r0 = pg.UnitScalar(r0)
    dr = pg.UnitScalar(dr)
    inner = pg.BallMask(r0)
    outer = pg.BallMask(r0 + dr)
    shell = outer & ~inner
    return shell


def covering_fractions(snap, qty, r_start, r_end, dr, treshold):
    """Calculate covering fractions for qty from r_start to
    r_end in steps of dr counting values over treshold"""
    r_start = pg.UnitScalar(r_start)
    r_end = pg.UnitScalar(r_end)
    dr = pg.UnitScalar(dr)

    bins = np.arange(r_start, r_end + dr, dr)
    qty_binned = np.zeros(bins.shape)
    total_mass_binned = np.zeros(bins.shape)
    for i, r in tqdm(enumerate(bins), total=bins.size):
        shell = radial_shell(r, dr)
        total_mass_binned[i] = snap.gas[shell]["mass"].sum()
        qty_binned[i] = snap.gas[shell][qty].sum()
    return qty_binned / total_mass_binned
