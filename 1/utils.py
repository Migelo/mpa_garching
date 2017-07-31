import numpy as np
import os
from random import shuffle
import math

figsize = np.array([8.268, 8.268*2**.5])

def prepare_step(y):
    '''
    Prepare data for a nice step plot.

    Args:
        y (list):   List containing bin values.

    Returns:
        y (list):   List containing corrected bin values.
    '''
    y = np.insert(y, 0, y[0])
    
    return y

'''Generate combinations of halos with all tracing types.'''
halos = ('M0408', 'M0501', 'M0616', 'M0664', 'M0858', 'M0959',
    'M0977', 'M1192', 'M1196', 'M1646', 'M1859', 'M2283')[::-1]
#types =  ('disc-Uebler', 'disc', 'ball', 'ism')
types =  ('ism', )
combinations = []
for halo in halos:
    for type in types:
        combinations.append((halo, type))
#shuffle(combinations)

def save(path, filename):
    '''
    Create directory if it does not exist.

    Args:
        path (str): Path for the new directory.

    Returns:
        None
    '''
    if not os.path.isdir(path):
        os.mkdir(path)
    plt.savefig('%s/%s' % (path, filename))

def finite(data):
    return data[np.isfinite(data)]

def norm_hist(data, maximum):
    return np.array(data) * float(maximum) / max(data) 

def round_to_n(x, n):
    " Round x to n significant figures "
    return round(x, -int(math.floor(np.sign(x) * np.log10(abs(x)))) + n)

def str_fmt(x, n=1):
    " Format x into nice Latex rounding to n"
    if x == 0:
        return 0
    power = int(np.log10(round_to_n(x, 0)))
    f_SF = round_to_n(x, n) * pow(10, -power)
    return r"%s\cdot 10^{%s}" % (f_SF, power)

def flatten_list(array):
    return [item for sublist in array for item in sublist]
