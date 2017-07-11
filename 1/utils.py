import numpy as np
import os
from random import shuffle

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
halos =  ('M0125', 'M0204', 'M0290', 'M0616', 'M0858', 'M0977', 'M1196', 'M1859', 'M4349',
         'M0175', 'M0408', 'M0664', 'M0959', 'M1646', 'M4323', 'M6782')
#halos =  ('M6782')
#halos =  ('M1196',)
types =  ('disc-Uebler', 'disc', 'ball', 'ism')
combinations = []
for halo in halos:
    for type in types:
        combinations.append((halo, type))
shuffle(combinations)

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
