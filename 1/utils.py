import numpy as np

def prepare_step(y):
    '''
    Prepare data for a nice step plot.

    Args:
        y (list):   List containing bin values.

    Returns:
        y (list):   List containing corrected bin values.
    '''
    y = np.insert(y, 0, y[0])
    
    return x

halos =  ('M0125', 'M0204', 'M0290', 'M0616', 'M0858', 'M0977', 'M1196', 'M1859', 'M4349',
         'M0175', 'M0408', 'M0664', 'M0959', 'M1192', 'M1646', 'M4323', 'M6782')
types =  ('disc-Uebler', 'disc', 'ball', 'ism')
combinations = []
for halo in halos:
    for type in types:
        combinations.append((halo, type))
