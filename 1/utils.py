def prepare_step(x, y):
    '''
    Prepare data for a nice step plot.

    Args:
        x (list):   List concaining bin edges.
        y (list):   List containing bin values.

    Returns:
        x (list):   List containing corrected bin edges.
        y (list):   List containing corrected bin values.
    '''
    if type(x) is not list: 
        x = x.tolist()
    if type(y) is not list:
        y = y.tolist()

    y.insert(0, y[0])
#    x.insert(0, x[0])
#    y.insert(0, 0)
#    x.append(x[-1])
#    y.append(0)
    
    return x, y

halos =  ('M0125', 'M0204', 'M0290', 'M0616', 'M0858', 'M0977', 'M1196', 'M1859', 'M4349',
         'M0175', 'M0408', 'M0664', 'M0959', 'M1192', 'M1646', 'M4323', 'M6782')
types =  ('disc-Uebler', 'disc', 'ball', 'ism')
combinations = []
for halo in halos:
    for type in types:
        combinations.append((halo, type))
