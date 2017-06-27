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
