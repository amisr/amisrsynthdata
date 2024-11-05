import numpy as np


def output_shape(ut, x):
    # determine the appropriate output shape for the given time and position
    # inputs
    if (np.isscalar(ut) and np.isscalar(x)):
        s = None
    elif np.isscalar(ut):
        s = x.shape
    elif np.isscalar(x):
        s = ut.shape
    else:
        s = ut.shape + x.shape
    return s
