import numpy as np
try:
    from gemini3d.grid.gridmodeldata import model2pointsgeogcoords
    import gemini3d.read as read
except ImportError:
    print('WARNING: pygemini is not installed.  GEMINI functionality will '
          'not be available.')


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


class gemini_helper(object):
    """
    This is small class to assist with the GEMINI ionospheric state functions.
    There are a few intricacies with properly calling the gemini3d functions,
    and this class serves to limit excessive copy/paste code blocks.
    """

    def __init__(self, gemini_output_dir, glat, glon, galt):
        self.gemini_output_dir = gemini_output_dir
        self.out_shape = galt.shape
        # gemini functions can't handle nans in coordinate arrays, so remove
        #   and replace in the output array
        self.glatf, self.glonf, self.galtf = self.remove_nans(glat, glon, galt)
        # read in GEMINI grid
        self.xg = read.grid(self.gemini_output_dir)

    def remove_nans(self, glat, glon, galt):
        # find/save the locations of NaNs and remove them from coordinate
        # arrays
        glatf = glat.flatten()
        glonf = glon.flatten()
        galtf = galt.flatten()
        finite_coords = np.any(np.isfinite([glatf, glonf, galtf]), axis=0)
        # find indices where NaNs will be removed and should be inserted in new
        # arrays
        self.nan_idx = np.array(
            [r - i for i, r in
             enumerate(np.argwhere(~finite_coords).flatten())])
        # return coordinate arrays without NaNs
        return glatf[finite_coords], glonf[finite_coords], galtf[finite_coords]

    def replace_nans(self, P):
        # replace NaNs in output arrays
        if np.any(self.nan_idx):
            P = np.insert(P, self.nan_idx, np.nan)
        return P

    def query_model(self, time, param):
        # Use gemini3d fuctions to get interpolated parameters and reform the
        #   arrays correctly
        dat = read.frame(self.gemini_output_dir, time, var=param)
        P = model2pointsgeogcoords(
            self.xg,
            dat[param],
            self.galtf,
            self.glonf,
            self.glatf)
        P = self.replace_nans(P)
        P0 = P.reshape(self.out_shape)
        return P0
