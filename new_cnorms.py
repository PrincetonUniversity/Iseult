import numpy as np
import matplotlib.colors as colors

class PowerNormWithNeg(colors.Normalize):

    ''' Custom color Norm: This one normalizes the negative data differently
    from the positive, and extends the power norm to negative values.  The main
    idea is that it plots the power_norm of the data.  If stretch_colors is
    true, then for diverging cmaps, the entire cmap is used, otherwise it only
    uses the amount so that -b and b are the same distance from the midpoint.'''

    def __init__(self, gamma = 1.0, vmin=None, vmax=None, clip=False, div_cmap = True, midpoint = 0.0, stretch_colors = True):
        colors.Normalize.__init__(self, vmin, vmax, clip)
        self.gamma = gamma
        self.div_cmap = div_cmap
        self.midpoint = midpoint
        self.stretch_colors = stretch_colors

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        # First see if there is a sign change:
        ans = PowerNormFunc(value, gamma = self.gamma, vmin=self.vmin, vmax=self.vmax, div_cmap = self.div_cmap,  midpoint = self.midpoint, stretch_colors = self.stretch_colors)
        if type(value) == np.ma.core.MaskedArray:
            ans.mask = value.mask
        return ans

def PowerNormFunc(data, gamma = 1.0, vmin=None, vmax=None, div_cmap = True,  midpoint = 0.0, stretch_colors = True):

    ''' Helper function for the PowerNorm  The main idea is that it norms data
    using np.sign(data-midpoint)*np.abs(data-midpoint)**gamma.  If
    stretch_colors is true, then for diverging cmaps, the entire cmap is used.
    If stretch colors is false then it so that -b+midpoint and b-midpoint are
    the same distance from the midpoint in the color space.'''

    vmin -= midpoint
    vmax -= midpoint
    left_clip = 0.0
    right_clip = 1.0
    if not stretch_colors:
        if np.sign(vmin) < 0 and np.sign(vmax) > 0:
            v_absmax = max(np.abs(vmin),np.abs(vmax))
            left_clip = 0.5*(1 - np.abs(vmin)**gamma/np.abs(v_absmax)**gamma)
            right_clip = 0.5*(1 + np.abs(vmax)**gamma/np.abs(v_absmax)**gamma)

    if div_cmap == True:
        if np.sign(vmin) != np.sign(vmax) and np.sign(vmin) != 0 and np.sign(vmax) != 0:
            x, y = [np.sign(vmin)*np.abs(vmin)**gamma,
                    0,
                    np.sign(vmax)*np.abs(vmax)**gamma],[left_clip, 0.5, right_clip]
        elif  np.sign(vmin) >= 0:
            # The data must be totally positive, extrapolate from midpoint
            x, y = [np.sign(vmin)*np.abs(vmin)**gamma, np.sign(vmax)*np.abs(vmax)**gamma], [0.5, right_clip]
        elif  np.sign(vmax) <= 0:
            # The data must be totally negative
            x, y = [np.sign(vmin)*np.abs(vmin)**gamma, np.sign(vmax)*np.abs(vmax)**gamma], [left_clip, 0.5]
    else:
        x, y = [np.sign(vmin)*np.abs(vmin)**gamma, np.sign(vmax)*np.abs(vmax)**gamma], [0, 1]
    if np.abs(midpoint)<=1E-8 and np.abs(gamma-1.0)<=1E-8:
        ans = np.ma.masked_array(np.interp(data, x, y))
    elif np.abs(gamma-1.0)<=1E-8:
        ans = np.ma.masked_array(np.interp(data-midpoint, x, y))
    else:
        ans = np.ma.masked_array(np.interp(np.sign(data-midpoint)*np.abs(data-midpoint)**gamma, x, y))
    return ans

class SymLogNorm(colors.Normalize):

    ''' Custom color Norm: This one normalizes the negative data differently
    from the positive, and extends the log norm to negative values.  The main
    idea is that it plots the log of the data, and stitches up around the midpoint, using a linear mapping.  If stretch_colors is
    true, then for diverging cmaps, the entire cmap is used, otherwise it only
    uses the amount so that -b and b are the same distance from the midpoint in color space.'''

    def __init__(self, vmin=None, vmax=None, clip=False, div_cmap = True, midpoint = 0.0, linthresh = 1E-3, linscale = 0.2, stretch_colors = True):
        colors.Normalize.__init__(self, vmin, vmax, clip)
        self.div_cmap = div_cmap
        self.thresh = linthresh
        self.linscale = linscale
        self.midpoint = midpoint
        self.stretch_colors = stretch_colors

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        # First see if there is a sign change:
        ans = SymLogNormFunc(value, vmin=self.vmin, vmax=self.vmax, div_cmap = self.div_cmap,  linthresh = self.thresh, linscale = self.linscale, stretch_colors = self.stretch_colors)
        if type(value) == np.ma.core.MaskedArray:
            ans.mask = value.mask
        return ans


def SymLogNormFunc(data,  vmin=None, vmax=None, div_cmap = True,  midpoint = 0.0, linthresh = 1E-3, linscale = 0.2, stretch_colors = True):

    ''' Helper function for the SymLogNorm  The main idea is that it norms data
    using np.sign(data-midpoint)*np.abs(data-midpoint)**gamma.  If
    stretch_colors is true, then for diverging cmaps, the entire cmap is used.
    If stretch colors is false then it so that -b+midpoint and b-midpoint are
    the same distance from the midpoint in the color space.'''

    vmin -= midpoint
    vmax -= midpoint
    left_clip = 0.0
    right_clip = 1.0
    left_lin = 0.5*(1-linscale)
    right_lin = 0.5*(1 + linscale)

    if not stretch_colors:
        if np.sign(vmin) < 0 and np.sign(vmax) > 0:
            v_absmax = max(np.abs(vmin),np.abs(vmax))
            left_clip = 0.5*(1-linscale)*(1-np.log10(np.abs(vmin/linthresh))/np.log10(v_absmax/linthresh))
            right_clip = .5+.5*linscale+0.5*(1-linscale)*np.log10(np.abs(vmax)/linthresh)/np.log10(v_absmax/linthresh)

    ans = np.copy(data)
    left_log=np.where(data<midpoint-linthresh)
    #####
    ans[left_log] = np.interp(np.sign(data[left_log]-midpoint)*np.abs(np.log10(np.abs(data[left_log]-midpoint)/linthresh)),
                                [np.sign(vmin)*np.log10(np.abs(vmin/linthresh)), np.sign(midpoint-linthresh)*np.log10(np.abs(midpoint-linthresh)/linthresh)],
                                [left_clip, left_lin])

    #Now the linear middle
    linear_section = np.where(np.abs(data-midpoint)<=linthresh)
    ans[linear_section] = np.interp(data[linear_section]-midpoint,[midpoint-linthresh,midpoint+linthresh],[left_lin, right_lin])

    # Now the right part
    right_log=np.where(data>=midpoint+linthresh)
    ans[right_log] = np.interp(np.sign(data[right_log]-midpoint)*np.abs(np.log10(np.abs(data[right_log]-midpoint)/linthresh)),
                                [np.sign(midpoint+linthresh)*np.log10(np.abs(midpoint+linthresh)/linthresh), np.sign(vmax)*np.log10(np.abs(vmax/linthresh))],
                                [right_lin, right_clip])



    return np.ma.masked_array(ans)

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matplotlib.mlab import bivariate_normal
    import new_cmaps
    N = 1000
    X, Y = np.mgrid[-3:3:complex(0, N), -2:2:complex(0, N)]
    Z1 = (bivariate_normal(X, Y, 1., 1., 1.0, 1.0))**2  \
    - 0.4 * (bivariate_normal(X, Y, 1.0, 1.0, -1.0, 0.0))**2
    Z1 = Z1/0.03
    print(Z1.size)

    fig, ax = plt.subplots(2, 1)


    pcm = ax[0].pcolormesh(X, Y, Z1,
                       norm=MySymLogNorm(linthresh= 1E-1, linscale = .1, stretch_colors = False),
                       cmap=new_cmaps.cmaps['plasma'])
    fig.colorbar(pcm, ax=ax[0], extend='both')

    pcm = ax[1].pcolormesh(X, Y, Z1, cmap=new_cmaps.cmaps['plasma'], vmin=-np.max(Z1))
    fig.colorbar(pcm, ax=ax[1], extend='both')
    plt.show()
