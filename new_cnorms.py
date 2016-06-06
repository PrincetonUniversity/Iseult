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
                    np.sign(midpoint)*np.abs(midpoint)**gamma,
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
        ans = np.ma.masked_array(np.interp(data, x, y))
    else:
        ans = np.ma.masked_array(np.interp(np.sign(data)*np.abs(data)**gamma, x, y))
    return ans

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matplotlib.mlab import bivariate_normal
    import new_cmaps
    N = 1000
    X, Y = np.mgrid[-3:3:complex(0, N), -2:2:complex(0, N)]
    Z1 = (bivariate_normal(X, Y, 1., 1., 1.0, 1.0))**2  \
    - 0.4 * (bivariate_normal(X, Y, 1.0, 1.0, -1.0, 0.0))**2
    Z1 = Z1/0.03


    fig, ax = plt.subplots(2, 1)

    pcm = ax[0].pcolormesh(X, Y, Z1,
                       norm=PowerNormWithNeg(gamma = 1.0),
                       cmap=new_cmaps.cmaps['vort'])
    fig.colorbar(pcm, ax=ax[0], extend='both')

    pcm = ax[1].pcolormesh(X, Y, Z1, cmap=new_cmaps.cmaps['vort'], vmin=-np.max(Z1))
    fig.colorbar(pcm, ax=ax[1], extend='both')
    plt.show()
