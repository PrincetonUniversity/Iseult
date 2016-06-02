import numpy as np
import matplotlib.colors as colors

class PowerNormWithNeg(colors.Normalize):
    ''' Custom color Norm: An example with a customized normalization.
    This one normalizes the negative data differently
    from the positive, and extends the power norm to negative values. '''

    def __init__(self, gamma = 1.0, vmin=None, vmax=None, clip=False):
        colors.Normalize.__init__(self, vmin, vmax, clip)
        self.gamma = gamma
    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        # First see if there is a sign change:
        if np.sign(self.vmin) != np.sign(self.vmax):
            x, y = [np.sign(self.vmin)*np.abs(self.vmin)**self.gamma, 0.0, np.sign(self.vmax)*np.abs(self.vmax)**self.gamma], [0, 0.5, 1]
        else:
            x, y = [np.sign(self.vmin)*np.abs(self.vmin)**self.gamma, np.sign(self.vmax)*np.abs(self.vmax)**self.gamma], [0, 1]
        ans = np.ma.masked_array(np.interp(np.sign(value)*np.abs(value)**self.gamma, x, y))
        if type(value) == np.ma.core.MaskedArray:
            ans.mask = value.mask
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
