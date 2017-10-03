import numpy as np
from numba import jit, guvectorize, float64, int32
from math import sqrt

@jit(nopython=True, cache = True)
def stepify(bins, hist):
    # make it a step, there probably is a much better way to do this
    tmp_hist = np.ones(2*len(hist))
    tmp_bin = np.ones(2*len(hist))
    for j in range(len(hist)):
        tmp_hist[2*j] = hist[j]
        tmp_hist[2*j+1] = hist[j]

        tmp_bin[2*j] = bins[j]
        if j != 0:
            tmp_bin[2*j-1] = bins[j]
        if j == len(hist)-1:
            tmp_bin[2*j+1] = bins[j+1]
    return tmp_bin, tmp_hist

@jit(nopython=True, cache = True)
def TwiceArr(xarr):
    # double the array, there probably is a much better way to do this
    tmp = np.empty(2*len(xarr))
    for j in range(len(xarr)):
        tmp[2*j] = xarr[j]
        tmp[2*j+1] = xarr[j]
    return tmp

@guvectorize([(float64[:], float64[:], float64[:],float64[:])], '(n),(n),(n)->(n)', nopython = True, cache = True)
def LorentzFactor(u,v,w, ans):
    # Iterate over the array, find the min and max and divide by a number
    # Calculate xmin and xmax
    for i in xrange(len(u)):
        ans[i] = u[i]*u[i]+v[i]*v[i]+ w[i]*w[i] + 1
        ans[i] = sqrt(ans[i])


@guvectorize([(float64[:], # x
               float64[:], # u
               float64[:], # v
               float64[:], # w
               float64[:], # bin_width
               float64[:], # xmin
               float64[:], # gamma (unbinned)
               float64[:], # vx binned
               float64[:], # E binned
               float64[:] # counts
               )],
               '(n),(n),(n),(n),(),(), (n), (m), (m), (m)', nopython = True, cache = True)
def CalcVxEHists(x, u, v, w, bin_width, xmin, g, vx, E, counts):
    bn = bin_width[0]
    minx = xmin[0]
    maxl = len(vx)-1
    for i in xrange(len(x)):
        l = int((x[i]-minx)//bn)
        if l > maxl:
            l = maxl
        g[i] = u[i]*u[i]+v[i]*v[i]+ w[i]*w[i] + 1
        g[i] = sqrt(g[i])
        vx[l] += u[i]*g[i]**-1
        E[l] += g[i]-1
        counts[l] += 1
    for l in xrange(len(vx)):
        c = counts[l]**-1
        if c != 0:
            vx[l] *= c
            E[l] *= c

@guvectorize([(float64[:], # x
               float64[:], # u
               float64[:], # v
               float64[:], # w
               float64[:], # weights
               float64[:], # bin_width
               float64[:], # xmin
               float64[:], # gamma (unbinned)
               float64[:], # vx binned
               float64[:], # E binned
               float64[:] # counts
               )],
               '(n),(n),(n),(n),(n), (),(),(n), (m), (m), (m)', nopython = True, cache = True)
def CalcVxEWeightedHists(x, u, v, w, weights, bin_width, xmin, g, vx,  E, counts):
    bn = bin_width[0]
    minx = xmin[0]
    maxl = len(vx)-1
    for i in xrange(len(x)):
        c1 = weights[i]
        l = int((x[i]-minx)//bn)
        if l > maxl:
            l = maxl
        g[i] = u[i]*u[i]+v[i]*v[i]+ w[i]*w[i] + 1
        g[i] = sqrt(g[i])
        vx[l] += u[i]*g[i]**-1*c1
        E[l] += (g[i]-1)*c1
        counts[l] += c1
    for l in xrange(len(vx)):
        c = counts[l]**-1
        if c != 0:
            vx[l] *= c
            E[l] *= c

@guvectorize([(float64[:], # x
               float64[:], # u
               float64[:], # v
               float64[:], # w
               float64[:], # bin_width
               float64[:], # xmin
               float64[:], # vx binned
               float64[:], # vy binned
               float64[:], # vz binned
               float64[:] # counts
               )],
               '(n),(n),(n),(n),(),(), (m),(m), (m), (m)', nopython = True, cache = True)
def CalcVHists(x, u, v, w, bin_width, xmin, vx, vy, vz, counts):
    bn = bin_width[0]
    minx = xmin[0]
    maxl = len(vx)-1
    for i in xrange(len(x)):
        l = int((x[i]-minx)//bn)
        if l > maxl:
            l = maxl
        g = u[i]*u[i]+v[i]*v[i]+ w[i]*w[i] + 1
        g = sqrt(g)**-1
        vx[l] += u[i]*g
        vy[l] += v[i]*g
        vz[l] += w[i]*g
        counts[l] += 1
    # Normalize
    for l in xrange(len(vx)):
        c = counts[l]**-1
        if c != 0:
            vx[l] *= c
            vy[l] *= c
            vz[l] *= c

@guvectorize([(float64[:], # x
               float64[:], # u
               float64[:], # v
               float64[:], # w
               float64[:], # weights
               float64[:], # bin_width
               float64[:], # xmin
               float64[:], # vx binned
               float64[:], # vy binned
               float64[:], # vz binned
               float64[:] # counts
               )],
               '(n),(n),(n),(n),(n),(),(), (m),(m), (m), (m)', nopython = True, cache = True)
def CalcVWeightedHists(x, u, v, w, weights, bin_width, xmin, vx, vy, vz, counts):
    bn = bin_width[0]
    minx = xmin[0]
    maxl = len(vx)-1
    for i in xrange(len(x)):
        l = int((x[i]-minx)//bn)
        if l > maxl:
            l = maxl
        g = u[i]*u[i]+v[i]*v[i]+ w[i]*w[i] + 1
        g = sqrt(g)**-1*weights[i]
        vx[l] += u[i]*g
        vy[l] += v[i]*g
        vz[l] += w[i]*g
        counts[l] += weights[i]
    # Normalize
    for l in xrange(len(vx)):
        c = counts[l]
        if c != 0:
            c = c**-1
            vx[l] *= c
            vy[l] *= c
            vz[l] *= c



@guvectorize([(float64[:], # x
               float64[:], # u
               float64[:], # v
               float64[:], # w
               float64[:], # bin_width
               float64[:], # xmin
               float64[:], # px binned
               float64[:], # py binned
               float64[:], # pz binned
               float64[:] # counts
               )],
               '(n),(n),(n),(n),(),(),(m),(m), (m), (m)', nopython = True, cache = True)
def CalcPHists(x, u, v, w, bin_width, xmin, px, py, pz, counts):
    bn = bin_width[0]
    minx =xmin[0]
    maxl = len(px)-1
    for i in xrange(len(x)):
        l = int(x[i]//bn)
        if l > maxl:
            l = maxl
        px[l] += u[i]
        py[l] += v[i]
        pz[l] += w[i]
        counts[l] += 1
    for l in xrange(len(px)):
        c = counts[l]
        if c != 0:
            c = c**-1
            px[l] *= c
            py[l] *= c
            pz[l] *= c

@guvectorize([(float64[:], # x
               float64[:], # u
               float64[:], # v
               float64[:], # w
               float64[:], # weights
               float64[:], # bin_width
               float64[:], # xmin
               float64[:], # px binned
               float64[:], # py binned
               float64[:], # pz binned
               float64[:] # counts
               )],
               '(n),(n),(n),(n),(n),(),(),(m),(m), (m), (m)', nopython = True)
def CalcPWeightedHists(x, u, v, w, weights, bin_width, xmin, px, py, pz, counts):
    bn = bin_width[0]
    minx = xmin[0]
    maxl = len(px)-1
    for i in xrange(len(x)):
        c1 = weights[i]
        l = int(x[i]//bn)
        if l > maxl:
            l = maxl
        px[l] += u[i]*c1
        py[l] += v[i]*c1
        pz[l] += w[i]*c1
        counts[l] += c1
    # normalize
    for l in xrange(len(px)):
        c = counts[l]
        if c != 0:
            c = c**-1
            px[l] *= c
            py[l] *= c
            pz[l] *= c

@guvectorize([(float64[:], float64[:], float64[:], float64[:],float64[:],float64[:])], '(m),(m),(m),(m),(m),(m)', cache = True)
def RestFrameBoost(vx_e, ecounts, vx_i, icounts, vx_avg, boost_g):
    for i in xrange(len(vx_e)):
        vx_avg[i] = (vx_e[i]*ecounts[i]+vx_i[i]*icounts[i])/(icounts[i]+ecounts[i])
        boost_g[i] = 1/sqrt(1+vx_avg[i]*vx_avg[i])

@guvectorize([(float64[:], float64[:], float64[:], float64[:], float64[:])], '(m),(m),(m),(m) ->(m)', nopython =True, cache = True)
def Total(e_arr,  ecounts, i_arr, icounts, ans):
    for i in xrange(len(e_arr)):
        ans[i] = (e_arr[i]*ecounts[i]+i_arr[i]*icounts[i])/(ecounts[i]+icounts[i])

@guvectorize([(float64[:], # x
               float64[:], # u
               float64[:], # v
               float64[:], # w
               float64[:], # g
               float64[:], # vx_avg
               float64[:], # boost_g,
               float64[:], # bin_width
               float64[:], # xmin
               float64[:], # Counts
               float64[:] # T
               )],
               '(n),(n),(n),(n),(n),(m),(m),(),(),(m),(m)', nopython = True, cache = True)
def CalcDelGamHists(x, u, v, w, g, vx_avg, boost_g, bin_width, xmin, counts, T):
    bn = bin_width[0]
    lmax = len(vx_avg)
    minx = xmin[0]
    for i in xrange(len(x)):
        l = int((x[i]-minx)//bn)
        if l >=lmax:
            l = lmax - 1
        c2 = boost_g[l]*(u[i]-g[i]*vx_avg[l]) # boosted
        T[l] += sqrt(c2*c2+v[i]*v[i]+w[i]*w[i]+1)-1
    for l in xrange(len(counts)):
        c = counts[l]
        if c != 0:
            T[l] *= c**-1

@guvectorize([(float64[:], # x
               float64[:], # u
               float64[:], # v
               float64[:], # w
               float64[:], # g
               float64[:], # weights
               float64[:], # vx_avg
               float64[:], # boost_g,
               float64[:], # bin_width
               float64[:], # xmin
               float64[:], # Counts
               float64[:] # T
               )],
               '(n),(n),(n),(n),(n),(n),(m),(m),(),(),(m),(m)', nopython = True, cache = True)
def CalcDelGamWeightedHists(x, u, v, w, g, weights, vx_avg, boost_g, bin_width, xmin, counts, T):
    bn = bin_width[0]
    lmax = len(vx_avg)
    minx = xmin[0]
    for i in xrange(len(x)):
        l = int((x[i]-minx)//bn)
        if l >=lmax:
            l = lmax - 1
        c2 = boost_g[l]*(u[i]-g[i]*vx_avg[l]) # boosted
        T[l] += (sqrt(c2*c2+v[i]*v[i]+w[i]*w[i]+1)-1)*weights[i]
    for l in xrange(len(counts)):
        c = counts[l]
        if c != 0:
            T[l] *= c**-1
