import numpy as np
from numba import jit, guvectorize, float64, int32
from math import sqrt
#import os
#os.environ['NUMBA_WARNINGS'] = '1'
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

@guvectorize([(float64[:], float64[:], float64[:],float64[:])], '(),(),(),()', nopython = True, cache = True, target='parallel')
def LorentzFactor(u,v,w, ans):
    # Iterate over the array, find the min and max and divide by a number
    # Calculate xmin and xmax
    ans[0] = u[0]*u[0]+v[0]*v[0]+ w[0]*w[0] + 1
    ans[0] = sqrt(ans[0])


@jit(nopython = True, cache = True)#, target='parallel')
def CalcVxEHists(x, u, g, bin_width, xmin, vx, E, counts):
    maxl = len(vx)
    for i in xrange(len(x)):
        l = int((x[i]-xmin)//bin_width)
        if l>=0 and l < maxl:
            vx[l] += u[i]*g[i]**-1
            E[l] += g[i]-1
            counts[l] += 1
    for l in xrange(len(vx)):
        if counts[l] != 0:
            c = counts[l]**-1
            vx[l] *= c
            E[l] *= c

@jit(nopython = True, cache = True)
def CalcVxEWeightedHists(x, u, g, weights, bin_width, xmin, vx,  E, counts):
    maxl = len(vx)
    for i in xrange(len(x)):
        c1 = weights[i]
        l = int((x[i]-xmin)//bin_width)
        if 0<= l and l < maxl:
            vx[l] += u[i]*g[i]**-1*c1
            E[l] += (g[i]-1)*c1
            counts[l] += c1
    for l in xrange(len(vx)):
        if counts[l] != 0:
            c = counts[l]**-1
            vx[l] *= c
            E[l] *= c

@jit(nopython = True, cache = True)
def CalcVHists(x, u, v, w, g, bin_width, xmin, vx, vy, vz, counts):
    maxl = len(vx)
    for i in xrange(len(x)):
        l = int((x[i]-xmin)//bin_width)
        if 0<=l and l < maxl:
            tmp = g[i]**-1
            vx[l] += u[i]*tmp
            vy[l] += v[i]*tmp
            vz[l] += w[i]*tmp
            counts[l] += 1
    # Normalize
    for l in xrange(len(vx)):
        if counts[l] != 0:
            c = counts[l]**-1
            vx[l] *= c
            vy[l] *= c
            vz[l] *= c

@jit(nopython = True, cache = True)
def CalcVWeightedHists(x, u, v, w, g, weights, bin_width, xmin, vx, vy, vz, counts):
    maxl = len(vx)
    for i in xrange(len(x)):
        l = int((x[i]-xmin)//bin_width)
        if 0<= l and l < maxl:
            tmp = g[i]**-1*weights[i]
            vx[l] += u[i]*tmp
            vy[l] += v[i]*tmp
            vz[l] += w[i]*tmp
            counts[l] += weights[i]
    # Normalize
    for l in xrange(len(vx)):
        c = counts[l]
        if c != 0:
            c = c**-1
            vx[l] *= c
            vy[l] *= c
            vz[l] *= c



@jit(nopython = True, cache = True)
def CalcPHists(x, u, v, w, bin_width, xmin, px, py, pz, counts):
    maxl = len(px)
    for i in xrange(len(x)):
        l = int((x[i]-xmin)//bin_width)
        if 0<= l  and l < maxl:
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
@jit(nopython = True, cache = True)
def CalcPWeightedHists(x, u, v, w, weights, bin_width, xmin, px, py, pz, counts):
    maxl = len(px)
    for i in xrange(len(x)):
        c1 = weights[i]
        l = int((x[i]-xmin)//bin_width)
        if 0<=l and l < maxl:
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

@guvectorize([(float64[:], float64[:], float64[:], float64[:],float64[:],float64[:])], '(),(),(),(),(),()', cache = True,target='parallel')
def RestFrameBoost(vx_e, ecounts, vx_i, icounts, vx_avg, boost_g):
    if ecounts[0] != 0 or icounts[0] != 0:
        vx_avg[0] = (vx_e[0]*ecounts[0]+vx_i[0]*icounts[0])/(icounts[0]+ecounts[0])
        boost_g[0] = 1/sqrt(1-vx_avg[0]*vx_avg[0])

@guvectorize([(float64[:], float64[:], float64[:], float64[:], float64[:])], '(),(),(),() ->()', nopython =True, cache = True,target='parallel')
def Total(e_arr,  ecounts, i_arr, icounts, ans):
    if ecounts[0] != 0 or icounts[0] != 0:
        ans[0] = (e_arr[0]*ecounts[0]+i_arr[0]*icounts[0])/(ecounts[0]+icounts[0])

@jit(nopython = True, cache = True)
def CalcDelGamHists(x, u, v, w, g, vx_avg, boost_g, bin_width, xmin, counts, T):
    lmax = len(vx_avg)
    for i in xrange(len(x)):
        l = int((x[i]-xmin)//bin_width)
        if 0 <= l and l <lmax:
            c2 = boost_g[l]*(u[i]-g[i]*vx_avg[l]) # boosted
            T[l] += sqrt(c2*c2+v[i]*v[i]+w[i]*w[i]+1)-1
    for l in xrange(len(counts)):
        c = counts[l]
        if c != 0:
            T[l] *= c**-1

@jit(nopython = True, cache = True)
def CalcDelGamWeightedHists(x, u, v, w, g, weights, vx_avg, boost_g, bin_width, xmin, counts, T):
    lmax = len(vx_avg)
    for i in xrange(len(x)):
        l = int((x[i]-xmin)//bin_width)
        if 0<= l  and l < lmax:
            c2 = boost_g[l]*(u[i]-g[i]*vx_avg[l]) # boosted
            T[l] += (sqrt(c2*c2+v[i]*v[i]+w[i]*w[i]+1)-1)*weights[i]
    for l in xrange(len(counts)):
        c = counts[l]
        if c != 0:
            T[l] *= c**-1
