from __future__ import print_function, division
from numpy import zeros, nan
from numba import jit, guvectorize, float64, int32
from math import sqrt, log10
#import os
#os.environ['NUMBA_WARNINGS'] = '1'
@guvectorize([(float64[:], float64[:],float64[:])], '(),()->()', nopython = True, target='parallel')
def vecLog10Norm(xArr,norm, ans):
    if xArr[0]>0:
        ans[0] = log10(xArr[0]/norm[0])
    else:
        ans[0]=nan

@jit(nopython=True)#
def Fast2DHist(x1, x2, min1, max1, bnum1, min2,max2, bnum2):
    hist= zeros((bnum1, bnum2))
    b1_w = ((max1-min1)/bnum1)**-1
    b2_w =((max2-min2)/bnum2)**-1
    for i in range(len(x1)):
        if x1[i]>=min1:
            if x1[i] == max1:
                j = bnum1-1
            else:
                j = (x1[i]-min1)*b1_w
            if j<bnum1:
                if x2[i]>=min2:
                    if x2[i]==max2:
                        k =bnum2-1
                    else:
                        k = (x2[i]-min2)*b2_w
                    if k<bnum2:
                        hist[int(j),int(k)] += 1
    return hist
@jit(nopython=True)#
def Fast2DWeightedHist(x1, x2, weights, min1, max1, bnum1, min2,max2, bnum2):
    hist= zeros((bnum1, bnum2))
    b1_w = ((max1-min1)/bnum1)**-1
    b2_w =((max2-min2)/bnum2)**-1
    for i in range(len(x1)):
        if x1[i]>=min1:
            j = int((x1[i]-min1)*b1_w)
            if j<bnum1:
                if x2[i]>=min2:
                    k = int((x2[i]-min2)*b2_w)
                    if k<bnum2:
                        hist[j,k] += weights[i]
    return hist


if __name__ =='__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    import data_loading

    filepath = 'output/prtl.tot.003'
    x = data_loading.load_dataset(filepath, 'xi', slice(None))
    px = data_loading.load_dataset(filepath, 'ui', slice(None))

    #for k in range(10):
    #    x = np.append(x, x)
    #    px = np.append(px, px)
    #print(x.min())
    hist = Fast2DHist(px, x, px.min(), px.max(),200, x.min(), x.max(), 200)
    numpyHist = np.histogram2d(px, x,
        bins = [200, 200],
        range = [[px.min(),px.max()],[x.min(),x.max()]])[0]
    #print(numpyHist)
    plt.subplot(211)
    plt.imshow(hist)
    plt.subplot(212)
    plt.imshow(numpyHist)
    plt.show()
    print(np.all(hist ==numpyHist))
