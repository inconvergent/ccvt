#!/usr/bin/python
# -*- coding: utf-8 -*-


from __future__ import print_function
from __future__ import division

from collections import defaultdict
from itertools import product
from itertools import repeat
from time import time

from numpy.random import randint
from numpy.linalg import norm

from scipy.spatial.distance import cdist
from numpy import mean
from numpy import ceil
from numpy import square
from numpy import column_stack
from numpy import row_stack
from numpy import zeros
from numpy import array
from numpy import reshape
from numpy import argsort



def __get_inv_tessellation(tessellation):

  inv = defaultdict(set)
  for i,t in enumerate(tessellation):
    inv[t].add(i)

  return inv

def __capacity_randint(m, n):

  tessellation = zeros(m, 'int')
  tessellation_count = {k:0 for k in xrange(n)}
  cap = m/n

  for i in xrange(m):

    while True:

      r = randint(n)
      if tessellation_count[r]>=cap:
        continue
      else:
        tessellation[i] = r
        tessellation_count[r] += 1
        break

  return tessellation

def __get_h(dd, xii, si, sj):

  #TODO: heap. this is too slow.

  # xii = (tessellation == si).nonzero()[0]

  hw = dd[xii,si] - dd[xii,sj]
  hi = array(xii, 'int')
  s = argsort(hw, kind='mergesort')

  return hi[s],hw[s]

def __remap(tes, inv, xi,xj,si,sj):
  tes[xi] = sj
  tes[xj] = si
  inv[si].remove(xi)
  inv[si].add(xj)
  inv[sj].remove(xj)
  inv[sj].add(xi)
  return


def Ccvt(
  domain,
  sites,
  tol=1.e-7,
  maxitt=10000
):

  m = len(domain)
  n = len(sites)
  cap = m/n

  sites_prev = zeros((n,2),'float')
  sites_prev[:] = sites[:]

  print('point cloud')
  print('points (m): {:d}, sites (n): {:d}, cap: {:f}'.format(m,n,cap))


  now = time()

  for k in xrange(maxitt):

    dd = cdist(domain, sites, 'euclidean')
    tessellation = __capacity_randint(m,n) #  x → s
    inv_tessellation = __get_inv_tessellation(tessellation) # s → x

    i = -1
    while True:

      i += 1

      print('itt: ', k, i)

      stable = True

      for si,sj in product(range(n), repeat=2):

        Hii, Hiw = __get_h(dd, list(inv_tessellation[si]), si, sj)
        Hjj, Hjw = __get_h(dd, list(inv_tessellation[sj]), sj, si)

        ml = max(len(Hii),len(Hjj))
        h = -1

        while h>-ml-1:

          if Hiw[h] + Hjw[h]<=0:
            break

          __remap(tessellation,inv_tessellation,Hii[h],Hjj[h],si,sj)

          stable = False
          h -= 1

      if stable:
        break

    agg = [[] for i in repeat(None, n)]
    for t,xy in zip(tessellation, domain):
      agg[t].append(xy)

    for k, v in enumerate(agg):
      if v:
        sites[k,:] = mean(v, axis=0)

    cap_count = array([len(v) for v in inv_tessellation.values()],'float')
    cap_err = square(cap_count/float(cap)-1.0).sum()/n
    diff_err = norm(sites - sites_prev)/n
    sites_prev[:] = sites[:]

    if abs(diff_err)<tol:
      print('terminating, reached diff error: {:0.5f}'.format(diff_err))
      print('capacity error: {:0.5f}'.format(cap_err))
      break
    else:
      print('diff error {:0.5f}, going again ...'.format(diff_err))
      print('capacity error: {:0.5f}'.format(cap_err))

  print('time: {:0.5f}'.format(time()-now))

  return sites, {k:v for k, v in inv_tessellation.iteritems() if v}

