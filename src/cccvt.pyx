#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division

cimport cython
cimport numpy as np
from cpython cimport bool

# from libc.math cimport sqrt

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



@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef dict __get_inv_tessellation(long[:] tessellation):

  cdef dict inv = {}
  cdef long i
  cdef long t

  cdef long num = len(tessellation)

  for i in range(num):

    t = tessellation[i]

    if t in inv:
      inv[t].add(i)
    else:
      inv[t] = set([i])

  return inv

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef void __capacity_randint(long[:] tessellation, long m, long n):

  cdef np.ndarray['long', ndim=1] tessellation_count = zeros(n, 'int')
  cdef long cap = <long>(m/n)

  cdef long i
  cdef long r

  for i in xrange(m):

    while True:

      r = randint(n)
      with nogil:
        if tessellation_count[r]>cap:
          continue
        else:
          tessellation[i] = r
          tessellation_count[r] += 1
          break

  return

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef void __simple_sort(long[:,:] a, double[:,:] b, long n, long ind) nogil:

  # im sorry ...

  cdef long i
  cdef long j
  cdef long maxind

  cdef long tmpa
  cdef double tmpb

  for i in xrange(n-1):

    maxind = i

    for j in xrange(i,n):

      if b[j,ind]>b[maxind,ind]:
        maxind = j

    if maxind != i:

      tmpb = b[i,ind]
      b[i,ind] = b[maxind,ind]
      b[maxind,ind] = tmpb
      tmpa = a[i,ind]
      a[i,ind] = a[maxind,ind]
      a[maxind,ind] = tmpa

  return


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef long __get_h(
  double[:,:] dd,
  long[:,:] Hi,
  double[:,:] Hw,
  dict inv,
  long si,
  long sj
):

  cdef np.ndarray['long', ndim=1] s

  cdef set ii = inv[si]
  cdef set jj = inv[sj]
  cdef long numi = len(ii)
  cdef long numj = len(jj)
  cdef long ml = min(numi,numj)
  cdef long k
  cdef long x

  k = 0
  for x in ii:
    with nogil:
      Hw[k,0] = dd[x,si] - dd[x,sj]
      Hi[k,0] = x
      k += 1

  k = 0
  for x in jj:
    with nogil:
      Hw[k,1] = dd[x,sj] - dd[x,si]
      Hi[k,1] = x
      k += 1

  __simple_sort(Hi,Hw,numi,0)
  __simple_sort(Hi,Hw,numj,1)

  return min(numi,numj)

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef void __remap(long[:] tes, dict inv, long xi, long xj, long si, long sj):
  tes[xi] = sj
  tes[xj] = si
  inv[si].remove(xi)
  inv[si].add(xj)
  inv[sj].remove(xj)
  inv[sj].add(xi)
  return

cpdef Ccvt(
  np.ndarray['double', ndim=2] domain,
  np.ndarray['double', ndim=2] sites,
  double tol=1.e-7,
  long maxitt=10000
):

  cdef long m = len(domain)
  cdef long n = len(sites)
  cdef double cap = m/n

  cdef np.ndarray['double', ndim=2] sites_prev = zeros((n,2),'float')
  sites_prev[:,:] = sites[:,:]

  print('point cloud')
  print('points (m): {:d}, sites (n): {:d}, cap: {:f}'.format(m,n,cap))

  cdef double now = time()

  cdef long k
  cdef long i
  cdef long h
  cdef long si
  cdef long sj
  cdef long p
  cdef set points
  cdef bool stable

  cdef long hsize = <long>(6*cap)

  cdef np.ndarray['long', ndim=2] Hi = zeros((hsize,2),'int')
  cdef np.ndarray['long', ndim=1] tessellation = zeros(m,'int')
  cdef np.ndarray['double', ndim=2] Hw = zeros((hsize,2),'float')
  cdef double cap_err
  cdef double c

  for k in xrange(maxitt):

    dd = cdist(domain, sites, 'euclidean')
    __capacity_randint(tessellation,m,n) #  x → s
    inv_tessellation = __get_inv_tessellation(tessellation) # s → x

    i = -1
    while True:

      i += 1

      print('itt: ', k, i)

      stable = True

      for si in xrange(n):
        for sj in xrange(n):

          ml = __get_h(
            dd,
            Hi,
            Hw,
            inv_tessellation,
            si,
            sj
          )

          h = 0

          while h<ml:

            if Hw[h,0] + Hw[h,1]<=0:
              break

            __remap(tessellation,inv_tessellation,Hi[h,0],Hi[h,1],si,sj)
            stable = False
            h += 1

      if stable:
        break

    for k in inv_tessellation:
      points = inv_tessellation[k]
      if points:

        c = 0.0
        sites[k,0] = 0.0
        sites[k,1] = 0.0
        for p in points:
          sites[k,0] += domain[p,0]
          sites[k,1] += domain[p,1]
          c += 1.

        if c>0.0:
          sites[k,0] /= c
          sites[k,1] /= c

    for k in xrange(n):
      cap_err += square(<double>len(inv_tessellation[k])/cap-1.0)
    cap_err /= <double>n

    diff_err = norm(sites - sites_prev)/<double>n
    sites_prev[:] = sites[:]

    if abs(diff_err)<tol:
      print('terminating, reached diff error: {:0.8f}'.format(diff_err))
      print('capacity error: {:0.8f}'.format(cap_err))
      break
    else:
      print('diff error {:0.8f}, going again ...'.format(diff_err))
      print('capacity error: {:0.8f}'.format(cap_err))

  print('time: {:0.8f}'.format(time()-now))

  return sites, {k:points for k,points in inv_tessellation.iteritems() if points}

