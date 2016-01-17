#!python
# -*- coding: utf-8 -*-
#cython: language_level=3, boundscheck=False, wraparound=False
#cython: nonecheck=False, cdivision=True

from __future__ import print_function
from __future__ import division

from libc.stdlib cimport malloc, free

cimport cython
cimport numpy as np
from cpython cimport bool

from libc.math cimport sqrt
from libc.math cimport pow

from time import time

from numpy.random import randint
from numpy.linalg import norm

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
cpdef double __capacity_randint(long[:] tes, long m, long n):

  cdef long i
  cdef long r

  cdef long* count = <long*>malloc(sizeof(long)*n)
  cdef double cap = <double>m/<double>n
  cdef long icap = <long>cap

  if m % n != 0:
    icap += 1

  for i in xrange(n):
    count[i] = 0

  print('cap', cap)
  print('icap', icap)

  for i in xrange(m):

    while True:

      r = randint(n)
      with nogil:
        if count[r]>icap-1:
          continue
        else:
          tes[i] = r
          count[r] += 1
          break

  free(count)

  return cap

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef double __get_cap_error(dict inv, long n, double cap):

  cdef long k
  cdef set s
  cdef double err = 0.0

  for k in xrange(n):
    s = inv[k]
    err += square(<double>len(s)/cap-1.0)

  return err/<double>n

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
cdef void __distance(double* dd, double[:,:] domain, double[:,:] sites, long m, long n) nogil:

  cdef long i
  cdef long j

  for i in xrange(m):
    for j in xrange(n):
      dd[i*n+j] = sqrt(pow(domain[i,0] - sites[j,0],2) +
                       pow(domain[i,1] - sites[j,1],2))

  return


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef long __get_h(
  double* dd,
  long m,
  long n,
  long[:,:] Hi,
  double[:,:] Hw,
  long[:] ii,
  long[:] jj,
  long numi,
  long numj,
  long si,
  long sj
) nogil:

  cdef long ml = min(numi,numj)
  cdef long k
  cdef long x

  cdef double tmp1
  cdef double tmp2

  for k in xrange(numi):
    x = ii[k]
    with nogil:
      Hi[k,0] = x
      Hw[k,0] = dd[x*n+si] - dd[x*n+sj]

  for k in xrange(numj):
    x = jj[k]
    with nogil:
      Hi[k,1] = x
      Hw[k,1] = dd[x*n+sj] - dd[x*n+si]

  __simple_sort(Hi,Hw,numi,0)
  __simple_sort(Hi,Hw,numj,1)

  return ml

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

  cdef long k
  cdef long i
  cdef long h
  cdef long si
  cdef long sj
  cdef long p
  cdef long m = len(domain)
  cdef long n = len(sites)
  cdef double cap = <double>m/<double>n
  print('points (m): {:d}, sites (n): {:d}, cap: {:f}'.format(m,n,cap))

  cdef double* dd = <double*>malloc(sizeof(double)*m*n)
  # for i in xrange(m):
    # dd[i] = <double*>malloc(sizeof(double)*n)

  cdef double now = time()

  cdef set points
  cdef bool stable

  cdef double cap_err
  cdef double count

  cdef np.ndarray['double', ndim=2] sites_prev = zeros((n,2),'float')
  for i in xrange(n):
    sites_prev[i,0] = sites[i,0]
    sites_prev[i,1] = sites[i,1]

  cdef long hsize = <long>(6*cap)
  cdef np.ndarray['long', ndim=2] Hi = zeros((hsize,2),'int')
  cdef np.ndarray['double', ndim=2] Hw = zeros((hsize,2),'float')
  cdef np.ndarray['long', ndim=1] tessellation = zeros(m,'int')

  for k in xrange(maxitt):

    print('distance')
    __distance(dd, domain, sites, m, n)
    __capacity_randint(tessellation,m,n) #  x → s
    inv_tessellation = __get_inv_tessellation(tessellation) # s → x

    cap_err = __get_cap_error(inv_tessellation, n, cap)
    print('capacity error: {:0.8f}'.format(cap_err))

    i = -1
    while True:

      i += 1

      print('itt: ', k, i)

      stable = True

      for si in xrange(n):
        for sj in xrange(n):

          ml = __get_h(
            dd,
            m,
            n,
            Hi,
            Hw,
            array(list(inv_tessellation[si]),'int'),
            array(list(inv_tessellation[sj]),'int'),
            len(inv_tessellation[si]),
            len(inv_tessellation[sj]),
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

        count = 0.0
        sites[k,0] = 0.0
        sites[k,1] = 0.0
        for p in points:
          sites[k,0] += domain[p,0]
          sites[k,1] += domain[p,1]
          count += 1.0

        if count>0.0:
          sites[k,0] /= count
          sites[k,1] /= count


    cap_err = __get_cap_error(inv_tessellation, n, cap)
    diff_err = norm(sites - sites_prev)/<double>n
    for i in xrange(n):
      sites_prev[i,0] = sites[i,0]
      sites_prev[i,1] = sites[i,1]

    print('capacity error: {:0.8f}'.format(cap_err))
    print('diff error {:0.8f} ... '.format(diff_err))

    if abs(diff_err)<tol:
      print('terminating')
      break
    else:
      print('going again')


  # for i in xrange(m):
    # free(dd[i])
  free(dd)

  print('time: {:0.8f}'.format(time()-now))

  return sites, {k:points for k,points in inv_tessellation.iteritems() if points}

