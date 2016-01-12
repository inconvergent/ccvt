#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function


from numpy import linspace
from numpy import meshgrid
from numpy import floor
from numpy import zeros
from numpy.random import random
from scipy.ndimage import imread


def grid(m, n):

  x = linspace(0, 1, m)
  y = linspace(0, 1, n)

  return meshgrid(x, y)

def get_dens_from_img(fn):

  return 1.0-imread(fn)/255.

def get_dens_example(n):

  print('making {:d} × {:d} grid, {:d} elements'.format(m,n,m*n))

  x,y = grid(n,n)

  return x*y

def sample_from_dens(dens, n):

  m = dens.shape[0]
  res = zeros((n,2),'float')
  k = 0

  while k<n:

    xy = random(2)
    ij = floor(xy*m)
    d = dens[ij[0],ij[1]]
    if random()<d:
      res[k,:] = xy
      k += 1

  return res

