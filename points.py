#!/usr/bin/python
# -*- coding: utf-8 -*-


from __future__ import division
from __future__ import print_function

# from numpy.random import seed
# seed(3)


BACK = [1,1,1,1]
FRONT = [0,0.7,0.7,0.7]
BLACK = [0,0,0,0.9]
LIGHT = [0,0,0,0.3]

SIZE = 1000
ONE = 1.0/SIZE

RAD = 0.45


def main():

  from numpy.random import random
  from numpy import zeros
  from dddUtils.pointCloud import point_cloud
  from dddUtils.ioOBJ import export_2d as export

  from modules.utils import get_dens_example
  from modules.utils import get_dens_from_img
  from modules.utils import sample_from_dens

  from ccvt import Ccvt as ccvt

  fn = './data/kelp.png'
  n = 10000
  m = 100000

  print('get density')
  dens = get_dens_from_img(fn)
  # dens = get_dens_example(100)

  print('sample domain')
  domain = sample_from_dens(dens, m)
  print('sample dens')
  org_sites = sample_from_dens(dens, n)

  sites, inv_tesselation = ccvt(domain, org_sites, maxitt=5)
  export('voronoi','./res/out.2obj', sites)



if __name__ == '__main__':

  if False:
    import pstats
    import cProfile
    pfilename = './profile/profile'
    cProfile.run('main()',pfilename)
    # p = pstats.Stats(pfilename)
    # p.strip_dirs().sort_stats('cumulative').print_stats()
  else:
    main()

