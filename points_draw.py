#!/usr/bin/python
# -*- coding: utf-8 -*-


from __future__ import division
from __future__ import print_function


BACK = [1,1,1,1]
FRONT = [0,0.7,0.7,0.7]
BLACK = [0,0,0,0.9]
LIGHT = [0,0,0,0.3]

SIZE = 1000
ONE = 1.0/SIZE

RAD = 0.45


def main():

  import gtk
  from render.render import Animate
  from numpy.random import random
  from numpy import zeros
  from dddUtils.pointCloud import point_cloud

  from modules.utils import get_dens_example
  from modules.utils import get_dens_from_img
  from modules.utils import sample_from_dens

  from ccvt import Ccvt as ccvt

  fn = './data/kelp.png'
  n = 100
  m = 1000

  print('get density')
  dens = get_dens_from_img(fn)
  # dens = get_dens_example(100)

  print('sample domain')
  domain = sample_from_dens(dens, m)
  print('sample dens')
  org_sites = sample_from_dens(dens, n)

  sites, inv_tesselation = ccvt(domain, org_sites, maxitt=5)

  def show(render):
    render.clear_canvas()

    render.set_front(LIGHT)
    for i, s in enumerate(domain):
      render.circle(*s, r=ONE, fill=True)

    render.set_front(BLACK)
    for s, sxy in enumerate(sites):
      render.circle(*sxy, r=3*ONE, fill=True)

    for s,xx in inv_tesselation.iteritems():

      sx, sy = sites[s]

      render.set_front(FRONT)
      for x in xx:
        render.line(sx, sy, domain[x,0], domain[x,1])

      render.set_front(BLACK)
      render.line(sx, sy, *org_sites[s,:])

    # render.set_front(BLACK)
    # for i, s in enumerate(org_sites):
      # render.circle(*s, r=3*ONE, fill=False)

  def wrap(render):
    show(render)
    return False

  render = Animate(SIZE, BACK, FRONT, wrap)
  render.set_line_width(ONE)

  gtk.main()



if __name__ == '__main__':
  main()

