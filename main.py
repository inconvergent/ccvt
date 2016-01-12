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
  from dddUtils.ioOBJ import export_2d as export

  from ccvt.utils import get_dens_example
  from ccvt.utils import get_dens_from_img
  from ccvt.utils import sample_from_dens
  from ccvt.ccvt import Ccvt

  fn = './data/mountain2.png'
  n = 200
  m = 4000

  dens = get_dens_from_img(fn)
  # dens = get_dens_example(100)

  domain = sample_from_dens(dens, m)
  org_sites = sample_from_dens(dens, n)

  sites, inv_tesselation = Ccvt(domain, org_sites, maxitt=4)
  export('voronoi','./res/exit.2obj', sites)

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
  show(render)
  render.write_to_png('./res/exit.png')

  def __write_svg_and_exit(*args):

    gtk.main_quit(*args)
    show(render)
  render.window.connect("destroy", __write_svg_and_exit)

  # gtk.main()



if __name__ == '__main__':

  if True:
    import pstats
    import cProfile
    pfilename = './profile/profile'
    cProfile.run('main()',pfilename)
    # p = pstats.Stats(pfilename)
    # p.strip_dirs().sort_stats('cumulative').print_stats()
  else:
    main()

