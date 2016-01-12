#!/usr/bin/python
# -*- coding: utf-8 -*-


from __future__ import division
from __future__ import print_function

BACK = [1,1,1,1]
FRONT = [0,0,0,1]


SIZE = 1000
ONE = 1.0/SIZE
RAD = 0.45


def draw(render, vertices):

  render.ctx.set_source_rgba(*FRONT)
  render.ctx.set_line_width(ONE)


  for vv in vertices:
    render.circle(vv[0], vv[1], 2*ONE, fill=True)



def main(args):

  from dddUtils.ioOBJ import load_2d as load
  from render.render import Render

  fn = args.fn

  data = load(fn)
  vertices = data['vertices']

  render = Render(SIZE, BACK, FRONT)
  draw(render, vertices)

  out = fn + '.png'
  render.write_to_png(out)



if __name__ == '__main__':

  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument(
    '--fn',
    type=str,
    required=True
  )
  parser.add_argument(
    '--size',
    type=int,
    default=SIZE
  )

  args = parser.parse_args()

  main(args)

