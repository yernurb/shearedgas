import cairo
import numpy as np 

ims = cairo.ImageSurface(cairo.FORMAT_ARGB32, 400, 400)
ctx = cairo.Context(ims)

ctx.set_source_rgb(0.2, 0.2, 0.2)
ctx.fill()

ctx.set_source_rgb(0.7, 0.2, 0.0)
ctx.translate(200, 200)
ctx.arc(0, 0, 50, 0, 2*np.pi)
ctx.stroke_preserve()
ctx.fill()

ctx.set_source_rgb(0.1, 0.2, 0.7)
ctx.arc(20, 10, 10, 0, 2*np.pi)
ctx.stroke_preserve()
ctx.fill()


ims.write_to_png("example.png")
