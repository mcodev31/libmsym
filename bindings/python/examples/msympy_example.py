#!/usr/bin/python

import msympy

e1 = msympy.Element(mass = 1.0, charge = 1, name = "H", coordinates = [1, 0, 0])
e2 = msympy.Element(mass = 1.0, charge = 1, name = "H", coordinates = [0, 1, 0])
e3 = msympy.Element(mass = 1.0, charge = 1, name = "H", coordinates = [0, 0, 1])

elements = [e1, e2, e3]

ctx = msympy.Context()

ctx.elements = elements

name = ctx.findSymmetry()

print "Python has received point group %s." % name

e1.name = "C"
e2.name = "C"
e1.coordinates = (0.0,0.0,1.0)
e2.coordinates = (0.0,0.0,-1.0)

ctx.elements = [e1,e2]

name = ctx.findSymmetry()

print "Python has received point group %s." % name
