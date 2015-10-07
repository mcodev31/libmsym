#!/usr/bin/python

import msympy

elements = [msympy.Element(mass = 1.0, charge = 1, name = "H", coordinates = [1, 0, 0]),
		msympy.Element(mass = 1.0, charge = 1, name = "H", coordinates = [0, 1, 0]),
		msympy.Element(mass = 1.0, charge = 1, name = "H", coordinates = [0, 0, 1.0000001])]

ctx = msympy.Context()

ctx.elements = elements

name = ctx.findSymmetry()

print "Python has received point group %s." % name

symelements = ctx.symmetrizeElements()
for element in symelements:
	print element.coordinates

elements[0].name = "C"
elements[1].name = "C"
elements[0].coordinates = (0.0,0.0,1.0)
elements[1].coordinates = (0.0,0.0,-1.0)

ctx.elements = [elements[0],elements[1]]

name = ctx.findSymmetry()

print "Python has received point group %s." % name

