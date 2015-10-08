#!/usr/bin/python

import msympy, argparse



def read_xyz(fin):
    length = int(fin.readline())
    comment = fin.readline()[:-1]
    elements = []
    for i in range(0,length):
        line = fin.readline().split()
        elements.append(msympy.Element(name = line[0], coordinates = map(float, line[1:4])))

    return (elements, comment)

def write_xyz(fout, elements, title):
    fout.write("%d\n%s\n" % (len(elements), title))
    for e in elements:
        v = e.coordinates
        fout.write("%s %.15g %.15g %.15g\n" % (e.name, v[0], v[1], v[2]))

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=argparse.FileType('r'))
parser.add_argument('outfile', type=argparse.FileType('w'))
args = parser.parse_args()

(elements, comment) = read_xyz(args.infile)
ctx = msympy.Context()
ctx.elements = elements
point_group_name = ctx.findSymmetry()
elements = ctx.symmetrizeElements()
write_xyz(args.outfile, elements, comment + " symmetrized by libmsympy according to point group " + point_group_name)




