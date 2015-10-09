#!/usr/bin/python

import libmsym as msym, argparse



def read_xyz(fin):
    length = int(fin.readline())
    comment = fin.readline()[:-1]
    elements = []
    for i in range(0,length):
        line = fin.readline().split()
        elements.append(msym.Element(name = line[0], coordinates = map(float, line[1:4])))

    return (elements, comment)

def write_xyz(fout, elements, comment):
    fout.write("%d\n%s\n" % (len(elements), comment))
    for e in elements:
        v = e.coordinates
        fout.write("%s %.14f %.14f %.14f\n" % (e.name, *v))

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=argparse.FileType('r'))
parser.add_argument('outfile', type=argparse.FileType('w'))
args = parser.parse_args()

(elements, comment) = read_xyz(args.infile)

with msym.Context() as ctx:
    ctx.elements = elements
    point_group = ctx.findSymmetry()
    selements = ctx.symmetrizeElements()
    write_xyz(args.outfile, selements, comment + " symmetrized by libmsympy according to point group " + point_group)




