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
 
def set_basis(element):
    basis_function = msym.RealSphericalHarmonic(element = element, name = "1s")
    element.basis_functions = [basis_function]
    return basis_function

basis_functions = [set_basis(e) for e in elements]

#with msym.Context() as ctx:
#    ctx.elements = elements
#    point_group = ctx.find_symmetry()
#    selements = ctx.symmetrize_elements()
#    write_xyz(args.outfile, selements, comment + " symmetrized by libmsym according to point group " + point_group)
 
with msym.Context(elements = elements, basis_functions = basis_functions) as ctx:
    point_group = ctx.find_symmetry()
    selements = ctx.symmetrize_elements()
    write_xyz(args.outfile, selements, comment + " symmetrized by libmsym according to point group " + point_group)
    ctx.generate_salc_subspaces()

#with msym.Context(elements = elements, point_group = "T") as ctx:
#    point_group = ctx.point_group
#    selements = ctx.symmetrize_elements()
#    write_xyz(args.outfile, selements, comment + " symmetrized by libmsym according to point group " + point_group)




