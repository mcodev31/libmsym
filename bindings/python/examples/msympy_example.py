#!/usr/bin/python

import libmsym as msym, numpy as np, argparse, random, sys



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


#msym.init(library_location='/path_to_libmsym_library/libmsym.dylib') # e.g. for OS X without libmsym in dload path

#with msym.Context() as ctx:
#    ctx.elements = elements
#    point_group = ctx.find_symmetry()
#    selements = ctx.symmetrize_elements()
#    write_xyz(args.outfile, selements, comment + " symmetrized by libmsym according to point group " + point_group)
 
with msym.Context(elements = elements, basis_functions = basis_functions) as ctx:
    point_group = ctx.find_symmetry()
    selements = ctx.symmetrize_elements()
    write_xyz(args.outfile, selements, comment + " symmetrized by libmsym according to point group " + point_group)
    ctx.symmetry_operations #symmetry operations
    ctx.subrepresentation_spaces # subspace
    ctx.subrepresentation_spaces[0].symmetry_species #symmetry species of space (index into character table)
    ctx.subrepresentation_spaces[0].salcs # salcs that span space
    ctx.subrepresentation_spaces[0].salcs[0].basis_functions # basis functions for salc    
    ctx.subrepresentation_spaces[0].salcs[0].partner_functions # numpy array of partner functions expressed in terms of basis_functions coefficients

    ctx.character_table.table # table as numpy array
    ctx.character_table.symmetry_operations # representative symmetry operations
    ctx.character_table.symmetry_species # symmetry species
    ctx.character_table.symmetry_species[0].dim  # dimensionality of symmetry species
    ctx.character_table.symmetry_species[0].name # name of symmetry species e.g. A2g

    somefunc = np.zeros((len(basis_functions)),dtype=np.float64)
    for i in range(0,len(somefunc)):
        somefunc[i] = i

    species_components = ctx.symmetry_species_components(somefunc)
    
    species = ctx.character_table.symmetry_species

    print(somefunc)
    for i, c in enumerate(species_components):
        print(str(c) + species[i].name)
    
    #matrix version of the above, as well as wave function symmetrization
    (matrix, species, partners) = ctx.salcs
    (d,d) = matrix.shape
    print(matrix)
    print(species)
    print([(p.index, p.dim) for p in partners])
    indexes = [x for x in range(0,d)]
    random.shuffle(indexes)
    matrix = matrix[indexes,:]
    matrix += 0.01
    print(matrix)
    (_same_matrix, species,partners) = ctx.symmetrize_wavefunctions(matrix)
    print(matrix)
    print(species)
    print([(p.index, p.dim) for p in partners])

    #generating elements
    ctx.point_group = "D6h"
    gen_elements = [msym.Element(name = "C", coordinates = [1.443524, 0.0,0.0]), msym.Element(name = "H", coordinates = [2.568381, 0.0, 0.0])]
    benzene = ctx.generate_elements(gen_elements)
    maxcomp = max([max(e.coordinates) for e in benzene])
    print(len(benzene),"\nbenzene")
    for e in benzene:
        vec = np.asarray(e.coordinates)
        vec[vec < maxcomp*sys.float_info.epsilon] = 0
        print(e.name, vec[0],vec[1],vec[2])
    

#with msym.Context(elements = elements, point_group = "T") as ctx:
#    point_group = ctx.point_group
#    selements = ctx.symmetrize_elements()
#    write_xyz(args.outfile, selements, comment + " symmetrized by libmsym according to point group " + point_group)




