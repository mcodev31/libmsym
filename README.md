# libmsym

libmsym is a C library dealing with point group symmetry in molecules. It can determine, symmetrize and generate molecules of any point group. It can also generate symmetry adapted linear combinations of atomic orbitals for a subset of all point groups and orbital angular momentum (l), and project orbitals into the irreducible representation with the larges component.

## Performance

120 (Ih) symmetry operations in C-720 fullerene found in ~30ms (2012 macbook air).
Detection, generation of permutation information and symmetrisation of above in ~70ms.
T point group protein with 15k+ elements detected and symmetrized in ~6s.
300 symmetry adapted orbitals of minimal basis C-60 buckminster fullerene generated in 0.18s

## Additional information

This was developed as a masters project. Thesis with descriptions of algorithms will be made available soon.
