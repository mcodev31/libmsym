# libmsym v0.2.0
libmsym is a C library dealing with point group symmetry in molecules.

## molecules
Determine, symmetrize and generate molecules of any point group as well determine/select subgroups.
Get symmetrically equivalent elements and symmetrize translation of any any element in a molecule.

## character tables
Generate character tables for any point group (complex characters form reducible representations)

## wave functions
Generate SALCs of real spherical harmonics for non icosahedral point groups with real characters (Ci, Cs, Cnv, Dn, Dnh, Dnd, Td, O and Oh) for any angular momentum, as well as symmetrize orbitals, determine partner functions etc.

## installing

```shell
git clone https://github.com/mcodev31/libmsym.git
cd libmsym
mkdir build
cd build
# build as shared library; build examples (built in ./examples,  not installed)
cmake -DBUILD_SHARED_LIBS:BOOL=ON -DMSYM_BUILD_EXAMPLES:BOOL=ON ../.
make
# sudo only required if installing in directory not owned by user
# use -DCMAKE_INSTALL_PREFIX:PATH=<libmsym installation path> to change
sudo make install
# run examples
./examples/msym_example <input xyz-file>
./examples/msym_tex D13h D13h.tex
```

### python

assumes python 3 is installed

```shell
cd ../bindings/python
python3 setup.py build
# install in user directory
python3 setup.py install --user
# run example
python3 ./examples/msympy_example.py <input xyz-file> <output xyz-file>
```
python requires that libmsym is built as a shared library and either installed or initialized before use e.g. on os x:

```py
import libmsym as msym
msym.init(library_location='/<libmsym installation path>/libmsym.dylib')
```

methods dealing with SALCs etc. require numpy to be installed

## notes

v0.1.0 is not compatible with v0.2.0
