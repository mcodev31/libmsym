# libmsym v0.2.2
libmsym is a C library dealing with point group symmetry in molecules.

## molecules
Determine, symmetrize and generate molecules of any point group as well determine/select subgroups.
Get symmetrically equivalent elements and symmetrize translation of any any element in a molecule.

## character tables
Generate character tables for any point group (complex characters form reducible representations)

## wave functions
Generate SALCs of real spherical harmonics with any angular momentum for point groups with real characters (Ci, Cs, Cnv, Dn, Dnh, Dnd, Td, O, Oh, I and Ih), as well as symmetrize orbitals, determine partner functions etc.

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

The libmsym module requires python 3.

If you have installed the libmsym library in a location that can be found by your loader (e.g. ldconfig):
```shell
cd ../bindings/python
# install libmsym module in user site
python setup.py install --user
# run example
python ./examples/msympy_example.py <input xyz-file> <output xyz-file>
```

If you want to install libmsym in a custom directory, the easies way it to use cmake:
```shell
# install libmsym shared library in $HOME/lib and the python module in the user site
cmake -DMSYM_BUILD_PYTHON:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=ON -DCMAKE_INSTALL_PREFIX=$HOME/lib -DMSYM_PYTHON_INSTALL_OPTS=--user ../.
make install
# run example
python ../bindings/python/examples/msympy_example.py <input xyz-file> <output xyz-file>
```

methods dealing with SALCs etc. require numpy to be installed

## notes

v0.1.0 is not compatible with v0.2.0
