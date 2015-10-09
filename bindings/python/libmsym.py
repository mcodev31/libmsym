from ctypes import *
from ctypes.util import find_library

class MSymError(Exception):
    def __init__(self, value, details):
         super().__init__(value)
         self.details = details
    def __str__(self):
        return repr(super().__str__() + ": " + self.details)

libmsym = CDLL(find_library('libmsym'))
 
class Element(Structure):
    _fields_ = [("_id", c_void_p),
                ("mass", c_double),
                ("_v", c_double*3),
                ("charge", c_int),
                ("_name",c_char*4)]
    
    @property
    def coordinates(self):
        return self._v
    @coordinates.setter
    def coordinates(self, coordinates):
        self._v = (c_double*3)(*coordinates)
    @property
    def name(self):
        return self._name.decode()
    @name.setter
    def name(self, name):
        self._name = name.encode('ascii')

class _RealSphericalHarmonic(Structure):
    _fields_ = [("n", c_int),
                ("l", c_int),
                ("m", c_int)]

class _BasisFunctionUnion(Union):
    _fields_ = [("_sh", _RealSphericalHarmonic)]

class BasisFunction(Structure):
    _fields_ = [("_id", c_void_p),
                ("_type", c_int),
                ("_element", POINTER(Element)),
                ("_f", _BasisFunctionUnion),
                ("_name",c_char*8)]

    def __init__(self, element=None):
        if element == None:
            raise MSymError("Basis funtion requires an element")
        super().__init__()
        self.element = element

    def _set_element_pointer(self, element):
        self._element = pointer(element)

    @property
    def name(self):
        return self._name.decode()
    @name.setter
    def name(self, name):
        self._name = name.encode('ascii')

class RealSphericalHarmonic(BasisFunction):
    def __init__(self, element=None, n=0, l=0, m=0, name=""):
        super().__init__(element=element)
        self._type = 0
        self._f._sh.n = n
        self._f._sh.l = l
        self._f._sh.m = m
        self.name = name

    @property
    def n(self):
        return self._f._sh.n
    @n.setter
    def n(self, n):
        self._f._sh.n = n

    @property
    def l(self):
        return self._f._sh.l
    @l.setter
    def l(self, n):
        self._f._sh.n = l

    @property
    def m(self):
        return self._f._sh.m
    @m.setter
    def m(self, n):
        self._f._sh.n = m
    

class Context(object):
    _Context = POINTER(type('msym_context', (Structure,), {}))

    libmsym.msymCreateContext.restype = _Context
    libmsym.msymCreateContext.argtypes = []

    libmsym.msymReleaseContext.restype = c_int
    libmsym.msymReleaseContext.argtypes = [_Context]
    
    libmsym.msymErrorString.argtypes = [c_int]
    libmsym.msymErrorString.restype = c_char_p

    libmsym.msymGetErrorDetails.argtypes = []
    libmsym.msymGetErrorDetails.restype = c_char_p
    
    libmsym.msymFindSymmetry.restype = c_int
    libmsym.msymFindSymmetry.argtypes = [_Context]

    libmsym.msymSetPointGroupByName.restype = c_int
    libmsym.msymSetPointGroupByName.argtypes = [_Context, c_char_p]

    libmsym.msymGetPointGroupName.restype = c_int
    libmsym.msymGetPointGroupName.argtypes = [_Context, c_int, c_char_p]

    libmsym.msymSetElements.restype = c_int
    libmsym.msymSetElements.argtypes = [_Context, c_int, POINTER(Element)]

    libmsym.msymGetElements.restype = c_int
    libmsym.msymGetElements.argtypes = [_Context, POINTER(c_int), POINTER(POINTER(Element))]

    libmsym.msymSymmetrizeElements.restype = c_int
    libmsym.msymSymmetrizeElements.argtypes = [_Context]

    libmsym.msymSetBasisFunctions.restype = c_int
    libmsym.msymSetBasisFunctions.argtypes = [_Context, c_int, POINTER(BasisFunction)]

    libmsym.msymGenerateSALCSubspaces.restype = c_int
    libmsym.msymGenerateSALCSubspaces.argtypes = [_Context]


    def __init__(self, elements=[], basis_functions=[], point_group="", _func=libmsym.msymCreateContext):
        self._elements = []
        self._basis_functions = []
        self._point_group = None
        self._ctx = _func()
        if not self._ctx:
            raise RuntimeError
        if len(elements) > 0:
            self._set_elements(elements)
        if len(basis_functions) > 0:
            self._set_basis_functions(basis_functions)
        if len(point_group) > 0:
            self._set_point_group(point_group)
            self.find_symmetry()
            
    def __del__(self):
        self._destruct()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._destruct()

    def _destruct(self, _func=libmsym.msymReleaseContext): 
        if self._ctx:
            _func(self._ctx)
        self._ctx = None

    @staticmethod
    def _assert_success(error, _error_string=libmsym.msymErrorString, _error_details=libmsym.msymGetErrorDetails):
        if not error == 0:
            raise MSymError(_error_string(error).decode(), _error_details().decode())
        
    def _update_elements(self, _func=libmsym.msymGetElements):
        celements = POINTER(Element)()
        csize = c_int(0)
        self._assert_success(_func(self._ctx,byref(csize),byref(celements)))
        self._elements = celements[0:csize.value]

    def _update_point_group(self, _func=libmsym.msymGetPointGroupName):
        cname = (c_char*8)()
        self._assert_success(_func(self._ctx,sizeof(cname),cname))
        self._point_group = cname.value.decode()

    def _set_elements(self, elements, _func=libmsym.msymSetElements):
        if not self._ctx:
            raise RuntimeError
        size = len(elements)
        element_array = (Element*size)(*elements)
        self._assert_success(_func(self._ctx, size, element_array))
        self._element_array = element_array
        self._elements = elements

    def _set_point_group(self, point_group, _func=libmsym.msymSetPointGroupByName):
        cname = c_char_p(point_group.encode('ascii'))
        self._assert_success(_func(self._ctx, cname))

    def _set_basis_functions(self, basis_functions,_func=libmsym.msymSetBasisFunctions):
        if not self._ctx:
            raise RuntimeError
        size = len(basis_functions)
        for bf in basis_functions:
            bf._set_element_pointer(self._element_array[basis_functions.index(bf)])
            
        self._assert_success(_func(self._ctx, size, (BasisFunction*size)(*basis_functions)))
        self._basis_functions = basis_functions
        
    @property
    def elements(self):
        return self._elements

    @elements.setter
    def elements(self, elements):
        self._set_elements(elements)

    @property
    def basis_functions(self):
        return self._basis_functions

    @basis_functions.setter
    def basis_functions(self, basis_functions):
        self._set_basis_functions(basis_functions)

    @property
    def point_group(self):
        return self._point_group

    @point_group.setter
    def point_group(self, point_group):
        self._set_point_group(point_group)

    def find_symmetry(self, _func=libmsym.msymFindSymmetry):
        if not self._ctx:
            raise RuntimeError
        self._assert_success(_func(self._ctx))
        self._update_point_group()
        return self._point_group

    def symmetrize_elements(self, _func=libmsym.msymSymmetrizeElements):
        if not self._ctx:
            raise RuntimeError
        self._assert_success(_func(self._ctx))
        self._update_elements()
        return self._elements

    def generate_salc_subspaces(self, _func=libmsym.msymGenerateSALCSubspaces):
        if not self._ctx:
            raise RuntimeError
        self._assert_success(_func(self._ctx))





