from ctypes import *
from ctypes.util import find_library

class MSymError(Exception):
     pass

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

class Context(object):
    _Context = POINTER(type('msym_context', (Structure,), {}))

    libmsym.msymCreateContext.restype = _Context
    libmsym.msymCreateContext.argtypes = []

    libmsym.msymReleaseContext.restype = c_int
    libmsym.msymReleaseContext.argtypes = [_Context]
    
    libmsym.msymErrorString.argtypes = [c_int]
    libmsym.msymErrorString.restype = c_char_p
    
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
    

    def __init__(self, func=libmsym.msymCreateContext):
        self.point_group = None
        self._ctx = func()
        
    def __del__(self):
        self.destruct()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.destruct()

    def destruct(self, func=libmsym.msymReleaseContext): 
        if self._ctx:
            func(self._ctx)
        self._ctx = None

    @staticmethod
    def _assertSuccess(success, func=libmsym.msymErrorString):
        if not success == 0:
            raise Error(func(success).decode())
        
    def _updateElements(self, func=libmsym.msymGetElements):
        celements = POINTER(Element)()
        csize = c_int(0)
        self._assertSuccess(func(self._ctx,byref(csize),byref(celements)))
        self._elements = celements[0:csize.value]

    def _updatePointGroup(self, func=libmsym.msymGetPointGroupName):
        buf = (c_char*6)()
        self._assertSuccess(func(self._ctx,6,buf))
        self.point_group = buf.value.decode()
        
    @property
    def elements(self):
        return self._elements
    @elements.setter
    def elements(self, elements, func=libmsym.msymSetElements):
        if not self._ctx:
            raise RuntimeError
        size = len(elements)
        self._assertSuccess(func(self._ctx, size, (Element*size)(*elements)))
        self._elements = elements

    def findSymmetry(self, func=libmsym.msymFindSymmetry):
        if not self._ctx:
            raise RuntimeError
        self._assertSuccess(func(self._ctx))
        self._updatePointGroup()
        return self.point_group

    def symmetrizeElements(self, func=libmsym.msymSymmetrizeElements):
        if not self._ctx:
            raise RuntimeError
        self._assertSuccess(func(self._ctx))
        self._updateElements()
        return self._elements





