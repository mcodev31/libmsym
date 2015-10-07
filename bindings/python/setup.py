#
#  setup.py
#  libmsympy
#
#  Created by Marcus Johansson on 07/10/15.
#  Copyright (c) 2015 Marcus Johansson.
#
#  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
#

from distutils.core import setup, Extension

msympymodule = Extension('msympy',
		sources = ['msympymodule.c',
			'msympy_context.c',
			'msympy_element.c'
			],
#		include_dirs = ['/usr/local/include'],
#		library_dirs = ['/usr/local/lib'],
		libraries=["msym"]
)

setup (name = 'libmsympy',
       version = '0.2',
       description = 'libmsym python binding',
       ext_modules = [msympymodule])
