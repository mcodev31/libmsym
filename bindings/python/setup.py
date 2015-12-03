#
#  setup.py
#  libmsym
#
#  Created by Marcus Johansson on 07/10/15.
#  Copyright (c) 2015 Marcus Johansson.
#
#  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
#

from distutils.core import setup
import sys

if (not sys.version_info[0] is 3):
    sys.exit('libmsym module requires python 3')
    
setup(name='libmsym',
      version='0.2.4',
      description = 'libmsym python binding',
      license='MIT',
      author='Marcus Johansson',
      author_email='mcodev31@gmail.com',
      url='https://github.com/mcodev31/libmsym',
      packages=['libmsym']
      )
