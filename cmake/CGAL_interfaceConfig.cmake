# MIT License
# Copyright (c) 2017 Benjamin Bercovici

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#



if(${IS_FORTUNA})
	set(CGAL_interface_INCLUDE_DIR /home/bebe0705/libs/local/include/CGAL_interface)
else()
	set(CGAL_interface_INCLUDE_DIR /usr/local/include/CGAL_interface/)
endif()

if (APPLE)
	set(CGAL_interface_LIBRARY /usr/local/lib/libCGAL_interface.dylib)
elseif(UNIX AND NOT APPLE)
	if(${IS_FORTUNA})
		set(CGAL_interface_LIBRARY /home/bebe0705/libs/local/lib/libCGAL_interface.so)
	else()
		set(CGAL_interface_LIBRARY /usr/local/lib/libCGAL_interface.so)
	endif()

else()
	message(FATAL_ERROR "Unsupported platform")
endif()

message("-- Found CGAL_interface: " ${CGAL_interface_LIBRARY})
