# MIT License
#
# Copyright (c) 2020 Robert Grupp
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required(VERSION 3.13.0 FATAL_ERROR)

# Version format is YYYY.MM.DD for a release and YYYY.MM.DD.1 for devel after the release
project(xreg_ipcai_regi VERSION 2020.0.0.1)

###########################################
# find all of the dependencies:

# Note: Let find_package(xReg) determine whether OpenCL needs to be found
# or if the ICD loader is to be used
find_package(xreg REQUIRED)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake" ${CMAKE_MODULE_PATH})
if (XREG_USE_LEGACY_FIND_TBB)
  set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/tbb" ${CMAKE_MODULE_PATH})
endif ()

if (XREG_USE_SYSTEM_OPENCL)
  find_package(OpenCL REQUIRED)
endif ()

find_path(VIENNACL_INC_DIR viennacl/matrix.hpp PATHS ${ViennaCL_DIR} ${ViennaCL_DIR}/..)
if (${VIENNACL_INC_DIR} STREQUAL "VIENNACL_INC_DIR-NOTFOUND")
  message(SEND_ERROR "ViennaCL not found!")
endif ()

find_package(TBB REQUIRED)

find_package(Boost REQUIRED)

find_package(Eigen3 REQUIRED)

find_package(fmt REQUIRED)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(OpenCV REQUIRED)

find_package(NLopt REQUIRED)

###########################################

# build this source tree:

add_subdirectory(lib)

add_subdirectory(apps)

