#!/bin/bash

#
#  Simple build script for seiscope optimization toolbox.
#
#  Requires: FoBiS
#

MODCODE='seiscope_optimization_toolbox.f90'    # module file name
LIBOUT='libOPTIM.a'           # name of library
DOCDIR='./doc/'                 # build directory for documentation
SRCDIR='./src/'                 # library source directory
TESTSRCDIR='./src/tests/'       # unit test source directory
BINDIR='./bin/'                 # build directory for unit tests
LIBDIR='./lib/'                 # build directory for library

#compiler flags:

FCOMPILER='gnu' #Set compiler to gfortran
FCOMPILERFLAGS='-c -O2 -std=f2008'
#FCOMPILER='intel' #Set compiler to intel
#FCOMPILERFLAGS='-c -O2 -warn -stand f08 -traceback'

#build using FoBiS:

if hash FoBiS.py 2>/dev/null; then

	echo "Building library..."

	FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${LIBDIR} -s ${SRCDIR} -dmod ./ -dobj ./ -t ${MODCODE} -o ${LIBOUT} -mklib static -colors

	echo "Building test programs..."

	FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${BINDIR} -s ${TESTSRCDIR} -dmod ./ -dobj ./ -colors -libs ${LIBDIR}${LIBOUT} --include ${LIBDIR}

else
	echo "FoBiS.py not found! Cannot build library. Install using: sudo pip install FoBiS.py"
fi
