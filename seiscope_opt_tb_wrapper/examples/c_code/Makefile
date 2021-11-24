#----------------------------------------------------------------------
#
# This makefile will build the C example program whose source is contained
# in this directory.  

# define some Makefile variables for the compiler and compiler flags
# to use Makefile variables later in the Makefile: $()
CC  = gcc
FC  = gfortran
RM	= rm -f

# Compile with gcc and link with the Fortran-compiler.
#-----------------------------------------------------------------
# Use the following definition to link the sample program with
# the shared OPTIM library dynamically.  
LIBS = -L../../seiscope_opt_toolbox_w_ctypes/lib -lOPTIM
#-----------------------------------------------------------------
# The -Wl, option is used to avoid the need to set LD_LIBRARY_PATH before execution
LFLAGS = -Wl,-rpath ../../seiscope_opt_toolbox_w_ctypes/lib $(LIBS)

# typing 'make' will invoke the first target entry in the file 
# (in this case the default target entry)

default: main

# To create the executable file main we need the objec file
# main.o:
#
main: main.o
	$(FC) main.o -o main $(LFLAGS)


# To create the object file main.o, we need the source
# files main.c and optim_bc.h:
#
main.o: main.c optim_bc.h
	$(CC) -c main.c

# To start over from scratch, type 'make clean'.  This
# removes the executable file, as well as old .o object
# files and *~ backup files:
#	
clean:
	$(RM) main *.o *~
