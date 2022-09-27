#
# User adjustable config options
#
# For a better overview, collect all variables here, which the user may want to override to
# influence the build behaviour
#


# Turn this on, if the libraries should be built as shared libraries
option(BUILD_SHARED_LIBS "Whether the libraries built should be shared" TRUE)

if (${BUILD_SHARED_LIBS})
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    message(STATUS "Darwin specific RPATH configuration")
    if(SKBUILD)
      set(CMAKE_MACOSX_RPATH FALSE)
      set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/sotb_wrapper")
    else()
      set(CMAKE_MACOSX_RPATH TRUE)
      set(CMAKE_SKIP_BUILD_RPATH FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
      set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
      set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
      set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    endif()
  else()
    file(RELATIVE_PATH relativeRpath
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
    )
    set(CMAKE_INSTALL_RPATH $ORIGIN $ORIGIN/${relativeRpath})
  endif()
endif()
#
# Fortran compiler dependent config options
#

if("GNU" STREQUAL "${CMAKE_Fortran_COMPILER_ID}")

    # Specific settings for the GNU compiler

    set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2018"
        CACHE STRING "General Fortran compiler flags")

    set(Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops"
        CACHE STRING "Extra Fortran compiler flags for Release build")

    set(Fortran_FLAGS_DEBUG "-g -Wall -pedantic -fbounds-check"
        CACHE STRING "Extra Fortran compiler flags for Debug build")

elseif("Intel" STREQUAL "${CMAKE_Fortran_COMPILER_ID}")

    # Specific settings for the Intel compiler

    set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand f18"
        CACHE STRING "General Fortran compiler flags")

    set(Fortran_FLAGS_RELEASE "-O3 -ip"
        CACHE STRING "Extra Fortran compiler flags for Release build")

    set(Fortran_FLAGS_DEBUG "-g -warn all -check -traceback"
        CACHE STRING "Extra Fortran compiler flags for Debug build")

elseif("NAG" STREQUAL "${CMAKE_Fortran_COMPILER_ID}")

    # Specific settings for the NAG compiler

    set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand f18"
        CACHE STRING "General Fortran compiler flags")

    set(Fortran_FLAGS_RELEASE "-O3"
        CACHE STRING "Extra Fortran compiler flags for Release build")

    set(Fortran_FLAGS_DEBUG "-g -f2008 -nan -C=all"
        CACHE STRING "Extra Fortran compiler flags for Debug build")

else()

    # Generic compiler settings (using CMake's default values)

    set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
        CACHE STRING "General Fortran compiler flags")

    set(Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
        CACHE STRING "Extra Fortran compiler flags for Release build")

    set(Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
        CACHE STRING "Extra compiler flags for Debug build")

endif()
