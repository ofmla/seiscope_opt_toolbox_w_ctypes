cmake_minimum_required(VERSION 3.16)

# Set basic project properties
project(sotb
    VERSION 0.1
    DESCRIPTION "Modern Fortran Edition of SEISCOPE Optimization Toolbox."
    LANGUAGES Fortran
)

# Set up install directories
include(GNUInstallDirs)
set(INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/sotb")
set(INSTALL_MODULEDIR "${INSTALL_INCLUDEDIR}/fmod")

# Include project specific CMake functions
include(cmake/sotb.cmake)

# Include user customizable config settings
include(config.cmake)

# Configure the build type (see cmake/sotb.cmake)
configure_build_type()

# Configure the compiler flags (see cmake/sotb.cmake)
configure_compiler_flags()

# Build library
add_subdirectory(src)

# Build applications
list(APPEND apps PSTD PNLCG PLBFGS LBFGS PTRN TRN)
foreach(_app IN LISTS apps)
  add_subdirectory(app/${_app})
endforeach()

if(NOT SKBUILD)
  #
  # Install package config find, so that other CMake project can find this project
  #
  include(CMakePackageConfigHelpers)

  # If project uses customized finders, they should be installed with it
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
  install(
      DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sotb
  )
  endif()

  # Install project, add namespace
  install(
      EXPORT sotb-targets
      FILE "sotb-targets.cmake"
      NAMESPACE "sotb::"
      DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/sotb"
  )

  # Create and install CMake package config file
  configure_package_config_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/export/sotb-config.cmake.in"
      "${CMAKE_CURRENT_BINARY_DIR}/cmake/sotb-config.cmake"
      INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/sotb"
  )
  write_basic_package_version_file(
      "${CMAKE_CURRENT_BINARY_DIR}/cmake/sotb-config-version.cmake"
      VERSION "${PROJECT_VERSION}"
      COMPATIBILITY SameMajorVersion
  )
  install(
      FILES "${CMAKE_CURRENT_BINARY_DIR}/cmake/sotb-config.cmake"
      "${CMAKE_CURRENT_BINARY_DIR}/cmake/sotb-config-version.cmake"
      DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/sotb"
  )
endif()
