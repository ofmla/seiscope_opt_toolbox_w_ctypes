# Collect source files for library
set(f90-sources
    seiscope_optimization_toolbox.f90
)

# Collect source files from subdirectories
add_subdirectory(sotb)

# Build library
add_library(sotb ${f90-sources})

# Store generated mod-files in a separate folder
set(moddir "${CMAKE_CURRENT_BINARY_DIR}/fmod")
set_target_properties(sotb PROPERTIES
    Fortran_MODULE_DIRECTORY "${moddir}"
)

# Add mod-file directory to include search path for targets depending on the library
target_include_directories(sotb PUBLIC
    $<BUILD_INTERFACE:${moddir}>
    $<INSTALL_INTERFACE:${INSTALL_MODULEDIR}>
)

# Install the library
if(SKBUILD)
  message(STATUS "The project is built using scikit-build")
  install(TARGETS sotb LIBRARY DESTINATION sotb_wrapper ARCHIVE DESTINATION sotb_wrapper RUNTIME DESTINATION sotb_wrapper)
else()
  install(TARGETS sotb
      EXPORT sotb-targets
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_INCLUDEDIR}"
)
endif()

# Install the folder with the modfiles (needed by projects using the library)
install(DIRECTORY "${moddir}/" DESTINATION "${INSTALL_MODULEDIR}")
