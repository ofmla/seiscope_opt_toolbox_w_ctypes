# Set up build type
function(configure_build_type)

    set(_build_types "Release;Debug")
    set(_default_build_type "Release")

    get_property(_multiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    if(_multiConfig)
        set(CMAKE_CONFIGURATION_TYPES "${_build_types}")
        message(STATUS "Build type: Multi-Config (build type selected at the build step)")
    else()
        if(NOT CMAKE_BUILD_TYPE)
            message(STATUS "Build type: ${_default_build_type} (default single-config)")
            set(CMAKE_BUILD_TYPE "${_default_build_type}" CACHE STRING "Build type" FORCE)
            set_property(CACHE CMAKE_BUILD_TYPE PROPERTY HELPSTRING "Choose the type of build")
            set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "${_build_types}")
        else()
            message(STATUS "Build type: ${CMAKE_BUILD_TYPE} (manually selected single-config)")
        endif()
    endif()

endfunction(configure_build_type)


# Copy user configured compiler flags into global compiler flags
macro(configure_compiler_flags)

    set(_languages "Fortran")

    if(CMAKE_BUILD_TYPE)
        set(_buildtypes ${CMAKE_BUILD_TYPE})
    else()
        set(_buildtypes ${CMAKE_CONFIGURATION_TYPES})
    endif()
    foreach(_buildtype IN LISTS _buildtypes)
        foreach (lang IN LISTS _languages)
            string(TOUPPER "${_buildtype}" _buildtype_upper)
            set(CMAKE_${lang}_FLAGS " ${${lang}_FLAGS}")
            set(CMAKE_${lang}_FLAGS_${_buildtype_upper} " ${${lang}_FLAGS_${_buildtype_upper}}")
            message(STATUS "Flags for ${lang}-compiler (build type: ${_buildtype}): "
                "${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${_buildtype_upper}}")
        endforeach()
    endforeach()
    unset(_buildtypes)
    unset(_buildtype)
    unset(_buildtype_upper)

endmacro(configure_compiler_flags)
