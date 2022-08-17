set(f90-sources
    test_PSTD.f90
)

add_executable(test_PSTD ${f90-sources})
target_link_libraries(test_PSTD sotb)

if(SKBUILD)
  install(TARGETS test_PSTD DESTINATION sotb_wrapper)
else()
  install(TARGETS test_PSTD DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
