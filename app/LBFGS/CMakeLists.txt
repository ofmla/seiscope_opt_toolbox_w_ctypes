set(f90-sources
    test_LBFGS.f90
)

add_executable(test_LBFGS ${f90-sources})
target_link_libraries(test_LBFGS sotb)

if(SKBUILD)
  install(TARGETS test_LBFGS DESTINATION sotb_wrapper)
else()
  install(TARGETS test_LBFGS DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
