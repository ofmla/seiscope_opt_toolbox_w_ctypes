set(f90-sources
    test_PLBFGS.f90
)

add_executable(test_PLBFGS ${f90-sources})
target_link_libraries(test_PLBFGS sotb)

if(SKBUILD)
  install(TARGETS test_PLBFGS DESTINATION sotb_wrapper)
else()
  install(TARGETS test_PLBFGS DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
