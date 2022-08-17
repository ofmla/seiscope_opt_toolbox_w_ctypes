set(f90-sources
    test_PTRN.f90
)

add_executable(test_PTRN ${f90-sources})
target_link_libraries(test_PTRN sotb)

if(SKBUILD)
  install(TARGETS test_PTRN DESTINATION sotb_wrapper)
else()
  install(TARGETS test_PTRN DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
