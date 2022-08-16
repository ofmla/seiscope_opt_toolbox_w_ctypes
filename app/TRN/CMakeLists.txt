set(f90-sources
    test_TRN.f90
)

add_executable(test_TRN ${f90-sources})
target_link_libraries(test_TRN sotb)

if(SKBUILD)
  install(TARGETS test_TRN DESTINATION sotb_wrapper)
else()
  install(TARGETS test_TRN DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
