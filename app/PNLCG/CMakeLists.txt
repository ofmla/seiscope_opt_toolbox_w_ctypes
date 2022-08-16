set(f90-sources
    test_PNLCG.f90
)

add_executable(test_PNLCG ${f90-sources})
target_link_libraries(test_PNLCG sotb)

if(SKBUILD)
  install(TARGETS test_PNLCG DESTINATION sotb_wrapper)
else()
  install(TARGETS test_PNLCG DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
