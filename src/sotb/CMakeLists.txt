set(curdir "${CMAKE_CURRENT_SOURCE_DIR}")

# Add further source from this directory to the library
list(APPEND f90-sources
    "${curdir}/mdle_typedef.f90"
    "${curdir}/merged-common.f90"
    "${curdir}/merged-LBFGS.f90"
    "${curdir}/merged-PLBFGS.f90"
    "${curdir}/merged-PNLCG.f90"
    "${curdir}/merged-PSTD.f90"
    "${curdir}/merged-PTRN.f90"
    "${curdir}/merged-TRN.f90"
)

set(f90-sources "${f90-sources}" PARENT_SCOPE)
