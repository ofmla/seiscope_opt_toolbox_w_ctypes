module typedef
use, intrinsic :: iso_c_binding, only: c_float, c_bool, c_char, c_int, c_ptr, c_loc, c_f_pointer
implicit none

type, bind(c) :: optim_type
                                                                      
!sequence

!#### DEBUG OPTION FOR USER
!# when debug is false: no information printed
!# when debug is true: additional information printed
logical(c_bool) :: debug

!#### BOUND CONSTRAINTS OPTION FOR USER
!# when bound is 0: no bound constraints
!# when bound is 1: bound constraints activated
!integer :: bound
real(c_float) :: threshold !#tolerance on bound constraints satisfaction

!#### PRINTING FLAG FOR MPI APPLICATION
integer(c_int) :: print_flag

!#### LINESEARCH PARAMETERS
logical(c_bool) :: first_ls
integer(c_int)   :: task !(0 -NEW_STEP, 1 -NEW_GRAD, 2 -FAILURE!)
integer(c_int)   :: nls_max,cpt_ls,nfwd_pb,cpt_iter,niter_max
real(c_float)    :: f0,fk,conv
real(c_float)    :: m1,m2,mult_factor,alpha_L,alpha_R,alpha
real(c_float)    :: q0, q

!### LBFGS PARAMETERS 
integer(c_int) :: cpt_lbfgs,l

end type optim_type

end module typedef
