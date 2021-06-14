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
logical(c_bool)  :: first_ls
integer(c_int)   :: task !(0 -NEW_STEP, 1 -NEW_GRAD, 2 -FAILURE!)
integer(c_int)   :: nls_max,cpt_ls,nfwd_pb,cpt_iter,niter_max
real(c_float)    :: f0,fk,conv
real(c_float)    :: m1,m2,mult_factor,alpha_L,alpha_R,alpha
real(c_float)    :: q0, q

!### LBFGS PARAMETERS 
integer(c_int) :: cpt_lbfgs,l

!### TRUNCATED NEWTON PARAMETERS
logical(c_bool)  :: conv_CG
integer(c_int)   :: cpt_iter_CG,niter_max_CG,nhess
integer(c_int)   :: CG_phase !(0 -INIT, 1 -IRUN)
integer(c_int)   :: comm !(1 -DES1, 2 -DES2, 3 -NSTE)
real(c_float)    :: qk_CG,qkm1_CG,hessian_term,eta,norm_grad,norm_grad_m1,norm_residual

!### PRECONDITIONED TRUNCATED NEWTON PARAMETERS
real(c_float)    :: dHd,res_scal_respreco,alpha_CG

end type optim_type

end module typedef
