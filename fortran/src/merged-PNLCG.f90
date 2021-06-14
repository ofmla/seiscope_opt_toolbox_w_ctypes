module opt_PNLCG
use, intrinsic :: iso_c_binding
use typedef
use miscellaneous
use opt_PSTD, only : xk, descent

real, allocatable, dimension(:), public :: grad_prev
real, allocatable, dimension(:), public :: descent_prev

public descent_PNLCG, finalize_PNLCG, init_PNLCG, PNLCG

!Module procedures have access to module variables through host association.

contains

subroutine pnlcg_centry(n,x,fcost,grad,grad_preco,optim,flag,lb,ub) bind(c, name='PNLCG')
  !IN
  integer(c_int)  :: n                          !dimension of the problem
  real(c_float)   :: fcost                      !cost associated with x
  real(c_float),dimension(n) :: grad,grad_preco !gradient and preconditioned gradient at x 
  !IN/OUT  
  integer(c_int)  :: flag
  real(c_float),dimension(n) :: x               !current point
  type(optim_type) :: optim                     !data structure 
  type(c_ptr), value  :: lb,ub
  real(c_float), pointer, dimension(:) :: lb_pass, ub_pass

  nullify(lb_pass)
  nullify(ub_pass)
  if (c_associated(ub)) call c_f_pointer(ub, ub_pass, (/n/))
  if (c_associated(lb)) call c_f_pointer(lb, lb_pass, (/n/))
  call pnlcg(n,x,fcost,grad,grad_preco,optim,flag,lb_pass,ub_pass)
end subroutine pnlcg_centry

!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! Computation of the descent direction given by the 
! preconditioned nonlinear conjugate gradient algorithm 
! of Dai and Yuan 
! Y. DAI AND Y. YUAN, A nonlinear conjugate gradient  !
! method with a strong global convergence property,   !
! SIAM Journal on Optimization, 10 (1999), pp. 177–182!
!                                                     !
! See also Nocedal, Numerical optimization,           !
! 2nd edition p.132                                   !
!-----------------------------------------------------!
! INPUT  : integer :: n (dimension)                   ! 
!          real,dimension(n) :: grad,grad_preco       !
! INPUT/OUTPUT : optim_typ optim (data structure)     !
!-----------------------------------------------------!
subroutine descent_PNLCG(n,grad,grad_preco)
  
  implicit none
  !IN
  integer :: n
  real,dimension(n) :: grad,grad_preco   
  !Local variables
  real :: gkpgk,skpk,beta
  real,dimension(:),allocatable :: sk
  
  !------------------------------------------------------------!
  ! Storing old descent direction                              !
  !------------------------------------------------------------!
  descent_prev(:)=descent(:)
  
  !------------------------------------------------------------!
  ! Computation of beta                                        !
  !------------------------------------------------------------!  ! 
  gkpgk =dot_product(grad,grad_preco)  
  allocate(sk(n))
  sk(:)=grad(:)-grad_prev(:)
  skpk = dot_product(sk,descent_prev)
  beta=gkpgk/skpk
  
  !------------------------------------------------------------!
  ! Safeguard (may be useful in some cases)                    !
  !------------------------------------------------------------! 
  if((beta.ge.1e5).or.(beta.le.-1e5)) then     
     beta=0.
  endif
  
  !------------------------------------------------------------!
  ! Computation of the descent direction                       !
  !------------------------------------------------------------! 
  descent(:)=-1.*grad_preco(:)+beta*descent_prev(:)
  
  !------------------------------------------------------------!
  ! Deallocation                                               ! 
  !------------------------------------------------------------!
  deallocate(sk)

end subroutine descent_PNLCG
!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! This routine is used to deallocate all the  !
! arrays that have been used during the       !
! PSTD optimization                           !
!---------------------------------------------!
! INPUT/OUT:  optim_type optim                !
!---------------------------------------------!
subroutine finalize_PNLCG
  
  implicit none
  
  deallocate(xk)
  deallocate(grad_prev)
  deallocate(descent)
  deallocate(descent_prev)
  
end subroutine finalize_PNLCG
!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routine is used for the initialization of the  !
! reverse communication mode preconditioned nonlinear !
! conjugate gradient algorithm from the TOOLBOX.      !
!                                                     !
! The parameters for the linesearch are set here, as  !
! well as the data structure allocations required for !
! the use of the algorithm                            !
!-----------------------------------------------------!
! INPUT : integer n, dimension of the problem         !
!       : real fcost                                  !
!       : real,dimension(n) x,grad,grad_preco         !
!                           first iterate, gradient   !
!                           preconditioned gradient   !
! INPUT/OUTPUT : optim_type optim data structure      ! 
!-----------------------------------------------------!
subroutine init_PNLCG(n,x,fcost,grad_preco,optim)
  
  implicit none

  !IN
  integer :: n
  real :: fcost
  real,dimension(n) :: x,grad_preco
  !IN/OUT
  type(optim_type) :: optim !data structure   

  !---------------------------------------!
  ! set counters                          !
  !---------------------------------------!
  optim%cpt_iter=0
  optim%f0=fcost
  optim%nfwd_pb=0

  !---------------------------------------!
  ! initialize linesearch parameters      !
  ! by default, the max number of         !
  ! linesearch iteration is set to 20     !
  ! and the initial steplength is set to 1!
  !---------------------------------------! 
  optim%m1=1e-4 ! Wolfe conditions parameter 1 (Nocedal value)
  optim%m2=0.9  ! Wolfe conditions parameter 2 (Nocedal value)
  optim%mult_factor=10 ! Bracketting parameter (Gilbert value)
  optim%fk=fcost
  optim%nls_max=30 ! max number of linesearch
  optim%cpt_ls=0
  optim%first_ls=.true.
  optim%alpha=1. ! first value for the linesearch steplength 
  
  !---------------------------------------!
  ! memory allocations                    !
  !---------------------------------------!
  allocate(grad_prev(n))
  allocate(descent_prev(n))
  allocate(xk(n))
  xk(:)=x(:)
  allocate(descent(n))
  
  !---------------------------------------!
  ! first descent direction               !
  !---------------------------------------!
  descent(:)=-1.*grad_preco(:)
  

end subroutine init_PNLCG
!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routine is the reverse communication mode      !
! preconditioned nonlinear conjugate gradient         !
! algorithm of Dai and Yuan                           !
! Y. DAI AND Y. YUAN, A nonlinear conjugate gradient  !
! method with a strong global convergence property,   !
! SIAM Journal on Optimization, 10 (1999), pp. 177–182!
!                                                     !
! See also Nocedal, Numerical optimization,           !
! 2nd edition p.132                                   !
!                                                     !
! This routine performs an iterative                  !
! minimization of a function f following the          !
! recurrence                                          !
!                                                     !
! x_{0}=x0                                            !
! x_{k+1}=x_k+\alpha_k d_k                            !
!                                                     !
! where the descent direction d_k is                  !
!                                                     !
! d_k=-Q_k \nabla f_k + \beta_k d_{k-1}               !
!                                                     !
! where Q_k       : preconditioner at iteration k     !
!      \nabla f_k : gradient of f in x_k              !
!      \beta_k    : scalar computed through the       !
!                   Dai and Yuan formula              !
!                                                     !
! and alpha_k is the steplength computed through the  !
! common linesearch algorithm of the TOOLBOX          !
!                                                     !
! The first call to the algorithm must be done with   !
! FLAG='INIT'. For this first call, the initial point !
! x0 is given through the variable x. The input       !
! variables fcost and grad, grad_preco must correspond!
! respectively to the misfit, gradient and            ! 
! preconditioned gradient at x0.                      !
!                                                     !
! The reverse communication with the user is          !
! performed through the variable FLAG. This           !
! variable indicates to the user on return what action! 
! he has to do, or the state of the algorithm.        !
! Possible values are                                 !
! - FLAG='GRAD' => the user must compute the cost,    !
!                  gradient and preconditioned        !
!                  gradient at current point x        ! 
! - FLAG='CONV' => a minimizer has been found         !
! - FLAG='NSTE' => a new step is performed            !
! - FLAG='FAIL' => the linesearch has failed          !
!-----------------------------------------------------!
! INPUT  : integer :: n (dimension)                   ! 
!          real fcost (current cost)                  !
!          real,dimension(n) grad                     !
!          real,dimension(n) grad_preco               !
! INPUT/OUTPUT : real,dimension(n) x                  !
!                optim_typ optim (data structure)     !
!                character*4 FLAG (communication)     !
!-----------------------------------------------------!
subroutine PNLCG(n,x,fcost,grad,grad_preco,optim,flag,lb,ub)

  implicit none
  
  !IN
  integer  :: n                           !dimension of the problem
  real   :: fcost                       !cost associated with x
  real,dimension(n) :: grad,grad_preco !gradient and preconditioned gradient at x 
  !IN/OUT  
  integer  :: flag
  real,dimension(n) :: x               !current point
  real(c_float),optional, dimension(n) :: lb,ub
  type(optim_type) :: optim            !data structure   
  !Local variable
  logical :: test_conv
  
  if(FLAG.eq. 0) then
     !-----------------------------------------------------!
     ! if FLAG is INIT, call the dedicated initialization  !
     ! subroutine to allocate data structure optim and     !
     ! initialize the linesearch process                   !
     !-----------------------------------------------------!
     call init_PNLCG(n,x,fcost,grad_preco,optim)
     call std_linesearch(n,x,fcost,grad,xk,descent,optim,lb,ub) !lb,ub,optim
     call print_info(n,'CG',optim,fcost,grad,FLAG)
     FLAG=1     
     optim%nfwd_pb=optim%nfwd_pb+1
     !Store current gradient before the user compute the new one
     grad_prev(:)=grad(:)      
  else
     !-----------------------------------------------------!
     ! else call the linesearch process                    !  
     !-----------------------------------------------------!
     call std_linesearch(n,x,fcost,grad,xk,descent,optim,lb,ub) !lb,ub,optim    
     if(optim%task.eq. 0) then 
        optim%cpt_iter=optim%cpt_iter+1
        !-----------------------------------------------------!
        ! test for convergence                                !
        !-----------------------------------------------------!
        call std_test_conv(optim,fcost,test_conv)
        if(test_conv) then
           FLAG=2
           !print info on current iteration        
           call print_info(n,'CG',optim,fcost,grad,FLAG)
           call finalize_PNLCG
        else
           FLAG=3
           !-----------------------------------------------------!
           ! if a NEW_STEP is taken, compute a new descent       !
           ! direction using current descent, gradient and       !
           ! preconditioned gradient                             !
           !-----------------------------------------------------!
           call descent_PNLCG(n,grad,grad_preco)                    
           !print info on current iteration        
           call print_info(n,'CG',optim,fcost,grad,FLAG)
        endif
     elseif(optim%task.eq. 1) then
        !-----------------------------------------------------!
        ! if the linesearch needs a new gradient then ask the !  
        ! user to provide it
        !-----------------------------------------------------!
        FLAG=1         
        optim%nfwd_pb=optim%nfwd_pb+1
        !Store current gradient before the user compute the new one
        grad_prev(:)=grad(:)        
     elseif(optim%task.eq. 2) then        
        !-----------------------------------------------------!
        ! if the linesearch has failed, inform the user       !
        !-----------------------------------------------------!
        FLAG=4
        !print info on current iteration
        call print_info(n,'CG',optim,fcost,grad,FLAG)
        call finalize_PNLCG
     endif
  endif
  
end subroutine PNLCG

end module opt_PNLCG

