!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This program provides an example for using the      !
! preconditioned steepest-descent algorithm within the! 
! TOOLBOX                                             !
!                                                     !
! The 2D Rosenbrock function is minimized             !
! f(x,y)  = (1-x)**2+100.*(y-x**2)**2                 !
!                                                     !
! The computation of this function and its gradient   !
! is implemented in                                   !
! ../../../COMMON/test/rosenbrock.f90                 !
!                                                     !
! The communication with the PSTD optimizer is        ! 
! achieved through the variable FLAG.                 !
!                                                     !
! Here are the actions requested for the user         !
! 1. Before first call, the user must set the         !
!    optimization parameters, namely                  !
!    - maximum number of iterations optim%niter_max   !
!    - tolerance for the stopping criterion           !
!      (this stopping criterion is the default one    !
!        in the TOOLBOX and is based on the decrease  !
!        of the misfit function scaled to one at the  !
!       first iteration                               !
!     - level of details in the output                !
!       * if optim%debug=.true. then the information  !
!         on the linesearch process is provided       !
!       * if optim%debug=.false. a more compact output!
!         file is provided                            !
!   - initialize x with an arbitrary initial guess    !
!   - initialize fcost with the cost function         !
!     associated with the initial guess               !
!   - initialize grad with the gradient associated    !
!     with the initial guess                          !
! 2. On first call, FLAG must be set to INIT          !
! 3. When the FLAG returned by PSTD is 'GRAD', then   !
!    the user must compute the cost and the gradient  !
!    associated with the current iterate x            !
!-----------------------------------------------------!
program test_PSTD
  !use typedef
  !use rosenbrock_module, only: Rosenbrock
  !use opt_PSTD
  use seiscope_optimization_toolbox
  implicit none  
  
  
  integer :: n                                ! dimension of the problem
  real :: fcost                               ! cost function value
  real,dimension(:),allocatable :: x          ! current point
  real,dimension(:),allocatable :: grad       ! current gradient
  real,dimension(:),allocatable :: grad_preco ! preconditioned gradient
  type(optim_type) :: optim                   ! data structure for the optimizer
  !character*4 :: FLAG                         ! communication FLAG
  integer :: FLAG                             ! communication FLAG
   
  !----------------------------------------------------!
  ! parameter initialization                           !
  !----------------------------------------------------!
  n=2                     ! dimension
  !FLAG='INIT'             ! first flag
  FLAG=0                  ! first flag
  optim%niter_max=10000   ! maximum iteration number 
  optim%conv=1e-8         ! tolerance for the stopping criterion
  optim%print_flag=1      ! print info in output files 
  optim%debug=.false.     ! level of details for output files


  !----------------------------------------------------!
  ! intial guess                                       !
  !----------------------------------------------------!
  allocate(x(n),grad(n),grad_preco(n)) 
  x(1)=1.5 
  x(2)=1.5
  
  !----------------------------------------------------!
  ! computation of the cost and gradient associated    !
  ! with the initial guess                             !
  !----------------------------------------------------!
  call Rosenbrock(x,fcost,grad)
  
  !----------------------------------------------------!
  ! optimization loop: while convergence not reached or!
  ! linesearch not failed, iterate                     !
  !----------------------------------------------------!
  !do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))  
  do while ((FLAG.ne.2).and.(FLAG.ne.4)) 
     call PSTD(n,x,fcost,grad,grad_preco,optim,FLAG)
     !if(FLAG.eq.'GRAD') then 
     if(FLAG.eq. 1) then       
        !compute cost and gradient at point x
        call Rosenbrock(x,fcost,grad)
        grad_preco(:)=grad(:)
     endif
  enddo
  
  !Helpful console writings
  write(*,*) 'END OF TEST'
  write(*,*) 'FINAL iterate is : ', x(:)
  write(*,*) 'See the convergence history in iterate_ST.dat'
  
end program test_PSTD


