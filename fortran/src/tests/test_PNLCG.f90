!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This program provides an example for using the      !
! preconditioned nonlinear conjugate gradient         !
! algorithm within the TOOLBOX .                      !
!                                                     !
! The version implemented is the one by               !
! Y. DAI AND Y. YUAN, A nonlinear conjugate gradient  !
! method with a strong global convergence property,   !
! SIAM Journal on Optimization, 10 (1999), pp. 177–182!
!                                                     !
! See also Nocedal, Numerical optimization,           !
! 2nd edition p.132                                   !
!                                                     !
! The 2D Rosenbrock function is minimized             !
! f(x,y)  = (1-x)**2+100.*(y-x**2)**2                 !
!                                                     !
! The computation of this function and its gradient   !
! is implemented in                                   !
! ../../../COMMON/test/rosenbrock.f90                 !
!                                                     !
! The communication with the PNLCG optimizer is       ! 
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
!   - initialize grad_preco with the preconditioned   !
!     version of the gradient.                        !
!    %%%%%%%%%%%                                      !
!    %IMPORTANT%                                      !
!    %%%%%%%%%%%                                      !  
!    If the user wants to use the PNLCG algorithm     !
!    WITHOUT preconditioning he has to COPY grad      !
!    in grad_preco                                    !
! 2. On first call, FLAG must be set to INIT          !
! 3. When the FLAG returned by PNLCG is 'GRAD', then  !
!    the user must compute the cost, the gradient and !
!    the preconditioned gradient at current iterate x !
!    %%%%%%%%%%%                                      !
!    %IMPORTANT%                                      !
!    %%%%%%%%%%%                                      !
!    If the user wants to use the PNLCG algorithm     !
!    WITHOUT preconditioning he has to COPY grad in   !
!    grad_preco                                       !
!-----------------------------------------------------!
program test_PNLCG
  !use typedef
  !use rosenbrock_module, only: Rosenbrock
  !use opt_PNLCG
  use seiscope_optimization_toolbox
  implicit none  
  
  integer :: n                                ! dimension of the problem
  real :: fcost                               ! cost function value
  real,dimension(:),allocatable :: x          ! current point
  real,dimension(:),allocatable :: grad       ! current gradient
  real,dimension(:),allocatable :: grad_preco ! preconditioned current gradient
  type(optim_type) :: optim                   ! data structure for the optimizer
  !character*4 :: FLAG                         ! communication 
  integer :: FLAG                             ! communication 
  
  !----------------------------------------------------!
  ! parameter initialization                           !
  !----------------------------------------------------!
  n=2                     ! dimension
  !FLAG='INIT'             ! first flag
  FLAG=0            ! first flag
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
  ! copy of grad in grad_preco: no preconditioning in  !
  ! this test                                          !
  !----------------------------------------------------!
  grad_preco(:)=grad(:) ! 
  
  !----------------------------------------------------!
  ! optimization loop: while convergence not reached or!
  ! linesearch not failed, iterate                     !
  !----------------------------------------------------!
  !call PNLCG(n,x,fcost,grad,grad_preco,optim,FLAG)
  !do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))
  do while ((FLAG.ne.2).and.(FLAG.ne.4))
     call PNLCG(n,x,fcost,grad,grad_preco,optim,FLAG)
     !if(FLAG.eq.'GRAD') then 
     !print*,FLAG    
     if(FLAG.eq.1) then        
        !compute cost and gradient at point x
        call Rosenbrock(x,fcost,grad)
        ! no preconditioning in this test: simply copy grad in 
        ! grad_preco
        grad_preco(:)=grad(:) 
     endif
  enddo
  
  !Helpful console writings
  write(*,*) 'END OF TEST'
  write(*,*) 'FINAL iterate is : ', x(:)
  write(*,*) 'See the convergence history in iterate_CG.dat'
  
end program test_PNLCG
