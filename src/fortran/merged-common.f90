module miscellaneous
use typedef

public print_info, project, std_linesearch, std_test_conv

contains

subroutine rosenbrock(x,fcost,grad) bind(c, name="rosenbrock")
  
  implicit none
  !IN
  real(c_float),dimension(2) :: x
  !IN/OUT
  real(c_float) :: fcost
  real(c_float),dimension(2) :: grad
 
  fcost=(1.-x(1))**2+100.*(x(2)-x(1)**2)**2
  grad(1)=2.*(x(1)-1.)-400.*x(1)*(x(2)-x(1)**2)
  grad(2)=200.*(x(2)-x(1)**2)

end subroutine rosenbrock

!*******************************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX                    *!
!*******************************************************!
! This routine is used to generate formatted            !
! ascii output files containing information             !
! on the convergence history.                           !
! The name of these files depend on the                 !
! optimization routine used:                            !
! - steepest descent                : iterate_ST.dat    !
! - nonlinear conjugate gradient    : iterate_CG.dat    !
! - l-BFGS                          : iterate_LB.dat    !
! - preconditioned l-BFGS           : iterate_PLB.dat   !
!                                                       !
! The truncated Newton method and the preconditioned    ! 
! truncated Newont method use their own printing        !
! subroutines to generate their output files since      !
! their format slightly differ.                         !
!                                                       !
! This printing routine is called basically at the      !
! initialization of the optimization and then at each   !
! iteration of the solver                               !
!-------------------------------------------------------!
! INPUT  : character*2 routine_type (solver)            !
!          character*4 FLAG (phase of the solver:       !
!                            initialization,            !
!                            iteration,etc...           !
!          integer n (dimension of the problem)         !
!          optim_type optim (data structure)            !
!-------------------------------------------------------!
subroutine print_info(n,routine_type,optim,fcost,grad,FLAG)

  implicit none

  !IN
  character*2 :: routine_type
  !character*4 :: FLAG
  integer :: FLAG
  integer :: n
  real :: fcost
  real,dimension(n) :: grad
  type(optim_type) :: optim !data structure   
  !Local variables
  real :: ng

  if(optim%print_flag.eq.1) then     
     ng= norm2(grad)
     if(FLAG.eq.0) then
        if(routine_type.eq.'ST') then
           open(10,file='iterate_ST.dat')     
           write(10,'(A70)') '**********************************************************************'
           write(10,'(A50)') '         STEEEPEST DESCENT ALGORITHM         '
           write(10,'(A70)') '**********************************************************************'
        elseif(routine_type.eq.'CG') then
           open(10,file='iterate_CG.dat')     
           write(10,'(A70)') '**********************************************************************'
           write(10,'(A50)') '     NONLINEAR CONJUGATE GRADIENT ALGORITHM  '
           write(10,'(A70)') '**********************************************************************'
        elseif(routine_type.eq.'LB') then
           open(10,file='iterate_LB.dat')     
           write(10,'(A70)') '**********************************************************************'
           write(10,'(A50)') '             l-BFGS ALGORITHM                '
           write(10,'(A70)') '**********************************************************************'
        elseif(routine_type.eq.'PL') then
           open(10,file='iterate_PLB.dat')     
           write(10,'(A70)') '**********************************************************************'
           write(10,'(A50)') '             PRECONDITIONED l-BFGS ALGORITHM                '
           write(10,'(A70)') '**********************************************************************'
        endif
        write(10,'(A30,ES10.2)') '     Convergence criterion  : ',optim%conv
        write(10,'(A30,I7)')     '     Niter_max              : ',optim%niter_max
        write(10,'(A30,ES10.2)') '     Initial cost is        : ',optim%f0
        write(10,'(A30,ES10.2)') '     Initial norm_grad is   : ',ng
        write(10,'(A70)') '**********************************************************************'
        write(10,'(A10,A10,A10,A10,A12,A12,A12)')&
             '   Niter   ' ,&
             '   fk      ' ,&
             '   ||gk||  ' ,&
             '     fk/f0   ' ,&
             '       alpha   ' ,&
             '       nls     ',&
             '   ngrad    '
        write(10,'(I6,ES12.2,ES12.2,ES12.2,ES12.2,I8,I8)') &
             optim%cpt_iter,&          
             fcost,&
             ng,&
             fcost/optim%f0,&
             optim%alpha,&
             optim%cpt_ls,&
             optim%nfwd_pb
     elseif(FLAG.eq.2) then
        write(10,'(I6,ES12.2,ES12.2,ES12.2,ES12.2,I8,I8)') &
             optim%cpt_iter,&
             fcost,&
             ng,&
             fcost/optim%f0,&
             optim%alpha,&
             optim%cpt_ls,&
             optim%nfwd_pb
        write(10,'(A70)') '**********************************************************************'
        if(optim%cpt_iter.ge.optim%niter_max) then
           write(10,'(A50)') '  STOP: MAXIMUM NUMBER OF ITERATION REACHED    '
        else
           write(10,'(A50)') '  STOP: CONVERGENCE CRITERION SATISFIED        '
        endif
        write(10,'(A70)') '**********************************************************************'
        close(10)
     elseif(FLAG.eq.4) then
        write(10,'(A70)') '**********************************************************************'
        write(10,'(A50)') '  STOP: LINESEARCH FAILURE    '
        write(10,'(A70)') '**********************************************************************'
        close(10)     
     else
        write(10,'(I6,ES12.2,ES12.2,ES12.2,ES12.2,I8,I8)') &
             optim%cpt_iter,&
             fcost,&
             ng,&
             fcost/optim%f0,&
             optim%alpha,&
             optim%cpt_ls,&
             optim%nfwd_pb    
     endif
  endif
  
end subroutine print_info

!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! The routine project compute the projection  !
! of the vector x into the box defined by     !
! the bound constraints optim%lb and optim%ub !
!---------------------------------------------!
! INPUT  : integer n                          !
!          optim_type optim                   !
! IN/OUT : real,dimension(n) x                !
!---------------------------------------------!
subroutine project(n,optim,x,lb,ub)

  implicit none
  
  !IN
  integer :: n
  type(optim_type) :: optim
  !IN/OUT
  real,dimension(n) :: x,ub,lb
  !Local variables
  integer :: i
  
  do i=1,n
     if(x(i).gt.ub(i)) then
        x(i)=ub(i)-optim%threshold
     endif
     if(x(i).lt.lb(i)) then
        x(i)=lb(i)+optim%threshold
     endif
  enddo
  
end subroutine project

!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! This routine is the reverse communication   ! 
! mode linesearch shared by all the           ! 
! optimization routines of the toolbox.       !
!                                             !                            
! This linesearch enforces the Wolfe          !
! conditions: sufficient decrease             !
!             sufficient curvature            !
! The Wolfe conditions can be found in        !
! Nocedal, Numerical Optimization, 2nd        !
! edition, p.33                               !
!                                             !
! The linesearch method implemented here is   !
! based first on a bracketing strategy, then  !
! on a dichotomy algorithm. A full description!
! of this strategy can be found in            !
! Numerical Optimizationn Theoretical and     !
! Practical Aspects, J.F.Bonnans, J.C.Gilbert,!
! C. Lemaréchal, C.A. Sagastizábal,           !
! Springer-Verlag, Universitext               !
!                                             !
!---------------------------------------------!
! INPUT  : integer :: n (dimension)           !
!          real fcost (current cost)          !
!          real,dimension(n) grad             !
! OUTPUT : real,dimension(n) x                !
!          optim_typ optim (data structure)   !
!---------------------------------------------!

subroutine std_linesearch(n,x,fcost,grad,xk,descent,optim,lb,ub) !lb,ub,optim
  implicit none

  !IN
  integer :: n
  real :: fcost
  real,dimension(n) :: grad,descent,xk
  !IN/OUT
  real,dimension(n) :: x
  real,optional, dimension(:) :: lb,ub
  type(optim_type) :: optim
  !Local variables
  real :: q0,new_alpha
    
  if(optim%first_ls) then
     !---------------------------------------!
     ! FIRST LINESEARCH: initialization step !
     !---------------------------------------!
     optim%fk=fcost
     q0 = dot_product(grad,descent)
     !print*,q0
     optim%q0=q0
     !set the search interval bounds to 0
     optim%alpha_L=0.
     optim%alpha_R=0.          
     optim%task=1
     optim%first_ls=.false.
     xk(:)=x(:)
     x(:)=xk(:)+optim%alpha*descent(:)
     !IF BOUNDS ACTIVATED, PROJECT x(:) TO THE FEASIBLE ENSEMBLE
     if(present(lb)) then
     	print*,'marica'
        call project(n,optim,x,lb,ub)
     endif
     optim%cpt_ls=0
  elseif( (optim%cpt_ls.ge.optim%nls_max) .and. (fcost<optim%fk)) then     
     !-----------------------------------------------------------------------!
     ! if the number of linesearch iteration outreaches the maximum allowed  !
     ! but a decrease of the misfit is produced then accept the steplength   !
     !-----------------------------------------------------------------------!
     optim%task=0        
     optim%first_ls=.true.
     !Compute new x in the descent direction     
     x(:)=xk(:)+optim%alpha*descent(:)     
     !IF BOUNDS ACTIVATED, PROJECT x(:) INTO TO THE FEASIBLE ENSEMBLE
     if(present(lb)) then
        call project(n,optim,x,lb,ub)
     endif
  elseif(optim%cpt_ls.ge.optim%nls_max) then     
     !-----------------------------------------------------------------------!
     ! if the number of linesearch iteration outreaches the maximum allowed  !
     ! without decreasing the misfit then the linesearch has failed          !
     !-----------------------------------------------------------------------!
     optim%task=2
  else
     !-----------------------------------------------------------------------!
     ! If not initialization step and number of linesearch iteration ok      !
     ! then perform one linesearch iteration                                 !
     !-----------------------------------------------------------------------!
     optim%q = dot_product(grad,descent)
     if( (fcost.le.(optim%fk+(optim%m1*optim%alpha*optim%q0)))&
          .and.(optim%q.ge.(optim%m2*optim%q0)) ) then
        !--------------------------------------------------------------------!
        ! First test if the Wolfe conditions are satisfied with              !     
        ! current steplength, if this is the case, linesearch                ! 
        ! ends here with success                                             !
        !--------------------------------------------------------------------!
        optim%task=0
        optim%first_ls=.true.
        if(optim%debug) then
           if(optim%print_flag.eq.1) then
              write(10,*) 'fcost :',fcost
              write(10,*) 'optim%f0 :',optim%f0
              write(10,*) 'optim%fk :',optim%fk
              write(10,*) 'optim%alpha :',optim%alpha
              write(10,*) 'optim%q :', optim%q
              write(10,*) 'optim%q0 :',optim%q0
              write(10,*) 'm1 :',optim%m1
              write(10,*) 'cpt_ls is : ',optim%cpt_ls
           endif
        endif
     elseif (fcost>(optim%fk+(optim%m1*optim%alpha*optim%q0))) then
        !--------------------------------------------------------------------!
        ! If the first condition is not satisfied then shrink the            !
        ! search interval                                                    !
        !--------------------------------------------------------------------!
        if(optim%debug) then
           if(optim%print_flag.eq.1) then
              write(10,*) 'failure 1'
              write(10,*) 'fcost :',fcost
              write(10,*) 'optim%fk :',optim%fk
              write(10,*) 'optim%alpha :',optim%alpha
              write(10,*) 'optim%q0 :',optim%q0
              write(10,*) 'm1 :',optim%m1
              write(10,*) 'cpt_ls is : ',optim%cpt_ls
           endif
        endif
        optim%alpha_R=optim%alpha
        new_alpha=(optim%alpha_L+optim%alpha_R)/2.
        optim%alpha=new_alpha        
        optim%task=1
        optim%cpt_ls=optim%cpt_ls+1
     elseif( (fcost.le. (optim%fk+(optim%m1*optim%alpha*optim%q0)))&
          .and.(optim%q<(optim%m2*optim%q0) ) ) then
        !--------------------------------------------------------------------!
        ! If the second condition is not satisfied then shrink the           !
        ! search interval unless the right bound of the search interval      !
        ! as not yet been defined                                            !
        !--------------------------------------------------------------------!
        if(optim%debug) then
           if(optim%print_flag.eq.1) then
              write(10,*) 'failure 2'
              write(10,*) 'fcost :',fcost
              write(10,*) 'optim%fk :',optim%fk
              write(10,*) 'optim%alpha :',optim%alpha
              write(10,*) 'optim%q0 :',optim%q0
              write(10,*) 'optim%q :',optim%q
              write(10,*) 'm1 :',optim%m1
              write(10,*) 'm2 :',optim%m2
              write(10,*) 'cpt_ls is : ',optim%cpt_ls
           endif
        endif
        optim%alpha_L=optim%alpha
        if(optim%alpha_R.ne.0.) then
           new_alpha=(optim%alpha_L+optim%alpha_R)/2.
        else
           new_alpha=optim%mult_factor*optim%alpha
        endif
        optim%alpha=new_alpha        
        optim%task=1
        optim%cpt_ls=optim%cpt_ls+1                        
     endif
     !Compute new x in the descent direction
     x(:)=xk(:)+optim%alpha*descent(:)     
     ! IF BOUNDS ACTIVATED, PROJECT x(:) TO THE FEASIBLE ENSEMBLE
     if(present(lb)) then
        call project(n,optim,x,lb,ub)
     endif
  endif
end subroutine std_linesearch



!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! This routine implements a simple            !
! convergence test based on the relative      !
! cost  decrease.                             !
! If the current relative cost is lower than  !
! a value set by the user in optim%conv       !
! then test_conv is returned equal to .true.  !
! otherwise, it is returned equal to .false.  ! 
!---------------------------------------------!
! INPUT:  optim_type optim                    !
! OUTPUT: logical test_conv                   !
!---------------------------------------------!
subroutine std_test_conv(optim,fcost,test_conv)

  implicit none

  !IN
  real :: fcost
  type(optim_type) :: optim !data structure   
  !INT/OUT
  logical :: test_conv
    
  test_conv=((fcost/optim%f0<optim%conv).or.&
       (optim%cpt_iter.ge.optim%niter_max))
  
end subroutine std_test_conv

end module miscellaneous
