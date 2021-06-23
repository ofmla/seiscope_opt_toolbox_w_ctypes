!******************************************************************
!Copyright 2013-2016 SEISCOPE II project, All rights reserved.
!
!Redistribution and use in source and binary forms, with or
!without modification, are permitted provided that the following
!conditions are met:
!
!    *  Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!    *  Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution.
!    *  Neither the name of the SEISCOPE project nor the names of
!       its contributors may be used to endorse or promote products
!       derived from this software without specific prior written permission.
!
!Warranty Disclaimer:
!THIS SOFTWARE IS PROVIDED BY THE SEISCOPE PROJECT AND CONTRIBUTORS
!"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!SEISCOPE PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
!STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
!IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!POSSIBILITY OF SUCH DAMAGE.

module opt_PSTD
use typedef
use miscellaneous

private
real, allocatable, dimension(:) :: xk, descent
public PSTD, xk, descent

contains

subroutine pstd_centry(n,x,fcost,grad,grad_preco,optim,flag,lb,ub) bind(c, name='PSTD')
  use, intrinsic :: iso_c_binding
  !IN
  integer(c_int)  :: n                           !dimension of the problem
  real(c_float)   :: fcost                       !cost associated with x
  real(c_float),dimension(n) :: grad,grad_preco !gradient and preconditioned gradient at x 
  !IN/OUT  
  integer(c_int)  :: flag
  real(c_float),dimension(n) :: x               !current point
  type(optim_type) :: optim            !data structure 
  type(c_ptr), value  :: lb,ub
  real(c_float), pointer, dimension(:) :: lb_pass, ub_pass

  nullify(lb_pass)
  nullify(ub_pass)
  if (c_associated(ub)) call c_f_pointer(ub, ub_pass, (/n/))
  if (c_associated(lb)) call c_f_pointer(lb, lb_pass, (/n/))
  call pstd(n,x,fcost,grad,grad_preco,optim,flag,lb_pass,ub_pass)
end subroutine pstd_centry

!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! This routine is used to deallocate all the  !
! arrays that have been used during the       !
! PSTD optimization                           !
!---------------------------------------------!
! INPUT/OUT:  optim_type optim                !
!---------------------------------------------!
subroutine finalize_PSTD()
  
  implicit none 
  
  deallocate(xk)
  deallocate(descent)
  
end subroutine finalize_PSTD
!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routine is used for the initialization of the  !
! reverse communication mode preconditioned steepest  !
! descent algorithm from the TOOLBOX.                 !
!                                                     !
! The parameters for the linesearch are set here, as  !
! well as the data structure allocations required for !
! the use of the algorithm                            !
!-----------------------------------------------------!
! INPUT : integer n, dimension of the problem         !
!       : real fcost                                  !
!       : real,dimension(n) x,grad,grad_preco         !
!                           initial guess, gradient,  !
!                           preconditioned gradient   !
! INPUT/OUTPUT : optim_type optim data structure      ! 
!-----------------------------------------------------!
subroutine init_PSTD(n,x,fcost,grad_preco,optim)
  
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
  optim%nls_max=20 ! max number of linesearch
  optim%cpt_ls=0
  optim%first_ls=.true.
  optim%alpha=1. ! first value for the linesearch steplength
  
  !---------------------------------------!
  ! memory allocations                    !
  !---------------------------------------!
  allocate(xk(n))
  xk(:)=x(:)
  allocate(descent(n))

  !---------------------------------------!
  ! first descent direction               !
  !---------------------------------------!
  descent(:)=-1.*grad_preco(:)

end subroutine init_PSTD
!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routine is the reverse communication mode      !
! preconditioned steepest descent algorithm. This     ! 
! routine performs an iterative minimization of a     !
! function f following the recurrence                 !
!                                                     !
! x_{0}=x0                                            !
! x_{k+1}=x_k+\alpha_k d_k                            !
!                                                     !
! where the descent direction d_k is                  !
!                                                     !
! d_k=-Q_k \nabla f_k                                 !
!                                                     !
! with Q_k        : preconditioner at iteration k     !
!      \nabla f_k : gradient of f in x_k              !
!                                                     !
! and alpha_k is the steplength computed through the  !
! common linesearch algorithm of the TOOLBOX          !
!                                                     !
! The first call to the algorithm must be done with   !
! FLAG='INIT'. For this first call, the initial point !
! x0 is given through the variable x, and the input   !
! variable fcost and grad must correspond respectively!
! to the misfit and gradient at x0.                   !
!                                                     !
! The reverse communication with the user is          !
! performed through the variable FLAG. This           !
! variable indicates to the user on return what action! 
! he has to do, or the state of the algorithm.        !
! Possible values are                                 !
! - FLAG='GRAD' => the user must compute the cost and !
!                  (preconditioned) gradient at       !
!                  current point x                    !
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
subroutine PSTD(n,x,fcost,grad,grad_preco,optim,flag,lb,ub)
  use, intrinsic :: iso_c_binding
  implicit none
  
  !IN
  integer  :: n                        !dimension of the problem
  real   :: fcost                      !cost associated with x
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
     call init_PSTD(n,x,fcost,grad_preco,optim)
     call std_linesearch(n,x,fcost,grad,xk,descent,optim,lb,ub)
     call print_info(n,'ST',optim,fcost,grad,flag)
     FLAG=1     
     optim%nfwd_pb=optim%nfwd_pb+1
  else
     !-----------------------------------------------------!
     ! else call the linesearch process                    !  
     !-----------------------------------------------------!
     call std_linesearch(n,x,fcost,grad,xk,descent,optim,lb,ub)
     if(optim%task.eq.0) then !NEW STEP
        !-----------------------------------------------------!
        ! test for convergence                                !
        !-----------------------------------------------------!        
        optim%cpt_iter=optim%cpt_iter+1        
        call std_test_conv(optim,fcost,test_conv)
        if(test_conv) then
           FLAG=2
           call print_info(n,'ST',optim,fcost,grad,flag)
           call finalize_PSTD
        else
           FLAG=3
           !-----------------------------------------------------!
           ! if a NEW_STEP is taken, compute a new descent       !
           ! direction using current descent, gradient and       !
           ! preconditioned gradient                             !          
           !-----------------------------------------------------!
           descent(:)=-1.*grad_preco(:)                
           !print info on current iteration           
           call print_info(n,'ST',optim,fcost,grad,flag)
        endif
     elseif(optim%task.eq. 1) then
        !-----------------------------------------------------!
        ! if the linesearch needs a new gradient then ask the !  
        ! user to provide it
        !-----------------------------------------------------!
        FLAG=1         
        optim%nfwd_pb=optim%nfwd_pb+1
     elseif(optim%task.eq. 2) then
        !-----------------------------------------------------!
        ! if the linesearch has failed, inform the user       !
        !-----------------------------------------------------!
        FLAG=4
        !print info on current iteration
        call print_info(n,'ST',optim,fcost,grad,flag)
        call finalize_PSTD
     endif
  endif
  
end subroutine PSTD

end module opt_PSTD

