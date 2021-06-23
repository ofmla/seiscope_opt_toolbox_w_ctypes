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

module opt_PTRN
use, intrinsic :: iso_c_binding
use typedef
use miscellaneous
use opt_PSTD, only : xk, descent
use opt_PNLCG, only : descent_prev
use opt_TRN, only : forcing_term_TRN, eisenvect

private
public PTRN

contains

subroutine ptrn_centry(n,x,fcost,grad,grad_preco,residual,residual_preco, &
					   d,Hd,optim,flag,lb,ub) bind(c, name='PTRN')
  !IN
  integer(c_int)  :: n                 !dimension of the problem
  real(c_float)   :: fcost             !cost associated with x
  real(c_float),dimension(n) :: grad, grad_preco   !gradient and preconditioned gradient at x 
  real(c_float),dimension(n) :: residual, residual_preco  
  real(c_float),dimension(n) :: d, Hd 
  !IN/OUT  
  integer(c_int)  :: flag
  real(c_float),dimension(n) :: x      !current point
  type(optim_type) :: optim            !data structure 
  type(c_ptr), value  :: lb,ub
  real(c_float), pointer, dimension(:) :: lb_pass, ub_pass

  nullify(lb_pass)
  nullify(ub_pass)
  if (c_associated(ub)) call c_f_pointer(ub, ub_pass, (/n/))
  if (c_associated(lb)) call c_f_pointer(lb, lb_pass, (/n/))
  call ptrn(n,x,fcost,grad,grad_preco,residual,residual_preco, &
			d,Hd,optim,flag,lb_pass,ub_pass)
end subroutine ptrn_centry

!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This file contains the routines for the computation !
! of the Newton descent direction using a             !
! preconditioned matrix-free conjugate gradient solver!
! The algorithm is taken from                         !
! Numerical Optimization, Nocedal, 2nd edition, 2006  !
! Algorithm 5.3 p.119                                 !
! The stopping criterion is taken from                !
! Choosing the forcing term in an Inexact Newton      !
! method, Eisenstat and Walker, 1994,                 !
! SIAM Journal on Scientific Computing 17 (1), 16-32  !
! The routine is implemented in a reverse             !
! communication format. The user is requested         !
! - to compute  Hessian-vector products at each       !
!   iteration of the process                          !
! - to apply its preconditioner to the residuals      !
!   at each iteration of the process                  !
! Because of the preconditioning, the routine has been!
! split into three routines, compared to the TRN      !
! version                                             !
! The first routine is descent_PTRN0 for              !
! initialization                                      !   
! The second routine is descent_PTRN1 for the first   !
! part of the conjugate gradient iteration, before    !
! preconditioning.                                    !
! The third routine is descent_PTRN2 for the second   !
! part of the conjugate gradient iteration, after     !
! preconditioning.                                    !

!-----------------------------------------------------!
! descent_PTRN0
!-----------------------------------------------------!
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         real,dimension(n) grad_preco preconditioned !
!                           current gradient          !
!         character*4 FLAG communication flag         !
! INPUT/OUTPUT : optim_type optim (data structure)    !
!-----------------------------------------------------!
subroutine descent_PTRN0(n,grad,grad_preco,residual, residual_preco,d,Hd,optim,FLAG)
  
  implicit none
  
  !IN
  integer :: FLAG                                          !communication FLAG
  integer :: n                                             !dimension of the problem
  real,dimension(n) :: grad                                !gradient at current point
  real,dimension(n) :: grad_preco                          !preconditioned gradient at current point
  real,dimension(n) :: d, Hd 
  real,dimension(n) :: residual, residual_preco  
  !IN/OUT
  type(optim_type) :: optim
  
  !-----------------------------------------------------!
  ! Initialization of the conjugate gradient process    !
  !-----------------------------------------------------!
  residual(:)=grad(:)
  residual_preco(:)=grad_preco(:)
  d(:)=-1.*residual_preco(:)
  Hd(:)=0.
  descent(:)=0.     
  optim%qk_CG=0.
  optim%hessian_term=0.     
  optim%norm_residual = norm2(residual)
  optim%conv_CG=.false.    
  optim%cpt_iter_CG=0
  call print_info_PTRN(optim,0e0,FLAG)     

end subroutine descent_PTRN0

!-----------------------------------------------------!
! descent_PTRN1
!-----------------------------------------------------!
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         character*4 FLAG communication flag         !
! INPUT/OUTPUT : optim_type optim (data structure)    !
!-----------------------------------------------------!
subroutine descent_PTRN1(n,grad,residual,residual_preco,d,Hd,optim)
  
  implicit none
  
  !IN
  integer :: n                               !dimension of the problem
  real,dimension(n) :: grad                  !gradient at current point
  real,dimension(n) :: d, Hd 
  real,dimension(n) :: residual, residual_preco  
  !IN/OUT
  type(optim_type) :: optim                  !data structure
  !Local variables
  real :: grad_term
  real,dimension(:),allocatable :: mgrad
  
  !-----------------------------------------------------!
  ! start one conjugate gradient iteration       !
  !-----------------------------------------------------!    
  optim%dHd = dot_product(d,Hd)        
  if(optim%dHd<0.) then
     !-----------------------------------------------------!
     ! if optim%dHd < 0, detection of a negative eigenvalue of   !
     ! the Hessian operator, stop the process              !        
     !-----------------------------------------------------!  
     optim%conv_CG=.true.
     write(21,*) 'Negative curvature'
     if(optim%cpt_iter_CG.eq.0) then
        !-----------------------------------------------------!
        ! if this is the first iteration, thenreturn the      !
        ! opposite of the preconditioned gradient as descent  !
        ! direction: preconditioned steepest descent direction!
        !-----------------------------------------------------!   
        descent(:)=d(:) 
        !-----------------------------------------------------!
        ! if the debug option is activated, compute the       !
        ! quadratic function minimized during the conjugate   !
        ! gradient process (check is this value decresae      !
        ! throughout the CG iterations )                      !
        !-----------------------------------------------------!
        if(optim%debug) then
           optim%res_scal_respreco = dot_product(residual,residual_preco)
           optim%alpha_CG=optim%res_scal_respreco/optim%dHd         
           optim%qkm1_CG=optim%qk_CG
           allocate(mgrad(n))
           mgrad(:)=-1.*grad(:)           
           grad_term = dot_product(descent,mgrad)           
           optim%hessian_term=optim%hessian_term+(optim%alpha_CG**2)*optim%dHd
           optim%qk_CG=-grad_term+0.5*optim%hessian_term  
           deallocate(mgrad)
        endif
     endif
  else         
     !-----------------------------------------------------!
     ! if optim%dHd > 0, then start one conjugate gradient       !
     ! iteration                                           !
     !-----------------------------------------------------!   
     !Update descent direction
     optim%res_scal_respreco = dot_product(residual,residual_preco)
     optim%alpha_CG=optim%res_scal_respreco/optim%dHd
     descent_prev(:)=descent(:)
     descent(:)=descent(:)+optim%alpha_CG*d(:)
     residual(:)=residual(:)+optim%alpha_CG*Hd(:)        
     !STOP HERE and wait for preconditioning
  end if
  
end subroutine descent_PTRN1

!-----------------------------------------------------!
! descent_PTRN2
!-----------------------------------------------------!
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         character*4 FLAG communication flag         !
! INPUT/OUTPUT : optim_type optim (data structure)    !
!-----------------------------------------------------!
subroutine descent_PTRN2(n,grad,residual,residual_preco,d,Hd,optim,FLAG)
  
  implicit none
  
  !IN
  integer :: FLAG                                     !communication FLAG
  integer :: n                                        !dimension of the problem
  real,dimension(n) :: grad                           !gradient at current point
  real,dimension(n) :: d, Hd 
  real,dimension(n) :: residual, residual_preco  
  !IN/OUT
  type(optim_type) :: optim                           !data structure
  !Local variables
  real :: res_scal_respreco_prev,&
       beta,grad_term,descent_scal_Hd
  real,dimension(:),allocatable :: mgrad
  
  !-----------------------------------------------------!
  ! continue the current conjugate gradient iteration   !
  !-----------------------------------------------------! 
  !Update CG direction
  res_scal_respreco_prev=optim%res_scal_respreco
  optim%res_scal_respreco = dot_product(residual,residual_preco)
  beta=(optim%res_scal_respreco)/(res_scal_respreco_prev)
  d(:)=-1.*residual_preco(:)+beta*d(:)                  
  !Update iteration counter 
  optim%cpt_iter_CG=optim%cpt_iter_CG+1  
  !-----------------------------------------------------!
  ! if the debug option is activated, compute the       !
  ! quadratic function minimized during the conjugate   !
  ! gradient process (check is this value decresae      !
  ! throughout the CG iterations )                      !
  !-----------------------------------------------------!  
  if(optim%debug) then
     optim%qkm1_CG=optim%qk_CG
     allocate(mgrad(n))
     mgrad(:)=-1.*grad(:)
     grad_term = dot_product(descent,mgrad)
     descent_scal_Hd= dot_product(descent_prev,Hd)
     optim%hessian_term=optim%hessian_term+(optim%alpha_CG**2)*optim%dHd+&
          2.*optim%alpha_CG*descent_scal_Hd
     optim%qk_CG=-grad_term+0.5*optim%hessian_term
     deallocate(mgrad)
  endif
  !-----------------------------------------------------!
  ! Check if the Eisenstat stopping critertion is       !
  ! satisfied                                           !
  !-----------------------------------------------------! 
  optim%norm_residual = norm2(residual)
  optim%conv_CG=(&
       (optim%norm_residual.le.(optim%eta*optim%norm_grad)).or.&
       (optim%cpt_iter_CG.ge.optim%niter_max_CG))        
  !-----------------------------------------------------!
  ! Print information on the current CG iteration       !
  !-----------------------------------------------------!
  call print_info_PTRN(optim,0e0,FLAG)     
  
end subroutine descent_PTRN2

!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! This routine is used to deallocate all the  !
! arrays that have been used during the       !
! PSTD optimization                           !
!---------------------------------------------!
! INPUT/OUT:  optim_type optim                !
!---------------------------------------------!
subroutine finalize_PTRN()
  
  implicit none   
  
  deallocate(xk)
  deallocate(descent)
  deallocate(descent_prev)
  deallocate(eisenvect)
  
end subroutine finalize_PTRN
!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This file contains the routine for the              !
! initialization of the truncated Newton algorithm    !
!------------------------------------------------------
! INPUT : integer n dimension                         !
!         real,dimension(n) x current point           !
!         real,dimension(n) grad current gradient     !
!         character*4 FLAG communication flag         !
! INPUT/OUTPUT : optim_type optim (data structure)    !
!-----------------------------------------------------!
subroutine init_PTRN(n,x,fcost,grad,optim)
  
  implicit none

  !IN
  integer :: n 
  real :: fcost
  real,dimension(n) :: x,grad
  !IN/OUT
  type(optim_type) :: optim !data structure   
    
  !---------------------------------------!
  ! set counters                          !
  !---------------------------------------!
  optim%cpt_iter=0
  optim%f0=fcost
  optim%nfwd_pb=0
  optim%nhess=0
  optim%eta=0.9  
  !optim%eta=0.1
  
  !---------------------------------------!
  ! initialize linesearch parameters      !
  ! by default, the max number of         !
  ! linesearch iteration is set to 20     !
  ! and the initial steplength is set to 1!
  !---------------------------------------! 
  optim%m1=1e-4  ! Wolfe conditions parameter 1 (Nocedal value)
  optim%m2=0.9   ! Wolfe conditions parameter 2 (Nocedal value)
  optim%mult_factor=10 !Bracketting parameter (Gilbert value)
  optim%fk=fcost
  optim%nls_max=20 ! max number of linesearch
  optim%cpt_ls=0
  optim%first_ls=.true.
  optim%alpha=1.   ! first value for the linesearch steplength

  !---------------------------------------!
  ! memory allocations                    !
  !---------------------------------------!
  allocate(xk(n))
  xk(:)=x(:)     
  allocate(descent(n))
  allocate(descent_prev(n))
  allocate(eisenvect(n))

  !---------------------------------------!
  ! norm of the first gradient            !
  !---------------------------------------!
  optim%norm_grad = norm2(grad)
   
  
end subroutine init_PTRN
!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This file contains the routine for printing the     !
! convergence history in the files iterate_PTRN.dat   !
! and iterate_PTRN_CG.dat                             !
! The file iterate_PTRN.dat is similar to the         !
! convergence history files of the other optimization !
! routines of the toolbox (PSTD,PNLCG,LBFGS,PLBFGS)   !
! The file iterate_PTRN_CG.dat contains additional    !
! information on the convergence of the inner         !
! preconditioned conjugate gradient iteration for the !
! computation of the inexact Newton descent direction.!
!------------------------------------------------------
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         character*4 FLAG communication flag         !
! INPUT/OUTPUT : optim_type optim (data structure)    !
!-----------------------------------------------------!
subroutine print_info_PTRN(optim,fcost,FLAG)

  implicit none

  !IN
  integer :: FLAG
  real :: fcost
  type(optim_type) :: optim !data structure   

  if(optim%print_flag.eq.1) then

     if(FLAG.eq. 0) then          
        open(10,file='iterate_PTRN.dat')     
        write(10,'(A90)') '******************************************************************************************'
        write(10,'(A70)') '                      PRECONDITIONED TRUNCATED NEWTON ALGORITHM                           '
        write(10,'(A90)') '******************************************************************************************'
        write(10,'(A30,ES10.2)') '     Convergence criterion  : ',optim%conv
        write(10,'(A30,I7)')     '     Niter_max              : ',optim%niter_max
        write(10,'(A30,ES10.2)') '     Initial cost is        : ',optim%f0
        write(10,'(A30,ES10.2)') '     Initial norm_grad is   : ',optim%norm_grad
        write(10,'(A30,I7)')     '     Maximum CG iter        : ',optim%niter_max_CG
        write(10,'(A90)') '******************************************************************************************'
        write(10,'(A10,A10,A10,A10,A12,A12,A8,A10,A9,A9)')&
             '   Niter   ' ,&
             '   fk      ' ,&
             '   ||gk||  ' ,&
             '     fk/f0   ' ,&
             '       alpha   ' ,&
             '        nls     ',&
             '  nit_CG',&
             '    eta     ',&
             '  ngrad ',&
             '  nhess '
        write(10,'(I6,ES12.2,ES12.2,ES12.2,ES12.2,I8,I8,ES12.2,I8,I8)') &
             optim%cpt_iter,&
             fcost,&
             optim%norm_grad,&
             fcost/optim%f0,&
             optim%alpha,&
             optim%cpt_ls,&
             optim%cpt_iter_CG,&
             optim%eta,&
             optim%nfwd_pb,&
             optim%nhess
        open(21,file='iterate_PTRN_CG.dat')     
        write(21,'(A90)') '******************************************************************************************'
        write(21,'(A70)') '                        PRECONDITIONED TRUNCATED NEWTON ALGORITHM                         '
        write(21,'(A70)') '                                      INNER CG HISTORY                                    '
        write(21,'(A90)') '******************************************************************************************'
        write(21,'(A30,ES10.2)') '     Convergence criterion  : ',optim%conv
        write(21,'(A30,I7)')     '     Niter_max              : ',optim%niter_max
        write(21,'(A30,ES10.2)') '     Initial cost is        : ',optim%f0
        write(21,'(A30,ES10.2)') '     Initial norm_grad is   : ',optim%norm_grad
        write(21,'(A30,I7)')     '     Maximum CG iter        : ',optim%niter_max_CG
        write(21,'(A90)') '******************************************************************************************'
     elseif(FLAG.eq.2) then
        write(10,'(A70)') '**********************************************************************'
        if(optim%cpt_iter.eq.optim%niter_max) then
           write(10,'(A50)') '  STOP: MAXIMUM NUMBER OF ITERATION REACHED    '
        else
           write(10,'(A50)') '  STOP: CONVERGENCE CRITERION SATISFIED        '
        endif
        write(10,'(A70)') '**********************************************************************'
        close(10)
     elseif(FLAG.eq. 4) then
        write(10,'(A70)') '**********************************************************************'
        write(10,'(A50)') '  STOP: LINESEARCH FAILURE    '
        write(10,'(A70)') '**********************************************************************'
        close(10)     
     elseif(optim%comm.eq. 1) then     
        write(21,'(A90)') '-------------------------------------------------------------------------------------------------'
        write(21,'(A20,I4,A10,ES12.2)') ' NONLINEAR ITERATION ',optim%cpt_iter, ' ETA IS : ',optim%eta
        write(21,'(A90)') '-------------------------------------------------------------------------------------------------'
        write(21,'(A12,A10,A10,A20)')&
             '  Iter_CG      ',&
             '  qk           ',&
             ' norm_res     ',&
             ' norm_res/||gk||  '
        write(21,'(I8,ES12.2,ES12.2,ES12.2)') &
             optim%cpt_iter_CG,&
             optim%qk_CG,&
             optim%norm_residual,&
             optim%norm_residual/optim%norm_grad       
     elseif(optim%comm.eq. 2) then
        write(21,'(I8,ES12.2,ES12.2,ES12.2)') &
             optim%cpt_iter_CG,&
             optim%qk_CG,&
             optim%norm_residual,&
             optim%norm_residual/optim%norm_grad
     else     
        write(10,'(I6,ES12.2,ES12.2,ES12.2,ES12.2,I8,I8,ES12.2,I8,I8)') &
             optim%cpt_iter,&          
             fcost,&
             optim%norm_grad,&
             fcost/optim%f0,&
             optim%alpha,&
             optim%cpt_ls,&
             optim%cpt_iter_CG,&
             optim%eta,&
             optim%nfwd_pb,&
             optim%nhess
     endif
  endif
end subroutine print_info_PTRN

!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routine is the reverse communication mode      !
! preconditioned truncated Newton.                    !
! This algorithm is described in                      !
! L.Metivier, R. Brossier, J. Virieux, S.Operto,      !
! Truncated Newton and full waveform inversion, 2013  !
! SIAM Journal on Scientific Computing, Vol. 35,      !
! No. 2, pp. B401–B437,                               !         
!                                                     !
! This routine performs an iterative                  !
! minimization of a function f following the          !
! recurrence                                          !
!                                                     !
! x_{0}=x0                                            !
! x_{k+1}=x_k+\alpha_k d_k  (1)                       !
!                                                     !
! where the descent direction d_k is computed through !
! the resolution of the preconditioned linear system  !
!                                                     !
! P_k H_k d_k=- P_k\nabla f_k  (2)                    !
!                                                     !
! with H_k        : Hessian operator at iteration k   !
!                                                     !
! with P_k        : Preconditioner                    !
!                                                     !
!\nabla f_k       : gradient of f in x_k              !
!                                                     !
! and alpha_k is the steplength computed through the  !
! common linesearch algorithm of the TOOLBOX          !
!                                                     !
! The linear system (2) is solved through a matrix    !
! free conjugate gradient algorithm which requires the!
! user to perform multiplication of given vector by   !
! the Hessian operator and the preconditioner.        !                     
!                                                     !
! Because of these tow nested algorithms, the reverse !
! communication strategy requires additional          !
! communicators within the code to clearly track which!
! stage the optimizer has reached                     !
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
! - FLAG='GRAD' => the user must compute the cost,the !
!                  gradient and the preconditioned    !
!                  gradient at current point x in     !
!                  fcost, grad, grad_preco            !
! - FLAG='HESS' => the user must multiply the vector  !
!                  d by the Hessian operator and!
!                  set the result in Hd         !
! - FLAG='PREC' => the user must multiply the vector  !
!                  residual by the              !
!                  preconditioner and set the result  !
!                  in residual_preco            !
! - FLAG='CONV' => a minimizer has been found         !
! - FLAG='NSTE' => a new step is performed            !
! - FLAG='FAIL' => the linesearch has failed          !
!-----------------------------------------------------!
! INPUT  : integer :: n (dimension)                   ! 
!          real fcost (current cost)                  !
!          real,dimension(n) grad                     !
! INPUT/OUTPUT : real,dimension(n) x                  !
!                optim_typ optim (data structure)     !
!                character*4 FLAG (communication)     !
!-----------------------------------------------------!
subroutine PTRN(n,x,fcost,grad,grad_preco,residual,residual_preco, &
				d,Hd,optim,FLAG,lb,ub)
  
  implicit none
  
  !IN
  integer  :: n                                           !dimension of the problem
  real :: fcost                                           !cost associated with x
  real,dimension(n) :: grad                               !gradient at x
  real,dimension(n) :: grad_preco                         !preconditioned gradient at x
  real,dimension(n) :: d, Hd 
  real,dimension(n) :: residual, residual_preco  
  !IN/OUT  
  integer:: FLAG  
  real,dimension(n) :: x                                  !current point
  type(optim_type) :: optim                               !data structure 
  real(c_float),optional, dimension(n) :: lb,ub 
  !Local variable
  logical :: test_conv
  
  if(FLAG.eq. 0) then
     !-----------------------------------------------------!
     ! if FLAG is INIT, call the dedicated initialization  !
     ! subroutine to allocate data structure optim         !
     !-----------------------------------------------------!
     call init_PTRN(n,x,fcost,grad,optim)     
     call print_info_PTRN(optim,fcost,FLAG)     
     optim%comm=1     
     optim%CG_phase=0
     optim%nfwd_pb=optim%nfwd_pb+1
     optim%conv_CG=.false.
     FLAG=6
  endif
  if(optim%comm.eq. 1) then  
     !-----------------------------------------------------!
     ! if optim%comm is DES1, the optimizer starts the     !
     ! computation of a descent direction through the      !
     ! conjugate gradient                                  !
     !-----------------------------------------------------!
     if(optim%CG_phase.eq. 0) then  
        !-----------------------------------------------------!
        ! if optim%CG_phase is INIT, initialization of the    !
        ! conjugate gradient process
        !-----------------------------------------------------!
        call descent_PTRN0(n,grad,grad_preco,residual,residual_preco,&
    					   d,Hd,optim,FLAG)  
        optim%CG_phase=1
        optim%comm=1
        FLAG=7
        optim%nhess=optim%nhess+1
     elseif(optim%CG_phase.eq. 1) then
        !-----------------------------------------------------!
        ! if optim%CG_phase is IRUN, iterate the conjugate    !
        ! gradient process (first part)                       !
        !-----------------------------------------------------!
        call descent_PTRN1(n,grad,residual,residual_preco,d,Hd,optim)    
        if(optim%conv_CG) then 
           !-----------------------------------------------------!
           ! if the conjugate gradient has already converged     ! 
           ! (detection of a negative curvature), go to next     !
           ! phase: linesearch in the descent direction          !
           !-----------------------------------------------------!
           optim%comm=3
           optim%CG_phase=0        
           FLAG=6                   
        else
           !-----------------------------------------------------!
           ! if the conjugate gradient has not converged         ! 
           ! prematurly, then ask the user to apply the          !
           ! preconditioner                                      !
           !-----------------------------------------------------!
           FLAG=8
           optim%comm=2
        endif
     endif
  elseif(optim%comm.eq. 2) then
     !-----------------------------------------------------!
     ! if optim%comm is DES2, the optimizer finish the     !
     ! current iteration of the conjugate gradient process !
     !-----------------------------------------------------!
     call descent_PTRN2(n,grad,residual,residual_preco,d,Hd,optim,FLAG)            
     if(optim%conv_CG) then
        !-----------------------------------------------------!
        ! if the conjugate gradient has converged go to next  !
        ! phase: linesearch in the descent direction          !
        !-----------------------------------------------------!        
        optim%comm=3
        optim%CG_phase=0        
        FLAG=6                
     else
        !-----------------------------------------------------!
        ! else start a new iteration of conjugate gradient and! 
        ! ask the user to compute a Hessian-vector product    !
        !-----------------------------------------------------!        
        optim%comm=1
        FLAG=7
        optim%nhess=optim%nhess+1
     endif
  elseif(optim%comm.eq. 3) then           
     !-----------------------------------------------------!
     ! if optim%comm is NSTE, a descent direction has been !
     ! computed, and a linesearch must be performed in this!
     ! direction                                           !
     !-----------------------------------------------------!
     call std_linesearch(n,x,fcost,grad,xk,descent,optim,lb,ub)     
     if(optim%task.eq. 0) then !NEW STEP                  
        !-----------------------------------------------------!
        ! if optim%task is 'NEW_STEP, the linesearch process  !
        ! has found the new step                              !
        !-----------------------------------------------------!
        optim%cpt_iter=optim%cpt_iter+1
        !Save the previous gradient norm 
        optim%norm_grad_m1=optim%norm_grad
        !Computation of the new gradient norm 
        optim%norm_grad = norm2(grad)
        !Print infor on current nonlinear iteration
        call print_info_PTRN(optim,fcost,FLAG)        
        !Test for convergence
        call std_test_conv(optim,fcost,test_conv)
        if(test_conv) then
           FLAG=2
           call print_info_PTRN(optim,fcost,FLAG)        
           close(10)
           close(21)
           call finalize_PTRN()
        else
           FLAG=3
           !Flag for the computation of the new descent direction
           optim%comm=1       
           !Update forcing term optim%eta following the Eisenstat and Walker formula
           call forcing_term_TRN(n,grad,residual,optim)  
        endif
     elseif(optim%task.eq.1) then 
        !-----------------------------------------------------!
        ! if optim%task is 'NEW_GRAD, the linesearch process  !
        ! is continuing, the gradient at the current point is !
        ! required                                            ! 
        !-----------------------------------------------------!
        FLAG=1         
        optim%nfwd_pb=optim%nfwd_pb+1
     elseif(optim%task.eq. 2) then        
        !-----------------------------------------------------!
        ! if optim%task is 'FAILURE, the linesearch process   !
        ! has failed, the iterations are stopped              !
        !-----------------------------------------------------!
        FLAG=4
        call print_info_PTRN(optim,fcost,FLAG)
        close(10)
        close(21)
        call finalize_PTRN()
     endif
  endif
  
end subroutine PTRN

end module opt_PTRN

