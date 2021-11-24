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

module opt_LBFGS
use typedef
use miscellaneous
use opt_PSTD, only : xk, descent

private
real, allocatable, dimension(:,:) :: sk,yk
public LBFGS, save_LBFGS, update_LBFGS, sk, yk

contains

subroutine lbfgs_centry(n,x,fcost,grad,optim,flag,lb,ub) bind(c, name='LBFGS')
  use, intrinsic :: iso_c_binding
  !IN
  integer(c_int), value  :: n                    !dimension of the problem
  real(c_float)   :: fcost                       !cost associated with x
  real(c_float),dimension(n) :: grad             !gradient at x 
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
  call lbfgs(n,x,fcost,grad,optim,flag,lb_pass,ub_pass)
end subroutine lbfgs_centry

!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routines is used by the l-BFGS routine to      !
! compute l-BFGS descent direction following the two  !
! recursion loop algorithm                            ! 
! Nocedal, 2nd edition, 2006 Algorithm 7.5 p. 179     !
!-----------------------------------------------------!
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         real,dimension(n,l) sk yk pairs of vectors  !
!                             for l-BFGS approximation!
!         integer cpt_lbfgs current number of stored  !
!                           pairs                     !
!         integer l maximum number of stored pairs    !
! OUTPUT : real,dimension(n) descent                  !
!-----------------------------------------------------!
subroutine descent_LBFGS(n,grad,sk,yk,cpt_lbfgs,l,descent)
  
  implicit none
  
  !IN
  integer :: n,cpt_lbfgs,l
  real,dimension(n) :: grad
  real,dimension(n,l) :: sk,yk
  !IN/OUT
  real,dimension(n) :: descent
  !Local variables
  real,dimension(:),allocatable :: q
  real,dimension(:),allocatable :: alpha,rho
  real :: beta,gamma,gamma_num,gamma_den
  real :: norml2_yk,norml2_sk
  integer :: i,borne_i
  
  borne_i=cpt_lbfgs-1
  
  !------------------------------------------!
  ! SAFEGUARD                                !
  !------------------------------------------!
  norml2_sk = norm2(sk(:,borne_i))
  norml2_yk = norm2(yk(:,borne_i))
  if( (norml2_sk==0.).or.(norml2_yk==0.)) then
     descent(:)=-1.*grad(:)     
  else
     !------------------------------------------!
     ! First phase of the recursion loop        !
     !------------------------------------------!
     allocate(alpha(cpt_lbfgs))
     allocate(rho(cpt_lbfgs))
     allocate(q(n))
     q(:)=grad(:)
     do i=1,borne_i
        rho(borne_i-i+1) = dot_product(yk(:,borne_i-i+1),sk(:,borne_i-i+1))
        rho(borne_i-i+1)=1./rho(borne_i-i+1)
        alpha(borne_i-i+1) = dot_product(sk(:,borne_i-i+1),q(:))
        alpha(borne_i-i+1)=rho(borne_i-i+1)*alpha(borne_i-i+1)
        q(:)=q(:)-alpha(borne_i-i+1)*yk(:,borne_i-i+1)
     enddo
     gamma_num = dot_product(sk(:,borne_i),yk(:,borne_i))
     gamma_den = norm2(yk(:,borne_i))
     !------------------------------------------!
     ! Scaling by gamma                         !
     !------------------------------------------!
     gamma=gamma_num/(gamma_den**2)     
     descent(:)=gamma*q(:)
     !------------------------------------------!
     ! Second phase of the recursion loop       !
     !------------------------------------------!
     do i=1,borne_i
        beta = dot_product(yk(:,i),descent(:))
        beta=rho(i)*beta
        descent(:)=descent(:)+(alpha(i)-beta)*sk(:,i)
     enddo
     descent(:)=-1.*descent(:)
     deallocate(q,alpha,rho)
  endif
  
end subroutine descent_LBFGS
!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! This routine is used to deallocate all the  !
! arrays that have been used during the       !
! PSTD optimization                           !
!---------------------------------------------!
! INPUT/OUT:  optim_type optim                !
!---------------------------------------------!
subroutine finalize_LBFGS()
  
  implicit none  
  
  deallocate(sk)
  deallocate(yk)
  deallocate(xk)
  !deallocate(optim%grad)
  deallocate(descent)
  
end subroutine finalize_LBFGS
!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routine is used for the initialization of the  !
! reverse communication mode l-BFGS algorithm from    !
! the TOOLBOX.                                        !
!                                                     !
! The parameters for the linesearch are set here, as  !
! well as the data structure allocations required for !
! the use of the algorithm                            !
!-----------------------------------------------------!
! INPUT : integer n, dimension of the problem         !
!       : real fcost                                  !
!       : real,dimension(n) x,grad first iterate and  !
!                           gradient                  !
! INPUT/OUTPUT : optim_type optim data structure      ! 
!-----------------------------------------------------!
subroutine init_LBFGS(n,l,x,fcost,grad,optim)
  
  implicit none

  !IN
  integer :: n,l 
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
  allocate(sk(n,l),yk(n,l))
  sk(:,:)=0.
  yk(:,:)=0.
  optim%cpt_lbfgs=1
  
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
  !allocate(optim%grad(n))
  !optim%grad(:)=grad(:)
  allocate(descent(n))

  !---------------------------------------!
  ! first descent direction               !
  !---------------------------------------!
  descent(:)=-1.*grad(:)  
  call save_LBFGS(n,optim%cpt_lbfgs,l,x,grad,sk,yk)
  
end subroutine init_LBFGS
!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routine is the reverse communication mode      !
! l-BFGS algorithm. This algorithm is described in    !
! Numerical Optimization, Nocedal, 2nd edition, 2006  !
! Algorithm 7.4 p. 178, Algorithm 7.5 p. 179          !
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
! d_k=-Q_k \nabla f_k                                 !
!                                                     !
! with Q_k        : l-BFGS approximation of the       !
!                   inverse Hessian at iteration k    !
!\nabla f_k       : gradient of f in x_k              !
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
! INPUT/OUTPUT : real,dimension(n) x                  !
!                optim_type optim (data structure)    !
!                character*4 FLAG (communication)     !
!-----------------------------------------------------!
subroutine LBFGS(n,x,fcost,grad,optim,flag,lb,ub) 
  use, intrinsic :: iso_c_binding
  implicit none
  
  !IN
  integer  :: n                         !dimension of the problem
  integer  :: l                         !number of stored pairs for the 
                                               !l-BFGS approximation of the 
                                               !inverse Hessian
  real :: fcost                         !cost associated with x
  real,dimension(n) :: grad             !associated gradient 
  !IN/OUT  
  integer  :: flag
  real,dimension(n) :: x                !current point
  type(optim_type) :: optim             !data structure 
  real(c_float),optional, dimension(n) :: lb,ub 
  !Local variables
  logical :: test_conv

  if(FLAG.eq. 0) then
     !-----------------------------------------------------!
     ! if FLAG is INIT, call the dedicated initialization  !
     ! subroutine to allocate data structure optim and     !
     ! initialize the linesearch process                   !
     !-----------------------------------------------------!
     call init_LBFGS(n,optim%l,x,fcost,grad,optim)
     call std_linesearch(n,x,fcost,grad,xk,descent,optim,lb,ub)
     call print_info(n,'LB',optim,fcost,grad,FLAG)
     FLAG=1    
     optim%nfwd_pb=optim%nfwd_pb+1
  else
     !-----------------------------------------------------!
     ! else call the linesearch process                    !  
     !-----------------------------------------------------!
     call std_linesearch(n,x,fcost,grad,xk,descent,optim,lb,ub)     
     if(optim%task.eq. 0) then
        !-----------------------------------------------------!
        ! test for convergence                                !
        !-----------------------------------------------------!
        optim%cpt_iter=optim%cpt_iter+1                   
        call std_test_conv(optim,fcost,test_conv)
        if(test_conv) then
           FLAG=2           
           !print info on current iteration
           call print_info(n,'LB',optim,fcost,grad,FLAG)        
           call finalize_LBFGS()
        else
           FLAG=3
           !-----------------------------------------------------!
           ! if a NEW_STEP is taken, compute a new descent       !
           ! direction using current gradient and l-BFGS         !
           ! approximation of the inverse Hessian                !
           ! preconditioned gradient                             !
           !-----------------------------------------------------!
           !LBFGS update 
           call update_LBFGS(n,optim%cpt_lbfgs,optim%l,x,grad,sk,yk)
           !Computation of the new descent direction
           call descent_LBFGS(n,grad,sk,yk,optim%cpt_lbfgs,&
                l,descent)         
           !LBFGS store
           call save_LBFGS(n,optim%cpt_lbfgs,optim%l,x,grad,sk,yk)           
           !print info on current iteration
           call print_info(n,'LB',optim,fcost,grad,FLAG)        
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
        call print_info(n,'LB',optim,fcost,grad,FLAG)
        call finalize_LBFGS
     endif
  endif
  
end subroutine LBFGS



!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routines is used by the l-BFGS routine to      !
! save the pairs of models and gradient for        !
! the inverse Hessian approximation                   !
! Nocedal, 2nd edition, 2006 Algorithm 7.5 p. 179     !
!-----------------------------------------------------!
! INPUT : integer n dimension                         !
!         integer cpt_lbfgs counts the number of      ! 
!                           already stored pairs      !
!         integer l maximum number of stored pairs    !
!         real,dimension(n) x current point           !
!         real,dimension(n) x current gradient        !
! OUTPUT : sk,yk updated vector                       !
!-----------------------------------------------------!
subroutine save_LBFGS(n,cpt_lbfgs,l,x,grad,sk,yk)
  
  implicit none
  
  !IN
  integer :: n,cpt_lbfgs,l
  real,dimension(n) :: x,grad
  !IN/OUT
  real,dimension(n,l) :: sk,yk
  
  !Local variables
  integer :: i
  
  if(cpt_lbfgs.le.l) then
     !---------------------------------------------------!
     ! if the number of stored pairs does not exceed the !
     ! maximum value, then save x and grad               !
     !---------------------------------------------------!
     sk(:,cpt_lbfgs)=x(:)
     yk(:,cpt_lbfgs)=grad(:)   
  else
     !---------------------------------------------------!
     ! otherwise, erase the oldest pair and save the     !
     ! new one (shift)                                   !
     !---------------------------------------------------!
     do i=1,l-1
        sk(:,i)=sk(:,i+1)
        yk(:,i)=yk(:,i+1)                
     enddo
     sk(:,l)=x(:)
     yk(:,l)=grad(:)
  endif
  
end subroutine save_LBFGS
!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routines is used by the l-BFGS routine to      !
! compute the new pairs of models and gradient for    !
! the inverse Hessian approximation                   !
! Nocedal, 2nd edition, 2006 Algorithm 7.5 p. 179     !
!-----------------------------------------------------!
! INPUT : integer n dimension                         !
!         integer cpt_lbfgs counts the number of      ! 
!                           already stored pairs      !
!         integer l maximum number of stored pairs    !
!         real,dimension(n) x current point           !
!         real,dimension(n) x current gradient        !
! OUTPUT : sk,yk updated vectors                      !
!-----------------------------------------------------!
subroutine update_LBFGS(n,cpt_lbfgs,l,x,grad,sk,yk)

  implicit none

  !IN
  integer :: n,cpt_lbfgs,l
  real,dimension(n) :: x,grad
  !IN/OUT
  real,dimension(n,l) :: sk,yk
  
  if(cpt_lbfgs.le.l) then
     !---------------------------------------------------!
     ! if the number of stored pairs does not exceed the !
     ! maximum value, then compute a new pair sk yk and  !
     ! update the counter cpt_lbfgs                      !
     !---------------------------------------------------!
     sk(:,cpt_lbfgs)=x(:)-sk(:,cpt_lbfgs)
     yk(:,cpt_lbfgs)=grad(:)-yk(:,cpt_lbfgs)
     cpt_lbfgs=cpt_lbfgs+1
  else
     !---------------------------------------------------!
     ! otherwise, simply update the lth pair             !
     !---------------------------------------------------!
     sk(:,l)=x(:)-sk(:,l)
     yk(:,l)=grad(:)-yk(:,l)
  endif
  
end subroutine update_LBFGS

end module opt_LBFGS
