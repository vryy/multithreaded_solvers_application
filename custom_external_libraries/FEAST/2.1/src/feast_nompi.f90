  !=========================================================================================
  !Copyright (c) 2009-2013, The Regents of the University of Massachusetts, Amherst.
  !Developed by E. Polizzi
  !All rights reserved.
  !
  !Redistribution and use in source and binary forms, with or without modification, 
  !are permitted provided that the following conditions are met:
  !
  !1. Redistributions of source code must retain the above copyright notice, this list of conditions 
  !   and the following disclaimer.
  !2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
  !   and the following disclaimer in the documentation and/or other materials provided with the distribution.
  !3. Neither the name of the University nor the names of its contributors may be used to endorse or promote
  !    products derived from this software without specific prior written permission.
  !
  !THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
  !BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
  !ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
  !EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
  !SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
  !LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
  !IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  !==========================================================================================
   
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! FEAST KERNEL - REVERSE COMMUNICATION INTERFACES !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! List of routines:
  !-------------------
  ! feastinit
  ! feastinit_driver
  ! check_fpm_input
  ! scheck_rci_input
  ! dcheck_rci_input
  ! scheck_rci_input_z
  ! dcheck_rci_input_z
  ! dfeast_srci
  ! zfeast_hrci
  ! sfeast_srci
  ! cfeast_hrci


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  include 'point_gauss_legendre.f90' 


  subroutine feastinit(fpm)
    !  Purpose
    !  =======
    !
    !  Define the default values for the input FEAST parameters.
    !
    !  Arguments
    !  =========
    !
    !  fpm  (output) INTEGER(*): FEAST parameters (size array at least 64)
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================
    implicit none
!    !-------------------------------------
!#ifdef MPI
!    include 'mpif.h'
!#endif
!    !-------------------------------------
    integer,dimension(*) :: fpm
    integer :: i
    do i=1,64
    fpm(i)=0
    enddo 

    fpm(1)=0 !com
    fpm(2)=8 !nbe
    fpm(3)=12!tol double precision
    fpm(4)=20!maxloop
    fpm(5)=0 !IS
    fpm(6)=0 !Convergence criteria
    fpm(7)=5 !tol single precision
    fpm(11)=1 !used for feast_sparse (0: no transpose-<MKL v10.3, 1: transpose capabilities)--undocumented   
    fpm(14)=0 ! return only subspace Q size M0 ==> 1 contour 
!!---------------------------------------
!#ifdef MPI
!    fpm(9)=MPI_COMM_WORLD ! default value
!#endif
!!----------------------------------------
    fpm(12)=0 ! customize eigenvalue solver for rci interface (No,Yes) (--undocumented--)
    fpm(64)=0 ! Additional feast parameters for driver interfaces (i.e size fpm>64) (0,1)      

 end subroutine feastinit



  subroutine feastinit_driver(fpm,N)
    !  Purpose
    !  =======
    !
    !  Define the default values for the input FEAST parameters from 1-64
    !  and 65-N are optional user defined inputs that can be used for FEAST predefined interfaces,
    !  here fpm(65:N) is initialized at -9999.
    !
    !  Arguments
    !  =========
    !
    !  fpm  (output) INTEGER(*): FEAST parameters (size array at least 64)
    !  N    (input)  INTEGER: size of array fpm (>=64)
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================
    implicit none
!    !-------------------------------------
!#ifdef MPI
!    include 'mpif.h'
!#endif
!    !-------------------------------------
    integer,dimension(*) :: fpm
    integer :: N
    integer :: i
    do i=1,64
    fpm(i)=0
    enddo 

    fpm(1)=0 !com
    fpm(2)=8 !nbe
    fpm(3)=12!tol double precision
    fpm(4)=20!maxloop
    fpm(5)=0 !IS
    fpm(6)=0 !Convergence criteria
    fpm(7)=5 !tol single precision
    fpm(11)=1 !used for feast_sparse (0: no transpose-<MKL v10.3, 1: transpose capabilities)--undocumented   
    fpm(14)=0 ! return only subspace Q size M0 ==> 1 contour 
!!---------------------------------------
!#ifdef MPI
!    fpm(9)=MPI_COMM_WORLD ! default value
!#endif
!!----------------------------------------  
    fpm(12)=0 ! customize eigenvalue solver for rci interface (No,Yes) (--undocumented--)
    fpm(64)=0 ! Additional feast parameters for driver interfaces (i.e size fpm>64) (0,1)      
    if (N>64) then
    fpm(64)=1
    do i=65,N
    fpm(i)=-9999 ! default values for the drivers inputs
    enddo
   end if

 end subroutine feastinit_driver



  subroutine check_fpm_input(fp,info)
    !  Purpose 
    !  =======
    !  Error handling for input FEAST parameters.
    !  Check the values for the input FEAST parameters, and return 
    !  info code error /=0 if incorrect values are found
    !
    !  Arguments
    !  =========
    !
    !  fp   (input) INTEGER(*) : FEAST parameters
    !  info (input/output) INTEGER
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================
    implicit none
    integer,dimension(*) :: fp
    integer :: info
    integer:: i
    logical :: test
    integer,parameter :: max=13
    integer, dimension(max):: tnbe=(/3,4,5,6,8,10,12,16,20,24,32,40,48/)


    if ((fp(1)/=0).and.(fp(1)/=1)) info=101
    test=.true.
    do i=1,max
       if (fp(2)==tnbe(i)) test=.false. 
    enddo
    if (test) info=102
    if (fp(3)<0) info=103
    if (fp(4)<0) info=104
    if ((fp(5)/=0).and.(fp(5)/=1)) info=105
    if ((fp(6)/=0).and.(fp(6)/=1)) info=106
  ! if ((fp(11)/=1).and.(fp(11)/=2)) info=111
    if ((fp(12)/=0).and.(fp(12)/=1)) info=112

  end subroutine check_fpm_input




  subroutine dcheck_rci_input(Emin,Emax,M0,N,info)
    !  Purpose 
    !  =======
    !  Error handling for the FEAST RCI double precision interfaces input parameters.
    !  Check the values of Emin, Emax, M0,N, and return 
    !  info code error /=0 if incorrect values are found
    !
    !  Arguments
    !  =========
    !
    !  Emin,Emax   (input) REAL DOUBLE PRECISION: search interval
    !  M0          (input) INTEGER: Size subspace
    !  N           (input) INTEGER: Size system
    !  info (input/output) INTEGER
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================
    implicit none
    double precision :: Emin,Emax
    integer :: N,M0,info

    if (Emin>=Emax) info=200 ! problem with Emin, Emax
    if ((M0<=0).or.(M0>N)) info=201 ! problem with M0 
    if (N<=0) info=202 ! problem with N

  end subroutine dcheck_rci_input





  subroutine scheck_rci_input(Emin,Emax,M0,N,info)
    !  Purpose 
    !  =======
    !  Error handling for the FEAST RCI single precision interfaces input parameters.
    !  Check the values of Emin, Emax, M0,N, and return 
    !  info code error /=0 if incorrect values are found
    !
    !  Arguments
    !  =========
    !
    !  Emin,Emax   (input) REAL SINGLE PRECISION: search interval
    !  M0          (input) INTEGER: Size subspace
    !  N           (input) INTEGER: Size system
    !  info (input/output) INTEGER
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================
    implicit none
    real :: Emin,Emax
    integer :: N,M0,info

    if (Emin>=Emax) info=200 ! problem with Emin, Emax
    if ((M0<=0).or.(M0>N)) info=201 ! problem with M0 
    if (N<=0) info=202 ! problem with N

  end subroutine scheck_rci_input




  subroutine dcheck_rci_input_z(Emid,r,M0,N,info)
    !  Purpose 
    !  =======
    !  Error handling for the FEAST RCI double precision interfaces input parameters.
    !  Check the values of Emid, r, M0,N, and return 
    !  info code error /=0 if incorrect values are found
    !
    !  Arguments
    !  =========
    !
    !  Emid,r      (input) REAL DOUBLE PRECISION: search interval (center+radius)
    !  M0          (input) INTEGER: Size subspace
    !  N           (input) INTEGER: Size system
    !  info (input/output) INTEGER
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================
    implicit none
    double precision :: Emid,r
    integer :: N,M0,info

    if (r<=0.0d0) info=200 ! problem with r
    if ((M0<=0).or.(M0>N)) info=201 ! problem with M0 
    if (N<=0) info=202 ! problem with N

  end subroutine dcheck_rci_input_z


  subroutine scheck_rci_input_z(Emid,r,M0,N,info)
    !  Purpose 
    !  =======
    !  Error handling for the FEAST RCI real interfaces input parameters.
    !  Check the values of Emid, r, M0,N, and return 
    !  info code error /=0 if incorrect values are found
    !
    !  Arguments
    !  =========
    !
    !  Emid,r      (input) REAL DOUBLE PRECISION: search interval (center+radius)
    !  M0          (input) INTEGER: Size subspace
    !  N           (input) INTEGER: Size system
    !  info (input/output) INTEGER
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================
    implicit none
    real :: Emid,r
    integer :: N,M0,info

    if (r<=0.0d0) info=200 ! problem with r
    if ((M0<=0).or.(M0>N)) info=201 ! problem with M0 
    if (N<=0) info=202 ! problem with N

  end subroutine scheck_rci_input_z


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  subroutine dfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST RCI (Reverse Communication Interfaces) 
    !  Solve generalized Aq=lambda Bq eigenvalue problems
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE (or B Identity) 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  ijob       (input/output) INTEGER :: ID of the RCI
    !                            INPUT on first entry: ijob=-1 
    !                            OUTPUT Return values (0,10,20,21,30,40)-- see FEAST documentation
    !  N          (input)        INTEGER: Size system
    !  work       (input/output) REAL DOUBLE PRECISION (N,M0)   :  Workspace 
    !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):  Workspace 
    !  Aq,Sq      (input/output) REAL DOUBLE PRECISION (M0,M0)  :  Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm     (input/output) INTEGER(*) : FEAST parameters (see FEAST documentation)
    !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) REAL DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                               
    !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================

    implicit none
!    !-------------------------------------
!#ifdef MPI
!    include 'mpif.h'
!#endif
!    !-------------------------------------
    include "f90_noruntime_interface.fi"
    integer :: ijob,N,M0
    complex(kind=(kind(1.0d0))) :: Ze
    double precision, dimension(N,*) ::work
    complex(kind=(kind(1.0d0))),dimension(N,*):: workc
    double precision,dimension(M0,*):: Aq,Sq
    integer,dimension(*) :: fpm
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    double precision,dimension(*)  :: lambda
    double precision,dimension(N,*):: q
    integer :: mode
    double precision,dimension(*) :: res
    integer :: info
    !! parameters
    double precision, Parameter :: pi=3.1415926535897932d0
    double precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO), ZEROC=(DZERO,DZERO)
    double precision, parameter :: ba=-pi/2.0d0, ab=pi/2.0d0
    integer(8),parameter :: fout =6
    !! variable for FEAST
    integer :: i,e,j,k
    integer,dimension(4) :: iseed
    double precision :: theta,r,Emid
    complex(kind=(kind(1.0d0))) :: jac
    double precision ::xe,we ! Gauss-Legendre
    double precision, dimension(:,:),pointer :: Sqo
    integer, dimension(:),pointer :: fpm_default
    logical :: testconv
    double precision :: trace
    !! Lapack variable (reduced system)
    character(len=1) :: JOBZ,UPLO
    double precision, dimension(:),pointer :: work_loc,work_loc2
    integer :: lwork_loc,info_lap,infoloc
    !! MPI compatibility variables
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
   
    rank=0
    nb_procs=1
!    !----------------------------------------------
!#ifdef MPI
!    NEW_COMM_WORLD=fpm(9)
!    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
!    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
!    if (rank/=0) fpm(1)=0 ! comment only in rank 0 if any
!#endif
!    !---------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialization!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (ijob==-1) then 
       info=0 ! default value
       if (fpm(1)==1) then
          call wwrite_n(fout)
          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_s(fout, '*********** FEAST- BEGIN **********************')
          call wwrite_n(fout) 
          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout)
          call wwrite_s(fout, 'Routine DFEAST_S{}{}')
          call wwrite_n(fout)
!!!!!!!!!!!! Print the FEAST parameters which has been changed from default
          call wwrite_s(fout, 'List of input parameters fpm(1:64)-- if different from default')
          call wwrite_n(fout)

          call wallocate_1i(fpm_default,64,infoloc)
          call feastinit(fpm_default)
          do i=1,20
             if (fpm(i)/=fpm_default(i)) then
                call wwrite_s(fout, '   fpm(')
                call wwrite_i(fout, i)
                call wwrite_s(fout, ')=')
                call wwrite_i(fout, fpm(i))
                call wwrite_n(fout)
             endif
          enddo
          call wdeallocate_1i(fpm_default)

          call wwrite_s(fout, 'Search interval [')
          call wwrite_d(fout,Emin)
          call wwrite_s(fout, '; ')
          call wwrite_d(fout,Emax)
          call wwrite_s(fout, ']')
          call wwrite_n(fout)
       end if
       call check_fpm_input(fpm,info)
       call dcheck_rci_input(Emin,Emax,M0,N,info)

       if (info/=0) fpm(21)=100 ! The End

!!!!!!!!!!!!!!!!
       IF (info==0) then
          if (fpm(1)==1) then
             call wwrite_s(fout, 'Size subspace')  
             call wwrite_t(fout) 
             call wwrite_i(fout,M0)
             call wwrite_n(fout)
             call wwrite_s(fout, '#Loop | #Eig  |       Trace           |     Error-Trace       |     Max-Residual')  
             call wwrite_n(fout)
          endif
          fpm(23)=min(M0,N) ! 'current M0' size (global value)
          fpm(25)=fpm(23) !! 'current M0' size (by default)
          fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
!          !----------------------------------------------
!#ifdef MPI
!          if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
!             fpm(25)=fpm(23)/nb_procs ! local size of 'current M0'
!             if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
!             fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
!          end if
!#endif
!          !-----------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          fpm(22)=fpm(2) ! only one half contour necessary here
          loop=0
          fpm(21)=1 ! prepare reentry
          if (fpm(5)==0) then !!! random vectors
             iseed=(/56,890,3456,2333/)
             call DLARNV(3,iseed,N*fpm(23),work(1,1))
          elseif (fpm(5)==1) then !!!!!! q is the initial guess
!             !----------------------------------------
!#ifdef MPI
!             work(1:N,1:fpm(23))=DZERO 
!#endif
!             !------------------------------------------
             ijob=40 !! B*q=>work
             return
          end if
       end IF ! info=0
!!!!!!!!!!!!!!
    end if   !ijob=-1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    IF (fpm(21)==1) then !! we initialize a new contour integration
!       !------------------------------------------------------------------------
!#ifdef MPI
!       if ((loop>0).or.(fpm(5)==1)) then        
!          if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work,N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
!       end if
!#endif
!       !------------------------------------------------------------------------
       q(1:N,1:fpm(23))=DZERO
       fpm(20)=1
       fpm(21)=2
       ijob=-2 ! just initialization 
    end IF


!!!!!!!!!!!!
    IF (fpm(21)==2) then !! we start or pursue the contour integration

       IF (info==0) then !! will end up checking info errors returned by FEAST drivers
          do e=fpm(20)+rank,fpm(22),nb_procs !!!! loop over the contour 

             if (ijob==-2) then !!Factorize the linear system (complex) (zS-A)
                call dset_point_gauss_legendre(fpm(2),e,xe,we) !! Gauss-points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0d0
                Emid=Emin+r
                Ze=Emid*ONEC+r*ONEC*wdcos(theta)+r*(DZERO,DONE)*wdsin(theta)                               
                fpm(20)=e-rank
                ijob=10 ! for fact
                if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor 
             endif

             if (ijob==10) then !!Solve the linear system (complex) (zS-A)q=v 
                call ZLACP2( 'F', N, fpm(23),work , N, workc, N )
                ijob=11 ! for solve
                return
             endif

             if (ijob==11) then 
                !! summation 
                call dset_point_gauss_legendre(fpm(2),e,xe,we) !! Gauss-points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0d0
                jac=(r*(DZERO,DONE)*wdsin(theta)+ONEC*r*wdcos(theta))
                q(1:N,1:fpm(23))=q(1:N,1:fpm(23))-(DONE/2.0d0)*we*dble(jac*workc(1:N,1:fpm(23)))                  
                ijob=-2 ! just for identification
             end if
          end do
       end IF  !! info=0

!       !------------------------------------------------
!#ifdef MPI
!       call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!#endif
!       !-----------------------------------------------  
       if (info/=0) fpm(21)=100 ! the end

       if (info==0) then
          fpm(21)=4 
!          !------------------------------------------------
!#ifdef MPI
!          call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
!#endif
!          !-----------------------------------------------  
       end if
    end IF     ! fpm(21)==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ((fpm(21)==4).and.(fpm(14)==1)) then !! only q vectors has been computed and is returned
       info=4
       if (info/=0) fpm(21)=100 ! The End
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq=eq Sq xq     with Aq=Q^TAQ Sq=Q^TAQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!form Aq=> Aq=Q^T A Q 
    if (fpm(21)==4) then
       fpm(21)=5 ! preparing reentry
       ijob=30 
       return  ! mat-vec A*q => work
    endif

    if (fpm(21)==5) then 
!       !------------------------------------------------
!#ifdef MPI 
!       Aq(1:M0,1:fpm(23))=DZERO
!#endif
!       !-------------------------------------------------
       call DGEMM('T','N',fpm(23),fpm(25),N,DONE,q(1,1),N,work(1,fpm(24)),N,DZERO,Aq(1,fpm(24)),M0)
       fpm(21)=6
    endif

!!!!!!!!!form  Sq=> Sq=Q^T S Q
    if (fpm(21)==6) then 
       fpm(21)=7 ! preparing reenty
       ijob=40 
       return! mat-vec S*q => work
    end if
    if (fpm(21)==7) then
!       !------------------------------------------------
!#ifdef MPI
!       Sq(1:M0,1:fpm(23))=DZERO
!#endif
!       !-------------------------------------------------
       call DGEMM('T','N',fpm(23),fpm(25),N,DONE,q(1,1),N,work(1,fpm(24)),N,DZERO,Sq(1,fpm(24)),M0)
       fpm(21)=8
    endif

    if (fpm(21)==8) then
!       !---------------------------------------- !(Aq,Sq known to all processors) 
!#ifdef MPI
!       if (fpm(23)/nb_procs>=1) then
!          call MPI_ALLREDUCE(MPI_IN_PLACE,Aq(1,1),M0*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
!          call MPI_ALLREDUCE(MPI_IN_PLACE,Sq(1,1),M0*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
!       end if
!#endif
!       !---------------------------------------
       if (fpm(12)==1) then ! customize eigenvalue solver
          fpm(21)=9 ! preparing reentry - could return new value of M0 in fpm(23) if reduced subspace is needed
          ijob=50
          return
       endif
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==8) then

       if (rank==0) then ! if using FEAST-MPI ==> solve on a single proc
          JOBZ='V'
          UPLO='L'
          info_lap=1 ! initialization
          i=1
          LWORK_LOC=3*fpm(23)-1 !! for lapack eig reduced system
          call wallocate_1d(WORK_LOC,LWORK_LOC,infoloc)
          if (infoloc/=0) info=-1
          do while ((info_lap/=0).and.(info==0))
             i=i+1
             if (i==10) info=-3 ! arbitrary maximum
             call wallocate_2d(Sqo,fpm(23),fpm(23),infoloc)
             if (infoloc/=0) info=-1

             ! # eigenvalues inside search (wait second contour)
             mode=0
             if (loop>0) then
                call DLACPY( 'F', fpm(23), fpm(23),Sq , M0, Sqo, fpm(23) )
                if (infoloc/=0) info=-1
                call DSYEV( 'N', UPLO, fpm(23), Sqo, fpm(23), lambda, work_loc, lwork_loc, INFO_lap )
                do j=1,fpm(23)
                   if (sqrt(abs(lambda(j)))>=0.5d0) mode=mode+1
                enddo
             end if

             call DLACPY( 'F', fpm(23), fpm(23),Sq , M0, Sqo, fpm(23) )
             call DSYGV(1, JOBZ, UPLO, fpm(23),  Aq, M0,  Sqo, fpm(23), lambda,work_loc,Lwork_loc,INFO_lap)

             if ((info_lap<=fpm(23)).and.(info_lap/=0)) info=-3
             if (info_lap>fpm(23)) then !! zSqo is not spd (a posteriori resize subspace)
                fpm(23)=info_lap-fpm(23)-1
                if (fpm(1)==1) then
                   call wwrite_s(fout, 'Resize subspace')  
                   call wwrite_t(fout) 
                   call wwrite_i(fout,fpm(23))
                   call wwrite_n(fout)
                end if
             end if
             call wdeallocate_2d(Sqo)
          end do
          call wdeallocate_1d(work_loc)
       end if !(rank 0)
!       !-------------------------------- !(info common to all processors)
!#ifdef MPI
!       call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!#endif
!       !--------------------------------
       if (info/=0) fpm(21)=100 ! the end

       if (info==0) then
          fpm(25)=fpm(23) !! current M0 size (by default) -- global
          fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
!          !----------------------------------------!(Aq==> vectors, lambda and fpm(23), known to all processors) 
!#ifdef MPI 
!          call MPI_BCAST(mode,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!          call MPI_BCAST(fpm(23),1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!          call MPI_BCAST(Aq,M0*fpm(23),MPI_DOUBLE_PRECISION,0,NEW_COMM_WORLD,code)
!          call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_PRECISION,0,NEW_COMM_WORLD,code)
!          if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
!             fpm(25)=fpm(23)/nb_procs ! local size of current M0
!             if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
!             fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
!          end if
!#endif
!          !-----------------------------
          fpm(21)=9
       end if
    end if !! fpm(21)=8


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==9) then
       call DLACPY( 'F', N, fpm(23),q , N, work, N ) 
       !! option - non shifted
!#ifdef MPI 
!          q(1:N,1:fpm(23))=DZERO
!#endif
          call DGEMM('N','N',N,fpm(25),fpm(23),DONE,work(1,1),N,Aq(1,fpm(24)),M0,DZERO,q(1,fpm(24)),N)
!#ifdef MPI       
!          if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
!#endif
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residual ||AX-lambda*BX||
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==9) then
       fpm(21)=10 ! preparing reentry
       ijob=30 
       return  ! mat-vec A*q => work
    endif

    if (fpm(21)==10) then
       call ZLACP2( 'F', N, fpm(25),work(1,fpm(24)), N, workc(1,fpm(24)), N )
       fpm(21)=11 ! preparing reentry
       ijob=40 
       !----------------------------------------
!#ifdef MPI
!       work(1:N,1:fpm(23))=DZERO  !! work is also needed for outer-loop if any 
!#endif
       !------------------------------------------
       return  ! mat-vec S*q => work 
    endif

    if (fpm(21)==11) then
       !----------------------------------------
!#ifdef MPI
!       res(1:fpm(23))=DZERO
!#endif
       !-------- Absolute residual
       do i=fpm(24),fpm(24)+fpm(25)-1
          res(i)=sum(abs(dble(workc(1:N,i))-lambda(i)*work(1:N,i)))/sum(abs(max(abs(Emin),abs(Emax))*(work(1:N,i))))
       end do
       !----------------------------------------
!#ifdef MPI 
!       if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
!#endif
       !-----------------------------------------
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Trace / count eigenvalues (remove spurious if any)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (fpm(21)==11) then
       k=0
       trace=DZERO
       theta=DZERO
       do i=1,fpm(23)
          if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
             k=k+1
             trace=trace+lambda(i)
             if (res(i)>theta) theta=res(i) ! max residual
          endif
       enddo

!!!!!!!!! remove spurious if any       
       if (loop==0) mode=k ! wait second contour

       if (mode<k) then
          do j=1,k-mode
             theta=DZERO
             e=1
             do i=1,fpm(23)
                if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
                   if (res(i)>theta) then
                      e=i
                      theta=res(i) !max
                   endif
                end if
             enddo
             trace=trace-lambda(e) ! update trace
             res(e)=-DONE !spurious
          end do
!!! max residual
          theta=DZERO
          do i=1,fpm(23)
             if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
                if (res(i)>theta) theta=res(i) 
             end if
          enddo
       end if

       if (mode==0) info=1  ! no mode detected in the interval
if (loop>0) then
   if ((mode==M0).and.(mode/=N)) info=3 ! size subspace too small
endif
       if (info/=0) fpm(21)=100 ! The End
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check FEAST Convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

    if (fpm(21)==11) then
       testconv=.false. ! initialization

       !! trace convergence
       if (loop==0) epsout=DONE
       if (loop>0) then
          epsout=(abs(trace-epsout))/max(abs(Emin),abs(Emax))           
          if ((fpm(6)==0).and.(log10(epsout)<(-fpm(3)))) testconv=.true.
       end if

       !! residual convergence
       if ((fpm(6)==1).and.(log10(theta)<(-fpm(3)))) testconv=.true.

       !! Two loops minimum if spurious are found
       if ((loop<=1).and.(k>mode)) testconv=.false.

       if (fpm(1)==1) then
          call wwrite_i(fout,loop)
          call wwrite_t(fout) 
          call wwrite_i(fout,mode)
          call wwrite_t(fout)
          call wwrite_d(fout,trace)
          call wwrite_t(fout) 
          call wwrite_d(fout,epsout)
          call wwrite_t(fout) 
          call wwrite_d(fout,theta)
          call wwrite_n(fout) 
       end if

       if (.not.testconv) then
          epsout=trace
          if (loop==fpm(4)) then
             info=2 ! FEAST did not converge (#loop reaches maximum)
             testconv=.true. ! return final eigenvector anyway
          endif
       endif
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! FEAST exit IF Convergence - FEAST iteration IF NOT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

    if (fpm(21)==11) then

       if (.not.testconv) then !!! need FEAST refinement loop

!!!!!!!!!!! here work=S*q is known alreasy,  using FEAST-MPI work is already partionned 
          fpm(21)=1   ! prepare reentry- reloop (with contour integration)
          !fpm(21)=4 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
          loop=loop+1
          ijob=-2 ! do nothing
          return  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       else    !!!!!!! final eigenvectors/eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! reorder (shift) lambda, eigenvector and residual
   
          call DLACPY( 'F', N, fpm(23),q , N, work, N )
          call wallocate_1d(work_loc,fpm(23),infoloc)
          call wallocate_1d(work_loc2,fpm(23),infoloc)
          call DCOPY(fpm(23),lambda , 1, work_loc, 1 )
          call DCOPY(fpm(23),res , 1, work_loc2, 1 )
          q(1:N,1:fpm(23))=DZERO
          e=0
          j=0
          k=0
          do i=1,fpm(23)
             if ((work_loc(i)>Emin).and.(work_loc(i)<Emax)) then ! inside the search interval
                if (work_loc2(i)/=-DONE) then ! spurious 
                   e=e+1
                   q(1:N,e)=work(1:N,i)
                   lambda(e)=work_loc(i)
                   res(e)=work_loc2(i)
                endif
             else
                j=j+1
                q(1:N,mode+j)=work(1:N,i)
                lambda(mode+j)=work_loc(i)
                res(mode+j)=work_loc2(i)
             end if
             if (work_loc2(i)==-DONE) then ! spurious at the end
                k=k+1
                q(1:N,fpm(23)-k+1)=work(1:N,i)
                lambda(fpm(23)-k+1)=work_loc(i)
                res(fpm(23)-k+1)=work_loc2(i)
             endif
          enddo
          call wdeallocate_1d(work_loc)
          call wdeallocate_1d(work_loc2)

!!!!!!!!!!!!!!!!!!!!!!!!!!
          M0=fpm(23)  ! update value of M0 (new subspace)
          fpm(21)=100 ! The End

       end if ! test convergence
    end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==100) then !! THE END (ijob=0) 
       ijob=0 !! exit FEAST

       if (fpm(1)==1) then !! Print  Information

          if (info>=200) then
             call wwrite_s(fout, 'PROBLEM with input parameters')
             call wwrite_n(fout) 
          end if

          if ((info>100).and.(info<200)) then
             call wwrite_s(fout, 'PROBLEM with FEAST array parameters')
             call wwrite_n(fout) 
          end if

          if (info==-3) then
             call wwrite_s(fout, 'ERROR with reduced system')  
             call wwrite_n(fout) 
          end if

          if (info==-2) then
             call wwrite_s(fout, 'ERROR from Inner Linear System Solver in FEAST driver')  
             call wwrite_n(fout) 
          end if

          if (info==-1) then
             call wwrite_s(fout, 'ERROR with Internal memory allocation')  
             call wwrite_n(fout) 
          end if

          if (info==1) then
             call wwrite_s(fout, '==>WARNING: No eigenvalue has been found in the proposed search interval')
             call wwrite_n(fout)
          endif

          if (info==3) then
             call wwrite_s(fout, '==>WARNING: Size subspace M0 too small')  
             call wwrite_n(fout)
          end if

          if (info==4) then
             call wwrite_s(fout, '==>WARNING: Only the subspace has been returned')  
             call wwrite_n(fout)
          end if


          if (info==2) then
             call wwrite_s(fout, '==>WARNING: FEAST did not converge "yet" (#loop reaches maximum allowed)')  
             call wwrite_n(fout)
          end if

          if (info==0) then
             call wwrite_s(fout, '==>FEAST has successfully converged (to desired tolerance)')  
             call wwrite_n(fout) 
          else
             call wwrite_s(fout, '==>INFO code = ') 
             call wwrite_i(fout,info)
             call wwrite_n(fout)
          end if


          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_s(fout, '*********** FEAST- END*************************')
          call wwrite_n(fout) 
          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_n(fout)
       endif
    end if

  end subroutine dfeast_srci






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
  
  
  subroutine zfeast_hrci(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST RCI (Reverse Communication Interfaces) 
    !  Solve generalized Aq=lambda Bq eigenvalue problems
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE (or B Identity)  
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  ijob       (input/output) INTEGER :: ID of the RCI
    !                            INPUT on first entry: ijob=-1 
    !                            OUTPUT Return values (0,10,20,21,30,40)-- see FEAST documentation
    !  N          (input)        INTEGER: Size system
    !  work       (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
    !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
    !  zAq,zSq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm     (input/output) INTEGER(*) : FEAST parameters (see FEAST documentation)
    !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                          
    !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================

    implicit none
!    !-------------------------------------
!#ifdef MPI
!    include 'mpif.h'
!#endif
!    !-------------------------------------
    include "f90_noruntime_interface.fi"
    integer :: ijob,N,M0
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
    complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zSq
    integer,dimension(*) :: fpm
    double precision :: epsout 
    integer :: loop
    double precision :: Emin,Emax
    double precision,dimension(*)  :: lambda
    complex(kind=(kind(1.0d0))),dimension(N,*):: q
    integer :: mode
    double precision,dimension(*) :: res
    integer :: info
    !! parameters
    double precision, Parameter :: pi=3.1415926535897932d0
    double precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO), ZEROC=(DZERO,DZERO)
    double precision, parameter :: ba=-pi/2.0d0,ab=pi/2.0d0
    integer(8),parameter :: fout =6
    !! variable for FEAST
    integer :: i,e,j,k
    integer,dimension(4) :: iseed
    double precision :: theta,r,Emid
    complex(kind=(kind(1.0d0))) :: jac,aux
    double precision ::xe,we ! Gauss-Legendre
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer :: zSqo
    integer, dimension(:),pointer :: fpm_default
    logical :: testconv
    double precision :: trace
    !! Lapack variable (reduced system)
    character(len=1) :: JOBZ,UPLO
    double precision, dimension(:),pointer :: work_loc,work_loc2
    complex(kind=(kind(1.0d0))), dimension(:),pointer :: zwork_loc
    integer :: lwork_loc,info_lap,infoloc
    !! MPI compatibility variables
    integer :: rank,code,nb_procs,NEW_COMM_WORLD

    rank=0
    nb_procs=1
!    !----------------------------------------------
!#ifdef MPI
!    NEW_COMM_WORLD=fpm(9)
!    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
!    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
!    if (rank/=0) fpm(1)=0 ! comment only in rank 0 if any
!#endif
!    !---------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialization!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (ijob==-1) then 
       info=0 ! default value
       if (fpm(1)==1) then
          call wwrite_n(fout)
          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_s(fout, '*********** FEAST- BEGIN **********************')
          call wwrite_n(fout) 
          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout)
          call wwrite_s(fout, 'Routine ZFEAST_H{}{}')
          call wwrite_n(fout)
!!!!!!!!!!!! Print the FEAST parameters which has been changed from default
          call wwrite_s(fout, 'List of input parameters fpm(1:64)-- if different from default')
          call wwrite_n(fout)

          call wallocate_1i(fpm_default,64,infoloc)
          call feastinit(fpm_default)
          do i=1,20
             if (fpm(i)/=fpm_default(i)) then
                call wwrite_s(fout, '   fpm(')
                call wwrite_i(fout, i)
                call wwrite_s(fout, ')=')
                call wwrite_i(fout, fpm(i))
                call wwrite_n(fout)
             endif
          enddo
          call wdeallocate_1i(fpm_default)

          call wwrite_s(fout, 'Search interval [')
          call wwrite_d(fout,Emin)
          call wwrite_s(fout, '; ')
          call wwrite_d(fout,Emax)
          call wwrite_s(fout, ']')
          call wwrite_n(fout)
       end if
       call check_fpm_input(fpm,info)
       call dcheck_rci_input(Emin,Emax,M0,N,info)

       if (info/=0) fpm(21)=100 ! The End

!!!!!!!!!!!!!!!!
       IF (info==0) then
          if (fpm(1)==1) then
             call wwrite_s(fout, 'Size subspace')  
             call wwrite_t(fout) 
             call wwrite_i(fout,M0)
             call wwrite_n(fout)
             call wwrite_s(fout, '#Loop | #Eig  |       Trace           |     Error-Trace       |     Max-Residual')
             call wwrite_n(fout)
          endif
          fpm(23)=min(M0,N) ! 'current M0' size (global value)
          fpm(25)=fpm(23) !! 'current M0' size (by default)
          fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
!          !----------------------------------------------
!#ifdef MPI
!          if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
!             fpm(25)=fpm(23)/nb_procs ! local size of 'current M0'
!             if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
!             fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
!          end if
!#endif
!          !-----------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          fpm(22)=fpm(2) ! only one half contour necessary here
          loop=0
          fpm(21)=1 ! prepare reentry
          if (fpm(5)==0) then !!! random vectors
             iseed=(/56,890,3456,2333/)
             call ZLARNV(3,iseed,N*fpm(23),work(1,1))
          elseif (fpm(5)==1) then !!!!!! q is the initial guess
!             !----------------------------------------
!#ifdef MPI
!             work(1:N,1:fpm(23))=ZEROC 
!#endif
!             !------------------------------------------
             ijob=40 !! B*q=>work
             return
          end if
       end IF ! info=0
!!!!!!!!!!!!!!
    end if   !ijob=-1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    IF (fpm(21)==1) then !! we initialize a new contour integration
!       !------------------------------------------------------------------------
!#ifdef MPI
!       if ((loop>0).or.(fpm(5)==1)) then        
!          if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work,N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
!       end if
!#endif
!       !------------------------------------------------------------------------
       q(1:N,1:fpm(23))=ZEROC
       fpm(20)=1
       fpm(21)=2
       ijob=-2 ! just initialization 
    end IF


!!!!!!!!!!!!
    IF (fpm(21)==2) then !! we start or pursue the contour integration

       IF (info==0) then !! will end up checking info errors returned by FEAST drivers
          do e=fpm(20)+rank,fpm(22),nb_procs !!!! loop over the contour 

             if (ijob==-2) then !!Factorize the linear system (complex) (zS-A)
                call dset_point_gauss_legendre(fpm(2),e,xe,we) !! Gauss-points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0d0
                Emid=Emin+r
                Ze=Emid*ONEC+r*ONEC*wdcos(theta)+r*(DZERO,DONE)*wdsin(theta)                               
                fpm(20)=e-rank
                ijob=10 ! for fact
                if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor 
             endif

             if (ijob==10) then !!Solve the linear system (complex) (zS-A)q=v 
                call ZLACPY( 'F', N, fpm(23),work , N, workc, N )
                ijob=11 ! for solve
                return
             endif

             if (ijob==11) then 
                !! summation 
                call dset_point_gauss_legendre(fpm(2),e,xe,we) !! Gauss-points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0d0 
                jac=(r*(DZERO,DONE)*wdsin(theta)+ONEC*r*wdcos(theta))
                aux=-(ONEC/(4.0d0))*we*jac
                call ZAXPY(N*fpm(23),aux,workc,1,q,1)
                !!Explicit Factorization of the linear system (complex) (zS-A)^T 
                !!needed if driver not capable to exploit Factorization of (zS-A) for solving (zS-A)^Tq=v          
                ijob=20 ! for fact
                if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor
             endif

             if (ijob==20) then!!!! Solve the linear system (complex) (zS-A)^Tq=v  
                call ZLACPY( 'F', N, fpm(23),work , N, workc, N )
                ijob=21 ! for solve with transpose
                return
             end if

             if (ijob==21) then
                !! summation
                call dset_point_gauss_legendre(fpm(2),e,xe,we) !! Gauss-points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0d0 
                jac=(r*(DZERO,DONE)*wdsin(theta)+ONEC*r*wdcos(theta))
                aux=-(ONEC/(4.0d0))*we*conjg(jac)
                call ZAXPY(N*fpm(23),aux,workc,1,q,1)
                ijob=-2 ! just for identification
             end if

          end do
       end IF !! info=0

!       !------------------------------------------------
!#ifdef MPI
!       call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!#endif
!       !-----------------------------------------------  
       if (info/=0) fpm(21)=100 ! the end

       if (info==0) then
          fpm(21)=4 
!          !------------------------------------------------
!#ifdef MPI
!          call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
!#endif
!          !-----------------------------------------------  
       end if
    end IF   ! fpm(21)==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ((fpm(21)==4).and.(fpm(14)==1)) then !! only q vectors has been computed and is returned
info=4
       if (info/=0) fpm(21)=100 ! The End
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq=eq Sq xq     with Aq=Q^TAQ Sq=Q^TAQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!form Aq=> Aq=Q^T A Q 
    if (fpm(21)==4) then
       fpm(21)=5 ! preparing reentry
       ijob=30 
       return  ! mat-vec A*q => work
    endif

    if (fpm(21)==5) then 
!       !------------------------------------------------
!#ifdef MPI 
!       zAq(1:M0,1:fpm(23))=ZEROC
!#endif
!       !-------------------------------------------------
       call ZGEMM('C','N',fpm(23),fpm(25),N,ONEC,q(1,1),N,work(1,fpm(24)),N,ZEROC,zAq(1,fpm(24)),M0)
       fpm(21)=6
    endif

!!!!!!!!!form  Sq=> Sq=Q^T S Q
    if (fpm(21)==6) then 
       fpm(21)=7 ! preparing reenty
       ijob=40 
       return! mat-vec S*q => work
    end if
    if (fpm(21)==7) then
!       !------------------------------------------------
!#ifdef MPI
!       zSq(1:M0,1:fpm(23))=ZEROC
!#endif
!       !-------------------------------------------------
       call ZGEMM('C','N',fpm(23),fpm(25),N,ONEC,q(1,1),N,work(1,fpm(24)),N,ZEROC,zSq(1,fpm(24)),M0)
       fpm(21)=8
    endif

    if (fpm(21)==8) then
!       !----------------------------------------!(zAq,zSq known to all processors) 
!#ifdef MPI 
!       if (fpm(23)/nb_procs>=1) then
!          call MPI_ALLREDUCE(MPI_IN_PLACE,zAq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
!          call MPI_ALLREDUCE(MPI_IN_PLACE,zSq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
!       end if
!#endif
!       !---------------------------------------
       if (fpm(12)==1) then ! customize eigenvalue solver
          fpm(21)=9 ! preparing reentry - could return new value of M0 in fpm(23) if reduced subspace is needed
          ijob=50
          return
       endif
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==8) then
       if (rank==0) then
          JOBZ='V'
          UPLO='L'
          info_lap=1 ! initialization
          i=1
          LWORK_LOC=2*fpm(23)-1 !! for lapack eig reduced system
          call wallocate_1z(zWORK_LOC,LWORK_LOC,infoloc)
          call wallocate_1d(WORK_LOC,3*fpm(23)-2,infoloc)
          if (infoloc/=0) info=-1
          do while ((info_lap/=0).and.(info==0))
             i=i+1
             if (i==10) info=-3 ! arbitrary maximum
             call wallocate_2z(zSqo,fpm(23),fpm(23),infoloc)
             if (infoloc/=0) info=-1

             ! # eigenvalues inside search (wait second contour)
             mode=0
             if (loop>0) then
                call ZLACPY( 'F', fpm(23), fpm(23),zSq , M0, zSqo, fpm(23) )
                if (infoloc/=0) info=-1
                call ZHEEV( 'N', UPLO, fpm(23), zSqo, fpm(23), lambda, zwork_loc,lwork_loc,work_loc, INFO_lap )
                do j=1,fpm(23)
                   if (sqrt(abs(lambda(j)))>=0.5d0) mode=mode+1
                enddo
             end if

             call ZLACPY( 'F', fpm(23), fpm(23),zSq , M0, zSqo, fpm(23) )
             call ZHEGV(1, JOBZ, UPLO, fpm(23), zAq, M0, zSqo, fpm(23), lambda, zWORK_loc,Lwork_loc, WORK_loc, INFO_lap)

             if ((info_lap<=fpm(23)).and.(info_lap/=0)) info=-3
             if (info_lap>fpm(23)) then !! zSqo is not spd (a posteriori resize subspace)
                fpm(23)=info_lap-fpm(23)-1
                if (fpm(1)==1) then
                   call wwrite_s(fout, 'Resize subspace')  
                   call wwrite_t(fout) 
                   call wwrite_i(fout,fpm(23))
                   call wwrite_n(fout)
                end if
             end if
             call wdeallocate_2z(zSqo)
          end do
          call wdeallocate_1z(zwork_loc)
          call wdeallocate_1d(work_loc)
       end if !(rank 0)
!       !-------------------------------- !(info common to all processors) 
!#ifdef MPI
!       call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!#endif
!       !--------------------------------
       if (info/=0) fpm(21)=100 ! the end

       if (info==0) then
          fpm(25)=fpm(23) !! current M0 size (by default) -- global
          fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
!          !---------------------------------------- !(zAq==> vectors, lambda and fpm(23), known to all processors) 
!#ifdef MPI
!          call MPI_BCAST(mode,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!          call MPI_BCAST(fpm(23),1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!          call MPI_BCAST(zAq,M0*fpm(23),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
!          call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_PRECISION,0,NEW_COMM_WORLD,code)
!          if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
!             fpm(25)=fpm(23)/nb_procs ! local size of current M0
!             if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
!             fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
!          end if
!#endif
!          !-----------------------------
          fpm(21)=9
       end if
    end if !! fpm(21)=8




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==9) then
       call ZLACPY( 'F', N, fpm(23),q , N, work, N ) 
       !! option - non shifted
!#ifdef MPI 
!          q(1:N,1:fpm(23))=ZEROC
!#endif
          call ZGEMM('N','N',N,fpm(25),fpm(23),ONEC,work(1,1),N,zAq(1,fpm(24)),M0,ZEROC,q(1,fpm(24)),N)
!#ifdef MPI       
!          if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
!#endif

    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residual ||AX-lambda*BX||
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==9) then
       fpm(21)=10 ! preparing reentry
       ijob=30 
       return  ! mat-vec A*q => work
    endif

    if (fpm(21)==10) then
       call ZLACPY( 'F', N, fpm(25),work(1,fpm(24)), N, workc(1,fpm(24)), N )
       fpm(21)=11 ! preparing reentry
       ijob=40 
!       !----------------------------------------
!#ifdef MPI
!       work(1:N,1:fpm(23))=ZEROC  !! work is also needed for outer-loop if any 
!#endif
!       !------------------------------------------
       return  ! mat-vec S*q => work 
    endif

    if (fpm(21)==11) then
!       !----------------------------------------
!#ifdef MPI
!       res(1:fpm(23))=DZERO
!#endif
!       !-------- Absolute residual
       do i=fpm(24),fpm(24)+fpm(25)-1
          res(i)=sum(abs(workc(1:N,i)-lambda(i)*work(1:N,i)))/sum(abs(max(abs(Emin),abs(Emax))*(work(1:N,i))))
       end do
!       !----------------------------------------
!#ifdef MPI 
!       if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
!#endif
!       !-----------------------------------------
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Trace / count eigenvalues (remove spurious if any)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (fpm(21)==11) then
       k=0
       trace=DZERO
       theta=DZERO
       do i=1,fpm(23)
          if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
             k=k+1
             trace=trace+lambda(i)
             if (res(i)>theta) theta=res(i) ! max residual
          endif
       enddo

!!!!!!!!! remove spurious if any       
       if (loop==0) mode=k ! wait second contour

       if (mode<k) then
          do j=1,k-mode
             theta=DZERO
             e=1
             do i=1,fpm(23)
                if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
                   if (res(i)>theta) then
                      e=i
                      theta=res(i) !max
                   endif
                end if
             enddo
             trace=trace-lambda(e) ! update trace
             res(e)=-DONE !spurious
          end do
!!! max residual
          theta=DZERO
          do i=1,fpm(23)
             if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
                if (res(i)>theta) theta=res(i) 
             end if
          enddo
       end if

       if (mode==0) info=1  ! no mode detected in the interval
if (loop>0) then
       if ((mode==M0).and.(mode/=N)) info=3 !size subspace too small
endif
       if (info/=0) fpm(21)=100 ! The End
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check FEAST Convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

    if (fpm(21)==11) then
       testconv=.false. ! initialization

       !! trace convergence
       if (loop==0) epsout=DONE
       if (loop>0) then
          epsout=(abs(trace-epsout))/max(abs(Emin),abs(Emax))           
          if ((fpm(6)==0).and.(log10(epsout)<(-fpm(3)))) testconv=.true.
       end if

       !! residual convergence
       if ((fpm(6)==1).and.(log10(theta)<(-fpm(3)))) testconv=.true.

       !! Two loops minimum if spurious are found
       if ((loop<=1).and.(k>mode)) testconv=.false.

       if (fpm(1)==1) then
          call wwrite_i(fout,loop)
          call wwrite_t(fout) 
          call wwrite_i(fout,mode)
          call wwrite_t(fout)
          call wwrite_d(fout,trace)
          call wwrite_t(fout) 
          call wwrite_d(fout,epsout)
          call wwrite_t(fout) 
          call wwrite_d(fout,theta)
          call wwrite_n(fout) 
       end if

       if (.not.testconv) then
          epsout=trace
          if (loop==fpm(4)) then
             info=2 ! FEAST did not converge (#loop reaches maximum)
             testconv=.true. ! return final eigenvector anyway
          endif
       endif
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! FEAST exit IF Convergence - FEAST iteration IF NOT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

    if (fpm(21)==11) then

       if (.not.testconv) then !!! need FEAST refinement loop

!!!!!!!!!!! here work=S*q is known alreasy,  using FEAST-MPI work is already partionned 
          fpm(21)=1   ! prepare reentry- reloop (with contour integration)
          !fpm(21)=4 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
          loop=loop+1
          ijob=-2 ! do nothing
          return  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       else    !!!!!!! final eigenvectors/eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! reorder (shift) lambda, eigenvector and residual

   call ZLACPY( 'F', N, fpm(23),q , N, work, N )
          call wallocate_1d(work_loc,fpm(23),infoloc)
          call wallocate_1d(work_loc2,fpm(23),infoloc)
          call DCOPY(fpm(23),lambda , 1, work_loc, 1 )
          call DCOPY(fpm(23),res , 1, work_loc2, 1 )
          q(1:N,1:fpm(23))=ZEROC
          e=0
          j=0
          k=0
          do i=1,fpm(23)
             if ((work_loc(i)>Emin).and.(work_loc(i)<Emax)) then ! inside the search interval
                if (work_loc2(i)/=-DONE) then ! spurious 
                   e=e+1
                   q(1:N,e)=work(1:N,i)
                   lambda(e)=work_loc(i)
                   res(e)=work_loc2(i)
                endif
             else
                j=j+1
                q(1:N,mode+j)=work(1:N,i)
                lambda(mode+j)=work_loc(i)
                res(mode+j)=work_loc2(i)
             end if
             if (work_loc2(i)==-DONE) then ! spurious at the end
                k=k+1
                q(1:N,fpm(23)-k+1)=work(1:N,i)
                lambda(fpm(23)-k+1)=work_loc(i)
                res(fpm(23)-k+1)=work_loc2(i)
             endif
          enddo
          call wdeallocate_1d(work_loc)
          call wdeallocate_1d(work_loc2)

!!!!!!!!!!!!!!!!!!!!!!!!!!
          M0=fpm(23)  ! update value of M0 (new subspace)
          fpm(21)=100 ! The End

       end if ! test convergence
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==100) then !! THE END (ijob=0) 
       ijob=0 !! exit FEAST

       if (fpm(1)==1) then !! Print  Information

          if (info>=200) then
             call wwrite_s(fout, 'PROBLEM with input parameters')
             call wwrite_n(fout) 
          end if

          if ((info>100).and.(info<200)) then
             call wwrite_s(fout, 'PROBLEM with FEAST array parameters')
             call wwrite_n(fout) 
          end if

          if (info==-3) then
             call wwrite_s(fout, 'ERROR with reduced system')  
             call wwrite_n(fout) 
          end if

          if (info==-2) then
             call wwrite_s(fout, 'ERROR from Inner Linear System Solver in FEAST driver')  
             call wwrite_n(fout) 
          end if

          if (info==-1) then
             call wwrite_s(fout, 'ERROR with Internal memory allocation')  
             call wwrite_n(fout) 
          end if

          if (info==1) then
             call wwrite_s(fout, '==>WARNING: No eigenvalue has been found in the proposed search interval')
             call wwrite_n(fout)
          endif

          if (info==3) then
             call wwrite_s(fout, '==>WARNING: Size subspace M0 too small')  
             call wwrite_n(fout)
          end if

          if (info==4) then
             call wwrite_s(fout, '==>WARNING: Only the subspace has been returned')  
             call wwrite_n(fout)
          end if



          if (info==2) then
             call wwrite_s(fout, '==>WARNING: FEAST did not converge "yet" (#loop reaches maximum allowed)')  
             call wwrite_n(fout)
          end if

          if (info==0) then
             call wwrite_s(fout, '==>FEAST has successfully converged (to desired tolerance)')  
             call wwrite_n(fout) 
          else
             call wwrite_s(fout, '==>INFO code = ') 
             call wwrite_i(fout,info)
             call wwrite_n(fout)
          end if


          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_s(fout, '*********** FEAST- END*************************')
          call wwrite_n(fout) 
          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_n(fout)
       endif
    end if

  end subroutine zfeast_hrci






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
    
  subroutine sfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST RCI (Reverse Communication Interfaces) 
    !  Solve generalized Aq=lambda Bq eigenvalue problems
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE (or B Identity) 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  ijob       (input/output) INTEGER :: ID of the RCI
    !                            INPUT on first entry: ijob=-1 
    !                            OUTPUT Return values (0,10,20,21,30,40)-- see FEAST documentation
    !  N          (input)        INTEGER: Size system
    !  work       (input/output) REAL SINGLE PRECISION (N,M0)   :  Workspace 
    !  workc      (input/output) COMPLEX SINGLE PRECISION (N,M0):  Workspace 
    !  Aq,Sq      (input/output) REAL SINGLE PRECISION (M0,M0)  :  Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm     (input/output) INTEGER(*) : FEAST parameters (see FEAST documentation)
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  lambda     (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) REAL SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                 
    !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================

    implicit none
!    !-------------------------------------
!#ifdef MPI
!    include 'mpif.h'
!#endif
!    !-------------------------------------
    include "f90_noruntime_interface.fi"
    integer :: ijob,N,M0
    complex :: Ze
    real, dimension(N,*) ::work
    complex,dimension(N,*):: workc
    real,dimension(M0,*):: Aq,Sq
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    real,dimension(*)  :: lambda
    real,dimension(N,*):: q
    integer :: mode
    real,dimension(*) :: res
    integer :: info
    !! parameters
    real, Parameter :: pi=3.1415926535897932e0
    real, Parameter :: SONE=1.0E0, SZERO=0.0E0
    complex,parameter :: ONEC=(SONE,SZERO), ZEROC=(SZERO,SZERO)
    real, parameter :: ba=-pi/2.0e0, ab=pi/2.0e0
    integer(8),parameter :: fout =6
    !! variable for FEAST
    integer :: i,e,j,k
    integer,dimension(4) :: iseed
    real :: theta,r,Emid
    complex :: jac
    real ::xe,we ! Gauss-Legendre
    real, dimension(:,:),pointer :: Sqo
    integer, dimension(:),pointer :: fpm_default
    logical :: testconv
    real :: trace
    !! Lapack variable (reduced system)
    character(len=1) :: JOBZ,UPLO
    real, dimension(:),pointer :: work_loc,work_loc2
    integer :: lwork_loc,info_lap,infoloc
    !! MPI compatibility variables
    integer :: rank,code,nb_procs,NEW_COMM_WORLD

    rank=0
    nb_procs=1
!    !----------------------------------------------
!#ifdef MPI
!    NEW_COMM_WORLD=fpm(9)
!    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
!    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
!    if (rank/=0) fpm(1)=0 ! comment only in rank 0 if any
!#endif
!    !---------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialization!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (ijob==-1) then 
       info=0 ! default value
       if (fpm(1)==1) then
          call wwrite_n(fout)
          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_s(fout, '*********** FEAST- BEGIN **********************')
          call wwrite_n(fout) 
          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout)
          call wwrite_s(fout, 'Routine SFEAST_S{}{}')
          call wwrite_n(fout)
!!!!!!!!!!!! Print the FEAST parameters which has been changed from default
          call wwrite_s(fout, 'List of input parameters fpm(1:64)-- if different from default')
          call wwrite_n(fout)

          call wallocate_1i(fpm_default,64,infoloc)
          call feastinit(fpm_default)
          do i=1,20
             if (fpm(i)/=fpm_default(i)) then
                call wwrite_s(fout, '   fpm(')
                call wwrite_i(fout, i)
                call wwrite_s(fout, ')=')
                call wwrite_i(fout, fpm(i))
                call wwrite_n(fout)
             endif
          enddo
          call wdeallocate_1i(fpm_default)

          call wwrite_s(fout, 'Search interval [')
          call wwrite_f(fout,Emin)
          call wwrite_s(fout, '; ')
          call wwrite_f(fout,Emax)
          call wwrite_s(fout, ']')
          call wwrite_n(fout)
       end if
       call check_fpm_input(fpm,info)
       call scheck_rci_input(Emin,Emax,M0,N,info)

       if (info/=0) fpm(21)=100 ! The End

!!!!!!!!!!!!!!!!
       IF (info==0) then
          if (fpm(1)==1) then
             call wwrite_s(fout, 'Size subspace')  
             call wwrite_t(fout) 
             call wwrite_i(fout,M0)
             call wwrite_n(fout) 
             call wwrite_s(fout, '#Loop | #Eig  |    Trace     |   Error-Trace |  Max-Residual') 
             call wwrite_n(fout)
          endif
          fpm(23)=min(M0,N) ! 'current M0' size (global value)
          fpm(25)=fpm(23) !! 'current M0' size (by default)
          fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
!          !----------------------------------------------
!#ifdef MPI
!          if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
!             fpm(25)=fpm(23)/nb_procs ! local size of 'current M0'
!             if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
!             fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
!          end if
!#endif
!          !-----------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          fpm(22)=fpm(2) ! only one half contour necessary here
          loop=0
          fpm(21)=1 ! prepare reentry
          if (fpm(5)==0) then !!! random vectors
             iseed=(/56,890,3456,2333/)
             call SLARNV(3,iseed,N*fpm(23),work(1,1))
          elseif (fpm(5)==1) then !!!!!! q is the initial guess
!             !----------------------------------------
!#ifdef MPI
!             work(1:N,1:fpm(23))=SZERO 
!#endif
!             !------------------------------------------
             ijob=40 !! B*q=>work
             return
          end if
       end IF ! info=0
!!!!!!!!!!!!!!
    end if   !ijob=-1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    IF (fpm(21)==1) then !! we initialize a new contour integration
!       !------------------------------------------------------------------------
!#ifdef MPI
!       if ((loop>0).or.(fpm(5)==1)) then        
!          if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work,N*fpm(23),MPI_REAL,MPI_SUM,NEW_COMM_WORLD,code)
!       end if
!#endif
!       !------------------------------------------------------------------------
       q(1:N,1:fpm(23))=SZERO
       fpm(20)=1
       fpm(21)=2
       ijob=-2 ! just initialization 
    end IF


!!!!!!!!!!!!
    IF (fpm(21)==2) then !! we start or pursue the contour integration

       IF (info==0) then !! will end up checking info errors returned by FEAST drivers
          do e=fpm(20)+rank,fpm(22),nb_procs !!!! loop over the contour 

             if (ijob==-2) then !!Factorize the linear system (complex) (zS-A)
                call sset_point_gauss_legendre(fpm(2),e,xe,we) !! Gauss-points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0e0
                Emid=Emin+r
                Ze=Emid*ONEC+r*ONEC*wscos(theta)+r*(SZERO,SONE)*wssin(theta)                               
                fpm(20)=e-rank
                ijob=10 ! for fact
                if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor 
             endif

             if (ijob==10) then !!Solve the linear system (complex) (zS-A)q=v 
                call CLACP2( 'F', N, fpm(23),work , N, workc, N )
                ijob=11 ! for solve
                return
             endif

             if (ijob==11) then 
                !! summation 
                call sset_point_gauss_legendre(fpm(2),e,xe,we) !! Gauss-points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0e0 
                jac=(r*(SZERO,SONE)*wssin(theta)+ONEC*r*wscos(theta))
                q(1:N,1:fpm(23))=q(1:N,1:fpm(23))-(SONE/2.0e0)*we*real(jac*workc(1:N,1:fpm(23)))                  
                ijob=-2 ! just for identification
             end if
          end do
       end IF  !! info=0

!       !------------------------------------------------
!#ifdef MPI
!       call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!#endif
!       !-----------------------------------------------  
       if (info/=0) fpm(21)=100 ! the end

       if (info==0) then
          fpm(21)=4 
!          !------------------------------------------------
!#ifdef MPI
!          call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_REAL,MPI_SUM,NEW_COMM_WORLD,code)
!#endif
!          !-----------------------------------------------  
       end if
    end IF    ! fpm(21)==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ((fpm(21)==4).and.(fpm(14)==1)) then !! only q vectors has been computed and is returned
info=4
       if (info/=0) fpm(21)=100 ! The End
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq=eq Sq xq     with Aq=Q^TAQ Sq=Q^TAQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!form Aq=> Aq=Q^T A Q 
    if (fpm(21)==4) then
       fpm(21)=5 ! preparing reentry
       ijob=30 
       return  ! mat-vec A*q => work
    endif

    if (fpm(21)==5) then 
!       !------------------------------------------------
!#ifdef MPI 
!       Aq(1:M0,1:fpm(23))=SZERO
!#endif
!       !-------------------------------------------------
       call SGEMM('T','N',fpm(23),fpm(25),N,SONE,q(1,1),N,work(1,fpm(24)),N,SZERO,Aq(1,fpm(24)),M0)
       fpm(21)=6
    endif

!!!!!!!!!form  Sq=> Sq=Q^T S Q
    if (fpm(21)==6) then 
       fpm(21)=7 ! preparing reenty
       ijob=40 
       return! mat-vec S*q => work
    end if
    if (fpm(21)==7) then
!       !------------------------------------------------
!#ifdef MPI
!       Sq(1:M0,1:fpm(23))=SZERO
!#endif
!       !-------------------------------------------------
       call SGEMM('T','N',fpm(23),fpm(25),N,SONE,q(1,1),N,work(1,fpm(24)),N,SZERO,Sq(1,fpm(24)),M0)
       fpm(21)=8
    endif

    if (fpm(21)==8) then
!       !---------------------------------------- !(Aq,Sq known to all processors) 
!#ifdef MPI
!       if (fpm(23)/nb_procs>=1) then
!          call MPI_ALLREDUCE(MPI_IN_PLACE,Aq(1,1),M0*fpm(23),MPI_REAL,MPI_SUM,NEW_COMM_WORLD,code)
!          call MPI_ALLREDUCE(MPI_IN_PLACE,Sq(1,1),M0*fpm(23),MPI_REAL,MPI_SUM,NEW_COMM_WORLD,code)
!       end if
!#endif
!       !---------------------------------------
       if (fpm(12)==1) then ! customize eigenvalue solver
          fpm(21)=9 ! preparing reentry - could return new value of M0 in fpm(23) if reduced subspace is needed
          ijob=50
          return
       endif
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==8) then

       if (rank==0) then
          JOBZ='V'
          UPLO='L'
          info_lap=1 ! initialization
          i=1
          LWORK_LOC=3*fpm(23)-1 !! for lapack eig reduced system
          call wallocate_1s(WORK_LOC,LWORK_LOC,infoloc)
          if (infoloc/=0) info=-1
          do while ((info_lap/=0).and.(info==0))
             i=i+1
             if (i==10) info=-3 ! arbitrary maximum
             call wallocate_2s(Sqo,fpm(23),fpm(23),infoloc)
             if (infoloc/=0) info=-1

             ! # eigenvalues inside search (wait second contour)
             mode=0
             if (loop>0) then
                call SLACPY( 'F', fpm(23), fpm(23),Sq , M0, Sqo, fpm(23) )
                if (infoloc/=0) info=-1
                call SSYEV( 'N', UPLO, fpm(23), Sqo, fpm(23), lambda, work_loc, lwork_loc, INFO_lap )
                do j=1,fpm(23)
                   if (sqrt(abs(lambda(j)))>=0.5e0) mode=mode+1
                enddo
             end if

             call SLACPY( 'F', fpm(23), fpm(23),Sq , M0, Sqo, fpm(23) )
             call SSYGV(1, JOBZ, UPLO, fpm(23),  Aq, M0,  Sqo, fpm(23), lambda,work_loc,Lwork_loc,INFO_lap)

             if ((info_lap<=fpm(23)).and.(info_lap/=0)) info=-3
             if (info_lap>fpm(23)) then !! zSqo is not spd (a posteriori resize subspace)
                fpm(23)=info_lap-fpm(23)-1
                if (fpm(1)==1) then
                   call wwrite_s(fout, 'Resize subspace')  
                   call wwrite_t(fout) 
                   call wwrite_i(fout,fpm(23))
                   call wwrite_n(fout)
                end if
             end if
             call wdeallocate_2s(Sqo)
          end do
          call wdeallocate_1s(work_loc)
       end if !(rank 0)
!       !-------------------------------- !(info common to all processors)
!#ifdef MPI
!       call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!#endif
!       !--------------------------------
       if (info/=0) fpm(21)=100 ! the end

       if (info==0) then
          fpm(25)=fpm(23) !! current M0 size (by default) -- global
          fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
!          !----------------------------------------!(Aq==> vectors, lambda and fpm(23), known to all processors) 
!#ifdef MPI 
!          call MPI_BCAST(mode,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!          call MPI_BCAST(fpm(23),1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!          call MPI_BCAST(Aq,M0*fpm(23),MPI_REAL,0,NEW_COMM_WORLD,code)
!          call MPI_BCAST(lambda,fpm(23),MPI_REAL,0,NEW_COMM_WORLD,code)
!          if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
!             fpm(25)=fpm(23)/nb_procs ! local size of current M0
!             if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
!             fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
!          end if
!#endif
!          !-----------------------------
          fpm(21)=9
       end if
    end if !! fpm(21)=8



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==9) then
       call SLACPY( 'F', N, fpm(23),q , N, work, N ) 
       !! option - non shifted
!#ifdef MPI 
!          q(1:N,1:fpm(23))=SZERO
!#endif
          call SGEMM('N','N',N,fpm(25),fpm(23),SONE,work(1,1),N,Aq(1,fpm(24)),M0,SZERO,q(1,fpm(24)),N)
!#ifdef MPI       
!          if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_REAL,MPI_SUM,NEW_COMM_WORLD,code)
!#endif

    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residual ||AX-lambda*BX||
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==9) then
       fpm(21)=10 ! preparing reentry
       ijob=30 
       return  ! mat-vec A*q => work
    endif

    if (fpm(21)==10) then
       call CLACP2( 'F', N, fpm(25),work(1,fpm(24)), N, workc(1,fpm(24)), N )
       fpm(21)=11 ! preparing reentry
       ijob=40 
!       !----------------------------------------
!#ifdef MPI
!       work(1:N,1:fpm(23))=SZERO  !! work is also needed for outer-loop if any 
!#endif
!       !------------------------------------------
       return  ! mat-vec S*q => work 
    endif

    if (fpm(21)==11) then
!       !----------------------------------------
!#ifdef MPI
!       res(1:fpm(23))=SZERO
!#endif
!       !-------- Absolute residual
       do i=fpm(24),fpm(24)+fpm(25)-1
          res(i)=sum(abs(real(workc(1:N,i))-lambda(i)*work(1:N,i)))/sum(abs(max(abs(Emin),abs(Emax))*(work(1:N,i))))
       end do
!       !----------------------------------------
!#ifdef MPI 
!       if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_REAL,MPI_SUM,NEW_COMM_WORLD,code)
!#endif
!       !-----------------------------------------
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Trace / count eigenvalues (remove spurious if any)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (fpm(21)==11) then
       k=0
       trace=SZERO
       theta=SZERO
       do i=1,fpm(23)
          if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
             k=k+1
             trace=trace+lambda(i)
             if (res(i)>theta) theta=res(i) ! max residual
          endif
       enddo

!!!!!!!!! remove spurious if any       
       if (loop==0) mode=k ! wait second contour

       if (mode<k) then
          do j=1,k-mode
             theta=SZERO
             e=1
             do i=1,fpm(23)
                if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
                   if (res(i)>theta) then
                      e=i
                      theta=res(i) !max
                   endif
                end if
             enddo
             trace=trace-lambda(e) ! update trace
             res(e)=-SONE !spurious
          end do
!!! max residual
          theta=SZERO
          do i=1,fpm(23)
             if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
                if (res(i)>theta) theta=res(i) 
             end if
          enddo
       end if

       if (mode==0) info=1  ! no mode detected in the interval
if (loop>0) then
       if ((mode==M0).and.(mode/=N)) info=3 !size subspace too small
endif
       if (info/=0) fpm(21)=100 ! The End
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check FEAST Convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

    if (fpm(21)==11) then
       testconv=.false. ! initialization

       !! trace convergence
       if (loop==0) epsout=SONE
       if (loop>0) then
          epsout=(abs(trace-epsout))/max(abs(Emin),abs(Emax))         
          if ((fpm(6)==0).and.(log10(epsout)<(-fpm(7)))) testconv=.true.
       end if

       !! residual convergence
       if ((fpm(6)==1).and.(log10(theta)<(-fpm(7)))) testconv=.true.

       !! Two loops minimum if spurious are found
       if ((loop<=1).and.(k>mode)) testconv=.false.

       if (fpm(1)==1) then
          call wwrite_i(fout,loop)
          call wwrite_t(fout) 
          call wwrite_i(fout,mode)
          call wwrite_t(fout)
          call wwrite_f(fout,trace)
          call wwrite_t(fout) 
          call wwrite_f(fout,epsout)
          call wwrite_t(fout) 
          call wwrite_f(fout,theta)
          call wwrite_n(fout) 
       end if

       if (.not.testconv) then
          epsout=trace
          if (loop==fpm(4)) then
             info=2 ! FEAST did not converge (#loop reaches maximum)
             testconv=.true. ! return final eigenvector anyway
          endif
       endif
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! FEAST exit IF Convergence - FEAST iteration IF NOT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

    if (fpm(21)==11) then

       if (.not.testconv) then !!! need FEAST refinement loop

!!!!!!!!!!! here work=S*q is known alreasy,  using FEAST-MPI work is already partionned 
          fpm(21)=1   ! prepare reentry- reloop (with contour integration)
          !fpm(21)=4 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
          loop=loop+1
          ijob=-2 ! do nothing
          return  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       else    !!!!!!! final eigenvectors/eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! reorder (shift) lambda, eigenvector and residual

          call SLACPY( 'F', N, fpm(23),q , N, work, N )
          call wallocate_1s(work_loc,fpm(23),infoloc)
          call wallocate_1s(work_loc2,fpm(23),infoloc)
          call SCOPY(fpm(23),lambda , 1, work_loc, 1 )
          call SCOPY(fpm(23),res , 1, work_loc2, 1 )
          q(1:N,1:fpm(23))=SZERO
          e=0
          j=0
          k=0
          do i=1,fpm(23)
             if ((work_loc(i)>Emin).and.(work_loc(i)<Emax)) then ! inside the search interval
                if (work_loc2(i)/=-SONE) then ! spurious 
                   e=e+1
                   q(1:N,e)=work(1:N,i)
                   lambda(e)=work_loc(i)
                   res(e)=work_loc2(i)
                endif
             else
                j=j+1
                q(1:N,mode+j)=work(1:N,i)
                lambda(mode+j)=work_loc(i)
                res(mode+j)=work_loc2(i)
             end if
             if (work_loc2(i)==-SONE) then ! spurious at the end
                k=k+1
                q(1:N,fpm(23)-k+1)=work(1:N,i)
                lambda(fpm(23)-k+1)=work_loc(i)
                res(fpm(23)-k+1)=work_loc2(i)
             endif
          enddo
          call wdeallocate_1s(work_loc)
          call wdeallocate_1s(work_loc2)

!!!!!!!!!!!!!!!!!!!!!!!!!!
          M0=fpm(23)  ! update value of M0 (new subspace)
          fpm(21)=100 ! The End

       end if ! test convergence
    end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==100) then !! THE END (ijob=0) 
       ijob=0 !! exit FEAST

       if (fpm(1)==1) then !! Print  Information

          if (info>=200) then
             call wwrite_s(fout, 'PROBLEM with input parameters')
             call wwrite_n(fout) 
          end if

          if ((info>100).and.(info<200)) then
             call wwrite_s(fout, 'PROBLEM with FEAST array parameters')
             call wwrite_n(fout) 
          end if

          if (info==-3) then
             call wwrite_s(fout, 'ERROR with reduced system')  
             call wwrite_n(fout) 
          end if

          if (info==-2) then
             call wwrite_s(fout, 'ERROR from Inner Linear System Solver in FEAST driver')  
             call wwrite_n(fout) 
          end if

          if (info==-1) then
             call wwrite_s(fout, 'ERROR with Internal memory allocation')  
             call wwrite_n(fout) 
          end if

          if (info==1) then
             call wwrite_s(fout, '==>WARNING: No eigenvalue has been found in the proposed search interval')
             call wwrite_n(fout)
          endif

          if (info==3) then
             call wwrite_s(fout, '==>WARNING: Size subspace M0 too small')  
             call wwrite_n(fout)
          end if


          if (info==4) then
             call wwrite_s(fout, '==>WARNING: Only the subspace has been returned')  
             call wwrite_n(fout)
          end if


          if (info==2) then
             call wwrite_s(fout, '==>WARNING: FEAST did not converge "yet" (#loop reaches maximum allowed)')  
             call wwrite_n(fout)
          end if

          if (info==0) then
             call wwrite_s(fout, '==>FEAST has successfully converged (to desired tolerance)')  
             call wwrite_n(fout) 
          else
             call wwrite_s(fout, '==>INFO code = ') 
             call wwrite_i(fout,info)
             call wwrite_n(fout)
          end if


          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_s(fout, '*********** FEAST- END*************************')
          call wwrite_n(fout) 
          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_n(fout)
       endif
    end if

  end subroutine sfeast_srci






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
  
  
  subroutine cfeast_hrci(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST RCI (Reverse Communication Interfaces) 
    !  Solve generalized Aq=lambda Bq eigenvalue problems
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE (or B Identity)  
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  ijob       (input/output) INTEGER :: ID of the RCI
    !                            INPUT on first entry: ijob=-1 
    !                            OUTPUT Return values (0,10,20,21,30,40)-- see FEAST documentation
    !  N          (input)        INTEGER: Size system
    !  work       (input/output) COMPLEX SINGLE PRECISION (N,M0):   Workspace 
    !  workc      (input/output) COMPLEX SINGLE PRECISION (N,M0):   Workspace 
    !  zAq,zSq    (input/output) COMPLEX SINGLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm     (input/output) INTEGER(*) : FEAST parameters (see FEAST documentation)
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  lambda     (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                              
    !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
    !=====================================================================
    ! Eric Polizzi 2009-2013
    ! ====================================================================

    implicit none
!    !-------------------------------------
!#ifdef MPI
!    include 'mpif.h'
!#endif
!    !-------------------------------------
    include "f90_noruntime_interface.fi"
    integer :: ijob,N,M0
    complex :: Ze
    complex,dimension(N,*):: work,workc
    complex,dimension(M0,*):: zAq,zSq
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    real,dimension(*)  :: lambda
    complex,dimension(N,*):: q
    integer :: mode
    real,dimension(*) :: res
    integer :: info
    !! parameters
    real, Parameter :: pi=3.1415926535897932e0
    real, Parameter :: SONE=1.0e0, SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO), ZEROC=(SZERO,SZERO)
    real, parameter :: ba=-pi/2.0e0,ab=pi/2.0e0
    integer(8),parameter :: fout =6
    !! variable for FEAST
    integer :: i,e,j,k
    integer,dimension(4) :: iseed
    real :: theta,r,Emid
    complex :: jac,aux
    real ::xe,we ! Gauss-Legendre
    complex, dimension(:,:),pointer :: zSqo
    integer, dimension(:),pointer :: fpm_default
    logical :: testconv
    real :: trace
    !! Lapack variable (reduced system)
    character(len=1) :: JOBZ,UPLO
    real, dimension(:),pointer :: work_loc,work_loc2
    complex, dimension(:),pointer :: zwork_loc
    integer :: lwork_loc,info_lap,infoloc
    !! MPI compatibility variables
    integer :: rank,code,nb_procs,NEW_COMM_WORLD

    rank=0
    nb_procs=1
!    !----------------------------------------------
!#ifdef MPI
!    NEW_COMM_WORLD=fpm(9)
!    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
!    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
!    if (rank/=0) fpm(1)=0 ! comment only in rank 0 if any
!#endif
!    !---------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialization!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (ijob==-1) then 
       info=0 ! default value
       if (fpm(1)==1) then
          call wwrite_n(fout)
          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_s(fout, '*********** FEAST- BEGIN **********************')
          call wwrite_n(fout) 
          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout)
          call wwrite_s(fout, 'Routine CFEAST_H{}{}')
          call wwrite_n(fout)
!!!!!!!!!!!! Print the FEAST parameters which has been changed from default
          call wwrite_s(fout, 'List of input parameters fpm(1:64)-- if different from default')
          call wwrite_n(fout)

          call wallocate_1i(fpm_default,64,infoloc)
          call feastinit(fpm_default)
          do i=1,20
             if (fpm(i)/=fpm_default(i)) then
                call wwrite_s(fout, '   fpm(')
                call wwrite_i(fout, i)
                call wwrite_s(fout, ')=')
                call wwrite_i(fout, fpm(i))
                call wwrite_n(fout)
             endif
          enddo
          call wdeallocate_1i(fpm_default)

          call wwrite_s(fout, 'Search interval [')
          call wwrite_f(fout,Emin)
          call wwrite_s(fout, '; ')
          call wwrite_f(fout,Emax)
          call wwrite_s(fout, ']')
          call wwrite_n(fout)
       end if
       call check_fpm_input(fpm,info)
       call scheck_rci_input(Emin,Emax,M0,N,info)

       if (info/=0) fpm(21)=100 ! The End

!!!!!!!!!!!!!!!!
       IF (info==0) then
          if (fpm(1)==1) then
             call wwrite_s(fout, 'Size subspace')  
             call wwrite_t(fout) 
             call wwrite_i(fout,M0)
             call wwrite_n(fout)  
             call wwrite_s(fout, '#Loop | #Eig  |    Trace     |   Error-Trace |  Max-Residual') 
             call wwrite_n(fout)
          endif
          fpm(23)=min(M0,N) ! 'current M0' size (global value)
          fpm(25)=fpm(23) !! 'current M0' size (by default)
          fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
!          !----------------------------------------------
!#ifdef MPI
!          if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
!             fpm(25)=fpm(23)/nb_procs ! local size of 'current M0'
!             if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
!             fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
!          end if
!#endif
!          !-----------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          fpm(22)=fpm(2) ! only one half contour necessary here
          loop=0
          fpm(21)=1 ! prepare reentry
          if (fpm(5)==0) then !!! random vectors
             iseed=(/56,890,3456,2333/)
             call CLARNV(3,iseed,N*fpm(23),work(1,1))
          elseif (fpm(5)==1) then !!!!!! q is the initial guess
!             !----------------------------------------
!#ifdef MPI
!             work(1:N,1:fpm(23))=ZEROC 
!#endif
!             !------------------------------------------
             ijob=40 !! B*q=>work
             return
          end if
       end IF ! info=0
!!!!!!!!!!!!!!
    end if   !ijob=-1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    IF (fpm(21)==1) then !! we initialize a new contour integration
!       !------------------------------------------------------------------------
!#ifdef MPI
!       if ((loop>0).or.(fpm(5)==1)) then        
!          if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work,N*fpm(23),MPI_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
!       end if
!#endif
!       !------------------------------------------------------------------------
       q(1:N,1:fpm(23))=ZEROC
       fpm(20)=1
       fpm(21)=2
       ijob=-2 ! just initialization 
    end IF


!!!!!!!!!!!!
    IF (fpm(21)==2) then !! we start or pursue the contour integration

       IF (info==0) then !! will end up checking info errors returned by FEAST drivers
          do e=fpm(20)+rank,fpm(22),nb_procs !!!! loop over the contour 

             if (ijob==-2) then !!Factorize the linear system (complex) (zS-A)
                call sset_point_gauss_legendre(fpm(2),e,xe,we) !! Gauss-points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0e0
                Emid=Emin+r
                Ze=Emid*ONEC+r*ONEC*wscos(theta)+r*(SZERO,SONE)*wssin(theta)                               
                fpm(20)=e-rank
                ijob=10 ! for fact
                if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor 
             endif

             if (ijob==10) then !!Solve the linear system (complex) (zS-A)q=v 
                call CLACPY( 'F', N, fpm(23),work , N, workc, N )
                ijob=11 ! for solve
                return
             endif

             if (ijob==11) then 
                !! summation 
                call sset_point_gauss_legendre(fpm(2),e,xe,we) !! Gauss-points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0e0 
                jac=(r*(SZERO,SONE)*wssin(theta)+ONEC*r*wscos(theta))
                aux=-(ONEC/(4.0e0))*we*jac
                call CAXPY(N*fpm(23),aux,workc,1,q,1)
                !!Explicit Factorization of the linear system (complex) (zS-A)^T 
                !!needed if driver not capable to exploit Factorization of (zS-A) for solving (zS-A)^Tq=v          
                ijob=20 ! for fact
                if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor
             endif

             if (ijob==20) then!!!! Solve the linear system (complex) (zS-A)^Tq=v  
                call CLACPY( 'F', N, fpm(23),work , N, workc, N )
                ijob=21 ! for solve with transpose
                return
             end if

             if (ijob==21) then
                !! summation
                call sset_point_gauss_legendre(fpm(2),e,xe,we) !! Gauss-points 
                theta=ba*xe+ab
                r=(Emax-Emin)/2.0e0 
                jac=(r*(SZERO,SONE)*wssin(theta)+ONEC*r*wscos(theta))
                aux=-(ONEC/(4.0e0))*we*conjg(jac)
                call CAXPY(N*fpm(23),aux,workc,1,q,1)
                ijob=-2 ! just for identification
             end if

          end do
       end IF !! info=0

!       !------------------------------------------------
!#ifdef MPI
!       call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!#endif
!       !-----------------------------------------------  
       if (info/=0) fpm(21)=100 ! the end

       if (info==0) then
          fpm(21)=4 
!          !------------------------------------------------
!#ifdef MPI
!          call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
!#endif
!          !-----------------------------------------------  
       end if
    end IF   ! fpm(21)==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ((fpm(21)==4).and.(fpm(14)==1)) then !! only q vectors has been computed and is returned
info=4
       if (info/=0) fpm(21)=100 ! The End
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq=eq Sq xq     with Aq=Q^TAQ Sq=Q^TAQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!form Aq=> Aq=Q^T A Q 
    if (fpm(21)==4) then
       fpm(21)=5 ! preparing reentry
       ijob=30 
       return  ! mat-vec A*q => work
    endif

    if (fpm(21)==5) then 
!       !------------------------------------------------
!#ifdef MPI 
!       zAq(1:M0,1:fpm(23))=ZEROC
!#endif
!       !-------------------------------------------------
       call CGEMM('C','N',fpm(23),fpm(25),N,ONEC,q(1,1),N,work(1,fpm(24)),N,ZEROC,zAq(1,fpm(24)),M0)
       fpm(21)=6
    endif

!!!!!!!!!form  Sq=> Sq=Q^T S Q
    if (fpm(21)==6) then 
       fpm(21)=7 ! preparing reenty
       ijob=40 
       return! mat-vec S*q => work
    end if
    if (fpm(21)==7) then
!       !------------------------------------------------
!#ifdef MPI
!       zSq(1:M0,1:fpm(23))=ZEROC
!#endif
!       !-------------------------------------------------
       call CGEMM('C','N',fpm(23),fpm(25),N,ONEC,q(1,1),N,work(1,fpm(24)),N,ZEROC,zSq(1,fpm(24)),M0)
       fpm(21)=8
    endif

    if (fpm(21)==8) then
!       !----------------------------------------!(zAq,zSq known to all processors) 
!#ifdef MPI 
!       if (fpm(23)/nb_procs>=1) then
!          call MPI_ALLREDUCE(MPI_IN_PLACE,zAq(1,1),M0*fpm(23),MPI_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
!          call MPI_ALLREDUCE(MPI_IN_PLACE,zSq(1,1),M0*fpm(23),MPI_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
!       end if
!#endif
!       !---------------------------------------
       if (fpm(12)==1) then ! customize eigenvalue solver
          fpm(21)=9 ! preparing reentry - could return new value of M0 in fpm(23) if reduced subspace is needed
          ijob=50
          return
       endif
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==8) then
       if (rank==0) then
          JOBZ='V'
          UPLO='L'
          info_lap=1 ! initialization
          i=1
          LWORK_LOC=2*fpm(23)-1 !! for lapack eig reduced system
          call wallocate_1c(zWORK_LOC,LWORK_LOC,infoloc)
          call wallocate_1s(WORK_LOC,3*fpm(23)-2,infoloc)
          if (infoloc/=0) info=-1
          do while ((info_lap/=0).and.(info==0))
             i=i+1
             if (i==10) info=-3 ! arbitrary maximum
             call wallocate_2c(zSqo,fpm(23),fpm(23),infoloc)
             if (infoloc/=0) info=-1

             ! # eigenvalues inside search (wait second contour)
             mode=0
             if (loop>0) then
                call CLACPY( 'F', fpm(23), fpm(23),zSq , M0, zSqo, fpm(23) )
                if (infoloc/=0) info=-1
                call CHEEV( 'N', UPLO, fpm(23), zSqo, fpm(23), lambda, zwork_loc,lwork_loc,work_loc, INFO_lap )
                do j=1,fpm(23)
                   if (sqrt(abs(lambda(j)))>=0.5e0) mode=mode+1
                enddo
             end if

             call CLACPY( 'F', fpm(23), fpm(23),zSq , M0, zSqo, fpm(23) )
             call CHEGV(1, JOBZ, UPLO, fpm(23), zAq, M0, zSqo, fpm(23), lambda, zWORK_loc,Lwork_loc, WORK_loc, INFO_lap)

             if ((info_lap<=fpm(23)).and.(info_lap/=0)) info=-3
             if (info_lap>fpm(23)) then !! zSqo is not spd (a posteriori resize subspace)
                fpm(23)=info_lap-fpm(23)-1
                if (fpm(1)==1) then
                   call wwrite_s(fout, 'Resize subspace')  
                   call wwrite_t(fout) 
                   call wwrite_i(fout,fpm(23))
                   call wwrite_n(fout)
                end if
             end if
             call wdeallocate_2c(zSqo)
          end do
          call wdeallocate_1c(zwork_loc)
          call wdeallocate_1s(work_loc)
       end if !(rank 0)
!       !-------------------------------- !(info common to all processors) 
!#ifdef MPI
!       call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!#endif
!       !--------------------------------
       if (info/=0) fpm(21)=100 ! the end

       if (info==0) then
          fpm(25)=fpm(23) !! current M0 size (by default) -- global
          fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
!          !---------------------------------------- !(zAq==> vectors, lambda and fpm(23), known to all processors) 
!#ifdef MPI
!          call MPI_BCAST(mode,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!          call MPI_BCAST(fpm(23),1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
!          call MPI_BCAST(zAq,M0*fpm(23),MPI_COMPLEX,0,NEW_COMM_WORLD,code)
!          call MPI_BCAST(lambda,fpm(23),MPI_REAL,0,NEW_COMM_WORLD,code)
!          if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
!             fpm(25)=fpm(23)/nb_procs ! local size of current M0
!             if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
!             fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
!          end if
!#endif
!          !-----------------------------
          fpm(21)=9
       end if
    end if !! fpm(21)=8





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==9) then
       call CLACPY( 'F', N, fpm(23),q , N, work, N ) 
       !! option - non shifted
!#ifdef MPI 
!          q(1:N,1:fpm(23))=ZEROC
!#endif
          call CGEMM('N','N',N,fpm(25),fpm(23),ONEC,work(1,1),N,zAq(1,fpm(24)),M0,ZEROC,q(1,fpm(24)),N)
!#ifdef MPI       
!          if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
!#endif
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residual ||AX-lambda*BX||
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==9) then
       fpm(21)=10 ! preparing reentry
       ijob=30 
       return  ! mat-vec A*q => work
    endif

    if (fpm(21)==10) then
       call CLACPY( 'F', N, fpm(25),work(1,fpm(24)), N, workc(1,fpm(24)), N )
       fpm(21)=11 ! preparing reentry
       ijob=40 
!       !----------------------------------------
!#ifdef MPI
!       work(1:N,1:fpm(23))=ZEROC  !! work is also needed for outer-loop if any 
!#endif
!       !------------------------------------------
       return  ! mat-vec S*q => work 
    endif

    if (fpm(21)==11) then
!       !----------------------------------------
!#ifdef MPI
!       res(1:fpm(23))=SZERO
!#endif
!       !-------- Absolute residual
       do i=fpm(24),fpm(24)+fpm(25)-1
          res(i)=sum(abs(workc(1:N,i)-lambda(i)*work(1:N,i)))/sum(abs(max(abs(Emin),abs(Emax))*(work(1:N,i))))
       end do
!       !----------------------------------------
!#ifdef MPI 
!       if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_REAL,MPI_SUM,NEW_COMM_WORLD,code)
!#endif
!       !-----------------------------------------
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Trace / count eigenvalues (remove spurious if any)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (fpm(21)==11) then
       k=0
       trace=SZERO
       theta=SZERO
       do i=1,fpm(23)
          if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
             k=k+1
             trace=trace+lambda(i)
             if (res(i)>theta) theta=res(i) ! max residual
          endif
       enddo

!!!!!!!!! remove spurious if any       
       if (loop==0) mode=k ! wait second contour

       if (mode<k) then
          do j=1,k-mode
             theta=SZERO
             e=1
             do i=1,fpm(23)
                if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
                   if (res(i)>theta) then
                      e=i
                      theta=res(i) !max
                   endif
                end if
             enddo
             trace=trace-lambda(e) ! update trace
             res(e)=-SONE !spurious
          end do
!!! max residual
          theta=SZERO
          do i=1,fpm(23)
             if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
                if (res(i)>theta) theta=res(i) 
             end if
          enddo
       end if

       if (mode==0) info=1  ! no mode detected in the interval
if (loop>0) then
       if ((mode==M0).and.(mode/=N)) info=3 !size subspace too small
endif
       if (info/=0) fpm(21)=100 ! The End
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check FEAST Convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

    if (fpm(21)==11) then
       testconv=.false. ! initialization

       !! trace convergence
       if (loop==0) epsout=SONE
       if (loop>0) then
          epsout=(abs(trace-epsout))/max(abs(Emin),abs(Emax))           
          if ((fpm(6)==0).and.(log10(epsout)<(-fpm(7)))) testconv=.true.
       end if

       !! residual convergence
       if ((fpm(6)==1).and.(log10(theta)<(-fpm(7)))) testconv=.true.

       !! Two loops minimum if spurious are found
       if ((loop<=1).and.(k>mode)) testconv=.false.

       if (fpm(1)==1) then
          call wwrite_i(fout,loop)
          call wwrite_t(fout) 
          call wwrite_i(fout,mode)
          call wwrite_t(fout)
          call wwrite_f(fout,trace)
          call wwrite_t(fout) 
          call wwrite_f(fout,epsout)
          call wwrite_t(fout) 
          call wwrite_f(fout,theta)
          call wwrite_n(fout) 
       end if

       if (.not.testconv) then
          epsout=trace
          if (loop==fpm(4)) then
             info=2 ! FEAST did not converge (#loop reaches maximum)
             testconv=.true. ! return final eigenvector anyway
          endif
       endif
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! FEAST exit IF Convergence - FEAST iteration IF NOT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

    if (fpm(21)==11) then

       if (.not.testconv) then !!! need FEAST refinement loop

!!!!!!!!!!! here work=S*q is known alreasy,  using FEAST-MPI work is already partionned 
          fpm(21)=1   ! prepare reentry- reloop (with contour integration)
          !fpm(21)=4 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
          loop=loop+1
          ijob=-2 ! do nothing
          return  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       else    !!!!!!! final eigenvectors/eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! reorder (shift) lambda, eigenvector and residual
   call CLACPY( 'F', N, fpm(23),q , N, work, N )
          call wallocate_1s(work_loc,fpm(23),infoloc)
          call wallocate_1s(work_loc2,fpm(23),infoloc)
          call SCOPY(fpm(23),lambda , 1, work_loc, 1 )
          call SCOPY(fpm(23),res , 1, work_loc2, 1 )
          q(1:N,1:fpm(23))=ZEROC
          e=0
          j=0
          k=0
          do i=1,fpm(23)
             if ((work_loc(i)>Emin).and.(work_loc(i)<Emax)) then ! inside the search interval
                if (work_loc2(i)/=-SONE) then ! spurious 
                   e=e+1
                   q(1:N,e)=work(1:N,i)
                   lambda(e)=work_loc(i)
                   res(e)=work_loc2(i)
                endif
             else
                j=j+1
                q(1:N,mode+j)=work(1:N,i)
                lambda(mode+j)=work_loc(i)
                res(mode+j)=work_loc2(i)
             end if
             if (work_loc2(i)==-SONE) then ! spurious at the end
                k=k+1
                q(1:N,fpm(23)-k+1)=work(1:N,i)
                lambda(fpm(23)-k+1)=work_loc(i)
                res(fpm(23)-k+1)=work_loc2(i)
             endif
          enddo
          call wdeallocate_1s(work_loc)
          call wdeallocate_1s(work_loc2)

!!!!!!!!!!!!!!!!!!!!!!!!!!
          M0=fpm(23)  ! update value of M0 (new subspace)
          fpm(21)=100 ! The End

       end if ! test convergence
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==100) then !! THE END (ijob=0) 
       ijob=0 !! exit FEAST

       if (fpm(1)==1) then !! Print  Information

          if (info>=200) then
             call wwrite_s(fout, 'PROBLEM with input parameters')
             call wwrite_n(fout) 
          end if

          if ((info>100).and.(info<200)) then
             call wwrite_s(fout, 'PROBLEM with FEAST array parameters')
             call wwrite_n(fout) 
          end if

          if (info==-3) then
             call wwrite_s(fout, 'ERROR with reduced system')  
             call wwrite_n(fout) 
          end if

          if (info==-2) then
             call wwrite_s(fout, 'ERROR from Inner Linear System Solver in FEAST driver')  
             call wwrite_n(fout) 
          end if

          if (info==-1) then
             call wwrite_s(fout, 'ERROR with Internal memory allocation')  
             call wwrite_n(fout) 
          end if

          if (info==1) then
             call wwrite_s(fout, '==>WARNING: No eigenvalue has been found in the proposed search interval')
             call wwrite_n(fout)
          endif

          if (info==3) then
             call wwrite_s(fout, '==>WARNING: Size subspace M0 too small')  
             call wwrite_n(fout)
          end if


          if (info==4) then
             call wwrite_s(fout, '==>WARNING: Only the subspace has been returned')  
             call wwrite_n(fout)
          end if


          if (info==2) then
             call wwrite_s(fout, '==>WARNING: FEAST did not converge "yet" (#loop reaches maximum allowed)')  
             call wwrite_n(fout)
          end if

          if (info==0) then
             call wwrite_s(fout, '==>FEAST has successfully converged (to desired tolerance)')  
             call wwrite_n(fout) 
          else
             call wwrite_s(fout, '==>INFO code = ') 
             call wwrite_i(fout,info)
             call wwrite_n(fout)
          end if


          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_s(fout, '*********** FEAST- END*************************')
          call wwrite_n(fout) 
          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_n(fout)
       endif
    end if

  end subroutine cfeast_hrci






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
  










