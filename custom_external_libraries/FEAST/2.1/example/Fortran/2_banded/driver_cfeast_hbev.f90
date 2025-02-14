!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! FEAST Driver banded example !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! solving Ax=ex with A complex Hermitian --- A banded matrix!!!!!!!
!!!!!!! by Eric Polizzi- 2009-2011!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program driver 

  implicit none

!!!!!!!!!!!!!!!!! Feast declaration variable
  integer,dimension(64) :: feastparam 
  real :: epsout
  integer :: loop
  character(len=1) :: UPLO='F';  ! 'L' or 'U' also fine

!!!!!!!!!!!!!!!!! Matrix declaration variable
  character(len=100) :: name
  integer :: n,nnz,kl,ku
  complex,dimension(:,:),allocatable :: A

!!!!!!!!!!!!!!!!! Others
  integer :: t1,t2,tim
  integer :: i,j,k
  integer :: M0,M,info
  real :: Emin,Emax,rea,img
  complex,dimension(:,:),allocatable :: X ! eigenvectors
  real,dimension(:),allocatable :: E,res ! eigenvalue+residual

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! read input file in coordinate format!!!!!!!
!!!!!!!!!!!!!!!!transform  in banded format directly !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  name='../../system2'

  !! find kl,ku (kl==ku here)
  kl=0
  open(10,file=trim(name),status='old')
  read(10,*) n,nnz
  do k=1,nnz
     read(10,*) i,j
     if (abs(i-j)>kl) kl=abs(i-j)
  end do
  close(10)
  ku=kl

  !! form the banded matrix
  allocate(A(1+kl+ku,N))
  A=(0.0,0.0)
  open(10,file=trim(name),status='old')
  read(10,*) N,nnz
  do k=1,nnz
     read(10,*) i,j,rea,img
  A(ku+1+(i-j),j)=rea*(1.0,0.0)+img*(0.0,1.0)
  end do
  close(10)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,'banded matrix -system2- size',n
  print *,'bandwidth',kl+ku+1


  call system_clock(t1,tim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! FEAST in banded format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! search interval [Emin,Emax] including M eigenpairs
  Emin=-0.35
  Emax= 0.23 
  M0=40 !! M0>=M

!!!!!!!!!!!!! ALLOCATE VARIABLE 
  allocate(e(1:M0))     ! Eigenvalue
  allocate(X(1:n,1:M0)) ! Eigenvectors
  allocate(res(1:M0))   ! Residual 

!!!!!!!!!!!!!  FEAST
  call feastinit(feastparam)
  feastparam(1)=1
  call cfeast_hbev(UPLO,N,kl,A,kl+ku+1,feastparam,epsout,loop,Emin,Emax,M0,E,X,M,res,info)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! POST-PROCESSING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call system_clock(t2,tim)
  print *,'FEAST OUTPUT INFO',info
  if (info==0) then
     print *,'*************************************************'
     print *,'************** REPORT ***************************'
     print *,'*************************************************'
     print *,'SIMULATION TIME',(t2-t1)*1.0d0/tim
     print *,'# Search interval [Emin,Emax]',Emin,Emax
     print *,'# mode found/subspace',M,M0
     print *,'# iterations',loop
     print *,'TRACE',sum(E(1:M))
     print *,'Relative error on the Trace',epsout
     print *,'Eigenvalues/Residuals'
     do i=1,M
        print *,i,E(i),res(i)
     enddo
  endif

end program driver



