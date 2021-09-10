!> solve (K-z_i.M)x=b, where K and M are real matrices of size ndof,
!> z_i (i=1,...,nsfts) is a complex number. One may have multiple right-hand size.
module module_scocg
  use module_globals
  use module_random_numbers
  implicit none

  type LinearSystems
     integer :: ndof !< DoF
     integer :: nnzm !< Number of non-zero elements in K matrix
     integer :: nnzk !< Number of non-zero elements in M matrix
     integer :: nshfts !< number of shifts
     integer :: nrhs !< number of right-hand sides
     integer,allocatable :: ik(:), jk(:) !< pointers for K matrix
     integer,allocatable :: im(:), jm(:) !< pointers for M matrix
     real(8),allocatable :: kmat(:) !< K matrix
     real(8),allocatable :: mmat(:) !< M matrix
     real(8),allocatable :: fmmat(:) !< M^-1 matrix
     complex(8),allocatable :: amat(:) !< coefficient matrix (for validation)
     complex(8),allocatable :: shfts(:) !< shfts(j): j-th shift
     complex(8),allocatable :: xvec(:,:,:) !< xvec(:,i): i-th/nrhs and j-th/nshfts sol .
     complex(8),allocatable :: bvec(:,:) !< bvec(:,i): i-th right-hand side
  end type LinearSystems

  integer,pointer :: ndof !< DoF
  integer,pointer :: nnzk !< Number of non-zero elements in K matrix
  integer,pointer :: nnzm !< Number of non-zero elements in M matrix
  integer,pointer :: nshfts !< number of shifts
  integer,pointer :: nrhs !< number of right-hand sides
  integer,pointer :: ik(:), jk(:) !< pointers for K matrix
  integer,pointer :: im(:), jm(:) !< pointers for M matrix
  real(8),pointer :: kmat(:) !< K matrix
  real(8),pointer :: mmat(:) !< M matrix
  real(8),pointer :: fmmat(:) !< M^-1 matrix
  complex(8),pointer :: shfts(:) !< shfts(j): j-th shift
  complex(8),pointer :: xvec(:,:,:) !< xvec(:,i,j): i-th/nrhs and j-th/nshfts sol .
  complex(8),pointer :: bvec(:,:) !< bvec(:,i): i-th right-hand side

  ! pardiso
  integer(8) :: pt(64)
  integer :: maxfct, mnum, mtype, phase, error, msglvl
  integer :: iparm(64)
  integer :: i, idum(1)
  real(8) :: ddum(1)
  
contains

  !> quick sort for ia
  recursive subroutine quicksort(ia,ja,a,first,last)
    implicit none
    integer :: first, last
    integer :: i, j, k, t
    integer :: ia(*), ja(*)
    real(8) :: a(*)
    integer :: l
    real(8) :: cmp

    l=first+last
    k=ia(l/2)
    i=first
    j=last
    do
       do while (ia(i)<k)
          i=i+1
       end do
       do while (k<ia(j))
          j=j-1
       end do
       if (i >= j) exit
       t=ia(i)
       ia(i)=ia(j)
       ia(j)=t
       t=ja(i)
       ja(i)=ja(j)
       ja(j)=t
       cmp=a(i)
       a(i)=a(j)
       a(j)=cmp
       i=i+1
       j=j-1
    end do
    if(first < i-1) call quicksort(ia,ja,a,first,i-1)
    if(j+1 < last)  call quicksort(ia,ja,a,j+1,last)
  end subroutine quicksort

  !> quick sort for ja
  recursive subroutine quicksort2(ia,a,first,last)
    implicit none
    integer :: first, last
    integer :: i, j, k, t
    integer :: ia(*)
    real(8) :: a(*)
    integer :: l
    real(8) :: cmp

    l=first+last
    k=ia(l/2)
    i=first
    j=last
    do
       do while (ia(i)<k)
          i=i+1
       end do
       do while (k<ia(j))
          j=j-1
       end do
       if (i >= j) exit
       t=ia(i)
       ia(i)=ia(j)
       ia(j)=t
       cmp=a(i)
       a(i)=a(j)
       a(j)=cmp
       i=i+1
       j=j-1
    end do
    if(first < i-1) call quicksort2(ia,a,first,i-1)
    if(j+1 < last)  call quicksort2(ia,a,j+1,last)
  end subroutine quicksort2
  
  !> initialise
  subroutine init_scocg(fk,fm,fs,ls_)
    character(len=32),intent(in) :: fk, fm
    character(len=32),intent(in) :: fs
    type(LinearSystems),intent(inout),target :: ls_

    integer,allocatable :: iia(:)
    integer :: itmp, ix, iy, i
    real(8) :: tmp
    
    ndof=>ls_%ndof
    nnzk=>ls_%nnzk
    nnzm=>ls_%nnzm
    nshfts=>ls_%nshfts
    nrhs=>ls_%nrhs

    ! load matrix K
    open(1,file=fk)
    read(1,*); read(1,*) ! skip headers
    read(1,*) ndof, itmp, nnzk
    allocate(iia(nnzk))
    allocate(ls_%ik(ndof+1)); ik=>ls_%ik
    allocate(ls_%jk(nnzk)); jk=>ls_%jk
    allocate(ls_%kmat(nnzk)); kmat=>ls_%kmat
    do i=1,nnzk
       read(1,*) ix, iy, tmp
       iia(i)=ix+1
       jk(i)=iy+1
       kmat(i)=tmp
    end do
    close(1)

    call quicksort(iia,jk,kmat,1,nnzk)
    ik(1)=1
    itmp=1
    do i=1,ndof-1
       do while(iia(itmp).eq.i)
          itmp=itmp+1
       end do
       ik(i+1)=itmp
    end do
    ik(ndof+1)=nnzk+1
    do i=1,ndof
       call quicksort2(jk,kmat,ik(i),ik(i+1)-1)
    end do
    deallocate(iia)
    
    ! load mass matrix M
    open(1,file=fm)
    read(1,*); read(1,*) ! skip headers
    read(1,*) ndof, itmp, nnzm
    allocate(iia(nnzm))
    allocate(ls_%im(ndof+1)); im=>ls_%im
    allocate(ls_%jm(nnzm)); jm=>ls_%jm
    allocate(ls_%mmat(nnzm)); mmat=>ls_%mmat
    allocate(ls_%fmmat(nnzm)); fmmat=>ls_%fmmat
    do i=1,nnzm
       read(1,*) ix, iy, tmp
       iia(i)=ix+1
       jm(i)=iy+1
       mmat(i)=tmp
    end do
    close(1)

    call quicksort(iia,jm,mmat,1,nnzm)
    im(1)=1
    itmp=1
    do i=1,ndof-1
       do while(iia(itmp).eq.i)
          itmp=itmp+1
       end do
       im(i+1)=itmp
    end do
    im(ndof+1)=nnzm+1
    do i=1,ndof
       call quicksort2(jm,mmat,im(i),im(i+1)-1)
    end do
    fmmat=mmat
    
    deallocate(iia)

    ! load shifts
    open(1,file=fs)
    read(1,*) nshfts
    allocate(ls_%shfts(nshfts)); shfts=>ls_%shfts
    read(1,*) shfts
    close(1)

    ! right-hand sides
    nrhs=20
    allocate(ls_%bvec(ndof,nrhs)); bvec=>ls_%bvec
    allocate(ls_%xvec(ndof,nrhs,nshfts)); xvec=>ls_%xvec

    call gen_random_numbers_complex(ndof*nrhs,bvec)
    
  end subroutine init_scocg

  !> destroy
  subroutine uninit_scocg(ls_)
    type(LinearSystems),intent(inout),target :: ls_

    nullify(ndof,nnzk,nnzm,nshfts,nrhs)
    deallocate(ls_%ik,ls_%jk,ls_%kmat,ls_%im,ls_%jm,ls_%mmat,ls_%fmmat,ls_%shfts)
    nullify(ik,jk,im,jk,kmat,mmat,fmmat,shfts)
    deallocate(ls_%bvec,ls_%xvec)
    nullify(bvec,xvec)
    
  end subroutine uninit_scocg

  subroutine lapack_complex
    integer :: i, j

    integer :: k, icnt, info, ipiv(ndof)
    complex(8),allocatable :: am(:,:)
    
    allocate(am(ndof,ndof))

    do i=1,nshfts
       am(:,:)=zero
       icnt=1
       do j=1,ndof
          do k=ik(j),ik(j+1)-1
             am(j,jk(icnt))=kmat(icnt)
             icnt=icnt+1
          end do
       end do
       icnt=1
       do j=1,ndof
          do k=im(j),im(j+1)-1
             am(j,jm(icnt))=am(j,jm(icnt))-shfts(i)*mmat(icnt)
             icnt=icnt+1
          end do
       end do
       xvec(:,:,i)=bvec
       call zgesv(ndof,nrhs,am,ndof,ipiv,xvec(:,:,i),ndof,info)
       do j=1,ndof
          write(9+i,*) real(xvec(j,3,i)), aimag(xvec(j,3,i))
       end do
    end do

    
    deallocate(am)
    
  end subroutine lapack_complex

  !> the main
  subroutine scocg
    integer :: irhs, ishf, ist, i, j, k, icnt
    real(8) :: trvec(ndof), tivec(ndof), qrvec(ndof), qivec(ndof)
    complex(8) :: pvec(ndof,nshfts), rvec(ndof), tvec(ndof), qvec(ndof) !!!!!!!!!!!!!!!!
    complex(8) :: alp(2,nshfts), bet(2,nshfts), c(5,nshfts) !!!!!!!!!!!!

    complex(8) :: bunbo, bunsi, sig

    ! precon. invert M
    maxfct = 1 
    mnum = 1
    iparm(:)=0
    iparm(1) = 1 ! no solver default
    iparm(2) = 2 ! fill-in reordering from METIS
    iparm(4) = 0 ! no iterative-direct algorithm
    iparm(5) = 0 ! no user fill-in reducing permutation
    iparm(6) = 0 ! =0 solution on the first n components of x
    iparm(8) = 0 ! numbers of iterative refinement steps
    iparm(10) = 13 ! perturb the pivot elements with 1E-13
    iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
    iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
    iparm(14) = 0 ! Output: number of perturbed pivots
    iparm(18) = 1 ! Output: number of nonzeros in the factor LU
    iparm(19) = 1 ! Output: Mflops for LU factorization
    iparm(20) = 0 ! Output: Numbers of CG Iterations

    error  = 0 ! initialize error flag
    msglvl = 0 ! print statistical information
    mtype  = 11 ! real, general

    pt(:)=0

    phase=12 ! Analysis, numerical factorization
    call pardiso (pt, maxfct, mnum, mtype, phase, ndof, fmmat, im, jm, &
         idum, 1, iparm, msglvl, ddum, ddum, error)
    if (error /= 0) then
       write(*,*) 'The following ERROR was detected: ', error
       stop 
    end if

    do irhs=3,3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       trvec(:)=real(bvec(:,irhs))
       tivec(:)=aimag(bvec(:,irhs))
       qrvec=zero
       qivec=zero
       phase=33
       call pardiso (pt, maxfct, mnum, mtype, phase, ndof, fmmat, im, jm, &
            idum, 1, iparm, msglvl, trvec, qrvec, error)
       if (error /= 0) then
          write(*,*) 'the following error was detected (doko): ', error
          stop
       endif
       call pardiso (pt, maxfct, mnum, mtype, phase, ndof, fmmat, im, jm, &
            idum, 1, iparm, msglvl, tivec, qivec, error)
       if (error /= 0) then
          write(*,*) 'the following error was detected (koko): ', error
          stop
       endif
       qvec=cmplx(qrvec,qivec,kind(1.d0))
       rvec(:)=bvec(:,irhs)
       do ishf=1,nshfts
          pvec(:,ishf)=qvec(:)
          xvec(:,irhs,ishf)=zero
       end do
       c(1:3,:)=one
       alp(:,:)=one
       bet(:,:)=zero

       ! CG loop
       do ist=1,ndof

          ! bunsi=r^t.M^-1r
          trvec(:)=real(rvec(:))
          tivec(:)=aimag(rvec(:))
          qrvec=zero
          qivec=zero
          phase=33
          call pardiso (pt, maxfct, mnum, mtype, phase, ndof, fmmat, im, jm, &
               idum, 1, iparm, msglvl, trvec, qrvec, error)
          if (error /= 0) then
             write(*,*) 'the following error was detected (doko): ', error
             stop
          endif
          call pardiso (pt, maxfct, mnum, mtype, phase, ndof, fmmat, im, jm, &
               idum, 1, iparm, msglvl, tivec, qivec, error)
          if (error /= 0) then
             write(*,*) 'the following error was detected (koko): ', error
             stop
          endif
          qvec=cmplx(qrvec,qivec,kind(1.d0))
          
          bunsi=zero
          do j=1,ndof
             bunsi=bunsi+rvec(j)*qvec(j)
          end do

          ! t=Kp
          tvec(:)=zero
          icnt=1
          do i=1,ndof
             do k=ik(i),ik(i+1)-1
                tvec(i)=tvec(i)+kmat(icnt)*pvec(jk(icnt),1)
                icnt=icnt+1
             end do
          end do
          icnt=1
          do i=1,ndof
             do k=im(i),im(i+1)-1
                tvec(i)=tvec(i)-shfts(1)*mmat(icnt)*pvec(jm(icnt),1)
                icnt=icnt+1
             end do
          end do

          !qvec(:)=tvec(:)-shfts(1)*pvec(:,1)
          
          ! bunbo=p^t.(M^-1.K)p
          bunbo=zero
          do j=1,ndof
             bunbo=bunbo+pvec(j,1)*tvec(j)
             !write(101,*) pvec(j,1)
          end do

          
          ! update alpha
          alp(2,1)=alp(1,1) ! keep the previous one
          alp(1,1)=bunsi/bunbo ! and update
          xvec(:,irhs,1)=xvec(:,irhs,1)+alp(1,1)*pvec(:,1)
          
          ! update rvec
          rvec(:)=rvec(:)-alp(1,1)*tvec(:)

          ! update beta
          bunbo=bunsi !beta's bunbo = alpha's bunsi

          trvec(:)=real(rvec(:)) !ika, bunsi
          tivec(:)=aimag(rvec(:))
          qrvec=zero
          qivec=zero
          phase=33
          call pardiso (pt, maxfct, mnum, mtype, phase, ndof, fmmat, im, jm, &
               idum, 1, iparm, msglvl, trvec, qrvec, error)
          if (error /= 0) then
             write(*,*) 'the following error was detected (doko): ', error
             stop
          endif
          call pardiso (pt, maxfct, mnum, mtype, phase, ndof, fmmat, im, jm, &
               idum, 1, iparm, msglvl, tivec, qivec, error)
          if (error /= 0) then
             write(*,*) 'the following error was detected (koko): ', error
             stop
          endif
          qvec=cmplx(qrvec,qivec,kind(1.d0))
          bunsi=zero
          do j=1,ndof
             bunsi=bunsi+rvec(j)*qvec(j)
          end do
          bet(2,1)=bet(1,1)
          bet(1,1)=bunsi/bunbo

          pvec(:,1)=qvec(:)+bet(1,1)*pvec(:,1)
          
          do ishf=2,nshfts
             sig=shfts(1)-shfts(ishf)
             c(3,ishf)=c(2,ishf)
             c(2,ishf)=c(1,ishf)
             c(1,ishf)=(1.d0+alp(1,1)*sig)*c(2,ishf)+(c(2,ishf)-c(3,ishf))*bet(2,1)*alp(1,1)/alp(2,1)
             alp(1,ishf)=c(2,ishf)/c(1,ishf)*alp(1,1)
             xvec(:,irhs,ishf)=xvec(:,irhs,ishf)+alp(1,ishf)*pvec(:,ishf)
             bet(1,ishf)=(c(2,ishf)/c(1,ishf))**2*bet(1,1)
             pvec(:,ishf)=qvec(:)/c(1,ishf)+bet(1,ishf)*pvec(:,ishf)
          end do
          if(real(dot_product(rvec,rvec)).le.real(dot_product(bvec(:,irhs),bvec(:,irhs)))*1.d-32) exit
       end do

       do i=1,nshfts
          do j=1,ndof
             write(29+i,*) real(xvec(j,3,i)), aimag(xvec(j,3,i))
          end do
       end do
    end do
    
  end subroutine scocg
  
end module module_scocg
