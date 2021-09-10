module module_random_numbers
  implicit none
contains

  subroutine gen_random_numbers_complex(n,rand)
    integer,intent(in) :: n
    complex(8),intent(out) :: rand(n)

    integer :: i, seedsize
    integer,allocatable :: seed(:)
    integer :: zikoku
    real(8) :: ransu(2*n)
    
    call system_clock(count=zikoku)
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    seed=zikoku
    call random_seed(put=seed)
    call random_number(ransu)
    deallocate(seed)
    do i=1,n
       !rand(i)=cmplx(ransu(2*i-1),ransu(2*i),kind(1.d0))
       !rand(i)=1.d0
       rand(i)=dble(i)
    end do
    
  end subroutine gen_random_numbers_complex
  
end module module_random_numbers
