program mk_identity
  integer,parameter :: n=1470
  write(1,*)
  write(1,*)
  write(1,*) n,n,n
  do i=1,n
     write(1,*) i-1,i-1,1.d0
  end do
end program mk_identity
