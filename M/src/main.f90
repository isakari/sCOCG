program main
  use module_globals
  use module_scocg
  implicit none

  type(LinearSystems) :: ls

  character(len=32) :: fk, fm, fs
  fk="k10wp.txt"
  fm="m10wp.txt"
  fs="shifts.txt"  
  
  !k(m).txt contains the non-zero components of the matrix K(M) in the COO format
  call init_scocg(fk,fm,fs,ls)
  call lapack_complex
  call scocg
  call uninit_scocg(ls)
  
end program main
