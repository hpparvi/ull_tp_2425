PROGRAM Parallel_Hello_World
  USE OMP_LIB
  implicit none
  integer :: nt = 1, tid = 0
  logical :: compiled_with_openmp = .false.

  !$OMP PARALLEL PRIVATE(nt, tid)
   !$nt = omp_get_num_threads()
  !$tid = OMP_GET_THREAD_NUM()
  !$ compiled_with_openmp = .true.

  PRINT *, &
	& "Inside parallel, thread ", &
	& tid, " of ", nt

  !$OMP END PARALLEL

  if (compiled_with_openmp) then
     write(*,*) 'OpenMP used'
  else
     write(*,*) 'OpenMP NOT used'
  end if

END PROGRAM Parallel_Hello_World
