C -*- Mode: Fortran; -*-
C
C  (C) 2004 by Argonne National Laboratory.
C      See COPYRIGHT in top-level directory.
C
C  This tests the null copy function (returns flag false; thus the
C  attribute should not be propagated to a dup'ed communicator
C  This is must like the test in commattr2f
C
      program main
      implicit none
      include 'mpif.h'
      integer errs, ierr
      include 'attraints.h'
      integer comm1, comm2
      integer keyval
      logical flag
C
C The only difference between the MPI-2 and MPI-1 attribute caching
C routines in Fortran is that the take an address-sized integer
C instead of a simple integer.  These still are not pointers,
C so the values are still just integers.
C
      errs      = 0
      call mtest_init( ierr )
      call mpi_comm_dup( MPI_COMM_WORLD, comm1, ierr )
C
      extrastate = 1001
      call mpi_comm_create_keyval( MPI_COMM_NULL_COPY_FN,
     &                             MPI_COMM_NULL_DELETE_FN, keyval,
     &                             extrastate, ierr )
      flag = .true.
      call mpi_comm_get_attr( comm1, keyval, valout, flag, ierr )
      if (flag) then
         errs = errs + 1
         print *, ' get attr returned true when no attr set'
      endif

C Test the null copy function
      valin = 5001
      call mpi_comm_set_attr( comm1, keyval, valin, ierr )
      call mpi_comm_dup( comm1, comm2, ierr )
C Because we set NULL_COPY_FN, the attribute should not
C appear on the dup'ed communicator
      flag = .false.
      call mpi_comm_get_attr( comm1, keyval, valout, flag, ierr )
      if (valout .ne. 5001) then
         errs = errs + 1
         print *, 'Unexpected output value in comm ', valout
      endif
      flag = .true.
      call mpi_comm_get_attr( comm2, keyval, valout, flag, ierr )
      if (flag) then
         errs = errs + 1
         print *, ' Attribute incorrectly present on dup communicator'
      endif
C Test the delete function
      call mpi_comm_free( comm2, ierr )
C
C Test the attr delete function
      call mpi_comm_dup( comm1, comm2, ierr )
      valin      = 6001
      extrastate = 1001
      call mpi_comm_set_attr( comm2, keyval, valin, ierr )
      call mpi_comm_delete_attr( comm2, keyval, ierr )
      flag = .true.
      call mpi_comm_get_attr( comm2, keyval, valout, flag, ierr )
      if (flag) then
         errs = errs + 1
         print *, ' Delete_attr did not delete attribute'
      endif
      call mpi_comm_free( comm2, ierr )
C
      ierr = -1
      call mpi_comm_free_keyval( keyval, ierr )
      if (ierr .ne. MPI_SUCCESS) then
         errs = errs + 1
         call mtestprinterror( ierr )
      endif
      call mpi_comm_free( comm1, ierr )

      call mtest_finalize( errs )
      call mpi_finalize( ierr )
      end
