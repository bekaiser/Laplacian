! ****************************************************************************

! FILE:         mpi_scatter_gather_2D.f90
! AUTHOR:       Bryan Kaiser
! DATE:         6/30/17
! DESCRIPTION:  A 2D array is split into 4 2D arrays that are scattered to 
!               4 ranks, where 1 is added to each value of each array. Then 
!               the 2D arrays on each rank are gathered back to the master 
!               rank.
! REQUIRED:     OpenMPI/1.10.5
! CREDIT:       This script was written using mpi_scatter.f by Blaise Barney 
!               (04/02/05) as a template.

! ****************************************************************************
program mpi_scatter_gather_2D
include 'mpif.h'

integer, parameter :: Np = 4, Nx = 8, Ny = 5
integer :: i, j, k, Nprocs, rank, ierror
real*8, dimension(Ny,Nx) :: u, ug  ! double precision
real*8, dimension(Ny,2) :: u0      ! double precision
integer, dimension(Np) :: uc = int(2*Ny) 
integer, dimension(Np) :: us

! Np = the number of processes this script was written for (equivalent to the 
!      number of columns of u). This is just the number of columns.
! Nx = the number of grid values in x (horizontal)
! Ny = the number of grid values in y (vertical). FORTRAN is column major, so 
!      the address of each cell in the array starts at the top of the left-most 
!      column and proceeds down the column.
! u  = the input array, which is given values corresponding to the address of 
!      each cell in the array.
! ug = the output array, equivalent to the input array + 1. The addition is 
!      performed by giving each rank one column.
! u0 = the dummy vector each rank is given to then populate.
! uc = a vector with the number of grid points that each rank is given
! us = a vector that tells each rank which address in the u array to take as 
!      the first value in the vector on that rank.
! Nprocs = the number of ranks that MPI recognizes. 
! rank = the processor id, from 0 to 3 for 4 processors.


call MPI_INIT(ierror) ! Start MPI
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, Nprocs, ierror)

! Do this on rank 0
if ( rank == 0 ) then
  print *, "Scatter this array from master rank:"
  do j = 1,Ny 
    do i = 1,Nx
      u(j,i) = j + (i-1)*Ny
    end do
    print *, int(u(j,:))
  end do
  us(1) = 0 ! the first rank vector starts with the address: u(1) + us(1) = 1
  us(2) = 10 ! the second rank vector starts with the address: u(1) + us(2) = 6
  us(3) = 20 ! the third rank vector starts with the address: u(1) + us(3) = 11
  us(4) = 30 ! the fourth rank vector starts with the address: u(1) + us(4) = 16
end if

if ( Nprocs == Np ) then 
  ! if the number of mpirun processes matches the number of processes this 
  ! script is set up for, then go ahead and run:
  
  ! Send vectors of length 5 from the master rank (rank 0) to each of the 4 ranks 
  ! (either scatter or scatter vector will do):

  !call MPI_SCATTER( u, int(2*Ny), MPI_DOUBLE, u0, int(2*Ny), MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
  call MPI_SCATTERV( u, uc, us, MPI_DOUBLE, u0, int(2*Ny), MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)

  ! Each of the 4 processes gets a vector length 5 of the 5x4 array u, adds 1:
  do i = 0,Nprocs-1 ! looping over processes
    if ( rank == i ) then
      !print *, "rank =",i
      !print *, int(u0)
      do j = 1,Ny
        do k = 1,2
          u0(j,k) = u0(j,k) + 1.0 ! add 1 to each value
        end do
      end do
    end if 
  end do 

 ! Send vectors of length 5 from each of the 4 ranks back to the master rank (rank 0) 
 ! (again either gather or gather vector will do this):

 !call MPI_GATHER( u0, int(2*Ny), MPI_DOUBLE, ug, int(2*Ny), MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
 call MPI_GATHERV( u0, int(2*Ny), MPI_DOUBLE, ug, uc, us, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)

 ! Print out the resulting array ug from the master rank:
 if ( rank == 0 ) then
   print *, "Add 1 to each component of each of the 4 vectors scattered to 4 ranks..."
   print *, "...then gather back to master rank:"
   do j = 1,Ny 
    print *, int(ug(j,:))
  end do
 end if

end if ! if Nprocs = Np

call MPI_FINALIZE(ierr) ! Close MPI

end program mpi_scatter_gather_2D
