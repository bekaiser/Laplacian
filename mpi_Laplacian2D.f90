! ****************************************************************************

! FILE:         mpi_Laplacian2D.f90
! AUTHOR:       Bryan Kaiser
! DATE:         7/7/2017
! DESCRIPTION:  This script computes the Laplacian of the variable u using 
!               second-order accurate central finite differences on Np number 
!               of processors. The domain is divided up in the x direction, so
!               the number of cells in x (Nx) must be divisible by the number
!               of MPI ranks (Np). The is one layer of halo cells surrounding 
!               the domain and the domain chunks for parallel computation. 

! REQUIRED:     OpenMPI/1.10.5 (at least)
      

! ****************************************************************************
program mpi_Laplacian2D

implicit none
include 'mpif.h'

! grid parameters
integer, parameter :: Nx = 1000, Ny = 1000 ! grid index dimensions
integer, parameter :: ih = 1, Np = 4 ! number of ghost cell layers and number of ranks (processors)
integer, parameter :: N = int(Nx*Ny) ! total number of computed grid points without ghost cells
integer, parameter :: Nxg = int(Nx+ih*2), Nyg = int(Ny+ih*2) ! grid points with ghost cells
integer, parameter :: Ng = int(Nxg*Nyg) ! total number of grid points without ghost cells
real*8, dimension(Nx+ih*2) :: x1 = 0.0 ! 1d grid of the x locations of grid points
real*8, dimension(Ny+ih*2) :: y1 = 0.0 ! 1d grid of the y locations of grid points
real*8, dimension(Ny+2*ih,Nx+2*ih) :: x, y ! 2d grid of x,y locations of grid points
real*8 :: dx, dy, kx, ky ! uniform grid spacings / signal wavenumbers 
real*8 :: Lx = 1.0, Ly = 2.0 ! grid dimensions in x,y
integer :: i, j, k ! loop indices

! input signal parameters
real*8, dimension(Ny,Nx) :: lap_signal = 0.0, u_signal = 0.0 ! analytical signals without ghost cells
real*8, parameter :: pi = 3.14159265358979323846264338 

! global variables
real*8, dimension(Ny+ih*2,Nx+ih*2) :: u_global = 0.0 ! velocity u field, with ghost cells
real*8, dimension(Ny,Nx) :: lap_global ! laplacian of the u field, without ghost cells

! local variables and parameters
integer, parameter :: Nxp = int(Nx/Np) ! number of cells in x that each rank receives
real*8, dimension(Ny+ih*2,Nxp+ih*2) :: u_local ! 1/4 of domain without ghost cells
real*8, dimension(Ny,Nxp) :: lap_local ! 1/4 of domain without ghost cells
integer, parameter :: Nch = int((Ny+ih*2)*(Nxp+ih*2)) ! chunk size with ghost cells
integer, parameter :: Nc = int(Ny*Nxp) ! chunk size without ghost cells

! script output
real*8, dimension(Nx*Ny) :: sig_out = 0.0, soln_out = 0.0 ! output vectors without ghost cells

! MPI variables
integer, dimension(Np) :: request, statuses
integer :: rank, Nprocs, ierror
integer, dimension(Np) :: N_start, Nh_start ! chunk start points without/with ghost cells (halos)
integer, dimension(Np) :: N_count = Nc, Nh_count = Nch ! total cell counts of each chunk without/with ghost cells (halos)
real*8 :: seconds ! MPI timer 

call MPI_INIT(ierror) ! Start MPI
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, Nprocs, ierror)

if ( rank == 0 ) then ! do on the master rank (rank 0):

  ! 1) construct the grid
  dx = real(Lx/real(Nx))
  dy = real(Ly/real(Ny))
  !print *, dx, dy
  do i = 1,Nx+2*ih
    x1(i) = dx*real(i-1) - dx*0.5
  end do
  do j = 1,Ny+2*ih
    y1(j) = dy*real(j-1) - dy*0.5
  end do
  do i = 1,Nx+2*ih
    do j = 1,Ny+2*ih
      x(j,i) = x1(i)
      y(j,i) = y1(j)
    end do
  end do

  ! 2) construct the signal (u_signal) and the analytical laplacian (lap_signal)
  kx = (2.0*pi)/Lx
  ky = (2.0*pi)/Ly
  u_signal = cos(x(2:Ny+1,2:Nx+1)*kx)*cos(y(2:Ny+1,2:Nx+1)*ky);
  lap_signal = -(kx**2.+ky**2.)*cos(x(2:Ny+1,2:Nx+1)*kx)*cos(y(2:Ny+1,2:Nx+1)*ky);

  ! 3) create a global array of the signal with periodic boundary conditions:
  do i = 1,Nx ! NOTE: explicit loops are faster than implicit!
    do j = 1,Ny 
      u_global(j+1,i+1) = u_signal(j,i)
    end do
  end do

  ! Ghost cells: periodic copy values (for practical applications) 
  !do i = 1,Nx 
  !  u_global(1,i+1) = u_signal(Ny,i)
  !  u_global(Ny+2,i+1) = u_signal(1,i)
  !end do
  !do j = 1,Ny
  !  u_global(j+1,1) = u_signal(j,Nx)
  !  u_global(j+1,Nx+2) = u_signal(j,1)
  !end do

  ! Ghost cells: analytical values
  do i = 2,Nx+1
    u_global(1,i) = cos(x(1,i)*kx)*cos(y(1,i)*ky);
    u_global(Ny+2,i) = cos(x(Ny+2,i)*kx)*cos(y(Ny+2,i)*ky);
  end do
  do j = 2,Ny+1
    u_global(j,1) = cos(x(j,1)*kx)*cos(y(j,1)*ky);
    u_global(j,Nx+2) = cos(x(j,Nx+2)*kx)*cos(y(j,Nx+2)*ky);
  end do
 
  ! Ghost cells: corners
  !u_global(1,1) = u_signal(Ny,Nx)
  !u_global(Ny+2,Nx+2) = u_signal(1,1)
  !u_global(Ny+2,1) = u_signal(1,Nx)
  !u_global(1,Nx+2) = u_signal(Ny,1)
 
  ! 4) Specify the starting indices for each chunk of the signal passed to 
  !    each rank.

  ! MPI start indices, with/without halo cells
  do i = 1,Np
    Nh_start(i) = int((i-1)*Nxp*(Ny+2*ih))
    N_start(i) = int((i-1)*Nc)
  end do
  
end if ! do this on the master rank (rank 0)

if ( Nprocs == Np ) then ! if Np is the same number of ranks as you specified with mpirun -n

  seconds = MPI_Wtime ( ) ! start timer

  ! 5) Pass the signal, the empty solution arrays, and finite difference variables to each rank:
  call MPI_SCATTERV( u_global, Nh_count, Nh_start, MPI_DOUBLE, u_local, Nch, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
  call MPI_BCAST( lap_local , Nc, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror) 
  call MPI_BCAST( dx , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror) 
  call MPI_BCAST( dy , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)   

  ! 6) Compute the Laplacian on each chunk of the domain:
  do k = 0,Nprocs-1
    if ( rank == k ) then
      do i = (1+ih),(Nxp+ih) 
        do j = (1+ih),(Ny+ih) 
          lap_local(j-1,i-1) = (u_local(j,i-1)-2.0*u_local(j,i)+u_local(j,i+1))/(dx**2.0) + & 
                               (u_local(j-1,i)-2.0*u_local(j,i)+u_local(j+1,i))/(dy**2.0)
        end do
      end do
    end if
  end do 

 ! 7) Collect the computed Laplacian chunks on the master rank to construct the solution over the entire domain: 
 call MPI_WAIT(request, statuses, ierror)
 call MPI_GATHERV( lap_local, Nc, MPI_DOUBLE, lap_global, N_count, N_start, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
 call MPI_WAIT(request, statuses, ierror)

 seconds = MPI_Wtime ( ) - seconds ! end timer
 print *, seconds

 if ( rank == 0 ) then ! do on the master rank (rank 0):

  ! 8) Write the solution to an ASCII file for plotting

  ! write velocity signal to file:
  sig_out = reshape( u_signal, (/ (Nx*Ny) /) )
  open (unit = 1,file = "./output/signal_mpi.txt")
  do j = 1,(Nx*Ny)
    write(1,*) sig_out(j)  
  end do
  close (1)

  ! write the analytical solution to file:
  sig_out = reshape( lap_signal, (/ (Nx*Ny) /) )
  open (unit = 1,file = "./output/lap_signal_mpi.txt")
  do j = 1,(Nx*Ny)
    write(1,*) sig_out(j)  
  end do
  close (1)

  ! write the computed solution to file:
  soln_out = reshape( lap_global , (/ (Nx*Ny) /) )
  open (unit = 2,file = "./output/lap_solution_mpi.txt")
  do j = 1,(Nx*Ny)
    write(2,*) soln_out(j)  
  end do
  close (2)
  
  ! 9) ADD BINARY FILE OUTPUT OPTION!

  end if ! rank == 0

end if ! if Nprocs = Np

call MPI_FINALIZE(ierror)

end program mpi_Laplacian2D
