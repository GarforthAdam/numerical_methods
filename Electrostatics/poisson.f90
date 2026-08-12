! Used to solve the Poisson problem for a 2D 10 cm X 10 cm
! grounding conducting box that has electrical potential U = 0 everywhere on the boundary.
! At the cell (25,50), the boxes potential appears to converge to a value found to be
! about-3.18351102

program poisson
implicit none             ! Turn off implicit typing
                                         !Variable Dictionary  
integer, parameter :: NX=100              ! Set the number of cells in x-direction
integer, parameter :: NY=100              ! Set the number of cells in y-direction
real, parameter :: H=1                  ! Cell size
real, dimension(0:NY+1,0:NX+1) :: u       ! Define array holding old estimate of the potential
real, dimension(0:NY+1,0:NX+1) :: unew    ! Define array holding new estimate of the potential
real, dimension(0:NY+1,0:NX+1) :: q       ! Charge density
real :: max_change=1.0                    ! Set value for max change
integer :: niter=0                        ! Set number of iterations to zero
integer :: i,j                            ! Define variables for loop indices
real, parameter :: PI=3.141592654         ! Pi
real :: x, y                              ! Define variables for cell control Coordinates
integer :: lun1                           ! LUN for I/O

u = 0.0 ! Initialize u to zero everywhere
unew = u ! Initialize unew to u
q = 0.0 ! First initialize q to zero everywhere, same as u
q(25,25) = -4.0 ! Set charge density in cell (25,25)
q(75,75) = 4.0 ! Set charge density in cell (75,75)

iter_loop: do while(max_change > 1.0e-5)   

   do i=1,NY ! Loop over interior of box with nested do loops
      do j=1,NX
         unew(i,j) = (u(i+1,j)+u(i-1,j)+u(i,j+1)+ & !Calculate given equation
         u(i,j-1)+4.0*PI*(H**2)*q(i,j))/4.0

      enddo
   enddo

! Use maxval intrinsic function to find the maximum
! absolute value of change over the entire set of cells
   max_change = maxval(abs(u-unew)) ! Calculate maximum change value

   u = unew ! Make old guess the new guess

   niter = niter+1 ! Count number of iterations
enddo iter_loop

write(*,*) ' It took ',niter,' iterations to converge'
write(*,*) 'potential at (25,50):', u(25,50)  !Write value for potential at (25,50) in the box

!Write out unew in Gnuplot form
open(newunit=lun1,file='poisson.dat',status='replace') !Set to replace other file called poisson.dat
write(*,*) 'New file successfully created'

x = 0.5*H ! Initialize y coordinate of cell
do j=0,NX+1 ! Loop over columns

   y = 0.5*H ! Initialize x coordinate of cell
   do i=0,NY+1 ! Loop over each row
      write(LUN1,*) x,y,unew(i,j)
      y = y+h ! Increment x coordinate of cell
   enddo
   write(lun1,*) ' ' ! Write blank line after row
   x = x+h ! Increment y coordinate to next row
enddo
close(unit=lun1) ! Close new file

stop 0
end program poisson
