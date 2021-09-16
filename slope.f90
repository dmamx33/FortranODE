PROGRAM slope_field				!INICIO DEL PROGRAMA <<FORTRAN NO DISTINGUE ENTRE MAYUS - MIN>>
implicit none					!ALL VARIABLES MUST BE DECLARED

!VARIABLE DECLARATION
integer nx, ny, i, j
real*4	x0, x1, y0, y1, L, xi, yj, dxi, dyj, th, m, dx,dy		!real*4 -> stablishes single precision
real*4, allocatable :: x(:), y(:)								!Store selected xi and yj points

!INITIALIZATION OF VARIABLES
nx=50 ; ny=50
x0=-3.0 ; x1=3.0 !Domain
y0=-3.0 ; y1=3.0 !Domain
dx = (x1 -x0) / float(nx)
dy = (y1 - y0) / float(ny)
L = 0.6 * sqrt( dx**2 + dy**2 )

!DYNAMIC MEMORY ALLOCATION
allocate( x(0:nx) , y(0:ny) )

!GRID GENERATION
do i=0,nx
	x(i) = x0 + (float(i)*dx)
end do

do i=0,ny
	y(i) = y0 + (float(i)*dy)
end do

!Open a file to save the results
Open(1,file="results.txt",status="replace")


!CALCULATE OF THE SLOPE FIELD
do i=0, nx
	do j=0, ny
		!Calculation of the slope
		m  = y(j) - x(i)
		th = atan(m)
		!Calculation of the sides of the triangle
		dxi = L*cos(th)
		dyj = L*sin(th)
		!Calculation of the position of the initial point
		xi = x(i) - dxi/2.0
		yj = y(j) - dyj/2.0
		!Stores the results in the file created
		write(1,*)xi,yj,dxi,dyj
	end do
end do

!Close the file
Close(1)

end PROGRAM

