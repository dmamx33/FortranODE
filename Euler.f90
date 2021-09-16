program euler_method !THIS METHOD IS UESD TO SOLVE AN ODE USING EULER METHOD
implicit none

!Variable declaration
integer i, itmax
real*4 h,time,y0,t0
real*4,allocatable :: t(:),y(:)

!VARIABLE INITIALIZATION
h = 0.1
time = 5.0
t0 = 0.0
y0 = 1
itmax = int(time/h)

!DYNAMIC MEMORY ALLOCATION
allocate(t(0:itmax),y(0:itmax))

!INITIAL CONDITIONS
t(0) = t0
y(0) = y0

!NUMERICAL SOLUTION
do i = 0, itmax-1
	t(i+1) = t(i) + h
	y(i+1) = y(i) + h * f1(t(i),y(i))
end do

!WRITE RESULTS INTO A FILE
open(1,file="euler.txt",status="replace")
do i = 0, itmax
	write(1,*)t(i),y(i),0.25*exp(-2.*t(i))*(t(i)**4+4)
end do
close(1)
write(*,*) "Program finished succesfully"

!SECTION TO WRITE FUNCTIONS FOR THE PROGRAM
contains	

function f1(ti,yi) !!!!!!!!!!!!!!!!!!!!!!!!
	
	implicit none
	real*4 ti,yi,f1

		f1 = (ti**3.0*exp(-2.0*ti))-(2.0*yi)
	
	return	
end function!!!!!!!!!!!!!!!!!!!!!!!!

end program