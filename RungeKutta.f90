program euler_method !THIS METHOD IS UESD TO SOLVE AN ODE USING EULER METHOD
implicit none

!Variable declaration
character*50 hchar
integer i, itmax
real*4 h,time,y0,t0,k1,k2,k3,k4
real*4,allocatable :: t(:),y(:)

!VARIABLE INITIALIZATION
h = 0.02
time = 5.0
	!convert h from real to h string
	write(hchar,"(f5.3)")h
	hchar = adjustl(hchar)
t0 = 0.0
y0 = 1
itmax = int(time-t0/h)

!DYNAMIC MEMORY ALLOCATION
allocate(t(0:itmax),y(0:itmax))

!INITIAL CONDITIONS
t(0) = t0
y(0) = y0

!NUMERICAL SOLUTION
do i = 0, itmax-1
	t(i+1) = t(i) + h
	k1 = h*f1(t(i),y(i))
	k2 = h*f1(t(i)+h/2,y(i)+k1/2)
	k3 = h*f1(t(i)+h/2,y(i)+k2/2)
	k4 = h*f1(t(i)+h,y(i)+k3)
	!y(i+1) = y(i) + h * f1(t(i),y(i))
	y(i+1) = y(i) + ( k1 + 2*k2 + 2*k3 +k4 )/6
end do

!WRITE RESULTS INTO A FILE
open(1,file="RungeKutta_h_"//hchar(1:len_trim(hchar))//".txt",status="replace")
do i = 0, itmax
	write(1,*)t(i),y(i),0.25*exp(-2.*t(i))*(t(i)**4+4)!Analitical solution
end do
close(1)

!SECTION TO WRITE FUNCTIOS FOR THE PROGRAM
contains	

function f1(ti,yi) !!!!!!!!!!!!!!!!!!!!!!!!
	
	implicit none
	real*4 ti,yi,f1

		f1 = (ti**3.0*exp(-2.0*ti))-(2.0*yi)!differential equation
	
	return	
end function!!!!!!!!!!!!!!!!!!!!!!!!

end program