program SDE	!This program solves a 3x3 SDE
implicit none
!The equations to solve are
!dx/dt = sigma(y-x) = f(x,y)
!dx/dt = rho*x - y - x*z  = g(x,y)
!dz/dt = -beta*x + x*y

integer i, itmax
real*4 beta,rho,sigma,x0,y0,z0
real*4 kf1,kf2,kf3,kf4,kg1,kg2,kg3,kg4,kh1,kh2,kh3,kh4,dt,time
real*4 ,allocatable :: x(:),y(:),z(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!
time = 45.0
dt   = 1.e-2

!!!!!!!!!!!!!!!!!!!!!!!!!!
x0 = 0.0
y0 = 1.0
z0 = 0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!
beta  = (8/3)
sigma = 10.0
rho   = 28.0
!!!!!!!!!!!!!!!!!!!!!!!!!!
itmax = INT(time/dt)

!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(x(0:itmax),y(0:itmax),z(0:itmax))

!!!!!!!!!!!!!!!!!!!!!!!!!!
x(0)=x0
y(0)=y0
z(0)=z0

!!!!!!!!!!!!!!!!!!!!!!!!!!
open(1,file="results.txt",status="replace")
!write(1,"(4f12,5)")0.0,x0,y
!write(1,"(4f12,5)")0.0,x0,y0
!write(1,"(4f12,5)")0.0,x0,y0
write(1,*)0.0,x0,y0,z0
!!!!!!!!!!!!!!!!!!!!!!!!!!
DO i=0, itmax-1
	
	kf1 = dt*f( sigma, x(i), y(i) )
	kg1 = dt*g( rho  , x(i), y(i) , z(i))
	kh1 = dt*h(beta  , x(i), y(i) , z(i) )

	kf2 = dt*f( sigma, x(i)+(kf1/2) , y(i)+(kg1/2) )
	kg2 = dt*g( rho  , x(i)+(kf1/2) , y(i)+(kg1/2) , z(i)+(kh1/2) )
	kh2 = dt*h( beta , x(i)+(kf1/2) , y(i)+(kg1/2) , z(i)+(kh1/2))

	kf3 = dt*f( sigma, x(i)+(kf2/2), y(i)+(kg2/2) )
	kg3 = dt*g( rho  , x(i)+(kf2/2), y(i)+(kg2/2) , z(i)+(kh2/2))
	kh3 = dt*h( beta , x(i)+(kf2/2), y(i)+(kg2/2) , z(i)+(kh2/2) )
	
	kf4 = dt*f(  sigma, x(i)+kf3, y(i)+kg3 )
	kg4 = dt*g(  rho  , x(i)+kf3, y(i)+kg3 , z(i)+(kh3/2))
	kh4 = dt*h(  beta , x(i)+kf3, y(i)+kg3 , z(i)+(kh3/2))

	x(i+1) = x(i) + (kf1 + (kf2*2) + (kf3*2) + kf4)/6.0
	y(i+1) = y(i) + (kg1 + (kg2*2) + (kg3*2) + kg4)/6.0
	z(i+1) = z(i) + (kh1 + (kh2*2) + (kh3*2) + kh4)/6.0

	write(1,*)float(i+1)*dt,x(i+1),y(i+1),z(i+1)

END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
close(1)
write(*,*)"Program finished succesfully"
write(*,*)itmax

contains
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
function f(sigma,xi,yi)
implicit none
real*4 sigma,xi,yi,f
!LV_f = (a1*xi) - (b1*xi*yi)
f=sigma*(yi-xi)
return
end function
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
function g(rho,xi,yi,zi)
implicit none
real*4 rho,xi,yi,zi,g
!LV_g = (a2*xi*yi) - (b2*yi)
g = (rho*xi) - yi - (xi*zi)
return 
end function
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
function h(beta,xi,yi,zi)
implicit none
real*4 beta,xi,yi,zi,h
!LV_g = (a2*xi*yi) - (b2*yi)
h = (-beta*zi) + (xi*yi)
return 
end function
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
end program








