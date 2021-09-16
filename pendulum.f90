program pendulum
implicit none
!The equation of the modell of a pendulum
!dQ2/d2t = -g/L sinQ - (c/mL) dQ/dt 
!Which can be expressed as
!dw/dt = (-g/L)sinQ -(c/mL)w
!dQ/dt = w
!!!!!!!!!!!!!!!!!!dx/dt = alpha_1*x + beta_1*x*y = f(x,y)
!!!!!!!!!!!!!!!!!!dx/dt = alpha_2*x*y - beta_2*y = g(x,y)
!!!!!!!!!!!!!!!!!!Where x is the prey polpulation and y is the predator population

integer i, itmax
!real*4 alpha_1,alpha_2,beta_1,beta_2,x0,y0
real*4 m,g,L,c,t0,theta0,w0,multiplo,ce
real*4 kf1,kf2,kf3,kf4,kg1,kg2,kg3,kg4,h,time,Pi
real*4 ,allocatable :: t(:),theta(:),w(:)
character*50 hchar

!!!!!!!!!!!!!!!!!!!!!!!!!!
time = 40.0
h = 1.e-3

Pi=3.1415926
g=9.81

!!!!!!!!!!!!!!!!!!!!!!!!!!
t0=0.0
theta0=0.0
multiplo = 3.0 !!!3.5*Pi --  4.0*Pi  ---  5.0*Pi
w0=multiplo*Pi

!!!!!!!!!!!!!!!!!!!!!!!!!!
ce=0.05
L=0.5
m = 0.2
!!!!!!!!!!!!!!!!!!!!!!!!!!
itmax = INT((time-t0)/h)

!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(t(0:itmax),w(0:itmax),theta(0:itmax))

!!!!!!!!!!!!!!!!!!!!!!!!!!
t(0)=t0
w(0)=w0
theta(0)=theta0
!!!!!!!!!!!!!!!!!!!!!!!!!!
	write(hchar,"(f6.4)")multiplo
	hchar = adjustl(hchar)
open(1,file="pendulum_wo_"//hchar(1:len_trim(hchar))//".txt",status="replace")
!write(1,"(3f12,5)")0.0,x0,y
!write(1,"(3f12,5)")0.0,x0,y0
!write(1,"(3f12,5)")0.0,x0,y0
write(1,*)t0,w0,theta0
!!!!!!!!!!!!!!!!!!!!!!!!!!
DO i=0, itmax-1
	
	kf1 = 	  h*dot_w( ce,m,L,g, w(i), theta(i) )
	kg1 = h*dot_tetha( ce,m,L,g, w(i), theta(i) )

	kf2 =     h*dot_w( ce,m,L,g, w(i)+(kf1/2), theta(i)+(kg1/2) )
	kg2 = h*dot_tetha( ce,m,L,g, w(i)+(kf1/2), theta(i)+(kg1/2) )

	kf3 =     h*dot_w( ce,m,L,g, w(i)+(kf2/2), theta(i)+(kg2/2) )
	kg3 = h*dot_tetha( ce,m,L,g, w(i)+(kf2/2), theta(i)+(kg2/2) )

	kf4 =     h*dot_w( ce,m,L,g, w(i)+kf3, theta(i)+kg3 )
	kg4 = h*dot_tetha( ce,m,L,g, w(i)+kf3, theta(i)+kg3 )

	w(i+1) = w(i) + (kf1 + (kf2*2) + (kf3*2) + kf4)/6.0
	theta(i+1) = theta(i) + (kg1 + (kg2*2) + (kg3*2) + kg4)/6.0

	write(1,*)float(i+1)*h,w(i+1),theta(i+1)

END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
close(1)
write(*,*)"Program finished succesfully"

contains
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
function dot_w(c1,m1,L1,g1,wi,thetai)
implicit none
real*4 c1,m1,L1,g1,wi,thetai,dot_w
dot_w = -(g1/L1)*sin(thetai)-((c1/(m1*L1)*wi))
return
end function
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
function dot_tetha(c1,m1,L1,g1,wi,thetai)
implicit none
real*4 c1,m1,L1,g1,wi,thetai,dot_tetha
dot_tetha = wi
return
end function
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
end program

! !!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
! function LV_f(a1,b1,xi,yi)
! implicit none
! real*4 a1,b1,xi,yi,LV_f
! LV_f = (a1*xi) - (b1*xi*yi)
! return
! end function
! !!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
! function LV_g(a2,b2,xi,yi)
! implicit none
! real*4 a2,b2,xi,yi,LV_g
! LV_g = (a2*xi*yi) - (b2*yi)
! return
! end function
! !!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!







