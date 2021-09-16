program lotka_volterra
implicit none
!The equations to solve are
!dx/dt = alpha_1*x + beta_1*x*y = f(x,y)
!dx/dt = alpha_2*x*y - beta_2*y = g(x,y)
!Where x is the prey polpulation and y is the predator population

integer i, itmax
real*4 alpha_1,alpha_2,beta_1,beta_2,x0,y0
real*4 kf1,kf2,kf3,kf4,kg1,kg2,kg3,kg4,h,time
real*4 ,allocatable :: x(:),y(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!
time = 10.0
h = 1.e-2

!!!!!!!!!!!!!!!!!!!!!!!!!!
x0=2.0
y0=1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!
alpha_1 = 2
alpha_2 = 0.9
beta_1 = 1.2
beta_2 = 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!
itmax = INT(time/h)

!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(x(0:itmax),y(0:itmax))

!!!!!!!!!!!!!!!!!!!!!!!!!!
x(0)=x0
y(0)=y0

!!!!!!!!!!!!!!!!!!!!!!!!!!
open(1,file="results.txt",status="replace")
!write(1,"(3f12,5)")0.0,x0,y
!write(1,"(3f12,5)")0.0,x0,y0
!write(1,"(3f12,5)")0.0,x0,y0
write(1,*)0.0,x0,y0
!!!!!!!!!!!!!!!!!!!!!!!!!!
DO i=0, itmax-1
	
	kf1 = h*LV_f( alpha_1, beta_1, x(i), y(i) )
	kg1 = h*LV_g( alpha_2, beta_2, x(i), y(i) )

	kf2 = h*LV_f( alpha_1, beta_1, x(i)+(kf1/2), y(i)+(kg1/2) )
	kg2 = h*LV_g( alpha_2, beta_2, x(i)+(kf1/2), y(i)+(kg1/2) )

	kf3 = h*LV_f( alpha_1, beta_1, x(i)+(kf2/2), y(i)+(kg2/2) )
	kg3 = h*LV_g( alpha_2, beta_2, x(i)+(kf2/2), y(i)+(kg2/2) )

	kf4 = h*LV_f( alpha_1, beta_1, x(i)+kf3, y(i)+kg3 )
	kg4 = h*LV_g( alpha_2, beta_2, x(i)+kf3, y(i)+kg3 )

	x(i+1) = x(i) + (kf1 + (kf2*2) + (kf3*2) + kf4)/6.0
	y(i+1) = y(i) + (kg1 + (kg2*2) + (kg3*2) + kg4)/6.0

	write(1,*)float(i+1)*h,x(i+1),y(i+1)

END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
close(1)
write(*,*)"Program finished succesfully"

contains
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
function LV_f(a1,b1,xi,yi)
implicit none
real*4 a1,b1,xi,yi,LV_f
LV_f = (a1*xi) - (b1*xi*yi)
return
end function
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
function LV_g(a2,b2,xi,yi)
implicit none
real*4 a2,b2,xi,yi,LV_g
LV_g = (a2*xi*yi) - (b2*yi)
return
end function
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
end program








