program Zombie_Scenario_2
implicit none
!The equation of the modell of a Zombien Outbrake
!Which can be expressed as
!dS/dt = -aSZ
!dZ/dt = aSZ 


integer i, itmax
!real*4 alpha_1,alpha_2,beta_1,beta_2,x0,y0
real*4 a,t0,S0,Z0
real*4 kf1,kf2,kf3,kf4,kg1,kg2,kg3,kg4,h,time
real*4 ,allocatable :: t(:),S(:),Z(:)
character*50 hchar

!!!!!!!!!!!!!!!!!!!!!!!!!!
t0=0.0
time = 80.0
h = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!
S0=100000.0   !Initial population all of the susceptible individuals
Z0=10.0		  !Number initial of zombies introduced	
!!!!!!!!!!!!!!!!!!!!!!!!!!
a=0.0000025
!!!!!!!!!!!!!!!!!!!!!!!!!!
itmax = INT((time-t0)/h)

!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(t(0:itmax),S(0:itmax),Z(0:itmax))

!!!!!!!!!!!!!!!!!!!!!!!!!!
t(0)=t0
S(0)=S0
Z(0)=Z0
!!!!!!!!!!!!!!!!!!!!!!!!!!
open(1,file="scenarioA.txt",status="replace")
!write(1,"(3f12,5)")0.0,x0,y
!write(1,"(3f12,5)")0.0,x0,y0
!write(1,"(3f12,5)")0.0,x0,y0
write(1,*)t0,S0,Z0
!!!!!!!!!!!!!!!!!!!!!!!!!!
DO i=0, itmax-1
	
	kf1 = h*dot_S( a, S(i), Z(i) )
	kg1 = h*dot_Z( a, S(i), Z(i) )

	kf2 = h*dot_S( a, S(i)+(kf1/2), Z(i)+(kg1/2) )
	kg2 = h*dot_Z( a, S(i)+(kf1/2), Z(i)+(kg1/2) )

	kf3 = h*dot_S( a, S(i)+(kf2/2), Z(i)+(kg2/2) )
	kg3 = h*dot_Z( a, S(i)+(kf2/2), Z(i)+(kg2/2) )

	kf4 = h*dot_S( a, S(i)+kf3, Z(i)+kg3 )
	kg4 = h*dot_Z( a, S(i)+kf3, Z(i)+kg3 )

	S(i+1) = S(i) + (kf1 + (kf2*2) + (kf3*2) + kf4)/6.0
	Z(i+1) = Z(i) + (kg1 + (kg2*2) + (kg3*2) + kg4)/6.0

	write(1,*)float(i+1)*h,S(i+1),Z(i+1)

IF(float(i+1)*h==30.0)THEN
	write(*,*)"Susceptible people (S) at day 30:"
	write(*,*)int(S(i))
	write(*,*)"Population of zombies (Z) at day 30:"
	write(*,*)int(Z(i))
ENDIF

END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
close(1)
write(*,*)"Program finished succesfully"

contains
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
function dot_S(ai,si,zi)
implicit none
real*4 ai,si,zi,dot_S
dot_S = -a*si*zi
return
end function
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
function dot_Z(ai,si,zi)
implicit none
real*4 ai,si,zi,dot_Z
dot_Z = a*si*zi
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







