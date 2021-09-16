program Zombie_Scenario_B
implicit none
!Model of a Zombien Outbrake
!
!dS/dt = -aSZ
!dZ/dt = aSZ 


integer i, itmax
!real*4 alpha_1,alpha_2,beta_1,beta_2,x0,y0
real*4 a,b,c,t0,S0,Z0,D0,R0
real*4 kf1,kf2,kf3,kf4,kg1,kg2,kg3,kg4,ki1,ki2,ki3,ki4,kh1,kh2,kh3,kh4,h,time
real*4 ,allocatable :: t(:),S(:),Z(:),D(:),R(:)
character*50 hchar
logical disease 

!!!!!!!!!!!!!!!!!!!!!!!!!!
t0=30.0
time = t0 + 100.0
h = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!
S0=87659   !Population at day 30
Z0=12350		  !Zombies at day 30
D0=0.0
R0=0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!
a=0.0000055
b=0.15
c=0.05
!!!!!!!!!!!!!!!!!!!!!!!!!!
itmax = INT((time-t0)/h)

!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(t(0:itmax),S(0:itmax),Z(0:itmax),D(0:itmax),R(0:itmax))

!!!!!!!!!!!!!!!!!!!!!!!!!!
t(0)=t0
S(0)=S0
Z(0)=Z0
D(0)=D0
R(0)=R0
!!!!!!!!!!!!!!!!!!!!!!!!!!
open(1,file="scenarioB.txt",status="replace")
!write(1,"(3f12,5)")0.0,x0,y
!write(1,"(3f12,5)")0.0,x0,y0
!write(1,"(3f12,5)")0.0,x0,y0
write(1,*)t0,S0,Z0,D0,R0
!!!!!!!!!!!!!!!!!!!!!!!!!!

disease = .TRUE.
DO i=0, itmax-1
	
	kf1 = h*dot_S( a, b, c, S(i), Z(i) )
	kg1 = h*dot_Z( a, b, c, S(i), Z(i) )
	kh1 = h*dot_D( a, b, c, S(i), Z(i) )
	ki1 = h*dot_R( a, b, c, S(i), Z(i) )

	kf2 = h*dot_S( a, b, c, S(i)+(kf1/2), Z(i)+(kg1/2) )
	kg2 = h*dot_Z( a, b, c, S(i)+(kf1/2), Z(i)+(kg1/2) )
	kh2 = h*dot_D( a, b, c, S(i)+(kf1/2), Z(i)+(kg1/2) )
	ki2 = h*dot_R( a, b, c, S(i)+(kf1/2), Z(i)+(kg1/2) )

	kf3 = h*dot_S( a, b, c, S(i)+(kf2/2), Z(i)+(kg2/2) )
	kg3 = h*dot_Z( a, b, c, S(i)+(kf2/2), Z(i)+(kg2/2) )
	kh3 = h*dot_D( a, b, c, S(i)+(kf2/2), Z(i)+(kg2/2) )
	ki3 = h*dot_R( a, b, c, S(i)+(kf2/2), Z(i)+(kg2/2) )


	kf4 = h*dot_S( a, b, c, S(i)+kf3, Z(i)+kg3 )
	kg4 = h*dot_Z( a, b, c, S(i)+kf3, Z(i)+kg3 )	
	kh4 = h*dot_D( a, b, c, S(i)+kf3, Z(i)+kg3 )
	ki4 = h*dot_R( a, b, c, S(i)+kf3, Z(i)+kg3 )

	S(i+1) = S(i) + (kf1 + (kf2*2) + (kf3*2) + kf4)/6.0
	Z(i+1) = Z(i) + (kg1 + (kg2*2) + (kg3*2) + kg4)/6.0
	D(i+1) = D(i) + (kh1 + (kh2*2) + (kh3*2) + kh4)/6.0
	R(i+1) = R(i) + (ki1 + (ki2*2) + (ki3*2) + ki4)/6.0

	write(1,*)t0+(float(i+1)*h),S(i+1),Z(i+1),D(i+1),R(i+1)

	IF ( (Z(i)< 1) .AND. (disease .EQV. .TRUE.)) THEN
	disease = .FALSE.
	write(*,*)"Zombie outbrake finish at day",int(t0+(float(i+1)*h)),"after first infection"
	write(*,*)"Deceased people:", int(D(i))
	write(*,*)"Recovered people :", int(R(i))
	write(*,*)"Never infected people :", int(S(i))
	ENDIF

END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
close(1)
!write(*,*)"Program finished succesfully"
!write(*,*)itmax

contains
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
function dot_S(ai,bi,ci,si,zi)
implicit none
real*4 ai,bi,ci,si,zi,dot_S
dot_S = -a*si*zi
return
end function
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
function dot_Z(ai,bi,ci,si,zi)
implicit none
real*4 ai,bi,ci,si,zi,dot_Z
dot_Z = a*si*zi - (bi*zi) - (ci*zi)
return
end function
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
function dot_D(ai,bi,ci,si,zi)
implicit none
real*4 ai,bi,ci,si,zi,dot_D
dot_D = (bi*zi) 
return
end function
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
function dot_R(ai,bi,ci,si,zi)
implicit none
real*4 ai,bi,ci,si,zi,dot_R
dot_R = (ci*zi) 
return
end function
!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!



end program


! contains
! !!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
! function f(sigma,xi,yi)
! implicit none
! real*4 sigma,xi,yi,f
! !LV_f = (a1*xi) - (b1*xi*yi)
! f=sigma*(yi-xi)
! return
! end function
! !!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
! function g(rho,xi,yi,zi)
! implicit none
! real*4 rho,xi,yi,zi,g
! !LV_g = (a2*xi*yi) - (b2*yi)
! g = (rho*xi) - yi - (xi*zi)
! return 
! end function
! !!!!!!!!!!!!&&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!
! function h(beta,xi,yi,zi)
! implicit none
! real*4 beta,xi,yi,zi,h
! !LV_g = (a2*xi*yi) - (b2*yi)
! h = (-beta*zi) + (xi*yi)
! return 
! end function