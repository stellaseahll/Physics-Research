Program fb3D
Implicit None
Type Particles
   DoublePrecision:: x, xold, y, yold, z, zold, VelX, VelXold, VelY, VelYold, VelZ, VelZold, Accelx, Accelxold, Accely, Accelyold, Accelz, Accelzold, &
   rotationx, rotationxold, rotationy, rotationyold, rotationz, rotationzold, AngularVelx, AngularVelxold, AngularVely, AngularVelyold, AngularVelz, AngularVelzold, &
   AngularAccelx, AngularAccelxold, AngularAccely, AngularAccelyold, AngularAccelz, AngularAccelzold, &
   Fx, BaseFx, Wall1Fx, Wall2Fx, ParticleFx, Fy, BaseFy, Wall1Fy, Wall2Fy, ParticleFy, &
   Fz, BaseFz, Wall1Fz, Wall2Fz, ParticleFz, Momentx, Momenty, Momentz, BaseMomentx, BaseMomenty, BaseMomentz, &
   Wall1Momentx, Wall1Momenty, Wall1Momentz, Wall2Momentx, Wall2Momenty, Wall2Momentz, &
   ParticleMomentx, ParticleMomenty, ParticleMomentz, &
   ContactFx, BaseContactFx, Wall1ContactFx, Wall2ContactFx, ParticleContactFx, &
   ContactFy, BaseContactFy, Wall1ContactFy, Wall2ContactFy, ParticleContactFy, &
   ContactFz, BaseContactFz, Wall1ContactFz, Wall2ContactFz, ParticleContactFz, &
   ContactMomentx, ContactMomenty, ContactMomentz, BaseContactMomentx, BaseContactMomenty, BaseContactMomentz, &
   Wall1ContactMomentx, Wall1ContactMomenty, Wall1ContactMomentz, &
   Wall2ContactMomentx, Wall2ContactMomenty, Wall2ContactMomentz, &
   ParticleContactMomentx, ParticleContactMomenty, ParticleContactMomentz, &
   group, ParticleF, drag, & 
   Fcap, FcapX, FcapY, FcapZ, & !FvX, FvY, FvZ, FddX, FddY, FddZ, FdlvoX, FdlvoY, FdlvoZ
   m, MomentInert, cnormp, cshearp, cnormw, cshearw
   Integer:: CellNum, R, C, L
   Logical:: BaseCollision, Wall1Collision, Wall2Collision, ParticleCollision
End Type
Integer, Parameter:: N=56000, RowTotal=240, ColTotal=30, LayTotal=12, MaxParticle=20, MaxIterate=1
Integer:: Row, Col, Lay, RowTemp, RowTemp2, ColTemp, ColTemp2, LayTemp, LayTemp2, CellNo, i, j, k, iterate, baserow !, Gs
Integer, Dimension(0:RowTotal+1,0:ColTotal+1,0:LayTotal+1,0:MaxParticle+1):: Cell
DoublePrecision, Parameter:: TotalTime=10.0, Diameter=0.15, Length=1.2, Depth=0.06, t=1.0e-6, &
timeframe=0.01, CFDtime=1.0e-6, d=3.0e-3, soliddensity1=1400, soliddensity2=2600, g=9.81, kn=1000.0, ks=1000.0, & 
phi=0.3, cohesion=0.0e-8, cnorm1=0.004, cshear1=0.004, cnorm2=0.006, cshear2=0.006, roll=5.0e-5, & 
vin=1.8, angle=90.0, gas_start=0.5, density=1.205, viscosity=1.8e-5, relax1=0.8, relax2=0.2, &
wetpercent=0.0, surften=0.073, contact=0.0, &
!Ha=1.0e-23, Boltz=1.381e-23, Temperature=298.0, &
!mu=1.2566e-6, magnetized=8.0e3, magneticfield=8.0e5, & !4.46e5
!lambdaB=0.7e-9, Zcharge=1.0, ions=1.0e3
ConvScheme=0.0, & !1.0 Central Differencing Scheme, 0.0 Upwind Differencing Scheme
TimeScheme=1.0 !1.0 Three Time Levels Scheme, 0.0 Euler Implicit Scheme
DoublePrecision:: a, b, c, CollisionTime, NormalOverlap, TangentOverlap1, TangentOverlap2, TangentOverlap3, &
NormalForce, ShearForce1, ShearForce2, ShearForce3, MaxShear, NormalDamp, ShearDamp1, ShearDamp2, ShearDamp3, &
Fn, Fs1, Fs2, Fs3, MomentInert, vx, vy, vz, e1, e2, e3, t1, t2, t3, t4, t5, t6, ndot, sdot1, sdot2, sdot3, &
xnew, ynew, znew, rotationxnew, rotationynew, rotationznew, RollFriction1, RollFriction2, RollFriction3, &
Velxnew, Velynew, Velznew, AngularVelxnew, AngularVelynew, AngularVelznew, &
ContactRoll1, ContactRoll2, ContactRoll3, time, pi, num, base, baseVel, vs, xpos, ypos, zpos, E, &
velYtotal, velXtotal, velZtotal, velYavg, velXavg, velZavg, vYtotal, vXtotal, vZtotal, vYavg, vXavg, vZavg, count, Pavg, &
avgcount, sample, group1, sigma2, sigma02, sigmaR2, mu, Lacey, h, &
!h, Fv, random1, random2, random3, FbX, FbY, FbZ, magnetization, dipolemoment, kappa, Fdlvo
Fcap, Vcap, App, Bpp, Cpp, Apw, Bpw, Cpw, hcpp, hcpw, &
dx, dy, dz, CE, CN, CP, FluxUUDS, FluxVUDS, FluxUCDS, FluxVCDS, &
PN, PS, PE, PW, APT, AWC, dPxMid, dPyMid, dPxE, dPyN, uMid, vMid, APUMid, APVMid, &
uE, vN, MassSum, PrN, PrS, PrE, PrW, &
ReynoldsX, ReynoldsY, ReynoldsZ, dragcoeffX, dragcoeffY, dragcoeffZ, chiX, chiY, chiZ, fxo, fyo, fzo, dragX, dragY, dragZ, Cd, Cdprime
DoublePrecision, Dimension(0:LayTotal+1):: MassFlow, OutFlow
DoublePrecision, Dimension(0:ColTotal+1,0:RowTotal+1,0:LayTotal+1):: u, uo, uoo, v, vo, voo, P, &
FlowX, FlowY, SU, SV, APU, APV, AN, AS, AE, AW, AP, dPx, dPy, Pr, &
porosity, DragForceX, DragForceY
DoublePrecision, Dimension(0:RowTotal):: factorX, factorY
logical:: overlap
Type(Particles), Dimension(1:N):: Particle
Open(1,File='fbseg1.8e.95.txt')
!Open(2,File='fb3d1.4vel1.txt')
!Open(3,File='fb3d1.4gas1.txt')
!Open(4,File='fb3d1.4vs1.txt')
!Open(5,File='fb3d1.4Gs.txt')
Open(6,File='fbseg1.8e.95Pr.txt')
!Open(7,File='fb3d1.4vel2.txt')
!Open(8,File='fb3d1.4gas2.txt')
!Open(9,File='fb3d1.4vs2.txt')
!Open(10,File='fb3d1.4T.txt')
Open(11,File='fbseg1.8e.95gvec.txt')
Open(12,File='fbseg1.8e.95F.txt')
Open(13,File='fbseg1.8e.95mix.txt')
Open(14,File='fbseg1.8e.95Lacey.txt')
time=0
Cell=0
pi=ATAN(1.0)*4.0
!m=4.0/3.0*pi*(d/2.0)**3*soliddensity
!MomentInert=0.4*m*(d/2.0)**2
!magnetization=0.0 !(N*4.0/3.0*pi*(d/2.0)**3)/(Length*Diameter*d)*magnetized
!dipolemoment=(mu*magnetization*pi*d**3)/6.0
!kappa=SQRT(4.0*pi*lambdaB*ions)
Vcap=4.0/3.0*pi*(d/2.0)**3*0.5*(soliddensity1+soliddensity2)*wetpercent/(100.0*1000.0) !m*wetpercent/(100.0*1000.0)
App=-1.1*(Vcap/(pi*(d/2.0)**3))**-0.53
Bpp=(-0.34*LOG(Vcap/(pi*(d/2.0)**3))-0.96)*contact**2-0.019*LOG(Vcap/(pi*(d/2.0)**3))+0.48
Cpp=0.0042*LOG(Vcap/(pi*(d/2.0)**3))+0.078
hcpp=(d/2.0)*((0.62*contact+0.99)*(Vcap/(pi*(d/2.0)**3))**0.34)
Apw=-1.9*(Vcap/(pi*(d/2.0)**3))**-0.51
Bpw=(-0.016*LOG(Vcap/(pi*(d/2.0)**3))-0.76)*contact**2-0.12*LOG(Vcap/(pi*(d/2.0)**3))+1.2
Cpw=0.013*LOG(Vcap/(pi*(d/2.0)**3))+0.18
hcpw=(d/2.0)*((0.22*contact+0.95)*(Vcap/(pi*(d/2.0)**3))**0.32)
dx=Diameter/ColTotal
dy=Length/RowTotal
dz=Depth/LayTotal
Call Random_Seed()
!Gs=0
!DEM Initialization
i=1
xpos=d
ypos=d+0.01
zpos=d
Do While (i<=N)
   Do While (i<=N .AND. zpos<Depth-d)
   Do While (i<=N .AND. xpos<Diameter-d)
   !overlap=.True.
   !Do While (overlap)
   !   Call Random_Number(xpos)
   !   xpos=xpos*Diameter*0.9+d
   !   Call Random_Number(ypos)
   !   ypos=ypos*Length*0.8+5.0e-3
   !   Call Random_Number(zpos)
   !   zpos=zpos*Depth*0.9+d !!!
   !  Particle(i)%x=xpos
   !  Particle(i)%y=ypos
   !  Particle(i)%z=zpos
   !  overlap=.False.
   !  If (i>1) then
   !     Do j=1,i-1
   !	    If ((Particle(i)%x-Particle(j)%x)**2+(Particle(i)%y-Particle(j)%y)**2+(Particle(i)%z-Particle(j)%z)**2<(d+1.0e-9)**2) overlap=.True.
   !     EndDo
   !  EndIf
   !EndDo
   Row=INT(ypos/(Length/RowTotal))+1
   Col=INT(xpos/(Diameter/ColTotal))+1
   Lay=INT(zpos/(Depth/LayTotal))+1
   k=1
   Do While (Cell(Row,Col,Lay,k)/=0)
      k=k+1
   EndDo
   Cell(Row,Col,Lay,k)=i
   Particle(i)%VelX=0.01
   Particle(i)%VelY=0.0
   Particle(i)%VelZ=0.01
   Particle(i)%VelXold=0.0
   Particle(i)%VelYold=0.0
   Particle(i)%VelZold=0.0
   Particle(i)%Accelx=0.0
   Particle(i)%Accely=0.0
   Particle(i)%Accelz=0.0
   Particle(i)%Accelxold=0.0
   Particle(i)%Accelyold=0.0
   Particle(i)%Accelzold=0.0
   Particle(i)%AngularVelx=0.0
   Particle(i)%AngularVely=0.0
   Particle(i)%AngularVelz=0.0
   Particle(i)%AngularVelxold=0.0
   Particle(i)%AngularVelyold=0.0
   Particle(i)%AngularVelzold=0.0
   Particle(i)%AngularAccelx=0.0
   Particle(i)%AngularAccely=0.0
   Particle(i)%AngularAccelz=0.0
   Particle(i)%AngularAccelxold=0.0
   Particle(i)%AngularAccelyold=0.0
   Particle(i)%AngularAccelzold=0.0
   Particle(i)%x=xpos
   Particle(i)%xold=Particle(i)%x-Particle(i)%VelX*t
   Particle(i)%y=ypos
   Particle(i)%yold=Particle(i)%y-Particle(i)%VelY*t
   Particle(i)%z=zpos
   Particle(i)%zold=Particle(i)%z-Particle(i)%VelZ*t
   Particle(i)%rotationx=0.0
   Particle(i)%rotationxold=Particle(i)%rotationx-Particle(i)%AngularVelx*t
   Particle(i)%rotationy=0.0
   Particle(i)%rotationyold=Particle(i)%rotationy-Particle(i)%AngularVely*t
   Particle(i)%rotationz=0.0
   Particle(i)%rotationzold=Particle(i)%rotationz-Particle(i)%AngularVelz*t
   Particle(i)%Fx=0.0
   Particle(i)%BaseFx=0.0
   Particle(i)%Wall1Fx=0.0
   Particle(i)%Wall2Fx=0.0
   Particle(i)%ParticleFx=0.0
   Particle(i)%Fy=0.0
   Particle(i)%BaseFy=0.0
   Particle(i)%Wall1Fy=0.0
   Particle(i)%Wall2Fy=0.0
   Particle(i)%ParticleFy=0.0
   Particle(i)%Fz=0.0
   Particle(i)%BaseFz=0.0
   Particle(i)%Wall1Fz=0.0
   Particle(i)%Wall2Fz=0.0
   Particle(i)%ParticleFz=0.0
   Particle(i)%Momentx=0.0
   Particle(i)%BaseMomentx=0.0
   Particle(i)%Wall1Momentx=0.0
   Particle(i)%Wall2Momentx=0.0
   Particle(i)%ParticleMomentx=0.0
   Particle(i)%Momenty=0.0
   Particle(i)%BaseMomenty=0.0
   Particle(i)%Wall1Momenty=0.0
   Particle(i)%Wall2Momenty=0.0
   Particle(i)%ParticleMomenty=0.0
   Particle(i)%Momentz=0.0
   Particle(i)%BaseMomentz=0.0
   Particle(i)%Wall1Momentz=0.0
   Particle(i)%Wall2Momentz=0.0
   Particle(i)%ParticleMomentz=0.0
   Particle(i)%ContactFx=0.0
   Particle(i)%Wall1ContactFx=0.0
   Particle(i)%Wall2ContactFx=0.0
   Particle(i)%BaseContactFx=0.0
   Particle(i)%ParticleContactFx=0.0
   Particle(i)%ContactFy=0.0
   Particle(i)%Wall1ContactFy=0.0
   Particle(i)%Wall2ContactFy=0.0
   Particle(i)%BaseContactFy=0.0
   Particle(i)%ParticleContactFy=0.0
   Particle(i)%ContactFz=0.0
   Particle(i)%Wall1ContactFz=0.0
   Particle(i)%Wall2ContactFz=0.0
   Particle(i)%BaseContactFz=0.0
   Particle(i)%ParticleContactFz=0.0
   Particle(i)%ContactMomentx=0.0
   Particle(i)%BaseContactMomentx=0.0
   Particle(i)%Wall1ContactMomentx=0.0
   Particle(i)%Wall2ContactMomentx=0.0
   Particle(i)%ParticleContactMomentx=0.0
   Particle(i)%ContactMomenty=0.0
   Particle(i)%BaseContactMomenty=0.0
   Particle(i)%Wall1ContactMomenty=0.0
   Particle(i)%Wall2ContactMomenty=0.0
   Particle(i)%ParticleContactMomenty=0.0
   Particle(i)%ContactMomentz=0.0
   Particle(i)%BaseContactMomentz=0.0
   Particle(i)%Wall1ContactMomentz=0.0
   Particle(i)%Wall2ContactMomentz=0.0
   Particle(i)%ParticleContactMomentz=0.0
   Particle(i)%CellNum=(Lay-1)*RowTotal*ColTotal+(Row-1)*ColTotal+Col
   Particle(i)%R=Row
   Particle(i)%C=Col
   Particle(i)%L=Lay
   Particle(i)%BaseCollision=.False.
   Particle(i)%Wall1Collision=.False.
   Particle(i)%Wall2Collision=.False.
   Particle(i)%ParticleCollision=.False.
   Particle(i)%group=0.0
   Particle(i)%ParticleF=0.0
   Particle(i)%drag=0.0
!   Particle(i)%FvX=0.0
!   Particle(i)%FvY=0.0
!   Particle(i)%FvZ=0.0
!   Particle(i)%FddX=0.0
!   Particle(i)%FddY=0.0
!   Particle(i)%FddZ=0.0
!   Particle(i)%FdlvoX=0.0
!   Particle(i)%FdlvoY=0.0
!   Particle(i)%FdlvoZ=0.0
    Particle(i)%FcapX=0.0
	Particle(i)%FcapY=0.0
	Particle(i)%FcapZ=0.0
	Particle(i)%Fcap=0.0
Particle(i)%m=4.0/3.0*pi*(d/2.0)**3*soliddensity1
Particle(i)%MomentInert=0.4*Particle(i)%m*(d/2.0)**2
Particle(i)%cnormp=cnorm1
Particle(i)%cshearp=cshear1
Particle(i)%cnormw=cnorm1
Particle(i)%cshearw=cshear1
Particle(i)%group=1.0
   i=i+1
   xpos=xpos+d
   EndDo
   xpos=d
   zpos=zpos+d
   EndDo
   xpos=d
   zpos=d
   ypos=ypos+d
EndDo
Do i=10002,50000,2
   Particle(i)%m=4.0/3.0*pi*(d/2.0)**3*soliddensity2
   Particle(i)%MomentInert=0.4*Particle(i)%m*(d/2.0)**2
   Particle(i)%cnormp=cnorm2
   Particle(i)%cshearp=cshear2
   Particle(i)%cnormw=cnorm2
   Particle(i)%cshearw=cshear2
   Particle(i)%group=2.0
EndDo
!CFD Initialization
FlowX=0.0
FlowY=0.0
u=0.0
v=0.0
P=1.0e5
uo=0.0
vo=0.0
AN=0.0
AS=0.0
AE=0.0
AW=0.0
AP=0.0
Pr=0.0
factorX=0.5
factorY=0.5
factorX(0)=0.0
factorY(0)=0.0
factorX(ColTotal)=1.0
factorY(RowTotal)=1.0
porosity=1.0
mu=0.5
sigma02=0.25
Do While (time<=TotalTime)
   Do i=1,N
	  Particle(i)%BaseCollision=.False.
	  Particle(i)%Wall1Collision=.False.
	  Particle(i)%Wall2Collision=.False.
	  Particle(i)%ParticleCollision=.False.
   EndDo
   MassFlow=0.0
   Do Lay=1,LayTotal
   Do Col=1,ColTotal
	  If (time<gas_start) then
	     v(Col,0,Lay)=0.0
		 !group1=0.0
		 !Do i=1,N
			!If (Particle(i)%y>0.0 .AND. Particle(i)%y<=0.044) then
			   !Particle(i)%group=1.0
			   !group1=group1+1.0
			!Else
			   !Particle(i)%group=2.0
			!EndIf
		 !EndDo
		 !mu=group1/N
		 !sigma02=(group1/N)*(1.0-group1/N)
	  Else
         v(Col,0,Lay)=vin
	  EndIf
	  FlowY(Col,0,Lay)=density*dx*v(Col,0,Lay) !*porosity(Col,0)
      MassFlow(Lay)=MassFlow(Lay)+FlowY(Col,0,Lay)
   EndDo
   EndDo
   uoo=uo
   voo=vo
   uo=u
   vo=v
   porosity=1.0
   DragForceX=0.0
   DragForceY=0.0
   base=0.0
   baserow=1
   baseVel=0.0
   !$OMP parallel Do
   Do i=1,N !DEM
      Do Row=(Particle(i)%R-1),(Particle(i)%R+1)
	     Do Col=(Particle(i)%C-1),(Particle(i)%C+1)
		    Do Lay=(Particle(i)%L-1),(Particle(i)%L+1)
		    !If (Row>RowTotal) Exit
	        If (Row<baserow .AND. Particle(i)%C==Col .AND. Particle(i)%L==Lay .AND. Particle(i)%y-d/2.0<base) then 
			!Base Collision
			   Particle(i)%BaseCollision=.True.
			   CollisionTime=(d/2.0+base-Particle(i)%y)/(ABS(baseVel-Particle(i)%VelY)+1.0e-30)
			   If (CollisionTime>t) CollisionTime=t
			   NormalOverlap=(Particle(i)%VelY-baseVel)*CollisionTime
		       TangentOverlap1=(Particle(i)%VelX+Particle(i)%AngularVelz*d/2.0)*CollisionTime
			   TangentOverlap2=(Particle(i)%VelZ+Particle(i)%AngularVelx*d/2.0)*CollisionTime
		       NormalForce=-NormalOverlap*kn
		       ShearForce1=-TangentOverlap1*ks
			   ShearForce2=-TangentOverlap2*ks
		       MaxShear=NormalForce*TAN(phi)+cohesion !Coulomb Friction Law
		       NormalDamp=-(Particle(i)%VelY-baseVel)*Particle(i)%cnormw
		       ShearDamp1=-(Particle(i)%VelX+Particle(i)%AngularVelz*d/2.0)*Particle(i)%cshearw
			   ShearDamp2=-(Particle(i)%VelZ+Particle(i)%AngularVelx*d/2.0)*Particle(i)%cshearw
		       Fn=NormalForce+NormalDamp
		       Fs1=ShearForce1+ShearDamp1
			   Fs2=ShearForce2+ShearDamp2
			   If (ABS(NormalForce*roll)<ABS(ShearForce1*d/2.0)) then
			      ContactRoll1=ShearForce1/(ABS(ShearForce1)+1.0e-30)*ABS(NormalForce*roll)
			   Else
			      ContactRoll1=ShearForce1*d/2.0
			   EndIf
			   If (ABS(NormalForce*roll)<ABS(ShearForce2*d/2.0)) then
			      ContactRoll2=ShearForce2/(ABS(ShearForce2)+1.0e-30)*ABS(NormalForce*roll)
			   Else
			      ContactRoll2=ShearForce2*d/2.0
			   EndIf
			   If (ABS(Fs1)>ABS(MaxShear)) then
			      Fs1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(MaxShear)
				  ShearForce1=0.0
				  ContactRoll1=0.0
			   EndIf
			   If (ABS(Fs2)>ABS(MaxShear)) then
			      Fs2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(MaxShear)
				  ShearForce2=0.0
				  ContactRoll2=0.0
			   EndIf
			   If (ABS(Fn*roll)>ABS(Fs1*d/2.0)) then
			      RollFriction1=Fs1*d/2.0
			   Else
			      RollFriction1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(Fn*roll) !Rolling Friction
			   EndIf
			   If (ABS(Fn*roll)>ABS(Fs2*d/2.0)) then
			      RollFriction2=Fs2*d/2.0
			   Else
			      RollFriction2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(Fn*roll) !Rolling Friction
			   EndIf
		       Particle(i)%BaseFx=Particle(i)%BaseContactFx+Fs1
			   Particle(i)%BaseContactFx=Particle(i)%BaseContactFx+ShearForce1
			   Particle(i)%BaseFy=Particle(i)%BaseContactFy+Fn
			   Particle(i)%BaseContactFy=Particle(i)%BaseContactFy+NormalForce
			   Particle(i)%BaseFz=Particle(i)%BaseContactFz+Fs2
			   Particle(i)%BaseContactFz=Particle(i)%BaseContactFz+ShearForce2
		       Particle(i)%BaseMomentz=Particle(i)%BaseContactMomentz+Fs1*d/2.0-RollFriction1
			   Particle(i)%BaseContactMomentz=Particle(i)%BaseContactMomentz+ShearForce1*d/2.0-ContactRoll1
			   Particle(i)%BaseMomentx=Particle(i)%BaseContactMomentx+Fs2*d/2.0-RollFriction2
			   Particle(i)%BaseContactMomentx=Particle(i)%BaseContactMomentx+ShearForce2*d/2.0-ContactRoll2
			ElseIf (Row>RowTotal .AND. Particle(i)%C==Col .AND. Particle(i)%L==Lay .AND. Particle(i)%y>Length-d/2.0) then 
			!Top Collision
			   Particle(i)%BaseCollision=.True.
			   CollisionTime=(Particle(i)%y-Length+d/2.0)/(ABS(Particle(i)%VelY)+1.0e-30)
			   If (CollisionTime>t) CollisionTime=t
			   NormalOverlap=Particle(i)%VelY*CollisionTime
		       TangentOverlap1=(Particle(i)%VelX-Particle(i)%AngularVelz*d/2.0)*CollisionTime
			   TangentOverlap2=(Particle(i)%VelZ-Particle(i)%AngularVelx*d/2.0)*CollisionTime
		       NormalForce=-NormalOverlap*kn
		       ShearForce1=-TangentOverlap1*ks
			   ShearForce2=-TangentOverlap2*ks
		       MaxShear=NormalForce*TAN(phi)+cohesion
		       NormalDamp=-Particle(i)%VelY*Particle(i)%cnormw
		       ShearDamp1=-(Particle(i)%VelX-Particle(i)%AngularVelz*d/2.0)*Particle(i)%cshearw
			   ShearDamp2=-(Particle(i)%VelZ-Particle(i)%AngularVelx*d/2.0)*Particle(i)%cshearw
		       Fn=NormalForce+NormalDamp
		       Fs1=ShearForce1+ShearDamp1
			   Fs2=ShearForce2+ShearDamp2
			   If (ABS(NormalForce*roll)<ABS(ShearForce1*d/2.0)) then
			      ContactRoll1=ShearForce1/(ABS(ShearForce1)+1.0e-30)*ABS(NormalForce*roll)
			   Else
			      ContactRoll1=ShearForce1*d/2.0
			   EndIf
			   If (ABS(NormalForce*roll)<ABS(ShearForce2*d/2.0)) then
			      ContactRoll2=ShearForce2/(ABS(ShearForce2)+1.0e-30)*ABS(NormalForce*roll)
			   Else
			      ContactRoll2=ShearForce2*d/2.0
			   EndIf
			   If (ABS(Fs1)>ABS(MaxShear)) then
			      Fs1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(MaxShear)
			      ShearForce1=0.0
			      ContactRoll1=0.0
			   EndIf
			   If (ABS(Fs2)>ABS(MaxShear)) then
			      Fs2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(MaxShear)
			      ShearForce2=0.0
			      ContactRoll2=0.0
			   EndIf
			   If (ABS(Fn*roll)>ABS(Fs1*d/2.0)) then
			      RollFriction1=Fs1*d/2.0
			   Else
			      RollFriction1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(Fn*roll)
			   EndIf
			   If (ABS(Fn*roll)>ABS(Fs2*d/2.0)) then
			      RollFriction2=Fs2*d/2.0
			   Else
			      RollFriction2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(Fn*roll)
			   EndIf
		       Particle(i)%BaseFx=Particle(i)%BaseContactFx+Fs1
			   Particle(i)%BaseContactFx=Particle(i)%BaseContactFx+ShearForce1
			   Particle(i)%BaseFy=Particle(i)%BaseContactFy+Fn
			   Particle(i)%BaseContactFy=Particle(i)%BaseContactFy+NormalForce
			   Particle(i)%BaseFz=Particle(i)%BaseContactFz+Fs2
			   Particle(i)%BaseContactFz=Particle(i)%BaseContactFz+ShearForce2
		       Particle(i)%BaseMomentz=Particle(i)%BaseContactMomentz-Fs1*d/2.0+RollFriction1
			   Particle(i)%BaseContactMomentz=Particle(i)%BaseContactMomentz-ShearForce1*d/2.0+ContactRoll1
			   Particle(i)%BaseMomentx=Particle(i)%BaseContactMomentx-Fs2*d/2.0+RollFriction2
			   Particle(i)%BaseContactMomentx=Particle(i)%BaseContactMomentx-ShearForce2*d/2.0+ContactRoll2
	        EndIf
		    If (Col<1 .AND. Particle(i)%R==Row .AND. Particle(i)%L==Lay .AND. Particle(i)%x<d/2.0) then
			!Left Wall Collision
			   Particle(i)%Wall1Collision=.True.
		       CollisionTime=(d/2.0-Particle(i)%x)/(ABS(Particle(i)%VelX)+1.0e-30)
			   If (CollisionTime>t) CollisionTime=t
		       NormalOverlap=Particle(i)%VelX*CollisionTime
		       TangentOverlap1=(Particle(i)%VelY-Particle(i)%AngularVelz*d/2.0)*CollisionTime
			   TangentOverlap2=(Particle(i)%VelZ-Particle(i)%AngularVely*d/2.0)*CollisionTime
		       NormalForce=-NormalOverlap*kn
		       ShearForce1=-TangentOverlap1*ks
			   ShearForce2=-TangentOverlap2*ks
		       MaxShear=NormalForce*TAN(phi)+cohesion 
		       NormalDamp=-Particle(i)%VelX*Particle(i)%cnormw
		       ShearDamp1=-(Particle(i)%VelY-Particle(i)%AngularVelz*d/2.0)*Particle(i)%cshearw
			   ShearDamp2=-(Particle(i)%VelZ-Particle(i)%AngularVely*d/2.0)*Particle(i)%cshearw
		       Fn=NormalForce+NormalDamp
		       Fs1=ShearForce1+ShearDamp1
			   Fs2=ShearForce2+ShearDamp2
			   If (ABS(NormalForce*roll)<ABS(ShearForce1*d/2.0)) then
			      ContactRoll1=ShearForce1/(ABS(ShearForce1)+1.0e-30)*ABS(NormalForce*roll)
			   Else
			      ContactRoll1=ShearForce1*d/2.0
			   EndIf
			   If (ABS(NormalForce*roll)<ABS(ShearForce2*d/2.0)) then
			      ContactRoll2=ShearForce2/(ABS(ShearForce2)+1.0e-30)*ABS(NormalForce*roll)
			   Else
			      ContactRoll2=ShearForce2*d/2.0
			   EndIf
			   If (ABS(Fs1)>ABS(MaxShear)) then
			      Fs1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(MaxShear)
				  ShearForce1=0.0
				  ContactRoll1=0.0
			   EndIf
			   If (ABS(Fs2)>ABS(MaxShear)) then
			      Fs2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(MaxShear)
				  ShearForce2=0.0
				  ContactRoll2=0.0
			   EndIf
			   If (ABS(Fn*roll)>ABS(Fs1*d/2.0)) then
			      RollFriction1=Fs1*d/2.0
			   Else
			      RollFriction1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(Fn*roll)
			   EndIf
			   If (ABS(Fn*roll)>ABS(Fs2*d/2.0)) then
			      RollFriction2=Fs2*d/2.0
			   Else
			      RollFriction2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(Fn*roll)
			   EndIf
		       Particle(i)%Wall1Fx=Particle(i)%Wall1ContactFx+Fn
			   Particle(i)%Wall1ContactFx=Particle(i)%Wall1ContactFx+NormalForce
		       Particle(i)%Wall1Fy=Particle(i)%Wall1ContactFy+Fs1
			   Particle(i)%Wall1ContactFy=Particle(i)%Wall1ContactFy+ShearForce1
			   Particle(i)%Wall1Fz=Particle(i)%Wall1ContactFz+Fs2
			   Particle(i)%Wall1ContactFz=Particle(i)%Wall1ContactFz+ShearForce2
		       Particle(i)%Wall1Momentz=Particle(i)%Wall1ContactMomentz-Fs1*d/2.0+RollFriction1
			   Particle(i)%Wall1ContactMomentz=Particle(i)%Wall1ContactMomentz-ShearForce1*d/2.0+ContactRoll1
			   Particle(i)%Wall1Momenty=Particle(i)%Wall1ContactMomenty-Fs2*d/2.0+RollFriction2
			   Particle(i)%Wall1ContactMomenty=Particle(i)%Wall1ContactMomenty-ShearForce2*d/2.0+ContactRoll2
		    ElseIf (Col>ColTotal .AND. Particle(i)%R==Row .AND. Particle(i)%L==Lay .AND. Particle(i)%x>Diameter-d/2.0) then
			!Right Wall Collision
			   Particle(i)%Wall1Collision=.True.
			   CollisionTime=(d/2.0-Diameter+Particle(i)%x)/(ABS(Particle(i)%VelX)+1.0e-30)
			   If (CollisionTime>t) CollisionTime=t
			   NormalOverlap=Particle(i)%VelX*CollisionTime
		       TangentOverlap1=(Particle(i)%VelY+Particle(i)%AngularVelz*d/2.0)*CollisionTime
			   TangentOverlap2=(Particle(i)%VelZ+Particle(i)%AngularVely*d/2.0)*CollisionTime
		       NormalForce=-NormalOverlap*kn
		       ShearForce1=-TangentOverlap1*ks
			   ShearForce2=-TangentOverlap2*ks
		       MaxShear=NormalForce*TAN(phi)+cohesion 
		       NormalDamp=-Particle(i)%VelX*Particle(i)%cnormw
			   ShearDamp1=-(Particle(i)%VelY+Particle(i)%AngularVelz*d/2.0)*Particle(i)%cshearw
			   ShearDamp2=-(Particle(i)%VelZ+Particle(i)%AngularVely*d/2.0)*Particle(i)%cshearw
		       Fn=NormalForce+NormalDamp
		       Fs1=ShearForce1+ShearDamp1
			   Fs2=ShearForce2+ShearDamp2
  			   If (ABS(NormalForce*roll)<ABS(ShearForce1*d/2.0)) then
			      ContactRoll1=ShearForce1/(ABS(ShearForce1)+1.0e-30)*ABS(NormalForce*roll)
			   Else
			      ContactRoll1=ShearForce1*d/2.0
			   EndIf
			   If (ABS(NormalForce*roll)<ABS(ShearForce2*d/2.0)) then
			      ContactRoll2=ShearForce2/(ABS(ShearForce2)+1.0e-30)*ABS(NormalForce*roll)
			   Else
			      ContactRoll2=ShearForce2*d/2.0
			   EndIf
			   If (ABS(Fs1)>ABS(MaxShear)) then
			      Fs1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(MaxShear)
				  ShearForce1=0.0
				  ContactRoll1=0.0
			   EndIf
			   If (ABS(Fs2)>ABS(MaxShear)) then
			      Fs2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(MaxShear)
				  ShearForce2=0.0
				  ContactRoll2=0.0
			   EndIf
			   If (ABS(Fn*roll)>ABS(Fs1*d/2.0)) then
			      RollFriction1=Fs1*d/2.0
			   Else
			      RollFriction1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(Fn*roll)
			   EndIf
			   If (ABS(Fn*roll)>ABS(Fs2*d/2.0)) then
			      RollFriction2=Fs2*d/2.0
			   Else
			      RollFriction2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(Fn*roll)
			   EndIf
		       Particle(i)%Wall1Fx=Particle(i)%Wall1ContactFx+Fn
			   Particle(i)%Wall1ContactFx=Particle(i)%Wall1ContactFx+NormalForce
		       Particle(i)%Wall1Fy=Particle(i)%Wall1ContactFy+Fs1
			   Particle(i)%Wall1ContactFy=Particle(i)%Wall1ContactFy+ShearForce1
			   Particle(i)%Wall1Fz=Particle(i)%Wall1ContactFz+Fs2
			   Particle(i)%Wall1ContactFz=Particle(i)%Wall1ContactFz+ShearForce2
			   Particle(i)%Wall1Momentz=Particle(i)%Wall1ContactMomentz+Fs1*d/2.0-RollFriction1
			   Particle(i)%Wall1ContactMomentz=Particle(i)%Wall1ContactMomentz+ShearForce1*d/2.0-ContactRoll1
			   Particle(i)%Wall1Momenty=Particle(i)%Wall1ContactMomenty+Fs2*d/2.0-RollFriction2
			   Particle(i)%Wall1ContactMomenty=Particle(i)%Wall1ContactMomenty+ShearForce2*d/2.0-ContactRoll2
		    EndIf
			If (Lay<1 .AND. Particle(i)%R==Row .AND. Particle(i)%C==Col .AND. Particle(i)%z<d/2.0) then
			!Front Wall Collision
			   Particle(i)%Wall2Collision=.True.
		       CollisionTime=(d/2.0-Particle(i)%z)/(ABS(Particle(i)%VelZ)+1.0e-30)
			   If (CollisionTime>t) CollisionTime=t
		       NormalOverlap=Particle(i)%VelZ*CollisionTime
		       TangentOverlap1=(Particle(i)%VelY-Particle(i)%AngularVelx*d/2.0)*CollisionTime
			   TangentOverlap2=(Particle(i)%VelX+Particle(i)%AngularVely*d/2.0)*CollisionTime
		       NormalForce=-NormalOverlap*kn
		       ShearForce1=-TangentOverlap1*ks
			   ShearForce2=-TangentOverlap2*ks
		       MaxShear=NormalForce*TAN(phi)+cohesion 
		       NormalDamp=-Particle(i)%VelZ*Particle(i)%cnormw
		       ShearDamp1=-(Particle(i)%VelY-Particle(i)%AngularVelx*d/2.0)*Particle(i)%cshearw
			   ShearDamp2=-(Particle(i)%VelX+Particle(i)%AngularVely*d/2.0)*Particle(i)%cshearw
		       Fn=NormalForce+NormalDamp
		       Fs1=ShearForce1+ShearDamp1
			   Fs2=ShearForce2+ShearDamp2
			   If (ABS(NormalForce*roll)<ABS(ShearForce1*d/2.0)) then
			      ContactRoll1=ShearForce1/(ABS(ShearForce1)+1.0e-30)*ABS(NormalForce*roll)
			   Else
			      ContactRoll1=ShearForce1*d/2.0
			   EndIf
			   If (ABS(NormalForce*roll)<ABS(ShearForce2*d/2.0)) then
			      ContactRoll2=ShearForce1/(ABS(ShearForce2)+1.0e-30)*ABS(NormalForce*roll)
			   Else
			      ContactRoll2=ShearForce2*d/2.0
			   EndIf
			   If (ABS(Fs1)>ABS(MaxShear)) then
			      Fs1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(MaxShear)
				  ShearForce1=0.0
				  ContactRoll1=0.0
			   EndIf
			   If (ABS(Fs2)>ABS(MaxShear)) then
			      Fs2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(MaxShear)
				  ShearForce2=0.0
				  ContactRoll2=0.0
			   EndIf
			   If (ABS(Fn*roll)>ABS(Fs1*d/2.0)) then
			      RollFriction1=Fs1*d/2.0
			   Else
			      RollFriction1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(Fn*roll)
			   EndIf
			   If (ABS(Fn*roll)>ABS(Fs2*d/2.0)) then
			      RollFriction2=Fs2*d/2.0
			   Else
			      RollFriction2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(Fn*roll)
			   EndIf
		       Particle(i)%Wall2Fz=Particle(i)%Wall2ContactFz+Fn
			   Particle(i)%Wall2ContactFz=Particle(i)%Wall2ContactFz+NormalForce
		       Particle(i)%Wall2Fy=Particle(i)%Wall2ContactFy+Fs1
			   Particle(i)%Wall2ContactFy=Particle(i)%Wall2ContactFy+ShearForce1
			   Particle(i)%Wall2Fx=Particle(i)%Wall2ContactFx+Fs2
			   Particle(i)%Wall2ContactFx=Particle(i)%Wall2ContactFx+ShearForce2
		       Particle(i)%Wall2Momentx=Particle(i)%Wall2ContactMomentx-Fs1*d/2.0+RollFriction1
			   Particle(i)%Wall2ContactMomentx=Particle(i)%Wall2ContactMomentx-ShearForce1*d/2.0+ContactRoll1
			   Particle(i)%Wall2Momenty=Particle(i)%Wall2ContactMomenty+Fs2*d/2.0-RollFriction2
			   Particle(i)%Wall2ContactMomenty=Particle(i)%Wall2ContactMomenty+ShearForce2*d/2.0-ContactRoll2
		    ElseIf (Lay>LayTotal .AND. Particle(i)%R==Row .AND. Particle(i)%C==Col .AND. Particle(i)%z>Depth-d/2.0) then
			!Rear Wall Collision
			   Particle(i)%Wall2Collision=.True.
			   CollisionTime=(d/2.0-Depth+Particle(i)%z)/(ABS(Particle(i)%VelZ)+1.0e-30)
			   If (CollisionTime>t) CollisionTime=t
			   NormalOverlap=Particle(i)%VelZ*CollisionTime
		       TangentOverlap1=(Particle(i)%VelY+Particle(i)%AngularVelx*d/2.0)*CollisionTime
			   TangentOverlap2=(Particle(i)%VelX-Particle(i)%AngularVely*d/2.0)*CollisionTime
		       NormalForce=-NormalOverlap*kn
		       ShearForce1=-TangentOverlap1*ks
			   ShearForce2=-TangentOverlap2*ks
		       MaxShear=NormalForce*TAN(phi)+cohesion 
		       NormalDamp=-Particle(i)%VelZ*Particle(i)%cnormw
			   ShearDamp1=-(Particle(i)%VelY+Particle(i)%AngularVelx*d/2.0)*Particle(i)%cshearw
			   ShearDamp2=-(Particle(i)%VelX-Particle(i)%AngularVely*d/2.0)*Particle(i)%cshearw
		       Fn=NormalForce+NormalDamp
		       Fs1=ShearForce1+ShearDamp1
			   Fs2=ShearForce2+ShearDamp2
  			   If (ABS(NormalForce*roll)<ABS(ShearForce1*d/2.0)) then
			      ContactRoll1=ShearForce1/(ABS(ShearForce1)+1.0e-30)*ABS(NormalForce*roll)
			   Else
			      ContactRoll1=ShearForce1*d/2.0
			   EndIf
			   If (ABS(NormalForce*roll)<ABS(ShearForce2*d/2.0)) then
			      ContactRoll2=ShearForce2/(ABS(ShearForce2)+1.0e-30)*ABS(NormalForce*roll)
			   Else
			      ContactRoll2=ShearForce2*d/2.0
			   EndIf
			   If (ABS(Fs1)>ABS(MaxShear)) then
			      Fs1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(MaxShear)
				  ShearForce1=0.0
				  ContactRoll1=0.0
			   EndIf
			   If (ABS(Fs2)>ABS(MaxShear)) then
			      Fs2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(MaxShear)
				  ShearForce2=0.0
				  ContactRoll2=0.0
			   EndIf
			   If (ABS(Fn*roll)>ABS(Fs1*d/2.0)) then
			      RollFriction1=Fs1*d/2.0
			   Else
			      RollFriction1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(Fn*roll)
			   EndIf
			   If (ABS(Fn*roll)>ABS(Fs2*d/2.0)) then
			      RollFriction2=Fs2*d/2.0
			   Else
			      RollFriction2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(Fn*roll)
			   EndIf
		       Particle(i)%Wall2Fz=Particle(i)%Wall2ContactFz+Fn
			   Particle(i)%Wall2ContactFz=Particle(i)%Wall2ContactFz+NormalForce
		       Particle(i)%Wall2Fy=Particle(i)%Wall2ContactFy+Fs1
			   Particle(i)%Wall2ContactFy=Particle(i)%Wall2ContactFy+ShearForce1
			   Particle(i)%Wall2Fx=Particle(i)%Wall2ContactFx+Fs2
			   Particle(i)%Wall2ContactFx=Particle(i)%Wall2ContactFx+ShearForce2
			   Particle(i)%Wall2Momentx=Particle(i)%Wall2ContactMomentx+Fs1*d/2.0-RollFriction1
			   Particle(i)%Wall2ContactMomentx=Particle(i)%Wall2ContactMomentx+ShearForce1*d/2.0-ContactRoll1
			   Particle(i)%Wall2Momenty=Particle(i)%Wall2ContactMomenty-Fs2*d/2.0+RollFriction2
			   Particle(i)%Wall2ContactMomenty=Particle(i)%Wall2ContactMomenty-ShearForce2*d/2.0+ContactRoll2
		    EndIf
		    k=1
		    Do While (k<=MaxParticle)
			   If (Cell(Row,Col,Lay,k)==0 .OR. Cell(Row,Col,Lay,k)==i .OR. Cell(Row,Col,Lay,k)<=i) then
			      k=k+1
			   Else
			      j=Cell(Row,Col,Lay,k)
				  If ((Particle(i)%x-Particle(j)%x)**2+(Particle(i)%y-Particle(j)%y)**2+(Particle(i)%z-Particle(j)%z)**2<d**2) then
				  !Particle-Particle Collision
			   		 Particle(i)%ParticleCollision=.True.
					 Particle(j)%ParticleCollision=.True.
					 vx=Particle(j)%VelX-Particle(i)%VelX
				     vy=Particle(j)%VelY-Particle(i)%VelY
					 vz=Particle(j)%VelZ-Particle(i)%VelZ
				     a=vx**2+vy**2+vz**2+1.0e-30
				     b=-2.0*(vx*(Particle(j)%x-Particle(i)%x)+vy*(Particle(j)%y-Particle(i)%y)+vz*(Particle(j)%z-Particle(i)%z))
				     c=(Particle(j)%x-Particle(i)%x)**2+(Particle(j)%y-Particle(i)%y)**2+(Particle(j)%z-Particle(i)%z)**2-d**2
				     CollisionTime=(-b+SQRT(b**2-4.0*a*c))/(2.0*a)
					 If (CollisionTime>t) CollisionTime=t
				     e1=(Particle(j)%x-Particle(i)%x)/d
				     e2=(Particle(j)%y-Particle(i)%y)/d
					 e3=(Particle(j)%z-Particle(i)%z)/d
				     ndot=(Particle(i)%VelX-Particle(j)%VelX)*e1+(Particle(i)%VelY-Particle(j)%VelY)*e2+(Particle(i)%VelZ-Particle(j)%VelZ)*e3	  
                     t1=e2
				     t2=-e1
				     sdot1=(Particle(i)%VelX-Particle(j)%VelX)*t1+(Particle(i)%VelY-Particle(j)%VelY)*t2- &
				     (Particle(i)%AngularVelz+Particle(j)%AngularVelz)*d/2.0 !xy-plane
					 t3=e3
				     t4=-e1
					 sdot2=(Particle(i)%VelX-Particle(j)%VelX)*t3+(Particle(i)%VelZ-Particle(j)%VelZ)*t4- &
				     (Particle(i)%AngularVely+Particle(j)%AngularVely)*d/2.0 !xz-plane
					 t5=e3
				     t6=-e2
					 sdot3=(Particle(i)%VelY-Particle(j)%VelY)*t5+(Particle(i)%VelZ-Particle(j)%VelZ)*t6- &
				     (Particle(i)%AngularVelx+Particle(j)%AngularVelx)*d/2.0 !yz-plane
				     NormalOverlap=ndot*CollisionTime
				     TangentOverlap1=sdot1*CollisionTime
					 TangentOverlap2=sdot2*CollisionTime
					 TangentOverlap3=sdot3*CollisionTime
     				 NormalForce=NormalOverlap*kn
				     ShearForce1=TangentOverlap1*ks
					 ShearForce2=TangentOverlap2*ks
					 ShearForce3=TangentOverlap3*ks
					 MaxShear=NormalForce*TAN(phi)+cohesion
				     NormalDamp=ndot*Particle(i)%cnormp
		             ShearDamp1=sdot1*Particle(i)%cshearp
					 ShearDamp2=sdot2*Particle(i)%cshearp
					 ShearDamp3=sdot3*Particle(i)%cshearp
					 Fn=NormalForce+NormalDamp
		             Fs1=ShearForce1+ShearDamp1
					 Fs2=ShearForce2+ShearDamp2
					 Fs3=ShearForce3+ShearDamp3
					 If (ABS(NormalForce*roll)<ABS(ShearForce1*d/2.0)) then
			            ContactRoll1=ShearForce1/(ABS(ShearForce1)+1.0e-30)*ABS(NormalForce*roll)
			         Else
			            ContactRoll1=ShearForce1*d/2.0
			         EndIf
					 If (ABS(NormalForce*roll)<ABS(ShearForce2*d/2.0)) then
			            ContactRoll2=ShearForce2/(ABS(ShearForce2)+1.0e-30)*ABS(NormalForce*roll)
			         Else
			            ContactRoll2=ShearForce2*d/2.0
			         EndIf
					 If (ABS(NormalForce*roll)<ABS(ShearForce3*d/2.0)) then
			            ContactRoll3=ShearForce3/(ABS(ShearForce3)+1.0e-30)*ABS(NormalForce*roll)
			         Else
			            ContactRoll3=ShearForce3*d/2.0
			         EndIf
					 If (ABS(Fs1)>ABS(MaxShear)) then
					    Fs1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(MaxShear)
						ShearForce1=0.0
						ContactRoll1=0.0
					 EndIf
					 If (ABS(Fs2)>ABS(MaxShear)) then
					    Fs2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(MaxShear)
						ShearForce2=0.0
						ContactRoll2=0.0
					 EndIf
					 If (ABS(Fs3)>ABS(MaxShear)) then
					    Fs3=Fs3/(ABS(Fs3)+1.0e-30)*ABS(MaxShear)
						ShearForce3=0.0
						ContactRoll3=0.0
					 EndIf
					 If (ABS(Fn*roll)>ABS(Fs1*d/2.0)) then
			            RollFriction1=Fs1*d/2.0
			         Else
			            RollFriction1=Fs1/(ABS(Fs1)+1.0e-30)*ABS(Fn*roll)
			         EndIf
					 If (ABS(Fn*roll)>ABS(Fs2*d/2.0)) then
			            RollFriction2=Fs2*d/2.0
			         Else
			            RollFriction2=Fs2/(ABS(Fs2)+1.0e-30)*ABS(Fn*roll)
			         EndIf
					 If (ABS(Fn*roll)>ABS(Fs3*d/2.0)) then
			            RollFriction3=Fs3*d/2.0
			         Else
			            RollFriction3=Fs3/(ABS(Fs3)+1.0e-30)*ABS(Fn*roll)
			         EndIf
				     Particle(i)%ParticleFx=Particle(i)%ParticleContactFx-Fn*e1-Fs1*t1-Fs2*t3
					 Particle(i)%ParticleContactFx=Particle(i)%ParticleContactFx-NormalForce*e1-ShearForce1*t1-ShearForce2*t3
		             Particle(i)%ParticleFy=Particle(i)%ParticleContactFy-Fn*e2-Fs1*t2-Fs3*t5
					 Particle(i)%ParticleContactFy=Particle(i)%ParticleContactFy-NormalForce*e2-ShearForce1*t2-ShearForce3*t5
					 Particle(i)%ParticleFz=Particle(i)%ParticleContactFz-Fn*e3-Fs2*t4-Fs3*t6
					 Particle(i)%ParticleContactFz=Particle(i)%ParticleContactFz-NormalForce*e3-ShearForce2*t4-ShearForce3*t6
		             
					 Particle(i)%ParticleMomentx=Particle(i)%ParticleContactMomentx+Fs3*d/2.0-RollFriction3
					 Particle(i)%ParticleContactMomentx=Particle(i)%ParticleContactMomentx+ShearForce3*d/2.0-ContactRoll3
					 Particle(i)%ParticleMomenty=Particle(i)%ParticleContactMomenty+Fs2*d/2.0-RollFriction2
					 Particle(i)%ParticleContactMomenty=Particle(i)%ParticleContactMomenty+ShearForce2*d/2.0-ContactRoll2
					 Particle(i)%ParticleMomentz=Particle(i)%ParticleContactMomentz+Fs1*d/2.0-RollFriction1
					 Particle(i)%ParticleContactMomentz=Particle(i)%ParticleContactMomentz+ShearForce1*d/2.0-ContactRoll1
					 
					 Particle(j)%ParticleFx=Particle(j)%ParticleContactFx+Fn*e1+Fs1*t1+Fs2*t3
					 Particle(j)%ParticleContactFx=Particle(j)%ParticleContactFx+NormalForce*e1+ShearForce1*t1+ShearForce2*t3
		             Particle(j)%ParticleFy=Particle(j)%ParticleContactFy+Fn*e2+Fs1*t2+Fs3*t5
					 Particle(j)%ParticleContactFy=Particle(j)%ParticleContactFy+NormalForce*e2+ShearForce1*t2+ShearForce3*t5
					 Particle(j)%ParticleFz=Particle(j)%ParticleContactFz+Fn*e3+Fs2*t4+Fs3*t6
					 Particle(j)%ParticleContactFz=Particle(j)%ParticleContactFz+NormalForce*e3+ShearForce2*t4+ShearForce3*t6
				     
					 Particle(j)%ParticleMomentx=Particle(j)%ParticleContactMomentx+Fs3*d/2.0-RollFriction3
					 Particle(j)%ParticleContactMomentx=Particle(j)%ParticleContactMomentx+ShearForce3*d/2.0-ContactRoll3
					 Particle(j)%ParticleMomenty=Particle(j)%ParticleContactMomenty+Fs2*d/2.0-RollFriction2
					 Particle(j)%ParticleContactMomenty=Particle(j)%ParticleContactMomenty+ShearForce2*d/2.0-ContactRoll2
					 Particle(j)%ParticleMomentz=Particle(j)%ParticleContactMomentz+Fs1*d/2.0-RollFriction1
					 Particle(j)%ParticleContactMomentz=Particle(j)%ParticleContactMomentz+ShearForce1*d/2.0-ContactRoll1

					 Particle(i)%ParticleF=SQRT(Particle(i)%ParticleFx**2+Particle(i)%ParticleFy**2+Particle(i)%ParticleFz**2)
					 Particle(j)%ParticleF=SQRT(Particle(j)%ParticleFx**2+Particle(j)%ParticleFy**2+Particle(j)%ParticleFz**2)
			      EndIf
				  k=k+1
			   EndIf
		    EndDo
	     EndDo
	  EndDo
	EndDo
	  !Drag Force
	  dragX=0.0
	  dragY=0.0
      num=0
      Do Row=(Particle(i)%R-1),(Particle(i)%R+1)
	     Do Col=(Particle(i)%C-1),(Particle(i)%C+1)
		    Do Lay=(Particle(i)%L-1),(Particle(i)%L+1)
	  !Row=Particle(i)%R
	  !Col=Particle(i)%C
	           Do k=1,MaxParticle
                  If (Cell(Row,Col,Lay,k)/=0) num=num+1
			   EndDo
		    EndDo
		 EndDo
      EndDo
	  !porosity(Particle(i)%C,Particle(i)%R)=1.0-num*m/(soliddensity*dx*dy*d)
      If (Particle(i)%C==1 .OR. Particle(i)%C==ColTotal .OR. Particle(i)%R==1 .OR. Particle(i)%R==RowTotal .OR. Particle(i)%L==1 .OR. Particle(i)%L==LayTotal) then
	     porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L)=1.0-num*(4.0/3.0*pi*(d/2.0)**3)/(18.0*dx*dy*dz)
      Else
         porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L)=1.0-num*(4.0/3.0*pi*(d/2.0)**3)/(27.0*dx*dy*dz)
	  EndIf
	  If (porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L)<0.4) porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L)=0.4
	  ReynoldsX=density*d*porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L)*ABS(u(Particle(i)%C,Particle(i)%R,Particle(i)%L)-Particle(i)%VelX)/viscosity
	  ReynoldsY=density*d*porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L)*ABS(v(Particle(i)%C,Particle(i)%R,Particle(i)%L)-Particle(i)%VelY)/viscosity
	  If (ReynoldsX>0.0) then
	     dragcoeffX=(0.63+4.8/SQRT(ReynoldsX+1.0e-30))**2
         chiX=3.7-0.65*EXP(-(1.5-LOG10(ReynoldsX+1.0e-30))**2/2.0)
         fxo=0.5*dragcoeffX*density*pi*(d/2.0)**2*(porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L))**2* &
	     ABS(u(Particle(i)%C,Particle(i)%R,Particle(i)%L)-Particle(i)%VelX)*(u(Particle(i)%C,Particle(i)%R,Particle(i)%L)-Particle(i)%VelX)
         !dragX=6.0*pi*viscosity*(u(Particle(i)%C,Particle(i)%R)-Particle(i)%VelX)*d/2.0 !Stokes Law
	     dragX=fxo*porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L)**(-(chiX+1.0))

		 !Wen and Yu Correlation
!		 If (ReynoldsX<1000) then
!		    Cd=24.0/(ReynoldsX+1.0e-30)*(1.0+0.15*ReynoldsX**0.687)
!		 Else
!		    Cd=0.44
!		 EndIf
!		 Cdprime=Cd*(porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L))**(-4.7)
!		 dragX=(1.0/8.0)*pi*(d**2)*Cdprime*density*((porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L))**2)* &
!		 ABS(u(Particle(i)%C,Particle(i)%R,Particle(i)%L)-Particle(i)%VelX)*(u(Particle(i)%C,Particle(i)%R,Particle(i)%L)-Particle(i)%VelX)

          DragForceX(Particle(i)%C,Particle(i)%R,Particle(i)%L)=DragForceX(Particle(i)%C,Particle(i)%R,Particle(i)%L)+dragX
	  EndIf
	  If (ReynoldsY>0.0) then
	     dragcoeffY=(0.63+4.8/SQRT(ReynoldsY+1.0e-30))**2
	     chiY=3.7-0.65*EXP(-(1.5-LOG10(ReynoldsY+1.0e-30))**2/2.0)
	     fyo=0.5*dragcoeffY*density*pi*(d/2.0)**2*(porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L))**2* &
	     ABS(v(Particle(i)%C,Particle(i)%R,Particle(i)%L)-Particle(i)%VelY)*(v(Particle(i)%C,Particle(i)%R,Particle(i)%L)-Particle(i)%VelY)
		 !dragY=6.0*pi*viscosity*(v(Particle(i)%C,Particle(i)%R)-Particle(i)%VelY)*d/2.0 !Stokes Law
	     dragY=fyo*porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L)**(-(chiY+1.0))

		 !Wen and Yu Correlation
!		 If (ReynoldsY<1000) then
!		    Cd=24.0/(ReynoldsY+1.0e-30)*(1.0+0.15*ReynoldsY**0.687)
!		 Else
!		    Cd=0.44
!		 EndIf
!		 Cdprime=Cd*(porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L))**(-4.7)
!		 dragY=(1.0/8.0)*pi*(d**2)*Cdprime*density*((porosity(Particle(i)%C,Particle(i)%R,Particle(i)%L))**2)* &
!		 ABS(v(Particle(i)%C,Particle(i)%R,Particle(i)%L)-Particle(i)%VelY)*(v(Particle(i)%C,Particle(i)%R,Particle(i)%L)-Particle(i)%VelY)

	     DragForceY(Particle(i)%C,Particle(i)%R,Particle(i)%L)=DragForceY(Particle(i)%C,Particle(i)%R,Particle(i)%L)+dragY
	  EndIf
	  Particle(i)%drag=SQRT(dragX**2+dragY**2)

	  !van der Waals Force
!	  Do Row=(Particle(i)%R-1),(Particle(i)%R+1)
!	     Do Col=(Particle(i)%C-1),(Particle(i)%C+1)
!            Do k=1,MaxParticle
!               If (Cell(Row,Col,k)/=0 .AND. Cell(Row,Col,k)>i) then
!	              j=Cell(Row,Col,k)
!	              h=SQRT((Particle(i)%x-Particle(j)%x)**2+(Particle(i)%y-Particle(j)%y)**2)-d
!		          If (h>5.0e-9) then
!				     e1=(Particle(j)%x-Particle(i)%x)/(h+d)
!		             e2=(Particle(j)%y-Particle(i)%y)/(h+d)
!		             Fv=Ha/6.0*(d**6/((h**2+2.0*d*h)**2*(h+d)**3))
!		             Particle(i)%FvX=Particle(i)%FvX+Fv*e1
!		             Particle(i)%FvY=Particle(i)%FvY+Fv*e2
!					 Particle(j)%FvX=Particle(j)%FvX-Fv*e1
!		             Particle(j)%FvY=Particle(j)%FvY-Fv*e2
!				  EndIf
!	           EndIf
!			EndDo
!		 EndDo
!	  EndDo
	  !DLVO
!	  Do Row=(Particle(i)%R-1),(Particle(i)%R+1)
!	     Do Col=(Particle(i)%C-1),(Particle(i)%C+1)
!            Do k=1,MaxParticle
!               If (Cell(Row,Col,k)/=0 .AND. Cell(Row,Col,k)>i) then
!	              j=Cell(Row,Col,k)
!	              h=SQRT((Particle(i)%x-Particle(j)%x)**2+(Particle(i)%y-Particle(j)%y)**2)
!		          If (h>d+5.0e-9) then
!				     e1=(Particle(j)%x-Particle(i)%x)/h
!		             e2=(Particle(j)%y-Particle(i)%y)/h
!		             Fdlvo=-Zcharge**2*lambdaB*Boltz*Temperature*(exp(kappa*d/2.0)/(1.0+kappa*d/2.0))**2*exp(-kappa*h)/h*(-kappa-1.0/h)
!		             Particle(i)%FdlvoX=Particle(i)%FdlvoX-Fdlvo*e1
!		             Particle(i)%FdlvoY=Particle(i)%FdlvoY-Fdlvo*e2
!					 Particle(j)%FdlvoX=Particle(j)%FdlvoX+Fdlvo*e1
!		             Particle(j)%FdlvoY=Particle(j)%FdlvoY+Fdlvo*e2
!				  EndIf
!	           EndIf
!			EndDo
!		 EndDo
!	  EndDo
      !Dipole-dipole Force
!	  Do Row=(Particle(i)%R-1),(Particle(i)%R+1) 
!	     Do Col=(Particle(i)%C-1),(Particle(i)%C+1)
!            Do k=1,MaxParticle
!               If (Cell(Row,Col,k)/=0 .AND. Cell(Row,Col,k)>i) then
!	              j=Cell(Row,Col,k)
!	              h=SQRT((Particle(i)%x-Particle(j)%x)**2+(Particle(i)%y-Particle(j)%y)**2)
!				  If (h>d+5.0e-9) then
!		             e1=(Particle(j)%x-Particle(i)%x)
!		             e2=(Particle(j)%y-Particle(i)%y)
!		             Particle(i)%FddX=Particle(i)%FddX+1.0/(4.0*pi*mu)*(-3.0*dipolemoment**2*e1/h**5+15.0*dipolemoment**2*e2**2*e1/h**7)
!		             Particle(i)%FddY=Particle(i)%FddY+1.0/(4.0*pi*mu)*(-9.0*dipolemoment**2*e2/h**5+15.0*dipolemoment**2*e2**3/h**7)
!					 Particle(j)%FddX=Particle(j)%FddX-1.0/(4.0*pi*mu)*(-3.0*dipolemoment**2*e1/h**5+15.0*dipolemoment**2*e2**2*e1/h**7)
!		             Particle(j)%FddY=Particle(j)%FddY-1.0/(4.0*pi*mu)*(-9.0*dipolemoment**2*e2/h**5+15.0*dipolemoment**2*e2**3/h**7)
!				  EndIf
!	           EndIf
!			EndDo
!		 EndDo
!	  EndDo

	  !Long Range Forces
	  Do Row=(Particle(i)%R-1),(Particle(i)%R+1)
	     Do Col=(Particle(i)%C-1),(Particle(i)%C+1)
		 Do Lay=(Particle(i)%L-1),(Particle(i)%L+1)
           Do k=1,MaxParticle
               If (Cell(Row,Col,Lay,k)/=0 .AND. Cell(Row,Col,Lay,k)/=i .AND. Cell(Row,Col,Lay,k)>i) then
	              j=Cell(Row,Col,Lay,k)
	              h=SQRT((Particle(i)%x-Particle(j)%x)**2+(Particle(i)%y-Particle(j)%y)**2+(Particle(i)%z-Particle(j)%z)**2)-d
		          If (h>0.0 .AND. h/2.0<=hcpp) then
				     e1=(Particle(j)%x-Particle(i)%x)/(h+d) !Particle-Particle Capillary Force
		             e2=(Particle(j)%y-Particle(i)%y)/(h+d)
					 e3=(Particle(j)%z-Particle(i)%z)/(h+d)
		             Fcap=(pi*(d/2.0)*surften)*(EXP(App*h/(2.0*(d/2.0))+Bpp)+Cpp)
		             Particle(i)%FcapX=Particle(i)%FcapX+Fcap*e1
		             Particle(i)%FcapY=Particle(i)%FcapY+Fcap*e2
					 Particle(i)%FcapZ=Particle(i)%FcapZ+Fcap*e3
					 Particle(j)%FcapX=Particle(j)%FcapX-Fcap*e1
		             Particle(j)%FcapY=Particle(j)%FcapY-Fcap*e2
					 Particle(j)%FcapZ=Particle(j)%FcapZ-Fcap*e3

!				     e1=(Particle(j)%x-Particle(i)%x)/(h+d) !van der Waals Force
!		             e2=(Particle(j)%y-Particle(i)%y)/(h+d)
!					 e3=(Particle(j)%z-Particle(i)%z)/(h+d)
!		             Fv=Ha/6.0*(d**6/((h**2+2.0*d*h)**2*(h+d)**3))
!		             Particle(i)%FvX=Particle(i)%FvX+Fv*e1
!		             Particle(i)%FvY=Particle(i)%FvY+Fv*e2
!					 Particle(i)%FvZ=Particle(i)%FvZ+Fv*e3
!					 Particle(j)%FvX=Particle(j)%FvX-Fv*e1
!		             Particle(j)%FvY=Particle(j)%FvY-Fv*e2
!					 Particle(j)%FvZ=Particle(j)%FvZ-Fv*e3
!					 h=h+d
!					 Fdlvo=-Zcharge**2*lambdaB*Boltz*Temperature*(exp(kappa*d/2.0)/(1.0+kappa*d/2.0))**2*exp(-kappa*h)/h*(-kappa-1.0/h) !DLVO
!		             Particle(i)%FdlvoX=Particle(i)%FdlvoX-Fdlvo*e1
!		             Particle(i)%FdlvoY=Particle(i)%FdlvoY-Fdlvo*e2
!					 Particle(i)%FdlvoZ=Particle(i)%FdlvoZ-Fdlvo*e3
!					 Particle(j)%FdlvoX=Particle(j)%FdlvoX+Fdlvo*e1
!		             Particle(j)%FdlvoY=Particle(j)%FdlvoY+Fdlvo*e2
!					 Particle(j)%FdlvoZ=Particle(j)%FdlvoZ+Fdlvo*e3
!					 e1=(Particle(j)%x-Particle(i)%x) !Dipole-dipole Force
!		             e2=(Particle(j)%y-Particle(i)%y)
!					 e3=(Particle(j)%z-Particle(i)%z)
!		             Particle(i)%FddX=Particle(i)%FddX+1.0/(4.0*pi*mu)*(-3.0*dipolemoment**2*e1/h**5+15.0*dipolemoment**2*e2**2*e1/h**7)
!		             Particle(i)%FddY=Particle(i)%FddY+1.0/(4.0*pi*mu)*(-9.0*dipolemoment**2*e2/h**5+15.0*dipolemoment**2*e2**3/h**7)
!					 Particle(i)%FddZ=Particle(i)%FddZ+1.0/(4.0*pi*mu)*(-3.0*dipolemoment**2*e3/h**5+15.0*dipolemoment**2*e2**2*e3/h**7)
!					 Particle(j)%FddX=Particle(j)%FddX-1.0/(4.0*pi*mu)*(-3.0*dipolemoment**2*e1/h**5+15.0*dipolemoment**2*e2**2*e1/h**7)
!		             Particle(j)%FddY=Particle(j)%FddY-1.0/(4.0*pi*mu)*(-9.0*dipolemoment**2*e2/h**5+15.0*dipolemoment**2*e2**3/h**7)
!					 Particle(j)%FddZ=Particle(j)%FddZ-1.0/(4.0*pi*mu)*(-3.0*dipolemoment**2*e3/h**5+15.0*dipolemoment**2*e2**2*e3/h**7)
				  EndIf
	           EndIf
			EndDo
		 EndDo
		 EndDo
	  EndDo
	  
	  !Particle-Wall Capillary Force
	  If (Particle(i)%x>d/2.0 .AND. Particle(i)%x<hcpw+d/2.0) then
	     h=Particle(i)%x-d/2.0
	     Particle(i)%FcapX=Particle(i)%FcapX-(pi*(d/2.0)*surften)*(EXP(Apw*h/(d/2.0)+Bpw)+Cpw)
	  ElseIf (Particle(i)%x<Diameter-d/2.0 .AND. Particle(i)%x>Diameter-hcpw-d/2.0) then
	     h=Diameter-Particle(i)%x-d/2.0
	     Particle(i)%FcapX=Particle(i)%FcapX+(pi*(d/2.0)*surften)*(EXP(Apw*h/(d/2.0)+Bpw)+Cpw)
	  EndIf
	  If (Particle(i)%y>d/2.0 .AND. Particle(i)%y<hcpw+d/2.0) then
	     h=Particle(i)%y-d/2.0
	     Particle(i)%FcapY=Particle(i)%FcapY-(pi*(d/2.0)*surften)*(EXP(Apw*h/(d/2.0)+Bpw)+Cpw)
	  EndIf
	  If (Particle(i)%z>d/2.0 .AND. Particle(i)%z<hcpw+d/2.0) then
	     h=Particle(i)%z-d/2.0
	     Particle(i)%FcapZ=Particle(i)%FcapZ-(pi*(d/2.0)*surften)*(EXP(Apw*h/(d/2.0)+Bpw)+Cpw)
	  ElseIf (Particle(i)%z<Depth-d/2.0 .AND. Particle(i)%z>Depth-hcpw-d/2.0) then
	     h=Depth-Particle(i)%z-d/2.0
	     Particle(i)%FcapZ=Particle(i)%FcapZ+(pi*(d/2.0)*surften)*(EXP(Apw*h/(d/2.0)+Bpw)+Cpw)
	  EndIf
	  Particle(i)%Fcap=SQRT(Particle(i)%FcapX**2+Particle(i)%FcapY**2+Particle(i)%FcapZ**2)

	  !Brownian Force
!	  Call Random_Number(random1)
!	  random1=random1-0.5
!	  Call Random_Number(random2)
!	  random2=random2-0.5
!	  Call Random_Number(random3)
!	  random3=random3-0.5
!	  FbX=random1*SQRT(12.0*pi*viscosity*(d/2.0)*Boltz*Temperature*t)
!	  FbY=random2*SQRT(12.0*pi*viscosity*(d/2.0)*Boltz*Temperature*t)
!	  FbZ=random3*SQRT(12.0*pi*viscosity*(d/2.0)*Boltz*Temperature*t)

  	  If (.NOT. Particle(i)%BaseCollision) then
	     Particle(i)%BaseFx=0.0
	     Particle(i)%BaseFy=0.0
		 Particle(i)%BaseFz=0.0
	     Particle(i)%BaseMomentx=0.0
		 Particle(i)%BaseMomenty=0.0
		 Particle(i)%BaseMomentz=0.0
		 Particle(i)%BaseContactFx=0.0
	     Particle(i)%BaseContactFy=0.0
		 Particle(i)%BaseContactFz=0.0
	     Particle(i)%BaseContactMomentx=0.0
		 Particle(i)%BaseContactMomenty=0.0
		 Particle(i)%BaseContactMomentz=0.0
	  EndIf
  	  If (.NOT. Particle(i)%Wall1Collision) then
	     Particle(i)%Wall1Fx=0.0
	     Particle(i)%Wall1Fy=0.0
		 Particle(i)%Wall1Fz=0.0
	     Particle(i)%Wall1Momentx=0.0
		 Particle(i)%Wall1Momenty=0.0
		 Particle(i)%Wall1Momentz=0.0
		 Particle(i)%Wall1ContactFx=0.0
	     Particle(i)%Wall1ContactFy=0.0
		 Particle(i)%Wall1ContactFz=0.0
	     Particle(i)%Wall1ContactMomentx=0.0
		 Particle(i)%Wall1ContactMomenty=0.0
		 Particle(i)%Wall1ContactMomentz=0.0
	  EndIf
	  If (.NOT. Particle(i)%Wall2Collision) then
	     Particle(i)%Wall2Fx=0.0
	     Particle(i)%Wall2Fy=0.0
		 Particle(i)%Wall2Fz=0.0
	     Particle(i)%Wall2Momentx=0.0
		 Particle(i)%Wall2Momenty=0.0
		 Particle(i)%Wall2Momentz=0.0
		 Particle(i)%Wall2ContactFx=0.0
	     Particle(i)%Wall2ContactFy=0.0
		 Particle(i)%Wall2ContactFz=0.0
	     Particle(i)%Wall2ContactMomentx=0.0
		 Particle(i)%Wall2ContactMomenty=0.0
		 Particle(i)%Wall2ContactMomentz=0.0
	  EndIf
  	  If (.NOT. Particle(i)%ParticleCollision) then
	     Particle(i)%ParticleFx=0.0
	     Particle(i)%ParticleFy=0.0
		 Particle(i)%ParticleFz=0.0
	     Particle(i)%ParticleMomentx=0.0
		 Particle(i)%ParticleMomenty=0.0
		 Particle(i)%ParticleMomentz=0.0
		 Particle(i)%ParticleContactFx=0.0
	     Particle(i)%ParticleContactFy=0.0
		 Particle(i)%ParticleContactFz=0.0
	     Particle(i)%ParticleContactMomentx=0.0
		 Particle(i)%ParticleContactMomenty=0.0
		 Particle(i)%ParticleContactMomentz=0.0
	  EndIf
	  Particle(i)%Fx=Particle(i)%BaseFx+Particle(i)%Wall1Fx+Particle(i)%Wall2Fx+Particle(i)%ParticleFx-Particle(i)%m*g*COS(angle*pi/180.0)+dragX+Particle(i)%FcapX !+Particle(i)%FvX+FbX+Particle(i)%FddX+Particle(i)%FdlvoX
	  Particle(i)%Fy=Particle(i)%BaseFy+Particle(i)%Wall1Fy+Particle(i)%Wall2Fy+Particle(i)%ParticleFy-Particle(i)%m*g*SIN(angle*pi/180.0)+dragY+Particle(i)%FcapY !+Particle(i)%FvY+FbY+Particle(i)%FddY+Particle(i)%FdlvoY
	  Particle(i)%Fz=Particle(i)%BaseFz+Particle(i)%Wall1Fz+Particle(i)%Wall2Fz+Particle(i)%ParticleFz+Particle(i)%FcapZ !+Particle(i)%FvZ+FbZ+Particle(i)%FddZ+Particle(i)%FdlvoZ !+dragZ
      Particle(i)%Momentx=Particle(i)%BaseMomentx+Particle(i)%Wall1Momentx+Particle(i)%Wall2Momentx+Particle(i)%ParticleMomentx
	  Particle(i)%Momenty=Particle(i)%BaseMomenty+Particle(i)%Wall1Momenty+Particle(i)%Wall2Momenty+Particle(i)%ParticleMomenty
	  Particle(i)%Momentz=Particle(i)%BaseMomentz+Particle(i)%Wall1Momentz+Particle(i)%Wall2Momentz+Particle(i)%ParticleMomentz
!	  Particle(i)%FvX=0.0
!	  Particle(i)%FvY=0.0
!	  Particle(i)%FvZ=0.0
!	  Particle(i)%FddX=0.0
!	  Particle(i)%FddY=0.0
!	  Particle(i)%FddZ=0.0
!	  Particle(i)%FdlvoX=0.0
!	  Particle(i)%FdlvoY=0.0
!	  Particle(i)%FdlvoZ=0.0
      Particle(i)%FcapX=0.0
	  Particle(i)%FcapY=0.0
	  Particle(i)%FcapZ=0.0

   !$OMP End parallel Do	  
	  !Verlet Algorithm
      !xnew=2.0*Particle(i)%x-Particle(i)%xold+(Particle(i)%Fx/m)*t**2
      !Particle(i)%VelX=(xnew-Particle(i)%xold)/(2.0*t)
      !Particle(i)%xold=Particle(i)%x
      !Particle(i)%x=xnew
      !ynew=2.0*Particle(i)%y-Particle(i)%yold+(Particle(i)%Fy/m)*t**2
      !Particle(i)%VelY=(ynew-Particle(i)%yold)/(2.0*t)
      !Particle(i)%yold=Particle(i)%y
      !Particle(i)%y=ynew
	  !znew=2.0*Particle(i)%z-Particle(i)%zold+(Particle(i)%Fz/m)*t**2
      !Particle(i)%VelZ=(znew-Particle(i)%zold)/(2.0*t)
      !Particle(i)%zold=Particle(i)%z
      !Particle(i)%z=znew
	  !rotationxnew=2.0*Particle(i)%rotationx-Particle(i)%rotationxold+(Particle(i)%Momentx/MomentInert)*t**2
      !Particle(i)%AngularVelx=(rotationxnew-Particle(i)%rotationxold)/(2.0*t)
      !Particle(i)%rotationxold=Particle(i)%rotationx
      !Particle(i)%rotationx=rotationxnew
	  !rotationynew=2.0*Particle(i)%rotationy-Particle(i)%rotationyold+(Particle(i)%Momenty/MomentInert)*t**2
      !Particle(i)%AngularVely=(rotationynew-Particle(i)%rotationyold)/(2.0*t)
      !Particle(i)%rotationyold=Particle(i)%rotationy
      !Particle(i)%rotationy=rotationynew
	  !rotationznew=2.0*Particle(i)%rotationz-Particle(i)%rotationzold+(Particle(i)%Momentz/MomentInert)*t**2
      !Particle(i)%AngularVelz=(rotationznew-Particle(i)%rotationzold)/(2.0*t)
      !Particle(i)%rotationzold=Particle(i)%rotationz
      !Particle(i)%rotationz=rotationznew
   
	  !2nd Order Adams-Bashforth Scheme
	  Particle(i)%x=Particle(i)%x+t*(1.5*Particle(i)%Velx-0.5*Particle(i)%Velxold)
	  Velxnew=Particle(i)%Velx+t*(1.5*Particle(i)%Fx/Particle(i)%m-0.5*Particle(i)%Accelxold)
	  Particle(i)%Velxold=Particle(i)%Velx
	  Particle(i)%Velx=Velxnew
	  Particle(i)%Accelxold=Particle(i)%Fx/Particle(i)%m
	  Particle(i)%y=Particle(i)%y+t*(1.5*Particle(i)%Vely-0.5*Particle(i)%Velyold)
	  Velynew=Particle(i)%Vely+t*(1.5*Particle(i)%Fy/Particle(i)%m-0.5*Particle(i)%Accelyold)
	  Particle(i)%Velyold=Particle(i)%Vely
	  Particle(i)%Vely=Velynew
	  Particle(i)%Accelyold=Particle(i)%Fy/Particle(i)%m
	  Particle(i)%z=Particle(i)%z+t*(1.5*Particle(i)%Velz-0.5*Particle(i)%Velzold)
	  Velznew=Particle(i)%Velz+t*(1.5*Particle(i)%Fz/Particle(i)%m-0.5*Particle(i)%Accelzold)
	  Particle(i)%Velzold=Particle(i)%Velz
	  Particle(i)%Velz=Velznew
	  Particle(i)%Accelzold=Particle(i)%Fz/Particle(i)%m
	  Particle(i)%rotationx=Particle(i)%rotationx+t*(1.5*Particle(i)%AngularVelx-0.5*Particle(i)%AngularVelxold)
	  AngularVelxnew=Particle(i)%AngularVelx+t*(1.5*Particle(i)%Momentx/Particle(i)%MomentInert-0.5*Particle(i)%AngularAccelxold)
	  Particle(i)%AngularVelxold=Particle(i)%AngularVelx
	  Particle(i)%AngularVelx=AngularVelxnew
	  Particle(i)%AngularAccelxold=Particle(i)%Momentx/Particle(i)%MomentInert
	  Particle(i)%rotationy=Particle(i)%rotationy+t*(1.5*Particle(i)%AngularVely-0.5*Particle(i)%AngularVelyold)
	  AngularVelynew=Particle(i)%AngularVely+t*(1.5*Particle(i)%Momenty/Particle(i)%MomentInert-0.5*Particle(i)%AngularAccelyold)
	  Particle(i)%AngularVelyold=Particle(i)%AngularVely
	  Particle(i)%AngularVely=AngularVelynew
	  Particle(i)%AngularAccelyold=Particle(i)%Momenty/Particle(i)%MomentInert
	  Particle(i)%rotationz=Particle(i)%rotationz+t*(1.5*Particle(i)%AngularVelz-0.5*Particle(i)%AngularVelzold)
	  AngularVelznew=Particle(i)%AngularVelz+t*(1.5*Particle(i)%Momentz/Particle(i)%MomentInert-0.5*Particle(i)%AngularAccelzold)
	  Particle(i)%AngularVelzold=Particle(i)%AngularVelz
	  Particle(i)%AngularVelz=AngularVelznew
	  Particle(i)%AngularAccelzold=Particle(i)%Momentz/Particle(i)%MomentInert

	  !Periodic Boundary Conditions
!	  If (Particle(i)%y>Length-d/2.0 .AND. Particle(i)%VelY>0.0) then
!	     Particle(i)%y=d/2.0
!	  	 overlap=.True.
!	  	 Do While (overlap)
!	  	    RowTemp=INT(Particle(i)%y/(Length/RowTotal))+1
!	        ColTemp=INT(Particle(i)%x/(Diameter/ColTotal))+1
!			LayTemp=INT(Particle(i)%z/(Depth/LayTotal))+1
!	  	 	overlap=.False.
!	  		Do RowTemp2=MAX(RowTemp-1,1),MIN(RowTemp+1,RowTotal)
!	  		   Do ColTemp2=MAX(ColTemp-1,1),MIN(ColTemp+1,ColTotal)
!			   Do LayTemp2=MAX(LayTemp-1,1),MIN(LayTemp+1,LayTotal)
!	  	 	      Do k=1,MaxParticle
!	  	 	         j=Cell(RowTemp2,ColTemp2,LayTemp2,k)
!	  	 	         If (j/=0 .AND. j/=i .AND. (Particle(i)%x-Particle(j)%x)**2+(Particle(i)%y-Particle(j)%y)**2+(Particle(i)%z-Particle(j)%z)**2<d**2) then
!	  	 	            overlap=.True.
!	  	 		        Particle(i)%y=Particle(i)%y+d
!	  	 		        Exit
!	  	 	         EndIf
!	  		      EndDo
!	  		      If (overlap) Exit
!	  		   EndDo
!			   If (overlap) Exit
!	  		   EndDo
!	  		   If (overlap) Exit
!	  	 	EndDo
!	  	 EndDo
!	  	 Particle(i)%yold=Particle(i)%y-Particle(i)%VelY*t
	  	 !Gs=Gs+1
!	  EndIf
!	  If (Particle(i)%y<d/2.0 .AND. Particle(i)%VelY<0.0) then
!	     Particle(i)%y=Length-d/2.0
!	  	 overlap=.True.
!	  	 Do While (overlap)
!	  	    RowTemp=INT(Particle(i)%y/(Length/RowTotal))+1
!	        ColTemp=INT(Particle(i)%x/(Diameter/ColTotal))+1
!			LayTemp=INT(Particle(i)%z/(Depth/LayTotal))+1
!	  	 	overlap=.False.
!	  		Do RowTemp2=MAX(RowTemp-1,1),MIN(RowTemp+1,RowTotal)
!	  		   Do ColTemp2=MAX(ColTemp-1,1),MIN(ColTemp+1,ColTotal)
!			   Do LayTemp2=MAX(LayTemp-1,1),MIN(LayTemp+1,LayTotal)
!	  	 	      Do k=1,MaxParticle
!	  	 	         j=Cell(RowTemp2,ColTemp2,LayTemp2,k)
!	  	 	         If (j/=0 .AND. j/=i .AND. (Particle(i)%x-Particle(j)%x)**2+(Particle(i)%y-Particle(j)%y)**2+(Particle(i)%z-Particle(j)%z)**2<d**2) then
!	  	 	            overlap=.True.
!	  	 		        Particle(i)%y=Particle(i)%y-d
!	  	 		        Exit
!	  	 	         EndIf
!	  		      EndDo
!	  		      If (overlap) Exit
!	  		   EndDo
!			   If (overlap) Exit
!	  		   EndDo
!	  		   If (overlap) Exit
!	  	 	EndDo
!	  	 EndDo
!	  	 Particle(i)%yold=Particle(i)%y-Particle(i)%VelY*t
	  	 !Gs=Gs+1
!	  EndIf
!	  If (Particle(i)%x>Diameter-d/2.0 .AND. Particle(i)%VelX>0.0) then
!	     Particle(i)%x=d/2.0
!	  	 overlap=.True.
!	  	 Do While (overlap)
!	  	    RowTemp=INT(Particle(i)%y/(Length/RowTotal))+1
!	        ColTemp=INT(Particle(i)%x/(Diameter/ColTotal))+1
!			LayTemp=INT(Particle(i)%z/(Depth/LayTotal))+1
!	  	 	overlap=.False.
!	  		Do RowTemp2=MAX(RowTemp-1,1),MIN(RowTemp+1,RowTotal)
!	  		   Do ColTemp2=MAX(ColTemp-1,1),MIN(ColTemp+1,ColTotal)
!			   Do LayTemp2=MAX(LayTemp-1,1),MIN(LayTemp+1,LayTotal)
!	  	 	      Do k=1,MaxParticle
!	  	 	         j=Cell(RowTemp2,ColTemp2,LayTemp2,k)
!	  	 	         If (j/=0 .AND. j/=i .AND. (Particle(i)%x-Particle(j)%x)**2+(Particle(i)%y-Particle(j)%y)**2+(Particle(i)%z-Particle(j)%z)**2<d**2) then
!	  	 	            overlap=.True.
!	  	 		        Particle(i)%x=Particle(i)%x+d
!	  	 		        Exit
!	  	 	         EndIf
!	  		      EndDo
!	  		      If (overlap) Exit
!	  		   EndDo
!			   If (overlap) Exit
!	  		   EndDo
!	  		   If (overlap) Exit
!	  	 	EndDo
!	  	 EndDo
!	  	 Particle(i)%xold=Particle(i)%x-Particle(i)%VelX*t
	  	 !Gs=Gs+1
!	  EndIf
!	  If (Particle(i)%x<d/2.0 .AND. Particle(i)%VelX<0.0) then
!	     Particle(i)%x=Diameter-d/2.0
!	  	 overlap=.True.
!	  	 Do While (overlap)
!	  	    RowTemp=INT(Particle(i)%y/(Length/RowTotal))+1
!	        ColTemp=INT(Particle(i)%x/(Diameter/ColTotal))+1
!			LayTemp=INT(Particle(i)%z/(Depth/LayTotal))+1
!	  	 	overlap=.False.
!	  		Do RowTemp2=MAX(RowTemp-1,1),MIN(RowTemp+1,RowTotal)
!	  		   Do ColTemp2=MAX(ColTemp-1,1),MIN(ColTemp+1,ColTotal)
!			   Do LayTemp2=MAX(LayTemp-1,1),MIN(LayTemp+1,LayTotal)
!	  	 	      Do k=1,MaxParticle
!	  	 	         j=Cell(RowTemp2,ColTemp2,LayTemp2,k)
!	  	 	         If (j/=0 .AND. j/=i .AND. (Particle(i)%x-Particle(j)%x)**2+(Particle(i)%y-Particle(j)%y)**2+(Particle(i)%z-Particle(j)%z)**2<d**2) then
!	  	 	            overlap=.True.
!	  	 		        Particle(i)%x=Particle(i)%x-d
!	  	 		        Exit
!	  	 	         EndIf
!	  		      EndDo
!	  		      If (overlap) Exit
!	  		   EndDo
!			   If (overlap) Exit
!	  		   EndDo
!	  		   If (overlap) Exit
!	  	 	EndDo
!	  	 EndDo
!	  	 Particle(i)%xold=Particle(i)%x-Particle(i)%VelX*t
	  	 !Gs=Gs+1
!	  EndIf
!	  If (Particle(i)%z>Depth-d/2.0 .AND. Particle(i)%VelZ>0.0) then
!	     Particle(i)%z=d/2.0
!	  	 overlap=.True.
!	  	 Do While (overlap)
!	  	    RowTemp=INT(Particle(i)%y/(Length/RowTotal))+1
!	        ColTemp=INT(Particle(i)%x/(Diameter/ColTotal))+1
!			LayTemp=INT(Particle(i)%z/(Depth/LayTotal))+1
!	  	 	overlap=.False.
!	  		Do RowTemp2=MAX(RowTemp-1,1),MIN(RowTemp+1,RowTotal)
!	  		   Do ColTemp2=MAX(ColTemp-1,1),MIN(ColTemp+1,ColTotal)
!			   Do LayTemp2=MAX(LayTemp-1,1),MIN(LayTemp+1,LayTotal)
!	  	 	      Do k=1,MaxParticle
!	  	 	         j=Cell(RowTemp2,ColTemp2,LayTemp2,k)
!	  	 	         If (j/=0 .AND. j/=i .AND. (Particle(i)%x-Particle(j)%x)**2+(Particle(i)%y-Particle(j)%y)**2+(Particle(i)%z-Particle(j)%z)**2<d**2) then
!	  	 	            overlap=.True.
!	  	 		        Particle(i)%z=Particle(i)%z+d
!	  	 		        Exit
!	  	 	         EndIf
!	  		      EndDo
!	  		      If (overlap) Exit
!	  		   EndDo
!			   If (overlap) Exit
!	  		   EndDo
!	  		   If (overlap) Exit
!	  	 	EndDo
!	  	 EndDo
!	  	 Particle(i)%zold=Particle(i)%z-Particle(i)%VelZ*t
	  	 !Gs=Gs+1
!	  EndIf
!	  If (Particle(i)%z<d/2.0 .AND. Particle(i)%VelZ<0.0) then
!	     Particle(i)%z=Depth-d/2.0
!	  	 overlap=.True.
!	  	 Do While (overlap)
!	  	    RowTemp=INT(Particle(i)%y/(Length/RowTotal))+1
!	        ColTemp=INT(Particle(i)%x/(Diameter/ColTotal))+1
!			LayTemp=INT(Particle(i)%z/(Depth/LayTotal))+1
!	  	 	overlap=.False.
!	  		Do RowTemp2=MAX(RowTemp-1,1),MIN(RowTemp+1,RowTotal)
!	  		   Do ColTemp2=MAX(ColTemp-1,1),MIN(ColTemp+1,ColTotal)
!			   Do LayTemp2=MAX(LayTemp-1,1),MIN(LayTemp+1,LayTotal)
!	  	 	      Do k=1,MaxParticle
!	  	 	         j=Cell(RowTemp2,ColTemp2,LayTemp2,k)
!	  	 	         If (j/=0 .AND. j/=i .AND. (Particle(i)%x-Particle(j)%x)**2+(Particle(i)%y-Particle(j)%y)**2+(Particle(i)%z-Particle(j)%z)**2<d**2) then
!	  	 	            overlap=.True.
!	  	 		        Particle(i)%z=Particle(i)%z-d
!	  	 		        Exit
!	  	 	         EndIf
!	  		      EndDo
!	  		      If (overlap) Exit
!	  		   EndDo
!			   If (overlap) Exit
!	  		   EndDo
!	  		   If (overlap) Exit
!	  	 	EndDo
!	  	 EndDo
!	  	 Particle(i)%zold=Particle(i)%z-Particle(i)%VelZ*t
	  	 !Gs=Gs+1
!	  EndIf
	  RowTemp=INT(Particle(i)%y/(Length/RowTotal))+1
	  ColTemp=INT(Particle(i)%x/(Diameter/ColTotal))+1
	  LayTemp=INT(Particle(i)%z/(Depth/LayTotal))+1
	  CellNo=(LayTemp-1)*RowTotal*ColTotal+(RowTemp-1)*ColTotal+ColTemp
      If (Particle(i)%CellNum/=CellNo) then
         k=1
	     Do While (Cell(Particle(i)%R,Particle(i)%C,Particle(i)%L,k)/=i)
	        k=k+1
	     EndDo
	     Cell(Particle(i)%R,Particle(i)%C,Particle(i)%L,k)=0
	     k=1
	     Do While (Cell(RowTemp,ColTemp,LayTemp,k)/=0)
	        k=k+1
	     EndDo
	     Cell(RowTemp,ColTemp,LayTemp,k)=i
	     Particle(i)%CellNum=CellNo
	     Particle(i)%R=RowTemp
	     Particle(i)%C=ColTemp
		 Particle(i)%L=LayTemp
      EndIf
   EndDo
!If (MOD(INT(time/t),INT(CFDtime/t))==0) then
Do Lay=1,LayTotal
   Do iterate=1,MaxIterate !CFD
      !Pressure Boundary Values
      Do Col=1,ColTotal
         P(Col,0,Lay)=P(Col,1,Lay)+(P(Col,1,Lay)-P(Col,2,Lay))/2.0
	     P(Col,RowTotal+1,Lay)=P(Col,RowTotal,Lay)+(P(Col,RowTotal,Lay)-P(Col,RowTotal-1,Lay))/2.0
      EndDo
      Do Row=1,RowTotal
         P(0,Row,Lay)=P(1,Row,Lay)+(P(1,Row,Lay)-P(2,Row,Lay))/2.0
	     P(ColTotal+1,Row,Lay)=P(ColTotal,Row,Lay)+(P(ColTotal,Row,Lay)-P(ColTotal-1,Row,Lay))/2.0
      EndDo
      SU=0.0
      SV=0.0
      APU=0.0
      APV=0.0
      !East Fluxes
      Do Col=1,ColTotal-1
         Do Row=1,RowTotal
		    CE=MIN(REAL(FlowX(Col,Row,Lay)),0.0)
		    CP=MAX(REAL(FlowX(Col,Row,Lay)),0.0)
		    FluxUUDS=CP*u(Col,Row,Lay)+CE*u(Col+1,Row,Lay)
		    FluxVUDS=CP*v(Col,Row,Lay)+CE*v(Col+1,Row,Lay)
		    FluxUCDS=FlowX(Col,Row,Lay)*(u(Col+1,Row,Lay)+u(Col,Row,Lay))/2.0
		    FluxVCDS=FlowX(Col,Row,Lay)*(v(Col+1,Row,Lay)+v(Col,Row,Lay))/2.0
		    AE(Col,Row,Lay)=(CE-viscosity)*porosity(Col,Row,Lay) 
		    AW(Col+1,Row,Lay)=(-CP-viscosity)*porosity(Col+1,Row,Lay) 
		    SU(Col,Row,Lay)=SU(Col,Row,Lay)+ConvScheme*(FluxUUDS-FluxUCDS)*porosity(Col,Row,Lay) 
		    SU(Col+1,Row,Lay)=SU(Col+1,Row,Lay)-ConvScheme*(FluxUUDS-FluxUCDS)*porosity(Col+1,Row,Lay) 
		    SV(Col,Row,Lay)=SV(Col,Row,Lay)+ConvScheme*(FluxVUDS-FluxVCDS)*porosity(Col,Row,Lay) 
		    SV(Col+1,Row,Lay)=SV(Col+1,Row,Lay)-ConvScheme*(FluxVUDS-FluxVCDS)*porosity(Col+1,Row,Lay) 
	     EndDo
      EndDo
      !North Fluxes
      Do Row=1,RowTotal-1
         Do Col=1,ColTotal
		    CN=MIN(REAL(FlowY(Col,Row,Lay)),0.0)
		    CP=MAX(REAL(FlowY(Col,Row,Lay)),0.0)
		    FluxUUDS=CP*u(Col,Row,Lay)+CN*u(Col,Row+1,Lay)
		    FluxVUDS=CP*v(Col,Row,Lay)+CN*v(Col,Row+1,Lay)
		    FluxUCDS=FlowY(Col,Row,Lay)*(u(Col,Row+1,Lay)+u(Col,Row,Lay))/2.0
		    FluxVCDS=FlowY(Col,Row,Lay)*(v(Col,Row+1,Lay)+v(Col,Row,Lay))/2.0
		    AN(Col,Row,Lay)=(CN-viscosity)*porosity(Col,Row,Lay) 
		    AS(Col,Row+1,Lay)=(-CP-viscosity)*porosity(Col,Row+1,Lay) 
		    SU(Col,Row,Lay)=SU(Col,Row,Lay)+ConvScheme*(FluxUUDS-FluxUCDS)*porosity(Col,Row,Lay) 
		    SU(Col,Row+1,Lay)=SU(Col,Row+1,Lay)-ConvScheme*(FluxUUDS-FluxUCDS)*porosity(Col,Row+1,Lay) 
		    SV(Col,Row,Lay)=SV(Col,Row,Lay)+ConvScheme*(FluxVUDS-FluxVCDS)*porosity(Col,Row,Lay) 
		    SV(Col,Row+1,Lay)=SV(Col,Row+1,Lay)-ConvScheme*(FluxVUDS-FluxVCDS)*porosity(Col,Row+1,Lay) 
	     EndDo
      EndDo
      !Source Terms
      Do Col=1,ColTotal
         Do Row=1,RowTotal
	        PE=P(Col+1,Row,Lay)*factorX(Col)+P(Col,Row,Lay)*(1.0-factorX(Col))
		    PW=P(Col,Row,Lay)*factorX(Col-1)+P(Col-1,Row,Lay)*(1.0-factorX(Col-1))
		    PN=P(Col,Row+1,Lay)*factorY(Row)+P(Col,Row,Lay)*(1.0-factorY(Row))
		    PS=P(Col,Row,Lay)*factorY(Row-1)+P(Col,Row-1,Lay)*(1.0-factorY(Row-1))
		    dPx(Col,Row,Lay)=(PE-PW)/dx !*porosity(Col,Row) 
		    dPy(Col,Row,Lay)=(PN-PS)/dy !*porosity(Col,Row) 
		    SU(Col,Row,Lay)=SU(Col,Row,Lay)-dPx(Col,Row,Lay)*dx*dy-DragForceX(Col,Row,Lay)-density*dx*dy*g*porosity(Col,Row,Lay)*COS(angle*pi/180.0)
		    SV(Col,Row,Lay)=SV(Col,Row,Lay)-dPy(Col,Row,Lay)*dx*dy-DragForceY(Col,Row,Lay)-density*dx*dy*g*porosity(Col,Row,Lay)*SIN(angle*pi/180.0)
		    APT=density*dx*dy/t
		    SU(Col,Row,Lay)=SU(Col,Row,Lay)+(1.0+TimeScheme)*APT*uo(Col,Row,Lay)*porosity(Col,Row,Lay)- &
			0.5*TimeScheme*APT*uoo(Col,Row,Lay)*porosity(Col,Row,Lay) 
		    SV(Col,Row,Lay)=SV(Col,Row,Lay)+(1.0+TimeScheme)*APT*vo(Col,Row,Lay)*porosity(Col,Row,Lay)- &
			0.5*TimeScheme*APT*voo(Col,Row,Lay)*porosity(Col,Row,Lay) 
		    APV(Col,Row,Lay)=APV(Col,Row,Lay)+(1.0+0.5*TimeScheme)*APT*porosity(Col,Row,Lay) 
		    APU(Col,Row,Lay)=APU(Col,Row,Lay)+(1.0+0.5*TimeScheme)*APT*porosity(Col,Row,Lay) 
	     EndDo
      EndDo
      !East Boundary
      Col=ColTotal
      Do Row=1,RowTotal
         APV(Col,Row,Lay)=APV(Col,Row,Lay)+viscosity*porosity(Col,Row,Lay) 
	     SV(Col,Row,Lay)=SV(Col,Row,Lay)+viscosity*v(Col+1,Row,Lay)*porosity(Col+1,Row,Lay) 
      EndDo
      !West Boundary
      Col=1
      Do Row=1,RowTotal
         APV(Col,Row,Lay)=APV(Col,Row,Lay)+viscosity*porosity(Col,Row,Lay) 
	     SV(Col,Row,Lay)=SV(Col,Row,Lay)+viscosity*v(Col-1,Row,Lay)*porosity(Col-1,Row,Lay) 
      EndDo
      !South Boundary
      Row=1
      Do Col=1,ColTotal
         AWC=viscosity*porosity(Col,Row-1,Lay)+FlowY(Col,Row-1,Lay) 
	     APV(Col,Row,Lay)=APV(Col,Row,Lay)+AWC
	     SU(Col,Row,Lay)=SU(Col,Row,Lay)+AWC*u(Col,Row-1,Lay)
	     SV(Col,Row,Lay)=SV(Col,Row,Lay)+AWC*v(Col,Row-1,Lay)
      EndDo
      !North Boundary
      Row=RowTotal
      Do Col=1,ColTotal
         AN(Col,Row,Lay)=0.0
      EndDo
      !Solving u
      Do Col=1,ColTotal
         Do Row=1,RowTotal
	        AP(Col,Row,Lay)=(-AE(Col,Row,Lay)-AW(Col,Row,Lay)-AN(Col,Row,Lay)-AS(Col,Row,Lay)+APU(Col,Row,Lay))/relax1
		    SU(Col,Row,Lay)=SU(Col,Row,Lay)+(1.0-relax1)*AP(Col,Row,Lay)*u(Col,Row,Lay)
		    APU(Col,Row,Lay)=1.0/(AP(Col,Row,Lay)+1.0e-30)
	     EndDo
      EndDo
      Call StoneSIP(u,Lay)
      !Solving v
      Do Col=1,ColTotal
         Do Row=1,RowTotal
	        AP(Col,Row,Lay)=(-AE(Col,Row,Lay)-AW(Col,Row,Lay)-AN(Col,Row,Lay)-AS(Col,Row,Lay)+APV(Col,Row,Lay))/relax1
		    SU(Col,Row,Lay)=SV(Col,Row,Lay)+(1.0-relax1)*AP(Col,Row,Lay)*v(Col,Row,Lay)
		    APV(Col,Row,Lay)=1.0/(AP(Col,Row,Lay)+1.0e-30)
	     EndDo
      EndDo
      Call StoneSIP(v,Lay)
      !Outlet Boundary Condition
      OutFlow=0.0
      Row=RowTotal
      Do Col=1,ColTotal
         FlowY(Col,Row,Lay)=density*dx*v(Col,Row,Lay) !*porosity(Col,Row) 
	     OutFlow(Lay)=OutFlow(Lay)+FlowY(Col,Row,Lay)
      EndDo
      Do Col=1,ColTotal
         FlowY(Col,Row,Lay)=FlowY(Col,Row,Lay)*MassFlow(Lay)/(OutFlow(Lay)+1.0e-30)
	     v(Col,Row+1,Lay)=v(Col,Row,Lay)*MassFlow(Lay)/(OutFlow(Lay)+1.0e-30)
      EndDo
      !Pressure Correction
      Do Col=1,ColTotal-1
         Do Row=1,RowTotal
	        dPxMid=(dPx(Col+1,Row,Lay)+dPx(Col,Row,Lay))/2.0
		    uMid=(u(Col+1,Row,Lay)*porosity(Col+1,Row,Lay)+u(Col,Row,Lay)*porosity(Col,Row,Lay))/2.0 
		    APUMid=(APU(Col+1,Row,Lay)+APU(Col,Row,Lay))/2.0
		    !dPxE=(P(Col+1,Row)*porosity(Col+1,Row)-P(Col,Row)*porosity(Col,Row))/dx 
			dPxE=(P(Col+1,Row,Lay)-P(Col,Row,Lay))/dx
		    uE=uMid-APUMid*dx*dy*(dPxE-dPxMid)
		    FlowX(Col,Row,Lay)=density*dy*uE
		    AE(Col,Row,Lay)=-density*dy*APUMid*dy
		    AW(Col+1,Row,Lay)=AE(Col,Row,Lay)
	     EndDo
      EndDo
      Do Row=1,RowTotal-1
         Do Col=1,ColTotal
	        dPyMid=(dPy(Col,Row+1,Lay)+dPy(Col,Row,Lay))/2.0
		    vMid=(v(Col,Row+1,Lay)*porosity(Col,Row+1,Lay)+v(Col,Row,Lay)*porosity(Col,Row,Lay))/2.0 
		    APVMid=(APV(Col,Row+1,Lay)+APV(Col,Row,Lay))/2.0
		    !dPyN=(P(Col,Row+1)*porosity(Col,Row+1)-P(Col,Row)*porosity(Col,Row))/dy 
			dPyN=(P(Col,Row+1,Lay)-P(Col,Row,Lay))/dy
		    vN=vMid-APVMid*dx*dy*(dPyN-dPyMid)
		    FlowY(Col,Row,Lay)=density*dx*vN
		    AN(Col,Row,Lay)=-density*dx*APVMid*dx
		    AS(Col,Row+1,Lay)=AN(Col,Row,Lay)
	     EndDo
      EndDo
      !Source Term
      MassSum=0.0
      Do Col=1,ColTotal
         Do Row=1,RowTotal
	        SU(Col,Row,Lay)=FlowX(Col-1,Row,Lay)-FlowX(Col,Row,Lay)+FlowY(Col,Row-1,Lay)-FlowY(Col,Row,Lay)
		    AP(Col,Row,Lay)=-(AE(Col,Row,Lay)+AW(Col,Row,Lay)+AN(Col,Row,Lay)+AS(Col,Row,Lay))
		    MassSum=MassSum+SU(Col,Row,Lay)
		    Pr(Col,Row,Lay)=0.0
	     EndDo
      EndDo
      Call StoneSIP(Pr,Lay)
      Do Col=1,ColTotal
         Pr(Col,0,Lay)=Pr(Col,1,Lay)+(Pr(Col,1,Lay)-Pr(Col,2,Lay))/2.0
	     Pr(Col,RowTotal+1,Lay)=Pr(Col,RowTotal,Lay)+(Pr(Col,RowTotal,Lay)-Pr(Col,RowTotal-1,Lay))/2.0
      EndDo
      Do Row=1,RowTotal
         Pr(0,Row,Lay)=Pr(1,Row,Lay)+(Pr(1,Row,Lay)-Pr(2,Row,Lay))/2.0
	     Pr(ColTotal+1,Row,Lay)=Pr(ColTotal,Row,Lay)+(Pr(ColTotal,Row,Lay)-Pr(ColTotal-1,Row,Lay))/2.0
      EndDo
      Do Col=1,ColTotal-1
         Do Row=1,RowTotal
	        FlowX(Col,Row,Lay)=FlowX(Col,Row,Lay)+AE(Col,Row,Lay)*(Pr(Col+1,Row,Lay)-Pr(Col,Row,Lay))
	     EndDo
      EndDo
      Do Col=1,ColTotal
         Do Row=1,RowTotal-1
	        FlowY(Col,Row,Lay)=FlowY(Col,Row,Lay)+AN(Col,Row,Lay)*(Pr(Col,Row+1,Lay)-Pr(Col,Row,Lay))
	     EndDo
      EndDo
      Do Col=1,ColTotal
         Do Row=1,RowTotal
	        PrE=Pr(Col+1,Row,Lay)*factorX(Col)+Pr(Col,Row,Lay)*(1.0-factorX(Col))
		    PrW=Pr(Col,Row,Lay)*factorX(Col-1)+Pr(Col-1,Row,Lay)*(1.0-factorX(Col-1))
		    PrN=Pr(Col,Row+1,Lay)*factorY(Row)+Pr(Col,Row,Lay)*(1.0-factorY(Row))
		    PrS=Pr(Col,Row,Lay)*factorY(Row-1)+Pr(Col,Row-1,Lay)*(1.0-factorY(Row-1))
		    u(Col,Row,Lay)=u(Col,Row,Lay)-(PrE-PrW)*dy*APU(Col,Row,Lay)
		    v(Col,Row,Lay)=v(Col,Row,Lay)-(PrN-PrS)*dx*APV(Col,Row,Lay)
		    P(Col,Row,Lay)=P(Col,Row,Lay)+relax2*(Pr(Col,Row,Lay)-Pr(Col,RowTotal,Lay))
	     EndDo
      EndDo
   EndDo
EndDo
!DragForceX=0.0
!DragForceY=0.0
!EndIf
   If (MOD(INT(time/t),INT(timeframe/t))==0) then
      Write(1,*) 'ZONE'
	  !Write(2,*) time !'ZONE'
	  !Write(3,*) time !'ZONE'
	  !Write(4,*) time !'ZONE'
	  Write(6,*) time !'ZONE'
	  !Write(7,*) time !'ZONE'
	  !Write(8,*) time !'ZONE'
	  !Write(9,*) time !'ZONE'
	  !Write(10,*) time !'ZONE'
	  Write(11,*) 'ZONE'
	  Write(12,*) 'ZONE'
	  Write(13,*) time
	  Print *, 'Wee Chuan''s Program running...'
	  Print *, 'Time= ',time, 's'
      Do i=1,N
	     Write(1,"(7E12.4)") Particle(i)%x,Particle(i)%y,Particle(i)%z,Particle(i)%group,Particle(i)%VelX,Particle(i)%VelY,Particle(i)%VelZ
		 Write(12,"(7E12.4)") Particle(i)%x,Particle(i)%y,Particle(i)%z,Particle(i)%group,Particle(i)%ParticleF,Particle(i)%drag,Particle(i)%Fcap
	  EndDo
!	  Do Col=1,ColTotal
!	     velYtotal=0.0
!	     velXtotal=0.0
!	     ugastotal=0.0
!	     vgastotal=0.0
!	     count=0.0
!		 xpos=(2.0*Col-1)*Diameter/(2.0*ColTotal)
!		 Do Row=1,RowTotal
!		    Do k=1,MaxParticle
!			   If (Cell(Row,Col,k)/=0) then
!			      count=count+1.0
!				  i=Cell(Row,Col,k)
!				  velYtotal=velYtotal+Particle(i)%VelY
!				  velXtotal=velXtotal+Particle(i)%VelX
!			   EndIf
!			EndDo
!			vgastotal=vgastotal+v(Col,Row)
!			ugastotal=ugastotal+u(Col,Row)
!		 EndDo
!		 If (count>=1.0) then
!		    velYavg=velYtotal/count
!		    velXavg=velXtotal/count
!		 EndIf
!		 vgasavg=vgastotal/REAL(RowTotal)
!		 ugasavg=ugastotal/REAL(RowTotal)
!        Write(2,*) xpos,velYavg,velXavg
!		 Write(3,*) xpos,vgasavg,ugasavg
!		 vs=(count*m*ColTotal)/(soliddensity*Length*Diameter*d)
!		 Write(4,*) xpos,vs
!	  EndDo
	  !If (time>gas_start) Write(5,*) time-gas_start,Gs*m/(time-gas_start)
	  Do Row=1,RowTotal
!	     velYtotal=0.0
!	     velXtotal=0.0
!	     ugastotal=0.0
!	     vgastotal=0.0
!		 vYtotal=0.0
!	     vXtotal=0.0
	     count=0.0
	     Pavg=0.0
	     ypos=(2.0*Row-1)*Length/(2.0*RowTotal)
		 Do Col=1,ColTotal
		    Pavg=Pavg+P(Col,Row,LayTotal/2)/REAL(ColTotal)
!		    Do k=1,MaxParticle
!			   If (Cell(Row,Col,LayTotal/2,k)/=0) then
!			      count=count+1.0
!				  i=Cell(Row,Col,LayTotal/2,k)
!				  velYtotal=velYtotal+Particle(i)%VelY
!				  velXtotal=velXtotal+Particle(i)%VelX
!			   EndIf
!			EndDo
!			vgastotal=vgastotal+v(Col,Row,LayTotal/2)
!			ugastotal=ugastotal+u(Col,Row,LayTotal/2)
		 EndDo
		 Write(6,*) ypos,Pavg
!		 If (count>=1.0) then
!		    velYavg=velYtotal/count
!		    velXavg=velXtotal/count
!		 EndIf
!		 vgasavg=vgastotal/REAL(ColTotal)
!		 ugasavg=ugastotal/REAL(ColTotal)
!        Write(7,*) ypos,velYavg,velXavg
!		 Write(8,*) ypos,vgasavg,ugasavg
!		 vs=(count*m*RowTotal)/(soliddensity*Length*Diameter*d)
!		 Write(9,*) ypos,vs
!		 Do Col=1,ColTotal
!		    Do k=1,MaxParticle
!		       If (Cell(Row,Col,k)/=0) then
!			      i=Cell(Row,Col,k)
!			      vYtotal=vYtotal+(Particle(i)%VelY-velYavg)**2
!			      vXtotal=vXtotal+(Particle(i)%VelX-velXavg)**2
!			   EndIf
!			EndDo
!		 EndDo
!		 If (count>=1.0) then
!		    vYavg=0.5*vYtotal/count
!		    vXavg=0.5*vXtotal/count
!		 EndIf
!		 Write(10,*) ypos,vYavg,vXavg
         count=0.0
		 group1=0.0
         Do Col=1,ColTotal
		    Do Lay=1,LayTotal
			   Do k=1,MaxParticle
			      If (Cell(Row,Col,Lay,k)/=0) then
				     count=count+1.0
					 i=Cell(Row,Col,Lay,k)
					 If (Particle(i)%group==1.0) group1=group1+1.0
				  EndIf
			   EndDo
			EndDo
		 EndDo
		 If (count>=1.0) Write(13,*) ypos,(group1*soliddensity1)/((group1*soliddensity1)+(count-group1)*soliddensity2)
	  EndDo
	  Do Row=1,RowTotal
	     Do Col=1,ColTotal
		    xpos=(2.0*Col-1)*Diameter/(2.0*ColTotal)
			ypos=(2.0*Row-1)*Length/(2.0*RowTotal)
			Write(11,"(4E14.6)") xpos,ypos,u(Col,Row,LayTotal/2),v(Col,Row,LayTotal/2)
		 EndDo
	  EndDo
	  count=0.0
	  avgcount=0.0
	  group1=0.0
	  sigma2=0.0
	  sample=0.0
	  Do Lay=1,LayTotal-1,2
	     Do Col=1,ColTotal-1,2
		    Do Row=1,RowTotal-1,2
			   Do LayTemp2=Lay,Lay+1
			      Do ColTemp2=Col,Col+1
				     Do RowTemp2=Row,Row+1
			            Do k=1,MaxParticle
			               If (Cell(RowTemp2,ColTemp2,LayTemp2,k)/=0) then
				              count=count+1.0
							  i=Cell(RowTemp2,ColTemp2,LayTemp2,k)
					          If (Particle(i)%group==1.0) group1=group1+1.0
					       EndIf
						EndDo
					 EndDo
				  EndDo
			   EndDo
			   If (count>40.0) then
			      sigma2=sigma2+(group1/count-mu)**2
				  sample=sample+1.0
				  avgcount=avgcount+count
			   EndIf
			   count=0.0
			   group1=0.0
			EndDo
		 EndDo
	  EndDo
	  If (sample>20.0) then
	     sigma2=sigma2/sample
		 avgcount=avgcount/sample
		 sigmaR2=sigma02/avgcount
         Lacey=(sigma02-sigma2)/(sigma02-sigmaR2)
		 Write(14,"(4E14.6)") time,sample,avgcount,Lacey
	  EndIf
   EndIf
   time=time+t
EndDo
Close(1)
!Close(2)
!Close(3)
!Close(4)
!Close(5)
Close(6)
!Close(7)
!Close(8)
!Close(9)
!Close(10)
Close(11)
Close(12)
Close(13)
Close(14)
Stop 'CFD-DEM Simulation Completed.'
CONTAINS
Subroutine StoneSIP(F,Lay)
Implicit None
DoublePrecision, Parameter:: alpha=0.92
DoublePrecision:: P1, P2, ResTotal
DoublePrecision, Dimension(0:ColTotal+1,0:RowTotal+1,0:LayTotal+1):: F, LW, LS, LPR, UE, UN, residual
Integer:: Lay
UE=0.0
UN=0.0
residual=0.0
   Do Col=1,ColTotal
      Do Row=1,RowTotal
         LW(Col,Row,Lay)=AW(Col,Row,Lay)/(1.0+alpha*UN(Col-1,Row,Lay)+1.0e-30)
	     LS(Col,Row,Lay)=AS(Col,Row,Lay)/(1.0+alpha*UE(Col,Row-1,Lay)+1.0e-30)
	     P1=alpha*LW(Col,Row,Lay)*UN(Col-1,Row,Lay)
	     P2=alpha*LS(Col,Row,Lay)*UE(Col,Row-1,Lay)
	     LPR(Col,Row,Lay)=1.0/(AP(Col,Row,Lay)+P1+P2-LW(Col,Row,Lay)*UE(Col-1,Row,Lay)-LS(Col,Row,Lay)*UN(Col,Row-1,Lay)+1.0e-30)
	     UN(Col,Row,Lay)=(AN(Col,Row,Lay)-P1)*LPR(Col,Row,Lay)
	     UE(Col,Row,Lay)=(AE(Col,Row,Lay)-P2)*LPR(Col,Row,Lay)
      EndDo
   EndDo
   ResTotal=1.0
   Do While (ResTotal>0.0001)
      ResTotal=0.0
	  Do Col=1,ColTotal
	     Do Row=1,RowTotal
		    residual(Col,Row,Lay)=SU(Col,Row,Lay)-AN(Col,Row,Lay)*F(Col,Row+1,Lay)-AS(Col,Row,Lay)*F(Col,Row-1,Lay)- &
			AE(Col,Row,Lay)*F(Col+1,Row,Lay)-AW(Col,Row,Lay)*F(Col-1,Row,Lay)-AP(Col,Row,Lay)*F(Col,Row,Lay)
			ResTotal=ResTotal+ABS(residual(Col,Row,Lay))
			residual(Col,Row,Lay)=(residual(Col,Row,Lay)-LS(Col,Row,Lay)*residual(Col,Row-1,Lay)- &
			LW(Col,Row,Lay)*residual(Col-1,Row,Lay))*LPR(Col,Row,Lay)
		 EndDo
	  EndDo
	  !Back Substitution and Correction
	  Do Col=ColTotal,1,-1
	     Do Row=RowTotal,1,-1
		    residual(Col,Row,Lay)=residual(Col,Row,Lay)-UN(Col,Row,Lay)*residual(Col,Row+1,Lay)-UE(Col,Row,Lay)*residual(Col+1,Row,Lay)
			F(Col,Row,Lay)=F(Col,Row,Lay)+residual(Col,Row,Lay)
		 EndDo
	  EndDo
   EndDo
End Subroutine
End Program
