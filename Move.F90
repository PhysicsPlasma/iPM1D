Module MoveModule
     Use TypeModule
     Use ParticleModule
     Implicit none
     !Real(8),parameter,Private :: Theta=1.d0
     !Real(8),parameter,Private :: FactorA=0.5d0*Theta
     !Real(8),parameter,Private :: FactorB=1.d0-0.5d0*Theta
     contains
            Subroutine   ParticleMove(PB,FG)
               Type(ParticleBundle),intent(inout) :: PB
               Type(Field),intent(in) :: FG
               Real(8) :: Ex(FG%Nx),Ey(FG%Nx),Bx(FG%Nx),By(FG%Nx)
               Integer(4) :: i,N
               Real(8) :: EFactor,BFactor
               Real(8) :: S1,S2
               Real(8) :: E1,E2,E3,B1,B2,B3
               Real(8):: SVelocity(1:3)=0.d0,Ep(1:3)=0.d0,Bp(1:3)=0.d0
               Real(8) :: Omega(1:3)=0.d0,OmegaT,TransB(3,3)=0.d0
               EFactor=PB%Charge/PB%Mass*FG%dt/(FG%dx/FG%dt)
               BFactor=PB%Charge/PB%Mass*FG%dt
               Ex=FG%Ex*EFactor
               Ey=FG%Ey*EFactor
               Bx=FG%Bx*BFactor
               By=FG%By*BFactor

               do i=1,PB%NPar
                      
                       N=Ceiling(PB%PhaseSpace(i)%X)
                       S1=Dble(N)-PB%PhaseSpace(i)%X
                       S2=1.d0-S1

                       Ep(1)=Ex(N)*S1+Ex(N+1)*S2
                       Ep(2)=Ey(N)*S1+Ey(N+1)*S2
                       Ep(3)=0.0
                       
                       Bp(1)=Bx(N)*S1+Bx(N+1)*S2
                       Bp(2)=By(N)*S1+By(N+1)*S2
                       Bp(3)=0.0
                       
                       Omega(1)=Bp(1)*0.5d0
                       Omega(2)=Bp(2)*0.5d0
                       Omega(3)=Bp(3)*0.5d0
                       OmegaT=1.d0+Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3)
                       
                       TransB(1,1)=1.d0+Omega(1)*Omega(1)
                       TransB(2,1)=Omega(1)*Omega(2)+Omega(3)
                       TransB(3,1)=Omega(1)*Omega(3)-Omega(2)
                       TransB(1,2)=Omega(1)*Omega(2)-Omega(3)
                       TransB(2,2)=1.d0+Omega(2)*Omega(2)
                       TransB(3,2)=Omega(2)*Omega(3)+Omega(1)
                       TransB(1,3)=Omega(1)*Omega(3)+Omega(2)
                       TransB(2,3)=Omega(2)*Omega(3)-Omega(1)
                       TransB(3,3)=1.d0+Omega(3)*Omega(3)
                       TransB(1:3,1:3)=TransB(1:3,1:3)/OmegaT

                       SVelocity(1)=(TransB(1,1)*Ep(1)+TransB(2,1)*Ep(2)+TransB(3,1)*Ep(3))*0.5d0
                       SVelocity(2)=(TransB(1,2)*Ep(1)+TransB(2,2)*Ep(2)+TransB(3,2)*Ep(3))*0.5d0
                       SVelocity(3)=(TransB(1,3)*Ep(1)+TransB(2,3)*Ep(2)+TransB(3,3)*Ep(3))*0.5d0

                       PB%PhaseSpace(i)%Vx=PB%PhaseSpace(i)%Vx+SVelocity(1)
                       PB%PhaseSpace(i)%Vy=PB%PhaseSpace(i)%Vy+SVelocity(2)
                       PB%PhaseSpace(i)%Vz=PB%PhaseSpace(i)%Vz+SVelocity(3)
                       
                       PB%PhaseSpace(i)%X=PB%PhaseSpace(i)%X+SVelocity(1)

                       PB%PhaseSpace(i)%Ax=0.5d0*(PB%PhaseSpace(i)%Ax+Ep(1)) 
                       PB%PhaseSpace(i)%Ay=0.5d0*(PB%PhaseSpace(i)%Ay+Ep(2))
                       PB%PhaseSpace(i)%Az=0.5d0*(PB%PhaseSpace(i)%Az+Ep(3))
                      
                       SVelocity(1)=PB%PhaseSpace(i)%Vx+0.5d0*PB%PhaseSpace(i)%Ax+&
                         PB%PhaseSpace(i)%Vy*Omega(3)-PB%PhaseSpace(i)%Vz*Omega(2)
                       SVelocity(2)=PB%PhaseSpace(i)%Vy+0.5d0*PB%PhaseSpace(i)%Ay+&
                         PB%PhaseSpace(i)%Vz*Omega(1)-PB%PhaseSpace(i)%Vx*Omega(3)
                       SVelocity(3)=PB%PhaseSpace(i)%Vz+0.5d0*PB%PhaseSpace(i)%Az+&
                         PB%PhaseSpace(i)%Vx*Omega(2)-PB%PhaseSpace(i)%Vy*Omega(1)
                       
                       PB%PhaseSpace(i)%Vx=TransB(1,1)*SVelocity(1)+TransB(2,1)*SVelocity(2)+TransB(3,1)*SVelocity(3)
                       PB%PhaseSpace(i)%Vy=TransB(1,2)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(3,2)*SVelocity(3)
                       PB%PhaseSpace(i)%Vz=TransB(1,3)*SVelocity(1)+TransB(2,3)*SVelocity(2)+TransB(3,3)*SVelocity(3)
        !---------------------定义了一个速度，是Vy和Vz方向上的速度平方根
                       
                       PB%PhaseSpace(i)%X=PB%PhaseSpace(i)%X+PB%PhaseSpace(i)%Vx
                      
               End do
               PB%Timer=PB%Timer+1
              return
             End Subroutine  ParticleMove
    End Module MoveModule
    
    
Module ParticleBoundaryModule
  Use TypeModule
  Use Constants
  Use ParticleModule
  Implicit none
       Type ParticleBoundary
                        !Integer(4) :: XStart=0,XEnd=1
                         Integer(4) :: ParticleBoundaryModel=12
                         Real(8) :: XMin=0.d0,XMax=dble(NxMax-1)!,XLength=99.d0
                         
                         !Integer(4) :: 
                         !Real(8) :: 
                         
       EndType ParticleBoundary
       
       Type ParticleBoundaryOne
                        !Integer(4) :: XStart=0,XEnd=1
                         Integer(4) :: ParticleBoundaryModel=11
                         Real(8) :: XMin=0.d0,XMax=dble(NxMax-1)!,XLength=99.d0 
                         Integer :: CountMinOne=0,CountmaxOne=0
                         Integer :: CountMin=0,Countmax=0
                         Real(8) :: Gamma=0.2d0
                         Integer(4) :: SEENparMin,SEENparMax
                         Type(ParticleBundle) :: PBLower,PBUpper
       EndType ParticleBoundaryOne
       
   contains 
   subroutine ParticleAborption(PB,PBDO)
        implicit none
        Type(ParticleBundle),intent(inout) :: PB
        Type(ParticleBoundaryOne),intent(inout) :: PBDO
        Integer :: i
        Select case (PBDO%ParticleBoundaryModel)
            case(10)
                do i=PB%NPar,1,-1
                     If (PB%PhaseSpace(i)%X<=PBDO%XMin.or.PB%PhaseSpace(i)%X>=PBDO%XMax) then
                        !Write(*,*) "Xmin"
                         Call DelParticle(i,PB)
                     end if
                end do
            case(11)
                PBDO%CountMinOne=0
                PBDO%CountMaxOne=0
                do i=PB%NPar,1,-1
                     If (PB%PhaseSpace(i)%X<=PBDO%XMin) then
                         PBDO%CountMinOne=PBDO%CountMinOne+1
                         !Write(*,*) "Xmin"
                         Call DelParticle(i,PB)
                     else If(PB%PhaseSpace(i)%X>=PBDO%XMax) then
                         PBDO%CountMaxOne=PBDO%CountMaxOne+1
                         !Write(*,*) "Xmax"
                         Call DelParticle(i,PB)                         
                     end if
                end do
              PBDO%CountMin=PBDO%CountMin+PBDO%CountMinOne
              PBDO%CountMax=PBDO%CountMax+PBDO%CountMaxOne
            case(12)
                PBDO%CountMinOne=0
                PBDO%CountMaxOne=0
                PBDO%PBLower%NPar=0
                PBDO%PBUpper%NPar=0
                do i=PB%NPar,1,-1
                     If (PB%PhaseSpace(i)%X<=PBDO%XMin) then
                         PBDO%CountMinOne=PBDO%CountMinOne+1
                         Call AddParticle(PB%PhaseSpace(i),PBDO%PBLower)
                         Call DelParticle(i,PB)
                     else If(PB%PhaseSpace(i)%X>=PBDO%XMax) then
                         PBDO%CountMaxOne=PBDO%CountMaxOne+1
                         Call AddParticle(PB%PhaseSpace(i),PBDO%PBUpper)
                         Call DelParticle(i,PB)                         
                     end if
                end do
                PBDO%CountMin=PBDO%CountMin+PBDO%CountMinOne
                PBDO%CountMax=PBDO%CountMax+PBDO%CountMaxOne
        End Select 
        return
  end subroutine ParticleAborption
  !
  !subroutine AborptionCount(PB,XMin,XMax,CountMin,CountMax)
  !      implicit none
  !      Type(ParticleBundle),intent(inout) :: PB
  !      Real(8),intent(in)  :: XMin,XMax
  !     ! Integer,intent(in) :: Sign
  !      Integer,intent(inout) :: CountMin,CountMax
  !      Integer :: i
  !      CountMin=0
  !      CountMax=0
  !      do i=PB%NPar,1,-1
  !         If (PB%PhaseSpace(i)%X<=XMin) then
  !            CountMin=CountMin+1
  !            Call DelParticle(i,PB)
  !         else If (PB%PhaseSpace(i)%X>=XMax) then
  !             CountMax=CountMax+1
  !            Call DelParticle(i,PB)
  !         end if
  !      end do
  !      return
  ! end subroutine AborptionCount
    end  Module ParticleBoundaryModule
    
Module SEEModule
    Use TypeModule
    Use ParticleBoundaryModule
    Use ParticleModule
    Use MCCModule
    Use Constants
    Implicit none
  
     contains
 Subroutine  Selectron(PB,PBDO)
               Implicit none
               Type(ParticleBundle),intent(inout) :: PB
               Type(ParticleBoundaryOne),intent(inout) :: PBDO
               Type(ParticleOne) :: ParticleTemp
               Integer(4) :: i
               Type(Gas),parameter :: SEEGas=(Gas(1,0,9.1095d-31,0.026d0*11605.d0,0.d0,0.d0)) 
               Real(8) :: XRandom=0.01,VFactor,Residue
               
               VFactor=1.d0/PB%VFactor
                
               !If (PBDO%CountMinOne/=0.or.PBDO%CountMaxOne/=0) Write(*,*) "SEE!"
               
               PBDO%SEENparMin=Int(PBDO%CountMinOne*PBDO%Gamma)
               Residue=PBDO%CountMinOne*PBDO%Gamma-Dble(PBDO%SEENparMin)
               Call DRandom(R)
               If (R<Residue) Then
                   PBDO%SEENparMin=PBDO%SEENparMin+1
                   !Write(*,*) "SEE2!"
               End If
               do i=1,PBDO%SEENparMin
                  Call DRandom(R)
                  ParticleTemp%X=PB%XMin+R*XRandom
                  ParticleTemp%Ax=0.0
                  ParticleTemp%Ay=0.0
                  ParticleTemp%Az=0.0
                Call Maxwellian(SEEGas,ParticleTemp)
                ParticleTemp%Vx=VFactor*ParticleTemp%Vx
                ParticleTemp%Vy=VFactor*ParticleTemp%Vy
                ParticleTemp%Vz=VFactor*ParticleTemp%Vz
                ParticleTemp%Vx=abs(ParticleTemp%Vx)
                Call AddParticle(ParticleTemp,PB)
               enddo
               
               PBDO%SEENparMax=Int(PBDO%CountMaxOne*PBDO%Gamma)
               Residue=PBDO%CountMaxOne*PBDO%Gamma-Dble(PBDO%SEENparMax)
               Call DRandom(R)
               If (R<Residue) Then 
                   PBDO%SEENparMax=PBDO%SEENparMax+1
                    !Write(*,*) "SEE3!"
               End If
               do i=1,PBDO%SEENparMax
                  Call DRandom(R)
                  ParticleTemp%X=PB%XMax-R*XRandom
                  ParticleTemp%Ax=0.0
                  ParticleTemp%Ay=0.0
                  ParticleTemp%Az=0.0
                Call Maxwellian(SEEGas,ParticleTemp)
                ParticleTemp%Vx=VFactor*ParticleTemp%Vx
                ParticleTemp%Vy=VFactor*ParticleTemp%Vy
                ParticleTemp%Vz=VFactor*ParticleTemp%Vz
                ParticleTemp%Vx=-abs(ParticleTemp%Vx)
                Call AddParticle(ParticleTemp,PB)
               enddo
            return
       end subroutine Selectron
  End Module SEEModule     