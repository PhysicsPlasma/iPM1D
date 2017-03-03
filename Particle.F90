Module ParticleModule
      ! This module is the operations for the particles.
      Use TypeModule
      Use Constants 
     Implicit none
     Contains
             Subroutine  AddParticle(OneParticle,InputParticle)
             Implicit none
                 Type(ParticleOne),intent(inout) ::  OneParticle
                 Type(ParticleBundle),intent(inout) ::  InputParticle
                 InputParticle%NPar=InputParticle%NPar+1
                 InputParticle%PhaseSpace(InputParticle%NPar)=OneParticle
             return   
             End  Subroutine  AddParticle
     
             Subroutine  DelParticle(NDel,InputParticle)
             Implicit none
                 integer(4),intent(in) ::  NDel
                 Type(ParticleBundle),intent(inout) ::  InputParticle
                 InputParticle%PhaseSpace(NDel)=InputParticle%PhaseSpace(InputParticle%NPar)
                 InputParticle%NPar=InputParticle%NPar-1
             return   
             End  Subroutine  DelParticle
             
            Subroutine  UpdatePositionBundle(InputParticle,Sign)
                 Implicit none
                    Type(ParticleBundle), intent(inout) ::  InputParticle
                    Integer(4),intent(in) :: Sign
                    Real(8) :: XFactor
                    Integer(4) :: i
                    Select case (Sign)
                           Case(1_4)
                                    XFactor=InputParticle%XFactor
                          Case(-1_4)
                                    XFactor=1.d0/InputParticle%XFactor
                    End Select       
                   do i=1,InputParticle%NPar
                         InputParticle%PhaseSpace(i)%X=XFactor*InputParticle%PhaseSpace(i)%X
                    end do 
                    return 
            End Subroutine  UpdatePositionBundle
             
            Subroutine  UpdateVelocityBundle(InputParticle,Sign)
                 Implicit none
                    Type(ParticleBundle), intent(inout) ::  InputParticle
                    Integer(4),intent(in) :: Sign
                    Real(8) :: VFactor
                    Integer(4) :: i
                    Select case (Sign)
                           Case(1_4)
                                    VFactor=InputParticle%VFactor
                            Case(-1_4)
                                    VFactor=1.d0/InputParticle%VFactor
                         End Select 
                     do i=1,InputParticle%NPar
                                     InputParticle%PhaseSpace(i)%Vx=VFactor*InputParticle%PhaseSpace(i)%Vx
                                     InputParticle%PhaseSpace(i)%Vy=VFactor*InputParticle%PhaseSpace(i)%Vy
                                     InputParticle%PhaseSpace(i)%Vz=VFactor*InputParticle%PhaseSpace(i)%Vz
                     end do 
                    return 
            End Subroutine  UpdateVelocityBundle
            
      subroutine AddParticleBundle(SmallParticle,MainParticle)
            Implicit none
            Type(ParticleBundle),intent(inout) :: MainParticle
            Type(ParticleBundle),intent(inout) :: SmallParticle
            Integer(4) :: i
            Do i=1,SmallParticle%NPar
                  Call  AddParticle(SmallParticle%PhaseSpace(i),MainParticle)
            End do
            return 
        end subroutine AddParticleBundle
        
        Subroutine PhaseSpaceNormalization(PB) 
              Implicit none
              Type(ParticleBundle),intent(inout) :: PB
              Type(ParticleOne) :: ParticleTemp
              Integer(4) :: i,Nold,Nnew,NMax,NMin
              Integer(4) :: Ndiff,Index
              NOld=PB%NPar
              NMin=0.4*PB%NParNormal
              NMax=2*PB%NParNormal
              Write(*,*) "PhaseSpaceNormalizationA",PB%NPar,Nmin,NMax
              If(NOld<NMin) then
                 NNew=2*NOld
                 Ndiff=NNew-NOld
                 do  i=1,Ndiff
                       Call DRandom(R)
                       Index=Ceiling(R*NOld)
                       ParticleTemp=PB%PhaseSpace(i)
                       Call AddParticle(ParticleTemp,PB)
                 end do
              else if (NOld>NMax) then
                 NNew=0.5*NOld 
                 Ndiff=NOld-NNew
                 do i=1,NDiff
                       Call DRandom(R)
                       Index=Ceiling(R*PB%NPar)
                       Call DelParticle(Index,PB)
                 end do
             end if
             PB%Weight=PB%Weight*dble(NOld)/dble(NNew)
             Write(*,*) "PhaseSpaceNormalizationB",PB%NPar,NOld
             return 
        end  Subroutine  PhaseSpaceNormalization
        
        Subroutine PhaseSpaceNormalization2(NParNew,PB) 
              Implicit none
              Integer(4),intent(in) :: NParNew
              Type(ParticleBundle),intent(inout) :: PB
              !Type(PICParticlePara),intent(inout) :: OutputParticlePara
              Type(ParticleOne) :: ParticleTemp
              Integer(4) :: i,NTemp,Ndiff,Index
              NTemp=PB%NPar
              PB%Weight=PB%Weight*dble(NTemp)/dble(NParNew)
              If(NParNew>NTemp) then
                 Ndiff=NParNew-NTemp
                 do  i=1,Ndiff
                       Call DRandom(R)
                       Index=Ceiling(R*NTemp)
                       ParticleTemp=PB%PhaseSpace(i)
                       Call AddParticle(ParticleTemp,PB)
                 end do
             else
                 Ndiff=NTemp-NParNew
                 do i=1,NDiff
                       Call DRandom(R)
                       Index=Ceiling(R*PB%NPar)
                       Call DelParticle(Index,PB)
                 end do
             end if  
             return 
           end  Subroutine  PhaseSpaceNormalization2
          
           Subroutine CalEnergy(Mass,InputParticle,Energy) 
              Implicit none
              Real(8),intent(in):: Mass 
              Type(ParticleOne),intent(in) :: InputParticle
              Real(8),intent(inout) ::Energy
              Energy=0.5d0*Mass*(InputParticle%Vx*InputParticle%Vx+InputParticle%Vy*InputParticle%Vy +InputParticle%Vz*InputParticle%Vz) 
             return 
           end  Subroutine  CalEnergy
           
           
           !Subroutine ParticleBundleCopy(PBT,PB) 
           !   Implicit none
           !   Type(ParticleBundleTimer),intent(inout) :: PBT
           !   Type(ParticleBundle),intent(in) :: PB
           !    PBT%SpecyIndex=PB%SpecyIndex
           !    PBT%Timer=0
           !    PBT%Period=PB%Period
           !    PBT%NPar=0
           !    PBT%NParNormal=PB%NParNormal
           !    PBT%Charge=PB%Charge
           !    PBT%Weight=PB%Weight
           !    PBT%XFactor=PB%XFactor
           !    PBT%VFactor=PB%VFactor
           !    PBT%dx=PB%dx
           !    PBT%dt=PB%dt
           !    PBT%XMin=PB%XMin
           !    PBT%XMax=PB%XMax
           !  return 
           !end  Subroutine  ParticleBundleCopy
           
End  Module ParticleModule

