Module MCCPublic
     Use TypeModule
     Use Constants
     Use ParticleModule
     Use  MCCModule
     Use ElectronMCCModule
    Use  IonMCC
     Implicit none
     Integer(4),parameter :: NSpecyMCCTemp=4_4
     Type(MCCParticleOne),private :: TempMCCParticle
     Type(ParticleBundle),private ::  MCCBefore,MCCAfter(0:NSpecyMCCTemp)

   Contains
   subroutine  MCCBundleInit(InputGasPhysical,InputGas,InputSpecy,InputReactionBundle,InputSigma,OutputMCCBundle,dt)
             Use SigmaModule
             Implicit none
             Type(GasPhysical),intent(in) :: InputGasPhysical
             Type(Gas),intent(in) ::  InputGas
             Type(OneSpecy),intent(in) :: InputSpecy(1:InputGas%NSpecy)
             Type(ReactionBundle),intent(in) :: InputReactionBundle(0:InputGas%NSpecy)
             Type(SigmaNormalized),intent(in) :: InputSigma(0:InputGas%NSpecy)
             Type(MCCBundle),intent(inout) ::  OutputMCCBundle(0:InputGas%NSpecy)
             Real(8),intent(in) :: dt 
             Integer(4) ::  i,j,NSpecy
             Real(8) :: Miu
    
             NSpecy=InputGasPhysical%NSpecy
             
             OutputMCCBundle(0)%Model=1_4
             OutputMCCBundle(0)%Mass=ElectronMass
             do i=1,NSpecy
                  OutputMCCBundle(i)%Mass=InputSpecy(i)%Mass
                  OutputMCCBundle(i)%Model=InputGasPhysical%Model
             End do
             
             do i=0,NSpecy
                  OutputMCCBundle(i)%EnergyInterval=InputSigma(i)%EnergyInterval
                  OutputMCCBundle(i)%NReaction=InputSigma(i)%NReaction
                  OutputMCCBundle(i)%NSigma=InputSigma(i)%NSigma
                  OutputMCCBundle(i)%Reaction(0)%Reactant=InputReactionBundle(i)%Reaction(1)%Reactant
                  do j=1,InputReactionBundle(i)%NReaction
                        OutputMCCBundle(i)%Reaction(j)=InputReactionBundle(i)%Reaction(j)
                  End do
                  OutputMCCBundle(i)%Ng=(InputGas%PGas)/(kB*InputGas%TGas)
             End do

             Call  ProbilityNonReactive(InputSigma(0)%NReaction,InputSigma(0)%NSigma,InputSigma(0)%Value,OutputMCCBundle(0)%Probility,OutputMCCBundle(0)%EnergyInterval,ElectronMass,OutputMCCBundle(0)%SigmaMax)
             OutputMCCBundle(0)%CollisionRatio=1.d0-DExp(-OutputMCCBundle(0)%SigmaMax*OutputMCCBundle(0)%Ng*dt)
             
             do i=1, NSpecy
                      Select case (OutputMCCBundle(i)%Model)
                          Case(2_4)
                               Miu=InputSpecy(i)%Mass*InputGas%MGas/(InputGas%MGas+InputSpecy(i)%Mass)
                              Call  ProbilityNonReactive(InputSigma(i)%NReaction,InputSigma(i)%NSigma,InputSigma(i)%Value,OutputMCCBundle(i)%Probility,OutputMCCBundle(i)%EnergyInterval,Miu,OutputMCCBundle(0)%SigmaMax)
                              OutputMCCBundle(i)%CollisionRatio=1.d0-DExp(-OutputMCCBundle(0)%SigmaMax*OutputMCCBundle(i)%Ng*dt)
                              
                          Case(3_4)
                              Call ProbilityReactive(InputSigma(i)%NReaction,InputSigma(i)%NSigma,InputSigma(i)%Value,OutputMCCBundle(i)%Probility)
                              Call CollisionRatioReactive(InputSpecy(i),InputGasPhysical,OutputMCCBundle(i)%Ng,dt,OutputMCCBundle(i)%CollisionRatio)
                      End Select
               enddo 
              Write(*,*) OutputMCCBundle%CollisionRatio,'ratio'  
           return
         End subroutine MCCBundleInit 
         
          subroutine  ProbilityNonReactive(NReaction,NSigma,Sigma,Probility,EnergyInterval,Mass,SigmaT)
                implicit none
                Integer(4),intent(in) ::  NReaction,NSigma
                Real(8),intent(in) ::  Sigma(NSigma,NReaction),EnergyInterval,Mass
                Real(8),intent(inout) ::  Probility(NReaction,NSigma),SigmaT
                Real(8) ::Energy,V,Max
                Integer(4) :: i,j
                do i=1,NSigma
                     do j=1,NReaction
                          Probility(j,i)=Sigma(i,j)
                     Enddo
                Enddo
                   
                do i=1,NSigma
                      do j=2,NReaction
                            Probility(j,i)=Probility(j,i)+Probility(j-1,i)
                      end do
                end do
                
                Max=0.d0
                do i=1,NSigma
                      Energy=dble(i)*EnergyInterval
                      V=DSqrt(2.d0*Energy*JtoeV/Mass)  
                      do j=1,NReaction
                            Probility(j,i)=Probility(j,i)*V
                            If  (Probility(j,i)>Max) Then
                                  Max=Probility(j,i)
                            End if
                      end do
                end do
                
                do i=1,NSigma
                     do j=1,NReaction
                          Probility(j,i)=Probility(j,i)/Max
                     Enddo
                Enddo
                SigmaT=Max
             return 
          End subroutine  ProbilityNonReactive
          
         subroutine  CollisionRatioReactive(InputSpecy,InputGas,NGas,dt,CollisionRatio)
             implicit none
             Type(GasPhysical),intent(in) :: InputGas
             Type(OneSpecy),intent(in) :: InputSpecy
             Real(8),intent(in) :: NGas,dt
             Real(8),intent(inout) :: CollisionRatio
             Real(8) :: Miu,BetaMax,Radius
             Miu=InputSpecy%Mass*InputGas%MGas/(InputGas%MGas+InputSpecy%Mass)
             BetaMax=InputGas%BetaMax
             Radius=InputGas%Radius
             CollisionRatio=DSQRT(PI*Radius*a0*a0*a0*ElectronCharge*ElectronCharge/(Epsilon*Miu))*BetaMax*BetaMax*Ngas*dt
             return 
          End subroutine  CollisionRatioReactive
          
          subroutine  ProbilityReactive(NReaction,NSigma,Sigma,Probility)
                implicit none
                Integer(4),intent(in) ::  NReaction,NSigma
                Real(8),intent(in) ::  Sigma(NReaction,NSigma)
                Real(8),intent(inout) ::  Probility(NReaction,NSigma)
                Integer(4) :: i,j
                do i=1,NSigma
                     do j=1,NReaction
                          Probility(j,i)=Sigma(j,i)
                     Enddo
                Enddo
                
               do i=1,NSigma
                        do j=2,NReaction
                            Probility(j,i)=Probility(j,i)+Probility(j-1,i)
                         end do
                end do
              return
     End subroutine  ProbilityReactive

     Subroutine MCC(j,NSpecyTotal,MainParticle,InputGas,InputMCCBundle) 
        Implicit none
        ! j is the Particle Index;
        Integer(4),intent(in) :: j,NSpecyTotal
        Type(ParticleBundle),intent(inout) :: MainParticle(0:NSpecyTotal)
        Type(Gas),intent(in) :: InputGas
        Type(MCCBundle),intent(in) :: InputMCCBundle
        Integer(4) :: k,m,NSpecy,Index,Nbegin
        NSpecy=InputGas%NSpecy
        Nbegin=InputGas%Shift
        If(InputMCCBundle%Model/=0_4) then
               Call SelectParticle(MainParticle(j),MCCBefore,InputMCCBundle%CollisionRatio)
               MCCBefore%VFactor=MainParticle(j)%VFactor

               MCCAfter(0:NSpecy)%NPar=0
                         
               MCCAfter(0)%VFactor=MainParticle(0)%VFactor 
               do k=1,NSpecy
                       MCCAfter(k)%VFactor=MainParticle(Nbegin+k)%VFactor
               End do
               
               Call UpdateVelocityBundle(MCCBefore,1_4)
                
               do k=1,MCCBefore%NPar
                      select case (InputMCCBundle%Model)
                         Case(1_4)
                             Call UpdateParticleMCCElectron(MCCBefore%PhaseSpace(k),TempMCCParticle,InputMCCBundle%Mass)
                             Call SelectProbility(TempMCCParticle,InputMCCBundle)
                             Index=TempMCCParticle%Index
                             Call  SelectCollisionElectron(TempMCCParticle,InputGas,InputMCCBundle%Reaction(Index),MCCAfter(0:NSpecy))
                         Case(2_4)
                              Call UpdateParticleMCCIon(MCCBefore%PhaseSpace(k),TempMCCParticle,InputMCCBundle%Mass,InputGas) 
                              Call SelectProbility(TempMCCParticle,InputMCCBundle)
                              Index=TempMCCParticle%Index
                              Call  SelectCollisionIon(TempMCCParticle,InputGas,InputMCCBundle%Reaction(Index),MCCAfter(0:NSpecy))
                         Case(3_4)
                              Call UpdateParticleMCCIon(MCCBefore%PhaseSpace(k),TempMCCParticle,InputMCCBundle%Mass,InputGas)   
                              Call  SelectProbilityReactive(TempMCCParticle,InputMCCBundle,InputGas)
                              Index=TempMCCParticle%Index
                              Call  SelectCollisionIon(TempMCCParticle,InputGas,InputMCCBundle%Reaction(Index),MCCAfter(0:NSpecy))
                         end select    
               End do

              do k=0, NSpecy
                      Call UpdateVelocityBundle(MCCAfter(k),-1_4)
                      !do m=1,MCCAfter(k)%NPar
                                !MCCAfter(k)%PhaseSpace(m)%A1=0.d0
                               !MCCAfter(k)%PhaseSpace(m)%A2=0.d0
                              !MCCAfter(k)%PhaseSpace(m)%Az=0.d0 
                      !End do
                      
               End do 
               
               Call AddParticleBundle(MCCAfter(0),MainParticle(0))
               do k=1,NSpecy
                       Call AddParticleBundle(MCCAfter(k),MainParticle(Nbegin+k))
               end do
           Endif      
        return
     End  Subroutine MCC

      subroutine SelectProbility(InputMCCParticle,InputMCCBundle)
           Implicit None 
           Type(MCCParticleOne),intent(inout) :: InputMCCParticle
           Type(MCCBundle),intent(in) :: InputMCCBundle
           Integer(4) :: i,N,Index,Upper,Center,Lower,NReaction
           Real(8) :: EnergyInterval,Energy,S1,S2
           Real(8) :: TempProbility(InputMCCBundle%NReaction)
            TempProbility=0.d0 
            EnergyInterval=InputMCCBundle%EnergyInterval
            Energy=InputMCCParticle%Energy
            NReaction=InputMCCBundle%NReaction
            N=Int(Energy/EnergyInterval)
            If(N<2) then
                 TempProbility=InputMCCBundle%Probility(1:NReaction)
            else if (N<InputMCCBundle%NSigma) then
                 Index=N*NReaction
                 do i=1,NReaction
                       Center=Index+i
                       Upper=Center+NReaction
                       Lower=Center-NReaction
                       If(Energy<=InputMCCBundle%Reaction(i)%Threshold) Then
                           TempProbility(i:NReaction)=-1.d0
                           exit       
                      else
                           S1=(Energy-(dble(N)*EnergyInterval))/EnergyInterval
                           S2=1.d0-S1 
                           TempProbility(i)=S1*InputMCCBundle%Probility(Lower)+S2*InputMCCBundle%Probility(Upper)
                      end If
                  end do
           else
                  Index=(InputMCCBundle%NSigma-1)*NReaction 
                  TempProbility=InputMCCBundle%Probility(Index+1:Index+NReaction)
           End if
           Call DRandom(R)
           InputMCCParticle%Index=0 
           do i=1, NReaction
               If(TempProbility(i)<MinReal) then  
                   InputMCCParticle%Index=0
                   exit 
              else if(R<TempProbility(i)) then
                   InputMCCParticle%Index=i
                  exit
               End If
           end do 
  End  subroutine SelectProbility

   subroutine SelectProbilityReactive(InputMCCParticle,InputMCCBundle,InputGas)
          Implicit none
          Type(MCCParticleOne),intent(inout) :: InputMCCParticle
          Type(MCCBundle),intent(in) :: InputMCCBundle
          Type(Gas),intent(in) :: InputGas 
          Call DRandom(R)
                  !InputMCCParticle%Beta=InputGas%BetaMax*Dsqrt(dble(R))
                  If (InputMCCParticle%Beta<=1.d0) then   !  Reactive collisions
                           If(InputMCCParticle%Energy<=InputMCCBundle%Reaction(1)%Threshold) then
                                !  Reactive elastic collision with isotropic scattering
                                InputMCCParticle%Index=101_4
                           Else
                                Call SelectProbility(InputMCCParticle,InputMCCBundle)
                           End If
                 Else        ! Nonreactive Elastic collision with anisotropic scattering
                          InputMCCParticle%Index=102_4
                 End If
           return
        end subroutine  SelectProbilityReactive 
   
         subroutine SelectParticle(InputParticle,OuputParticle,CollisionRatio)
                    Implicit none
                    Type(ParticleBundle),intent(inout) :: InputParticle
                    Type(ParticleBundle),intent(inout) :: OuputParticle
                    Real(8),intent(in) :: CollisionRatio
                    Integer(4) :: i,N,NCollision
                    NCollision=Int(InputParticle%NPar*CollisionRatio)
                    OuputParticle%NPar=0 
                    do i=1, NCollision
                          Call DRandom(R)
                          N=Ceiling(R*InputParticle%NPar)
                          Call AddParticle(InputParticle%PhaseSpace(N),OuputParticle)
                          Call DelParticle(N,InputParticle)
                    end do
                    return
         end subroutine   SelectParticle  
           
        Subroutine UpdateParticleMCCElectron(InputParticle,OutputParticle,Mass)
                     Implicit none
                     Type(ParticleOne),intent(in) :: InputParticle
                     Type(MCCParticleOne),intent(inout) :: OutputParticle
                     Real(8),intent(in) :: Mass
                     Real(8) :: Vx,Vy,Vz,V2,V,Energy
                     OutputParticle%PhaseSpace=InputParticle
                     OutputParticle%Mass=Mass
                     Vx=InputParticle%Vx
                     Vy=InputParticle%Vy
                     Vz=InputParticle%Vz
                     V2=Vx*Vx+Vy*Vy+Vz*Vz
                     V=Dsqrt(V2)
                     Energy=0.5d0*Mass*V2/JtoeV
                     OutputParticle%V2=V2
                     OutputParticle%V=V
                     OutputParticle%Energy=Energy
                    return  
            End subroutine UpdateParticleMCCElectron
                
            Subroutine UpdateParticleMCCIon(InputParticle,OutputParticle,Mass,InputGas)
                     Implicit none
                     Type(ParticleOne),intent(in) :: InputParticle
                     Type(MCCParticleOne),intent(inout) :: OutputParticle
                     Real(8),intent(in) :: Mass
                     Type(Gas),intent(in) :: InputGas
                     Real(8) :: Vx,Vy,Vz,V2,V,Energy,Miu
                       OutputParticle%GasParticle=InputParticle
                       OutputParticle%PhaseSpace=InputParticle
                       OutputParticle%Mass=Mass 
                       Miu=Mass*InputGas%MGas/(InputGas%MGas+Mass)
                       Call Maxwellian(InputGas,OutputParticle%GasParticle)
                       Vx=InputParticle%Vx-OutputParticle%GasParticle%Vx
                       Vy=InputParticle%Vy-OutputParticle%GasParticle%Vy
                       Vz=InputParticle%Vz-OutputParticle%GasParticle%Vz
                       V2=Vx*Vx+Vy*Vy+Vz*Vz
                       V=Dsqrt(V2)
                       Energy=0.5d0*Miu*V2/JtoeV
                       OutputParticle%V2=V2
                       OutputParticle%V=V
                       OutputParticle%Energy=Energy
                    return  
                End subroutine UpdateParticleMCCIon 
end module MCCPublic