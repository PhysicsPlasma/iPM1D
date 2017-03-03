Module InitilalizationModule
Use FileIO
Use Constants
Implicit none
           contains
           Subroutine  ParticleInit(InputGas,NSpecy,InputSpecy,OutputParticleBundle,Period)  
            Implicit none
            Type(Gas), intent(in) :: InputGas 
            Integer(4), intent(in) :: NSpecy,Period
            Type(OneSpecy), intent(in) ::  InputSpecy(1:NSpecy)
            Type(ParticleBundle),intent(inout) :: OutputParticleBundle(1:NSpecy)
             Logical :: Status
             Integer(4) :: i
               do i=1,NSpecy
                    Call LoadParticle(InputSpecy(i),OutputParticleBundle(i),Status)
                     If (Status) then
                         OutputParticleBundle%Charge=ElectronCharge*InputSpecy%Charge
                         OutputParticleBundle%Mass=InputSpecy%Mass
                         Call PhaseSpaceInit1D(InputGas,InputSpecy(i),OutputParticleBundle(i))
                     else
                         Call UpdatePositionBundle(OutputParticleBundle(i),-1_4)
                      End If
                      Call UpdateVelocityBundle(OutputParticleBundle(i),-1_4) 
               End do
               OutputParticleBundle%Period=Period
                return
           End Subroutine  ParticleInit
 
             Subroutine PhaseSpaceInit1D(InputGas,InputSpecy,PB)
                      Use MCCModule
                      Implicit none
                      Type(Gas), intent(in) :: InputGas
                      Type(OneSpecy),intent(in) :: InputSpecy 
                      Type(ParticleBundle),intent(inout) :: PB
                      Type(ParticleOne) :: ParticleTemp
                      Integer(4) :: i,j,k
                      Real(8) :: XMin,XMax
                      PB%Name=InputSpecy%Name
                      XMin=PB%XMin
                      XMax=PB%XMax
                      do i=1,PB%NParNormal
                                 Call RandomPosition(ParticleTemp%X,XMin,XMax)
                                 Call Maxwellian(InputGas,ParticleTemp) 
                                 Call AddParticle(ParticleTemp,PB)
                          end do
 
                      return
                    end subroutine PhaseSpaceInit1D
                   
                 Subroutine RandomPosition(X,MinPosition,MaxPosition)
                              Implicit none
                              Real(8),intent(in) ::  MaxPosition,MinPosition 
                              Real(8),intent(inout) ::  X
                              Call DRandom(R)
                              X=MinPosition+R*(MaxPosition-MinPosition)
                              return 
                  end  Subroutine RandomPosition   
End Module InitilalizationModule