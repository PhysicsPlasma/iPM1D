 Module WeightingModule
         Use TypeModule 
         Use Constants
         implicit none
         contains 
         subroutine WeightingOne(PB,FO)
                implicit none
                Type(ParticleBundle),intent(in) :: PB
                Type(FieldOne),intent(inout) :: FO
                Real(8) :: S1,S2
                Integer(4) :: i,N
                Real(8) :: RhoFactor,ChiFactor!,JFactor,TFactor
                FO%RhoOne=0.d0
                do i=1,PB%Npar
                   N=Ceiling(PB%PhaseSpace(i)%X)
                   S1=Dble(N)-PB%PhaseSpace(i)%X
                   S2=1.d0-S1
                   FO%RhoOne(N)=FO%RhoOne(N)+S1
                   FO%RhoOne(N+1)=FO%RhoOne(N+1)+S2
                end do
                RhoFactor=PB%Charge*PB%Weight
                FO%RhoOne=FO%RhoOne*RhoFactor
                ChiFactor=0.5d0*PB%Charge/PB%Mass*PB%dt*PB%dt/Epsilon
                FO%ChiOne=FO%RhoOne*ChiFactor
                
             return
         end subroutine WeightingOne
         
         
End Module WeightingModule