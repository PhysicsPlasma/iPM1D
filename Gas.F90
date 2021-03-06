Module GasModule
     Use TypeModule
     Use Constants
     Implicit none
     contains  
              subroutine  GasInit(GasPressure,GasTemperature,InputGasPhysical,OutputGas,Shift)
                  Implicit none
                  Integer(4), intent(inout) :: Shift
                  Real(8), intent(in) :: GasPressure,GasTemperature
                  Type(GasPhysical), intent(in) :: InputGasPhysical 
                  Type(Gas), intent(inout) ::  OutputGas
                  OutputGas%Shift=Shift
                  OutputGas%NSpecy=InputGasPhysical%NSpecy
                  OutputGas%MGas=InputGasPhysical%MGas
                  OutputGas%TGas=GasTemperature
                  OutputGas%PGas=GasPressure*mTorrtoPa
                  !OutputGas%BetaMax=InputGasPhysical%BetaMax
                  return
             End  subroutine  GasInit
 End Module  GasModule

Module SigmaModule
     Use TypeModule
     Use Constants
     Implicit none
     Real(8),parameter,private :: ElectronEnergyInterval=0.1d0,IonEnergyInterval=1.d0
     Integer(4),parameter,private :: NSigmaElectron=5000,NSigmaIon=1000
     !  The points of the Sigma data is simply devided into two kinds, 0-larger for Electron, 1-smaller for ion.
     contains
              subroutine  SigmaInit(Model,OutputSigma)
                  Implicit none
                  Integer(4), intent(in):: Model
                  Type(SigmaNormalized), intent(inout) :: OutputSigma
                  Real(8),external :: ExFunc
                  
                  Integer(4),save :: NSigma
                  Real(8),save ::  EnergyInterval
                   Integer(4) :: i,j,Index
                  Real(8) :: Energy,Alfa,LEnergy,UEnergy,S1,S2
                              Select case (Model)
                                  Case(0_4)
                                      EnergyInterval=ElectronEnergyInterval
                                      OutputSigma%EnergyInterval=EnergyInterval
                                      NSigma=NSigmaElectron 
                                      OutputSigma%NSigma=NSigma
                                  Case(1_4)
                                      EnergyInterval=IonEnergyInterval
                                      OutputSigma%EnergyInterval=EnergyInterval
                                      NSigma=NSigmaIon
                                      OutputSigma%NSigma=NSigma
                              End Select
                           return 
                   Entry  FuncSigmaNormalization(OutputSigma,ExFunc)
                            Index=OutputSigma%NReaction*NSigma
                             do i=0,NSigma-1
                                  Energy=dble(i)*EnergyInterval
                                  OutputSigma%Value(Index+i+1)=ExFunc(Energy)
                              end do
                              OutputSigma%NReaction=OutputSigma%NReaction+1
                             return
               End  subroutine  SigmaInit
             
             Subroutine  ReactiveIonSigmaNormalization(NComplex,InputReactionBundle,OutputSigma)
                      Implicit none
                      Integer(4),intent(in) ::  NComplex
                     Type(ReactionBundle),intent(in) :: InputReactionBundle
                     Type(SigmaNormalized), intent(inout) :: OutputSigma  
                     Real(8) :: EnergyInterval,Energy,S,PTemp(InputReactionBundle%NReaction),Threshold,PSum
                     Integer(4) :: i,j,Index,NSigma,NReaction
                     
                     EnergyInterval=OutputSigma%EnergyInterval
                     NSigma=OutputSigma%NSigma 
                     
                     NReaction=InputReactionBundle%NReaction
                     OutputSigma%NReaction=NReaction
                      
                     S=(3.d0*dble(NComplex)-6.d0)/2.d0-1.d0
                     
                     Index=0 
                     do i=0,NSigma-1
                             Energy=EnergyInterval*dble(i)
                             do j=1, NReaction
                                   Threshold=InputReactionBundle%Reaction(j)%Threshold 
                                   If(Energy>Threshold)  then
                                      PTemp(j)=(Energy-Threshold)**S
                                  else
                                      PTemp(j)=0.d0
                                  end if
                             end do
                                     
                            PSum=0.d0
                            do  j=1, NReaction
                                   Psum=Psum+PTemp(j)
                            end do
                            
                            do  j=1, NReaction
                                  OutputSigma%Value(Index+j)=PTemp(j)/PSum
                             end do
                             Index=Index+NReaction              
                    end do
        return
     end subroutine  ReactiveIonSigmaNormalization  
 End module SigmaModule

Module EnergyKai
     Use Constants
     Implicit none
     contains
     Function IsotropicCosKai(Energy) 
       Implicit none
       Real(8) :: IsotropicCosKai,Energy
       Call DRandom(R)
       IsotropicCosKai=1.d0-2.d0*R
       If(abs(IsotropicCosKai)<MinReal) then
                 If(IsotropicCosKai>0.d0)   then
                     IsotropicCosKai=MinReal
                else
                     IsotropicCosKai=-MinReal
                end if
        end  if     
       return
    end Function IsotropicCosKai
    
     Function EArCosKai(energy)
       !!!!!Warinig: ArcosTheta cant not be equal to +-1.0.
	    implicit none
    	real(8) :: energy
        real(8) :: EArCosKai
        if(energy < 1.d-30) then
		     EArCosKai = 1.d0
          else
		     Call DRandom(R)  
		     EArCosKai= (energy +2.0d0 -2.0d0*(energy+1.0d0)**R)/energy
		end if
          If(abs(EArCosKai)<MinReal) then
                 If(EArCosKai>0.d0)   then
                     EArCosKai=MinReal
                else
                     EArCosKai=-MinReal
                end if
         end  if   
 		return 
	  end function EArCosKai  
   
    Function IsotropicEnergy(Energy)
        implicit none
    	Real(8) :: Energy
        Real(8) :: IsotropicEnergy
        Call DRandom(R)
        IsotropicEnergy=R*Energy
        If(abs(IsotropicEnergy)<MinReal) then
                   IsotropicEnergy=MinReal
         end  if   
        return
    end Function  IsotropicEnergy
   
     Function EArEnergy(Energy)
        implicit none
    	Real(8) :: Energy
        Real(8) :: EArEnergy
        Call DRandom(R)
        EArEnergy=10.d0*Dtan(R*DAtan(Energy)/20.d0)
        return
    end Function  EArEnergy
End Module  EnergyKai    