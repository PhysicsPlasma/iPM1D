Module Input
   implicit none
    !  All in SI unit except gas pressure in mTorr.
    !  This section defines the physical dimensions
           !Integer(4),parameter :: Dims=2
           !Integer(4),parameter :: Nz=65,Nr=1
           !Real(8),parameter,private :: ZLength=0.02d0,RLength=0.08d0
           !Real(8),parameter :: Inputdx=ZLength/dble(Nz-1)
           !Real(8),parameter :: Inputdt=1.d-10  !4.d-12
           !Integer(4),parameter :: InputGap=Int(0.02d0/Inputdx)
           
    !  This  section defines the numrical setting.
           !Integer(4),parameter :: ParticlePerGrid=100  !100 
           !Real(8),parameter :: InitDensity=1.d16  !1.d18 !5.d17  !1.d16
           !Real(8),parameter :: Weighting=InitDensity/ParticlePerGrid
       
    !   This section defines the external rf sources.
           !Integer(4) :: RFKind=11
           !Real(8) :: RFFrequency(2)=(/13.56d6,60.d6/),RfVol(2)=(/-200.d0,-100.d0/)
  ! Real(8) :: RFFrequency(2)=(/1.1d5,60.d6/),RfVol(2)=(/600.d0,-100.d0/)
   
    !   This section defines the controlable properties of backgroud gas.
            Integer(4) :: InputNGas=1
            Integer(4) :: GasDefinition(1:3)=(/1,0,0/)
            Real(8) :: GasPressure(1:3)=(/50.d0,0.d0,0.d0/),GasTemperature(1:3)=(/300.d0,300.d0,300.d0/)
end Module Input




 
