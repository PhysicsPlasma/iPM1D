Module Constants
    Implicit none
    Real(8),parameter :: PI=3.141592653589793238D0
    Real(8),parameter :: a0=5.2918d-11
    Real(8),parameter:: lightspeed=3.0d8 
    Real(8),parameter :: kB=1.3807d-23
    Real(8),parameter :: exp=2.71828d0
    Real(8),parameter :: ElectronCharge=1.6022d-19
    Real(8),parameter :: ElectronMass=9.1095d-31
    Real(8),Parameter :: Epsilon=8.8542d-12
    Real(8),Parameter :: JtoeV=1.6022d-19
    Real(8),Parameter :: eVtoK=11605.d0
    Real(8),Parameter :: mTorrtoPa=0.13332d0
    Real(8),Parameter :: MinReal=1.d-15
    Real(8)  :: R
    contains
           Subroutine DRandomInt(RR)
                USE IFPORT
                Implicit none
                 Real(8),intent(inout) :: RR
                 Integer :: AA
                 AA=Int(RR*100.d0)
                 RR=DRAND(AA)
                 Return 
           Entry  DRandom(RR)
                   RR=DRAND(0)
                   !CALL RANDOM_NUMBER(RR)
                   !Write(*,*) RR
                Return
          End  subroutine DRandomInt 
End  Module Constants

Module TypeModule
    Implicit none
    !  This Section is the definitions of the physicla parameters.

              !
              Integer(4),Parameter :: NxMax=113
              Real(8),parameter,private :: ZLength=0.035d0
              Real(8),parameter :: Inputdx=ZLength/dble(NxMax-1)
              Real(8),parameter:: Inputdt=1.d-10  !4.d-12
              
              ! Type :: Test2
              !          Real(8), Allocatable :: Test2(:)
              !End Type

              Type :: Grid2D(Nx,Ny,Ns)
                        Integer(4),Len :: Nx=NxMax,Ny=1000,Ns=1_4
                        Real(8) ::  Value(1:Nx,0:Ny,1:Ns)
                         Real(8) :: Dx=0.d0,Dy=0.d0
                         Integer(4) :: Timer=0,Period=1
                        !Integer(4) :: NxAva=1,NyAva=1
                End Type Grid2D
     
              Type :: Grid1D(Nx,Ns)
                        Integer(4),Len :: Nx=NxMax,Ns=1_4
                        Real(8) ::  Value(1:Nx,1:Ns)
                        Integer(4) :: Timer=0,Period=1
                        !Integer(4) :: NxAva=1
                        Real(8) :: Dx=0.d0
              End Type Grid1D

              Type  OneSpecy      ! Define a physical particle.
                       Character(len=16) Name
                       Integer(4) :: NAtom,Charge
                       Real(8) :: Mass
                End  Type  OneSpecy
                

                
                Type GasPhysical    ! Define a physical gas.
                     Integer(4) ::  Model,NAtom,NSpecy
                     Real(8) ::  MGas,Radius,BetaMax 
                End Type  GasPhysical
                
                   !  ReactionType is definition of the reactions:
                   !   -1-The injected particle are removed because the generated particle is no longer included in simulations . 
                   !   0-no collision
                   !   1~99 for electrons  (1~10 Isotropic : 1-Elastic; 2-Exitation; 3- Ionization; 4-Attachment; 5-Dissociation)
                   !                                    (11-20 Anisotropic Ar: 11-Elastic; 12-Exitation; 13- Ionization;)
                   !   100~199 for  ions  (101~110 non-reactive : 101-isotropic; 102-Anisotropic; 103-ChargeExchange)
                   !                                    (111~120 reactive : 111-original particle; 112-new particle)    
                Type  ReactionOne    ! Define a collosion process.
                     Integer(4) :: Reactant,ReactionType,Resultant
                     Real(8) ::  Threshold
                End Type  ReactionOne
    
                
                Type ParticleOne
                    Real(8) ::  X,Vx,Vy,Vz,Ax,Ay,Az
                EndType ParticleOne
                
                !Type,extends(ParticleOne) :: ParticleOneTimer
                !    Real(8) ::  Timer
                !End Type ParticleOneTimer

                Integer(4),parameter :: ParticlePerGrid=400  !100 
                Real(8),parameter :: InitDensity=1.d15  !1.d18 !5.d17  !1.d16
                Real(8),parameter :: Weighting=InitDensity/ParticlePerGrid
                Integer(4),Parameter :: NParMax=ParticlePerGrid*NxMax*5.0d0
                !Integer(4),Parameter :: NSpecy=2
    
                Type  :: ParticleBundle !(NParMax)
                    !Integer(4),Len :: NParMax=0
                    Integer(4) :: SpecyIndex
                    Character(len=99) Name
                    Integer(4) :: Timer=0,Period=1
                    Integer(4) :: NPar=0,NParNormal=(NxMax-1)*ParticlePerGrid
                    Real(8) :: Charge,Mass,Weight=Weighting
                    Real(8) :: XFactor=Inputdx,VFactor=Inputdx/Inputdt
                    Real(8) :: dx=Inputdx,dt=Inputdt
                    Real(8) :: XMin=0.d0,XMax=dble(NxMax-1)!,XLength=99.d0
                    Type(ParticleOne) :: PhaseSpace(NParMax)
                EndType ParticleBundle
                
                !Type  :: ParticleBundleTimer !(NParMax)
                !    !Integer(4),Len :: NParMax=0
                !    Integer(4) :: SpecyIndex
                !    Character(len=99) Name
                !    Integer(4) :: Timer=0,Period=1
                !    Integer(4) :: NPar=0,NParNormal=(NxMax-1)*ParticlePerGrid
                !    Real(8) :: Charge,Mass,Weight=Weighting
                !    Real(8) :: XFactor=Inputdx,VFactor=Inputdx/Inputdt
                !    Real(8) :: dx=Inputdx,dt=Inputdt
                !    Real(8) :: XMin=0.d0,XMax=dble(NxMax-1)!,XLength=99.d0
                !    Type(ParticleOneTimer) :: PhaseSpace(NParMax)
                !EndType ParticleBundleTimer

                Type Field !(Nx)
                        !Integer(4),Len :: Nx=100
                        !Integer(4) :: FieldModel=101
                        Integer(4) :: Nx=NxMax
                        Integer(4) :: Timer=0,Period=1
                        !Integer(4) :: XStart=0,XEnd=NxMax-1
                        Real(8) :: Dx=Inputdx,Dt=Inputdt
                        Real(8) ::  Ex(1:NxMax)=0.d0,Ey(1:NxMax)=0.d0!,Ey(0:Nx+1),Ez(0:Nx+1)
                        Real(8) ::  Bx(1:NxMax)=0.d0,By(1:NxMax)=0.d0!,Bz(0:Nx+1)
                        Real(8) ::  Rho(1:NxMax),Phi(1:NxMax)
                        Real(8) ::  Chi(1:NxMax)
                   EndType Field
                   
                    Type FieldSolver
                        !Integer(4) :: XStart=0,XEnd=1
                        Integer(4) :: Ns=NxMax-2
                        Real(8) :: Dx=Inputdx,Dt=Inputdt
                        Real(8) :: Source(1:NxMax-2)
                        Real(8) :: Solve(1:NxMax-2)
                        Real(8) :: CoeA(1:NxMax-2),CoeB(1:NxMax-2),CoeC(1:NxMax-2)
                    EndType FieldSolver

                    Type FieldOne !(Nx)
                        !Integer(4),Len :: Nx=100
                        Integer(4) :: Nx=NxMax
                        !Integer(4) :: XStart=0,XEnd=NxMax-1
                        Real(8) :: Dx=Inputdx,Dt=Inputdt
                        !Real(8) :: QdM
                        Real(8) ::  RhoOne(1:NxMax),ChiOne(1:NxMax)
                        !Real(8) ::  JxOne(1:NxMax),JyOne(1:NxMax),JzOne(1:NxMax)
                        !Real(8) ::  TOne(1:NxMax)
                        !Real(8) ::  Ex(1:Nx)!,Ey(0:Nx+1),Ez(0:Nx+1)
                        !Real(8) ::  Bx(1:Nx),By(1:Nx)
                    EndType FieldOne

                        !  This tpye defines the values for the collision praticle.
                         !  Index- ReactionType, same as the ReactionType defined below.
                         !  PhaseSpace- the position and velocity of the indection particle.
                         !  GasParticle- the position and velocity of the Gas particle .
                   Type MCCParticleOne
                          Integer(4) :: Index                            
                          Type(ParticleOne) :: PhaseSpace
                          Real(8) :: Mass,V2,V,Energy,Beta
                          Type(ParticleOne) ::  GasParticle
                   End Type  MCCParticleOne
      
   !!  This Sections Defines the grids  used in the simulations.
   !                   Integer(4),parameter,private :: GridMax=5*1024
   !                   Integer(4),parameter,private :: GridMaxSmall=100000_4 

    
     ! This  Sections Defines the gas and collision properties used For PIC simulations. 
                       Type Gas
                             Integer(4) ::  NSpecy,Shift
                             Real(8) :: MGas,TGas=300.d0,PGas=50.d0,BetaMax 
                       End Type  Gas
                           
                       Integer(4),parameter ::  NReactionMax=3_4
                       Integer(4),parameter,private :: NSigmaMax=150000_4 
                         
                       Type  ReactionBundle        !  Reactive bundle for all collision process of one particle withe one gas.
                             Integer(4) ::  NReaction 
                             Type(ReactionOne) ::  Reaction(1:NReactionMax)
                        End Type ReactionBundle
                       
                       !  Mcc bundle for all collision process of one particle withe one gas. 
                       ! Model is the index for the different MCC models:
                       ! 0-no collision; 1-elelctrons,the gas moculars are stationary; 2-nonreactive ions; 3-reactive ions;
                         
                        Type MCCBundle
                            Integer(4) ::  Model,NReaction=0,NSigma=0
                            Real(8) :: Mass,CollisionRatio,EnergyInterval
                            Real(8) :: Ng,SigmaMax
                            Real(8) :: Probility(1:NSigmaMax)
                            Type(ReactionOne) :: Reaction(0:NReactionMax)
                        End Type MCCBundle
                        
                        Integer(4),parameter,private :: NSigmaRaw=100_4
                              Type  SigmaRaw
                                   Real(8)  :: Threshold,MaxEnergy
                                   Integer(4) :: NSigma
                                   Real(8) :: Energy(NSigmaRaw),Sigma(NSigmaRaw)
                              end Type  SigmaRaw
                              
                              Type SigmaNormalized
                                  Real(8) :: EnergyInterval=0.d0
                                  Integer(4) :: NReaction=0,NSigma=0
                                  Real(8) ::  Value(1:NSigmaMax)
                              End Type SigmaNormalized
End Module TypeModule