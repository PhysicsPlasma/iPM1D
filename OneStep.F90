Module OneStepModule
      Use FieldModule
      Use MCCPublic
      Use ParticleBoundaryModule
      Use Diagnostics
!      Use Diagnostics
      Implicit none
      ! This section defines the particles.
                  Integer(4),parameter :: NSpecyMax=1_4
                  Integer(4),save :: NParticle=0  
                  Type(ParticleBundle),save ::  ParticleGlobal(0:NSpecyMax)
                  Type(ParticleBoundaryOne),save  :: ParticleBDOneGlobal(0:NSpecyMax)
                  Type(ParticleBoundary),save  :: ParticleBDGlobal
      
      !  This section defines the MCC and gas properties.
                 Integer(4),parameter :: NSpecyGasMax=1_4
                 Integer(4),parameter :: NGasMax=1_4
                 Integer(4),save ::  Ngas=1
                 Type(Gas),save  ::  GasGlobal(1:NGasMax)
                 Type(MCCBundle),save ::  MCCBundleGlobal(0:NSpecyGasMax,1:NGasMax)
       
       !  This section defines the PIC properties.
                 Integer(4),save :: Period
                 Real(8),save :: dx,dt
                 Type(FieldOne),save :: FieldOneGlobal(0:NSPecyMax)
                 Type(Field),save :: FieldGlobal
                 Type(FieldSolver),save  :: FieldSolverGlobal
                 Type(FieldBoundary),save  :: FieldBoundaryGlobal
                 
                 
                  
    contains
    Subroutine Initilalization()
             Use GasModule
             Use InitilalizationModule
             Use MCCPublic
             Use ElectronModule
             Use ArGasModule
             Use Input
             Implicit none
             Integer(4) ::  i,NSpecy,Shift=0
             Ngas=InputNGas

             !Call BoundaryLoadInit(RFFrequency(1),dt,Period)
             Call FieldBoundayInitilalization(FieldBoundaryGlobal,Period,Inputdt)
             FieldGlobal%Period=Period
            ! ParticleInit(InputGas,NSpecy,InputSpecy,OutputParticleBundle)  
             Call ParticleInit(ElectronGas,1_4,Electron,ParticleGlobal(0),Period)
             ParticleGlobal(0)%SpecyIndex=0
             do i=1,Ngas
                       select case (GasDefinition(i))
                          Case(1_4) ! Ar
                                Call ArSigmaInitialization()
                                Call GasInit(GasPressure(i),GasTemperature(i),ArGas,GasGlobal(i),Shift)
                                Call MCCBundleInit(ArGas,GasGlobal(i),ArSpecy,ArReaction,ArSigma,MCCBundleGlobal(0:,i),Inputdt)  
                                NSpecy=ArGas%NSpecy
                               Call ParticleInit(GasGlobal(i),NSpecy,ArSpecy,ParticleGlobal(Shift+1:Shift+NSpecy),Period)
                                Shift=Shift+NSpecy
                          Case(2_4) !CF4
                        End Select
             End do

             ParticleGlobal(1:NSpecyMax)%SpecyIndex=1
             ParticleBDOneGlobal(0:NSpecyMax)%PBLower=ParticleGlobal(0:NSpecyMax)
             ParticleBDOneGlobal(0:NSpecyMax)%PBUpper=ParticleGlobal(0:NSpecyMax)
             ParticleBDOneGlobal(0:NSpecyMax)%PBLower%NPar=0
             ParticleBDOneGlobal(0:NSpecyMax)%PBUpper%NPar=0
                  NParticle=Shift
                  Write(*,*) NParticle
                  return  
    End Subroutine Initilalization

    Subroutine OneStep()
               Use MoveModule 
               Use WeightingModule
               Use FieldModule 
               Use SEEModule
             !  Use ,only : R
               Implicit none
            !   Real(8)::GamaE=0.2d0,GamaI=0.2d0
               Integer(4) :: i,j
               do i=0,NParticle
                     Call ParticleMove(ParticleGlobal(i),FieldGlobal)
                     Call ParticleAborption(ParticleGlobal(i),ParticleBDOneGlobal(i))
                   !  Call Selectron(ParticleGlobal(0),ParticleBDOneGlobal(i))
                     Call WeightingOne(ParticleGlobal(i), FieldOneGlobal(i))
               end do

                Call FieldOneStep(FieldGlobal,NParticle,FieldOneGlobal,FieldSolverGlobal,FieldBoundaryGlobal)
                               !Write(*,*) FieldBoundaryGlobal%V1,FieldBoundaryGlobal%V2
               !Write(*,*)  ParticleGlobal%NPar,'Before'  
                do i=1, Ngas
                   do j=0,NParticle
                       Call MCC(j,NParticle,ParticleGlobal,GasGlobal(i),MCCBundleGlobal(j,i)) 
                   end do
               end do
              return
    End  subroutine OneStep
    
        Subroutine OneStepRestart()
               Use FileIO 
               Implicit none
               Integer(4) :: i
               Write(*,*) Period,ParticleGlobal%NPar,"Period"!,FieldBoundaryGlobal%Qmin,FieldBoundaryGlobal%Qmax
               Call DumpFieldSolver(FieldSolverGlobal,0)
               Call DumpField(FieldGlobal,0)
               do i=0,1
                     Call DumpParticle(ParticleGlobal(i),0)
               End do
               Call FieldBoundayFinalization(FieldBoundaryGlobal)
              return
        End  subroutine OneStepRestart
 
    End Module OneStepModule
    
    
  Module DiagOneStepModule
      Use OneStepModule
      Use Diagnostics
      Implicit none
      Type(Grid1D(Nx=NxMax,Ns=11)),save :: G1DDiagParticleField
      Type(Grid2D(Nx=NxMax,Ny=100,Ns=11)),save :: G2DDiagParticleField
      
      Type(Grid1D(Nx=NxMax,Ns=3)),save :: G1DDiagParticleCR
      Type(Grid2D(Nx=NxMax,Ny=100,Ns=3)),save :: G2DDiagParticleCR

      Type(Grid1D(Nx=500,Ns=2)),save :: G1DDiagParticleEDF
      Type(Grid2D(Nx=500,Ny=100,Ns=2)),save :: G2DDiagParticleEDF
      
      Type(ParticleElectrode),save :: ParticleElectrodeGlobal(0:NSpecyMax)
      
      contains
      Subroutine DiagInitilalization()
            Implicit none
                 Call DiagParticleFieldPeriod(G1DDiagParticleField,1,ParticleGlobal,FieldGlobal,-1)
                 Call DiagParticleFieldPeriod(G2DDiagParticleField,1,ParticleGlobal,FieldGlobal,-1)
                 Call DiagParticleEDFOne(G1DDiagParticleEDF,ParticleGlobal(0),-1)
                 Call DiagParticleEDFOne(G2DDiagParticleEDF,ParticleGlobal(0),-1)
                 Call DiagParticleCollisionRateOne(G1DDiagParticleCR,ParticleGlobal(0),MCCBundleGlobal(0,1),-1)
                 Call DiagParticleCollisionRateOne(G2DDiagParticleCR,ParticleGlobal(0),MCCBundleGlobal(0,1),-1)
            return  
      End Subroutine DiagInitilalization
      
      Subroutine  DiagOneStep()
          Implicit none   
             Call DiagParticleFieldPeriod(G1DDiagParticleField,1,ParticleGlobal,FieldGlobal,0)
             Call DiagParticleFieldPeriod(G2DDiagParticleField,1,ParticleGlobal,FieldGlobal,0)
             Call DiagParticleEDFOne(G1DDiagParticleEDF,ParticleGlobal(0),0)
             Call DiagParticleEDFOne(G2DDiagParticleEDF,ParticleGlobal(0),0)
             Call DiagParticleCollisionRateOne(G1DDiagParticleCR,ParticleGlobal(0),MCCBundleGlobal(0,1),0)
             Call DiagParticleCollisionRateOne(G2DDiagParticleCR,ParticleGlobal(0),MCCBundleGlobal(0,1),0)
        return  
      End Subroutine DiagOneStep
      
       Subroutine DiagOneStepFinal()
       Implicit none  
             Call DiagParticleFieldPeriod(G1DDiagParticleField,1,ParticleGlobal,FieldGlobal,1)
             Call DiagParticleFieldPeriod(G2DDiagParticleField,1,ParticleGlobal,FieldGlobal,1)
             Call DiagParticleEDFOne(G1DDiagParticleEDF,ParticleGlobal(0),1)
             Call DiagParticleEDFOne(G2DDiagParticleEDF,ParticleGlobal(0),1)
             Call DiagParticleCollisionRateOne(G1DDiagParticleCR,ParticleGlobal(0),MCCBundleGlobal(0,1),1)
             Call DiagParticleCollisionRateOne(G2DDiagParticleCR,ParticleGlobal(0),MCCBundleGlobal(0,1),1)
          return  
        End Subroutine DiagOneStepFinal

  End Module DiagOneStepModule      

    !

    !
    !

        
        
        
        

