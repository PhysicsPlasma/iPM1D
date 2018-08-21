! This version is finished by 2016/10/12, the code can run well and it can be a produciable code, but still improvements are possible, please refer next version.
    
    
    
Program PIC_MCC_for_CCP
   Use OneStepModule
   Use FileIO 
   USE Input 
   Use WeightingModule
   Use DiagOneStepModule
   Implicit none
   Integer(4) :: i,j,k
   real(8) Cpu1,Cpu2
   !Integer(4) :: NRun=0,NDiagShort=1,NDiagLong=0
   !Integer(4) :: NRun=20,NDiagShort=1,NDiagLong=0
    Integer(4) :: NRun=5000,NDiagShort=200,NDiagLong=200
   Call Initilalization()
   ParticleBDOneGlobal(0)%Gamma=0.2d0
   ParticleBDOneGlobal(1)%Gamma=0.2d0
   !Call InitilalizationDiag()
  Write(*,*) 'aa',NParticle,Period
  do j=1,NRun!00  !10  !200
           do i=1,Period
                 Call OneStep()
                 If (ParticleGlobal(0)%Npar>ParticleGlobal(0)%NParNormal) then
                     do k=0,1                
                             Write(*,*) ParticleGlobal(k)%Npar,k,"before"
                             Call PhaseSpaceNormalization2(ParticleGlobal(k)%Npar/2,ParticleGlobal(k))
                             Write(*,*) ParticleGlobal(k)%Npar,k,"after" 
                     end do
                 End If
            ENDDO
            !Call UpdateFieldBounday(FieldBoundaryGlobal,NParticle,ParticleBDOneGlobal)
            Call OneStepRestart()
            Write(*,*) i,Period,ParticleGlobal%NPar,j,"Period",FieldBoundaryGlobal%Vdc
    End do
  

     Call DiagInitilalization()
     do j=1,NDiagShort
         do i=1,Period
             Call OneStep()
             Call DiagOneStep()
         End do
     ENd Do
     Call DiagOneStepFinal()
    
    Do k=0,1
           Call DiagParticleElectrode(ParticleElectrodeGlobal(k),ParticleGlobal(k),FieldGlobal,-1)
     End do
     !Call DiagParticleTestParticle(ParticleGlobal(0),FieldGlobal,-1)
     
     do j=1,NDiagLong
         do i=1,Period
             Call OneStep()
              DO k=0,NParticle
                !Write(*,*) ParticleBDOneGlobal(k)%CountMinOne,ParticleBDOneGlobal(k)%CountMaxOne,ParticleBDOneGlobal(k)%CountMin,ParticleBDOneGlobal(k)%CountMax,k,i
               Call DiagParticleElectrode(ParticleElectrodeGlobal(k),ParticleGlobal(k),FieldGlobal,0)
              ENd do
              !Call DiagParticleTestParticle(ParticleGlobal(0),FieldGlobal,0)

         End do
         If (Mod(j,50)==0) then
         Do k=0,1
           Call DiagParticleElectrode(ParticleElectrodeGlobal(k),ParticleGlobal(k),FieldGlobal,1)
         End do
         ENd if
         Write(*,*) ParticleElectrodeGlobal(0)%PBLower%Npar,ParticleElectrodeGlobal(0)%PBUpper%Npar,ParticleBDOneGlobal(0)%CountMin,ParticleBDOneGlobal(0)%CountMax
     ENd do
     Do k=0,1
           Call DiagParticleElectrode(ParticleElectrodeGlobal(k),ParticleGlobal(k),FieldGlobal,1)
     End do
     !Call DiagParticleTestParticle(ParticleGlobal(0),FieldGlobal,1)
stop
end  Program  

