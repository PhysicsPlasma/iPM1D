    Module FieldBoundaryModule
  Use TypeModule
  Use Constants
  Use ParticleBoundaryModule
  Implicit none
     Type FieldBoundary
                        !Integer(4) :: XStart=0,XEnd=1
                         Integer(4) :: FieldBoundaryModel=11
                         Integer(4) :: Timer=0,Period
                         Real(8) :: Dt=Inputdt
                         Real(8) :: Frequency(2)=(/13.56d6,2.d6/),Voltage(2)=(/300.d0,100.d0/)
                         Real(8) :: V1=0.d0,V2=0.d0
                         Real(8) :: Vdc=0.d0!-20.d0
     EndType FieldBoundary
     
     
     
     !Type FieldBoundary
     !                   !Integer(4) :: XStart=0,XEnd=1
     !                    Integer(4) :: FieldBoundaryModel=11
     !                    Integer(4) :: Timer=0,Period
     !                    Real(8) :: Dt
     !                    Real(8) :: Frequency(2)=(/13.56d6,2.d6/),Voltage(2)=(/200.d0,100.d0/)
     !                    Real(8) :: Capacitor=1.d-9,Area=PI*0.1d0*0.1d0
     !                    Real(8) :: V1,V2
     !                    Real(8) :: PhiNMin=0.d0,PhiNmax=0.d0
     !                    Real(8) :: Qcap=0.d0,Qconv=0.d0
     !                    Real(8) :: Sigma=0.d0
     !EndType FieldBoundary
  !Real(8),private :: Bottom(100000),Top(100000)
       contains
     Subroutine FieldBoundayInitilalization(FB,Period,dt)
               Implicit none
               Type(FieldBoundary),intent(inout) :: FB
               Integer(4), intent(inout) ::  Period
               Real(8) :: Dt
               Character(len=99) :: Filename
               Logical :: alive
               Write(filename,*) "10000","FieldBoundary",".dat"
               Inquire(file=filename,exist=alive)
               If(alive) then
                   Open (10,file=filename)
                   Read(10,*) FB%Timer,FB%V1,FB%V2,FB%Vdc
                   Close(10)
               Else
                   Write(*,*) 'Can not find the file for ', filename,' the FieldBoundary will be set to Zero.' 
                   FB%Timer=0
                   FB%Vdc=0.d0
                   FB%V1=0.d0
                   FB%V2=0.d0
               End If
               Period=Int(1.d0/FB%Frequency(1)/dt)
               FB%Dt=dt
               FB%Period=Period
               return
     End Subroutine FieldBoundayInitilalization
     
     Subroutine FieldBoundayOneStep(FS,FG,FB)
            implicit none
            Type(FieldSolver),intent(inout) :: FS
            Type(Field),intent(inout) :: FG
            Type(FieldBoundary),intent(inout) :: FB
            Integer i,Nx
            Real(8),save ::  GeoFactor
            
            Select case (FB%FieldBoundaryModel)
                      case (11)
                               FB%V1=FB%Vdc+FB%Voltage(1)*DSin(2*PI*FB%Frequency(1)*FB%dt*Dble(FB%Timer))
                               FB%V2=0.d0
                      case (21)
                               FB%V1=FB%Vdc+FB%Voltage(1)*DCos(2*PI*FB%Frequency(1)*FB%dt*Dble(FB%Timer))+FB%Voltage(2)*DCos(2*PI*FB%Frequency(2)*FB%dt*Dble(FB%Timer))
                               FB%V2=0.d0
                      case(22)
                               FB%V1=FB%Vdc+FB%Voltage(2)*DCos(2*PI*FB%Frequency(2)*FB%dt*Dble(FB%Timer))
                               FB%V2=FB%Vdc+FB%Voltage(1)*DCos(2*PI*FB%Frequency(1)*FB%dt*Dble(FB%Timer))  
              end  Select
            
            GeoFactor=-FG%dx*FG%dx/Epsilon
            Nx=FG%Nx
            FS%Source(1:Nx-2)=FG%Rho(2:Nx)*GeoFactor
            FS%Source(1)=FS%Source(1)-FB%V1*FS%CoeA(1)
            FS%Source(Nx-2)=FS%Source(Nx-2)-FB%V2*FS%CoeC(Nx-2)

            FB%Timer=FB%Timer+1
            return
     end subroutine  FieldBoundayOneStep
     
     Subroutine FieldBoundayFinalization(FB)
               Implicit none
               Type(FieldBoundary),intent(inout) :: FB
               Character(len=99) :: Filename
               Write(filename,*) "10000","FieldBoundary",".dat"
               Open (10,file=filename)
               Write(10,*) FB%Timer,FB%V1,FB%V2,FB%Vdc
               Close(10)
               Write(*,*) "Saving ",filename," Please wait..."
               return
     End Subroutine FieldBoundayFinalization
     
     Subroutine UpdateFieldBounday(FB,NSpecy,PBDO)
               Implicit none
               Type(FieldBoundary),intent(inout) :: FB
               Integer(4), intent(in) :: NSpecy
               Type(ParticleBoundaryOne),intent(inout) :: PBDO(0:NSpecy)
               Integer(4) :: i
               Real(8) :: NetCharge
               NetCharge=0.d0
               do i=0,NSpecy
                    NetCharge=NetCharge+PBDO(i)%CountMin*PBDO(i)%PBLower%Charge*PBDO(i)%PBLower%Weight
                    PBDO(i)%CountMin=0
               End Do
               If (NetCharge<0.d0) then
                   FB%Vdc=FB%Vdc-1.d0
               Else if (NetCharge>0.d0) Then
                   FB%Vdc=FB%Vdc+1.d0
                ENd IF
                return
     End Subroutine UpdateFieldBounday
                  
    !   Subroutine FieldBoundayOneStep(FS,FG,FB)
    !        implicit none
    !        Type(FieldSolver),intent(inout) :: FS
    !        Type(Field),intent(inout) :: FG
    !        Type(FieldBoundary),intent(inout) :: FB
    !        Integer i,Nx
    !        Real(8),save ::  GeoFactor
    !        
    !        Select case (FB%FieldBoundaryModel)
    !                  case (11)
    !                           FB%V1=FB%Voltage(1)*DCOs(2*PI*FB%Frequency(1)*FB%dt*Dble(FB%Timer))
    !                           FB%V2=0.d0
    !                  case (21)
    !                           FB%V1=FB%Voltage(1)*DCos(2*PI*FB%Frequency(1)*FB%dt*Dble(FB%Timer))+FB%Voltage(2)*DCos(2*PI*FB%Frequency(2)*FB%dt*Dble(FB%Timer))
    !                           FB%V2=0.d0
    !                  case(22)
    !                           FB%V1=FB%Voltage(2)*DCos(2*PI*FB%Frequency(2)*FB%dt*Dble(FB%Timer))
    !                           FB%V2=FB%Voltage(1)*DCos(2*PI*FB%Frequency(1)*FB%dt*Dble(FB%Timer))  
    !          end  Select
    !        
    !        GeoFactor=-FG%dx*FG%dx/Epsilon
    !        Nx=FG%Nx
    !        FS%Source(1:Nx-2)=FG%Rho(2:Nx)*GeoFactor
    !        FS%Source(1)=FS%Source(1)-FB%V1*FS%CoeA(1)
    !        FS%Source(Nx-2)=FS%Source(Nx-2)-FB%V2*FS%CoeC(Nx-2)
    !
    !        FB%Timer=FB%Timer+1
    !        return
    !end subroutine  FieldBoundayOneStep
end  Module FieldBoundaryModule

    
    
    Module FieldModule
  Use TypeModule
  Use FileIO
  Use Constants
  Use FieldBoundaryModule
  Implicit none
  !Real(8),private :: Bottom(100000),Top(100000)
  contains

                !GeoFactor=-dx*dx/Epsilon

            subroutine  FieldOneStep(FG,NSpecy,FO,FS,FB)
                Implicit none
                Type(Field),intent(inout) :: FG
                Integer(4),intent(in) :: NSpecy
                Type(FieldOne),intent(in) :: FO(0:NSPecy)
                Type(FieldSolver),intent(inout) :: FS
                Type(FieldBoundary),intent(inout) :: FB
                
                Call AddRhoChi(FG,NSpecy,FO)
                !FG%Rho=0.d0
                !FG%Chi=0.d0
                Call CalCoeff(FS,FG)
                Call FieldBoundayOneStep(FS,FG,FB)
                Call Tridag(FS%CoeA,FS%CoeB,FS%CoeC,FS%Source,FS%Solve,Fs%Ns)
                Call FaiToE(FS,FG,FB)
                FG%Timer=FG%Timer+1
                return
        End subroutine FieldOneStep
         
         subroutine AddRhoChi(FG,NSpecy,FO)
            Implicit none
            Type(FieldOne),intent(in) :: FO(0:NSPecy)
            Type(Field),intent(inout) :: FG
            Integer(4),intent(in) :: NSpecy
            Integer(4) :: i,Nx
            FG%Rho=0.d0
            FG%Chi=0.d0
            Nx=FG%Nx
            do i=0,NSpecy
                    FG%Rho=FG%Rho+FO(i)%RhoOne
                    FG%Chi=FG%Chi+FO(i)%ChiOne
            End do
            return
        End subroutine AddRhoChi
        
       subroutine  CalCoeff(FS,FG)
            implicit none 
            Type(FieldSolver),intent(inout) :: FS
            Type(Field),intent(in) :: FG
            Integer(4) :: i
            Real(8) :: ATemp, BTemp,CTemp
            do i=1,FS%Ns
                  ATemp=1.d0+0.5d0*(FG%Chi(i)+FG%Chi(i+1))
                  !ATemp=1.d0+Max(Chi(i),Chi(i+1))
                 FS% CoeA(i)=ATemp
                  !CTemp=1.d0+Max(Chi(i),Chi(i+1))
                  CTemp=1.d0+0.5d0*(FG%Chi(i+1)+FG%Chi(i+2))
                  FS%CoeC(i)=CTemp
                  
                  BTemp=-(ATemp+CTemp)
                  FS%CoeB(i)=BTemp
             End do
           return
       end subroutine  CalCoeff
     
       subroutine FaiToE(FS,FG,FB)
            implicit none
            Type(FieldSolver),intent(inout) :: FS
            Type(Field),intent(inout) :: FG
            Type(FieldBoundary),intent(in) :: FB
            Integer i,Nx
            !Real(8) :: PhiMin,PhiMax
            Nx=FG%Nx
            FG%Phi(2:Nx-1)=FS%Solve(1:Nx-2)
            FG%Phi(1)=FB%V1
            FG%Phi(Nx)=FB%V2
            
            !PhiMin=2.d0*FG%Phi(2)-FG%Phi(1)
            FG%Ex(1)=(3.d0*FG%Phi(1)-4.d0*FG%Phi(2)+FG%Phi(3))/(2.d0*FG%dx)
            do i=2,Nx-1
                FG%Ex(i)=(FG%Phi(i-1)-FG%Phi(i+1))/(2.d0*FG%dx)
            end do
            !PhiMax=2.d0*FG%Phi(Nx-1)-FG%Phi(Nx)
            FG%Ex(Nx)=-1.d0*(FG%Phi(Nx-2)-4.d0*FG%Phi(Nx-1)+3.d0*FG%Phi(Nx))/(2.d0*FG%dx)
          return
       end subroutine FaiToE
       
        !subroutine FaiToE(FS,FG,FB)
        !    implicit none
        !    Type(FieldSolver),intent(inout) :: FS
        !    Type(Field),intent(inout) :: FG
        !    Type(FieldBoundary),intent(in) :: FB
        !    Integer i,Nx
        !    !Real(8) :: PhiMin,PhiMax
        !    Nx=FG%Nx
        !    FG%Phi(2:Nx-1)=FS%Solve(1:Nx-2)
        !    FG%Phi(1)=FB%V1
        !    FG%Phi(Nx)=FB%V2
        !    !PhiMin=2.d0*FG%Phi(2)-FG%Phi(1)
        !    FG%Ex(1)=(3.d0*FG%Phi(1)-4.d0*FG%Phi(2)+FG%Phi(3))/(2.d0*FG%dx)
        !    do i=2,Nx-1
        !        FG%Ex(i)=(FG%Phi(i-1)-FG%Phi(i+1))/(2.d0*FG%dx)
        !    end do
        !    !PhiMax=2.d0*FG%Phi(Nx-1)-FG%Phi(Nx)
        !    FG%Ex(Nx)=-1.d0*(FG%Phi(Nx-2)-4.d0*FG%Phi(Nx-1)+3.d0*FG%Phi(Nx))/(2.d0*FG%dx)
        !  return
        !  end subroutine FaiToE

SUBROUTINE tridag(a,b,c,r,u,n)
    integer,PARAMETER ::nmax=10000
    Integer,intent(in) :: n
    real(8), intent(in) :: a(n),b(n),c(n),r(n)
    real(8), intent(inout) :: u(n)
    INTEGER j
    REAL(8) gam(nmax)
    real(8) bet
        !A=1.d0
        !B=-2.d0
        !C=1.d0
        bet=b(1)
        u(1)=r(1)/bet
    do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        u(j)=(r(j)-a(j)*u(j-1))/bet
    end do
    do j=n-1,1,-1
       u(j)=u(j)-gam(j+1)*u(j+1)
    end do
    return
   END SUBROUTINE tridag
    end  Module FieldModule
    


              
              