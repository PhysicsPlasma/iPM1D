 Module FileIO
    Use TypeModule
    Use ParticleModule  
    Implicit none
    Integer(4),Parameter  ::  DefaultNameIndex=10000
    Integer(4),Parameter  ::  DefaultNameIndexInit=20000
    Integer(4),Parameter  ::  ModeMultiplier=1000
    
    contains 
    subroutine DumpParticle(PB,Mode)
        implicit none
        Type(ParticleBundle),intent(inout) :: PB
        Type(ParticleBundle) :: TempParticleBundle
        Integer(4),intent(in)  :: Mode
        Character(len=99) :: Filename
        Integer(4) :: i,Dim,NameIndex
        Real(8) :: XFactor,YFactor
        Integer(4),save :: NameTimer=1
        If (Mode==0) Then
                   NameIndex=DefaultNameIndex
        else
                   NameIndex=DefaultNameIndexInit+ModeMultiplier*Mode+NameTimer
                   NameTimer=NameTimer+1
        End If 
        
        TempParticleBundle=PB
        Call UpdatePositionBundle(TempParticleBundle,1_4)
        Call UpdateVelocityBundle(TempParticleBundle,1_4)
        Write(filename,*) NameIndex,trim(TempParticleBundle%Name),".dat"
        Write(*,*) 'Saving ', trim(TempParticleBundle%Name), TempParticleBundle%NPar,' Please Wait...' 
        open (10,file=filename)
        Write(10,*) TempParticleBundle%NPar
        Write(10,*) TempParticleBundle%Charge,TempParticleBundle%Mass,TempParticleBundle%Weight
 !       Write(10,*) 
        Dim=Sizeof(TempParticleBundle%PhaseSpace(1))/8_4
        do i=1,TempParticleBundle%NPar
             Write(10,FMt="(<Dim>D23.15)")  TempParticleBundle%PhaseSpace(i)
        end do
        close(10)
        Write(*,*) 'Saving ', trim(TempParticleBundle%Name),' Complete!'  
        return
    end subroutine  DumpParticle
   
    subroutine LoadParticle(InputSpecy,PB,Status)
        implicit none
        Type(OneSpecy),intent(in) ::  InputSpecy
        Type(ParticleBundle),intent(inout) :: PB
        Logical,intent(inout) ::  Status
        Logical :: alive
        Character(len=99) :: Filename
        Integer(4) :: i,NameIndex
        NameIndex=DefaultNameIndex
        Write(filename,*) NameIndex,trim(InputSpecy%Name),".dat"
        Inquire(file=filename,exist=alive)
        If(alive) then
               Open (10,file=filename)
               PB%Name=InputSpecy%Name
               Read(10,*) PB%NPar
               Read(10,*)  PB%Charge, PB%Mass, PB%Weight
               Write(*,*) 'Loading ', trim(PB%Name), PB%NPar,' Please Wait...' 
               Write(filename,*) trim(InputSpecy%Name),".dat"
               do i=1,PB%NPar
                     Read(10,*)  PB%PhaseSpace(i)
               end do
               Status=.False.
               Write(*,*) 'Loading ', trim(PB%Name),' Complete!'   
         else
               Write(*,*) 'Can not find the file for ', trim(InputSpecy%Name),' the particle will be randomly initilalized.'    
               Status=.True. 
         End if        
        return
    end subroutine  LoadParticle
         
    subroutine DumpField(FG,Mode)
        Implicit none
        Type(Field),intent(inout) :: FG
        Integer(4) :: Mode
        Character(len=99) :: Filename
        Integer(4):: i,j,NameIndex
        Integer(4),save :: NameTimer=1
        If (Mode==0) Then
                   NameIndex=DefaultNameIndex
        else
                   NameIndex=DefaultNameIndexInit+ModeMultiplier*Mode+NameTimer
                   NameTimer=NameTimer+1
        End If 
        
        Write(filename,*) NameIndex,"FieldGlobal",".dat"
        Write(*,*) "Saving ",filename," Please wait..."
        open (10,file=filename)
        Write(10,*)  FG%Nx,FG%Dx,FG%Dt
        do i=1,FG%Nx
               Write(10,FMt="(5D23.15)") dble(i-1)*FG%Dx,FG%Ex(i),FG%Phi(i),FG%Rho(i),FG%Chi(i)!,FG%Bx(i),FG%By(i)
        end do
        close(10)
        Write(*,*) "Save ",filename,"Complete!"  
        return
     end subroutine DumpField
    ! 
     subroutine LoadField(FG,Status)
        Implicit none
        Type(Field),intent(inout) :: FG
        Logical,intent(inout) ::  Status
        Character(len=99) :: Filename
        logical :: alive
        Integer(4) :: i,j,NameIndex
        Real(8) :: TempDx
        NameIndex=DefaultNameIndex
        Write(filename,*) NameIndex,"FieldGlobal",".dat"
        
        Inquire(file=filename,exist=alive)
        If(alive) then
                    Write(*,*) "Loading ",filename," Please wait..."
                    open (10,file=filename)
                    Read(10,*)  FG%Nx,FG%Dx,FG%Dt
                    do i=1,FG%Nx
                           Read(10,*) TempDx,FG%Ex(i),FG%Phi(i),FG%Rho(i),FG%Chi(i)!,FG%Bx(i),FG%By(i)
                    end do
                    close(10)
                    Status=.False. 
                    Write(*,*) "Load ",filename,"Complete!" 
         else
                    Write(*,*) 'Can not find the file for ', filename,' set for zero.'    
                    FG%Ex=0.d0
                    FG%Phi=0.d0
                    FG%Rho=0.d0
                    FG%Chi=0.d0
                    Status=.True. 
         End if              
        return
     end subroutine LoadField
     
     subroutine GridDump(GD,Mode)
        Implicit none
        Class(*),intent(inout)  :: GD
        Integer(4),intent(in) :: Mode
        Character(len=99) :: Filename
        Integer(4) :: i,j,k,NameIndex,Ns
        Integer(4),save :: NameTimer=1
        If (Mode==0) Then
                   NameIndex=DefaultNameIndex
        else
                   NameIndex=DefaultNameIndexInit+ModeMultiplier*Mode+NameTimer
                   NameTimer=NameTimer+1
        End If 
        Select Type (GD)
            Type is (Grid1D(*,*))
                Write(filename,*) "Grid1D",NameIndex,".dat"
                Write(*,*) "Saving ",Filename," Please wait..."
                open (10,file=filename)
                Ns=GD%Ns
                do i=1,GD%Nx
                       Write(10,FMt="(<1+Ns>D23.15)")  dble(i-1)*GD%dx,(GD%Value(i,k),k=1,Ns)
                end do
                close(10)
                Write(*,*) "Save ",Filename,"Complete!"
            Type is (Grid2D(*,*,*))
                Write(filename,*) "Grid2D",NameIndex,".dat"
                Filename=Trim(filename)
                Write(*,*) "Saving ",Filename," Please wait..."
                open (10,file=filename)
                Ns=GD%Ns
                do i=0,GD%Ny
                     do j=1,GD%Nx
                       Write(10,FMt="(<2+Ns>D23.15)")  dble(j-1)*GD%dx,dble(i-1)*GD%dy,(GD%Value(j,i,k),k=1,Ns)
                     End do
                end do
                close(10)
                Write(*,*) "Save ",Filename,"Complete!"
            End Select
        return
     end subroutine GridDump
     
     subroutine GridInitialization(GD,Period,Dx,Dy)
        Implicit none
        Class(*),intent(inout)  :: GD
        Integer(4),intent(in) :: Period
        Real(8),intent(in) :: Dx,Dy
        Select Type (GD)
            Type is (Grid1D(*,*))
                    GD%Timer=0
                    GD%Period=Period
                    GD%Dx=Dx
                    GD%Value=0.d0
            Type is (Grid2D(*,*,*))
                    GD%Timer=0
                    GD%Period=Period
                    GD%Dx=Dx
                    GD%Dy=2.d0/Dble(GD%Ny)!FG%Dt
                    GD%Value=0.d0
            End Select
        return
     end subroutine GridInitialization
     
              subroutine DumpFieldSolver(FS,Mode)
        Implicit none
        Type(FieldSolver),intent(inout) :: FS
        Integer(4) :: Mode
        Character(len=99) :: Filename
        Integer(4):: i,j,NameIndex
        Integer(4),save :: NameTimer=1
        If (Mode==0) Then
                   NameIndex=DefaultNameIndex
        else
                   NameIndex=DefaultNameIndexInit+ModeMultiplier*Mode+NameTimer
                   NameTimer=NameTimer+1
        End If 
        
        Write(filename,*) NameIndex,"FieldSolver",".dat"
        Write(*,*) "Saving ",filename," Please wait..."
        open (10,file=filename)
        Write(10,*)  FS%Ns,FS%Dx,FS%Dt
        do i=1,FS%Ns
               Write(10,FMt="(6D23.15)") dble(i-1)*FS%Dx,FS%Solve(i),FS%Source(i),FS%CoeA(i),FS%CoeB(i),FS%CoeC(i)!,FG%Bx(i),FG%By(i)
        end do
        close(10)
        Write(*,*) "Save ",filename,"Complete!"  
        return
     end subroutine DumpFieldSolver
     
   END  Module FileIO
        !subroutine DumpFieldOne(SpecyIndex,FO)
    !    Implicit none
    !    Type(Field),intent(in) :: FO
    !    Integer(4),intent(in) :: SpecyIndex
    !    Character(len=99) :: Filename
    !    Integer(4):: i,j
    !    write(filename,*) SpecyIndex,"Particle.dat"
    !    Write(*,*) "Saving ",filename," Please wait..."
    !    open (10,file=filename)
    !    do i=1,FO%Nx
    !           Write(10,FMt="(4D23.15)")  i*FO%dx,FO%RhoOne(i),FO%JxOne(i),FO%TOne(i)
    !    end do
    !    close(10)
    !    Write(*,*) "Save ",filename,"Complete!"  
    !    return
    !end subroutine DumpFieldOne