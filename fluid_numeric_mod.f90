!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fluid_numeric_mod.f90   Bengt Eliasson     2017-07-22     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains numerical objects that can act on    !
! physical objects from the class TwoDimFluidPlasma, or on !
! parts of the object.                                      !
!                                                           !
! The module contains the following numerical objects:      !
!                                                           !
! Domain_Initializer - Calculates the initial condition for !
!   the problem, calculates various constants, opens data   !
!   files, et.c.                                            !
!                                                           !
! Domain_Finalizer - Closes data files, et.c.               !
!                                                           !
! Domain_StatisticsWriter - Writes various statistics of    !
!   the solution to file and screen.                        !
!                                                           !
! Domain_StatisticsCalculator - Calculates statistics from  !
!   the solution.                                           !
!                                                           !
! Fluid_rhs - Calculates the righthand side of the time-   !
!   dependent partial differential eqation.                 !
!                                                           !
! ElectrStaticField - Calculates the electric field by      !
!   integrating the Poisson equation (Maxwell's equations). !
!                                                           !
! x1_Differentiator - Calculates d/dx of a two-dimensional  !
!   function f(x,eta), using the fourth order difference    ! 
!   approximation.                                          !
!                                                           !
! x1_FourierDifferentiator - Calculates d/dx of a           !
!   two-dimensional function f(x,eta), using the            !
!   Fourier transform.                                      !
!                                                           !
! eta1_Differetiator - Calculates d/deta of a two-          !
!   dimensional function f(x,eta), using the Pade'          ! 
!   approximation.                                          !
!                                                           !
! SimpsonIntegrator1D - integrates a function f(x),         !
!   using the Simpson formula.                              !
!                                                           !
! RungeKutta_TimeStepper - Advances the solution one time   !
!   step by using the fouth order Runge-Kutta method.       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module fluid_numeric_mod
  USE fluid_domain_mod

  IMPLICIT NONE


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the initial condition, various constants, et.c.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Domain_Initializer(D1,t)
  IMPLICIT NONE
  TYPE (TwoDimFluidPlasma),INTENT(inout) :: D1
  REAL(8),INTENT(inout) :: t

  INTEGER :: i1,i2
  REAL(8) :: x,y,ni02

  ! The grid size in x direction
  D1%dx = (x_end-x_start)/DBLE(TotalNx)
  D1%dy = (y_end-y_start)/DBLE(TotalNy)

  ! Time starts at zero
  t=zero

  ! Initialize MPI
  CALL MPI_INIT(errcode)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,errcode)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,errcode)

  cartsize(1)=dim_x
  cartsize(2)=dim_y
  
  period=(/ .TRUE., .TRUE./)
  call MPI_CART_CREATE(MPI_COMM_WORLD,CART_DIM,cartsize, &
     period,reorder,cart_comm,errcode)

  call MPI_CART_COORDS(cart_comm,rank,cart_dim,mycoord,errcode)
  call MPI_CART_SHIFT(cart_comm,SIDEWAYS,1,link(LEFT),link(RIGHT),errcode)
  call MPI_CART_SHIFT(cart_comm,UPDOWN  ,1,link(DOWN)  ,link(UP) ,errcode)

  my_x=mycoord(1)
  my_y=mycoord(2)

  ! Create a communicator for my row.
  CALL MPI_COMM_SPLIT(cart_comm, my_y, my_x, row_comm, errcode)

  ! Create a communicator for my column.
  CALL MPI_COMM_SPLIT(cart_comm, my_x, my_y, col_comm, errcode)

  IF (my_y.EQ.0) THEN
    my_down=dim_y-1
  ELSE
    my_down=my_y-1
  END IF

  IF (my_y.EQ.dim_y-1) THEN
    my_up=0
  ELSE
    my_up=my_y+1
  END IF

  IF (my_x.EQ.0) THEN
    my_left=dim_x-1
  ELSE
    my_left=my_x-1
  END IF

  IF (my_x.EQ.dim_x -1) THEN
    my_right=0
  ELSE
    my_right=my_x+1
  END IF

  ! The known variables x and eta
  DO i1=1,Nx
    D1%x(i1) = x_start+DBLE(my_x*Nx+i1-1)*D1%dx
  END DO

  DO i2=1,Ny
    D1%y(i2) = y_start+DBLE(my_y*Ny+i2-1)*D1%dy
  END DO

  ! Mass and charge
  D1%ions1%mass=Ion1_mass
  D1%ions1%charge=Ion1_charge
  D1%ions1%temperature=Ion1_temperature
  D1%ions2%mass=Ion2_mass
  D1%ions2%charge=Ion2_charge
  D1%ions2%temperature=Ion2_temperature

  ! Initial values
  DO i2 = 1,Ny
    y = D1%y(i2)
    DO i1 = 1,Nx
      x = D1%x(i1)

      D1%ions1%ni(i1,i2)=ni1_initial(x,y)
      D1%ions1%vix(i1,i2)=vi1x_initial(x,y)
      D1%ions1%viy(i1,i2)=vi1y_initial(x,y)
      D1%ions2%ni(i1,i2)=ni2_initial(x,y)
      D1%ions2%vix(i1,i2)=vi2x_initial(x,y)
      D1%ions2%viy(i1,i2)=vi2y_initial(x,y)
      D1%Bz(i1,i2)=Bz_initial(x,y)
    END DO
  END DO


  CALL PeriodicIntegrator2D(D1%ions2%ni,D1%dx,D1%dy,Nx,Ny,ni02)
  Ni02=NParticles/(Ni02*two*R0)
  DO i2 = 1,Ny
  DO i1 = 1,Nx
    D1%ions2%ni(i1,i2)=D1%ions2%ni(i1,i2)*Ni02
  END DO
  END DO


  D1%v_parallel=v_parallel
  D1%acc_time=acc_time
  D1%R0=R0
  
  DO i2 = 1,Ny
    y = D1%y(i2)
    DO i1 = 1,Nx
      x = D1%x(i1)

      D1%v0x(i1,i2)=v_radial*x/DSQRT(R0**2+(x**2+y**2)/four)+v_satellite
      D1%v0y(i1,i2)=v_radial*y/DSQRT(R0**2+(x**2+y**2)/four)
    END DO
  END DO

  ! Maximum frequency, used to calculate dt.
  D1%nu_x_max=zero
  D1%nu_y_max=zero

END SUBROUTINE Domain_Initializer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Close files et.c.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Domain_Finalizer(domain)
  IMPLICIT NONE
  TYPE (TwoDimFluidPlasma),INTENT(inout) :: domain

END SUBROUTINE Domain_Finalizer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write results to screen and to data files 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Domain_StatisticsWriter(domain,k,t,dt)
  IMPLICIT NONE

  ! Arguments
  TYPE (TwoDimFluidPlasma),INTENT(in) :: domain
  INTEGER, INTENT(in) :: k
  REAL(8), INTENT(in) :: t,dt
  CHARACTER(len=4) :: str
  INTEGER :: iocheck
  REAL(8) :: start,finish1,finish2

  !!!! MPI
  REAL(8),DIMENSION(TotalNx*TotalNy) :: DataGlobal
  !INTEGER :: NProcX,NProcY
  !REAL(8),DIMENSION(:,:),ALLOCATABLE :: BzList
  !INTEGER :: Nproc,NRowsLocal,NColsLocal
  INTEGER :: LocArraySize

  !NProc=NProcX*NProcY
  !NRowsLocal=NRowsGlobal/NProcY
  !NColsLocal=NColsGlobal/NProcX

  LocArraySize=Nx*Ny
  !ALLOCATE(BzList(LocArraySize,NProc))


  ! Write dumpnr to str
  WRITE(str,1000)domain%dumpnr

  1000 FORMAT(I4)

  IF (NP.GE.7) THEN
    IF (write_Bz) THEN
      CALL MPI_GATHER(domain%Bz,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)
    END IF

    IF (write_ni1) THEN
      CALL MPI_GATHER(domain%ions1%ni,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,1,MPI_COMM_WORLD,errcode)
    END IF

    IF (write_vi1x) THEN  
      CALL MPI_GATHER(domain%ions1%vix,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,2,MPI_COMM_WORLD,errcode)
    END IF

    IF (write_vi1y) THEN
      CALL MPI_GATHER(domain%ions1%viy,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,3,MPI_COMM_WORLD,errcode)
    END IF

    IF (write_ni2) THEN
      CALL MPI_GATHER(domain%ions2%ni,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,4,MPI_COMM_WORLD,errcode)
    END IF

    IF (write_vi2x) THEN
      CALL MPI_GATHER(domain%ions2%vix,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,5,MPI_COMM_WORLD,errcode)
    END IF

    IF (write_vi2y) THEN
      CALL MPI_GATHER(domain%ions2%viy,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,6,MPI_COMM_WORLD,errcode)
    END IF


    IF (write_Bz) THEN
      IF (rank.EQ.root) THEN
        OPEN (UNIT=out_Bz,FILE=FileDir//'Bz_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)

        
        WRITE (out_Bz,*) DataGlobal

        CLOSE(out_Bz)
      END IF
    END IF

    IF (write_ni1) THEN

      IF (rank.EQ.1) THEN
        OPEN (UNIT=out_ni1,FILE=FileDir//'ni1_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)
        WRITE (out_ni1,*) DataGlobal

        CLOSE(out_ni1)
      END IF
    END IF

    IF (write_vi1x) THEN
            
      IF (rank.EQ.2) THEN
        OPEN (UNIT=out_vi1x,FILE=FileDir//'vi1x_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)
 
        WRITE (out_vi1x,*) DataGlobal

        CLOSE(out_vi1x)
      END IF
    END IF

    IF (write_vi1y) THEN

      IF (rank.EQ.3) THEN
        OPEN (UNIT=out_vi1y,FILE=FileDir//'vi1y_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)

        WRITE (out_vi1y,*) DataGlobal

        CLOSE(out_vi1y)
      END IF
    END IF

    IF (write_ni2) THEN

      IF (rank.EQ.4) THEN
        OPEN (UNIT=out_ni2,FILE=FileDir//'ni2_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)

        WRITE (out_ni2,*) DataGlobal

        CLOSE(out_ni2)
      END IF
    END IF

    IF (write_vi2x) THEN

      IF (rank.EQ.5) THEN
        OPEN (UNIT=out_vi2x,FILE=FileDir//'vi2x_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)

        WRITE (out_vi2x,*) DataGlobal

        CLOSE(out_vi2x)
      END IF
    END IF

    IF (write_vi2y) THEN

      IF (rank.EQ.6) THEN
        OPEN (UNIT=out_vi2y,FILE=FileDir//'vi2y_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)

        WRITE (out_vi2y,*) DataGlobal

        CLOSE(out_vi2y)
      END IF
    END IF
  ELSE

    IF (write_Bz) THEN
      CALL MPI_GATHER(domain%Bz,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)
      IF (rank.EQ.root) THEN
        OPEN (UNIT=out_Bz,FILE=FileDir//'Bz_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)

        
        WRITE (out_Bz,*) DataGlobal

        CLOSE(out_Bz)
      END IF
    END IF

    IF (write_ni1) THEN
      CALL MPI_GATHER(domain%ions1%ni,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)
      IF (rank.EQ.root) THEN
        OPEN (UNIT=out_ni1,FILE=FileDir//'ni1_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)
        WRITE (out_ni1,*) DataGlobal

        CLOSE(out_ni1)
      END IF
    END IF

    IF (write_vi1x) THEN
      CALL MPI_GATHER(domain%ions1%vix,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)            
      IF (rank.EQ.root) THEN
        OPEN (UNIT=out_vi1x,FILE=FileDir//'vi1x_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)
 
        WRITE (out_vi1x,*) DataGlobal

        CLOSE(out_vi1x)
      END IF
    END IF

    IF (write_vi1y) THEN
      CALL MPI_GATHER(domain%ions1%viy,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)
      IF (rank.EQ.root) THEN
        OPEN (UNIT=out_vi1y,FILE=FileDir//'vi1y_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)

        WRITE (out_vi1y,*) DataGlobal

        CLOSE(out_vi1y)
      END IF
    END IF

    IF (write_ni2) THEN
      CALL MPI_GATHER(domain%ions2%ni,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)
      IF (rank.EQ.root) THEN
        OPEN (UNIT=out_ni2,FILE=FileDir//'ni2_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)

        WRITE (out_ni2,*) DataGlobal

        CLOSE(out_ni2)
      END IF
    END IF

    IF (write_vi2x) THEN
      CALL MPI_GATHER(domain%ions2%vix,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)
      IF (rank.EQ.root) THEN
        OPEN (UNIT=out_vi2x,FILE=FileDir//'vi2x_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)

        WRITE (out_vi2x,*) DataGlobal

        CLOSE(out_vi2x)
      END IF
    END IF

    IF (write_vi2y) THEN
      CALL MPI_GATHER(domain%ions2%viy,LocArraySize,MPI_DOUBLE_PRECISION,DataGlobal, &
          LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)

      IF (rank.EQ.root) THEN
        OPEN (UNIT=out_vi2y,FILE=FileDir//'vi2y_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)

        WRITE (out_vi2y,*) DataGlobal

        CLOSE(out_vi2y)
      END IF
    END IF

  END IF

  IF (rank.EQ.root) THEN
    WRITE(*,906)k,' dt=',dt,' t=',t
  END IF


906 FORMAT (I7,A,F15.10,A,F15.10)

END SUBROUTINE Domain_StatisticsWriter


SUBROUTINE TimeStepCalculator(domain,dt)
  IMPLICIT NONE
  TYPE (TwoDimFluidPlasma),INTENT(inout) :: domain
  REAL(8),INTENT(out) :: dt
  REAL(8),DIMENSION(Nx,Ny) :: dt_arr
  REAL(8),DIMENSION(Nx,Ny) :: ne,inv_ne,vi_max,d1x_inv_ne,d1y_inv_ne
  REAL(8),DIMENSION(Nx,Ny) :: vdy_max,vdx_max,vdy_max_minus,vdx_max_minus
  REAL(8) :: inv_mu0_e,mu0_mi,inv_dx,inv_dy,temp
  INTEGER :: i1,i2

  ! MPI
  INTEGER :: errcode

  IF (IsAdaptive) THEN
    IF ((domain%nu_x_max.GT.zero).AND.(domain%nu_y_max.GT.zero)) THEN
            
      temp=domain%nu_x_max
      CALL MPI_REDUCE(temp,domain%nu_x_max,1,  &
        MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,errcode)

      temp=domain%nu_y_max
      CALL MPI_REDUCE(temp,domain%nu_y_max,1,  &
        MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,errcode)

      IF (rank.EQ.root) THEN
        dt=CFL*SqrtEight/(domain%nu_x_max/domain%dx+domain%nu_y_max/domain%dy)
      END IF
      
      ! Broadcast the time step to all procs.
      CALL MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION,root, MPI_COMM_WORLD, errcode)
      
      !dt=CFL/domain%nu_max
    ELSE
      inv_mu0_e=one/(mu0*domain%ions1%charge)
      mu0_mi=mu0*domain%ions1%mass
      inv_dx=one/domain%dx
      inv_dy=one/domain%dy

      DO i2=1,Ny
      DO i1=1,Nx
        ne(i1,i2)=MAX(domain%ions1%ni(i1,i2)+domain%ions2%ni(i1,i2),2.0d11)
        inv_ne(i1,i2)=one/ne(i1,i2)
        vi_max(i1,i2)=max(sqrt(domain%ions1%vix(i1,i2)**2+domain%ions1%viy(i1,i2)**2),  &
              sqrt(domain%ions2%vix(i1,i2)**2+domain%ions2%viy(i1,i2)**2));
      END DO
      END DO

      CALL d1x_plus(inv_ne,d1x_inv_ne,domain%dx)    
      CALL d1y_plus(inv_ne,d1y_inv_ne,domain%dy)    

      vdy_max(:,:)=abs(domain%Bz(:,:)*d1x_inv_ne(:,:))*inv_mu0_e
      vdx_max(:,:)=abs(domain%Bz(:,:)*d1y_inv_ne(:,:))*inv_mu0_e

      temp=(inv_dx+inv_dy)
      dt_arr(:,:)=(vi_max(:,:)+ABS(domain%Bz(:,:))*SQRT(inv_ne(:,:))/SQRT(mu0_mi))*temp &
            +vdy_max(:,:)*inv_dy+vdx_max(:,:)*inv_dx

      dt_arr(:,:)=CFL*sqrt(eight)/dt_arr(:,:);
      !dt_arr(:,:)=CFL/dt_arr(:,:);
      dt=MINVAL(dt_arr)

      ! Calculate the minimum dt among processors
      temp=dt
      CALL MPI_REDUCE(temp,dt,1,  &
        MPI_DOUBLE_PRECISION,MPI_MIN,root,MPI_COMM_WORLD,errcode)

      ! Broadcast the time step to all procs.
      CALL MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION,root, MPI_COMM_WORLD, errcode)

    END IF


  ELSE
    dt=TimeStep
    dt_arr(:,:)=TimeStep
  END IF
END SUBROUTINE TimeStepCalculator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Runge-Kutta time stepper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RungeKutta_TimeStepper(domain,t,dt,direction,vx_frame_global)
  IMPLICIT NONE

  ! Arguments
  TYPE(TwoDimFluidPlasma),INTENT(inout) :: domain
  REAL(8),INTENT(inout) :: t
  REAL(8),INTENT(in) :: dt,vx_frame_global
  CHARACTER, INTENT(in) :: direction

  !
  TYPE(TwoDimFluidPlasma) :: D1,D2,RHS_D

  REAL(8) :: frac_dt_2, frac_dt_3, frac_dt_6

  frac_dt_2=dt/two
  frac_dt_3=dt/three
  frac_dt_6=dt/six


  !----

  CALL Copy_Domain(domain,D1)

  !----
  
  CALL RHS(domain,RHS_D,t,direction,vx_frame_global)

  !----

  CALL Add_domain(RHS_D,domain,frac_dt_6)

  CALL Copy_domain(D1,D2)
  CALL Add_domain(RHS_D,D2,frac_dt_2)

  !----

  t=t+frac_dt_2
  
  CALL RHS(D2,RHS_D,t,direction,vx_frame_global)

  !----

  CALL Add_domain(RHS_D,domain,frac_dt_3)
 
  CALL Copy_domain(D1,D2)
  CALL Add_domain(RHS_D,D2,frac_dt_2)

  !----

  CALL RHS(D2,RHS_D,t,direction,vx_frame_global)

  !----

  CALL Add_domain(RHS_D,domain,frac_dt_3)
 
  CALL Copy_domain(D1,D2)
  CALL Add_domain(RHS_D,D2,dt)

  !----

  t=t+frac_dt_2

  CALL RHS(D2,RHS_D,t,direction,vx_frame_global)

  !----

  CALL Add_domain(RHS_D,domain,frac_dt_6)

  IF (direction.EQ.'x') THEN
    domain%nu_x_max=RHS_D%nu_x_max
  ELSE
    domain%nu_y_max=RHS_D%nu_y_max
  END IF
END SUBROUTINE RungeKutta_TimeStepper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Euler time stepper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Euler_TimeStepper(domain,t,dt,direction,vx_frame_global)
  IMPLICIT NONE

  ! Arguments
  TYPE(TwoDimFluidPlasma),INTENT(inout) :: domain
  REAL(8),INTENT(inout) :: t
  REAL(8),INTENT(in) :: dt,vx_frame_global
  CHARACTER, INTENT(in) :: direction

  !
  TYPE(TwoDimFluidPlasma) :: RHS_D

!  REAL(8) :: frac_dt_2, frac_dt_3, frac_dt_6

!  frac_dt_2=dt/two
!  frac_dt_3=dt/three
!  frac_dt_6=dt/six


  !----

!  CALL Copy_Domain(domain,D1)

  !----
  
  CALL RHS(domain,RHS_D,t,direction,vx_frame_global)

  !----

  CALL Add_domain(RHS_D,domain,dt)

  IF (direction.EQ.'x') THEN
    domain%nu_x_max=RHS_D%nu_x_max
  ELSE
    domain%nu_y_max=RHS_D%nu_y_max
  END IF

END SUBROUTINE Euler_TimeStepper


SUBROUTINE Copy_domain(D1,D2)
  TYPE(TwoDimFluidPlasma), INTENT(in) :: D1
  TYPE(TwoDimFluidPlasma), INTENT(out) :: D2

  INTEGER i1,i2


  ! The x variable
  DO i1=1,Nx
    D2%x(i1)=D1%x(i1)
  END DO
  DO i2=1,Ny
    D2%y(i2)=D1%y(i2)
  END DO

    
  ! The step sizes in x direction
  D2%dx=D1%dx
  D2%dy=D1%dy

  DO i2=1,Ny
  DO i1=1,Nx
    D2%Bz(i1,i2)=D1%Bz(i1,i2)
  END DO
  END DO

  ! Copy the particle distribution functions
  CALL Copy_species(D1%ions1,D2%ions1)
  CALL Copy_species(D1%ions2,D2%ions2)

  ! Copy initial velocities etc.
  D2%acc_time=D1%acc_time
  D2%R0=D1%R0
  D2%v_parallel=D1%v_parallel

  ! Copy maximum frequency, used to calculate dt
  D2%nu_x_max=D1%nu_x_max
  D2%nu_y_max=D1%nu_y_max
  
  DO i2=1,Ny
  DO i1=1,Nx
    D2%v0x(i1,i2)=D1%v0x(i1,i2)
    D2%v0y(i1,i2)=D1%v0y(i1,i2)
  END DO
  END DO


END SUBROUTINE Copy_domain


SUBROUTINE Copy_species(S1,S2)
  TYPE(TwoDimFluidSpecies), INTENT(in) :: S1
  TYPE(TwoDimFluidSpecies), INTENT(out) :: S2

  INTEGER :: i1,i2

  ! Copy the particle density and particle current density
  DO i2=1,Ny
  DO i1=1,Nx
    S2%ni(i1,i2)=S1%ni(i1,i2)
    S2%vix(i1,i2)=S1%vix(i1,i2)
    S2%viy(i1,i2)=S1%viy(i1,i2)
  END DO
  END DO

  ! Copy the charge, mass and temperature
  S2%charge=S1%charge
  S2%mass=S1%mass
  S2%temperature=S1%temperature

END SUBROUTINE Copy_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Add some of the fields in two domains, D2 <-- D2+alpha*D1   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Add_domain(D1,D2,alpha)
  TYPE(TwoDimFluidPlasma), INTENT(in) :: D1
  TYPE(TwoDimFluidPlasma), INTENT(inout) :: D2
  REAL(8), INTENT(in) :: alpha

  INTEGER i1,i2

  ! Add the particle densities and velocities
  CALL Add_species(D1%ions1,D2%ions1,alpha)
  CALL Add_species(D1%ions2,D2%ions2,alpha)

  ! Add magnetic field
  DO i2=1,Ny
  DO i1=1,Nx
    D2%Bz(i1,i2)=D2%Bz(i1,i2)+D1%Bz(i1,i2)*alpha
  END DO
  END DO
  
  
END SUBROUTINE Add_domain


SUBROUTINE Add_species(S1,S2,alpha)
  TYPE(TwoDimFluidSpecies), INTENT(in) :: S1
  TYPE(TwoDimFluidSpecies), INTENT(inout) :: S2
  REAL(8), INTENT(in) :: alpha

  INTEGER :: i1,i2

  DO i2=1,Ny
  DO i1=1,Nx
    S2%ni(i1,i2)=S2%ni(i1,i2)+S1%ni(i1,i2)*alpha
    S2%vix(i1,i2)=S2%vix(i1,i2)+S1%vix(i1,i2)*alpha
    S2%viy(i1,i2)=S2%viy(i1,i2)+S1%viy(i1,i2)*alpha
  END DO
  END DO
END SUBROUTINE Add_species


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the righthand side P of the differential 
! equation df/dt=P
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RHS (domain,RHS_D,t,direction,vx_frame_global)
  IMPLICIT NONE

  ! Arguments
  TYPE (TwoDimFluidPlasma),INTENT(in) :: domain
  TYPE (TwoDimFluidPlasma),INTENT(out) :: RHS_D
  REAL(8),INTENT(in) :: t,vx_frame_global
  CHARACTER, INTENT(in) :: direction

  REAL(8), DIMENSION(1:Nx,1:Ny) :: vi_sq,dvi_sq,dvi,dd,dd2,dissip
  REAL(8), DIMENSION(1:Nx,1:Ny) :: ni_vix,d1x_ni_vix,ni_viy,d1y_ni_viy
  REAL(8), DIMENSION(1:Nx,1:Ny) :: inv_ne,ni1_Bz_ne,ni2_Bz_ne,ve_Bz
  REAL(8), DIMENSION(1:Nx,1:Ny) :: diff_ve_Bz,B2_2mu0,d1x_B2_2mu0,d1y_B2_2mu0
  REAL(8) :: charge_mass, inv_m, inv_e, temp, diss_dd,dx,dy

  REAL(8), DIMENSION(1:Nx,1:Ny) :: dni1_d_plus,dvi1x_d_plus,dvi1y_d_plus
  REAL(8), DIMENSION(1:Nx,1:Ny) :: dni2_d_plus,dvi2x_d_plus,dvi2y_d_plus
  REAL(8), DIMENSION(1:Nx,1:Ny) :: dBz_d_plus
  REAL(8), DIMENSION(1:Nx,1:Ny) :: dni1_d_minus,dvi1x_d_minus,dvi1y_d_minus
  REAL(8), DIMENSION(1:Nx,1:Ny) :: dni2_d_minus,dvi2x_d_minus,dvi2y_d_minus
  REAL(8), DIMENSION(1:Nx,1:Ny) :: dBz_d_minus
  REAL(8), DIMENSION(1:Nx,1:Ny) :: dinv_ne_d
  REAL(8), DIMENSION(1:Nx,1:Ny) :: dni1_d2,dvi1x_d2,dvi1y_d2
  REAL(8), DIMENSION(1:Nx,1:Ny) :: dni2_d2,dvi2x_d2,dvi2y_d2
  REAL(8), DIMENSION(1:Nx,1:Ny) :: dBz_d2  
  REAL(8) :: nu_x_max,nu_y_max

  !!! Upstream !!!!!
  INTEGER, PARAMETER :: lwork=4384

  !!! 5x5 system !!!
  INTEGER, PARAMETER :: lda=5
  REAL(8), DIMENSION(1:lda) :: d1_plus,d1_minus,wr,wi,d2
  REAL(8), DIMENSION(1:lda) :: d1,d1_out
  COMPLEX(8), DIMENSION(1:lda) :: d1_C,d1_out_C
  REAL(8), DIMENSION(1:lda,1:3) :: d1_minus_plus
  COMPLEX(8), DIMENSION(1:lda,1:3) :: d1_minus_plus_C
  REAL(8), DIMENSION(1:lda,1:lda) :: matrix_A
  REAL(8), DIMENSION(1:lda,1:lda) :: vr,vr_copy
  COMPLEX(8), DIMENSION(1:lda,1:lda) :: vr_C,vr_C_copy
  INTEGER, DIMENSION(1:lda) :: ipvt
  !!! 2x2 system !!!
  INTEGER, PARAMETER :: lda_2x2=2
  REAL(8), DIMENSION(1:lda_2x2) :: d1_plus_2x2,d1_minus_2x2,wr_2x2,wi_2x2,d2_2x2
  REAL(8), DIMENSION(1:lda_2x2) :: d1_2x2,d1_out_2x2
  COMPLEX(8), DIMENSION(1:lda_2x2) :: d1_C_2x2,d1_out_C_2x2
  REAL(8), DIMENSION(1:lda_2x2,1:3) :: d1_minus_plus_2x2
  COMPLEX(8), DIMENSION(1:lda_2x2,1:3) :: d1_minus_plus_C_2x2
  REAL(8), DIMENSION(1:lda_2x2,1:lda_2x2) :: matrix_A_2x2
  REAL(8), DIMENSION(1:lda_2x2,1:lda_2x2) :: vr_2x2,vr_copy_2x2
  COMPLEX(8), DIMENSION(1:lda_2x2,1:lda_2x2) :: vr_C_2x2,vr_C_copy_2x2
  INTEGER, DIMENSION(1:lda_2x2) :: ipvt_2x2
   !!! 3x3 system !!!
  INTEGER, PARAMETER :: lda_3x3=3
  REAL(8), DIMENSION(1:lda_3x3) :: d1_plus_3x3,d1_minus_3x3,wr_3x3,wi_3x3,d2_3x3
  REAL(8), DIMENSION(1:lda_3x3) :: d1_3x3,d1_out_3x3
  COMPLEX(8), DIMENSION(1:lda_3x3) :: d1_C_3x3,d1_out_C_3x3
  REAL(8), DIMENSION(1:lda_3x3,1:3) :: d1_minus_plus_3x3
  COMPLEX(8), DIMENSION(1:lda_3x3,1:3) :: d1_minus_plus_C_3x3
  REAL(8), DIMENSION(1:lda_3x3,1:lda_3x3) :: matrix_A_3x3
  REAL(8), DIMENSION(1:lda_3x3,1:lda_3x3) :: vr_3x3,vr_copy_3x3
  COMPLEX(8), DIMENSION(1:lda_3x3,1:lda_3x3) :: vr_C_3x3,vr_C_copy_3x3
  INTEGER, DIMENSION(1:lda_3x3) :: ipvt_3x3
  !!!!!!!!!!!!!!!!!

  REAL(8) :: vl
  INTEGER :: i1,i2,i3,i4,ierr,ldvr,ldvl
  REAL(8), DIMENSION(1:lwork) :: work

  REAL(8) :: v_min,n_min,vi1x,vi2x

  ! MPI
  INTEGER :: root
  INTEGER :: request_send,request_recv
  INTEGER :: status(MPI_STATUS_SIZE),ierror
  INTEGER :: sendtag,recvtag
  REAL(8) :: sendbuf(1:Ny),sendbuf_all_r(1:Ny,1:7),sendbuf_all_u(1:Nx,1:7)
  REAL(8) :: lc(1:Ny),lc_all(1:Ny,1:7),dc_all(1:Nx,1:7)
  !!!!!!!!!!!!!!!!!!
  REAL(8), PARAMETER :: alpha=SQRT(two)-one
  


  ! Treat velocities as zero when below v_min
  v_min=1.0d0

  ! Treat densities as zero when below v_min
  n_min=1.0d5

  dx=domain%dx
  dy=domain%dy


  nu_y_max=zero
  nu_x_max=zero

  inv_ne(:,:)= one/(domain%ions1%ni(:,:)+domain%ions2%ni(:,:)) 
  ni1_Bz_ne(:,:)=domain%ions1%ni(:,:)*domain%Bz(:,:)*inv_ne(:,:)
  ni2_Bz_ne(:,:)=domain%ions2%ni(:,:)*domain%Bz(:,:)*inv_ne(:,:)

  inv_e=one/domain%ions1%charge

  IF (direction.EQ.'x') THEN

    !!! x-direction
    CALL d1x_plus(domain%ions1%ni,dni1_d_plus,domain%dx)
    CALL d1x_plus(domain%ions1%vix,dvi1x_d_plus,domain%dx)
    CALL d1x_plus(domain%ions1%viy,dvi1y_d_plus,domain%dx)
    CALL d1x_plus(domain%ions2%ni,dni2_d_plus,domain%dx)
    CALL d1x_plus(domain%ions2%vix,dvi2x_d_plus,domain%dx)
    CALL d1x_plus(domain%ions2%viy,dvi2y_d_plus,domain%dx)
    CALL d1x_plus(domain%Bz,dBz_d_plus,domain%dx)

    sendbuf_all_r(:,1)=dni1_d_plus(Nx,:)
    sendbuf_all_r(:,2)=dvi1x_d_plus(Nx,:)
    sendbuf_all_r(:,3)=dvi1y_d_plus(Nx,:)
    sendbuf_all_r(:,4)=dni2_d_plus(Nx,:)
    sendbuf_all_r(:,5)=dvi2x_d_plus(Nx,:)
    sendbuf_all_r(:,6)=dvi2y_d_plus(Nx,:)
    sendbuf_all_r(:,7)=dBz_d_plus(Nx,:)

    sendtag=100+my_x
    recvtag=100+my_left
    CALL MPI_IRECV(lc_all,Ny*7,MPI_DOUBLE_PRECISION,  &
      my_left,recvtag,row_comm,request_recv,ierror)
    CALL MPI_ISEND(sendbuf_all_r,Ny*7,MPI_DOUBLE_PRECISION,  &
      my_right,sendtag,row_comm,request_send,ierror)

    dni1_d_minus(2:Nx,:)=dni1_d_plus(1:Nx-1,:)
    dvi1x_d_minus(2:Nx,:)=dvi1x_d_plus(1:Nx-1,:)
    dvi1y_d_minus(2:Nx,:)=dvi1y_d_plus(1:Nx-1,:)
    dni2_d_minus(2:Nx,:)=dni2_d_plus(1:Nx-1,:)
    dvi2x_d_minus(2:Nx,:)=dvi2x_d_plus(1:Nx-1,:)
    dvi2y_d_minus(2:Nx,:)=dvi2y_d_plus(1:Nx-1,:)
    dBz_d_minus(2:Nx,:)=dBz_d_plus(1:Nx-1,:)

    CALL MPI_WAIT(request_send,status,ierror)
    CALL MPI_WAIT(request_recv,status,ierror)

    dni1_d_minus(1,:)=lc_all(:,1)
    dvi1x_d_minus(1,:)=lc_all(:,2)
    dvi1y_d_minus(1,:)=lc_all(:,3)
    dni2_d_minus(1,:)=lc_all(:,4)
    dvi2x_d_minus(1,:)=lc_all(:,5)
    dvi2y_d_minus(1,:)=lc_all(:,6)
    dBz_d_minus(1,:)=lc_all(:,7)


    CALL d2y(domain%ions1%ni,dni1_d2,dy)
    CALL d2y(domain%ions1%vix,dvi1x_d2,dy)
    CALL d2y(domain%ions1%viy,dvi1y_d2,dy)
    CALL d2y(domain%ions2%ni,dni2_d2,dy)
    CALL d2y(domain%ions2%vix,dvi2x_d2,dy)
    CALL d2y(domain%ions2%viy,dvi2y_d2,dy)
    CALL d2y(domain%Bz,dBz_d2,dy)

    CALL d1y(inv_ne,dinv_ne_d,domain%dy)

    DO i2=1,Ny
    DO i1=1,Nx
      RHS_D%ions1%vix(i1,i2) = zero
      RHS_D%ions2%vix(i1,i2) = zero
      RHS_D%ions1%viy(i1,i2) = zero
      RHS_D%ions2%viy(i1,i2) = zero
      RHS_D%ions1%ni(i1,i2) = zero
      RHS_D%ions2%ni(i1,i2) = zero
      RHS_D%Bz(i1,i2)        = zero

      vi1x=domain%ions1%vix(i1,i2)-vx_frame_global
      vi2x=domain%ions2%vix(i1,i2)-vx_frame_global

      ! Small density ni2
      IF (ABS(domain%ions2%ni(i1,i2)).LT.n_min) THEN
        matrix_A_2x2(:,:)=zero

        matrix_A_2x2(1,1)=vi1x
        matrix_A_2x2(1,2)=inv_ne(i1,i2)*domain%Bz(i1,i2)/(mu0*domain%ions1%mass)

        matrix_A_2x2(2,1)=ni1_Bz_ne(i1,i2)
        matrix_A_2x2(2,2)=(domain%ions1%ni(i1,i2)*vi1x+  &
          domain%ions2%ni(i1,i2)*vi2x)*inv_ne(i1,i2)+  &
          inv_e*dinv_ne_d(i1,i2)*domain%Bz(i1,i2)/mu0

        d1_plus_2x2(1)=dvi1x_d_plus(i1,i2)
        d1_plus_2x2(2)=dBz_d_plus(i1,i2)

        d1_minus_2x2(1)=dvi1x_d_minus(i1,i2)
        d1_minus_2x2(2)=dBz_d_minus(i1,i2)

        d2_2x2(1)=dvi1x_d2(i1,i2)
        d2_2x2(2)=dBz_d2(i1,i2)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ldvl=1
        ldvr=lda_2x2
      
        CALL DGEEV('N','V',lda_2x2,matrix_A_2x2,lda_2x2,wr_2x2,wi_2x2,vl,ldvl,vr_2x2,ldvr,work,lwork,ierr)

        IF (MAXVAL(wi_2x2).EQ.zero) THEN
          d1_minus_plus_2x2(:,1)=d1_minus_2x2
          d1_minus_plus_2x2(:,2)=d1_plus_2x2
          d1_minus_plus_2x2(:,3)=d2_2x2

          vr_copy_2x2(:,:)=vr_2x2(:,:)
          CALL DGETRF(lda_2x2,lda_2x2,vr_2x2,lda_2x2,ipvt_2x2,ierr)

          CALL DGETRS('N', lda_2x2, 3, vr_2x2,lda_2x2,ipvt_2x2,d1_minus_plus_2x2,lda_2x2, ierr)

          DO i3=1,lda_2x2
            IF (wr_2x2(i3).GT.zero) THEN
              d1_2x2(i3)=d1_minus_plus_2x2(i3,1)*wr_2x2(i3)-alpha*half*dy*ABS(wr_2x2(i3))*d1_minus_plus_2x2(i3,3)
            ELSE
              d1_2x2(i3)=d1_minus_plus_2x2(i3,2)*wr_2x2(i3)-alpha*half*dy*ABS(wr_2x2(i3))*d1_minus_plus_2x2(i3,3)
            END IF
            nu_x_max=MAX(ABS(wr_2x2(i3)),nu_x_max)
          END DO

          CALL DGEMV('N',lda_2x2,lda_2x2,one,vr_copy_2x2,lda_2x2,d1_2x2,1,zero,d1_out_2x2,1)

        ELSE

          DO i3=1,lda_2x2
            IF (wi_2x2(i3).EQ.zero) THEN
              DO i4=1,lda_2x2
                vr_C_2x2(i4,i3)=DCMPLX(vr_2x2(i4,i3))            
              END DO
            ELSE IF (wi_2x2(i3).GT.zero) THEN
              DO i4=1,lda_2x2
                vr_C_2x2(i4,i3)=DCMPLX(vr_2x2(i4,i3),vr_2x2(i4,i3+1))            
              END DO
            ELSE
              DO i4=1,lda_2x2
                vr_C_2x2(i4,i3)=DCMPLX(vr_2x2(i4,i3-1),-vr_2x2(i4,i3))        
              END DO
            END IF
          END DO

          d1_minus_plus_C_2x2(:,1)=DCMPLX(d1_minus_2x2(:))
          d1_minus_plus_C_2x2(:,2)=DCMPLX(d1_plus_2x2(:))
          d1_minus_plus_C_2x2(:,3)=DCMPLX(d2_2x2(:))

          vr_C_copy_2x2(:,:)=vr_C_2x2(:,:)
          CALL ZGETRF(lda_2x2,lda_2x2,vr_C_2x2,lda_2x2,ipvt_2x2,ierr)
          CALL ZGETRS('N', lda_2x2, 3, vr_C_2x2, lda_2x2, ipvt_2x2, d1_minus_plus_C_2x2, lda_2x2, ierr)

          DO i3=1,lda_2x2
            IF (wr_2x2(i3).GT.zero) THEN
              d1_C_2x2(i3)=d1_minus_plus_C_2x2(i3,1)*DCMPLX(wr_2x2(i3),wi_2x2(i3))  &
                      -alpha*half*dy*ABS(DCMPLX(wr_2x2(i3),wi_2x2(i3)))*d1_minus_plus_C_2x2(i3,3)
            ELSE
              d1_C_2x2(i3)=d1_minus_plus_C_2x2(i3,2)*DCMPLX(wr_2x2(i3),wi_2x2(i3))  &
                      -alpha*half*dy*ABS(DCMPLX(wr_2x2(i3),wi_2x2(i3)))*d1_minus_plus_C_2x2(i3,3)
            END IF
            nu_x_max=MAX(ABS(wr_2x2(i3)),nu_x_max)
          END DO

          CALL ZGEMV('N',lda_2x2,lda_2x2,DCMPLX(one),vr_C_copy_2x2,lda_2x2,d1_C_2x2,1,DCMPLX(zero),d1_out_C_2x2,1)

          d1_out_2x2(:)=DREAL(d1_out_C_2x2(:))
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        RHS_D%ions1%vix(i1,i2) = RHS_D%ions1%vix(i1,i2)-d1_out_2x2(1)
        RHS_D%Bz(i1,i2)        = RHS_D%Bz(i1,i2)-d1_out_2x2(2)

        RHS_D%ions1%ni(i1,i2) = RHS_D%ions1%ni(i1,i2) -  &
          domain%ions1%ni(i1,i2)*(dvi1x_d_minus(i1,i2)+dvi1x_d_plus(i1,i2))*half

        RHS_D%ions2%ni(i1,i2) = RHS_D%ions2%ni(i1,i2) -  &
          domain%ions2%ni(i1,i2)*(dvi2x_d_minus(i1,i2)+dvi2x_d_plus(i1,i2))*half

        RHS_D%ions2%vix(i1,i2) = RHS_D%ions2%vix(i1,i2) - &
          (inv_ne(i1,i2)*domain%Bz(i1,i2)/(mu0*domain%ions2%mass))*(dBz_d_minus(i1,i2)+dBz_d_plus(i1,i2))*half

        IF (vi1x.GT.zero) THEN
          RHS_D%ions1%viy(i1,i2) = RHS_D%ions1%viy(i1,i2)-vi1x*dvi1y_d_minus(i1,i2)
          RHS_D%ions1%ni(i1,i2)  = RHS_D%ions1%ni(i1,i2)-vi1x*dni1_d_minus(i1,i2)
        ELSE
          RHS_D%ions1%viy(i1,i2) = RHS_D%ions1%viy(i1,i2)-vi1x*dvi1y_d_plus(i1,i2)
          RHS_D%ions1%ni(i1,i2)  = RHS_D%ions1%ni(i1,i2)-vi1x*dni1_d_plus(i1,i2)
        END IF

        IF (vi2x.GT.zero) THEN
          RHS_D%ions2%viy(i1,i2) = RHS_D%ions2%viy(i1,i2)-vi2x*dvi2y_d_minus(i1,i2)
          RHS_D%ions2%ni(i1,i2)  = RHS_D%ions2%ni(i1,i2)-vi2x*dni2_d_minus(i1,i2)
          RHS_D%ions2%vix(i1,i2) = RHS_D%ions2%vix(i1,i2)-vi2x*dvi2x_d_minus(i1,i2)
        ELSE
          RHS_D%ions2%viy(i1,i2) = RHS_D%ions2%viy(i1,i2)-vi2x*dvi2y_d_plus(i1,i2)
          RHS_D%ions2%ni(i1,i2)  = RHS_D%ions2%ni(i1,i2)-vi2x*dni2_d_plus(i1,i2)
          RHS_D%ions2%vix(i1,i2) = RHS_D%ions2%vix(i1,i2)-vi2x*dvi2x_d_plus(i1,i2)
        END IF

        RHS_D%Bz(i1,i2)=RHS_D%Bz(i1,i2)-ni1_Bz_ne(i1,i2)*inv_ne(i1,i2)*  &
          (vi2x-vi1x)*  &
          (dni2_d_minus(i1,i2)+dni2_d_plus(i1,i2))*half
      !!! Small differential velocity vi1y-vi2y
      ELSE IF (ABS(vi1x-vi2x).LT.v_min) THEN
        matrix_A_3x3(:,:)=zero

        matrix_A_3x3(1,1)=vi1x
        matrix_A_3x3(1,3)=inv_ne(i1,i2)*domain%Bz(i1,i2)/(mu0*domain%ions1%mass)

        matrix_A_3x3(2,2)=vi2x
        matrix_A_3x3(2,3)=inv_ne(i1,i2)*domain%Bz(i1,i2)/(mu0*domain%ions2%mass)

        matrix_A_3x3(3,1)=ni1_Bz_ne(i1,i2)
        matrix_A_3x3(3,2)=ni2_Bz_ne(i1,i2)
        matrix_A_3x3(3,3)=(domain%ions1%ni(i1,i2)*vi1x+ &
          domain%ions2%ni(i1,i2)*vi2x)*inv_ne(i1,i2)+ &
          inv_e*dinv_ne_d(i1,i2)*domain%Bz(i1,i2)/mu0

        d1_plus_3x3(1)=dvi1x_d_plus(i1,i2)
        d1_plus_3x3(2)=dvi2x_d_plus(i1,i2)
        d1_plus_3x3(3)=dBz_d_plus(i1,i2)

        d1_minus_3x3(1)=dvi1x_d_minus(i1,i2)
        d1_minus_3x3(2)=dvi2x_d_minus(i1,i2)
        d1_minus_3x3(3)=dBz_d_minus(i1,i2)

        d2_3x3(1)=dvi1x_d2(i1,i2)
        d2_3x3(2)=dvi2x_d2(i1,i2)
        d2_3x3(3)=dBz_d2(i1,i2)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ldvl=1
        ldvr=lda_3x3
      
        CALL DGEEV('N','V',lda_3x3,matrix_A_3x3,lda_3x3,wr_3x3,wi_3x3,vl,ldvl,vr_3x3,ldvr,work,lwork,ierr)

        IF (MAXVAL(wi_3x3).EQ.zero) THEN
          d1_minus_plus_3x3(:,1)=d1_minus_3x3
          d1_minus_plus_3x3(:,2)=d1_plus_3x3
          d1_minus_plus_3x3(:,3)=d2_3x3          

          vr_copy_3x3(:,:)=vr_3x3(:,:)
          CALL DGETRF(lda_3x3,lda_3x3,vr_3x3,lda_3x3,ipvt_3x3,ierr)

          CALL DGETRS('N', lda_3x3, 3, vr_3x3,lda_3x3,ipvt_3x3,d1_minus_plus_3x3,lda_3x3, ierr)

          DO i3=1,lda_3x3
            IF (wr_3x3(i3).GT.zero) THEN
              d1_3x3(i3)=d1_minus_plus_3x3(i3,1)*wr_3x3(i3)-alpha*half*dy*ABS(wr_3x3(i3))*d1_minus_plus_3x3(i3,3)
            ELSE
              d1_3x3(i3)=d1_minus_plus_3x3(i3,2)*wr_3x3(i3)-alpha*half*dy*ABS(wr_3x3(i3))*d1_minus_plus_3x3(i3,3)
            END IF
            nu_x_max=MAX(ABS(wr_3x3(i3)),nu_x_max)
          END DO

          CALL DGEMV('N',lda_3x3,lda_3x3,one,vr_copy_3x3,lda_3x3,d1_3x3,1,zero,d1_out_3x3,1)

        ELSE

          DO i3=1,lda_3x3
            IF (wi_3x3(i3).EQ.zero) THEN
              DO i4=1,lda_3x3
                vr_C_3x3(i4,i3)=DCMPLX(vr_3x3(i4,i3))            
              END DO
            ELSE IF (wi_3x3(i3).GT.zero) THEN
              DO i4=1,lda_3x3
                vr_C_3x3(i4,i3)=DCMPLX(vr_3x3(i4,i3),vr_3x3(i4,i3+1))            
              END DO
            ELSE
              DO i4=1,lda_3x3
                vr_C_3x3(i4,i3)=DCMPLX(vr_3x3(i4,i3-1),-vr_3x3(i4,i3))        
              END DO
            END IF
          END DO

          d1_minus_plus_C_3x3(:,1)=DCMPLX(d1_minus_3x3(:))
          d1_minus_plus_C_3x3(:,2)=DCMPLX(d1_plus_3x3(:))
          d1_minus_plus_C_3x3(:,3)=DCMPLX(d2_3x3(:))

          vr_C_copy_3x3(:,:)=vr_C_3x3(:,:)
          CALL ZGETRF(lda_3x3,lda_3x3,vr_C_3x3,lda_3x3,ipvt_3x3,ierr)
          CALL ZGETRS('N', lda_3x3, 3, vr_C_3x3, lda_3x3, ipvt_3x3, d1_minus_plus_C_3x3, lda_3x3, ierr)

          DO i3=1,lda_3x3
            IF (wr_3x3(i3).GT.zero) THEN
              d1_C_3x3(i3)=d1_minus_plus_C_3x3(i3,1)*DCMPLX(wr_3x3(i3),wi_3x3(i3))  &
                          -alpha*half*dy*ABS(DCMPLX(wr_3x3(i3),wi_3x3(i3)))*d1_minus_plus_C_3x3(i3,3)
              ELSE
              d1_C_3x3(i3)=d1_minus_plus_C_3x3(i3,2)*DCMPLX(wr_3x3(i3),wi_3x3(i3))  &
                          -alpha*half*dy*ABS(DCMPLX(wr_3x3(i3),wi_3x3(i3)))*d1_minus_plus_C_3x3(i3,3)
            END IF
            nu_x_max=MAX(ABS(wr_3x3(i3)),nu_x_max)
          END DO

          CALL ZGEMV('N',lda_3x3,lda_3x3,DCMPLX(one),vr_C_copy_3x3,lda_3x3,d1_C_3x3,1,DCMPLX(zero),d1_out_C_3x3,1)

          d1_out_3x3(:)=DREAL(d1_out_C_3x3(:))
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        RHS_D%ions1%vix(i1,i2) = RHS_D%ions1%vix(i1,i2)-d1_out_3x3(1)
        RHS_D%ions2%vix(i1,i2) = RHS_D%ions2%vix(i1,i2)-d1_out_3x3(2)
        RHS_D%Bz(i1,i2)        = RHS_D%Bz(i1,i2)-d1_out_3x3(3)

        RHS_D%ions1%ni(i1,i2) = RHS_D%ions1%ni(i1,i2)-  &
          domain%ions1%ni(i1,i2)*(dvi1x_d_minus(i1,i2)+dvi1x_d_plus(i1,i2))*half

        RHS_D%ions2%ni(i1,i2) = RHS_D%ions2%ni(i1,i2)-  &
          domain%ions2%ni(i1,i2)*(dvi2x_d_minus(i1,i2)+dvi2x_d_plus(i1,i2))*half

        IF (vi1x.GT.zero) THEN
          RHS_D%ions1%viy(i1,i2) = RHS_D%ions1%viy(i1,i2)-vi1x*dvi1y_d_minus(i1,i2)
          RHS_D%ions1%ni(i1,i2) = RHS_D%ions1%ni(i1,i2)-vi1x*dni1_d_minus(i1,i2)
        ELSE
          RHS_D%ions1%viy(i1,i2) = RHS_D%ions1%viy(i1,i2)-vi1x*dvi1y_d_plus(i1,i2)
          RHS_D%ions1%ni(i1,i2) = RHS_D%ions1%ni(i1,i2)-vi1x*dni1_d_plus(i1,i2)
        END IF

        IF (vi2x.GT.zero) THEN
          RHS_D%ions2%viy(i1,i2) = RHS_D%ions2%viy(i1,i2)-vi2x*dvi2y_d_minus(i1,i2)
          RHS_D%ions2%ni(i1,i2) = RHS_D%ions2%ni(i1,i2)-vi2x*dni2_d_minus(i1,i2)
        ELSE
          RHS_D%ions2%viy(i1,i2) = RHS_D%ions2%viy(i1,i2)-vi2x*dvi2y_d_plus(i1,i2)
          RHS_D%ions2%ni(i1,i2) = RHS_D%ions2%ni(i1,i2)-vi2x*dni2_d_plus(i1,i2)
        END IF

      ELSE
        matrix_A(:,:)=zero

        matrix_A(1,1)=vi1x
        matrix_A(1,2)=domain%ions1%ni(i1,i2) 

        matrix_A(2,2)=vi1x 
        matrix_A(2,5)=inv_ne(i1,i2)*domain%Bz(i1,i2)/(mu0*domain%ions1%mass)

        matrix_A(3,3)=vi2x
        matrix_A(3,4)=domain%ions2%ni(i1,i2)
      
        matrix_A(4,4)=vi2x
        matrix_A(4,5)=inv_ne(i1,i2)*domain%Bz(i1,i2)/(mu0*domain%ions2%mass)
      
        matrix_A(5,1)=ni2_Bz_ne(i1,i2)*inv_ne(i1,i2)*(vi1x-vi2x)
        matrix_A(5,2)=ni1_Bz_ne(i1,i2)
        matrix_A(5,3)=ni1_Bz_ne(i1,i2)*inv_ne(i1,i2)*(vi2x-vi1x)
        matrix_A(5,4)=ni2_Bz_ne(i1,i2)
        matrix_A(5,5)=(domain%ions1%ni(i1,i2)*vi1x+  &
          domain%ions2%ni(i1,i2)*vi2x)*inv_ne(i1,i2)+  &
          inv_e*dinv_ne_d(i1,i2)*domain%Bz(i1,i2)/mu0

        d1_plus(1)=dni1_d_plus(i1,i2)
        d1_plus(2)=dvi1x_d_plus(i1,i2)
        d1_plus(3)=dni2_d_plus(i1,i2)
        d1_plus(4)=dvi2x_d_plus(i1,i2)
        d1_plus(5)=dBz_d_plus(i1,i2)

        d1_minus(1)=dni1_d_minus(i1,i2)
        d1_minus(2)=dvi1x_d_minus(i1,i2)
        d1_minus(3)=dni2_d_minus(i1,i2)
        d1_minus(4)=dvi2x_d_minus(i1,i2)
        d1_minus(5)=dBz_d_minus(i1,i2)

        d2(1)=dni1_d2(i1,i2)
        d2(2)=dvi1x_d2(i1,i2)
        d2(3)=dni2_d2(i1,i2)
        d2(4)=dvi2x_d2(i1,i2)
        d2(5)=dBz_d2(i1,i2)

        ldvl=1
        ldvr=lda
        CALL DGEEV('N','V',lda,matrix_A,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,ierr)

        IF (MAXVAL(wi).EQ.zero) THEN
          d1_minus_plus(:,1)=d1_minus
          d1_minus_plus(:,2)=d1_plus
          d1_minus_plus(:,3)=d2

          vr_copy(:,:)=vr(:,:)
          CALL DGETRF(lda,lda,vr,lda,ipvt,ierr)

          CALL DGETRS('N', lda, 3, vr, lda, ipvt, d1_minus_plus, lda, ierr)

          DO i3=1,lda
            IF (wr(i3).GT.zero) THEN
              d1(i3)=d1_minus_plus(i3,1)*wr(i3)-alpha*half*dy*ABS(wr(i3))*d1_minus_plus(i3,3)
            ELSE
              d1(i3)=d1_minus_plus(i3,2)*wr(i3)-alpha*half*dy*ABS(wr(i3))*d1_minus_plus(i3,3)
            END IF
            nu_x_max=MAX(ABS(wr(i3)),nu_x_max)
          END DO

          CALL DGEMV('N',lda,lda,one,vr_copy,lda,d1,1,zero,d1_out,1)

        ELSE

          DO i3=1,lda
            IF (wi(i3).EQ.zero) THEN
              DO i4=1,lda
                vr_C(i4,i3)=DCMPLX(vr(i4,i3))            
              END DO
            ELSE IF (wi(i3).GT.zero) THEN
              DO i4=1,lda
                vr_C(i4,i3)=DCMPLX(vr(i4,i3),vr(i4,i3+1))            
              END DO
            ELSE
              DO i4=1,lda
                vr_C(i4,i3)=DCMPLX(vr(i4,i3-1),-vr(i4,i3))        
              END DO
            END IF
          END DO

          d1_minus_plus_C(:,1)=DCMPLX(d1_minus(:))
          d1_minus_plus_C(:,2)=DCMPLX(d1_plus(:))
          d1_minus_plus_C(:,3)=DCMPLX(d2(:))
   
          vr_C_copy(:,:)=vr_C(:,:)
          CALL ZGETRF(lda,lda,vr_C,lda,ipvt,ierr)
          CALL ZGETRS('N', lda, 3, vr_C, lda, ipvt, d1_minus_plus_C, lda, ierr)

          DO i3=1,lda
            IF (wr(i3).GT.zero) THEN
              d1_C(i3)=d1_minus_plus_C(i3,1)*DCMPLX(wr(i3),wi(i3))  &
                       -alpha*half*dy*ABS(DCMPLX(wr(i3),wi(i3)))*d1_minus_plus_C(i3,3)
            ELSE
              d1_C(i3)=d1_minus_plus_C(i3,2)*DCMPLX(wr(i3),wi(i3))  &
                       -alpha*half*dy*ABS(DCMPLX(wr(i3),wi(i3)))*d1_minus_plus_C(i3,3)
            END IF
            nu_x_max=MAX(ABS(wr(i3)),nu_x_max)
          END DO

          CALL ZGEMV('N',lda,lda,DCMPLX(one),vr_C_copy,lda,d1_C,1,DCMPLX(zero),d1_out_C,1)

          d1_out(:)=DREAL(d1_out_C(:))
        END IF

        RHS_D%ions1%ni(i1,i2)  = RHS_D%ions1%ni(i1,i2)-d1_out(1)
        RHS_D%ions1%vix(i1,i2) = RHS_D%ions1%vix(i1,i2)-d1_out(2)
        RHS_D%ions2%ni(i1,i2)  = RHS_D%ions2%ni(i1,i2)-d1_out(3)
        RHS_D%ions2%vix(i1,i2) = RHS_D%ions2%vix(i1,i2)-d1_out(4)
        RHS_D%Bz(i1,i2)        = RHS_D%Bz(i1,i2)-d1_out(5)

        IF (vi1x.GT.zero) THEN
          RHS_D%ions1%viy(i1,i2) = RHS_D%ions1%viy(i1,i2)-vi1x*dvi1y_d_minus(i1,i2)
        ELSE
          RHS_D%ions1%viy(i1,i2) = RHS_D%ions1%viy(i1,i2)-vi1x*dvi1y_d_plus(i1,i2)
        END IF

        IF (vi2x.GT.zero) THEN
          RHS_D%ions2%viy(i1,i2) = RHS_D%ions2%viy(i1,i2)-vi2x*dvi2y_d_minus(i1,i2)
        ELSE
          RHS_D%ions2%viy(i1,i2) = RHS_D%ions2%viy(i1,i2)-vi2x*dvi2y_d_plus(i1,i2)
        END IF

      END IF

    END DO
    END DO

    RHS_D%nu_x_max=nu_x_max

  ELSE IF(direction.EQ.'y') THEN
    !!! y-direction

    CALL d1y_plus(domain%ions1%ni,dni1_d_plus,domain%dy)
    CALL d1y_plus(domain%ions1%vix,dvi1x_d_plus,domain%dy)
    CALL d1y_plus(domain%ions1%viy,dvi1y_d_plus,domain%dy)
    CALL d1y_plus(domain%ions2%ni,dni2_d_plus,domain%dy)
    CALL d1y_plus(domain%ions2%vix,dvi2x_d_plus,domain%dy)
    CALL d1y_plus(domain%ions2%viy,dvi2y_d_plus,domain%dy)
    CALL d1y_plus(domain%Bz,dBz_d_plus,domain%dx)


    sendbuf_all_u(:,1)=dni1_d_plus(:,Ny)
    sendbuf_all_u(:,2)=dvi1x_d_plus(:,Ny)
    sendbuf_all_u(:,3)=dvi1y_d_plus(:,Ny)
    sendbuf_all_u(:,4)=dni2_d_plus(:,Ny)
    sendbuf_all_u(:,5)=dvi2x_d_plus(:,Ny)
    sendbuf_all_u(:,6)=dvi2y_d_plus(:,Ny)
    sendbuf_all_u(:,7)=dBz_d_plus(:,Ny)

    sendtag=100+my_y
    recvtag=100+my_down
    CALL MPI_IRECV(dc_all,Nx*7,MPI_DOUBLE_PRECISION,  &
      my_down,recvtag,col_comm,request_recv,ierror)
    CALL MPI_ISEND(sendbuf_all_u,Nx*7,MPI_DOUBLE_PRECISION,  &
      my_up,sendtag,col_comm,request_send,ierror)

    dni1_d_minus(:,2:Ny)=dni1_d_plus(:,1:Ny-1)
    dvi1x_d_minus(:,2:Ny)=dvi1x_d_plus(:,1:Ny-1)
    dvi1y_d_minus(:,2:Ny)=dvi1y_d_plus(:,1:Ny-1)
    dni2_d_minus(:,2:Ny)=dni2_d_plus(:,1:Ny-1)
    dvi2x_d_minus(:,2:Ny)=dvi2x_d_plus(:,1:Ny-1)
    dvi2y_d_minus(:,2:Ny)=dvi2y_d_plus(:,1:Ny-1)
    dBz_d_minus(:,2:Ny)=dBz_d_plus(:,1:Ny-1)

    CALL MPI_WAIT(request_send,status,ierror)
    CALL MPI_WAIT(request_recv,status,ierror)

    dni1_d_minus(:,1)=dc_all(:,1)
    dvi1x_d_minus(:,1)=dc_all(:,2)
    dvi1y_d_minus(:,1)=dc_all(:,3)
    dni2_d_minus(:,1)=dc_all(:,4)
    dvi2x_d_minus(:,1)=dc_all(:,5)
    dvi2y_d_minus(:,1)=dc_all(:,6)
    dBz_d_minus(:,1)=dc_all(:,7)

    CALL d2x(domain%ions1%ni,dni1_d2,dx)
    CALL d2x(domain%ions1%vix,dvi1x_d2,dx)
    CALL d2x(domain%ions1%viy,dvi1y_d2,dx)
    CALL d2x(domain%ions2%ni,dni2_d2,dx)
    CALL d2x(domain%ions2%vix,dvi2x_d2,dx)
    CALL d2x(domain%ions2%viy,dvi2y_d2,dx)
    CALL d2x(domain%Bz,dBz_d2,dx)

    CALL d1x(inv_ne,dinv_ne_d,domain%dx)

    DO i2=1,Ny
    DO i1=1,Nx
      RHS_D%ions1%vix(i1,i2) = zero
      RHS_D%ions2%vix(i1,i2) = zero
      RHS_D%ions1%viy(i1,i2) = zero
      RHS_D%ions2%viy(i1,i2) = zero
      RHS_D%ions1%ni(i1,i2) = zero
      RHS_D%ions2%ni(i1,i2) = zero
      RHS_D%Bz(i1,i2)        = zero

      ! Small density ni2
      IF (ABS(domain%ions2%ni(i1,i2)).LT.n_min) THEN
        matrix_A_2x2(:,:)=zero

        matrix_A_2x2(1,1)=domain%ions1%viy(i1,i2)
        matrix_A_2x2(1,2)=inv_ne(i1,i2)*domain%Bz(i1,i2)/(mu0*domain%ions1%mass)

        matrix_A_2x2(2,1)=ni1_Bz_ne(i1,i2)
        matrix_A_2x2(2,2)=(domain%ions1%ni(i1,i2)*domain%ions1%viy(i1,i2)+ &
          domain%ions2%ni(i1,i2)*domain%ions2%viy(i1,i2))*inv_ne(i1,i2)- &
          inv_e*dinv_ne_d(i1,i2)*domain%Bz(i1,i2)/mu0

        d1_plus_2x2(1)=dvi1y_d_plus(i1,i2)
        d1_plus_2x2(2)=dBz_d_plus(i1,i2)

        d1_minus_2x2(1)=dvi1y_d_minus(i1,i2)
        d1_minus_2x2(2)=dBz_d_minus(i1,i2)

        d2_2x2(1)=dvi1y_d2(i1,i2)
        d2_2x2(2)=dBz_d2(i1,i2)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ldvl=1
        ldvr=lda_2x2
      
        CALL DGEEV('N','V',lda_2x2,matrix_A_2x2,lda_2x2,wr_2x2,wi_2x2,vl,ldvl,vr_2x2,ldvr,work,lwork,ierr)

        IF (MAXVAL(wi_2x2).EQ.zero) THEN
          d1_minus_plus_2x2(:,1)=d1_minus_2x2
          d1_minus_plus_2x2(:,2)=d1_plus_2x2
          d1_minus_plus_2x2(:,3)=d2_2x2

          vr_copy_2x2(:,:)=vr_2x2(:,:)
          CALL DGETRF(lda_2x2,lda_2x2,vr_2x2,lda_2x2,ipvt_2x2,ierr)

          CALL DGETRS('N', lda_2x2, 3, vr_2x2,lda_2x2,ipvt_2x2,d1_minus_plus_2x2,lda_2x2, ierr)

          DO i3=1,lda_2x2
            IF (wr_2x2(i3).GT.zero) THEN
              d1_2x2(i3)=d1_minus_plus_2x2(i3,1)*wr_2x2(i3)-alpha*half*dx*ABS(wr_2x2(i3))*d1_minus_plus_2x2(i3,3)
            ELSE
              d1_2x2(i3)=d1_minus_plus_2x2(i3,2)*wr_2x2(i3)-alpha*half*dx*ABS(wr_2x2(i3))*d1_minus_plus_2x2(i3,3)
            END IF
            nu_y_max=MAX(ABS(wr_2x2(i3)),nu_y_max)
          END DO

          CALL DGEMV('N',lda_2x2,lda_2x2,one,vr_copy_2x2,lda_2x2,d1_2x2,1,zero,d1_out_2x2,1)

        ELSE

          DO i3=1,lda_2x2
            IF (wi_2x2(i3).EQ.zero) THEN
              DO i4=1,lda_2x2
                vr_C_2x2(i4,i3)=DCMPLX(vr_2x2(i4,i3))            
              END DO
            ELSE IF (wi_2x2(i3).GT.zero) THEN
              DO i4=1,lda_2x2
                vr_C_2x2(i4,i3)=DCMPLX(vr_2x2(i4,i3),vr_2x2(i4,i3+1))            
              END DO
            ELSE
              DO i4=1,lda_2x2
                vr_C_2x2(i4,i3)=DCMPLX(vr_2x2(i4,i3-1),-vr_2x2(i4,i3))        
              END DO
            END IF
          END DO

          d1_minus_plus_C_2x2(:,1)=DCMPLX(d1_minus_2x2(:))
          d1_minus_plus_C_2x2(:,2)=DCMPLX(d1_plus_2x2(:))
          d1_minus_plus_C_2x2(:,3)=DCMPLX(d2_2x2(:))

          vr_C_copy_2x2(:,:)=vr_C_2x2(:,:)
          CALL ZGETRF(lda_2x2,lda_2x2,vr_C_2x2,lda_2x2,ipvt_2x2,ierr)
          CALL ZGETRS('N', lda_2x2, 3, vr_C_2x2, lda_2x2, ipvt_2x2, d1_minus_plus_C_2x2, lda_2x2, ierr)

          DO i3=1,lda_2x2
            IF (wr_2x2(i3).GT.zero) THEN
              d1_C_2x2(i3)=d1_minus_plus_C_2x2(i3,1)*DCMPLX(wr_2x2(i3),wi_2x2(i3))  &
                      -alpha*half*dx*ABS(DCMPLX(wr_2x2(i3),wi_2x2(i3)))*d1_minus_plus_C_2x2(i3,3)
            ELSE
              d1_C_2x2(i3)=d1_minus_plus_C_2x2(i3,2)*DCMPLX(wr_2x2(i3),wi_2x2(i3))  &
                      -alpha*half*dx*ABS(DCMPLX(wr_2x2(i3),wi_2x2(i3)))*d1_minus_plus_C_2x2(i3,3)
            END IF
            nu_y_max=MAX(ABS(wr_2x2(i3)),nu_y_max)
          END DO

          CALL ZGEMV('N',lda_2x2,lda_2x2,DCMPLX(one),vr_C_copy_2x2,lda_2x2,d1_C_2x2,1,DCMPLX(zero),d1_out_C_2x2,1)

          d1_out_2x2(:)=DREAL(d1_out_C_2x2(:))
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        RHS_D%ions1%viy(i1,i2) = RHS_D%ions1%viy(i1,i2)-d1_out_2x2(1)
        RHS_D%Bz(i1,i2)        = RHS_D%Bz(i1,i2)-d1_out_2x2(2)

        RHS_D%ions1%ni(i1,i2) = RHS_D%ions1%ni(i1,i2)-  &
          domain%ions1%ni(i1,i2)*(dvi1y_d_minus(i1,i2)+dvi1y_d_plus(i1,i2))*half

        RHS_D%ions2%ni(i1,i2) = RHS_D%ions2%ni(i1,i2)-  &
          domain%ions2%ni(i1,i2)*(dvi2y_d_minus(i1,i2)+dvi2y_d_plus(i1,i2))*half

        RHS_D%ions2%viy(i1,i2) = RHS_D%ions2%viy(i1,i2)-  &
          (inv_ne(i1,i2)*domain%Bz(i1,i2)/(mu0*domain%ions2%mass))*(dBz_d_minus(i1,i2)+dBz_d_plus(i1,i2))*half

        IF (domain%ions1%viy(i1,i2).GT.zero) THEN
          RHS_D%ions1%vix(i1,i2) = RHS_D%ions1%vix(i1,i2)-domain%ions1%viy(i1,i2)*dvi1x_d_minus(i1,i2)
          RHS_D%ions1%ni(i1,i2)  = RHS_D%ions1%ni(i1,i2)-domain%ions1%viy(i1,i2)*dni1_d_minus(i1,i2)
        ELSE
          RHS_D%ions1%vix(i1,i2) = RHS_D%ions1%vix(i1,i2)-domain%ions1%viy(i1,i2)*dvi1x_d_plus(i1,i2)
          RHS_D%ions1%ni(i1,i2)  = RHS_D%ions1%ni(i1,i2)-domain%ions1%viy(i1,i2)*dni1_d_plus(i1,i2)
        END IF

        IF (domain%ions2%viy(i1,i2).GT.zero) THEN
          RHS_D%ions2%vix(i1,i2) = RHS_D%ions2%vix(i1,i2)-domain%ions2%viy(i1,i2)*dvi2x_d_minus(i1,i2)
          RHS_D%ions2%ni(i1,i2)  = RHS_D%ions2%ni(i1,i2)-domain%ions2%viy(i1,i2)*dni2_d_minus(i1,i2)
          RHS_D%ions2%viy(i1,i2) = RHS_D%ions2%viy(i1,i2)-domain%ions2%viy(i1,i2)*dvi2y_d_minus(i1,i2)
        ELSE
          RHS_D%ions2%vix(i1,i2) = RHS_D%ions2%vix(i1,i2)-domain%ions2%viy(i1,i2)*dvi2x_d_plus(i1,i2)
          RHS_D%ions2%ni(i1,i2)  = RHS_D%ions2%ni(i1,i2)-domain%ions2%viy(i1,i2)*dni2_d_plus(i1,i2)
          RHS_D%ions2%viy(i1,i2) = RHS_D%ions2%viy(i1,i2)-domain%ions2%viy(i1,i2)*dvi2y_d_plus(i1,i2)
        END IF
        
        RHS_D%Bz(i1,i2)=RHS_D%Bz(i1,i2)-ni1_Bz_ne(i1,i2)*inv_ne(i1,i2)*  &
          (domain%ions2%viy(i1,i2)-domain%ions1%viy(i1,i2))*  &
          (dni2_d_minus(i1,i2)+dni2_d_plus(i1,i2))*half

      !!! Small differential velocity vi1y-vi2y
      ELSE IF (ABS(domain%ions1%viy(i1,i2)-domain%ions2%viy(i1,i2)).LT.v_min) THEN
        matrix_A_3x3(:,:)=zero

        matrix_A_3x3(1,1)=domain%ions1%viy(i1,i2)
        matrix_A_3x3(1,3)=inv_ne(i1,i2)*domain%Bz(i1,i2)/(mu0*domain%ions1%mass)

        matrix_A_3x3(2,2)=domain%ions2%viy(i1,i2)
        matrix_A_3x3(2,3)=inv_ne(i1,i2)*domain%Bz(i1,i2)/(mu0*domain%ions2%mass)

        matrix_A_3x3(3,1)=ni1_Bz_ne(i1,i2)
        matrix_A_3x3(3,2)=ni2_Bz_ne(i1,i2)
        matrix_A_3x3(3,3)=(domain%ions1%ni(i1,i2)*domain%ions1%viy(i1,i2)+ &
          domain%ions2%ni(i1,i2)*domain%ions2%viy(i1,i2))*inv_ne(i1,i2)-  &
          inv_e*dinv_ne_d(i1,i2)*domain%Bz(i1,i2)/mu0

        d1_plus_3x3(1)=dvi1y_d_plus(i1,i2)
        d1_plus_3x3(2)=dvi2y_d_plus(i1,i2)
        d1_plus_3x3(3)=dBz_d_plus(i1,i2)

        d1_minus_3x3(1)=dvi1y_d_minus(i1,i2)
        d1_minus_3x3(2)=dvi2y_d_minus(i1,i2)
        d1_minus_3x3(3)=dBz_d_minus(i1,i2)

        d2_3x3(1)=dvi1y_d2(i1,i2)
        d2_3x3(2)=dvi2y_d2(i1,i2)
        d2_3x3(3)=dBz_d2(i1,i2)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ldvl=1
        ldvr=lda_3x3
      
        CALL DGEEV('N','V',lda_3x3,matrix_A_3x3,lda_3x3,wr_3x3,wi_3x3,vl,ldvl,vr_3x3,ldvr,work,lwork,ierr)

        IF (MAXVAL(wi_3x3).EQ.zero) THEN
          d1_minus_plus_3x3(:,1)=d1_minus_3x3
          d1_minus_plus_3x3(:,2)=d1_plus_3x3
          d1_minus_plus_3x3(:,3)=d2_3x3

          vr_copy_3x3(:,:)=vr_3x3(:,:)
          CALL DGETRF(lda_3x3,lda_3x3,vr_3x3,lda_3x3,ipvt_3x3,ierr)

          CALL DGETRS('N', lda_3x3, 3, vr_3x3,lda_3x3,ipvt_3x3,d1_minus_plus_3x3,lda_3x3, ierr)

          DO i3=1,lda_3x3
            IF (wr_3x3(i3).GT.zero) THEN
              d1_3x3(i3)=d1_minus_plus_3x3(i3,1)*wr_3x3(i3)-alpha*half*dx*ABS(wr_3x3(i3))*d1_minus_plus_3x3(i3,3)
            ELSE
              d1_3x3(i3)=d1_minus_plus_3x3(i3,2)*wr_3x3(i3)-alpha*half*dx*ABS(wr_3x3(i3))*d1_minus_plus_3x3(i3,3)
            END IF
            nu_y_max=MAX(ABS(wr_3x3(i3)),nu_y_max)
          END DO

          CALL DGEMV('N',lda_3x3,lda_3x3,one,vr_copy_3x3,lda_3x3,d1_3x3,1,zero,d1_out_3x3,1)

        ELSE

          DO i3=1,lda_3x3
            IF (wi_3x3(i3).EQ.zero) THEN
              DO i4=1,lda_3x3
                vr_C_3x3(i4,i3)=DCMPLX(vr_3x3(i4,i3))            
              END DO
            ELSE IF (wi_3x3(i3).GT.zero) THEN
              DO i4=1,lda_3x3
                vr_C_3x3(i4,i3)=DCMPLX(vr_3x3(i4,i3),vr_3x3(i4,i3+1))            
              END DO
            ELSE
              DO i4=1,lda_3x3
                vr_C_3x3(i4,i3)=DCMPLX(vr_3x3(i4,i3-1),-vr_3x3(i4,i3))        
              END DO
            END IF
          END DO

          d1_minus_plus_C_3x3(:,1)=DCMPLX(d1_minus_3x3(:))
          d1_minus_plus_C_3x3(:,2)=DCMPLX(d1_plus_3x3(:))
          d1_minus_plus_C_3x3(:,3)=DCMPLX(d2_3x3(:))

          vr_C_copy_3x3(:,:)=vr_C_3x3(:,:)
          CALL ZGETRF(lda_3x3,lda_3x3,vr_C_3x3,lda_3x3,ipvt_3x3,ierr)
          CALL ZGETRS('N', lda_3x3, 3, vr_C_3x3, lda_3x3, ipvt_3x3, d1_minus_plus_C_3x3, lda_3x3, ierr)

          DO i3=1,lda_3x3
            IF (wr_3x3(i3).GT.zero) THEN
              d1_C_3x3(i3)=d1_minus_plus_C_3x3(i3,1)*DCMPLX(wr_3x3(i3),wi_3x3(i3))  &
                      -alpha*half*dx*ABS(DCMPLX(wr_3x3(i3),wi_3x3(i3)))*d1_minus_plus_C_3x3(i3,3)
            ELSE
              d1_C_3x3(i3)=d1_minus_plus_C_3x3(i3,2)*DCMPLX(wr_3x3(i3),wi_3x3(i3))  &
                      -alpha*half*dx*ABS(DCMPLX(wr_3x3(i3),wi_3x3(i3)))*d1_minus_plus_C_3x3(i3,3)
            END IF
            nu_y_max=MAX(ABS(wr_3x3(i3)),nu_y_max)
          END DO

          CALL ZGEMV('N',lda_3x3,lda_3x3,DCMPLX(one),vr_C_copy_3x3,lda_3x3,d1_C_3x3,1,DCMPLX(zero),d1_out_C_3x3,1)

          d1_out_3x3(:)=DREAL(d1_out_C_3x3(:))
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        RHS_D%ions1%viy(i1,i2) = RHS_D%ions1%viy(i1,i2)-d1_out_3x3(1)
        RHS_D%ions2%viy(i1,i2) = RHS_D%ions2%viy(i1,i2)-d1_out_3x3(2)
        RHS_D%Bz(i1,i2)        = RHS_D%Bz(i1,i2)-d1_out_3x3(3)

        RHS_D%ions1%ni(i1,i2) = RHS_D%ions1%ni(i1,i2)-  &
          domain%ions1%ni(i1,i2)*(dvi1y_d_minus(i1,i2)+dvi1y_d_plus(i1,i2))*half

        RHS_D%ions2%ni(i1,i2) = RHS_D%ions2%ni(i1,i2)-  &
          domain%ions2%ni(i1,i2)*(dvi2y_d_minus(i1,i2)+dvi2y_d_plus(i1,i2))*half

        IF (domain%ions1%viy(i1,i2).GT.zero) THEN
          RHS_D%ions1%vix(i1,i2) = RHS_D%ions1%vix(i1,i2)-domain%ions1%viy(i1,i2)*dvi1x_d_minus(i1,i2)
          RHS_D%ions1%ni(i1,i2) = RHS_D%ions1%ni(i1,i2)-domain%ions1%viy(i1,i2)*dni1_d_minus(i1,i2)
        ELSE
          RHS_D%ions1%vix(i1,i2) = RHS_D%ions1%vix(i1,i2)-domain%ions1%viy(i1,i2)*dvi1x_d_plus(i1,i2)
          RHS_D%ions1%ni(i1,i2) = RHS_D%ions1%ni(i1,i2)-domain%ions1%viy(i1,i2)*dni1_d_plus(i1,i2)
        END IF

        IF (domain%ions2%viy(i1,i2).GT.zero) THEN
          RHS_D%ions2%vix(i1,i2) = RHS_D%ions2%vix(i1,i2)-domain%ions2%viy(i1,i2)*dvi2x_d_minus(i1,i2)
          RHS_D%ions2%ni(i1,i2) = RHS_D%ions2%ni(i1,i2)-domain%ions2%viy(i1,i2)*dni2_d_minus(i1,i2)
        ELSE
          RHS_D%ions2%vix(i1,i2) = RHS_D%ions2%vix(i1,i2)-domain%ions2%viy(i1,i2)*dvi2x_d_plus(i1,i2)
          RHS_D%ions2%ni(i1,i2) = RHS_D%ions2%ni(i1,i2)-domain%ions2%viy(i1,i2)*dni2_d_plus(i1,i2)
        END IF

      ELSE
        matrix_A(:,:)=zero

        !Checked
        matrix_A(1,1)=domain%ions1%viy(i1,i2)
        matrix_A(1,2)=domain%ions1%ni(i1,i2)

        !Checked
        matrix_A(2,2)=domain%ions1%viy(i1,i2)
        matrix_A(2,5)=inv_ne(i1,i2)*domain%Bz(i1,i2)/(mu0*domain%ions1%mass)
      
        !Checked
        matrix_A(3,3)=domain%ions2%viy(i1,i2)
        matrix_A(3,4)=domain%ions2%ni(i1,i2)
      
        !Checked
        matrix_A(4,4)=domain%ions2%viy(i1,i2)
        matrix_A(4,5)=inv_ne(i1,i2)*domain%Bz(i1,i2)/(mu0*domain%ions2%mass)
      
        matrix_A(5,1)=ni2_Bz_ne(i1,i2)*inv_ne(i1,i2)*(domain%ions1%viy(i1,i2)-domain%ions2%viy(i1,i2))
        matrix_A(5,2)=ni1_Bz_ne(i1,i2)
        matrix_A(5,3)=ni1_Bz_ne(i1,i2)*inv_ne(i1,i2)*(domain%ions2%viy(i1,i2)-domain%ions1%viy(i1,i2))
        matrix_A(5,4)=ni2_Bz_ne(i1,i2)
        matrix_A(5,5)=(domain%ions1%ni(i1,i2)*domain%ions1%viy(i1,i2)+    &
          domain%ions2%ni(i1,i2)*domain%ions2%viy(i1,i2))*inv_ne(i1,i2)-  &
          inv_e*dinv_ne_d(i1,i2)*domain%Bz(i1,i2)/mu0


        d1_plus(1)=dni1_d_plus(i1,i2)
        d1_plus(2)=dvi1y_d_plus(i1,i2)
        d1_plus(3)=dni2_d_plus(i1,i2)
        d1_plus(4)=dvi2y_d_plus(i1,i2)
        d1_plus(5)=dBz_d_plus(i1,i2)

        d1_minus(1)=dni1_d_minus(i1,i2)
        d1_minus(2)=dvi1y_d_minus(i1,i2)
        d1_minus(3)=dni2_d_minus(i1,i2)
        d1_minus(4)=dvi2y_d_minus(i1,i2)
        d1_minus(5)=dBz_d_minus(i1,i2)

        d2(1)=dni1_d2(i1,i2)
        d2(2)=dvi1y_d2(i1,i2)
        d2(3)=dni2_d2(i1,i2)
        d2(4)=dvi2y_d2(i1,i2)
        d2(5)=dBz_d2(i1,i2)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ldvl=1
        ldvr=lda
      
        CALL DGEEV('N','V',lda,matrix_A,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,ierr)

        IF (MAXVAL(wi).EQ.zero) THEN
          d1_minus_plus(:,1)=d1_minus
          d1_minus_plus(:,2)=d1_plus
          d1_minus_plus(:,3)=d2

          vr_copy(:,:)=vr(:,:)
          CALL DGETRF(lda,lda,vr,lda,ipvt,ierr)

          CALL DGETRS('N', lda, 3, vr, lda, ipvt, d1_minus_plus, lda, ierr)

          DO i3=1,lda
            IF (wr(i3).GT.zero) THEN
              d1(i3)=d1_minus_plus(i3,1)*wr(i3)-alpha*half*dx*ABS(wr(i3))*d1_minus_plus(i3,3)
            ELSE
              d1(i3)=d1_minus_plus(i3,2)*wr(i3)-alpha*half*dx*ABS(wr(i3))*d1_minus_plus(i3,3)
            END IF
            nu_y_max=MAX(ABS(wr(i3)),nu_y_max)
          END DO

          CALL DGEMV('N',lda,lda,one,vr_copy,lda,d1,1,zero,d1_out,1)

        ELSE

          DO i3=1,lda
            IF (wi(i3).EQ.zero) THEN
              DO i4=1,lda
                vr_C(i4,i3)=DCMPLX(vr(i4,i3))            
              END DO
            ELSE IF (wi(i3).GT.zero) THEN
              DO i4=1,lda
                vr_C(i4,i3)=DCMPLX(vr(i4,i3),vr(i4,i3+1))            
              END DO
            ELSE
              DO i4=1,lda
                vr_C(i4,i3)=DCMPLX(vr(i4,i3-1),-vr(i4,i3))        
              END DO
            END IF
          END DO

          d1_minus_plus_C(:,1)=DCMPLX(d1_minus(:))
          d1_minus_plus_C(:,2)=DCMPLX(d1_plus(:))
          d1_minus_plus_C(:,3)=DCMPLX(d2(:))
   
          vr_C_copy(:,:)=vr_C(:,:)
          CALL ZGETRF(lda,lda,vr_C,lda,ipvt,ierr)
          CALL ZGETRS('N', lda, 3, vr_C, lda, ipvt, d1_minus_plus_C, lda, ierr)

          DO i3=1,lda
            IF (wr(i3).GT.zero) THEN
              d1_C(i3)=d1_minus_plus_C(i3,1)*DCMPLX(wr(i3),wi(i3))  &
                      -alpha*half*dx*ABS(DCMPLX(wr(i3),wi(i3)))*d1_minus_plus_C(i3,3)
            ELSE
              d1_C(i3)=d1_minus_plus_C(i3,2)*DCMPLX(wr(i3),wi(i3))  &
                      -alpha*half*dx*ABS(DCMPLX(wr(i3),wi(i3)))*d1_minus_plus_C(i3,3)
            END IF
            nu_y_max=MAX(ABS(wr(i3)),nu_y_max)
          END DO

          CALL ZGEMV('N',lda,lda,DCMPLX(one),vr_C_copy,lda,d1_C,1,DCMPLX(zero),d1_out_C,1)

          d1_out(:)=DREAL(d1_out_C(:))
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        RHS_D%ions1%ni(i1,i2)  = RHS_D%ions1%ni(i1,i2)-d1_out(1)
        RHS_D%ions1%viy(i1,i2) = RHS_D%ions1%viy(i1,i2)-d1_out(2)
        RHS_D%ions2%ni(i1,i2)  = RHS_D%ions2%ni(i1,i2)-d1_out(3)
        RHS_D%ions2%viy(i1,i2) = RHS_D%ions2%viy(i1,i2)-d1_out(4)
        RHS_D%Bz(i1,i2)        = RHS_D%Bz(i1,i2)-d1_out(5)

        IF (domain%ions1%viy(i1,i2).GT.zero) THEN
          RHS_D%ions1%vix(i1,i2) = RHS_D%ions1%vix(i1,i2)-domain%ions1%viy(i1,i2)*dvi1x_d_minus(i1,i2)
        ELSE
          RHS_D%ions1%vix(i1,i2) = RHS_D%ions1%vix(i1,i2)-domain%ions1%viy(i1,i2)*dvi1x_d_plus(i1,i2)
        END IF

        IF (domain%ions2%viy(i1,i2).GT.zero) THEN
          RHS_D%ions2%vix(i1,i2) = RHS_D%ions2%vix(i1,i2)-domain%ions2%viy(i1,i2)*dvi2x_d_minus(i1,i2)
        ELSE
          RHS_D%ions2%vix(i1,i2) = RHS_D%ions2%vix(i1,i2)-domain%ions2%viy(i1,i2)*dvi2x_d_plus(i1,i2)
        END IF

      END IF

    END DO
    END DO

    RHS_D%nu_y_max=nu_y_max

    END IF
    

    !!! Add force terms with a factor half to fit time-splitting !!!
    IF (t.LT.acc_time) THEN
      ! Initial acceleration of radial expansion and satellite speed
      temp=(pi/(two*domain%acc_time))*sin(pi*t/domain%acc_time)
      DO i2=1,Ny
      DO i1=1,Nx
        RHS_D%v0x(i1,i2)=v_radial*(domain%x(i1)+zero*v_satellite*t)  &
              /DSQRT(R0**2+((domain%x(i1)+zero*v_satellite*t)**2+domain%y(i2)**2)/four)+v_satellite
        RHS_D%v0y(i1,i2)=v_radial*domain%y(i2)  &
              /DSQRT(R0**2+((domain%x(i1)+zero*v_satellite*t)**2+domain%y(i2)**2)/four)
      END DO
      END DO
      RHS_D%ions2%vix(:,:) = RHS_D%ions2%vix(:,:)+half*RHS_D%v0x(:,:)*temp
      RHS_D%ions2%viy(:,:) = RHS_D%ions2%viy(:,:)+half*RHS_D%v0y(:,:)*temp
    END IF

    ni2_Bz_ne(:,:)=domain%ions2%ni(:,:)*domain%Bz(:,:)*inv_ne(:,:)
    charge_mass=domain%ions1%charge/domain%ions1%mass

    RHS_D%ions1%vix(:,:) = RHS_D%ions1%vix(:,:)  &
      +half*charge_mass*(domain%ions1%viy(:,:)-domain%ions2%viy(:,:))*Ni2_Bz_ne(:,:)

    RHS_D%ions1%viy(:,:) = RHS_D%ions1%viy(:,:)  &
            -half*charge_mass*(domain%ions1%vix(:,:)-domain%ions2%vix(:,:))*Ni2_Bz_ne(:,:)

    ni1_Bz_ne(:,:)=domain%ions1%ni(:,:)*domain%Bz(:,:)*inv_ne(:,:)
    charge_mass=domain%ions2%charge/domain%ions2%mass

    RHS_D%ions2%vix(:,:) = RHS_D%ions2%vix(:,:)  &
      +half*charge_mass*(domain%ions2%viy(:,:)-domain%ions1%viy(:,:))*Ni1_Bz_ne(:,:)

    RHS_D%ions2%viy(:,:) = RHS_D%ions2%viy(:,:)  &
            -half*charge_mass*(domain%ions2%vix(:,:)-domain%ions1%vix(:,:))*Ni1_Bz_ne(:,:)

    ! Parallel expansion term
    temp=domain%v_parallel/(domain%R0+domain%v_parallel*t)

    RHS_D%ions2%ni(:,:) = RHS_D%ions2%ni(:,:)  &
     -half*domain%ions2%ni(:,:)*temp

    !!!!!!!!!!!!!!!!!!!!!!!


    IF(.TRUE.) THEN
    CALL dissipation(domain%ions1%ni,dd,domain%dx,domain%dy)
    RHS_D%ions1%ni(:,:) = RHS_D%ions1%ni(:,:) + 600.0d0*dd(:,:)
    
    CALL dissipation(domain%ions2%ni,dd,domain%dx,domain%dy)
    RHS_D%ions2%ni(:,:) = RHS_D%ions2%ni(:,:) + 600.0d0*dd(:,:)

    CALL dissipation(domain%ions1%vix,dd,domain%dx,domain%dy)
    RHS_D%ions1%vix(:,:) = RHS_D%ions1%vix(:,:) + 600.0d0*dd(:,:)

    CALL dissipation(domain%ions2%vix,dd,domain%dx,domain%dy)
    RHS_D%ions2%vix(:,:) = RHS_D%ions2%vix(:,:) + 600.0d0*dd(:,:)

    CALL dissipation(domain%ions1%viy,dd,domain%dx,domain%dy)
    RHS_D%ions1%viy(:,:) = RHS_D%ions1%viy(:,:) + 600.0d0*dd(:,:)

    CALL dissipation(domain%ions2%viy,dd,domain%dx,domain%dy)
    RHS_D%ions2%viy(:,:) = RHS_D%ions2%viy(:,:) + 600.0d0*dd(:,:)

    CALL dissipation(domain%Bz,dd,domain%dx,domain%dy)
    RHS_D%Bz(:,:) = RHS_D%Bz(:,:) + 600.0d0*dd(:,:)
    END IF
!    WRITE(*,*)'nu_x_max',nu_x_max
!    WRITE(*,*)'nu_y_max',nu_y_max

   
    

  IF (.FALSE.) THEN
  dissip(:,:)=3.5d4;
  diss_dd=2.0e3

  inv_ne(:,:)= one/(domain%ions1%ni(:,:)+domain%ions2%ni(:,:)) 
  ni1_Bz_ne(:,:)=domain%ions1%ni(:,:)*domain%Bz(:,:)*inv_ne(:,:)
  ni2_Bz_ne(:,:)=domain%ions2%ni(:,:)*domain%Bz(:,:)*inv_ne(:,:)

  B2_2mu0(:,:)=domain%Bz(:,:)*domain%Bz(:,:)*one_2mu0
  CALL d1x(B2_2mu0,d1x_B2_2mu0,domain%dx)
  CALL d1y(B2_2mu0,d1y_B2_2mu0,domain%dy)

  ! Ions 1, continuity equation
  ni_vix(:,:)=domain%ions1%ni(:,:)*domain%ions1%vix(:,:)
  CALL d1x(ni_vix,d1x_ni_vix,domain%dx)

  ni_viy(:,:)=domain%ions1%ni(:,:)*domain%ions1%viy(:,:)
  CALL d1y(ni_viy,d1y_ni_viy,domain%dy)
  
  CALL dissipation(domain%ions1%ni,dd,domain%dx,domain%dy)
  CALL dissipation(dd,dd2,domain%dx,domain%dy)
  RHS_D%ions1%ni(:,:) =-(d1x_ni_vix(:,:)+d1y_ni_viy(:,:)) &
    +dissip(:,:)*(dd(:,:)-diss_dd*dd2(:,:))


  ! Ions 1, momentum equation
  inv_m=one/domain%ions1%mass
  charge_mass=domain%ions1%charge/domain%ions1%mass

  vi_sq(:,:)=domain%ions1%vix(:,:)*domain%ions1%vix(:,:)
  CALL d1x(vi_sq,dvi_sq,domain%dx)
  CALL d1y(domain%ions1%vix,dvi,domain%dy)

  CALL dissipation(domain%ions1%vix,dd,domain%dx,domain%dy)
  CALL dissipation(dd,dd2,domain%dx,domain%dy)
  RHS_D%ions1%vix(:,:) = -(half*dvi_sq(:,:)  &
    +domain%ions1%viy(:,:)*dvi(:,:))  &
    -inv_m*inv_ne(:,:)*d1x_B2_2mu0(:,:) &
    +charge_mass*(domain%ions1%viy(:,:)-domain%ions2%viy(:,:))*Ni2_Bz_ne(:,:)  &
    +dissip(:,:)*(dd(:,:)-diss_dd*dd2(:,:))


  vi_sq(:,:)=domain%ions1%viy(:,:)*domain%ions1%viy(:,:)
  CALL d1x(vi_sq,dvi_sq,domain%dy)
  CALL d1x(domain%ions1%viy,dvi,domain%dx)

  CALL dissipation(domain%ions1%viy,dd,domain%dx,domain%dy)
  CALL dissipation(dd,dd2,domain%dx,domain%dy)
  RHS_D%ions1%viy(:,:) = -(domain%ions1%vix(:,:)*dvi(:,:)  &
    +half*dvi_sq(:,:))  &
    -inv_m*inv_ne(:,:)*d1y_B2_2mu0(:,:)  &
    -charge_mass*(domain%ions1%vix(:,:)-domain%ions2%vix(:,:))*Ni2_Bz_ne(:,:)  &
    +dissip(:,:)*(dd(:,:)-diss_dd*dd2(:,:))


  ! Ions 2, continuity equation

  ni_vix(:,:)=domain%ions2%ni(:,:)*domain%ions2%vix(:,:)
  CALL d1x(ni_vix,d1x_ni_vix,domain%dx)

  ni_viy(:,:)=domain%ions2%ni(:,:)*domain%ions2%viy(:,:)
  CALL d1y(ni_viy,d1y_ni_viy,domain%dy)

  ! Parallel expansion term
  temp=domain%v_parallel/(domain%R0+domain%v_parallel*t)

  CALL dissipation(domain%ions2%ni,dd,domain%dx,domain%dy)
  CALL dissipation(dd,dd2,domain%dx,domain%dy)
  RHS_D%ions2%ni(:,:) = -(d1x_ni_vix(:,:)+d1y_ni_viy(:,:)) &
    -domain%ions2%ni(:,:)*temp  &
    +dissip(:,:)*(dd(:,:)-diss_dd*dd2(:,:))


  ! Ions 2, momentum equation
  inv_m=one/domain%ions2%mass
  charge_mass=domain%ions2%charge/domain%ions2%mass
  IF (t.LT.acc_time) THEN
    ! Initial acceleration of radial expansion and satellite speed
    temp=(pi/(two*domain%acc_time))*sin(pi*t/domain%acc_time)
  END IF
  
  vi_sq(:,:)=domain%ions2%vix(:,:)*domain%ions2%vix(:,:)
  CALL d1x(vi_sq,dvi_sq,domain%dx)
  CALL d1y(domain%ions2%vix,dvi,domain%dy)

  CALL dissipation(domain%ions2%vix,dd,domain%dx,domain%dy)
  CALL dissipation(dd,dd2,domain%dx,domain%dy)
  RHS_D%ions2%vix(:,:) = -(half*dvi_sq(:,:)  &
    +domain%ions2%viy(:,:)*dvi(:,:))  & 
    -inv_m*inv_ne(:,:)*d1x_B2_2mu0(:,:)  &
    +charge_mass*(domain%ions2%viy(:,:)-domain%ions1%viy(:,:))*Ni1_Bz_ne(:,:) &
    +domain%v0x(:,:)*temp  &
    +dissip(:,:)*(dd(:,:)-diss_dd*dd2(:,:))


  vi_sq(:,:)=domain%ions2%viy(:,:)*domain%ions2%viy(:,:)
  CALL d1y(vi_sq,dvi_sq,domain%dy)
  CALL d1x(domain%ions2%viy,dvi,domain%dx)

  CALL dissipation(domain%ions2%viy,dd,domain%dx,domain%dy)
  CALL dissipation(dd,dd2,domain%dx,domain%dy)
  RHS_D%ions2%viy(:,:) =  -(domain%ions2%vix(:,:)*dvi(:,:)  &
    +half*dvi_sq(:,:))  &
    -inv_m*inv_ne(:,:)*d1y_B2_2mu0(:,:)  &
    -charge_mass*(domain%ions2%vix(:,:)-domain%ions1%vix(:,:))*Ni1_Bz_ne(:,:) &
    +domain%v0y(:,:)*temp  &
    +dissip(:,:)*(dd(:,:)-diss_dd*dd2(:,:))

  ! Magnetic field
  inv_e=one/domain%ions1%charge
  
  DO i2 = 1,Ny
  DO i1 = 1,Nx
    ve_Bz(i1,i2) = ni1_Bz_ne(i1,i2)*domain%ions1%vix(i1,i2)+ni2_Bz_ne(i1,i2)*domain%ions2%vix(i1,i2) &
            -inv_e*inv_ne(i1,i2)*d1y_B2_2mu0(i1,i2)
  END DO
  END DO

  CALL d1x(ve_Bz,diff_ve_Bz,domain%dx)
  RHS_D%Bz(:,:) = -diff_ve_Bz(:,:)

  DO i2 = 1,Ny
  DO i1 = 1,Nx
    ve_Bz(i1,i2) = ni1_Bz_ne(i1,i2)*domain%ions1%viy(i1,i2)+ni2_Bz_ne(i1,i2)*domain%ions2%viy(i1,i2) &
            +inv_e*inv_ne(i1,i2)*d1x_B2_2mu0(i1,i2)
  END DO
  END DO

  CALL d1y(ve_Bz,diff_ve_Bz,domain%dy)
  CALL dissipation(domain%Bz,dd,domain%dx,domain%dy)
  CALL dissipation(dd,dd2,domain%dx,domain%dy)

  RHS_D%Bz(:,:) = RHS_D%Bz(:,:)-diff_ve_Bz(:,:)+dissip(:,:)*(dd(:,:)-diss_dd*dd2(:,:))

  END IF
END SUBROUTINE RHS



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The d1x function calculates d/dx of u_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE d1x(u_in,u_out,dx)
  IMPLICIT NONE

  REAL(8),DIMENSION(:,:), INTENT (in) :: u_in
  REAL(8),DIMENSION(:,:), INTENT (out) :: u_out
  REAL(8), INTENT (in) :: dx
  REAL(8) :: a1
  
  ! MPI
  INTEGER :: request_send_r,request_recv_r,request_send_l,request_recv_l
  INTEGER :: status(MPI_STATUS_SIZE),ierror
  INTEGER :: sendtag_r,recvtag_r,sendtag_l,recvtag_l
  REAL(8) :: sendbuf_r(1:Ny),sendbuf_l(1:Ny)
  REAL(8) :: rc(1:Ny),lc(1:Ny)


  ! Exchange the boundary cells with neighbouring processors
  sendbuf_l(:)=u_in(1,:)
  sendtag_l=1000+my_x
  recvtag_l=1000+my_right
  CALL MPI_IRECV(rc,Ny,MPI_DOUBLE_PRECISION,  &
    my_right,recvtag_l,row_comm,request_recv_l,ierror)
  CALL MPI_ISEND(sendbuf_l,Ny,MPI_DOUBLE_PRECISION,  &
    my_left,sendtag_l,row_comm,request_send_l,ierror)

  sendbuf_r(:)=u_in(Nx,:)
  sendtag_r=2000+my_x
  recvtag_r=2000+my_left
  CALL MPI_IRECV(lc,Ny,MPI_DOUBLE_PRECISION,  &
    my_left,recvtag_r,row_comm,request_recv_r,ierror)
  CALL MPI_ISEND(sendbuf_r,Ny,MPI_DOUBLE_PRECISION,  &
    my_right,sendtag_r,row_comm,request_send_r,ierror)

  a1=one/(two*dx);
  u_out(2:Nx-1,1:Ny)=(u_in(3:Nx,:)-u_in(1:Nx-2,:))*a1;

  CALL MPI_WAIT(request_send_l,status,ierror)
  CALL MPI_WAIT(request_send_r,status,ierror)
  CALL MPI_WAIT(request_recv_l,status,ierror)
  CALL MPI_WAIT(request_recv_r,status,ierror)


  u_out(1,1:Ny)=(u_in(2,:)-lc(:))*a1;
  u_out(Nx,1:Ny)=(rc(:)-u_in(Nx-1,:))*a1;

END SUBROUTINE d1x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The d2x function calculates d2/dx2 of u_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE d2x(u_in,u_out,dx)
  IMPLICIT NONE

  REAL(8),DIMENSION(:,:), INTENT (in) :: u_in
  REAL(8),DIMENSION(:,:), INTENT (out) :: u_out
  REAL(8), INTENT (in) :: dx
  REAL(8) :: a1
  
  ! MPI
  INTEGER :: request_send_r,request_recv_r,request_send_l,request_recv_l
  INTEGER :: status(MPI_STATUS_SIZE),ierror
  INTEGER :: sendtag_r,recvtag_r,sendtag_l,recvtag_l
  REAL(8) :: sendbuf_r(1:Ny),sendbuf_l(1:Ny)
  REAL(8) :: rc(1:Ny),lc(1:Ny)


  ! Exchange the boundary cells with neighbouring processors
  sendbuf_l(:)=u_in(1,:)
  sendtag_l=1000+my_x
  recvtag_l=1000+my_right
  CALL MPI_IRECV(rc,Ny,MPI_DOUBLE_PRECISION,  &
    my_right,recvtag_l,row_comm,request_recv_l,ierror)
  CALL MPI_ISEND(sendbuf_l,Ny,MPI_DOUBLE_PRECISION,  &
    my_left,sendtag_l,row_comm,request_send_l,ierror)

  sendbuf_r(:)=u_in(Nx,:)
  sendtag_r=2000+my_x
  recvtag_r=2000+my_left
  CALL MPI_IRECV(lc,Ny,MPI_DOUBLE_PRECISION,  &
    my_left,recvtag_r,row_comm,request_recv_r,ierror)
  CALL MPI_ISEND(sendbuf_r,Ny,MPI_DOUBLE_PRECISION,  &
    my_right,sendtag_r,row_comm,request_send_r,ierror)

  a1=one/(dx*dx);
  u_out(2:Nx-1,1:Ny)=(u_in(3:Nx,:)-two*u_in(2:Nx-1,:)+u_in(1:Nx-2,:))*a1;

  CALL MPI_WAIT(request_send_l,status,ierror)
  CALL MPI_WAIT(request_send_r,status,ierror)
  CALL MPI_WAIT(request_recv_l,status,ierror)
  CALL MPI_WAIT(request_recv_r,status,ierror)


  u_out(1,1:Ny)=(u_in(2,:)-two*u_in(1,:)+lc(:))*a1;
  u_out(Nx,1:Ny)=(rc(:)-two*u_in(Nx,:)+u_in(Nx-1,:))*a1;

END SUBROUTINE d2x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The d1x_plus function calculates one-sided d/dx of u_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE d1x_plus(u_in,u_out,dx)
  IMPLICIT NONE

  REAL(8),DIMENSION(:,:), INTENT (in) :: u_in
  REAL(8),DIMENSION(:,:), INTENT (out) :: u_out
  REAL(8), INTENT (in) :: dx
  REAL(8) :: a1

  ! MPI
  INTEGER :: request_send,request_recv
  INTEGER :: status(MPI_STATUS_SIZE),ierror
  INTEGER :: sendtag,recvtag
  REAL(8) :: sendbuf(1:Ny)
  REAL(8) :: rc(1:NY),lc(1:NY)


  ! Exchange the boundary cells with neighbouring processors
  sendbuf(:)=u_in(1,:)
  sendtag=100+my_x
  recvtag=100+my_right
  CALL MPI_IRECV(rc,Ny,MPI_DOUBLE_PRECISION,  &
    my_right,recvtag,row_comm,request_recv,ierror)
  CALL MPI_ISEND(sendbuf,Ny,MPI_DOUBLE_PRECISION,  &
    my_left,sendtag,row_comm,request_send,ierror)

  !!! Hide some communication
  a1=one/(dx);

  u_out(1:Nx-1,:)=(u_in(2:Nx,:)-u_in(1:Nx-1,:))*a1;
  !!!


  CALL MPI_WAIT(request_send,status,ierror)
  CALL MPI_WAIT(request_recv,status,ierror)
    

  u_out(Nx,:)=(rc(:)-u_in(Nx,:))*a1;



END SUBROUTINE d1x_plus


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The d1y function calculates d/dy of u_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE d1y(u_in,u_out,dy)
  IMPLICIT NONE

  REAL(8),DIMENSION(:,:), INTENT (in) :: u_in
  REAL(8),DIMENSION(:,:), INTENT (out) :: u_out
  REAL(8), INTENT (in) :: dy

  REAL(8) :: a1

  ! MPI
  INTEGER :: request_send_u,request_recv_u,request_send_d,request_recv_d
  INTEGER :: status(MPI_STATUS_SIZE),ierror
  INTEGER :: sendtag_u,recvtag_u,sendtag_d,recvtag_d
  REAL(8) :: sendbuf_u(1:Nx),sendbuf_d(1:Nx)
  REAL(8) :: uc(1:Nx),dc(1:Nx)


  ! Exchange the boundary cells with neighbouring processors
  sendbuf_d(:)=u_in(:,1)
  sendtag_d=1000+my_y
  recvtag_d=1000+my_up
  CALL MPI_IRECV(uc,Nx,MPI_DOUBLE_PRECISION,  &
    my_up,recvtag_d,col_comm,request_recv_d,ierror)
  CALL MPI_ISEND(sendbuf_d,Nx,MPI_DOUBLE_PRECISION,  &
    my_down,sendtag_d,col_comm,request_send_d,ierror)

  sendbuf_u(:)=u_in(:,Ny)
  sendtag_u=2000+my_y
  recvtag_u=2000+my_down
  CALL MPI_IRECV(dc,Nx,MPI_DOUBLE_PRECISION,  &
    my_down,recvtag_u,col_comm,request_recv_u,ierror)
  CALL MPI_ISEND(sendbuf_u,Nx,MPI_DOUBLE_PRECISION,  &
    my_up,sendtag_u,col_comm,request_send_u,ierror)

  a1=one/(two*dy)
  u_out(1:Nx,2:Ny-1)=(u_in(:,3:Ny)-u_in(:,1:Ny-2))*a1;

  CALL MPI_WAIT(request_send_d,status,ierror)
  CALL MPI_WAIT(request_send_u,status,ierror)
  CALL MPI_WAIT(request_recv_d,status,ierror)
  CALL MPI_WAIT(request_recv_u,status,ierror)

  u_out(1:Nx,1)=(u_in(:,2)-dc(:))*a1;
  u_out(1:Nx,Ny)=(uc(:)-u_in(:,Ny-1))*a1;

END SUBROUTINE d1y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The d2y function calculates d2/dy2 of u_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE d2y(u_in,u_out,dy)
  IMPLICIT NONE

  REAL(8),DIMENSION(:,:), INTENT (in) :: u_in
  REAL(8),DIMENSION(:,:), INTENT (out) :: u_out
  REAL(8), INTENT (in) :: dy

  REAL(8) :: a1

  ! MPI
  INTEGER :: request_send_u,request_recv_u,request_send_d,request_recv_d
  INTEGER :: status(MPI_STATUS_SIZE),ierror
  INTEGER :: sendtag_u,recvtag_u,sendtag_d,recvtag_d
  REAL(8) :: sendbuf_u(1:Nx),sendbuf_d(1:Nx)
  REAL(8) :: uc(1:Nx),dc(1:Nx)


  ! Exchange the boundary cells with neighbouring processors
  sendbuf_d(:)=u_in(:,1)
  sendtag_d=1000+my_y
  recvtag_d=1000+my_up
  CALL MPI_IRECV(uc,Nx,MPI_DOUBLE_PRECISION,  &
    my_up,recvtag_d,col_comm,request_recv_d,ierror)
  CALL MPI_ISEND(sendbuf_d,Nx,MPI_DOUBLE_PRECISION,  &
    my_down,sendtag_d,col_comm,request_send_d,ierror)

  sendbuf_u(:)=u_in(:,Ny)
  sendtag_u=2000+my_y
  recvtag_u=2000+my_down
  CALL MPI_IRECV(dc,Nx,MPI_DOUBLE_PRECISION,  &
    my_down,recvtag_u,col_comm,request_recv_u,ierror)
  CALL MPI_ISEND(sendbuf_u,Nx,MPI_DOUBLE_PRECISION,  &
    my_up,sendtag_u,col_comm,request_send_u,ierror)

  a1=one/(dy*dy)
  u_out(1:Nx,2:Ny-1)=(u_in(:,3:Ny)-two*u_in(:,2:Ny-1)+u_in(:,1:Ny-2))*a1;

  CALL MPI_WAIT(request_send_d,status,ierror)
  CALL MPI_WAIT(request_send_u,status,ierror)
  CALL MPI_WAIT(request_recv_d,status,ierror)
  CALL MPI_WAIT(request_recv_u,status,ierror)

  u_out(1:Nx,1)=(u_in(:,2)-two*u_in(:,1)+dc(:))*a1;
  u_out(1:Nx,Ny)=(uc(:)-two*u_in(:,Ny)+u_in(:,Ny-1))*a1;

END SUBROUTINE d2y


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The d1y_plus function calculates one-sided d/dy of u_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE d1y_plus(u_in,u_out,dy)
  IMPLICIT NONE

  REAL(8),DIMENSION(:,:), INTENT (in) :: u_in
  REAL(8),DIMENSION(:,:), INTENT (out) :: u_out
  REAL(8), INTENT (in) :: dy
  REAL(8) :: a1

  ! MPI
  INTEGER :: request_send,request_recv
  INTEGER :: status(MPI_STATUS_SIZE),ierror
  INTEGER :: sendtag,recvtag
  REAL(8) :: sendbuf(1:Nx)
  REAL(8) :: uc(1:Nx),dc(1:Nx)

  ! Exchange the boundary cells with neighbouring processors
  sendbuf(:)=u_in(:,1)
  sendtag=100+my_y
  recvtag=100+my_up
  CALL MPI_IRECV(uc,Nx,MPI_DOUBLE_PRECISION,  &
    my_up,recvtag,col_comm,request_recv,ierror)
  CALL MPI_ISEND(sendbuf,Nx,MPI_DOUBLE_PRECISION,  &
    my_down,sendtag,col_comm,request_send,ierror)

  !!! Hide some communication
  a1=one/dy;
  u_out(:,1:Ny-1)=(u_in(:,2:Ny)-u_in(:,1:Ny-1))*a1;
  !!!

  CALL MPI_WAIT(request_send,status,ierror)
  CALL MPI_WAIT(request_recv,status,ierror)

  u_out(:,Ny)=(uc(:)-u_in(:,Ny))*a1;

END SUBROUTINE d1y_plus



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The dissipation function calculates d^2/dx^2+d^2/dy^2 of u_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dissipation(u_in,u_out,dx,dy)
  IMPLICIT NONE

  REAL(8),DIMENSION(:,:), INTENT (in) :: u_in
  REAL(8),DIMENSION(:,:), INTENT (out) :: u_out
  REAL(8), INTENT (in) :: dx,dy

  REAL(8) :: a1x,a1y,a00

  ! MPI
  INTEGER :: request_send_r,request_recv_r,request_send_l,request_recv_l
  INTEGER :: request_send_u,request_recv_u,request_send_d,request_recv_d
  INTEGER :: status(MPI_STATUS_SIZE),ierror
  INTEGER :: sendtag_r,recvtag_r,sendtag_l,recvtag_l
  INTEGER :: sendtag_u,recvtag_u,sendtag_d,recvtag_d
  REAL(8) :: sendbuf_r(1:Ny),sendbuf_l(1:Ny)
  REAL(8) :: sendbuf_u(1:Nx),sendbuf_d(1:Nx)
  REAL(8) :: rc(1:Ny),lc(1:Ny)
  REAL(8) :: uc(1:Nx),dc(1:Nx)


  ! Exchange the boundary cells with neighbouring processors
  sendbuf_l(:)=u_in(1,:)
  sendtag_l=1000+my_x
  recvtag_l=1000+my_right
  CALL MPI_ISEND(sendbuf_l,Ny,MPI_DOUBLE_PRECISION,  &
    my_left,sendtag_l,row_comm,request_send_l,ierror)
  CALL MPI_IRECV(rc,Ny,MPI_DOUBLE_PRECISION,  &
    my_right,recvtag_l,row_comm,request_recv_l,ierror)

  sendbuf_r(:)=u_in(Nx,:)
  sendtag_r=2000+my_x
  recvtag_r=2000+my_left
  CALL MPI_IRECV(lc,Ny,MPI_DOUBLE_PRECISION,  &
    my_left,recvtag_r,row_comm,request_recv_r,ierror)
  CALL MPI_ISEND(sendbuf_r,Ny,MPI_DOUBLE_PRECISION,  &
    my_right,sendtag_r,row_comm,request_send_r,ierror)

  sendbuf_d(:)=u_in(:,1)
  sendtag_d=3000+my_y
  recvtag_d=3000+my_up
  CALL MPI_IRECV(uc,Nx,MPI_DOUBLE_PRECISION,  &
    my_up,recvtag_d,col_comm,request_recv_d,ierror)
  CALL MPI_ISEND(sendbuf_d,Nx,MPI_DOUBLE_PRECISION,  &
    my_down,sendtag_d,col_comm,request_send_d,ierror)

  sendbuf_u(:)=u_in(:,Ny)
  sendtag_u=4000+my_y
  recvtag_u=4000+my_down
  CALL MPI_IRECV(dc,Nx,MPI_DOUBLE_PRECISION,  &
    my_down,recvtag_u,col_comm,request_recv_u,ierror)
  CALL MPI_ISEND(sendbuf_u,Nx,MPI_DOUBLE_PRECISION,  &
    my_up,sendtag_u,col_comm,request_send_u,ierror)


  a1x=one/dx**2;
  a1y=one/dy**2;
  a00=two*(a1x+a1y);


  u_out(2:Nx-1,2:Ny-1)=(u_in(2:Nx-1,3:Ny)+u_in(2:Nx-1,1:Ny-2))*a1x  &
    +(u_in(3:Nx,2:Ny-1)+u_in(1:Nx-2,2:Ny-1))*a1y-a00*u_in(2:Nx-1,2:Ny-1);


  CALL MPI_WAIT(request_send_l,status,ierror)
  CALL MPI_WAIT(request_send_r,status,ierror)
  CALL MPI_WAIT(request_recv_l,status,ierror)
  CALL MPI_WAIT(request_recv_r,status,ierror)

  CALL MPI_WAIT(request_send_d,status,ierror)
  CALL MPI_WAIT(request_send_u,status,ierror)
  CALL MPI_WAIT(request_recv_d,status,ierror)
  CALL MPI_WAIT(request_recv_u,status,ierror)

  u_out(1,2:Ny-1)=(u_in(1,3:Ny)+u_in(1,1:Ny-2))*a1x+(u_in(2,2:Ny-1)+lc(2:Ny-1))*a1y-a00*u_in(1,2:Ny-1);
  u_out(Nx,2:Ny-1)=(u_in(Nx,3:Ny)+u_in(Nx,1:Ny-2))*a1x+(rc(2:Ny-1)+u_in(Nx-1,2:Ny-1))*a1y-a00*u_in(Nx,2:Ny-1);
  u_out(2:Nx-1,1)=(u_in(2:Nx-1,2)+dc(2:Nx-1))*a1x+(u_in(3:Nx,1)+u_in(1:Nx-2,1))*a1y-a00*u_in(2:Nx-1,1);
  u_out(2:Nx-1,Ny)=(uc(2:Nx-1)+u_in(2:Nx-1,Ny-1))*a1x+(u_in(3:Nx,Ny)+u_in(1:Nx-2,Ny))*a1y-a00*u_in(2:Nx-1,Ny);

  u_out(1,1)=(u_in(1,2)+dc(1))*a1x+(u_in(2,1)+lc(1))*a1y-a00*u_in(1,1);
  u_out(Nx,1)=(u_in(Nx,2)+dc(Nx))*a1x+(rc(1)+u_in(Nx-1,1))*a1y-a00*u_in(Nx,1);
  u_out(1,Ny)=(uc(1)+u_in(1,Ny-1))*a1x+(u_in(2,Ny)+lc(Ny))*a1y-a00*u_in(1,Ny);
  u_out(Nx,Ny)=(uc(Nx)+u_in(Nx,Ny-1))*a1x+(rc(Ny)+u_in(Nx-1,Ny))*a1y-a00*u_in(Nx,Ny);

END SUBROUTINE dissipation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The Heaviside function                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8) function Heaviside(x)
  IMPLICIT NONE

  ! Argument
  REAL(8),INTENT(in) :: x

  IF (x.GE.zero) THEN
    Heaviside=one
  ELSE
    Heaviside=zero
  END IF
end function Heaviside


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The PeriodicIntegrator2D integrates a periodic 
! function f(x) in two dimensions 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PeriodicIntegrator2D(f,d1,d2,N1,N2,result)
  IMPLICIT NONE

  REAL(8),DIMENSION(:,:),INTENT(in) :: f
  REAL(8),INTENT(in) :: d1,d2
  INTEGER, INTENT(in) :: N1,N2
  REAL(8),INTENT(out) :: result
  REAL(8) :: temp

  INTEGER :: i1,i2

  ! MPI
  INTEGER :: errcode

  result = zero
  DO i2 = 1,N2
  DO i1 = 1,N1
    result = result+f(i1,i2)
  END DO
  END DO
  result = result*d1*d2

  ! Sum all partial sums to proc 0.
  temp=result
  CALL MPI_REDUCE(temp,result,1,  &
    MPI_DOUBLE_PRECISION,MPI_SUM,root,MPI_COMM_WORLD,errcode)

  ! Broadcast the total sum to all procs.
  CALL MPI_BCAST(result, 1, MPI_DOUBLE_PRECISION,root, MPI_COMM_WORLD, errcode)

END SUBROUTINE PeriodicIntegrator2D


END MODULE fluid_numeric_mod
