!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fluid.f90    Bengt Eliasson    2017-07-22            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This main program simulates the evolution of         !
! electrostatic waves in a collisionless plasma,       !
! according to a two-fluid model.                      !
! The program uses a physical object from the class    !
! 'TwoDimFluidPlasma', and some numerical objects     !
! that can operate on objects from this class.         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM fluid

USE fluid_numeric_mod

IMPLICIT NONE

  INTEGER :: k
  REAL(8) :: t,dt,start,finish,next_print,ddt,vmax,t_copy,temp
  REAL(8) :: v_frame_local,dx_frame_local,ni2_frame_local,x_frame
  REAL(8) :: vx_frame_global,dx_frame_global,ni2_frame_global,dx_frame_global_old

  INTEGER :: i1,i2
  TYPE (TwoDimFluidPlasma) :: domain1
  CHARACTER(len=4) :: str
  INTEGER :: iocheck
  CHARACTER :: direction

  !!!! MPI
  REAL(8),DIMENSION(TotalNx*TotalNy) :: DataGlobal
  INTEGER :: LocArraySize


  ! Allocate memory for the domain
  CALL Constructor(domain1)


  ! Read from dump file
  IF (read_dump) THEN

    ! Initialize the initial conditions and various constants, et.c.
    CALL Domain_Initializer(domain1,t)

    IF (rank.EQ.root) THEN
      WRITE(*,*)'Dumpnr:'
      READ(*,*)domain1%dumpnr
      domain1%dumpnr=domain1%dumpnr+1000
    END IF

    ! Broadcast the dumpnr to all procs.
    CALL MPI_BCAST(domain1%dumpnr, 1, MPI_INTEGER,root, MPI_COMM_WORLD, errcode)

    ! Write dumpnr to str
    WRITE(str,1000)domain1%dumpnr
    1000 FORMAT(I4)

    LocArraySize=Nx*Ny

    IF (rank.EQ.root) THEN
      OPEN (UNIT=out_Bz,FILE=FileDir//'Bz_'//str//'.dat',ACCESS='SEQUENTIAL', &
      status='UNKNOWN',form='FORMATTED',ACTION='READ',IOSTAT=iocheck)
      READ(out_Bz,*,end=100)DataGlobal
100   CLOSE(out_Bz)
    END IF

    CALL MPI_SCATTER(DataGlobal,LocArraySize,MPI_DOUBLE_PRECISION,domain1%Bz, &
        LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)


    IF (rank.EQ.root) THEN
      OPEN (UNIT=out_ni1,FILE=FileDir//'ni1_'//str//'.dat',ACCESS='SEQUENTIAL', &
      status='UNKNOWN',form='FORMATTED',ACTION='READ',IOSTAT=iocheck)
      READ(out_ni1,*,end=101)DataGlobal
101   CLOSE(out_ni1)
    END IF

    CALL MPI_SCATTER(DataGlobal,LocArraySize,MPI_DOUBLE_PRECISION,domain1%ions1%ni, &
        LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)


    IF (rank.EQ.root) THEN
      OPEN (UNIT=out_ni2,FILE=FileDir//'ni2_'//str//'.dat',ACCESS='SEQUENTIAL', &
      status='UNKNOWN',form='FORMATTED',ACTION='READ',IOSTAT=iocheck)
      READ(out_ni2,*,end=102)DataGlobal
102   CLOSE(out_ni2)
    END IF

    CALL MPI_SCATTER(DataGlobal,LocArraySize,MPI_DOUBLE_PRECISION,domain1%ions2%ni, &
        LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)


    IF (rank.EQ.root) THEN
      OPEN (UNIT=out_vi1x,FILE=FileDir//'vi1x_'//str//'.dat',ACCESS='SEQUENTIAL', &
      status='UNKNOWN',form='FORMATTED',ACTION='READ',IOSTAT=iocheck)
      READ(out_vi1x,*,end=103)DataGlobal
103   CLOSE(out_vi1x)
    END IF

    CALL MPI_SCATTER(DataGlobal,LocArraySize,MPI_DOUBLE_PRECISION,domain1%ions1%vix, &
        LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)


    IF (rank.EQ.root) THEN
      OPEN (UNIT=out_vi1y,FILE=FileDir//'vi1y_'//str//'.dat',ACCESS='SEQUENTIAL', &
      status='UNKNOWN',form='FORMATTED',ACTION='READ',IOSTAT=iocheck)
      READ(out_vi1y,*,end=104)DataGlobal
104   CLOSE(out_vi1y)
    END IF

    CALL MPI_SCATTER(DataGlobal,LocArraySize,MPI_DOUBLE_PRECISION,domain1%ions1%viy, &
        LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)


    IF (rank.EQ.root) THEN
      OPEN (UNIT=out_vi2x,FILE=FileDir//'vi2x_'//str//'.dat',ACCESS='SEQUENTIAL', &
      status='UNKNOWN',form='FORMATTED',ACTION='READ',IOSTAT=iocheck)
      READ(out_vi2x,*,end=105)DataGlobal
105   CLOSE(out_vi2x)
    END IF

    CALL MPI_SCATTER(DataGlobal,LocArraySize,MPI_DOUBLE_PRECISION,domain1%ions2%vix, &
        LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)


    IF (rank.EQ.root) THEN
      OPEN (UNIT=out_vi2y,FILE=FileDir//'vi2y_'//str//'.dat',ACCESS='SEQUENTIAL', &
      status='UNKNOWN',form='FORMATTED',ACTION='READ',IOSTAT=iocheck)
      READ(out_vi2y,*,end=106)DataGlobal
106   CLOSE(out_vi2y)
    END IF

    CALL MPI_SCATTER(DataGlobal,LocArraySize,MPI_DOUBLE_PRECISION,domain1%ions2%viy, &
        LocArraySize,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,errcode)


    IF (rank.EQ.root) THEN
      OPEN (UNIT=out_x_vx,FILE=FileDir//'x_vx_'//str//'.dat',ACCESS='SEQUENTIAL', &
      status='UNKNOWN',form='FORMATTED',ACTION='READ',IOSTAT=iocheck)
      READ(out_x_vx,*,end=107)x_frame, vx_frame_global,dx_frame_global_old, domain1%nu_x_max,domain1%nu_y_max
107   CLOSE(out_x_vx)
    END IF

    ! Broadcast the frame velocity to all procs.
    CALL MPI_BCAST(vx_frame_global, 1, MPI_DOUBLE_PRECISION,root, MPI_COMM_WORLD, errcode)

    ! Broadcast the maximum eigenvalue to all procs for calculation of dt.
    CALL MPI_BCAST(domain1%nu_x_max, 1, MPI_DOUBLE_PRECISION,root, MPI_COMM_WORLD, errcode)
    CALL MPI_BCAST(domain1%nu_y_max, 1, MPI_DOUBLE_PRECISION,root, MPI_COMM_WORLD, errcode)

    domain1%t=0.01d0*DBLE(domain1%dumpnr-1000)

    t=domain1%t


    next_print=t+dt_print
    k=0

  ELSE

    domain1%dumpnr=1000

    ! Initialize the initial conditions and various constants, et.c.
    CALL Domain_Initializer(domain1,t)

    CALL TimeStepCalculator(domain1,dt)


    k=0
    CALL Domain_StatisticsWriter(domain1,k,t,dt)


    vx_frame_global=zero
    dx_frame_global=zero
    dx_frame_global_old=zero
    x_frame=zero
    IF (rank.EQ.root) THEN
      ! Write dumpnr to str
      WRITE(str,1000)domain1%dumpnr

      IF (write_x_vx) THEN
        OPEN (UNIT=out_x_vx,FILE=FileDir//'x_vx_'//str//'.dat',ACCESS='SEQUENTIAL',   &
        FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)

        WRITE (out_x_vx,*) x_frame, vx_frame_global,dx_frame_global_old, domain1%nu_x_max,domain1%nu_y_max

        CLOSE(out_x_vx)
!      CALL execute_command_line ('gzip -9 -f '//FileDir//'ni2_'//str//'_'//str2//'_'//str3//'.dat', wait=.false.)
      END IF
    END IF

    next_print=dt_print
  END IF



  DO k = 1,Nt
    IF (rank.EQ.root) THEN
      CALL cpu_time(start)
    END IF
    ! Calculate the time step
    CALL TimeStepCalculator(domain1,dt)

    IF (t+dt.LT.next_print) THEN
            
      ! Runge-Kutta steps, first in x-direction then in y-direction
      t_copy=t
      CALL RungeKutta_TimeStepper(domain1,t_copy,dt,'x',vx_frame_global)
      t_copy=t
      CALL RungeKutta_TimeStepper(domain1,t_copy,dt,'y',vx_frame_global)
      t=t+dt

      domain1%t=t

      !!! Calculate the center of mass velocity of the beam ions
      dx_frame_local=zero
      ni2_frame_local=zero
      DO i2=1,Ny
      DO i1=1,Nx
        dx_frame_local=dx_frame_local+domain1%x(i1)*domain1%ions2%ni(i1,i2)
        ni2_frame_local=ni2_frame_local+domain1%ions2%ni(i1,i2)
      END DO
      END DO
      

      ! Calculate the sum among processors
      CALL MPI_REDUCE(dx_frame_local,dx_frame_global,1,  &
        MPI_DOUBLE_PRECISION,MPI_SUM,root,MPI_COMM_WORLD,errcode)

      ! Calculate the sum among processors
      CALL MPI_REDUCE(ni2_frame_local,ni2_frame_global,1,  &
        MPI_DOUBLE_PRECISION,MPI_SUM,root,MPI_COMM_WORLD,errcode)

      IF (rank.EQ.root) THEN
        dx_frame_global=dx_frame_global/ni2_frame_global
        vx_frame_global=vx_frame_global+(dx_frame_global-dx_frame_global_old)/dt
        dx_frame_global_old=dx_frame_global
        x_frame=x_frame+vx_frame_global*dt
      END IF

      ! Broadcast the result to all procs.
      CALL MPI_BCAST(vx_frame_global, 1, MPI_DOUBLE_PRECISION,root, MPI_COMM_WORLD, errcode)
!      CALL MPI_BCAST(ni2_frame_global, 1, MPI_DOUBLE_PRECISION,root, MPI_COMM_WORLD, errcode)

    ELSE
      ddt=next_print-t
      t_copy=t
      CALL RungeKutta_TimeStepper(domain1,t_copy,ddt,'x',vx_frame_global)
      t_copy=t
      CALL RungeKutta_TimeStepper(domain1,t_copy,ddt,'y',vx_frame_global)
      t=next_print
      domain1%t=t
      next_print=next_print+dt_print

      IF (ddt.GT.1e-10) THEN

        !!! Calculate the center of mass velocity of the beam ions
        dx_frame_local=zero
        ni2_frame_local=zero
        DO i2=1,Ny
        DO i1=1,Nx
          dx_frame_local=dx_frame_local+domain1%x(i1)*domain1%ions2%ni(i1,i2)
          ni2_frame_local=ni2_frame_local+domain1%ions2%ni(i1,i2)
        END DO
        END DO
      

        ! Calculate the sum among processors
        CALL MPI_REDUCE(dx_frame_local,dx_frame_global,1,  &
          MPI_DOUBLE_PRECISION,MPI_SUM,root,MPI_COMM_WORLD,errcode)

        ! Calculate the sum among processors
        CALL MPI_REDUCE(ni2_frame_local,ni2_frame_global,1,  &
          MPI_DOUBLE_PRECISION,MPI_SUM,root,MPI_COMM_WORLD,errcode)

        IF (rank.EQ.root) THEN
          dx_frame_global=dx_frame_global/ni2_frame_global
          vx_frame_global=vx_frame_global+(dx_frame_global-dx_frame_global_old)/ddt
          dx_frame_global_old=dx_frame_global
          x_frame=x_frame+vx_frame_global*ddt
        END IF

        ! Broadcast the result to all procs.
        CALL MPI_BCAST(vx_frame_global, 1, MPI_DOUBLE_PRECISION,root, MPI_COMM_WORLD, errcode)
      END IF

      domain1%dumpnr=domain1%dumpnr+1
      CALL Domain_StatisticsWriter(domain1,k,t,dt)

      temp=domain1%nu_x_max
      CALL MPI_REDUCE(temp,domain1%nu_x_max,1,  &
        MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,errcode)

      temp=domain1%nu_y_max
      CALL MPI_REDUCE(temp,domain1%nu_y_max,1,  &
        MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,errcode)

      IF (rank.EQ.root) THEN
        ! Write dumpnr to str
        WRITE(str,1000)domain1%dumpnr

        IF (write_x_vx) THEN
          OPEN (UNIT=out_x_vx,FILE=FileDir//'x_vx_'//str//'.dat',ACCESS='SEQUENTIAL',   &
          FORM='FORMATTED',ACTION='WRITE',IOSTAT=iocheck)

          WRITE (out_x_vx,*) x_frame,vx_frame_global,dx_frame_global_old,  &
                  domain1%nu_x_max,domain1%nu_y_max

          CLOSE(out_x_vx)
!      CALL execute_command_line ('gzip -9 -f '//FileDir//'ni2_'//str//'_'//str2//'_'//str3//'.dat', wait=.false.)
        END IF

      END IF
    END IF

    IF (rank.EQ.root) THEN
      CALL cpu_time(finish)
      WRITE(*,906)'k=',k,' dt=',dt,' t=',t,' nextdump=',next_print,' cpu-time=',finish-start
    END IF

    IF (ABS(t-t_end).LT.1.0d-10) THEN
      STOP 0
    END IF
  END DO



  ! Clean up
  CALL Domain_Finalizer(domain1)

  ! Deallocate memory
  CALL Destructor(domain1)

  ! CALL MPI_FINALIZE(errcode)

  STOP 0
  906 FORMAT (A,I5,A,F11.8,A,F13.7,A,F13.7,A,F10.7)

END PROGRAM
