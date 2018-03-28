!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fluid_param_mod.f90    Bengt Eliasson   2014-07-22   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains mathematical constants and      !
! parameters used to set up a run.                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE fluid_param_mod
  IMPLICIT NONE

  !---- Include file for MPI, containing various
  !---- parameters and variables
  include 'mpif.h'
  INTEGER :: errcode,size,rank,cartsize(1:2)
  INTEGER :: cart_comm,row_comm,col_comm 
  INTEGER :: link(1:5),mycoord(1:2),leftcoord(1:2),my_x,my_y
  INTEGER :: my_up,my_down,my_right,my_left
  LOGICAL :: period(1:2)
  INTEGER,PARAMETER :: cart_dim=2         ! Dimension of cartesian map
  LOGICAL,PARAMETER:: reorder  = .true.   ! Reorder nodes if nessesary
  INTEGER,PARAMETER:: UPDOWN   = 1
  INTEGER,PARAMETER:: SIDEWAYS = 0
  INTEGER,PARAMETER:: LEFT     = 1
  INTEGER,PARAMETER:: RIGHT    = 2
  INTEGER,PARAMETER:: DOWN     = 4
  INTEGER,PARAMETER:: UP       = 3
  INTEGER,PARAMETER:: LEFTMOST = 5

  INTEGER,PARAMETER:: root = 0

  ! The number of processors
  INTEGER,PARAMETER :: dim_x=30
  INTEGER,PARAMETER :: dim_y=20
  INTEGER,PARAMETER :: NP=dim_x*dim_y


  !---- Static parameters
  REAL(8),PARAMETER :: pi=    3.141592653589793238463d0
  REAL(8),PARAMETER :: twopi= 6.283185307179586476925d0
  REAL(8),PARAMETER :: halfpi=1.570796326794896619232d0
  REAL(8),PARAMETER :: zero=0.0d0
  REAL(8),PARAMETER :: one=1.0d0
  REAL(8),PARAMETER :: two=2.0d0
  REAL(8),PARAMETER :: three=3.0d0
  REAL(8),PARAMETER :: four=4.0d0
  REAL(8),PARAMETER :: five=5.0d0
  REAL(8),PARAMETER :: six=6.0d0
  REAL(8),PARAMETER :: seven=7.0d0
  REAL(8),PARAMETER :: eight=8.0d0
  REAL(8),PARAMETER :: ten=10.0d0
  REAL(8),PARAMETER :: thirteen=13.0d0
  REAL(8),PARAMETER :: fourteen=14.0d0
  REAL(8),PARAMETER :: twentytwo=22.0d0
  REAL(8),PARAMETER :: half=0.5d0
  REAL(8),PARAMETER :: frac_1_6=0.166666666666666666666d0
  REAL(8),PARAMETER :: frac_5_2=2.5d0 
  REAL(8),PARAMETER :: frac_1_3=0.333333333333333333333d0 
  REAL(8),PARAMETER :: frac_2_3=0.666666666666666666667d0 
  REAL(8),PARAMETER :: frac_4_3=1.333333333333333333333d0 
  REAL(8),PARAMETER :: SqrtThree=1.732050807568877293527d0
  REAL(8),PARAMETER :: SqrtEight=2.828427124746190097604d0
  REAL(8),PARAMETER :: SqrtTen=3.16227766016838d0


  ! File units, beginning with 'out...'
  ! If data files is to be be created, beginning with 'write...':
  ! Yes: .TRUE.
  ! No: .FALSE.

  INTEGER,PARAMETER :: out_Bz=21
  LOGICAL,PARAMETER :: write_Bz=.TRUE.

  INTEGER,PARAMETER :: out_ni1=out_Bz+1
  LOGICAL,PARAMETER :: write_ni1=.TRUE.

  INTEGER,PARAMETER :: out_vi1x=out_ni1+1
  LOGICAL,PARAMETER :: write_vi1x=.TRUE.

  INTEGER,PARAMETER :: out_vi1y=out_vi1x+1
  LOGICAL,PARAMETER :: write_vi1y=.TRUE.

  INTEGER,PARAMETER :: out_ni2=out_vi1y+1
  LOGICAL,PARAMETER :: write_ni2=.TRUE.

  INTEGER,PARAMETER :: out_vi2x=out_ni2+1
  LOGICAL,PARAMETER :: write_vi2x=.TRUE.

  INTEGER,PARAMETER :: out_vi2y=out_vi2x+1
  LOGICAL,PARAMETER :: write_vi2y=.TRUE.

  INTEGER,PARAMETER :: out_x_vx = out_vi2y+1
  LOGICAL,PARAMETER :: write_x_vx=.TRUE.


  LOGICAL,PARAMETER :: write_dump=.TRUE.

  ! Dump file read
  LOGICAL,PARAMETER :: read_dump=.FALSE.

  ! File directory
  CHARACTER(len=*),PARAMETER :: FileDir='data/'

  !---- Numerical settings. ----

  ! The number of grid points
  INTEGER,PARAMETER :: TotalNx=600
  INTEGER,PARAMETER :: TotalNy=600

  ! The number of gridpoints on this processor (don't change).
  INTEGER,PARAMETER :: Nx=TotalNx/dim_x
  INTEGER,PARAMETER :: Ny=TotalNy/dim_y

  ! The number of time steps
  INTEGER,PARAMETER :: Nt=200000000

  ! The domain in x direction
  REAL(8),PARAMETER :: x_start=-6000.0d0
  REAL(8),PARAMETER :: x_end=6000.0d0
  REAL(8),PARAMETER :: y_start=-6000.0d0
  REAL(8),PARAMETER :: y_end=6000.0d0

  ! Numerical dissipation
  REAL(8),PARAMETER :: diss=0.005d0

  ! Random Seed
  REAL(8),PARAMETER :: seed=12345

  !---------- Physical parameters -----------

  ! Constants
  REAL(8),PARAMETER :: mu0 = 4.0d-7*pi
  REAL(8),PARAMETER :: one_2mu0 = one/(two*mu0)

  ! Ions
  REAL(8),PARAMETER :: unit_mass =1.660539040d-27
  REAL(8),PARAMETER :: Ion1_mass = 15.9994d0*unit_mass
  REAL(8),PARAMETER :: Ion1_temperature = 1.0d3
  REAL(8),PARAMETER :: Ion1_charge = 1.60217662d-19
  REAL(8),PARAMETER :: Ion2_mass = 26.981539d0*unit_mass
  REAL(8),PARAMETER :: Ion2_temperature = 1.0d3
  REAL(8),PARAMETER :: Ion2_charge = 1.60217662d-19

  ! Initial velocities, etc.
  REAL(8),PARAMETER :: acc_time = 0.1d0
  REAL(8),PARAMETER :: v_rad_init = 0.0e3
  REAL(8),PARAMETER :: v_radial = 5.0d3
  REAL(8),PARAMETER :: v_satellite = 0.0d3
  REAL(8),PARAMETER :: v_parallel = 1.0d3
  REAL(8),PARAMETER :: R0 = 500.0d0
  REAL(8),PARAMETER :: D0 = 300.0d0
  REAL(8),PARAMETER :: NParticles = 1.0d24

  ! Random Noise (noise level 0.01 means 1% fluctuations)
  LOGICAL,PARAMETER :: ni2Noise = .TRUE.
  LOGICAL,PARAMETER :: vi2Noise = .TRUE.
  LOGICAL,PARAMETER :: BzNoise = .TRUE.
  REAL(8),PARAMETER :: ni2NoiseLevel = 0.01d0
  REAL(8),PARAMETER :: vi2NoiseLevel = 0.01d0
  REAL(8),PARAMETER :: BzNoiseLevel = 0.01d0

  ! Asymmetry - initial y velocity scaled by aspectRatio
  LOGICAL,PARAMETER :: vi2Asym = .FALSE.
  REAL(8),PARAMETER :: aspectRatio = 0.4

  ! Hole - reduce injected plasma density within phi deg of vertical by
  ! a factor of densRed
  LOGICAL,PARAMETER :: ni2Hole = .FALSE.
  REAL(8),PARAMETER :: phiDeg = 15
  REAL(8),PARAMETER :: phiRad = phiDeg*pi/180
  REAL(8),PARAMETER :: densRed = 1.0d3

  ! Density Enhancement
  REAL(8),PARAMETER :: ni1_enhance = 10.0
  
  !------------------------------------------

  ! Stability parameters combining the stability condition for the
  ! Runge-Kutta and difference schemes.
  ! If the time should be calculated adaptively:
  ! IsAdaptive=.TRUE.   : Yes.
  ! IsAdaptive=.FALSE.  : No.
  LOGICAL, PARAMETER :: IsAdaptive = .TRUE.
  
  ! If IsAdaptive=.TRUE., give the CFL number, CFL<1 for stability
  ! to be used by TimeStepCalculator to calculate the time step adaptively.
  REAL(8),PARAMETER :: CFL=0.8d0

  ! If IsAdaptive=.FALSE., give the constant timestep
  REAL(8),PARAMETER :: TimeStep=0.0001d0

  ! Do File write or not
  LOGICAL,PARAMETER :: DoFileWrite=.TRUE.

  ! The number of times to write result
  INTEGER,PARAMETER :: NrPrints=Nt/100
  
  REAL(8),PARAMETER :: dt_print=0.01d0
  REAL(8),PARAMETER :: t_end=2.0d0;

  !------ Variables used in calculations ----

  REAL(8),DIMENSION(Nx,Ny) :: diss_boost

  !===================================================
  ! Initial function
  !===================================================

  ! Select initial function
  INTEGER,PARAMETER :: initial=0

  ! A seed for generating random numbers
  REAL(8) :: harvest1

  CONTAINS

  ! Calculate initial function, species 1
  REAL(8) FUNCTION ni1_initial(x,y)
    REAL(8),INTENT(in) :: x,y

    INTEGER :: i1,i2

    SELECT CASE(initial)
    CASE (0) 
      ni1_initial=ni1_enhance*5.63d11
    CASE (1)
      ni1_initial = one 
    END SELECT
  END FUNCTION ni1_initial

  REAL(8) FUNCTION vi1x_initial(x,y)
    REAL(8),INTENT(in) :: x,y

    INTEGER :: i1,i2

    SELECT CASE(initial)
    CASE (0) 
      vi1x_initial=zero
    CASE (1)
      vi1x_initial = one 
    END SELECT
  END FUNCTION vi1x_initial

  REAL(8) FUNCTION vi1y_initial(x,y)
    REAL(8),INTENT(in) :: x,y

    INTEGER :: i1,i2

    SELECT CASE(initial)
    CASE (0) 
      vi1y_initial=zero
    CASE (1)
      vi1y_initial = one 
    END SELECT
  END FUNCTION vi1y_initial


    ! Calculate initial function, species 2
  REAL(8) FUNCTION ni2_initial(x,y)
    REAL(8),INTENT(in) :: x,y

    REAL(8) ::r,D,t
    INTEGER :: i1,i2

    SELECT CASE(initial)
    CASE (0) 
      D=D0

      r=DSQRT(x**2+y**2)
      ni2_initial=EXP(-(r-R0)**2/D**2)+EXP(-(r+R0)**2/D**2)
    CASE (1)
      ni2_initial = one 
    END SELECT
    IF (ni2Noise) THEN
       ni2_initial=ni2_initial*(1+ni2NoiseLevel*(2*rand()-1))
    ENDIF
    IF (ni2Hole) THEN
       t=atan2(y,x)
       IF ((abs(t-halfpi)<phiRad) .OR. (abs(t+halfpi)<phiRad)) THEN
          ni2_initial=ni2_initial/densRed
       ENDIF
    ENDIF
      
  END FUNCTION ni2_initial

  REAL(8) FUNCTION vi2x_initial(x,y)
    REAL(8),INTENT(in) :: x,y
    REAL(8) :: theta
    INTEGER :: i1,i2

    theta=atan2(y,x)
    SELECT CASE(initial)
    CASE (0)
      vi2x_initial=v_rad_init*cos(theta)
    CASE (1)
      vi2x_initial = one 
    END SELECT
    IF (vi2Noise) THEN
       vi2x_initial=vi2x_initial*(1+ni2NoiseLevel*(2*rand()-1))
    ENDIF
  END FUNCTION vi2x_initial

  REAL(8) FUNCTION vi2y_initial(x,y)
    REAL(8),INTENT(in) :: x,y
    REAL(8) :: theta
    INTEGER :: i1,i2

    theta=atan2(y,x)
    SELECT CASE(initial)
    CASE (0) 
      vi2y_initial=v_rad_init*sin(theta)
    CASE (1)
      vi2y_initial = one 
    END SELECT
    IF (vi2Noise) THEN
      vi2y_initial=vi2y_initial*(1+vi2NoiseLevel*(2*rand()-1))
    ENDIF
    IF (vi2Asym) THEN
       vi2y_initial=vi2y_initial*aspectRatio
    ENDIF
  END FUNCTION vi2y_initial

  ! Calculate initial function, magnetic field
  REAL(8) FUNCTION Bz_initial(x,y)
    REAL(8),INTENT(in) :: x,y

    SELECT CASE(initial)
    CASE (0) 
      Bz_initial=4.0d-5
    CASE (1)
      Bz_initial = one 
    END SELECT
    IF (BzNoise) THEN
       Bz_initial=Bz_initial*(1+BzNoiseLevel*(2*rand()-1))
    ENDIF
  END FUNCTION Bz_initial

END MODULE fluid_param_mod
