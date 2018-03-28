!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fluid_domain_mod.f90   Bengt Eliasson   2014-07-22        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Conatains the class TwoDimFluidPlasma, and constructors/  !
! destructors for this class.                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fluid_domain_mod
  USE fluid_species_mod

  IMPLICIT NONE

  type TwoDimFluidPlasma 
    ! The time variable
    REAL(8) :: t
    ! The x and eta variables
    REAL(8),DIMENSION(1:Nx) :: x
    REAL(8),DIMENSION(1:Ny) :: y
    
    ! The grid sizes
    REAL(8) :: dx,dy

    ! The magnetic field
    REAL(8),DIMENSION(1:Nx,1:Ny) :: Bz

    ! Particle species
    TYPE (TwoDimFluidSpecies) :: ions1,ions2

    ! Initial velocities and acceleration time
    REAL(8) :: acc_time,R0,v_parallel
    REAL(8),DIMENSION(1:Nx,1:Ny) :: v0x,v0y

    ! Maximum frequency, used to calculate dt
    REAL(8) :: nu_x_max
    REAL(8) :: nu_y_max

    ! Dump number
    INTEGER :: dumpnr
  end type

  CONTAINS

  SUBROUTINE Constructor(domain)
    IMPLICIT NONE
    TYPE (TwoDimFluidPlasma) :: domain
  END SUBROUTINE Constructor


  SUBROUTINE Destructor(domain)
    IMPLICIT NONE
    TYPE (TwoDimFluidPlasma) :: domain
  END SUBROUTINE Destructor

end module fluid_domain_mod
