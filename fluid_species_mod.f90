!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fluid_species_mod.f90   Bengt Eliasson   2017-07-22       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Conatains the class TwoDimFluidSpecies, and constructors/  !
! destructors for this class.                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fluid_species_mod
  USE fluid_param_mod

  IMPLICIT NONE

  type TwoDimFluidSpecies
    ! The particle density and velocity
    REAL(8), DIMENSION(1:Nx,1:Ny) :: ni
    REAL(8), DIMENSION(1:Nx,1:Ny) :: vix,viy
    ! Charge, mass and temperature
    REAL(8) :: charge, mass, temperature
  end type

  CONTAINS

  ! Constructor and Destructor, used if dynamic memory is used ...
  SUBROUTINE SpeciesConstructor(species)
    IMPLICIT NONE
    TYPE (TwoDimFluidSpecies) :: species
  END SUBROUTINE SpeciesConstructor

  SUBROUTINE SpeciesDestructor(species)
    IMPLICIT NONE
    TYPE (TwoDimFluidSpecies) :: species
  END SUBROUTINE SpeciesDestructor
end module fluid_species_mod
