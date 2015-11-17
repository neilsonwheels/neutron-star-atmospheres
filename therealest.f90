
module constants

	integer, parameter :: dp = selected_real_kind(15, 307)

	real(dp), parameter :: fine_structure=1.0D0
	real(dp), parameter :: plancks_constant=1.0D0
	real(dp), parameter :: lightspeed=1.0D0
	real(dp), parameter :: thompson_crosssection=1.0D0
	real(dp), parameter :: electon_mass=1.0D0
	real(dp) :: compton_wavelength=plancks_constant/(electon_mass*lightspeed)

	real(dp), parameter::omega_min=1, omega_max=1.8 !maximum (1.8 for level 10)  and minimum energies of photon in MeV
	real(dp), parameter:: r_c=10 !critical radius in magnetic field
	real(dp), parameter:: electron_density=1000 !density of electrons

	INTEGER, PARAMETER:: montecarlo=100
	INTEGER, PARAMETER:: maxlevel=15 !maximum landau level considered in the simulation 

end module constants


program gettingtooreal
use constants

real(dp):: photon_ie,photon_ix,photon_iy,photon_iz,photon_ia,photon_imeu,photon_iphi
real(dp):: B

real(dp):: tscs !total scattering cross_section
real(dp):: mean_free_path !mean free path travelled by the photon

B=3

photon_ei=inject()

call angle_pos_init(photon_iy,photon_iz,photon_ia,photon_imeu,photon_iphi)

photon_ia=0  !setting photon initial angle to be 0 by hand
photon_imeu=1 !setting cos(photon_ia)=cos(0)=1

!part(a): getting the mean free path of the process

call totcross(omega_i,B,tscs)

mean_free_path= 1/abs(electron_density*tscs)

!propagate the photon by this path, along the magnetic field

photon_iz=photon_iz+mean_free_path



!part(b): comparing between cyclotron absorption and compton scattering

!part(c): lorentz transform to ERF





!-------------------------------------------!
real(dp) FUNCTION inject()
!-------------------------------------------!
!find the spectral weight associated with a photon
use constants
real(dp):: r1
CALL RANDOM_NUMBER(r1)
inject= omega_min+ (omega_max-omega_min)*r1

RETURN
END

!--------------------------------------------!
SUBROUTINE angle_pos_init(x_inj,y_inj,z_inj,theta_inj,phi_inj,meu_inj)
!--------------------------------------------!
!initialize the injection position and direction of propagation at the time of injection.
USE constants
real(dp), INTENT(OUT):: x_inj,y_inj,z_inj,theta_inj,phi_inj,meu_inj
real(dp):: r2, r3, r



CALL RANDOM_NUMBER(r2)
r=r_c*SQRT(r2)

CALL RANDOM_NUMBER(r3)
phi_inj=2*Pi*r3

x_inj=r*cos(phi_inj)
y_inj=r*sin(phi_inj)
z_inj=0.00

meu_inj=r3

CALL RANDOM_NUMBER(r3)

theta_inj=Pi*r3/2


RETURN
END 

!-------------------------------------------!
subroutine mean_free_path(local_freepath,freepath)
!-------------------------------------------!

use constants


real(dp), INTENT(IN):: local_freepath
real(dp), INTENT(OUT):: freepath

real(dp):: r2


CALL RANDOM_NUMBER(r2)

freepath=-local_freepath*log(r2)

RETURN
END


