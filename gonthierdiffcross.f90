!subroutine cyclorate is used from cylcorate5.f90


module constants

	integer, parameter :: dp = selected_real_kind(15, 307)

	real(dp), parameter :: fine_structure=1.0D0
	real(dp), parameter :: plancks_constant=1.0D0
	real(dp), parameter :: lightspeed=1.0D0
	real(dp), parameter :: thompson_crosssection=1.0D0
	real(dp), parameter :: electon_mass=1.0D0
	real(dp) :: compton_wavelength=plancks_constant/(electon_mass*lightspeed)

	REAL(dp), PARAMETER::omega_min=1, omega_max=1.8 !maximum (1.8 for level 10)  and minimum energies of photon in MeV
	real(dp),PARAMETER:: r_c=10 !critical radius in magnetic field

	INTEGER, PARAMETER:: montecarlo=100
	INTEGER, PARAMETER:: maxlevel=15 !maximum landau level considered in the simulation 
end module constants

program diffcross
	use constants

	!short hand for variables electron/photon a: angle e: energy p: momentum lev: energy level (integer). j & l are initial are final 		!electrons, i & f are initial and final photons.	
		
	real(kind=dp) :: B, photon_ai, photon_ei,photon_af
	integer :: electon_levj, electron_levl,s,loop1

	real(kind=dp) :: z, energy, linewidth,photon_ef,r,electron_ej,electron_el,electon_pl,electon_el
	real(kind=dp) :: temp, energy_perp, delta, Ns, Ds, F, diff_crosssection, Tave, diff_cross_peak, tot_crosssection
	real(kind=dp) :: temp1,temp2,temp3,temp4,temp5,temp6

		

	B=3
	photon_ei=1.002D0*B
	photon_ai=0
		
		
	electron_levj=0
	electron_levl=1

	open(unit=1,file='gonthierB3w1.002.dat')	
		
	do loop1=1,3000		
	
	
		call random_number(temp1)
		temp1=-4.0D0+temp1*5.0D0		
		
		photon_af=10**(temp1) !randomly choose photon final angle between 0 and 2pi?
		r=1.0D0/(1.0D0+photon_ei*(1-cos(photon_af))) 
		

		photon_ef=2.0D0*(photon_ei-real(electon_levl)*B)*r !final energy of photon from eqn 3
		temp1=-2.0D0*(photon_ei-(electon_levl)*B)
		temp2=(r*sin(photon_af))**2
		temp3=sqrt(1.0D0-temp1*temp2)
		photon_ef=photon_ef/(1.0D0+ temp3 )
	
		
		electron_pl=photon_ei-photon_ef*cos(photon_af) !final momentum/energy of electron from eqn 5
		electron_el=1+photon_ei-photon_ef

		temp=photon_ei*photon_ef*(1.0D0-cos(photon_af)) 
		energy_perp= sqrt(1.0D0+2.0D0*B)

		temp1=photon_ei*(2.0D0*photon_ei-photon_ef-temp) !calculating F from eqn 40 
		temp2= (photon_ef*sin(photon_af))**2
		temp2=-temp2/(2*B)
		F=photon_ef*photon_ef*exp(temp2)
		F=0.375D0*thompson_crosssection*F/temp1

		delta=2.0D0*(photon_ei-B)	

		z=1.0D0+1.0D0/B
		linewidth=0
		energy=sqrt(1+2*B)
	
		call CYCLORATE(photon_ai,electron_pl,B,electron_levl,linewidth) !calculating TAU_n0, from Gonthier (2005)	
	
!		linewidth=linewidth*finestructure*B/energy	


		Tave=2.0D0*((1.0D0+photon_ei)*(photon_ei-temp)-photon_ef) ! Calculating T_ave from eqn 77

!		write(*,*)"temp",temp		
!		write(*,*)"F",F
!		write(*,*)"delta",delta
!		write(*,*)"energy_perp",energy_perp
!		write(*,*)"Tave",Tave
!		write(*,*)"photon_ef",photon_ef
!		write(*,*)"line width",linewidth	
	
		temp4=2.0D0*delta*(1.0D0+cos(photon_af))*((1.0D0+photon_ei)*temp-photon_ei*(photon_ei+photon+photon_ef))
		temp6=4.0D0*(photon_ei-B)**2
			
		diff_crosssection=0.0D0	
		s=-1
		
		do s=-1,1,2
	
			temp1=energy_perp+real(s)  !Calculating N_s and D_s from eqn 81
			temp2=energy_perp-real(s)
			temp3=temp2*temp2*((2*energy_perp+real(s))*Tave-real(s)*(temp1*temp1)*(photon_ei-photon_ef))
			temp5=real(s)*temp4
		
			Ns=temp3-temp5
		
	
			Ds=temp6+(temp2*(1.0D0+B)*linewidth/energy_perp)**2
		
			!write(*,*)"Ns,Ds",Ns,Ds

			diff_crosssection=diff_crosssection+Ns/Ds
		end do
	
		diff_crosssection=diff_crosssection*F/(2*energy_perp**3)
	
		write(1,*)loop1,log10(photon_af),diff_crosssection	
	
		!write(*,*)"z,final",z,diff_crosssection
	
	end do
	close(1)
	
	diff_cross_peak=2*F*(Tave-(photon_ei-photon_ef))/(((1.0D0+B)*linewidth)**2)

end program diffcross



