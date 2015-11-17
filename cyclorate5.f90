!called from gonthierdiffcross.f90

!calculates the cyclotron transition rate, Tau_n0 using the following inputs:photon_ai,electron_pi,B,n.
!And stores it in the output variableTau_final

!CYCLORATE(photon_ai,electron_pi,B,n,Tau_final)



subroutine cyclorate(photon_ai,electron_pi,B,n,Tau_final)

	use constants

	real(dp), intent (IN) :: photon_ai, electron_pi,B
	real(dp), intent (OUT) :: Tau_final
	integer, intent (IN) :: n
	

        real(dp):: electron_ei, J(1:500),interpol(1:6),xvar(1:6)
	real(dp) :: gama, beta, meu, final,z,temp,Bn,final2,tempo,final3,error,B2
	integer::i,tempn

	
	open(unit=2,file='cyclorateIII5.dat')
	open(unit=144,file='cyclorateIIInb3.dat')
	

	call calcTau2(final3,n)	

	meu=cos(photon_ai)


	J=0

	do i=1,3
		call random_number(temp)
		tempo=-3.0D0+temp*1.0D0


		B2=10**(tempo)
	
			
		z=1.0D0+1.0D0/(real(n)*B2)

		
		Bn=B2**(n+1)
		final2=abs((fine_structure)/(compton_wavelength))*final3*Bn
		
		xvar(i)=log10(real(n)*B2)
		interpol(i)=log10(final2)
	end do
	
	do i=4,6
		call random_number(temp)
		tempo=temp*1.0D0


		B2=10**(tempo)
	
			
		z=1.0D0+1.0D0/(real(n)*B2)
		
		electron_ei=sqrt(1.0D0+2*real(n)*B2)

		call calcJ(z,J)	

		beta=-1.0D0*electon_pi/sqrt(1.0D0+2*B2*real(n)*electon_pi*electon_pi)
		gama=1/sqrt(1.0D0-beta*beta)

		call calcTau(z,final,n,J)
		final2=abs((fine_structure*B2)/(compton_wavelength*gama*electon_ei))*final

	
		xvar(i)=log10(real(n)*B2)
		interpol(i)=log10(final2)
	end do
	
	do i=1,6
	end do


	
	
				
		z=1.0D0+1.0D0/(real(n)*B)

		
		final=0.0D0	

		if (z.lt.1.0D0) then
			write(*,*)"ERROR Z<1"
		
		end if

		if(log10(real(n)*B).lt.-2.0D0)then
			write(2,*)"asymptotic"
			Bn=B**(n+1)
			final2=abs((fine_structure)/(compton_wavelength))*final3*Bn
			write(2,*)tempo,B,log10(B), final2, log10(final2)
			write(144,*)real(n)*B,log10(real(n)*B),final2, log10(final2)
			
			Tau_final=final2

		else if(log10(real(n)*B).gt.-2.0D0 .and. log10(real(n)*B).lt. 0.00D0)then
			write(2,*)"interpolated"
			Bn=log10(real(n)*B)
			tempn=6
			call AITKEN(tempn, xvar, interpol, Bn, final2, error)
			write(2,*)tempo,B,log10(B), final2, final2
			write(144,*)real(n)*B,log10(real(n)*B),final2, final2

			Tau_final=10**(final2)
		else
			write(2,*)"normal"
			electron_ei=sqrt(1.0D0+2*real(n)*B)

			call calcJ(z,J)	
	
			beta=-1.0D0*electon_pi/sqrt(1.0D0+2*B*real(n)*electon_pi*electon_pi)
			gama=1/sqrt(1.0D0-beta*beta)

			call calcTau(z,final,n,J)

			final2=abs((fine_structure*B)/(compton_wavelength*gama*electon_ei))*final

			write(2,*)tempo,B,log10(B), final2, log10(final2)
			write(144,*)real(n)*B,log10(real(n)*B),final2, log10(final2)

			Tau_final=final2
		end if


	close(2)
	close(144)
end subroutine cyclorate

subroutine calcTau(z,final,nomber,J)
	use constants

	real(dp), intent(IN):: z,J(1:500)
	real(dp), intent(OUT):: final
	integer, intent(IN):: nomber

	real(dp) :: J1, J2, Jnext, temp, facn ,pown, pownk, fack, maximum,final2
	integer :: k,l



	fack=1.0D0
	pownk=1.0D0

	final=J(nomber)	
	final2=0

	k=0

	
	maximum=10.0D253	
	do k=1,200
	

			fack=fack*real(k)
			pownk=pownk*real(-nomber)
			final2=final

		
			final=final+J(nomber+k)*(pownk/fack)

	
		if(abs(final-final2).lt.1.0D-6)then
			exit
		end if		

	end do
	
	pown=1.0D0
	facn=1.0D0

	do l=1,nomber	
		pown=pown*real(nomber)
		facn=facn*real(l)
	end do

	facn=facn/real(nomber)

	final=final*pown/facn

end subroutine calcTau 



subroutine calcj(z,J)
use constants

	real(dp), intent(IN) :: z
	real(dp), intent(OUT) :: J(1:500)

	real(dp) :: J1, J2, Jnext, temp,term,sum1,fac
	integer :: n,n2

	

	jnext=0.0D0

	temp=log(abs((z+1.0D0)/(z-1.0D0)))
	J1= z-0.5*(z*z-1.0D0)*temp
	J2= -1+ 1.5D0*z*z -0.75D0*z*(z*z-1)*temp



	J(1)=J1
	J(2)=J2

	do n=3,500

	Jnext=z*J2*(2*real(n)+1.0D0)/(real(n)+1.0D0)-J1*(real(n)-1.0D0)/real(n)	

	J(n)=Jnext

	
	j1=j2
	j2=jnext


	end do

	


end subroutine




subroutine calcTau2(final,nomber)
	use constants

	
	real(dp), intent(OUT):: final
	integer, intent(IN):: nomber


	real(dp):: temp1,temp2,fac1,fac2

	integer:: i,j,k

	fac1=1
	fac2=1

	temp1=2.0D0*real(nomber*nomber)
	temp2=1	


	do i=1,nomber+1
		fac1=fac1*i
		temp2=temp2*temp1	
	end do
	
	fac2=fac1
	temp2=temp2/temp1

	do i=nomber+2,2*nomber+1
		fac2=fac2*i
	end do

	
	final=fac1*temp2/fac2


end subroutine calcTau2


SUBROUTINE AITKEN (Num,XI,FI,X,F,DF)
use constants
!
! Subroutine to carry out the Aitken recursions.
!
  INTEGER, INTENT (IN) :: Num
  REAL(dp), INTENT (IN) :: X
  REAL(dp), INTENT (IN):: XI(1:6), FI(1:6)
  REAL(dp), INTENT (OUT) :: F, DF
  
 
 INTEGER, PARAMETER :: NMAX=21
 INTEGER :: I,J
 REAL(dp) :: X1, X2, F1, F2	
 REAL(dp):: FT(1:21)
!
  IF (Num.GT.NMAX) STOP 'Dimension of the data is too large.'
  DO I = 1, Num
    FT(I) = FI(I)
  END DO
!
  DO I = 1, Num-1  
    DO J = 1, Num-I
      X1 = XI(J)
      X2 = XI(J+I)
      F1 = FT(J)
      F2 = FT(J+1)
      FT(J) = (X-X1)/(X2-X1)*F2+(X-X2)/(X1-X2)*F1
    END DO
  END DO
  F = FT(1) 
  DF = (ABS(F-F1)+ABS(F-F2))/2.0
END SUBROUTINE AITKEN

