
subroutine totcross(omega_i,B,final)

use constants
	
	real(dp), intent(IN):: omega_i,B

	real(dp), intent(OUT):: final
	

 	real(dp) QM,QD,X
        DIMENSION QM(0:100,0:100),QD(0:100,0:100)

	INTEGER M,N,MM,K
	
 	real(dp) p,Y1,Y2,tau,FACT,temp1,temp2,temp3,temp4
	
	real(dp) tauk,pk,ptaum

	real(dp) H0,H1,G0,G1,G2,one,two,three,four,final
		

	open(unit=1,file='totalcrossdata.dat')
	open(unit=2,file='totalcrossdata2.dat')

	
	do i=1,2000

		call random_number(temp4)

	
		omega_i=temp4*B*2

		y1= (omega_i/(omega_i-B))**2
		y2=(omega_i+B)*(omega_i+B+2*omega_i*B)
		y2=omega_i**2.00/y2

		tau=2.00*B*y2
		tau=1/tau

		p=omega_i/B

		X=1.00+1.00/omega_i

		write(*,*)"y1,y2,p,x,tau",y1,y2,p,x,tau
	


		M=0
		N=50
		MM=100
	
		FACT=1.0D0
		pk=1

		call LQMN(MM,M,N,X,QM,QD)	
	
		H0=QM(M,0)
		H1=QM(M,1)	
	

		G0=QM(M,0)/tau
		G1=QM(M,1)/tau
		G2=QM(M,2)/tau

		tauk=tau

		do k=1,30

			tauk=tauk*tau
			FACT=FACT*k
			pk=pk*(-p)	
	
			H0=H0+(pk*QM(M,k)/FACT)
			H1=H1+(pk*QM(M,k+1)/FACT)
	
			write(1,*)"each term num, den , tot", ((-p)**k)*QM(M,k), FACT , (((-p)**k)*QM(M,k)/FACT)
		
			temp1=1.0D0	
			temp3=1.0D0		
			ptaum=1.0D0
	
			do m=1,k
		
				temp1=temp1*m
				ptaum=ptaum*(-p*tau)	
				temp3=temp3+(ptaum/temp1)

			end do

			G0=G0+QM(M,k)*temp3/tauk
			G1=G1+QM(M,k+1)*temp3/tauk
			G2=G2+QM(M,k+2)*temp3/tauk	

		end do
	
		write(*,*)H0,H1
	
		write(*,*)G0,G1,G2

		temp1=G0-2*X*G1+G2
		temp2=-2*B*Y2+X*G0-G1
		temp3=((X*B+1.0D0)**2+B*B)/(4*B**3)
		four=temp3*(temp2+temp1*omega_i/B)

		three=temp1/(2*B*B)
		
		two=G0-X*G1
		two=two*X/(2*B)

		one=H0-X*H1
		one=one*X*Y1/(1.0D0+2*omega_i)


		final=0.75*thompson_crosssection*(one+two-three-four)

		write(1,*)"one,two,three,four	",one,two,three,four
		write(2,*)omega_i/B,log(abs(final))

			
	end do	

	close(1)
	close(2)

end subroutine totcross
                                 



                        
                                                                  
        SUBROUTINE LQMN(MM,M,N,X,QM,QD) 
!                                                                       
!       ==========================================================      
!       Purpose: Compute the associated Legendre functions of the       
!                second kind, Qmn(x) and Qmn'(x)                        
!       Input :  x  --- Argument of Qmn(x)                              
!                m  --- Order of Qmn(x)  ( m = 0,1,2,úúú )           
!                n  --- Degree of Qmn(x) ( n = 0,1,2,úúú )           
!                mm --- Physical dimension of QM and QD                 
!       Output:  QM(m,n) --- Qmn(x)                                     
!                QD(m,n) --- Qmn'(x)                                    
!       ==========================================================      
        INTEGER MM,M,N                                                         
        DOUBLE PRECISION QM,QD,X
        DIMENSION QM(0:MM,0:N),QD(0:MM,0:N) 

	DOUBLE PRECISION XS,XQ,Q0,Q1,Q10,QF0,QF1,QF2,LS
	
	INTEGER I,J,K,KM



        IF (DABS(X).EQ.1.0D0) THEN 
           DO 10 I=0,M 
           DO 10 J=0,N 
              QM(I,J)=1.0D+300 
              QD(I,J)=1.0D+300 
   10      CONTINUE 
           RETURN 
        ENDIF 
        LS=1 
        IF (DABS(X).GT.1.0D0) LS=-1 
        XS=LS*(1.0D0-X*X) 
        XQ=DSQRT(XS) 
        Q0=0.5D0*DLOG(DABS((X+1.0D0)/(X-1.0D0))) 
        IF (DABS(X).LT.1.0001D0) THEN 
           QM(0,0)=Q0 
           QM(0,1)=X*Q0-1.0D0 
           QM(1,0)=-1.0D0/XQ 
           QM(1,1)=-XQ*(Q0+X/(1.0D0-X*X)) 
           DO 15 I=0,1 
           DO 15 J=2,N 
              QM(I,J)=((2.0D0*J-1.0D0)*X*QM(I,J-1)                      &
     &               -(J+I-1.0D0)*QM(I,J-2))/(J-I)                      
   15      CONTINUE 
           DO 20 J=0,N 
           DO 20 I=2,M 
              QM(I,J)=-2.0D0*(I-1.0D0)*X/XQ*QM(I-1,J)-LS*               &
     &                (J+I-1.0D0)*(J-I+2.0D0)*QM(I-2,J)                 
   20      CONTINUE 
        ELSE 
           IF (DABS(X).GT.1.1D0) THEN 
              KM=40+M+N 
           ELSE 
              KM=(40+M+N)*INT(-1.0-1.8*LOG(X-1.0)) 
           ENDIF 
           QF2=0.0D0 
           QF1=1.0D0 
           DO 25 K=KM,0,-1 
              QF0=((2*K+3.0D0)*X*QF1-(K+2.0D0)*QF2)/(K+1.0D0) 
              IF (K.LE.N) QM(0,K)=QF0 
              QF2=QF1 
   25         QF1=QF0 
           DO 30 K=0,N 
   30         QM(0,K)=Q0*QM(0,K)/QF0 
           QF2=0.0D0 
           QF1=1.0D0 
           DO 35 K=KM,0,-1 
              QF0=((2*K+3.0D0)*X*QF1-(K+1.0D0)*QF2)/(K+2.0D0) 
              IF (K.LE.N) QM(1,K)=QF0 
              QF2=QF1 
   35         QF1=QF0 
           Q10=-1.0D0/XQ 
           DO 40 K=0,N 
   40         QM(1,K)=Q10*QM(1,K)/QF0 
           DO 45 J=0,N 
              Q0=QM(0,J) 
              Q1=QM(1,J) 
              DO 45 I=0,M-2 
                 QF=-2.0D0*(I+1)*X/XQ*Q1+(J-I)*(J+I+1.0D0)*Q0 
                 QM(I+2,J)=QF 
                 Q0=Q1 
                 Q1=QF 
   45      CONTINUE 
        ENDIF 
        QD(0,0)=LS/XS 
        DO 50 J=1,N 
   50      QD(0,J)=LS*J*(QM(0,J-1)-X*QM(0,J))/XS 
        DO 55 J=0,N 
        DO 55 I=1,M 
           QD(I,J)=LS*I*X/XS*QM(I,J)+(I+J)*(J-I+1.0D0)/XQ*QM(I-1,J) 
   55   CONTINUE 
        RETURN 
      END                                           
