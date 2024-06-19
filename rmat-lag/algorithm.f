ccccccc
!algorithm and special functions
        module algorithm
            contains
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SUBROUTINE LEGZO(N,X,W)
C
C       =========================================================
C       Purpose : Compute the zeros of Legendre polynomial Pn(x)
C                 in the interval [0,1], and the corresponding
C                 weighting coefficients for Gauss-Legendre
C                 integration
C       Input :   n    --- Order of the Legendre polynomial
C       Output:   X(n) --- Zeros of the Legendre polynomial
C                 W(n) --- Corresponding weighting coefficients
C       =========================================================
C
C Author: J. M. Jin
C Downloaded from http://jin.ece.illinois.edu/routines/routines.html
          IMPLICIT REAL*8 (A-H,O-Z)
          DIMENSION X(N),W(N)
          data pi/3.1415926535898d0/,one/1/
          N0=(N+1)/2
          DO 45 NR=1,N0
c            Z=COS(pi*(NR-0.25D0)/N)
             Z=COS(pi*(NR-0.25D0)/(N+0.5d0))
  10         Z0=Z
             P=1
             DO 15 I=1,NR-1
  15            P=P*(Z-X(I))
             F0=1
             IF (NR.EQ.N0.AND.N.NE.2*INT(N/2)) Z=0
             F1=Z
             DO 20 K=2,N
                PF=(2-one/K)*Z*F1-(1-one/K)*F0
                PD=K*(F1-Z*PF)/(1-Z*Z)
                F0=F1
  20            F1=PF
             IF (Z.EQ.0) GO TO 40
             FD=PF/P
             Q=0
             DO 35 I=1,NR-1
                WP=1
                DO 30 J=1,NR-1
                   IF (J.NE.I) WP=WP*(Z-X(J))
  30            CONTINUE
  35            Q=Q+WP
             GD=(PD-Q*FD)/P
             Z=Z-FD/GD
             IF (ABS(Z-Z0).GT.ABS(Z)*1.0D-15) GO TO 10
  40         X(NR)=Z
             X(N+1-NR)=-Z
             W(NR)=2/((1-Z*Z)*PD*PD)
  45         W(N+1-NR)=W(NR)
          x(1:n)=(1+x(n:1:-1))/2
          w(1:n)=w(n:1:-1)/2
          RETURN
          END
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!gauss-legendre integral, N:integral mesh number; 
!x1,x2: integral interval;    
!x,w: mesh point and weight
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gauleg(N,x1,x2,X,W)
        IMPLICIT NONE
        INTEGER N
        REAL*8 x1,x2,X(N),W(N)
        REAL*8 z1,z,xm,xl,pp,p3,p2,p1,pi,tol
        INTEGER m,i,j

        pi=acos(-1.0)
        tol=1.E-12

        m=(n+1)/2
        xm=0.5*(x2+x1)
        xl=0.5*(x2-x1)

         DO 10 i=1,m
         z=cos(pi*(i-0.25)/(N+0.5))

 20      CONTINUE
         p1=1.0E0
         p2=0.0E0
         DO 30 j=1,N
          p3=p2
          p2=p1
          p1=((2*j-1)*z*p2-(j-1)*p3)/j
 30      CONTINUE
         pp=N*(z*p1-p2)/(z*z-1.0E0)
         z1=z
         z=z1-p1/pp
         IF( abs(z1-z) .GT. tol) GOTO 20 ! Scheifenende

         X(i) = xm - xl*z
         X(n+1-i) = xm + xl*z
         W(i) = 2.E0*xl/((1.0-z*z)*pp*pp)
         W(n+1-i) = W(i)
 10     CONTINUE
        END SUBROUTINE gauleg 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc
c *** Calculate du(r)/dr using five points derivative formula
c     f(ndim)=function to make derivative
c     h      =step
c     j      =point for derivative
      complex*16 function deriv1(f,h,ndim,j)
        implicit none
        integer ndim,j
        complex*16::f(ndim),h

        if ((j.eq.1).or.(j.eq.2)) then
           deriv1=(-f(j+2)+4d0*f(j+1)-3d0*f(j))/2d0/h
        else if (j.eq.ndim-1) then
           deriv1=(3d0*f(j)-4d0*f(j-1)+f(j-2))/2d0/h
        else if (j.eq.ndim) then
           deriv1=0 !!!CHECK
        else ! five points formula
           deriv1=(f(j-2)-8*f(j-1)+8*f(j+1)-f(j+2))/h/12.
         end if
      end function deriv1
ccccccc        
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
! five points derivative formula for second derivative
! y: function value array
! d2y: second derivative array
! n:size of the array
! uniform grid 
!!complex type 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine second_derivative(y,d2y,n,dx)
       implicit none
       integer,intent(in)::n
       complex*16,intent(in)::dx
       complex*16,dimension(1:n),intent(in)::y
       complex*16,dimension(1:n),intent(out)::d2y
       integer::i 
       d2y(1)=(35.d0/12.d0*y(1)-26.d0/3.d0*y(2)+19.d0/2.d0*y(3)
     &   -14.d0/3.d0*y(4)+11.d0/12.d0*y(5))/(dx**2)
       d2y(2)=(11.d0/12.d0*y(1)-5.d0/3.d0*y(2)+1.d0/2.d0*y(3)
     &  +1.d0/3.d0*y(4)-1.d0/12.d0*y(5))/(dx**2)
       d2y(n-1)=(-1.d0/12.d0*y(N-4)+1.d0/3.d0*y(N-3)+1.d0/2.d0*y(N-2)
     &  -5.d0/3.d0*y(N-1)+11.d0/12.d0*y(N))/(dx**2)
       d2y(n)=(11.d0/12.d0*y(N-4)-14.d0/3.d0*y(N-3)+19.d0/2.d0*y(N-2)
     &  -26.d0/3.d0*y(N-1)+35.d0/12.d0*y(N))/(dx**2)
       do i=3,n-2
       d2y(i)=(-y(i-2)+16.d0*y(i-1)-30.d0*y(i)+ 
     & 16.d0*y(i+1)-y(i+2))/(12.d0*dx**2)
          !  write(*,*) d2y(i)
       end do
       end subroutine second_derivative
ccccccc
!complex interpolation function for uniform grids
      FUNCTION FFC(PP,F,N)
      COMPLEX*16 FFC,F(N)
      REAL*8 PP
      PARAMETER(X=.16666666666667)
      I=PP
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=PP-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFC=(-P2*F(I)+Q*F(I+3))*(P*P1*X)+(P1*F(I+1)-P*F(I+2))*(Q*P2*.5)
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFC=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFC=F(N)
      RETURN
      END function
ccccccc
        end module algorithm