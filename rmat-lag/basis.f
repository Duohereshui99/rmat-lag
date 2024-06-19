!!hobasis (complex)
        module basis
            use parameter
            implicit none

!hobasisindex: size(1:nbasis,2)
!phi: basis, (1:Nr,1:Nbasis), each channel (L) in line with Nbasis basis functions
!d2phi: 2nd derivative of basis functions, size the same as that of phi
!phi1: the value of basis at uniformed mesh, size (1:ndiff)
!d2phi1:: the 2nd derivative of basis at uniformed mesh,size (1:ndiff)
!phia: the value of basis at the boundary r=a, phi(a), size: (1:Nbasis)
!phipa: the derivative of basis at the boundary r=a, phi'(a), size: (1:Nbasis)
            contains
!!norm of ho basis
           function normho(nu,n,l)
            implicit none
            real*8::nu,normho
            integer::n,l
            normho=sqrt(sqrt(2*nu**3/pi)*2**(n+2*l+3)*fact(n)
     &      *nu**l/doublefact(2*n+2*l+1))
      !   if(abs(normho)<1e-6) then
      !      write(*,*)'nu,n,l,normho',nu,n,l,normho
      !      write(*,*)'fact(n)',fact(n)
      !      write(*,*)'fact(n+l)=',fact(n+l)
      !      write(*,*)'fact(2n+2l+1)=',fact(2*n+2*l+1)
      !      stop
      !   endif
            end function
!!hobasis
            function ho3d(n,l,nu,r)            !!3d hobasis
            implicit none                           !!nu=mu*omega/(2hbar)
            integer l,n
            real*8::norma,nu
            complex*16::r,ho3d
            norma=normho(nu,n,l)
      !        if (norma<1e-6) then
      !           write(*,*)'ho3d: Norm=0!!!for  nu,n,l',nu,n,l
      !        endif
            ho3d=norma*r**l*exp(-nu*r**2)*
     &       generalized_laguerre(n,l+0.5d0,cmplx(2d0*nu*r**2,kind=8)) !!convert argument x into complex 16 type
            end function ho3d
!! fact            
            function fact(n)           
                  implicit none
                  integer n
                  real*8 fact, dgamma,x
                  x=dfloat(n+1)
                  fact=dgamma(x)    
            end function fact
!!double fact
      function doublefact(n)       !double factorial
            implicit none
            integer::n,i
            real::s,doublefact
            s=1.0
            if(mod(n,2)==0) then
                  do i=n,2,-2
                        s=s*i
                  end do
            else
                  do i=n,1,-2
                        s=s*i
                  end do
            end if
            doublefact=s
      end function
!!generalized_laguerre, non-recursive,alpha=l+1/2
      function generalized_laguerre(n, alpha, x) result(Ln_alpha_x)
            implicit none
            integer, intent(in) :: n
            real*8, intent(in) :: alpha
            complex*16, intent(in) :: x  
            complex*16 :: Ln_alpha_x          
            integer :: i
            complex*16 :: L0, L1, L2
!initialize the first two polynomials
            L0=1.0d0
            L1=1.0d0+alpha-x
ccccccc
            if (n==0) then
                Ln_alpha_x=L0
                return
            endif
ccccccc
            if (n==1) then
                Ln_alpha_x=L1
                return
            endif
ccccccc
            do i=1,n-1
                L2=((2.0d0*i+1.0d0+alpha-x)*L1-(i+alpha)*L0)/(i+1.0d0)
                L0=L1
                L1=L2
            end do
ccccccc
            Ln_alpha_x=L2
        end function generalized_laguerre
ccccccc       
      function LSTFUN(gamma,m,r) !LST transformation s(r)
      implicit none
      real*8::gamma,m
      complex*16::r,LSTFUN
      LSTFUN=(1/((1.d0/r)**m+
     & (1.d0/gamma/sqrt(r))**m))**(1.d0/m)
      end function
ccccccc
      function D1LSTFUN(gamma,m,r)     !s'(r) analytic expression 
      implicit none
      real*8::gamma,m
      complex*16::x,y
      complex*16::r,D1LSTFUN
      x=(1.d0/r)**(m+1.d0)+1.d0/(2.0*r)*(1.d0/(gamma*sqrt(r)))**m
      y=(1.d0/r)**m+(1.d0/(gamma*sqrt(r)))**m
      D1LSTFUN=LSTFUN(gamma,m,r)*x/y
      end function
ccccccc
!!alpha = nu in this module 
      function THOFUNC(n,l,alpha,gamma,m,r)    !alpha: dimensionless parameter in hobasis
      implicit none
      integer::n,l
      real*8::alpha,gamma,m
      complex*16::r,THOFUNC
      THOFUNC=sqrt(D1LSTFUN(gamma,m,r))
     & *ho3d(n,l,alpha,LSTFUN(gamma,m,r))
     & *LSTFUN(gamma,m,r)/r
      end function
ccccccc
c *** ------------------------------------------------
c     Generalized Laguerre function L(n,l+1/2,x)
c *** -----------------------------------------------
      function laguerre(n,l,x)
        implicit none
        integer n,l,p
        real*8:: x,eps,aux,laguerre,fact,aux1,aux2,alpha
        parameter(eps=1e-6)

        alpha=l+0.5d0
        aux1=1d0
        aux2=-x+alpha+1d0


        if (n.eq.0) then
           aux=1d0
        else if (n.eq.1) then
           aux=aux2
        else
           do p=2,n
              aux=((2*p-1+alpha-x)*aux2-(p-1+alpha)*aux1)/p
              aux1=aux2
              aux2=aux
           enddo
        endif
         laguerre=aux


      end function laguerre
      end module