ccccccc
        module rmatmod
                use parameter
                use channels
                use matinv
                use coulfunc
                use deltaf
                use coulvar
                use wtkvar
                use rmatvar
                use mesh
                use system
                use potential
                use potvar
                implicit none
        contains
ccccccc
        subroutine rmat_int()
                implicit none
ccccccc
                integer::i,k             
ccccccc
                if(allocated(WTK)) deallocate(WTK)
                if(allocated(WTKP)) deallocate(WTKP)
                if(allocated(xle)) deallocate(xle)
                if(allocated(wle)) deallocate(wle)
                if(allocated(phia)) deallocate(phia)
                if(allocated(Cmat)) deallocate(Cmat)
                if(allocated(Vcouple)) deallocate(Vcouple)
                if(allocated(T)) deallocate(T)
                if(allocated(B_i)) deallocate(B_i)
                if(allocated(Ech)) deallocate(Ech)
                if(allocated(C)) deallocate(C)
                if(allocated(Rmat)) deallocate(Rmat)
                if(allocated(Z_O)) deallocate(Z_O)
                if(allocated(Z_I)) deallocate(Z_I)
                if(allocated(Smat)) deallocate(Smat)
ccccccc
                allocate(WTK(1:beta%nchmax+1),WTKP(1:beta%nchmax+1))
                allocate(xle(1:nr),wle(1:nr))
                allocate(phia(1:nr))
                allocate(Cmat(1:nr,1:nr,1:beta%nchmax,1:beta%nchmax))
                allocate(Vcouple(1:nr,1:nr,1:beta%nchmax,1:beta%nchmax))
                allocate(T(1:nr,1:nr,1:beta%nchmax))
                allocate(B_i(1:beta%nchmax))
                allocate(Ech(1:nr,1:nr,1:beta%nchmax))
                allocate(C(1:nr*beta%nchmax,1:nr*beta%nchmax))
                allocate(Rmat(1:beta%nchmax,1:beta%nchmax),Smat(1:beta%nchmax,1:beta%nchmax))
                allocate(Z_O(1:beta%nchmax,1:beta%nchmax),Z_I(1:beta%nchmax,1:beta%nchmax))
ccccccc
                call LEGZO(nr,xle,wle)
ccccccc
                do k=1,nr
                    phia(k)=(-1)**(nr+k)*sqrt(1d0/rmax/xle(k)/(1d0-xle(k)))
                end do
ccccccc              
                do i=1,beta%nchmax
                      if(E>Ec(i)) then
                                B_i(i)=0d0
                        else
                        ki=sqrt(2d0*mu*abs(E-Ec(i))/hbarc**2)
                        eta=z1*z2*e2*mu/hbarc**2/ki
                        call WHIT(eta,rmax,ki,E,int(lc(i),4),WTK,WTKP,0)
                        B_i(i)=2*ki*rmax*WTKP(i)/WTK(i)
                      end if
                end do
ccccccc
        end subroutine rmat_int
ccccccc
!this subroutine gives the potential
        subroutine getpot(str)
            implicit none
            integer::i
            character(len=*)::str
ccccccc
            real*8::xx,zz,vcen,vtens,vls
ccccccc
        if(allocated(Vc)) deallocate(Vc)
        allocate(Vc(1:nr,1:beta%nchmax,1:beta%nchmax))    
ccccccc
        select case(str)
ccccccc
!t: tensor force term included in the coupled pot for neutron-proton scattering
!!(only for 2 channels l=0,2)
            case('t') 
c Reid neutron-proton potential (T=1, soft core)
        do i=1,nr
            xx=0.7d0*xle(i)*rmax
            zz=exp(-xx)
            vcen=(-10.463d0*zz+105.468d0*zz**2-3187.8d0*zz**4+9924.3d0*zz**6)/xx
            vtens=-10.463d0*((1+3/xx+3/xx**2)*zz-(12/xx+3/xx**2)*zz**4)/xx+351.77d0*zz**4/xx-1673.5d0*zz**6/xx
            vls=708.91d0*zz**4/xx-2713.1d0*zz**6/xx
            Vc(i,1,1)=vcen-2*(beta%j_tot-1)*vtens/(2*beta%j_tot+1)+(beta%j_tot-1)*vls
            Vc(i,1,2)=6*vtens*sqrt(beta%j_tot*(beta%j_tot+1.0d0))/(2*beta%j_tot+1)
            Vc(i,2,1)=Vc(i,1,2)
            Vc(i,2,2)=vcen-2*(beta%j_tot+2)*vtens/(2*beta%j_tot+1)-(beta%j_tot+2)*vls
        end do
        end select
        end subroutine


ccccccc
        subroutine rmatrix()
ccccccc
        implicit none
        integer::i,j,mm,nn !sum variables
        integer::li,lj     !lc(i),lc(j)
ccccccc
        do mm=1,nr
            do i=1,beta%nchmax
                Ech(mm,mm,i)=Ec(i)-E
            end do
        end do        
ccccccc coupled potential matrix elements Vcouple_{im,jn}
        do i=1,beta%nchmax
            do j=1,beta%nchmax
                do mm=1,nr
                    Vcouple(mm,mm,i,j)=Vc(mm,i,j)
                end do
            end do
        end do
ccccccc
!T+L(B), kinetic energy and Bloch term
        do mm=1,nr 
            do nn=1,nr
                do i=1,beta%nchmax
                    if(mm==nn) then 
                        T(mm,nn,i)=hbarc**2/2/mu/rmax*phia(nn)**2
     &         *(((4*nr**2+4*nr+3d0)*xle(nn)*(1-xle(nn))-6d0*xle(nn)+1d0)/3d0/xle(nn)/(1-xle(nn))-B_i(i))
                    else
                        T(mm,nn,i)=hbarc**2/2/mu/rmax
     &         *phia(mm)*phia(nn)*(nr**2+nr+1d0+(xle(nn)+xle(mm)-2d0*xle(mm)*xle(nn))/(xle(nn)-xle(mm))**2-1d0/(1-xle(nn))-1d0/(1-xle(mm))-B_i(i))        
                    end if
                end do
            end do
        end do
!then we add the several matrix elements together to get Cmatrix (without centrifugal term)
        do mm=1,nr
            do nn=1,nr
                do i=1,beta%nchmax
                    do j=1,beta%nchmax
                            if(i==j) then 
                            Cmat(mm,nn,i,j)=Ech(mm,nn,i)+T(mm,nn,i)+Vcouple(mm,nn,i,j)
                            else
                            Cmat(mm,nn,i,j)=Vcouple(mm,nn,i,j)
                            end if                 
                    end do
                end do
            end do
        end do     
!centrifugal term added to Cmatrix, now Cmatrix is complete
        do i=1,beta%nchmax
            do mm=1,nr
                Cmat(mm,mm,i,i)=Cmat(mm,mm,i,i)+hbarc**2/2d0/mu*lc(i)*(lc(i)+1d0)/xle(mm)**2/rmax**2
            end do
        end do
ccccccc
!reconstruct the Cmatrix for inversion, the size of C is (nr*beta%nchmax,nr*beta%nchmax)
ccccccc
        do mm=1,nr
            do nn=1,nr
                do i=1,beta%nchmax
                    do j=1,beta%nchmax
                        C((i-1)*nr+mm,(j-1)*nr+nn)=Cmat(mm,nn,i,j)
                    end do
                end do
            end do
        end do
!get the inversion of Cmatrix, and the inversion is stored just in C.
        call mat_inv(C,nr*beta%nchmax,nr*beta%nchmax)
!Rmatrix , R_{ij}=hbar^2/(2mu a)*\sum_{mn}φ_n(a)(C^{-1})_{in,jm}φ_m(a)
        do i=1,beta%nchmax
            do j=1,beta%nchmax
                Rmat(i,j)=0d0
                do mm=1,nr
                    do nn=1,nr
                        Rmat(i,j) =Rmat(i,j)+hbarc**2/2/mu/rmax*phia(mm)*C((i-1)*nr+mm,(j-1)*nr+nn)*phia(nn)        
                    end do
                end do
            end do
        end do
ccccccc
!Zmatrix: Z_O,Z_I, Smatrix: S=(Z_O)^{-1}Z_I
ccccccc
        KFN=0
        do i=1,beta%nchmax
            do j=1,beta%nchmax
                k_i=sqrt(2d0*mu*abs(Ec(i)-E)/hbarc**2)
                k_j=sqrt(2d0*mu*abs(Ec(j)-E)/hbarc**2)
ccccccc
                li=int(lc(i),4)
                lj=int(lc(j),4)
ccccccc
                allocate(FC_i(0:li),GC_i(0:li),FCP_i(0:li),GCP_i(0:li))
                allocate(FC_j(0:lj),GC_j(0:lj),FCP_j(0:lj),GCP_j(0:lj))
ccccccc
                call COUL90(k_i*rmax,z1*z2*e2*mu/hbarc**2/k_i,0d0,li,FC_i,GC_i,FCP_i,GCP_i,KFN,IFAIL)
                call COUL90(k_j*rmax,z1*z2*e2*mu/hbarc**2/k_j,0d0,lj,FC_j,GC_j,FCP_j,GCP_j,KFN,IFAIL)
ccccccc H^{+}=G+iF, H^{-}=G-iF, and their derivatives
                hlp_i=cmplx(GC_i(li),FC_i(li),kind=8)
                hln_i=cmplx(GC_i(li),-FC_i(li),kind=8)
                dhlp_i=cmplx(GCP_i(li),FCP_i(li),kind=8)
                dhln_i=cmplx(GCP_i(li),-FCP_i(li),kind=8)
                hlp_j=cmplx(GC_j(lj),FC_j(lj),kind=8)
                hln_j=cmplx(GC_j(lj),-FC_j(lj),kind=8)
                dhlp_j=cmplx(GCP_j(lj),FCP_j(lj),kind=8)
                dhln_j=cmplx(GCP_j(lj),-FCP_j(lj),kind=8)
ccccccc
                Z_O(i,j)=(k_j*rmax)**(-0.5d0)*(hlp_i*delta(i,j)-k_j*rmax*Rmat(i,j)*dhlp_j)
                Z_I(i,j)=(k_j*rmax)**(-0.5d0)*(hln_i*delta(i,j)-k_j*rmax*Rmat(i,j)*dhln_j)
ccccccc
                deallocate(FC_i,GC_i,FCP_i,GCP_i)
                deallocate(FC_j,GC_j,FCP_j,GCP_j)
            end do
        end do
ccccccc
! Smatrix: S=(Z_O)^{-1}Z_I,first we get the inverse of Z_O
                call mat_inv(Z_O,beta%nchmax,beta%nchmax)
                Smat=matmul(Z_O,Z_I)
ccccccc
                write(*,*) 'Smatrix:'
                do i=1,beta%nchmax
                    write(*,*) Smat(i,:)
                end do
                write(*,*) 'the module of S'
                do i=1,beta%nchmax
                    write(*,*) abs(Smat(i,:))
                end do
                write(*,*) 'the phase of S'
                do i=1,beta%nchmax
                        write(*,*) atan2(aimag(Smat(i,:)),real(Smat(i,:)))/2
                end do

        end subroutine
ccccccc

        end module