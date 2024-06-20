ccccccc
        module potential
            implicit none
            contains                                
            real*8 function gausspot(r,v0,r0,a)
            real*8::r,v0,r0,a
            if (a.gt.1e-6) then
            gausspot=V0*exp(-(r-r0)**2/a**2)
            else
            write(*,*)'a too small in gausspot!'
            stop
            end if
            return
            end function gausspot
        end module