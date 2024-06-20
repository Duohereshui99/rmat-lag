ccccccc
      module input
        use rmatmod
        implicit none
ccccccc
        contains
ccccccc
        subroutine readinput()
ccccccc
            namelist /channelbeta/ beta
            namelist /meshs/ nr,rmax
            namelist /systems/  E,mass1,mass2,z1,z2
            namelist /potentials/ str,v0,r0,a
ccccccc
            read(5,nml=channelbeta)
            read(5,nml=meshs)
            read(5,nml=systems)
            read(5,nml=potentials)
ccccccc
            write(1,nml=channelbeta)
            write(1,nml=meshs)
            write(1,nml=systems)
            write(1,nml=potentials)
ccccccc
            z12=z1*z2
            mu=amu*(mass1*mass2/(mass1+mass2))
ccccccc
        end subroutine
ccccccc
!-----------------------------------------------------------------------
ccccccc
        subroutine get_info()
#ifdef BASE
        print *, 'Base directory: ', BASE
#endif

#ifdef VERDATE
        print *, 'Version date: ', VERDATE
#endif

#ifdef VERREV
        print *, 'Version revision: ', VERREV
#endif

#ifdef COMPDATE
        print *, 'Compilation date: ', COMPDATE
#endif
            end subroutine  
ccccccc
      end module