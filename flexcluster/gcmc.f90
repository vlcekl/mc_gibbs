program gcmc
    use global_mod
    use model_mod
    use config_mod
    use energy_mod
    use move_mod
    use measure_mod
    implicit none
    integer*8 :: ic, it, itprn, ipre, iequil, i
    real*8 :: g
    real*8, dimension(3, NMX) :: f, m        !approximate forces and torques acting on molecules used for SMC-MPM

    call random_seed
    call init

    call rntprn(0_8)

    ! main cycle
    do it = 1, ipre + iequil
        do ic = 1, itprn
            call random_number(g)
            if (g < p_exc .and. nt(ntype) < nmax(ntype)) then !particle deletion/creation
                call move_exchangeParticle(ntype, lx, nt, utot, p, ti, f, m, pd)
            else !SMC-MPM move (translation or rotation with equal probability 0.5)
                call move_AllParticlesSMC(0.5d0, nt, lx, ti, utot, p, f, m, pd)
            end if        
        end do
        call rntprn(it)
        if (it > ipre) then
            call measure_do
        end if
    end do

    call measure_output

    contains

    subroutine init
        use model_mod
        use config_mod
        use energy_mod
        use measure_mod
        implicit none
        character*40 :: str, filinp, filcfg, filfld
        integer*8 :: i, ix, is, tt, iaux, tx
 
        call getarg(1, filinp)
        open(1, file=filinp, status='old')
            read(1,1)str
            read(1,*)temper
            read(1,1)str
            read(1,*)vol
            read(1,1)str
            read(1,*)ntype
            read(1,1)str
            read(1,*)(chp(i), i=1,ntype)
            read(1,1)str
            read(1,*)(nmax(i), i=1,ntype)
            read(1,1)str
            read(1,*)rcut
            read(1,1)str
            read(1,*)(p_trn(i), i=1,ntype)
            read(1,1)str
            read(1,*)(p_rot(i), i=1,ntype)
            read(1,1)str
            read(1,*)(p_int(i), i=1,ntype)
            read(1,1)str
            read(1,*)(p_exc(i), i=1,ntype)
            read(1,1)str
            read(1,*)ipre, iequil, itprn
            read(1,1)str
            read(1,1)filfld
            read(1,1)str
            read(1,1)filcfg
            read(1,1)str
            read(1,1)filcrl
            read(1,1)str
            read(1,*)lbar
            if (lbar) then
                read(1,1)str
                read(1,*)ttry, nins, nrem
                read(1,1)filbar
            end if
        close(1, status='keep')

        ! Force field (call mol_input)
        call model_input(filfld)
        call config_input(filcfg)

        call energy_init

        call measure_init

        call analyzeNewConfiguration(nt, lx, ti, p, utot, f, m, pd)
 
        return
1       format(a40)
2       format(a20)
    end subroutine init
    
    subroutine rntprn(km)
        use energy_mod
        use move_mod
        implicit none
        integer*8 :: km
        real*8 :: cx
        cx = (-dGp + log(beta*PH_stdP*vol)/beta)
        cx = log(beta*PH_stdP*dexp(-beta*dGp))/beta
        cx = chp(2)

        if (ttrn > 0.0) actrn = atrn/ttrn
        if (trot > 0.0) acrot = arot/trot
        if (texc > 0.0) acexc = aexc/texc
        if (tit  > 0.0) acit  = ait/tit
        
        write(*,'(i8,2i5,1x,2(1x,1f9.3),4(1x,1f5.3))') km, nt, cx*PH_NA/1000.0, utot*PH_NA/1000.0, actrn, acrot, acexc, acit
        ttrn = 0 ; trot = 0 ; texc = 0 ; tit = 0
        atrn = 0 ; arot = 0 ; aexc = 0 ; ait = 0

        return
    end subroutine rntprn
 
end program gcmc
