!
!  File name: ewtest.f90
!  Date:      2005/03/29 15:06
!  Author:    Lukas Vlcek
! 


program ewtest
    use movemc
    implicit none
    integer*8 :: ik, ia, ir, i, j
    real*8 :: dr, da, usumx

    da = 1.0
    dr = 1.0

    call input         ! read input data (parameters, configuration, ...)
    call init          ! initialization of variables (normalization, ...)

    do ir = 14, 18, 1
        rcut = real(ir)*dr
        do ik = 4, 9, 1
            kmax = ik
            do ia = 3, 11, 1
                alfa = real(ia)*da/a
                deallocate(uij, usj)
                deallocate(uself)
                call init
                usumx = (sum(uij(1:nmol, 1:nmol))+ufour + sum(uself(1:nmol)))/real(nmol)
                write(*, 1) rcut, kmax, alfa, alfa*a, usumx, sum(uij(1:nmol, 1:nmol))/nmol, ufour/nmol, sum(uself(1:nmol))/nmol
            end do
            write(*, *)
        end do
    end do
    stop
1   format(f6.2,1x,i2,1x,f6.2,1x,f6.2,7(1x,f12.5))
end program ewtest  

    subroutine input
        use movemc
        use measmc
        implicit none
        character*40 :: str
        character*20 :: filinp, inpcfg, filmol
 
        call getarg(1, filinp)
        open(1, file=filinp, status='old')

        read(1,1)str
        read(1,*)nmol, volume, temper
        read(1,1)str
        read(1,*)preq, dvol
        read(1,1)str
        read(1,*)rcut, kmax
        read(1,1)str
        read(1,*)trnreq, rotreq
        read(1,1)str
        read(1,*)kpre, kconf, ksweep
        read(1,1)str
        read(1,*)kprn
        read(1,1)str
        read(1,*)dm
        read(1,1)str
        read(1,2)filmol
        read(1,1)str
        read(1,2)inpcfg
        read(1,1)str
        read(1,2)cfgout
        read(1,1)str
        read(1,2)filcrl
        close(1, status='keep')
 
        call mol_input(filmol)
        call cfg_input(inpcfg)

        return
1       format(a40)
2       format(a20)
    end subroutine input

    subroutine init
        use movemc
        implicit none
        real*8 :: start

        start = 3434908123.0
        call rnd_init(start)
        call mol_init
        call cfg_init
        call ene_init

        return
    end subroutine init
