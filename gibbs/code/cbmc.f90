program cbmc
    use movemc
    use measmc

    implicit none
    logical lequil
    integer*4 :: i, kg, ki, ks, kgmax, kpr, kcfg
    real*8 :: dum, du, fnrnd, t1, t2

    call input
    call init

    call cpu_time(t1)

    fnrnd = real(3*nmol + 3)   ! n*trans + n*rot + n*grow + 1*vol

    call rntprn(0)
    lequil = .false.
    kgmax = kpre/kprn
    ks = 1
10  continue

    do kg = 1, kgmax
        do kpr = 1, kprn
            do ki = 1, nmol*ks
         
                i = int(rnd(dum)*fnrnd) + 1
         
                if (i <= nmol) then                    ! translation
                    call trans(i)
                    call enesite(i, 1, nsite, 1, nsite)
                    call fouriere(i)
                    du = udif(i)
                    if (du > 0.0) then
                        if (rnd(dum) > exp(-beta*du)) cycle
                    end if
                    call accept(i, 0)
                elseif (i <= 2*nmol) then              ! rotation
                    i = i - nmol
                    call rot(i)
                    call enesite(i, 1, nsite, 1, nsite)
                    call fouriere(i)
                    du = udif(i)
                    if (du > 0.0) then
                        if (rnd(dum) > exp(-beta*du)) cycle
                    end if
                    call accept(i, 1)
                elseif (i <= 3*nmol) then              ! growing
                    i = i - 2*nmol
                    tgrow = tgrow + 1.0
                    if (rnd(dum) < grow(i)) then
                        agrow = agrow + 1.0
                        call accept(i, 2)
                    end if
                else          ! volume change
                    call volchng
                end if

            end do

            if (lequil) then
                call measure
                kcfg = kcfg + 1
                if (mod(kcfg, 2000) == 0) call cfg_dump(kcfg/2000)
            end if

        end do
        call rntprn(kg*kprn)
        if (.not.lequil) call accset

    end do

    call cpu_time(t2)
    print *, t1, t2, t2 - t1

    if (.not.lequil) then
        lequil = .true.
        kgmax = kconf/kprn
        ks = ksweep
        kcfg = 0
        goto 10
    end if

    call output

    stop
end program cbmc

    subroutine input
        use movemc
        use measmc
        implicit none
        character*40 :: str, filinp, inpcfg, filmol
 
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
        use measmc
        implicit none
        real*8 :: start

        start = 3434908123.0
        call rnd_init(start)
        call ene_init
        call move_init
        call meas_init

        return
    end subroutine init

    subroutine rntprn(km)
        use movemc
        implicit none
        integer*8 :: km

        if(ttrans > 0.0) acctrn = atrans/ttrans
        if(trot > 0.0) accrot = arot/trot
        if(tgrow > 0.0) accgr = agrow/tgrow
        if(tvo > 0.0) accvol = avo/tvo
        write(*, *) km, (sum(uij(1:nmol, 1:nmol))+ufour)/real(nmol), sum(uself(1:nmol))/real(nmol), volume, acctrn, accrot, accgr, accvol
        atrans = 0.0 ; ttrans = 0.0
        arot = 0.0 ; trot = 0.0
        agrow = 0.0 ; tgrow = 0.0
        avo = 0.0 ; tvo = 0.0

        return
    end subroutine rntprn

    subroutine output
        use measmc
        implicit none

        call cfg_output(cfgout)
        call meas_output

        return
    end subroutine output
