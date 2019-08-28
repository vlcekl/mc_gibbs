module cfgmc
    use molmc

    implicit none
    real*8, dimension(:), allocatable :: cx0, cy0, cz0, dd
    real*8, dimension(:), allocatable :: x0, y0, z0
    real*8, dimension(:), allocatable :: tx, ty, tz
    real*8 :: tcx, tcy, tcz, tdd
    real*8 :: volume, a, ah, lx, ly, lz, lxhi, lyhi, lzhi
    character*40 :: cfgout

    contains

    subroutine cfg_init
        implicit none
        integer*4 :: i, ismin, ismax

!        a = volume**(1.0/3.0)
        volume = a**3
        ah = 0.5*a
        lx = a
        ly = a
        lz = a
        lxhi = 1.0/(0.5*lx)
        lyhi = 1.0/(0.5*ly)
        lzhi = 1.0/(0.5*lz)

        ! cx, dx
        do i = 1, nmol
            ismin = (i - 1)*nsite + 1
            ismax = i*nsite
            cx0(i) = sum(x0(ismin:ismax))*fnsitei
            cy0(i) = sum(y0(ismin:ismax))*fnsitei
            cz0(i) = sum(z0(ismin:ismax))*fnsitei
            dd(i) = sqrt(maxval((x0(ismin:ismax)-cx0(i))**2 & 
                              + (y0(ismin:ismax)-cy0(i))**2 &
                              + (z0(ismin:ismax)-cz0(i))**2))
        end do

 
    end subroutine cfg_init

    subroutine cfg_allocate
        implicit none

        allocate(x0(n), y0(n), z0(n))
        allocate(cx0(nmol), cy0(nmol), cz0(nmol), dd(nmol))
        allocate(tx(nsite), ty(nsite), tz(nsite))

        return
    end subroutine cfg_allocate

    subroutine cfg_input(filcfg)
        implicit none
        character*20 :: filcfg
        integer*8 :: i, nt, nst, is, ii, ix
        real*8 :: at

        call cfg_allocate
 
        open(4, file = filcfg, status = 'old')
        read(4,*) nt, nst, a
        if((nt /= nmol) .or. (nst /= nsite)) then
            close(2, status = 'keep')
            print *, nmol, nt, a, at
            stop 'Wrong input cfg!'
            return
        end if
        is = 0
        do i = 1, nmol
            read(4, *) ix
            do ii = 1, nsite
                is = is + 1
                read(4, *) x0(is), y0(is), z0(is)
            end do
        end do
        close(4, status = 'keep')

        call cfg_init
 
        return
    end subroutine cfg_input

    subroutine cfg_output(filcfg)
        implicit none
        character*20 :: filcfg
        integer*8 :: i, is, ii
 
        open(4, file = filcfg, status = 'new')
        write(4,*) nmol, nsite, a
        is = 0
        do i = 1, nmol
            write(4, *) i
            do ii = 1, nsite
                is = is + 1
                write(4, *) x0(is), y0(is), z0(is)
            end do
        end do
        close(4, status = 'keep')
 
        return
    end subroutine cfg_output

    subroutine cfg_dump(kk)
        implicit none
        character*20 :: fil
        integer*4 :: kk, ii

        write(fil, '(i4)') kk
        fil = adjustl(fil)
        ii = len_trim(fil)
        fil(ii+1:ii+4) = '.cfg'
        open(4, file = fil, form = 'unformatted', status = 'new')
            write(4) n, nmol, a, x0, y0, z0
        close(4, status = 'keep')

        return
    end subroutine cfg_dump

end module cfgmc
