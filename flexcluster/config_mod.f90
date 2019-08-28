module config_mod
    use global_mod
    use model_mod
    implicit none

    real*8 :: temper, beta, vol, lx, mu_water, muwat, dGp, press  ! T, 1/kT, box size, ...
    integer*8, dimension(TMX) :: nt, nmax    ! number of particles of each type, number of interaction sites in each type
    real*8, dimension(TMX) :: chp, qkin      ! mass and chemical potential of each type
    integer*8, dimension(:), allocatable :: ti   !type of each molecule (1-na,2-cl,3-h2o)
    real*8, dimension(:,:,:), allocatable :: p, pd    !position of molecules and its interaction sites (pd=drude) in simulation box coordinates (0 in the center of the box)
    character*40 :: cfgout

    contains

    subroutine config_init
        implicit none
        real*8 :: lamb, qtrn, qrot, symm, inerts, p1bar, qk
        integer*8 :: tp

        call model_moments_of_inertia

        ! TDM
        beta = 1.d0/(PH_kb*temper)
        chp = chp*1d3/PH_NA

        if (ntype == 2) then
            tp = 2  ! molecular type of water
            !chp(2) = (-240.3008403029283d0+ 15 )*1d3/PH_NA  ! water ch.p. calculated for SPC/E bulk water, try values close to this
            mu_water = chp(tp)
            dGp = (-228.582*1d3/PH_NA - mu_water)    !-228.582 = JANAF gas at 298.15K    
   
            lamb = PH_h*sqrt(beta/(2.0*M_PI*mt(tp)))
            print *, 'lamb', lamb, log(lamb)/beta, mu_water, dGp, log(lamb**3)/beta*PH_NA/1000.0
            !qtrn =   (sqrt(3*mt(2)/beta)/PH_h)**3     ! from equipartition theorem
            symm = 2.0                                          ! water symmetry number
            inerts = (inx(tp)*iny(tp)*inz(tp))**(1./3.)             ! water moments of inertia ^(1/3)
            !qtrn = sqrt(3.0*mt(tp)/(beta*PH_h**2))**3                     ! translational p.f.
            !qrot = sqrt(8.0*M_PI**2*ine/(beta*PH_h**2))**3 * (M_PI)**0.5/symm  ! rotational p.f.
            !print *, 'qrot', qrot, log(qrot)/beta*PH_NA/1000.0
            qrot = sqrt(2.0*M_PI*inerts/(beta*PH_h**2))**3 * 8*M_PI**2/symm  ! rotational p.f.
            qtrn = sqrt(2.0*M_PI*mt(tp)/(beta*PH_h**2))**3                     ! translational p.f.
            qk = qtrn*qrot
            print *, 'qtrn', qtrn, log(qtrn)/beta*PH_NA/1000.0
            print *, 'qrot', qrot, log(qrot)/beta*PH_NA/1000.0
   
            p1bar = 1e5
            muwat = log(p1bar*beta/qk)/beta*PH_NA/1000.0 ! std chem pot of ig H2O
            print *, 'muwat, 1e5', muwat
            print *, 'Dch', chp(2)*PH_NA/1000.0 - muwat
            p1bar = 1e-1
            print *, 'muwat, 1e-1', log(p1bar*beta/qk)/beta*PH_NA/1000.0
            press = qk*exp(beta*chp(2))/beta
            print *, 'press', press
        else
            tp = 1
            qtrn = sqrt(2.0*M_PI*mt(tp)/(beta*PH_h**2))**3                     ! translational p.f.
            qrot = 1.0
            qk = qtrn*qrot
            print *, 'qtrn', qtrn, log(qtrn)/beta*PH_NA/1000.0
        end if

        qkin(ntype) = qk

        !nt = (/1_8,0_8,0_8/) !starting particle counts  - tady zmenit, kdyz se ma pouzit jiny iont
        vol = vol*1.0d-30
        lx = vol**(1./3.)   ! box size
!        ti(1) = 1                           !initial type array (with one particle)
!        p(:,:,1) = 0                        !set the first particle in the center

    end subroutine config_init

    subroutine config_allocate
        implicit none

        allocate(p(3, SMX, NMX), pd(3, SMX, NMX))
        allocate(ti(NMX))

        return
    end subroutine config_allocate

    subroutine config_input(filcfg)
        implicit none
        character*20 :: filcfg
        integer*8 :: i, ix, is, tt, iaux, tx

        call config_allocate

        open(1, file=filcfg, status='old')
            read(1, *) lx, lx, lx
            if (lx**3 /= vol) then
                write(*, *) 'Wrong box size', lx**3, vol
                stop
            end if
            read(1, *) ix  ! number of types in cfg file
            do tx = 1, ix
                read(1, *) tt, nt(tt)
                do i = 1, nt(tt)
                    ti(i) = tt
                    do is = 1, ns(tt)
                        read(1, *) iaux, iaux, p(1,is,i), p(2,is,i), p(3,is,i)
                    end do
                end do
            end do
        close(1, status='keep')

        call config_init
 
        return
    end subroutine config_input

    subroutine config_output(filcfg)
        implicit none
        character*20 :: filcfg
        integer*8 :: i, is, ii

        open(4, file=filcfg, status='unknown')
            write(4, *) lx, lx, lx
            write(4, *) ntype  ! number of types in cfg file
            do ii = 1, ntype
                write(1, *) ii, nt(ii)
                do i = 1, nt(ii)
                    do is = 1, ns(ii)
                        write(4, *) i, is, p(1,is,i), p(2,is,i), p(3,is,i)
                    end do
                end do
            end do
        close(4, status='keep')
 
        return
    end subroutine config_output

    subroutine config_dump(kk)
        implicit none
        character*20 :: fil
        integer*4 :: kk, ii

        write(fil, '(i4)') kk
        fil = adjustl(fil)
        ii = len_trim(fil)
        fil(ii+1:ii+4) = '.cfg'
        open(4, file = fil, form = 'unformatted', status = 'new')
            write(4) p, pd
        close(4, status = 'keep')

        return
    end subroutine config_dump

end module config_mod
