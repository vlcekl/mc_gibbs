module measure_mod
    use config_mod
    use energy_mod
    use move_mod

    logical :: lbar
    integer*8, dimension(:,:,:,:,:), allocatable :: ihist
    integer*8 :: lmax, icount
    integer*8 :: ttry, nins, nrem ! number of insertion and removal trials
    real*8 :: dm, dminv, ah, muid
    character*40 :: filcrl, filbar, frmti, frmto
    real*8, dimension(:), allocatable :: eins, erem

    contains

    subroutine measure_init
        implicit none

        dm = 0.02e-10
        dminv  = 1.0/dm
        ah = lx/2.0
        lmax   = int(ah*dminv)+1
        allocate(ihist(lmax, TMX, TMX, SMX, SMX))
        ihist = 0
        icount = 0
        if (lbar) then
            allocate(eins(nins), erem(nmax(ttry)))
            eins = 0.0 ; erem = 0.0
            open(7, file=filbar, status='unknown')
            muid = -log(qkin(ntype)*vol/float(nmax(ttry)+1))/beta*PH_NA/1000.0
            !muid = log(p1bar*beta/qkin)/beta*PH_NA/1000.0
            write(7, *) muid, muwat, muid - muwat ! ideal gas contriubtion to dF(n,n-1)
!            close(4, status='keep')
!            frmti = "(A1,1X,I3,1X,6f12.4)"
!            frmto = "(A1,1X,I3,1X,6f12.4)"
        end if

        return
    end subroutine measure_init

    function eachParticle(typ, nt, ti, inow)
        integer*8 :: eachParticle, inow
        integer*8 ,intent(in):: typ, nt(TMX), ti(NMX)
        !local variables
        integer*8 :: c1,c2,c3,nmol
        real*8 g
        nmol = sum(nt)
        if (nt(typ) < 1) then
            eachParticle = 0
            return
        end if
        c3 = 0
        do c1 = 1,nmol
            if (ti(c1) == typ) c3 = c3+1
            if (c3 == inow) then
                exit
            end if
        end do

        eachParticle = c1
    end function

    subroutine measure_do
        implicit none
        integer*8 :: i, j, l, is, js, iit, jjt, nmol
        real*8 :: rnx, rny, rnz, rn
        integer*8 :: n_nmol, n_nt(ntype), idel, n_ti(NMX), nrmv
        real*8 :: n_e, n_p(3,SMX,NMX), n_f(3,NMX), n_m(3,NMX), n_pd(3,SMX,NMX)
        real*8 :: fac
        fac = PH_NA/1000.0
 
        icount = icount + 1
        nmol = sum(nt(1:ntype))
        do i = 1, nmol-1
            iit = ti(i)
            do is = 1, ns(iit)
                do j = i + 1, nmol
                    jjt = ti(j)
                    do js = 1, ns(jjt)
                        rnx = ah - abs(ah - abs(p(1, is, i) - p(1, js, j)))
                        rny = ah - abs(ah - abs(p(2, is, i) - p(2, js, j)))
                        rnz = ah - abs(ah - abs(p(3, is, i) - p(3, js, j)))
                        rn = rnx*rnx + rny*rny + rnz*rnz
                        if(rn < ah*ah) then
                            l = int(sqrt(rn)*dminv) + 1
                            !print *, i, is, j, js, l, sqrt(rn)*1e10
                            ihist(l, iit, jjt, is, js) = ihist(l, iit, jjt, is, js) + 1
                        end if
                    end do
                end do
            end do
        end do

        if (lbar) then
            n_nt = nt
            n_p(:,:,1:nmol) = p(:,:,1:nmol)
            n_ti(1:nmol) = ti(1:nmol)

            n_nmol = nmol+1
            n_nt(ttry) = nt(ttry)+1
            n_ti(n_nmol) = ttry
            do i = 1, nins
                n_p(:,:,n_nmol) = randomPosition(ttry, lx) !generate new position
                call analyzeNewConfiguration(n_nt, lx, n_ti, n_p, n_e, n_f, n_m, n_pd)    
                eins(i) = n_e - utot
                !print *, 'ins', n_e, utot
            end do

            nrmv = min(nrem, nt(ttry))
            n_nmol = nmol-1
            n_nt(ttry) = nt(ttry)-1
            do i = 1, nrmv
                if (nrmv == nt(ttry)) then
                    idel = eachParticle(ttry, nt, ti, i) !choose random particle to be deleted
                else
                    idel = randomParticle(ttry, nt, ti) !choose random particle to be deleted
                end if
                !copy old positions while omitting idel-particle (move last particle to idel-index)
                n_p(:,:,1:nmol) = p(:,:,1:nmol)
                n_p(:,:,idel) = p(:,:,nmol)
                n_ti(1:nmol) = ti(1:nmol)
                n_ti(idel) = ti(nmol)
                call analyzeNewConfiguration(n_nt, lx, n_ti, n_p, n_e, n_f, n_m, n_pd)    
                erem(i) = n_e - utot
                !print *,'rem', n_e, utot, idel
            end do

!            open(4, file=filbar, status='old', position='append')
            write(7, *) eins(1:nins)*fac
            write(7, *) erem(1:nrmv)*fac
!            write(4, frmti) 'i', nins, eins(1:nins)*fac
!            write(4, frmto) 'r', nrmv, erem(1:nrmv)*fac
!            close(4, status='keep')

        end if
             
        return
    end subroutine measure_do

    subroutine measure_output
        implicit none
        integer*8 :: l, is, js, iit, jjt
        real*8 :: vx, dvxi, r, fhist

        open(4, file = filcrl, status = 'unknown')
        write(4, *) '#', nt(1:ntype)
        write(4, *) '#', ns(1:ntype)

        vx = 3.0*vol/(2.0*M_PI*real(icount))
        do l = 1, lmax
            r = dm*real(l - 1)
            dvxi = vx/((r+dm)**3-r**3)
            write(4, 1, advance='no') (r + 0.5*dm)*1d10
            do iit = 1, ntype
                do is = 1, ns(iit)
                    do jjt = iit, ntype
                        do js = 1, ns(jjt)
                            if (iit == 1 .and. is == 2) cycle
                            if (jjt == 1 .and. js == 2) cycle
                            if (iit == 2 .and. is == 2) cycle
                            if (jjt == 2 .and. js == 2) cycle
                            if (is == 5 .or. js == 5) cycle
                            if (iit == jjt) cycle
                            write(4, 2, advance='no') real(ihist(l, iit, jjt, is, js))*dvxi/real(nt(iit)*nt(jjt))
                            !write(4, 2, advance='no') real(ihist(l, iit, jjt, is, js))/real(nt(iit)*nt(jjt))
                        end do
                    end do
                end do
            end do
            write(4, *) ''
        end do
        
        close(4, status = 'keep')

        deallocate(ihist)
        if (lbar) then
            close(7, status='keep')
            deallocate(eins, erem)
        end if

        return
1       format(f7.4)
2       format(1x,f12.4)
    end subroutine measure_output

end module measure_mod
