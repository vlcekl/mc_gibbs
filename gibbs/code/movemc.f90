module movemc
    use rndmod
    use enemc

    implicit none
    ! translation + rotation
    real*8 :: trnreq, rotreq, atrans, ttrans, arot, trot, flp, dfi, acctrn, accrot
    ! volume change
    real*8 :: avo, tvo, accvol, dvol, oufour
    real*8, dimension(:), allocatable :: ouself
!    real*8, dimension(:, :), allocatable :: uijo
    real*8, dimension(:, :), pointer :: uijo, uija
    complex*16, dimension(NMX, KTMX) :: oelmn
    complex*16, dimension(KTMX) :: oQsum, oAk

    ! grow
    integer*4, parameter :: TRMX = 10
    real*8 :: agr, accgr, agrow, tgrow
    real*8, dimension(TRMX) :: eusi, tbx, tby, tbz, ux
    real*8, dimension(:, :), allocatable :: ust
    integer*4 :: ntry, isfirst, islast, istart
    
    contains

    subroutine move_allocate
        implicit none

        allocate(uijo(nmol, nmol))
        uijo = 0.0
        allocate(ust(TRMX, nmol))
        ust = 0.0
        allocate(ouself(nmol))
        ouself = 0.0

        return
    end subroutine move_allocate

    subroutine move_init
        implicit none

        call move_allocate

        if(trnreq /= 0.0) trnreq = 1.0/trnreq
        if(rotreq /= 0.0) rotreq = 1.0/rotreq
        atrans = 0.0
        ttrans = 0.0
        arot   = 0.0
        trot   = 0.0
        avo = 0.0
        tvo = 0.0
        agrow = 0.0
        tgrow = 0.0
        flp = 0.1
        dfi = 0.1
        uijo = 0.0

        return
    end subroutine move_init

    subroutine trans(i)
        implicit none
        integer*4 :: is, i, ismin
        real*8 :: dx, dy, dz, dum

        ttrans = ttrans + 1.0

        do
            dx = rnd(dum) - 0.5
            dy = rnd(dum) - 0.5
            dz = rnd(dum) - 0.5
            if ((dx*dx + dy*dy + dz*dz) < 0.25) exit
        end do
        dx = dx*flp
        dy = dy*flp
        dz = dz*flp

        ismin = (i - 1)*nsite

        forall (is = 1:nsite) tx(is) = x0(ismin + is) + dx
        tcx = sum(tx(1:nsite))*fnsitei
        if (tcx > ah) then
            tcx = tcx - lx
            tx(1:nsite) = tx(1:nsite) - lx
        elseif (tcx < -ah) then
            tcx = tcx + lx
            tx(1:nsite) = tx(1:nsite) + lx
        end if

        forall (is = 1:nsite) ty(is) = y0(ismin + is) + dy
        tcy = sum(ty(1:nsite))*fnsitei
        if (tcy > ah) then
            tcy = tcy - ly
            ty(1:nsite) = ty(1:nsite) - ly
        elseif (tcy < -ah) then
            tcy = tcy + ly
            ty(1:nsite) = ty(1:nsite) + ly
        end if

        forall (is = 1:nsite) tz(is) = z0(ismin + is) + dz
        tcz = sum(tz(1:nsite))*fnsitei
        if (tcz > ah) then
            tcz = tcz - lz
            tz(1:nsite) = tz(1:nsite) - lz
        elseif (tcz < -ah) then
            tcz = tcz + lz
            tz(1:nsite) = tz(1:nsite) + lz
        end if

        tdd = dd(i)

        return
    end subroutine trans


    subroutine rot(i)
        implicit none
        integer*4 :: is, i, ismin, ismax
        real*8 :: dx, dy, dz, px, py, pz, pr, ox, oy, oz, fi, alpha, sina, cosa, dummy
 
        trot = trot + 1.0

        fi = 2.0*M_PI*rnd(dummy)
        oz = 2.0*rnd(dummy) - 1.0
        oy = sqrt(1.0 - oz*oz)
        ox = cos(fi)*oy
        oy = sin(fi)*oy
 
        alpha = dfi*(rnd(dummy) - 0.5)
        sina = sin(alpha)
        cosa = cos(alpha) - 1.0
 
        ismin = (i - 1)*nsite + 1
        tcx = cx0(i)
        tcy = cy0(i)
        tcz = cz0(i)
            ismax = i*nsite
            tcx = sum(x0(ismin:ismax))*fnsitei
            tcy = sum(y0(ismin:ismax))*fnsitei
            tcz = sum(z0(ismin:ismax))*fnsitei
        ismin = ismin - 1
        do is = 1, nsite
            dx = x0(ismin + is) - tcx
            dy = y0(ismin + is) - tcy
            dz = z0(ismin + is) - tcz
            pr = ox*dx + oy*dy + oz*dz
            px = dx - ox*pr
            py = dy - oy*pr
            pz = dz - oz*pr
            tx(is) = tcx + dx + cosa*px + sina*(py*oz - pz*oy)
            ty(is) = tcy + dy + cosa*py + sina*(pz*ox - px*oz)
            tz(is) = tcz + dz + cosa*pz + sina*(px*oy - py*ox)
        end do

        tdd = dd(i)
 
        return
    end subroutine rot

    subroutine accept(i, mtyp)
        implicit none
        integer*4 :: i, is, j, ismin, mtyp, kc

        select case(mtyp)
            case(0)
                atrans = atrans + 1.0
                cx0(i) = tcx
                cy0(i) = tcy
                cz0(i) = tcz
            case(1)
                arot = arot + 1.0
            case(2)
                agr = agr + 1.0
                uself(i) = tuself
                cx0(i) = sum(tx(1:nsite))*fnsitei
                cy0(i) = sum(ty(1:nsite))*fnsitei
                cz0(i) = sum(tz(1:nsite))*fnsitei
                dd(i) = sqrt(maxval((tx(1:nsite)-cx0(i))**2    & 
                                  + (ty(1:nsite)-cy0(i))**2    &
                                  + (tz(1:nsite)-cz0(i))**2))
        end select

        ismin = (i - 1)*nsite
        forall (is = 1:nsite)
            x0(ismin + is) = tx(is) 
            y0(ismin + is) = ty(is) 
            z0(ismin + is) = tz(is) 
        end forall

        forall (j = 1:i-1) uij(j, i) = sum(usj(1:nsite, j))
        forall (j = i:nmol) uij(i, j) = sum(usj(1:nsite, j))

        forall (kc = 1:kcntmax)
            Qsum(kc) = Qsum(kc) + telmn(kc) - elmn(i, kc)
            elmn(i, kc) = telmn(kc)
        end forall

        ufour = dufour
 
        return
    end subroutine accept


    subroutine accset
        implicit none
 
        if(abs(acctrn - trnreq) > 0.001) then
            acctrn = acctrn*trnreq
            if(acctrn < 0.75) then
                acctrn = 0.75
            else
                if(acctrn > 1.25) acctrn = 1.25
            end if
            flp = flp*acctrn
            if(flp > ah) flp = ah
        end if
        if(abs(accrot - rotreq) > 0.001) then
            accrot = accrot*rotreq
            if(accrot < 0.75) then
                accrot = 0.75
            else
                if(accrot > 1.25) accrot = 1.25
            end if
            dfi = dfi*accrot
            if(dfi > M_PI) dfi = M_PI
        end if
 
        return
    end subroutine accset


    subroutine volchng
        implicit none
        integer*4 :: i, j, is, js
        real*8 :: aold, volold, scal, du, cond, ahold, dlj, dum, dx, dy, dz

        aold = a
        ahold = ah
        volold = volume
        uija => uij
        uij => uijo
        uijo => uija
        oufour = ufour
!        ouself = uself
        oelmn = elmn
        oQsum = Qsum
        oAk = Ak


      ! create new

        volume = volume + dvol*(rnd(dum) - 0.5)
        a = volume**(1.0/3.0)
        ah = 0.5*a
        lx = a ; lxhi = 2.0/lx ; pipix = M_PI*lxhi ; cpipix = (0.0, 1.0)*pipix
        ly = a ; lyhi = 2.0/ly ; pipiy = M_PI*lyhi ; cpipiy = (0.0, 1.0)*pipiy
        lz = a ; lzhi = 2.0/lz ; pipiz = M_PI*lzhi ; cpipiz = (0.0, 1.0)*pipiz
        fac = 2.0*PH_e*PH_e*PH_Na/(PH_eps0*PH_A2m*1000.0*lx*ly*lz)
        scal = a/aold
!        alfa = alfa/scal
!        alf4i = -1.0/(4.0*alfa**2)

        do i = 1, nmol
            dx = cx0(i)*(scal - 1.0)
            dy = cy0(i)*(scal - 1.0)
            dz = cz0(i)*(scal - 1.0)
            cx0(i) = cx0(i) + dx
            cy0(i) = cy0(i) + dy
            cz0(i) = cz0(i) + dz
            forall (is = (i-1)*nsite+1:i*nsite)
                x0(is) = x0(is) + dx
                y0(is) = y0(is) + dy
                z0(is) = z0(is) + dz
            end forall
        end do

        call ene_tot
        call four_tot
        call ene_self_tot

        
!        du = sum(uij) + ufour + sum(uself) - sum(uijo) - oufour  - sum(ouself)
        du = sum(uij) + ufour - sum(uijo) - oufour
!        print *, 'self', scal, alfa*lx
!        print *, sum(uself), sum(ouself)
!        print *, ufour, oufour
!        print *, sum(uij), sum(uijo)
!        print *, 'v', scal, du
!        print *, sum(uij)-sum(uijo), ufour - oufour

        ! LJ correction
        dlj = real(nmol)*(1.0/volume - 1.0/volold)*sum(eps4*sig6*(sig6/(6.0*rcut**6) - 1.0))/rcut**3

        du = du + 2.0*real(nmol)*M_PI*dlj/3.0

        cond = 3.0*real(nmol)*log(scal) - beta*(du + preq*(volume - volold))

        if (cond < 0.0) then
            if (exp(cond) < rnd(dum)) then   ! not accepted: restore old
                volume = volold
                a = aold
                ah = 0.5*a
                lx = a ; lxhi = 2.0/lx ; pipix = M_PI*lxhi ; cpipix = (0.0, 1.0)*pipix
                ly = a ; lyhi = 2.0/ly ; pipiy = M_PI*lyhi ; cpipiy = (0.0, 1.0)*pipiy
                lz = a ; lzhi = 2.0/lz ; pipiz = M_PI*lzhi ; cpipiz = (0.0, 1.0)*pipiz
                fac = 2.0*PH_e*PH_e*PH_Na/(PH_eps0*PH_A2m*1000.0*lx*ly*lz)
                scal = 1.0/scal
!                alfa = alfa/scal
!                alf4i = -1.0/(4.0*alfa**2)
                do i = 1, nmol
                    dx = cx0(i)*(scal - 1.0)
                    dy = cy0(i)*(scal - 1.0)
                    dz = cz0(i)*(scal - 1.0)
                    cx0(i) = cx0(i) + dx
                    cy0(i) = cy0(i) + dy
                    cz0(i) = cz0(i) + dz
                    forall (is = (i-1)*nsite+1:i*nsite)
                        x0(is) = x0(is) + dx
                        y0(is) = y0(is) + dy
                        z0(is) = z0(is) + dz
                    end forall
                end do
                uija => uij
                uij => uijo
                uijo => uija
                ufour = oufour
!                uself = ouself
                elmn = oelmn
                Qsum = oQsum
                Ak = oAk
                return
            end if
        end if

        avo = avo + 1

        return
    end subroutine volchng

    function grow(i)
        implicit none
        integer*4 :: op, i, ix, is, itry
        real*8 :: grow, ww, dum, rn, psum, rosn, roso, udi, udis, dufi, txo, tyo, tzo

        ix = (i-1)*nsite
        forall (is = 1:nsite)
            tx(is) = x0(ix + is)
            ty(is) = y0(ix + is)
            tz(is) = z0(ix + is)
        end forall

        istart = int(rnd(dum)*real(nsite)) + 1
        if (rnd(dum) < 0.5) then
            op = 1
            isfirst = 1
            islast = nsite
            call enesite(i, isfirst, istart, isfirst, istart)
        else
            op = -1
            isfirst = nsite
            islast = 1
            call enesite(i, istart, isfirst, isfirst, istart)
        end if

        udis = 0.0
        ntry = 7
        roso = 1.0
        do is = istart+op, islast, op 
            call enesite(i, is, is, isfirst, is)
            txo = tx(is)
            tyo = ty(is)
            tzo = tz(is)
            ux(ntry+1) = sum(usj(is, 1:nmol))
            select case (abs(is - op - isfirst))
                case (2:SMX)
                    call rot_bend_dih(i, is, op)
                case (1)
                    call rot_bend(i, is, op)
                case (0)
                    call rot_plain(i, is, op)
                case default
                    stop 'Error in grow()!'
            end select
            ux(1:ntry) = sum(ust(1:ntry, 1:nmol), dim=2)
            udi = minval(ux(1:ntry+1))
            roso = roso*sum(exp(-beta*(ux(1:ntry+1)-udi)))
            udis = udis - udi
            tx(is) = txo
            ty(is) = tyo
            tz(is) = tzo
        end do

        ntry = ntry + 1
        rosn = 1.0
        do is = istart+op, islast, op 
            select case (abs(is - op - isfirst))
                case (2:SMX)
                    call rot_bend_dih(i, is, op)
                case (1)
                    call rot_bend(i, is, op)
                case (0)
                    call rot_plain(i, is, op)
                case default
                    stop 'Error in grow()!'
            end select
            ux(1:ntry) = sum(ust(1:ntry,1:nmol), dim=2)
            udi = minval(ux(1:ntry))
            eusi(1:ntry) = exp(-beta*(ux(1:ntry)-udi))
            ww = sum(eusi(1:ntry))
            rosn = rosn*ww
            udis = udis + udi

            rn = rnd(dum)*ww
            psum = 0.0
            do itry = 1, ntry
                psum = psum + eusi(itry)
                if (psum > rn) then
                    tx(is) = tbx(itry)
                    ty(is) = tby(itry)
                    tz(is) = tbz(itry)
                    usj(is, 1:nmol) = ust(itry, 1:nmol)
                    exit
                end if
            end do
        end do

        ! compute the fouriere energy of the new molecule
        call fouriere(i)
        call eneself
        dufi = dufour - ufour + tuself - uself(i)

!        print *, 'i', i, istart, islast
!        print *, dufour, ufour
!        print *, tuself, uself(i)
 
        if (-beta*(udis + dufi) > log(huge(rosn))) then
            print *, 'huge'
            grow = 1.0
        else
            grow = rosn/roso*exp(-beta*(udis + dufi))
        end if

        return
    end function grow


    subroutine rot_plain(i, is, op)
        implicit none
        integer*4 :: itry, is, op, is1, i
        real*8 ::  oy, oz, lb, dum, fi

        is1 = is - op
        lb = lbond(is, is1)

        do itry = 1, ntry
            fi = 2.0*M_PI*rnd(dum)
            oz = 2.0*rnd(dum) - 1.0
            oy = sqrt(1.0 - oz*oz)*lb
            tx(is) = tx(is1) + cos(fi)*oy
            ty(is) = ty(is1) + sin(fi)*oy
            tz(is) = tz(is1) + oz*lb
            call enesite(i, is, is, isfirst, is)
            tbx(itry) = tx(is)
            tby(itry) = ty(is)
            tbz(itry) = tz(is)
            ust(itry, 1:nmol) = usj(is, 1:nmol)
        end do

        return
    end subroutine rot_plain


    subroutine rot_bend(i, is, op)
        implicit none
        integer*4 :: itry, is, op, is1, is2, i
        real*8 :: rx, ry, rz, lb, dum, pr, px, py, pz, fi, an1, an2, ox, oy, oz, ax, ay, az, cosa, sina, ur

        is1 = is - op
        is2 = is - 2*op
        lb = lbond(is, is1)
        an1 = ang1(is, is2)
        an2 = ang2(is, is2)

        ax = (tx(is2) - tx(is1))  ! unit -a
        ay = (ty(is2) - ty(is1))
        az = (tz(is2) - tz(is1))
        pr = 1.0/sqrt(ax**2 + ay**2 + az**2)
        ax = ax*pr
        ay = ay*pr
        az = az*pr
        ox = ay
        oy = -ax
        oz = 0.0
        pr = 1.0/sqrt(ox*ox + oy*oy)
        ox = ox*pr
        oy = oy*pr

!        an1 = sqrt(1.0/(an1*beta))
        do itry = 1, ntry
            do
                rz = 2.0*rnd(dum) - 1.0
                ur = an1*(acos(rz) - an2)**2
                if (rnd(dum) < exp(-beta*ur)) exit
            end do
            sina = sqrt(1.0 - rz*rz)
            cosa = rz - 1.0
!            do
!                rz = an1*rndgauss_old(dum) + an2
!                sina = sin(rz)
!                if (rnd(dum) < sina) exit 
!            end do
!            cosa = sqrt(1.0 - sina*sina) - 1.0
            pr = ox*ax + oy*ay + oz*az
            px = ax - ox*pr
            py = ay - oy*pr
            pz = az - oz*pr
            rx = lb*(ax + cosa*px + sina*(py*oz - pz*oy))    ! dih.ang = 0
            ry = lb*(ay + cosa*py + sina*(pz*ox - px*oz))
            rz = lb*(az + cosa*pz + sina*(px*oy - py*ox))


            fi = 2.0*M_PI*rnd(dum)
            sina = sin(fi)
            cosa = cos(fi) - 1.0
            pr = ax*rx + ay*ry + az*rz 
            px = rx - ax*pr
            py = ry - ay*pr
            pz = rz - az*pr
            tx(is) = tx(is1) + rx + cosa*px + sina*(py*az - pz*ay)   ! set dih.ang = fi
            ty(is) = ty(is1) + ry + cosa*py + sina*(pz*ax - px*az)
            tz(is) = tz(is1) + rz + cosa*pz + sina*(px*ay - py*ax)
            call enesite(i, is, is, isfirst, is)
            tbx(itry) = tx(is)
            tby(itry) = ty(is)
            tbz(itry) = tz(is)
            ust(itry, 1:nmol) = usj(is, 1:nmol)

        end do

        return
    end subroutine rot_bend


    subroutine rot_bend_dih(i, is, op)
        implicit none
        integer*4 :: itry, is, op, is1, is2, is3, i
        real*8 :: rx, ry, rz, lb, dum, ur, fi, an1, an2, dh1, dh2, dh3, dh4, ax, ay, az, ox, oy, oz
        real*8 :: pr, px, py, pz, sina, cosa

        is1 = is - op
        is2 = is - 2*op
        is3 = is - 3*op
        lb = lbond(is, is1)
        an1 = ang1(is, is2)
        an2 = ang2(is, is2)
        dh1 = dih1(is, is3)
        dh2 = dih2(is, is3)
        dh3 = dih3(is, is3)
        dh4 = dih4(is, is3)
        ax = (tx(is2) - tx(is1))  ! minus a
        ay = (ty(is2) - ty(is1))
        az = (tz(is2) - tz(is1))
        pr = 1.0/sqrt(ax**2 + ay**2 + az**2)
        ax = ax*pr
        ay = ay*pr
        az = az*pr
        rx = (tx(is3) - tx(is2))
        ry = (ty(is3) - ty(is2))
        rz = (tz(is3) - tz(is2))
        ox = az*ry - ay*rz           ! minus (reverted) to change -a to a
        oy = ax*rz - az*rx
        oz = ay*rx - ax*ry
        pr = 1.0/sqrt(ox**2 + oy**2 + oz**2)
        ox = ox*pr
        oy = oy*pr
        oz = oz*pr

!        an1 = sqrt(1.0/(an1*beta))
        do itry = 1, ntry
            do
                rz = 2.0*rnd(dum) - 1.0
                ur = an1*(acos(rz) - an2)**2
                if (rnd(dum) < exp(-beta*ur)) exit
            end do
            sina = sqrt(1.0 - rz*rz)
            cosa = rz - 1.0
!            do
!                rz = an1*rndgauss_old(dum) + an2
!                sina = sin(rz)
!                if (rnd(dum) < sina) exit 
!            end do
!            cosa = sqrt(1.0 - sina*sina) - 1.0
            pr = ox*ax + oy*ay + oz*az
            px = ax - ox*pr
            py = ay - oy*pr
            pz = az - oz*pr
            rx = lb*(ax + cosa*px + sina*(py*oz - pz*oy))    ! dih.ang = 0
            ry = lb*(ay + cosa*py + sina*(pz*ox - px*oz))
            rz = lb*(az + cosa*pz + sina*(px*oy - py*ox))
            do  ! must reflect position of site-3
                fi = 2.0*M_PI*rnd(dum)
                ur = dh1 + dh2*(1.0 + cos(fi)) + dh3*(1.0 - cos(2.0*fi)) + dh4*(1.0 + cos(3.0*fi))
                if (rnd(dum) < exp(-beta*ur)) exit
            end do
            sina = sin(fi)
            cosa = cos(fi) - 1.0
            pr = ax*rx + ay*ry + az*rz 
            px = rx - ax*pr
            py = ry - ay*pr
            pz = rz - az*pr
            tx(is) = tx(is1) + rx + cosa*px + sina*(py*az - pz*ay)   ! set dih.ang = fi
            ty(is) = ty(is1) + ry + cosa*py + sina*(pz*ax - px*az)
            tz(is) = tz(is1) + rz + cosa*pz + sina*(px*ay - py*ax)
            call enesite(i, is, is, isfirst, is)
            tbx(itry) = tx(is)
            tby(itry) = ty(is)
            tbz(itry) = tz(is)
            ust(itry, 1:nmol) = usj(is, 1:nmol)
        end do

        return
    end subroutine rot_bend_dih

end module movemc
