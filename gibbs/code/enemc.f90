module enemc
    use cfgmc

    implicit none
    integer*4, parameter :: QSMX = 3
    integer*4, parameter :: QMX = NMX*QSMX
    integer*4, parameter :: KTMX = (2*KMX + 1)*(2*KMX + 1)*(KMX + 1)*4/10 ! max.# of useful k-vectors

    integer*4 :: kmax, kcntmax, nqsit
    real*8 :: rcut, rcutsq, alfa, facq, pipix, pipiy, pipiz, fac, alf4i
    complex*16 :: cpipix, cpipiy, cpipiz

    real*8, dimension(QSMX) :: qi
    integer*4, dimension(QSMX) :: sqit

    complex*16, dimension(QMX, 0:KMX) :: el
    complex*16, dimension(QMX) :: elm
    complex*16, dimension(QMX, -KMX:KMX) :: em, en
    complex*16, dimension(NMX, KTMX) :: elmn
    complex*16, dimension(KTMX) :: Ak, Qsum

    complex*16, dimension(QSMX, 0:KMX) :: tel
    complex*16, dimension(QSMX, -KMX:KMX) :: tem, ten
    complex*16, dimension(QSMX) :: telm
    complex*16, dimension(KTMX) :: telmn


    real*8, dimension(:, :), pointer :: uij
!    real*8, dimension(:, :), allocatable :: uij
    real*8, dimension(:, :), allocatable :: usj
    real*8, dimension(:), allocatable :: uself
    real*8 :: tuself, dufour, ufour

    real*8, parameter :: a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741
    real*8, parameter :: a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911

    contains

    subroutine ene_allocate
        implicit none

        allocate(uij(nmol, nmol), usj(nsite, nmol))
        allocate(uself(nmol))
        uij = 0.0 ; usj = 0.0 ; uself = 0.0
!
!        allocate(qi(qsmax), sqit(qsmax))
!        allocate()
!
!        ktmax = (2*kmax + 1)**2*(kmax + 1)*4/10
!        allocate(Qsum(ktmax))
!
        return
    end subroutine ene_allocate

    subroutine ene_self_tot
        implicit none
        integer*4 :: i, ix, iy, is, js
        real*8 :: rij, t, erfv

        uself = 0.0
        do i = 1, nmol
            do is = 1, nsite
                uself(i) = uself(i) - facq*alfa*qs(is)**2/sqrt(M_PI)
                do js = is+1, nsite
                    ix = (i-1)*nsite + is
                    iy = (i-1)*nsite + js
                    rij = sqrt((x0(ix) - x0(iy))**2 + (y0(ix) - y0(iy))**2 + (z0(ix) - z0(iy))**2)
                    t = 1.0/(1.0 + p*alfa*rij)
                    erfv = 1.0 - ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp(-(alfa*rij)**2)
                    uself(i) = uself(i) - facq*qs(is)*qs(js)*erfv/rij
                end do
            end do
        end do

        return
    end subroutine ene_self_tot

    subroutine ene_init
        implicit none
        integer*4 :: i, is, iq, kl, km, kn, kmmin, knmin, kcount, ix, iqmax, iy, js, iimin, iimax
        real*8 :: rl, rm, rn, rksq, rij, erfv, t

        call ene_allocate

        alfa = 7/a
        rcutsq = rcut*rcut
        facq = PH_4epi*PH_e*PH_e*PH_Na/(1000.0*PH_A2m)

        call ene_self_tot
        ! correction for self-term - all molecules

        call ene_tot

        ! fouriere part
        pipix = 2.0*M_PI/lx
        pipiy = 2.0*M_PI/ly
        pipiz = 2.0*M_PI/lz
        cpipix = (0.0, 1.0)*pipix
        cpipiy = (0.0, 1.0)*pipiy
        cpipiz = (0.0, 1.0)*pipiz
        fac = 2.0*PH_e*PH_e*PH_Na/(PH_eps0*PH_A2m*1000.0*lx*ly*lz)
        alf4i = -1.0/(4.0*alfa**2)

        ! id number of all charges in the system for later use in Ewald summtion
        el(:, 0) = 1.0
        em(:, 0) = 1.0
        en(:, 0) = 1.0

        iq = 0
        sqit = 0
        do is = 1, nsite
            if (abs(qs(is)) > 1.e-5) then
                iq = iq + 1
                sqit(iq) = is
                qi(iq) = qs(is)
            end if
        end do
        nqsit = iq

        iq = 0
        ix = 0
        do i = 1, nmol
            do is = 1, nsite
                ix = ix + 1
                if (abs(qs(is)) > 1.e-5) then
                    iq = iq + 1
                    if (iq > QMX) stop "Too many charges!"
                    el(iq, 1) = exp(cpipix*x0((ix)))
                    em(iq, 1) = exp(cpipiy*y0((ix))) ; em(iq, -1) = conjg(em(iq, 1))
                    en(iq, 1) = exp(cpipiz*z0((ix))) ; en(iq, -1) = conjg(en(iq, 1))
                end if
            end do
        end do
        iqmax = iq
        do kl = 2, kmax
            forall (iq = 1:iqmax)
                el(iq, kl) = el(iq, kl-1)*el(iq, 1)
                em(iq, kl) = em(iq, kl-1)*em(iq, 1) ; em(iq, -kl) = conjg(em(iq, kl))
                en(iq, kl) = en(iq, kl-1)*en(iq, 1) ; en(iq, -kl) = conjg(en(iq, kl))
            end forall
        end do

        ufour = 0.0
        Qsum = (0.0, 0.0)
        elmn = (0.0, 0.0)
        kmmin = 0
        knmin = 1
        kcount = 0
        do kl = 0, kmax
            rl = pipix*real(kl)
            do km = kmmin, kmax
                rm = pipiy*real(km)
                forall(iq = 1:iqmax) elm(iq) = em(iq, km)*el(iq, kl)*qi(mod(iq-1, 3)+1)
                do kn = knmin, kmax
                    if(kl*kl + km*km + kn*kn <= kmax*kmax) then
                        kcount = kcount + 1
                        rn = pipiz*real(kn)
                        rksq = rl*rl + rm*rm + rn*rn
                        Ak(kcount) = fac*exp(rksq*alf4i)/rksq
                        do i = 1, nmol
                            iimin = (i-1)*nqsit + 1
                            iimax = i*nqsit
                            elmn(i, kcount) = sum(elm(iimin:iimax)*en(iimin:iimax, kn))
                        end do
                        Qsum(kcount) = sum(elmn(1:nmol, kcount))
                        ufour = ufour + Qsum(kcount)*conjg(Qsum(kcount))*Ak(kcount)
                    end if
                end do
                knmin = -kmax
            end do
            kmmin = -kmax
        end do
        kcntmax = kcount
        if (kcntmax > KTMX) stop 'too large kcntmax'

        call four_tot
        do i = 1, nmol
            ix = (i-1)*nsite
            tx(1:nsite) = x0(ix+1:ix+nsite)
            ty(1:nsite) = y0(ix+1:ix+nsite)
            tz(1:nsite) = z0(ix+1:ix+nsite)
            call fouriere(i)
            if (abs(ufour - dufour) > 1.e-10) print *, 'b', i, ufour, dufour
        end do

        return
    end subroutine ene_init

    subroutine four_tot 
        implicit none
        integer*4 :: i, is, ii, iq, kl, km, kn, kmmin, knmin, kcount, ix, iqmax, iimin, iimax
        real*8 :: rl, rm, rn, rksq

        iq = 0
        do i = 1, nmol
            ii = (i-1)*nsite
            do ix = 1, nqsit
                iq = iq + 1
                is = ii + sqit(ix)
                el(iq, 1) = exp(cpipix*x0(is))
                em(iq, 1) = exp(cpipiy*y0(is)) ; em(iq, -1) = conjg(em(iq, 1))
                en(iq, 1) = exp(cpipiz*z0(is)) ; en(iq, -1) = conjg(en(iq, 1))
            end do
        end do
        iqmax = iq
        do kl = 2, kmax
            forall (iq = 1:iqmax)
                el(iq, kl) = el(iq, kl-1)*el(iq, 1)
                em(iq, kl) = em(iq, kl-1)*em(iq, 1) ; em(iq, -kl) = conjg(em(iq, kl))
                en(iq, kl) = en(iq, kl-1)*en(iq, 1) ; en(iq, -kl) = conjg(en(iq, kl))
            end forall
        end do

        ufour = 0.0
        Qsum = (0.0, 0.0)
        elmn = (0.0, 0.0)
        kmmin = 0
        knmin = 1
        kcount = 0
        do kl = 0, kmax
            rl = pipix*real(kl)
            do km = kmmin, kmax
                rm = pipiy*real(km)
                forall(iq = 1:iqmax) elm(iq) = em(iq, km)*el(iq, kl)*qi(mod(iq-1, 3)+1)
                do kn = knmin, kmax
                    if(kl*kl + km*km + kn*kn <= kmax*kmax) then
                        kcount = kcount + 1
                        rn = pipiz*real(kn)
                        rksq = rl*rl + rm*rm + rn*rn
                        Ak(kcount) = fac*exp(rksq*alf4i)/rksq
                        do i = 1, nmol
                            iimin = (i-1)*nqsit + 1
                            iimax = i*nqsit
                            elmn(i, kcount) = sum(elm(iimin:iimax)*en(iimin:iimax, kn))
                        end do
                        Qsum(kcount) = sum(elmn(1:nmol, kcount))
                        ufour = ufour + Qsum(kcount)*conjg(Qsum(kcount))*Ak(kcount)
                    end if
                end do
                knmin = -kmax
            end do
            kmmin = -kmax
        end do
        kcntmax = kcount

        return
    end subroutine four_tot 

    subroutine ene_tot
        implicit none
        integer*4 :: i, j, is, js, ii, jj
        real*8 :: rij, alfar, t, erfx, fqq, rnx, rny, rnz, rn, r6, rnn, qx

        uij = 0.0
        do i = 1, nmol
            do j = i+1, nmol
                rnx = ah - abs(ah - abs(cx0(i) - cx0(j)))
                rny = ah - abs(ah - abs(cy0(i) - cy0(j)))
                rnz = ah - abs(ah - abs(cz0(i) - cz0(j)))
                rn = sqrt(rnx*rnx + rny*rny + rnz*rnz)
                if (rn - dd(i) - dd(j) > rcut) cycle
                do is = 1, nsite
                    ii = (i - 1)*nsite + is
                    qx = facq*qs(is)
                    do js = 1, nsite
                        jj = (j - 1)*nsite + js
                        rnx = ah - abs(ah - abs(x0(ii) - x0(jj)))
                        rny = ah - abs(ah - abs(y0(ii) - y0(jj)))
                        rnz = ah - abs(ah - abs(z0(ii) - z0(jj)))
                        rn = rnx*rnx + rny*rny + rnz*rnz
                        if (rn < rcutsq) then
                            r6 = sig6(is, js)/(rn*rn*rn)
                            uij(i, j) = uij(i, j) + eps4(is, js)*r6*(r6 - 1.0)
                            fqq = qx*qs(js)
                            if (fqq /= 0.0) then
                                rij = sqrt(rn)
                                alfar = alfa*rij
                                t = 1.0/(1.0 + p*alfar)
                                erfx = ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp(-alfar*alfar)
                                uij(i, j) = uij(i, j) + erfx*fqq/rij
                            end if
                        end if
                    end do
                end do
            end do
        end do

        do i = 1, nmol
            ii = (i - 1)*nsite
            do is = 1, nsite-4
                do js = is + 4, nsite
                    rnx = x0(ii + is) - x0(ii + js)
                    rny = y0(ii + is) - y0(ii + js)
                    rnz = z0(ii + is) - z0(ii + js)
                    rn = rnx*rnx + rny*rny + rnz*rnz
                    r6 = sig6(is, js)/(rn*rn*rn)
                    uij(i, i) = uij(i, i) + eps4(is, js)*r6*(r6 - 1.0)
                end do
            end do
        end do

        return
    end subroutine ene_tot

    subroutine enesite(i, is1, is2, js1, js2)
        implicit none
        integer*4 :: i, j, jm, is1, is2, is, js, op, js1, js2
        real*8 :: rij, alfar, t, erfx, fqq, rnx, rny, rnz, rn, r6, rnn, qx

        usj(is1:is2, 1:nmol) = 0.0
        do is = is1, is2
            qx = facq*qs(is)
            do jm = 1, nmol
                if (jm == i) cycle
                rnx = ah - abs(ah - abs(tx(is) - cx0(jm)))
                rny = ah - abs(ah - abs(ty(is) - cy0(jm)))
                rnz = ah - abs(ah - abs(tz(is) - cz0(jm)))
                rn = sqrt(rnx*rnx + rny*rny + rnz*rnz)
                if (rn - dd(jm) > rcut) cycle
                do js = 1, nsite
                    j = (jm - 1)*nsite + js
                    rnx = ah - abs(ah - abs(tx(is) - x0(j)))
                    rny = ah - abs(ah - abs(ty(is) - y0(j)))
                    rnz = ah - abs(ah - abs(tz(is) - z0(j)))
                    rn = rnx*rnx + rny*rny + rnz*rnz
                    if (rn < rcutsq) then
                        r6 = sig6(is, js)/(rn*rn*rn)
                        usj(is, jm) = usj(is, jm) + eps4(is, js)*r6*(r6 - 1.0)
                        fqq = qx*qs(js)
                        if (fqq /= 0.0) then
                            rij = sqrt(rn)
                            alfar = alfa*rij
                            t = 1.0/(1.0 + p*alfar)
                            erfx = ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp(-alfar*alfar)
                            usj(is, jm) = usj(is, jm) + erfx*fqq/rij
                        end if
                    end if
                end do
            end do
        end do

        op = 1
        if (js2 < js1) op = -1
        do is = is1, is2
            do js = js1, js2, op
                if (abs(is - js) < 4 .or. js*op > is*op) cycle
                rnx = tx(is) - tx(js)
                rny = ty(is) - ty(js)
                rnz = tz(is) - tz(js)
                rn = rnx*rnx + rny*rny + rnz*rnz
                r6 = sig6(is, js)/(rn*rn*rn)
                usj(is, i) = usj(is, i) + eps4(is, js)*r6*(r6 - 1.0)
            end do
        end do

        return
    end subroutine enesite

    subroutine eneself
        implicit none
        integer*4 :: is, js
        real*8 :: rij, t, erfv

        tuself = 0.0
        do is = 1, nqsit
            tuself = tuself - facq*alfa*qi(is)**2/sqrt(M_PI)
            do js = is+1, nqsit
                rij = sqrt((tx(sqit(is)) - tx(sqit(js)))**2 + (ty(sqit(is)) - ty(sqit(js)))**2 + (tz(sqit(is)) - tz(sqit(js)))**2)
                t = 1.0/(1.0 + p*alfa*rij)
                erfv = 1.0 - ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp(-(alfa*rij)**2)
                tuself = tuself - facq*qi(is)*qi(js)*erfv/rij
            end do
        end do
        
        return
    end subroutine eneself

    ! fouriere part
    subroutine fouriere(i)
        implicit none
        integer*4 :: i, iq, is, kl, km, kn, kmmin, knmin, kcount
        complex*16 :: delmn

!        forall (kcount = 1:kcntmax) telmn(kcount) = elmn(i, kcount)

        tel(:, 0) = 1.0
        tem(:, 0) = 1.0
        ten(:, 0) = 1.0

        do iq = 1, nqsit
            is = sqit(iq)
            tel(iq, 1) = exp(cpipix*tx(is))
            tem(iq, 1) = exp(cpipiy*ty(is)) ; tem(iq, -1) = conjg(tem(iq, 1))
            ten(iq, 1) = exp(cpipiz*tz(is)) ; ten(iq, -1) = conjg(ten(iq, 1))
        end do

        do kl = 2, kmax
            forall (iq = 1:nqsit)
                tel(iq, kl) = tel(iq, kl-1)*tel(iq, 1)
                tem(iq, kl) = tem(iq, kl-1)*tem(iq, 1) ; tem(iq, -kl) = conjg(tem(iq, kl))
                ten(iq, kl) = ten(iq, kl-1)*ten(iq, 1) ; ten(iq, -kl) = conjg(ten(iq, kl))
            end forall
        end do
  
        dufour = 0.0
        kmmin = 0
        knmin = 1
        kcount = 0
        do kl = 0, kmax
            tel(1:nqsit, kl) = qi(1:nqsit)*tel(1:nqsit, kl)
            do km = kmmin, kmax
                telm(1:nqsit) = tel(1:nqsit, kl)*tem(1:nqsit, km)
                do kn = knmin, kmax
                    if(kl*kl + km*km + kn*kn <= kmax*kmax) then
                        kcount = kcount + 1
                        telmn(kcount) = sum(telm(1:nqsit)*ten(1:nqsit, kn))
                        delmn = Qsum(kcount) + telmn(kcount) - elmn(i, kcount)
                        dufour = dufour + Ak(kcount)*delmn*conjg(delmn)
                    end if
                end do
                knmin = -kmax
            end do
            kmmin = -kmax
        end do

        return
    end subroutine fouriere

    function udif(i)
        implicit none
        integer*4 :: i
        real*8 :: udif

        udif = sum(usj) - sum(uij(1:i-1, i)) - sum(uij(i, i:nmol))
        udif = udif + dufour - ufour

        return
    end function udif

end module enemc
