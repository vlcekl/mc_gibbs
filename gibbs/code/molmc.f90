module molmc
    use globmc

    implicit none
    integer*8 :: n, nmol
    integer*8 :: nsite, nbond, nang, ndih
    real*8 :: fnsitei
    real*8 :: temper, beta, preq
    real*8, dimension(SMX) :: qs, si, ep
    real*8, dimension(SMX, SMX) :: sig6, eps4
    real*8, dimension(SMX, SMX) :: lbond, ang1, ang2, dih1, dih2, dih3, dih4

    contains


    subroutine mol_input(filmol)
        implicit none
        character*20 :: filmol
        character*80 :: str
        integer*8 :: is, js, ks, ls, l
    
        open(4, file = filmol, status = 'old')
            read(4,*) str, nsite
            do is = 1, nsite
                read(4, *) str, qs(is), si(is), ep(is)
            end do
        
            read(4,*) str, nbond
            do l = 1, nbond
                read(4, *) is, js, lbond(is, js)
                lbond(js, is) = lbond(is, js)
            end do
 
            read(4,*) str, nang
            do l = 1, nang
                read(4, *) is, js, ks
                read(4, *) ang1(is, ks), ang2(is, ks)
                ang1(ks, is) = ang1(is, ks)
                ang2(ks, is) = ang2(is, ks)
            end do
 
            read(4,*) str, ndih
            do l = 1, ndih
                read(4, *) is, js, ks, ls
                read(4, *) dih1(is, ls), dih2(is, ls), dih3(is, ls), dih4(is, ls)
                dih1(ls, is) = dih1(is, ls)
                dih2(ls, is) = dih2(is, ls)
                dih3(ls, is) = dih3(is, ls)
                dih4(ls, is) = dih4(is, ls)
            end do
        close(4)

        call mol_init

        return
    end subroutine mol_input

    subroutine mol_init
        implicit none
        integer*8 :: is, js
        real*8 :: k2kj

        n = nmol*nsite
        fnsitei = 1.0/real(nsite)
        beta = 1000.0/(temper*PH_R)               ! 1/RT (kJ/mol)^(-1)
        preq = preq*PH_Na*1.e-30*1.e-3*1.e5 ! p[bar]->PV[kJ/mol]

        k2kj = PH_R/1000.0
        ep = ep*k2kj
        ang1 = ang1*k2kj/2.0
        ang2 = ang2*M_PI/180.0
        dih1 = dih1*k2kj
        dih2 = dih2*k2kj
        dih3 = dih3*k2kj
        dih4 = dih4*k2kj
        forall (is = 1:nsite, js = 1:nsite)
            sig6(is, js) = (0.5*(si(is) + si(js)))**6
            eps4(is, js) = 4.0*sqrt(ep(is)*ep(js))
        end forall

        return
    end subroutine mol_init

end module molmc

