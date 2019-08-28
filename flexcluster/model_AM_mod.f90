module model_mod
    use global_mod
    implicit none
    character*2 :: fftyp
    integer*8 :: ntype
    integer*8, dimension(TMX) :: ns             ! number of particles of each type, number of interaction sites in each type
    integer*8, dimension(TMX) :: isfix_f, isfix_l, isdru_f, isdru_l ! first and last sites of fixed and drude sites
    real*8, dimension(TMX, TMX) :: sig, eps, rcin2  !LJ parameters, internal cutoff
    real*8, dimension(SMX, TMX) :: qs, ks      !partial charges, their widths, and spring force constant
    real*8, dimension(3, SMX, TMX) :: rs            !poloha sajtu vzhledem k pocatku molekuly (prvni sajt), prvni sajt je ve vode kyslik, stred otaceni a zaroven jadro NERA nebo LJ nebo EXP6
    real*8, dimension(TMX) :: mt, inx, iny, inz     ! mass and chemical potential of each type
    real*8, dimension(SMX, TMX) :: ms, si, ep, aa, bb, cc

    contains

    subroutine model_input(filmol)
        implicit none
        character*40 :: str, filmol
        integer*8 :: nmolt, tt, typ, iaux, is, iis, t1, t2, a1, a2
        real*8 :: faux
        character*8, dimension(TMX) :: molname
        integer*8, dimension(TMX) :: nst
        character*1, dimension(SMX, TMX) :: corshl
        character*8, dimension(SMX, TMX) :: sitname

        ns = 0        !site count
        rs = 0.0      !site positions
        qs = 0.0      !charge
        ks = 0.0      !spring constant

        open(1, file=filmol, status='old')
            read(1,1) str
            read(1,*) nmolt, fftyp
            if (fftyp .ne. 'LJ') then
                print *, 'Wrong force field type'
                stop
            end if
            do tt = 1, nmolt
                read(1,1) str
                read(1,*) typ, molname(typ)
                read(1,1) str
                read(1,*) ns(typ)
                isfix_f(typ) = 1
                isfix_l(typ) = 0
                do is = 1, ns(typ)
                    read(1,*) sitname(is, typ), corshl(is, typ)
                    read(1,*) ms(is, typ), qs(is, typ)
                    if (corshl(is, typ) == 'c') then
                        read(1,*) si(is, typ), ep(is, typ)
                        isfix_l(typ) = isfix_l(typ) + 1
                    else
                        read(1,*) ks(is, typ)    ! drude force constant
                    end if
                end do
                sig(typ, typ) = si(1, typ)
                eps(typ, typ) = ep(1, typ)
                mt(typ) = sum(ms(1:ns(typ), typ))
                isdru_f(typ) = isfix_l(typ) + 1
                isdru_l(typ) = ns(typ)
                read(1,1) str ! geometry
                do is = 1, ns(typ)
                    read(1,*) rs(1, is, typ), rs(2, is, typ), rs(3, is, typ)
                end do
            end do
            read(1,1) str ! internal cutoffs
            read(1,*) iaux
            do tt = 1, iaux
                read(1,*) typ, is, faux
                rcin2(typ, is) = faux**2
                rcin2(is, typ) = rcin2(typ, is)
            end do
        close(1, status='keep')

        !LB cross rules
        do t1 = 1, ntype  
            do t2 = 1, ntype  
                eps(t1, t2) = (eps(t1,t1)*eps(t2, t2))**0.5d0
                sig(t1, t2) = (sig(t1,t1)+sig(t2, t2))*0.5d0        
            end do
        end do

        eps   = eps*PH_kcal/PH_NA   !Parmeters in kcal/mol
        sig   = sig*1.0d-10
        ks    = ks*PH_kcal/PH_NA/1.0d-20
        rcin2 = rcin2*1.0d-20
        rs = rs*1.0d-10
        qs = qs*PH_e
        ms = ms*PH_Mu
        mt = mt*PH_Mu

!        print *, 'ns', ns
!        print *, 'mt', mt
!        print *, 'eps', eps
!        print *, 'sig', sig
!        print *, 'qs', qs
!        print *, 'ks', ks
!        print *, 'rcin2', rcin2
!        print *, 'rs', rs
        
        write(*,*)'# polarizability ion  =  ', PH_4epi*qs(2,1)**2/ks(2,1)
        write(*,*)'# polarizability H2O  =  ', PH_4epi*qs(5,2)**2/ks(5,2)
    
        return
1       format(a40)
    end subroutine model_input

    subroutine model_moments_of_inertia
        implicit none
        integer*8 :: it
        real*8, dimension(SMX) :: rx, ry, rz

        do it = 1, ntype
            if (it == 2) then  ! only for water
                rx(1:ns(it)) = rs(1,1:ns(it),it) - sum(rs(1,1:ns(it),it)*ms(1:ns(it),it))/mt(it)
                ry(1:ns(it)) = rs(2,1:ns(it),it) - sum(rs(2,1:ns(it),it)*ms(1:ns(it),it))/mt(it)
                rz(1:ns(it)) = rs(3,1:ns(it),it) - sum(rs(3,1:ns(it),it)*ms(1:ns(it),it))/mt(it)
                inx(it) = sum(ms(1:ns(it),it)*(ry(1:ns(it))**2 + rz(1:ns(it))**2))
                iny(it) = sum(ms(1:ns(it),it)*(rz(1:ns(it))**2 + rx(1:ns(it))**2))
                inz(it) = sum(ms(1:ns(it),it)*(rx(1:ns(it))**2 + ry(1:ns(it))**2))
            end if
        end do

        return
    end subroutine model_moments_of_inertia

end module model_mod

