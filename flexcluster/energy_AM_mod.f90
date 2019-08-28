module energy_mod
    use global_mod
    use model_mod
    use config_mod
    implicit none
    real*8 :: utot, rcut, rcutsq, tit, ait

    contains

    subroutine energy_init
        implicit none

        rcut = rcut*1d-10
        rcutsq = rcut**2

        return
    end subroutine energy_init

    ! analyzeNewConfiguration -> calculates new properties of a previously generated new configutation: n_f,n_m,n_e
    subroutine analyzeNewConfiguration(nt,lx, ti, p, utot,f,m, pd )
    !subroutine analyzeNewConfiguration(nt,lx, ti, p, utot,f,m, pd )
        integer*8,intent(in):: nt(ntype), ti(NMX)        !analyzed simulation box index (1 or 2), number of particles in the box
        real*8,intent(in):: lx, p(3,SMX,NMX)        !analyzed simulation box index (1 or 2)
        real*8,intent(out):: utot, f(3,NMX), m(3,NMX), pd(3,SMX,NMX)    !output properties (energy, forces, torques, drude site positions)
        !lokalni promenne
        integer*8 :: nmol, c1, c2, cit
        integer*8 :: a1, a2, t1, t2, max_c1,max_a1
        real*8 maxdrudepos2,drudepos2,drudepos(3),pomv2(3),pomv(3),pomm(3),pomf(3),rm1(3),rm2(3)
        real*8 pair_e,f_r(3),pair_f1(3),pair_f2(3),pair_m1(3),pair_m2(3)    !pair properties
        real*8 :: rsq,rv(3),r,k,pomk,rsqp,rvc(3),rc, pom_pd(3,SMX,NMX), e_r, lxh, aux
        
        lxh = 0.5*lx
        nmol = sum(nt)

        ! test hardcore cuttoffs -> fast infinite energy rejection
        e_r = 0
        do c1 = 1, nmol-1
            do c2 = c1+1, nmol
                !interanal cutoff
                rv = p(:,1,c2)-p(:,1,c1)
                rsq = sum((lxh-dabs(lxh-dabs(rv(:))))**2)
                if (rsq < rcin2(ti(c1),ti(c2))) then
                    e_r = e_r+1.0d30
                end if
            end do
        end do
        
        if (e_r > 1.0d20) then
            utot = 1.0d30
            return
        end if
        
        ! delete:    f, m  and set initial drude positions to the origin of the molecules: pd
        f(:,1:nmol) = 0
        m(:,1:nmol) = 0
        pd(:,:,1:nmol) = p(:,:,1:nmol) !nastaveni vsech drude poloh na fixed polohy, dale se pocita jen s drude polohami
         
        ! induced dipole iteration -> add contributions to: n_e, n_f, n_m
        ccit:do cit = 1,30
            maxdrudepos2 = 0
           
            !calculate new force acting on drude particles and new drude positions
            do c1 = 1, nmol
                t1 = ti(c1)
                do a1 = isdru_f(t1), isdru_l(t1)
                    !calculate real space contribution to forces acting on drude particles
                    !toto lze asi vytahnout z cyklu pres c1 a prepsat na cyklus pres pary + intramolekularni prispevek
                    f_r(:) = 0.0d0
                    do c2 = 1, nmol
                        t2 = ti(c2)
                        pomk = PH_4epi*qs(a1,t1)
                        
                        rvc = pd(:,1,c1)-pd(:,1,c2)
                        where(rvc>lxh)rvc = rvc-lx
                        where(rvc<-lxh)rvc = rvc+lx                
                        if (sum(rvc**2) > rcutsq) cycle
                                            
                        do a2 = isfix_f(t2),isdru_l(t2)    !fixed sites must be before drude sites in the model definition
                            rv = pd(:,a1,c1)-pd(:,a2,c2)
                            where(rv>lxh)rv = rv-lx
                            where(rv<-lxh)rv = rv+lx                
                            rsq = sum(rv**2)
                            !if(rsq<rcutsq.and.rsq>1.0d-30)then
                            if (rsq > 1.0d-30) then
                                r = dsqrt(rsq)
                                if (c1 /= c2) then    
                                    !drude with another particle
                                    f_r = f_r + pomk*qs(a2,t2)*rv/(r*rsq)
                                end if
                            end if
                        end do
                    end do
                                           
                    !calculate new drude positions and maximum drude square displacement
                    pomv(:) = f_r(:)/ks(a1,t1)
                    
                    drudepos(:) = p(:,a1,c1)+pomv(:)-pd(:,a1,c1)
                    pd(:,a1,c1) = p(:,a1,c1)+pomv(:)    
                    drudepos2 = sum(drudepos*drudepos)
                    if (drudepos2 > maxdrudepos2) maxdrudepos2 = drudepos2
                end do
            end do        
            
            !end iteration criterion
            if (maxdrudepos2 < 1.0d-26) exit ccit    
        end do ccit
        
        if (cit == 31) then
            write(*,*)'diverged',maxdrudepos2
            utot = 1.0d30
            return
        end if
        
        ait = ait + cit - 1
        tit = tit + 1
        
        ! analyze optimized configuration -> calculate: n_e,n_f,n_m
        ! real space pair contributions
        do c1 = 1, nmol-1
            t1 = ti(c1)
            do c2 = c1+1, nmol
                t2 = ti(c2)
                
                pair_e = 0
                pair_f1 = 0
                pair_m1 = 0
                pair_m2 = 0
        
                !first site - first site exp6 contribution
                rv = p(:,1,c1)-p(:,1,c2)
                where (rv >  lxh) rv = rv - lx
                where (rv < -lxh) rv = rv + lx
                rsq = sum(rv**2)
                if (rsq > rcutsq) cycle

                !non-coul
                pair_e = pair_e + 4*eps(t1,t2)*(((sig(t1,t2)**2)/rsq)**6-((sig(t1,t2)**2)/rsq)**3)
                pomf = rv/rsq*4*eps(t1,t2)*(12*((sig(t1,t2)**2)/rsq)**6 - 6*((sig(t1,t2)**2)/rsq)**3)
                pair_f1 = pair_f1 + pomf
                !zadny moment sil, protoze je exp6 ve stredu otaceni
            
                !charge-charge both fixed and drude p is copied to pd for fixed sites
                do a1 = isfix_f(t1),isdru_l(t1)
                    rm1 = pd(:,a1,c1)-pd(:,1,c1)
                    do a2 = isfix_f(t2),isdru_l(t2)
                        rm2 = pd(:,a2,c2)-pd(:,1,c2)
                        rv = pd(:,a1,c1)-pd(:,a2,c2)
                        where (rv >  lxh) rv = rv - lx
                        where (rv < -lxh) rv = rv + lx
                        rsq = sum(rv**2)
                        r = dsqrt(rsq)
                        
                        pair_e =  pair_e + PH_4epi*qs(a1,t1)*qs(a2,t2)/r !PC energy
                        pomf = PH_4epi*qs(a1,t1)*qs(a2,t2)* rv/(r*rsq)  !PC force
                        
                        pair_f1 = pair_f1+pomf
                        !GRF torque
                        pomm(1) = rm1(2)*pomf(3) - rm1(3)*pomf(2)
                        pomm(2) = rm1(3)*pomf(1) - rm1(1)*pomf(3)
                        pomm(3) = rm1(1)*pomf(2) - rm1(2)*pomf(1)
                        pair_m1 = pair_m1+pomm
                        pomm(1) = rm2(3)*pomf(2) - rm2(2)*pomf(3)
                        pomm(2) = rm2(1)*pomf(3) - rm2(3)*pomf(1)
                        pomm(3) = rm2(2)*pomf(1) - rm2(1)*pomf(2)
                        pair_m2 = pair_m2+pomm
                    end do
                end do        
                
                e_r = e_r+pair_e
                f(:,c1) = f(:,c1)+pair_f1(:)
                f(:,c2) = f(:,c2)-pair_f1(:)
                m(:,c1) = m(:,c1)+pair_m1(:)
                m(:,c2) = m(:,c2)+pair_m2(:)
            end do
        end do
        
        !real space self contributions -> spring energies, GRF self energies, forces, torques
        do c1 = 1,nmol
            t1 = ti(c1)
            pair_e = 0
                
            !drude spring energy
            do a1 = isdru_f(t1),isdru_l(t1)
                rv = pd(:,a1,c1)-p(:,a1,c1)
                rsq = sum(rv**2)
                pair_e = pair_e+0.5d0*ks(a1,t1)*rsq
            end do
                
            !linear charge scaling by the lambda coupling parameter
            e_r = e_r+pair_e
        end do
        utot  =  e_r
    end subroutine 

end module energy_mod
