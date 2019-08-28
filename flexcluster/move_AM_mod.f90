module move_mod
    use global_mod
    use model_mod
    use config_mod
    use energy_mod
    implicit none
    real*8, dimension(TMX) :: p_trn, p_rot, p_int, p_exc
    real*8 :: ttrn, atrn, trot, arot, tint, aint, aexc, texc, actrn, acrot, acexc, acit 
    real*8, parameter :: dposa = 1*(1.0d0/6.0d0)*(0.07d-10)**2   !SMC-MPM steps !translation step parameter
    real*8, parameter :: drota = 5*(1.0d0/6.0d0)*(0.11d0)**2     !SMC-MPM steps !rotation step parameter

    contains

    subroutine move_AllParticlesSMC(mez, nt, lx, ti, utot, p, f, m, pd)
        integer*8,intent(in) :: nt(ntype), ti(NMX)
        real*8,intent(in) :: mez, lx 
        real*8 :: utot, p(3,SMX,NMX), f(3,NMX), m(3,NMX), pd(3,SMX,NMX)    
        !lokalni
        logical rotate, intern
        integer*8 :: nmol, c1 ,c2, t1
        real*8 :: n_e, n_p(3,SMX,NMX), n_f(3,NMX), n_m(3,NMX), n_pd(3,SMX,NMX)    
        real*8 :: g,vg4(4,NMX), g3(3), pr,pom(3), posrot(3,NMX)
        real*8 :: log_alpha_01(NMX), log_alpha_tot, alp, rm(3,3), x, y, z, lxh

        ! generating new positions and calculation of alfa01
        lxh = 0.5d0*lx
        nmol = sum(nt)
        
        call random_number(g)
        call random_number(vg4(:,1:nmol))
        rotate = .false.
        intern = .false.
        if (g < p_trn) then
            !translate all the particles
            ttrn = ttrn+1
            do c1 = 1,nmol
                t1 = ti(c1)
                !generate translation and alpha01
                !gaussian contribution
                pom(1) = (2*dposa)**0.5d0*dsqrt(-2*dlog(vg4(1,c1)))*dcos(2*M_PI*vg4(2,c1))
                pom(2) = (2*dposa)**0.5d0*dsqrt(-2*dlog(vg4(1,c1)))*dsin(2*M_PI*vg4(2,c1))
                pom(3) = (2*dposa)**0.5d0*dsqrt(-2*dlog(vg4(3,c1)))*dsin(2*M_PI*vg4(4,c1))
                
                posrot(:,c1) = pom(:)+beta*dposa*f(:,c1)
                
                !alpha_mn matrix element - can be simplified
                log_alpha_01(c1) = -sum(pom(:)**2)/(4*dposa)
                    
                !geneate new position
                do c2 = 1,ns(t1)
                    n_p(:,c2,c1) = p(:,c2,c1)+posrot(:,c1)
                end do

                do c2 = 1,3
                    if (n_p(c2,1,c1) < -lxh) n_p(c2,1:ns(t1),c1) = n_p(c2,1:ns(t1),c1)+lx
                    if (n_p(c2,1,c1) >  lxh) n_p(c2,1:ns(t1),c1) = n_p(c2,1:ns(t1),c1)-lx                
                end do                
                            
            end do
        else if (g < p_rot) !rotate all particles
            trot = trot+1
            rotate = .true.
            do c1 = 1,nmol
                t1 = ti(c1)
                if (t1 == 2) then    !*** only water

                    !generate rotation and alpha01
                    !gaussian contribution
                    pom(1) = (2*drota)**0.5d0*dsqrt(-2*dlog(vg4(1,c1)))*dcos(2*M_PI*vg4(2,c1))
                    pom(2) = (2*drota)**0.5d0*dsqrt(-2*dlog(vg4(1,c1)))*dsin(2*M_PI*vg4(2,c1))
                    pom(3) = (2*drota)**0.5d0*dsqrt(-2*dlog(vg4(3,c1)))*dsin(2*M_PI*vg4(4,c1))
 
                    posrot(:,c1) = pom(:)+beta*drota*m(:,c1)
                    
                    !alpha_mn matrix element - can be simplified
                    log_alpha_01(c1) = -sum(pom(:)**2)/(4*drota)
                                    
                    alp = sum(posrot(:,c1)**2)**0.5
                    x = posrot(1,c1)/alp
                    y = posrot(2,c1)/alp
                    z = posrot(3,c1)/alp                
                    rm(1,1) =            1 + (1 - dcos(alp))*(x*x - 1)
                    rm(1,2) = -z*dsin(alp) + (1 - dcos(alp))*(x*y)
                    rm(1,3) =  y*dsin(alp) + (1 - dcos(alp))*(x*z)
                    rm(2,1) =  z*dsin(alp) + (1 - dcos(alp))*(y*x)
                    rm(2,2) =            1 + (1 - dcos(alp))*(y*y - 1)
                    rm(2,3) = -x*dsin(alp) + (1 - dcos(alp))*(y*z)
                    rm(3,1) = -y*dsin(alp) + (1 - dcos(alp))*(z*x)
                    rm(3,2) =  x*dsin(alp) + (1 - dcos(alp))*(z*y)
                    rm(3,3) =            1 + (1 - dcos(alp))*(z*z - 1)
                    !generate new position
                    n_p(:,1,c1) = p(:,1, c1)
                    do c2 = 2,ns(t1)
                        pom(:) = p(:,c2, c1)-p(:,1, c1)
                        n_p(:,c2,c1) = p(:,1, c1)+matmul(rm(:,:),pom(:))
                    end do
                else    !*** ions
                    !spherically symetrical particle doesn't rotate
                    n_p(:,1:2,c1) = p(:,1:2,c1)
                end if
                
            end do
        else ! intramolecular change
            tint = tint+1
            intern = .true.
            do c1 = 1,nmol
                t1 = ti(c1)
                if (t1 == 2) then    !*** only water
                    !gaussian contribution
                    pom(1) = (2*dposa)**0.5d0*dsqrt(-2*dlog(vg4(1,c1)))*dcos(2*M_PI*vg4(2,c1))
                    pom(2) = (2*dposa)**0.5d0*dsqrt(-2*dlog(vg4(1,c1)))*dsin(2*M_PI*vg4(2,c1))
                    pom(3) = (2*dposa)**0.5d0*dsqrt(-2*dlog(vg4(3,c1)))*dsin(2*M_PI*vg4(4,c1))
                    
                    posrot(:,c1) = pom(:)+beta*dposa*f(:,c1)
                    
                    !alpha_mn matrix element - can be simplified
                    log_alpha_01(c1) = -sum(pom(:)**2)/(4*dposa)
                        
                    !geneate new position
                    do c2 = 1,ns(t1)
                        n_p(:,c2,c1) = p(:,c2,c1)+posrot(:,c1)
                    end do
                else    !*** ions
                    !spherically symetrical particle doesn't rotate
                    n_p(:,1:2,c1) = p(:,1:2,c1)
                end if
            end do

        end if
            
        ! analyze new configuration -> calculate: n_f,n_m,n_e,n_pd(:,:,b)
        call analyzeNewConfiguration(nt, lx, ti, n_p, n_e, n_f, n_m, n_pd)
        
        ! calculate alpha10 and alpha_tot
        log_alpha_tot = 0.0d0
        if (rotate) then
            do c1 = 1, nmol
                t1 = ti(c1)
                if (t1 == 2) then
                    !alpha_mn matrix element - can be simplified
                    pom = posrot(:,c1)+beta*drota*n_m(:,c1)
                    g = -sum(pom(:)**2 )/(4*drota)
                          !(4*dpos*M_PI)**(-3.0d0/2.0d0)* dexp(-sum(pom(:)**2)/(4*dpos))
                    log_alpha_tot = log_alpha_tot-log_alpha_01(c1)+g
                end if
            end do
        else
            do c1 = 1, nmol
                !alpha_mn matrix element - can be simplified
                pom = posrot(:,c1)+beta*dposa*n_f(:,c1)
                g = -sum(pom(:)**2 )/(4*dposa)
                      !(4*dpos*M_PI)**(-3.0d0/2.0d0)* dexp(-sum(pom(:)**2)/(4*dpos))
                log_alpha_tot = log_alpha_tot-log_alpha_01(c1)+g
            end do
        end if
        
        ! acceptance test
        pr = dexp(-beta*(n_e-utot) +log_alpha_tot) 
 
        call random_number(g)
        if (g < pr) then !MPM MC step is accepted
            p(:,:,1:nmol) = n_p(:,:,1:nmol)
            pd(:,:,1:nmol) = n_pd(:,:,1:nmol)
            utot = n_e
            f(:,1:nmol) = n_f(:,1:nmol)
            m(:,1:nmol) = n_m(:,1:nmol)
            if (rotate) then
                arot = arot + 1.0d0
            else
                atrn = atrn + 1.0d0
            end if
        end if
    
    end subroutine move_AllParticlesSMC

    function randomParticle(typ, nt, ti)
        integer*8 :: randomParticle
        integer*8 ,intent(in):: typ, nt(TMX), ti(NMX)
        !local variables
        integer*8 :: c1,c2,c3,nmol
        real*8 g
        nmol = sum(nt)
        if (nt(typ) < 1) then
            randomParticle = 0
            return
        end if
        call random_number(g)
        c2 = floor(g*nt(typ))+1
        c3 = 0
        do c1 = 1,nmol
            if (ti(c1) == typ) c3 = c3+1
            if (c3 == c2) then
                c3 = c1
                exit
            end if
        end do

        randomParticle = c3
    end function

    !returns randomly uniformly generated (trial) position including random uniform rotation for molecular type typ
    function randomPosition(typ, lx)
        real*8 :: randomPosition(3,SMX)
        integer*8,intent(in) :: typ
        real*8,intent(in) :: lx
        !local variables
        integer*8 :: c1, c2
        real*8 :: newp(3,SMX), pos(3), teta, fi, alp, g, x, y, z, rm(3,3)
        
        call random_number(pos)
        pos = (pos-0.5d0)*lx
        newp(:,1) = pos(:)
        
        call random_number(fi)
        fi = fi*2*M_PI
        call random_number(teta)
        teta = dacos(2*teta-1)
        call random_number(alp)
        alp = alp*2.0*M_PI
        x = dsin(teta)*dcos(fi)
        y = dsin(teta)*dsin(fi)
        z = dcos(teta)
        
        rm(1,1) =            1 + (1 - dcos(alp))*(x*x - 1)
        rm(1,2) = -z*dsin(alp) + (1 - dcos(alp))*(x*y)
        rm(1,3) =  y*dsin(alp) + (1 - dcos(alp))*(x*z)
        rm(2,1) =  z*dsin(alp) + (1 - dcos(alp))*(y*x)
        rm(2,2) =            1 + (1 - dcos(alp))*(y*y - 1)
        rm(2,3) = -x*dsin(alp) + (1 - dcos(alp))*(y*z)
        rm(3,1) = -y*dsin(alp) + (1 - dcos(alp))*(z*x)
        rm(3,2) =  x*dsin(alp) + (1 - dcos(alp))*(z*y)
        rm(3,3) =            1 + (1 - dcos(alp))*(z*z - 1)
            
        do c1 = 2, ns(typ)
            newp(:,c1) = pos(:) + matmul(rm(:,:), rs(:, c1, typ))
        end do
        
        randomPosition = newp
    end function
    
    subroutine move_exchangeParticle(typ, lx, nt, utot, p, ti, f, m, pd)
        real*8, intent(in):: lx
        integer*8, intent(in):: typ 
        integer*8 :: ti(NMX),nt(ntype)
        real*8 :: utot, p(3,SMX,NMX), f(3,NMX), m(3,NMX), pd(3,SMX,NMX)
        !lokalni
        integer*8 :: nmol ,n_nmol, n_nt(ntype), idel, n_ti(NMX)
        real*8 :: n_e, n_p(3,SMX,NMX), n_f(3,NMX)
        real*8 :: n_m(3,NMX), n_pd(3,SMX,NMX), g, pr, pos(3)

        texc = texc+1
        
        nmol = sum(nt)
        n_nt = nt
        n_p(:,:,1:nmol) = p(:,:,1:nmol)
        n_ti(1:nmol) = ti(1:nmol)
        
        !generate random particle exchange
        call random_number(g)
        if (g < 0.5d0) then !create new particle
            if(nmol > NMX-1)then
                write(*,*)"too many particles in the box"
                stop
            end if
            n_nmol = nmol+1
            n_nt(typ) = nt(typ)+1
            n_ti(n_nmol) = typ
        
            n_p(:,:,n_nmol) = randomPosition(typ, lx) !generate new position
        
            call analyzeNewConfiguration(n_nt, lx, n_ti, n_p, n_e, n_f, n_m, n_pd)    
        
            !pr = (beta*press*vol)/n_nt(typ)*dexp(-beta*(n_e - utot)) !calculate acceptance criterion
            pr = qkin(typ)*vol/n_nt(typ)*dexp(-beta*(n_e - utot))*dexp(beta*chp(typ)) !calculate acceptance criterion

        else !delete random particle
            if (nt(typ) < 1) return        !no particles in the box
            n_nmol = nmol-1
            n_nt(typ) = nt(typ)-1
            
            idel = randomParticle(typ, nt, ti) !choose random particle to be deleted
            !copy old positions while omitting idel-particle (move last particle to idel-index)
            n_p(:,:,idel) = p(:,:,nmol)
            n_ti(idel) = ti(nmol)
            
            call analyzeNewConfiguration(n_nt, lx, n_ti, n_p, n_e, n_f, n_m, n_pd)    
            
            !pr = nt(typ)/(beta*press*vol)*dexp(-beta*(n_e - utot)) !calculate acceptance criterion
            pr = nt(typ)/(qkin(typ)*vol)*dexp(-beta*(n_e - utot))*dexp(-beta*chp(typ)) !calculate acceptance criterion
        end if
 
        !accept/reject new configuration
        call random_number(g)
        if (g < pr) then
            aexc = aexc + 1
            nt(typ) = n_nt(typ)
            utot = n_e
            p(:,:,1:n_nmol)  = n_p(:,:,1:n_nmol)
            pd(:,:,1:n_nmol) = n_pd(:,:,1:n_nmol)
            ti(1:n_nmol)     = n_ti(1:n_nmol)
            f(:,1:n_nmol)    = n_f(:,1:n_nmol)
            m(:,1:n_nmol)    = n_m(:,1:n_nmol)
        end if
 
        return
    end subroutine move_exchangeParticle

end module move_mod
