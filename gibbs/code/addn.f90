!
!  File name: addn.f90
!  Date:      2005/07/17 21:59
!  Author:    Lukas Vlcek
! 


program  addn
    implicit none
    character*20 :: fil
    integer*4, parameter :: NMX = 512, SMX = 8
    integer*4 :: n, nmol, ns, i
    real*8 :: a, vol, x0(NMX*SMX), y0(NMX*SMX), z0(NMX*SMX)

    ns = 4
    nmol = 512 
    vol = 51464.0
    n = nmol*ns
    a = (vol)**(1.0/3.0)

    call getarg(1, fil)

    open(1, file = fil, form='unformatted')
        read(1) x0(1:n), y0(1:n), z0(1:n)
    close(1)

    open(2, file = fil, form='unformatted')
        write(2) n, nmol, a, x0(1:n), y0(1:n), z0(1:n)
    close(2)

    stop
end program addn

!  end of addn.f90 
