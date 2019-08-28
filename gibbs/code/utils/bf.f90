program alvis

    implicit none
    integer, parameter :: NMX = 512, SMX = 8
    character*20 :: filin
    integer*8 :: i, n, ns, is, nmol, itot
    real*8, dimension(NMX*SMX) :: x0, y0, z0
    real*8 :: a

    call getarg(1, filin)

    open(1, file=filin, form='unformatted')
        read(1) n, nmol, a, x0(1:n), y0(1:n), z0(1:n)
    close(1)

    ns = n/nmol
    itot = 0
    write(*, 10) nmol, ns, a
    do i = 1, nmol
        write(*, 20) i
        do is = 1, ns
            itot = itot + 1
            write(*, 30) x0(itot), y0(itot), z0(itot)
        end do
    end do

    stop
10  format(I5,1X,I5,1X,F15.8)
20  format(I5)
30  format(3(1X, F12.6))
end program alvis
