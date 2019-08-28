module measmc
    use cfgmc

    integer*8, parameter :: LMX = 10000
    integer*8, dimension(SMX, SMX) :: ijtf
    integer*8, dimension(LMX, (SMX+1)*SMX/2) :: ihist, thehist, fihist
    integer*8 :: lmax, nijt
    real*8 :: dm, dminv
    character*40 :: filcrl

    contains

    subroutine meas_init
        implicit none
        integer*8 :: i, j, ijt

        dminv  = 1.0/dm
        lmax   = int(ah*dminv)+1
        ihist = 0
        fihist = 0
        thehist = 0
        ijt = 0
        do i = 1, nsite
            do j = i, nsite
                ijt = ijt + 1
                ijtf(i, j) = ijt
                ijtf(j, i) = ijt
            end do
        end do
        nijt = ijt

        return
    end subroutine meas_init

    subroutine measure
        implicit none
        integer*8 :: i, j, is, js, im, jm, ijt, l, ibend, idih, is1, is2, is3
        real*8 :: rnx, rny, rnz, rn, rx1, ry1, rz1, rx2, ry2, rz2, ox, oy, oz
        real*8 :: ax, ay, az, rx3, ry3, rz3, bx, by, bz, cx, fi, norm
 
        do is = 1, nsite
            do js = 1, nsite
                ijt = ijtf(is, js)
                do im = 1, nmol - 1
                    i = (im - 1)*nsite + is
                    do jm = im + 1, nmol
                        j = (jm - 1)*nsite + js
                        rnx = ah - abs(ah - abs(x0(i) - x0(j)))
                        rny = ah - abs(ah - abs(y0(i) - y0(j)))
                        rnz = ah - abs(ah - abs(z0(i) - z0(j)))
                        rn = rnx*rnx + rny*rny + rnz*rnz
                        if(rn < ah*ah) then
                            l = int(sqrt(rn)*dminv) + 1
                            ihist(l, ijt) = ihist(l, ijt) + 1
                        end if
                    end do
                end do
            end do
        end do


        do i = 1, nmol
            ! bending
            do ibend = 1, nsite-2
                is = (i-1)*nsite + ibend
                is1 = is + 1
                is2 = is + 2

                rx1 = x0(is) - x0(is1)
                ry1 = y0(is) - y0(is1)
                rz1 = z0(is) - z0(is1)
                ax = rx1*rx1 + ry1*ry1 + rz1*rz1
            
                rx2 = x0(is2) - x0(is1)
                ry2 = y0(is2) - y0(is1)
                rz2 = z0(is2) - z0(is1)
                bx = rx2*rx2 + ry2*ry2 + rz2*rz2
            
                rx3 = x0(is2) - x0(is)
                ry3 = y0(is2) - y0(is)
                rz3 = z0(is2) - z0(is)
                cx = rx3*rx3 + ry3*ry3 + rz3*rz3;
            
                fi = ax + bx - cx;
                fi = fi/(2.0*sqrt(ax*bx))
                fi = acos(fi);
            
                l = int(fi*1000.0)
                thehist(l, ibend) = thehist(l, ibend) + 1
            end do

            ! torsion
            do idih = 1, nsite-3
                is = (i-1)*nsite + idih
                is1 = is + 1
                is2 = is + 2
                is3 = is + 3

                rx1 = x0(is) - x0(is1)
                ry1 = y0(is) - y0(is1)
                rz1 = z0(is) - z0(is1)
                ox = x0(is1) - x0(is2)
                oy = y0(is1) - y0(is2)
                oz = z0(is1) - z0(is2)
                rx2 = x0(is3) - x0(is2)
                ry2 = y0(is3) - y0(is2)
                rz2 = z0(is3) - z0(is2)
             
                ax = ry1*oz - rz1*oy
                ay = rz1*ox - rx1*oz
                az = rx1*oy - ry1*ox
                norm = sqrt(ax*ax + ay*ay + az*az)
                bx = ry2*oz - rz2*oy
                by = rz2*ox - rx2*oz
                bz = rx2*oy - ry2*ox
                norm = norm*sqrt(bx*bx + by*by + bz*bz)
             
                l = int(acos((ax*bx + ay*by + az*bz)/norm)*200.0)
                fihist(l, idih) = fihist(l, idih) + 1
            end do
        end do
             
        return
    end subroutine measure

    subroutine meas_output
        implicit none
        integer*8 :: l
        real*8 :: vx, dvxi, r, fisum

        open(4, file = filcrl, status = 'unknown')

        vx = 3.0*volume/(2.0*M_PI*real(kconf)*real(nmol*nmol))
        do l = 1, lmax
            r = dm*real(l - 1)
            dvxi = vx/((r+dm)**3-r**3)
            write(4, 1) r + 0.5*dm, real(ihist(l, 1:nijt))*dvxi
        end do

        vx = 180.0/(M_PI*1000.0)
        fisum = sum(thehist(:, 1))*vx
        do l = 1, int(1000*M_PI)
            write(4, 2) vx*(real(l)+0.5), thehist(l, 1:nsite-3)/fisum
        end do

        vx = 180.0/(M_PI*200.0)
        fisum = sum(fihist(:, 1))*vx
        do l = 1, int(200*M_PI)
            write(4, 3) vx*(real(l)+0.5), fihist(l, 1:nsite-3)/fisum
        end do
        
        close(4, status = 'keep')

        return
1       format(f7.4,15(1x,f8.5))
2       format(f8.3,7(1x,e12.5))
3       format(f8.3,6(1x,e12.5))
    end subroutine meas_output

end module measmc
