module globmc
    implicit none

    real*8, parameter ::  M_PI      = 3.1415926535898          !/* pie */
    real*8, parameter ::  PH_e      = 1.6021892E-19            !/* elementary charge */
    real*8, parameter ::  PH_R      = 8.31441                  !/* universal gas constant */
    real*8, parameter ::  PH_NA     = 6.022045E23              !/* Avogadro's number */
    real*8, parameter ::  PH_eps0   = 8.85418782E-12           !/* vacuum permitivity */
    real*8, parameter ::  PH_4epi   = 8.987551784952e+9        !/* 1./(4.*PI*eps0) */
    real*8, parameter ::  PH_A2m    = 1.0E-10                  !/* Angstrom to SI: 1 A = 1.0E-10 m */

    integer*8, parameter :: NMX = 512, SMX = 8
    integer*8, parameter :: RMX = NMX*SMX
    integer*8, parameter :: KMX = 12

    integer*8 :: kpre, kconf, ksweep, kprn

end module globmc
