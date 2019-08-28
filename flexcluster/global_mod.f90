module global_mod
    implicit none
    real*8, parameter ::  M_PI      = 3.1415926535897433832795d0  !/* pie */
!    real*8, parameter ::  M_PI      = 3.1415926535898          !/* pie */
    real*8, parameter ::  PH_e      = 1.6021892E-19            !/* elementary charge */
    real*8, parameter ::  PH_R      = 8.31441                  !/* universal gas constant */
    real*8, parameter ::  PH_NA     = 6.022045E23              !/* Avogadro's number */
    real*8, parameter ::  PH_eps0   = 8.85418782E-12           !/* vacuum permitivity */
    real*8, parameter ::  PH_4epi   = 8.987551784952e+9        !/* 1./(4.*PI*eps0) */
    real*8, parameter ::  PH_A2m    = 1.0E-10                  !/* Angstrom to SI: 1 A = 1.0E-10 m */
    real*8, parameter ::  PH_kcal   = 4184                     !/* 1000 cal in J */
    real*8, parameter ::  PH_debye  = 3.335640952d-30
    real*8, parameter ::  PH_kb     = 1.380662d-23             !/* Boltzmann constant */
    real*8, parameter ::  PH_Mu     = 1.660538782d-27          !/* unit mass */
    real*8, parameter ::  PH_esu    = 3.335640951981521d-010
    real*8, parameter ::  PH_stdP   = 1.0d5                    !/* standard pressure 1 bar = 1e5 Pa */
    real*8, parameter ::  PH_h      = 6.62606957d-34                       !/* Planck constant */

    integer*8, parameter :: TMX = 2
    integer*8, parameter :: NMX = 4096, SMX = 5
    integer*8, parameter :: RMX = NMX*SMX
    integer*8, parameter :: KMX = 12

    integer*8 :: kpre, kconf, ksweep, kprn

end module global_mod
