MODULE Parameters
IMPLICIT NONE
SAVE

!GRID SIZES
integer, parameter :: Igrid   = 200            ! asset grid points
integer, parameter :: nstates = 4              ! (u,unemp) (u,empl) (s,unemp) (s,empl)
integer, parameter :: ntot    = Igrid*nstates  ! stacked state dimension

!HJB ITERATION
real(8), parameter :: Delta_hjb = 1.0e4
integer, parameter :: maxit_hjb = 2000
real(8), parameter :: crit_hjb  = 1.0e-11

!FPK (KFE) ITERATION
real(8), parameter :: Delta_fpk = 1.0e3
integer, parameter :: maxit_fpk = 40000
real(8), parameter :: tol_fpk   = 1.0e-13

!OUTER (EQUILIBRIUM) LOOP
integer, parameter :: maxit_prices = 200
real(8), parameter :: crit_prices    = 1.0e-6
real(8), parameter :: omega_prices   = 0.1     ! damped update for L
real(8), parameter :: omega_adopt    = 0.1     ! damped update for adoption I*
real(8), parameter :: omega_tau      = 0.5     ! damped update for tau (fast subloop)

!ASSET GRID BOUNDS
real(8), parameter :: amin = 0.0
real(8), parameter :: amax = 200.0

END MODULE Parameters
