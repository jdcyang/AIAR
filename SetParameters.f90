SUBROUTINE SetParameters
!China first-layer calibration (see set_parameters_AR.m)
USE Parameters
USE Globals
IMPLICIT NONE

!preferences
ggamma = 2.0
rrho   = 0.010        ! quarterly, ~4% annualized

!production
ddelta   = 0.025      ! quarterly, ~10% annualized
aalpha   = 0.45       ! capital share (China: Bai-Hsieh-Qian 0.45-0.55)
cces_eta = 0.3        ! CES parameter (elasticity 1/(1-eta) ~ 1.43)
N_task   = 1.0

!task space
A0       = 0.20
I_A0     = A0 * N_task
M0       = 0.65
Theta0_u = M0 - I_A0
Theta0_s = N_task - M0
m_task   = 0.06
nu_chi   = 0.30
iota_task = 1.0

!CES weights
lambda_base    = 0.43
lambda_u_coeff = 1.5
lambda_s_coeff = 0.2
sbtc_zeta      = 1.0
task_int_weights = 1

!Markov transitions (quarterly, China: surveyed unemployment ~5.1%)
lambda_eu0_u = 0.55
lambda_ue0_u = 0.033
lambda_eu0_s = 0.60
lambda_ue0_s = 0.022
lambda_us0   = 0.01
lambda_su0   = lambda_us0 * 0.75 / 0.25   ! 75/25 skill split
theta_auto = 0.35
theta_rein = 0.20
theta_stab = 0.10
phi_find   = 0.80

!labor / UI
lbar = 1.0
mmu  = 0.12           ! effective UI replacement rate (China: low coverage)

!token
p_token = 0.0
token_per_task = 1.0
token_markup_share = 0.0
token_tax = 0.0
token_tax_prog = 0.0
token_price_endog = 0
p_token0 = 1.0
Q_bar = 0.25
eps_supply = 1.0

!adoption margin (Race 2018)
adopt_margin = 0
rho_mach = 15.89
g_gamma  = 1.5
margin_m0 = 1.10
w_u_cal = 1.6684      ! baseline equilibrium w_u at A=0.05 (recompute on recalibration)
r_cal   = 0.00983     ! baseline equilibrium r at A=0.05

!experiment
A_bar = A0

END SUBROUTINE SetParameters
