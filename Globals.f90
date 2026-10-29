MODULE Globals
USE Parameters
IMPLICIT NONE
SAVE

!PREFERENCES
real(8) :: ggamma, rrho

!PRODUCTION
real(8) :: ddelta, aalpha, cces_eta, N_task

!TASK SPACE (Acemoglu-Restrepo)
real(8) :: A0, I_A0, M0, Theta0_u, Theta0_s, m_task, nu_chi, iota_task
real(8) :: lambda_base, lambda_u_coeff, lambda_s_coeff, sbtc_zeta
integer :: task_int_weights

!MARKOV TRANSITIONS
real(8) :: lambda_eu0_u, lambda_ue0_u, lambda_eu0_s, lambda_ue0_s
real(8) :: lambda_us0, lambda_su0
real(8) :: theta_auto, theta_rein, theta_stab, phi_find

!LABOR / UI
real(8) :: lbar, mmu

!TOKEN MECHANISM
real(8) :: p_token, token_per_task, token_markup_share, token_tax, token_tax_prog
integer :: token_price_endog
real(8) :: p_token0, Q_bar, eps_supply

!ADOPTION MARGIN (Race 2018 cost condition)
integer :: adopt_margin
real(8) :: rho_mach, g_gamma, margin_m0, w_u_cal, r_cal

!EXPERIMENT STATE
real(8) :: A_bar

!GRIDS
real(8) :: a(Igrid), da, aa(Igrid,nstates)

!EQUILIBRIUM OBJECTS
real(8) :: V(Igrid,nstates), ccons(Igrid,nstates), gmat(Igrid,nstates)
real(8) :: y_base(Igrid,nstates), savings(Igrid,nstates)
real(8) :: Q(nstates,nstates), pi_ss(nstates)
real(8) :: I_A, I_star, M_bnd, Theta_u, Theta_s, chi
real(8) :: lambda_u, lambda_s, lambda_u_eff, lambda_s_eff, sbtc_factor
real(8) :: TokenCost, token_leak, token_rebate, token_tax_eff, token_tax_revenue
real(8) :: fee_per_skilled
real(8) :: K, L_u, L_s, tau, r, w_u, w_s, P_agg, Y, C_agg, I_inv, labor_share
real(8) :: K_imp, L_u_imp, L_s_imp, tau_imp
real(8) :: mass_uu, mass_ue, mass_su, mass_se
real(8) :: unemp_rate_u, unemp_rate_s
real(8) :: wealth_gini, income_gini, cons_gini
real(8) :: best_err, err_total, err_I, omega_K, e_K_prev, K_init
integer :: converged, Display

!BEST-ITERATE SNAPSHOT (fallback when the outer loop ends off the fixed point)
real(8) :: bV(Igrid,nstates), bc(Igrid,nstates), bg(Igrid,nstates)
real(8) :: byb(Igrid,nstates)
real(8) :: bK, bLu, bLs, btau, br, bwu, bws, bPagg, bIstar, bchi
real(8) :: bLam_u, bLam_s, bTu, bTs, bTCost, bleak, brebate, btaxrev, bfee, bptok

!FPK WARM START
real(8) :: gg_prev(ntot)

END MODULE Globals
