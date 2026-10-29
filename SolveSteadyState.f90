SUBROUTINE SolveSteadyState

USE Parameters
USE Globals
IMPLICIT NONE
logical :: first
integer :: iter
real(8) :: err_K, err_Lu, err_Ls, err_tau, e_K
real(8) :: total_benefits, total_wages, benefit_bill
logical :: did_fallback

CALL TaskSpaceInit

IF (K_init > 0.0) THEN
    K = K_init
ELSE
    K = 15.0
END IF
L_u = pi_ss(2) * lbar
L_s = pi_ss(4) * lbar
tau = mmu * (pi_ss(1) + pi_ss(3)) / (pi_ss(2) + pi_ss(4))

best_err = 1.0e30
err_total = 1.0e30
err_I = 0.0
gg_prev = 0.0
first = .true.
omega_K = omega_prices
e_K_prev = 0.0
converged = 0

DO iter = 1, maxit_prices

    P_agg = (lambda_u_eff * L_u**cces_eta + &
             lambda_s_eff * (sbtc_factor * L_s)**cces_eta)**(1.0/cces_eta)

    r = aalpha * K**(aalpha - 1.0) * P_agg**(1.0 - aalpha) - ddelta

    w_s = (1.0 - aalpha) * K**aalpha * P_agg**(1.0 - aalpha - cces_eta) * &
          lambda_s_eff * sbtc_factor**cces_eta * max(L_s, 1e-10)**(cces_eta - 1.0)
    w_u = (1.0 - aalpha) * K**aalpha * P_agg**(1.0 - aalpha - cces_eta) * &
          lambda_u_eff * max(L_u, 1e-10)**(cces_eta - 1.0)
    w_s = max(w_s, 1e-8)
    w_u = max(w_u, 1e-8)

    IF (adopt_margin == 1) CALL AdoptionMargin

    y_base(:,1) = mmu * w_u * lbar + token_rebate
    y_base(:,2) = (1.0 - tau) * w_u * lbar + token_rebate
    y_base(:,3) = mmu * w_s * lbar + token_rebate
    y_base(:,4) = (1.0 - tau) * w_s * lbar - fee_per_skilled + token_rebate

    CALL IterateBellman(first)
    CALL StationaryDistribution(first)
    first = .false.

    K_imp = SUM(gmat * aa) * da
    L_u_imp = lbar * SUM(gmat(:,2)) * da
    L_s_imp = lbar * SUM(gmat(:,4)) * da
    mass_uu = SUM(gmat(:,1)) * da
    mass_ue = SUM(gmat(:,2)) * da
    mass_su = SUM(gmat(:,3)) * da
    mass_se = SUM(gmat(:,4)) * da

    total_benefits = mmu * lbar * (w_u * mass_uu + w_s * mass_su)
    total_wages = lbar * (w_u * mass_ue + w_s * mass_se)
    benefit_bill = total_benefits - token_tax_revenue
    IF (total_wages > 1e-10) THEN
        tau_imp = benefit_bill / total_wages
    ELSE
        tau_imp = tau
    END IF
    IF (tau_imp > 0.99 .or. tau_imp < -0.99) THEN
        WRITE(*,*) 'warning: tau hit the +-0.99 bound, iteration ', iter
    END IF
    tau_imp = min(max(tau_imp, -0.99), 0.99)

    err_K  = abs(K_imp - K) / max(K, 1e-8)
    err_Lu = abs(L_u_imp - L_u) / max(L_u, 1e-8)
    err_Ls = abs(L_s_imp - L_s) / max(L_s, 1e-8)
    err_tau = abs(tau_imp - tau) / max(tau, 1e-2)
    err_total = max(err_K, max(err_Lu, max(err_Ls, max(err_tau, err_I))))

    !snapshot the lowest-residual iterate
    IF (err_total < best_err) THEN
        best_err = err_total
        bK = K_imp; bLu = L_u_imp; bLs = L_s_imp; btau = tau_imp
        bV = V; bc = ccons; bg = gmat; byb = y_base
        br = r; bwu = w_u; bws = w_s; bPagg = P_agg
        bIstar = I_star; bchi = chi
        bLam_u = lambda_u; bLam_s = lambda_s; bTu = Theta_u; bTs = Theta_s
        bTCost = TokenCost; bleak = token_leak; brebate = token_rebate
        btaxrev = token_tax_revenue; bfee = fee_per_skilled; bptok = p_token
    END IF

    IF ((Display >= 2 .and. mod(iter, 10) == 0) .or. (Display >= 3 .and. iter <= 3)) THEN
        WRITE(*,'(A,I3,A,F9.3,A,F9.5,A,ES10.2)') '  it=', iter, ' K=', K, ' r=', r, ' err=', err_total
    END IF
    IF (Display >= 3 .and. iter == 1) THEN
        WRITE(*,'(A,F9.5,A,F9.5,A,F9.5,A,F9.5,A,F9.5)') '  Pagg=', P_agg, ' Lu=', L_u, ' Ls=', L_s, &
            ' lam_u=', lambda_u_eff, ' sbtc=', sbtc_factor
        WRITE(*,'(A,F9.3,A,F9.3,A,ES10.2,A,ES10.2,A,ES10.2,A,ES10.2)') &
            '  Kimp=', K_imp, ' Lu_imp=', L_u_imp, ' eK=', err_K, ' eLu=', err_Lu, &
            ' eLs=', err_Ls, ' etau=', err_tau
        WRITE(*,'(A,F9.5,A,F9.5,A,F9.5,A,F9.5,A,ES10.2)') &
            '  tau=', tau, ' tau_imp=', tau_imp, ' muu=', mass_uu, ' msu=', mass_su, &
            ' taxrev=', token_tax_revenue
    END IF
    IF (err_total < crit_prices) THEN
        converged = 1
        EXIT
    END IF

    !adaptive damping for the asset market
    e_K = K_imp - K
    !reset the FPK warm start after large price moves: at near-reducible
    !drift operators (extreme prices), a stale density is a wrong attractor
    !(path-dependent spurious high-asset equilibria)
    IF (abs(e_K) / max(K, 1e-8) > 0.2) gg_prev = 0.0
    IF (iter > 1 .and. sign(1.0_8, e_K) /= sign(1.0_8, e_K_prev) .and. e_K_prev /= 0.0) THEN
        omega_K = max(omega_K * 0.5, 1e-3)
    ELSE
        omega_K = min(omega_K * 1.1, 0.5)
    END IF
    e_K_prev = e_K
    K = max(K + omega_K * e_K, 0.1)
    L_u = omega_prices * L_u_imp + (1.0 - omega_prices) * L_u
    L_s = omega_prices * L_s_imp + (1.0 - omega_prices) * L_s
    tau = omega_tau * tau_imp + (1.0 - omega_tau) * tau

END DO

!fallback to the best iterate if the loop ended off the fixed point
did_fallback = .false.
IF (err_total >= crit_prices .and. best_err < err_total) THEN
    did_fallback = .true.
    K_imp = bK; L_u_imp = bLu; L_s_imp = bLs; tau_imp = btau
    V = bV; ccons = bc; gmat = bg; y_base = byb
    r = br; w_u = bwu; w_s = bws; P_agg = bPagg
    I_star = bIstar; chi = bchi
    lambda_u = bLam_u; lambda_s = bLam_s; Theta_u = bTu; Theta_s = bTs
    TokenCost = bTCost; token_leak = bleak; token_rebate = brebate
    token_tax_revenue = btaxrev; fee_per_skilled = bfee; p_token = bptok
END IF

K = K_imp
L_u = L_u_imp
L_s = L_s_imp
tau = tau_imp
IF (best_err < crit_prices) converged = 1

!final aggregates and inequality measures
Y = K**aalpha * P_agg**(1.0 - aalpha)
C_agg = SUM(ccons * gmat) * da
I_inv = ddelta * K
labor_share = (w_u * L_u + w_s * L_s) / Y
unemp_rate_u = mass_uu / (mass_uu + mass_ue)
unemp_rate_s = mass_su / (mass_su + mass_se)

CALL ComputeGinis

END SUBROUTINE SolveSteadyState
