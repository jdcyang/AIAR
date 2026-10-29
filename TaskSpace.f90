SUBROUTINE TaskSpaceInit
!Task-space quantities at the technological frontier I_A = A_bar*N,
!transition matrix Q and its stationary distribution pi_ss.
USE Parameters
USE Globals
IMPLICIT NONE

!technological frontier and skill boundary
I_A = I_A0 + iota_task * (A_bar - A0) * N_task
M_bnd = M0 - m_task * (A_bar - A0) * N_task
I_A = max(min(I_A, M_bnd - 0.05), 0.05)
M_bnd = max(min(M_bnd, N_task - 0.05), I_A + 0.05)

I_star = I_A

sbtc_factor = (A_bar / A0)**sbtc_zeta

!progressive token tax schedule (constant over the solve) — must precede
!DownstreamObjects, which consumes token_tax_eff
token_tax_eff = max(token_tax + token_tax_prog * (A_bar - A0), 0.0)

CALL DownstreamObjects

END SUBROUTINE TaskSpaceInit


SUBROUTINE AdoptionMargin
!Race (2018) cost condition: task i operated by capital+AI iff machine unit
!cost <= labor unit cost w_u/gamma(i). Recomputes all downstream objects at
!the adopted frontier I_star (called inside the outer loop when
!adopt_margin = 1).
USE Parameters
USE Globals
IMPLICIT NONE
real(8) :: gamma_bar, cKt, I_tilde, I_new

gamma_bar = w_u_cal / (margin_m0 * rho_mach * (r_cal + ddelta))
cKt = rho_mach * (r + ddelta) + p_token + token_tax_eff
I_tilde = I_A + log(w_u / (gamma_bar * cKt)) / g_gamma
I_new = max(min(I_A, I_tilde), 0.05)
I_new = min(I_new, M_bnd - 0.05)
err_I = abs(I_new - I_star) / max(I_star, 0.01)
I_star = omega_adopt * I_new + (1.0 - omega_adopt) * I_star

CALL DownstreamObjects

END SUBROUTINE AdoptionMargin


SUBROUTINE DownstreamObjects
!All objects that depend on the adopted frontier I_star:
!task shares, chi, CES weights, transition matrix, stationary distribution,
!token accounting.
USE Parameters
USE Globals
IMPLICIT NONE
real(8) :: A_use, leu_u, lee_u, leu_s, lee_s, lus, lsu
real(8) :: Mq(nstates,nstates), rhs(nstates), piv(nstates)
integer :: ipiv(nstates), info, s

Theta_u = M_bnd - I_star
Theta_s = N_task - M_bnd
chi = (Theta_u / Theta0_u)**nu_chi
chi = min(max(chi, 0.1), 2.0)

IF (task_int_weights == 1) THEN
    lambda_u = (1.0 - lambda_base) * (Theta_u / Theta0_u)**(1.0 - cces_eta)
    lambda_s = lambda_base * (Theta_s / Theta0_s)**(1.0 - cces_eta)
ELSE
    lambda_u = (1.0 - lambda_base) * (Theta_u / Theta0_u)**lambda_u_coeff
    lambda_s = lambda_base * (Theta_s / Theta0_s)**lambda_s_coeff
END IF
lambda_s_eff = lambda_s / (lambda_u + lambda_s)
lambda_u_eff = lambda_u / (lambda_u + lambda_s)

!transition matrix at automation in use
A_use = I_star / N_task
leu_u = lambda_eu0_u * chi**phi_find
lee_u = lambda_ue0_u * (1.0 + theta_auto * (A_use - A0))
leu_s = lambda_eu0_s
lee_s = lambda_ue0_s * (1.0 - theta_stab * (A_use - A0))
lus   = lambda_us0 * (1.0 + theta_rein * (A_use - A0))
lsu   = lambda_su0
Q = 0.0_8
Q(1,1) = -(leu_u + lus); Q(1,2) = leu_u;  Q(1,3) = lus
Q(2,1) = lee_u;          Q(2,2) = -(lee_u + lus); Q(2,4) = lus
Q(3,1) = lsu;            Q(3,3) = -(leu_s + lsu); Q(3,4) = leu_s
Q(4,2) = lsu;            Q(4,3) = lee_s;  Q(4,4) = -(lee_s + lsu)

!stationary distribution: solve Q' pi = 0 with sum(pi)=1
!(replace last equation of Q' with the normalization)
Mq = TRANSPOSE(Q)
Mq(nstates,:) = 1.0
rhs = 0.0
rhs(nstates) = 1.0
piv = rhs
CALL dgesv(nstates, 1, Mq, nstates, ipiv, piv, nstates, info)
IF (info /= 0) THEN
    WRITE(*,*) 'dgesv failed in pi_ss, info = ', info
    STOP
END IF
pi_ss = max(piv, 0.0)
pi_ss = pi_ss / SUM(pi_ss)

!token accounting (only adopted tasks consume tokens)
IF (token_price_endog == 1) THEN
    p_token = p_token0 * (I_star / Q_bar)**(1.0 / eps_supply)
END IF
TokenCost  = p_token * token_per_task * I_star
token_leak = (1.0 - token_markup_share) * TokenCost
token_rebate = token_markup_share * TokenCost
token_tax_revenue = token_tax_eff * token_per_task * I_star
fee_per_skilled = (TokenCost + token_tax_revenue) / max(pi_ss(4), 1e-10)

END SUBROUTINE DownstreamObjects
