PROGRAM Main
!AIAR model — Fortran port (architecture follows Kaplan-Moll-Violante 2018
!replication). Experiment modes via input.nml:
!  mode = 0: single steady state (A_bar, token_tax, adopt_margin, ...)
!  mode = 1: automation sweep (A_min:A_step:A_max) -> sweep.csv
!  mode = 2: token tax x A grid, adoption margin off/on -> adoption.csv
!  mode = 3: welfare CEV sweep (reference tau=0) -> welfare.csv
USE Parameters
USE Globals
IMPLICIT NONE

integer :: mode, token_tax_npts, ios
real(8) :: A_min, A_max, A_step, tax_min, tax_step

character(len=256) :: outdir
namelist /config/ mode, A_bar, token_tax, adopt_margin, task_int_weights, K_init, Display, &
    p_token, token_price_endog, A_min, A_max, A_step, tax_min, tax_step, &
    token_tax_npts, outdir

!defaults
mode = 0
K_init = -1.0
Display = 0
A_min = 0.05; A_max = 0.50; A_step = 0.025
tax_min = 0.0; tax_step = 0.25; token_tax_npts = 9
outdir = 'output'

CALL SetParameters

!apply namelist overrides if input.nml exists
OPEN(10, FILE='input.nml', STATUS='old', IOSTAT=ios)
IF (ios == 0) THEN
    READ(10, NML=config)
    CLOSE(10)
END IF
CALL execute_command_line('mkdir -p ' // trim(outdir))

CALL Grids

SELECT CASE (mode)
CASE (0)
    CALL RunSingle
CASE (1)
    CALL RunSweep
CASE (2)
    CALL RunAdoptionGrid
CASE (3)
    CALL RunWelfare
CASE DEFAULT
    WRITE(*,*) 'unknown mode = ', mode
    STOP 1
END SELECT

CONTAINS

SUBROUTINE RunSingle
    CALL SolveSteadyState
    CALL PrintSummary('single steady state')
    CALL SavePolicyCSV(trim(outdir) // '/policy.csv')
    WRITE(*,'(A)') 'wrote ' // trim(outdir) // '/policy.csv'
END SUBROUTINE RunSingle

SUBROUTINE RunSweep
    integer :: u, kk, nA
    real(8) :: Ak
    nA = nint((A_max - A_min) / A_step) + 1
    OPEN(NEWUNIT=u, FILE=trim(outdir) // '/sweep.csv', STATUS='replace', ACTION='write')
    WRITE(u,'(A)') 'A,premium,gini_wealth,gini_income,gini_cons,unemp_u,unemp_s,' // &
                   'r,w_u,w_s,K,Y,C,I_adopt,converged'
    DO kk = 1, nA
        Ak = A_min + (kk-1) * A_step
        CALL ResetForNewSolve
        A_bar = Ak
        CALL SolveSteadyState
        WRITE(*,'(A,F6.3,A,F8.4,A,F8.4,A,I2)') 'A = ', Ak, '  premium = ', w_s/w_u, &
            '  gini_w = ', wealth_gini, '  conv = ', converged
        WRITE(u,'(*(ES14.6,:,","))') Ak, w_s/w_u, wealth_gini, income_gini, cons_gini, &
            unemp_rate_u, unemp_rate_s, r, w_u, w_s, K, Y, C_agg, I_star, real(converged,8)
    END DO
    CLOSE(u)
END SUBROUTINE RunSweep

SUBROUTINE RunAdoptionGrid
    integer :: u, kk, j, m, nA
    real(8) :: Ak, tk
    real(8), allocatable :: uu_base(:,:), tau_base(:,:)
    integer :: margins(2)
    margins = [0, 1]
    nA = nint((A_max - A_min) / A_step) + 1
    ALLOCATE(uu_base(token_tax_npts, nA), tau_base(token_tax_npts, nA))
    OPEN(NEWUNIT=u, FILE=trim(outdir) // '/adoption.csv', STATUS='replace', ACTION='write')
    WRITE(u,'(A)') 'A,token_tax,unemp_u_margin,unemp_u_base,unemp_s_margin,' // &
        'I_adopt,tau_margin,tau_base,w_u,w_s,Y,C,gini_wealth,converged'
    DO m = 1, 2
        DO j = 1, token_tax_npts
            tk = tax_min + (j-1) * tax_step
            DO kk = 1, nA
                Ak = A_min + (kk-1) * A_step
                CALL ResetForNewSolve
                A_bar = Ak
                token_tax = tk
                adopt_margin = margins(m)
                CALL SolveSteadyState
                IF (margins(m) == 1) THEN
                    WRITE(u,'(*(ES14.6,:,","))') Ak, tk, unemp_rate_u, uu_base(j,kk), &
                        unemp_rate_s, I_star, tau, tau_base(j,kk), w_u, w_s, Y, C_agg, &
                        wealth_gini, real(converged,8)
                ELSE
                    uu_base(j,kk) = unemp_rate_u
                    tau_base(j,kk) = tau
                END IF
            END DO
            WRITE(*,'(A,I2,A,F6.2)') 'margin = ', margins(m), '  tax = ', tk
        END DO
    END DO
    CLOSE(u)
    DEALLOCATE(uu_base, tau_base)
END SUBROUTINE RunAdoptionGrid

SUBROUTINE RunWelfare
    !CEV relative to the zero-tax steady state at the same A.
    !lambda = (W1/W0)^{1/(1-gamma)} - 1, W = sum(V*g)*da (utilitarian);
    !group CEVs conditional on unskilled (states 1,2) / skilled (3,4).
    integer :: u, kk, j, m
    real(8) :: Ak, tk, pw
    real(8) :: V0(Igrid,nstates), g0(Igrid,nstates)
    real(8) :: W0, W0u, W0s, m0u, m0s, W1, W1u, W1s, m1u, m1s
    integer :: margins(2)
    margins = [0, 1]
    OPEN(NEWUNIT=u, FILE=trim(outdir) // '/welfare.csv', STATUS='replace', ACTION='write')
    WRITE(u,'(A)') 'A,token_tax,adopt_margin,CEV_total,CEV_unskilled,CEV_skilled,unemp_u,converged'
    DO kk = 1, nint((A_max - A_min) / A_step) + 1
        Ak = A_min + (kk-1) * A_step
        CALL ResetForNewSolve
        A_bar = Ak
        token_tax = 0.0      ! reference is always the zero-tax steady state
        adopt_margin = 0     ! (left over from the previous A's tax loop)
        CALL SolveSteadyState
        V0 = V; g0 = gmat
        W0 = SUM(V0 * g0) * da
        WRITE(*,'(A,F5.2,A,F10.4,A,F12.6,A,I2)') '  ref A=', Ak, ' K=', K, ' W0=', W0, ' conv=', converged
        m0u = SUM(g0(:,1:2)) * da; m0s = SUM(g0(:,3:4)) * da
        W0u = SUM(V0(:,1:2) * g0(:,1:2)) * da / m0u
        W0s = SUM(V0(:,3:4) * g0(:,3:4)) * da / m0s
        DO m = 1, 2
            DO j = 1, token_tax_npts
                tk = tax_min + (j-1) * tax_step
                CALL ResetForNewSolve
                A_bar = Ak
                token_tax = tk
                adopt_margin = margins(m)
                CALL SolveSteadyState
                W1 = SUM(V * gmat) * da
                m1u = SUM(gmat(:,1:2)) * da; m1s = SUM(gmat(:,3:4)) * da
                W1u = SUM(V(:,1:2) * gmat(:,1:2)) * da / m1u
                W1s = SUM(V(:,3:4) * gmat(:,3:4)) * da / m1s
                pw = 1.0 / (1.0 - ggamma)
                WRITE(u,'(*(ES14.6,:,","))') Ak, tk, real(margins(m),8), &
                    (W1/W0)**pw - 1.0, (W1u/W0u)**pw - 1.0, (W1s/W0s)**pw - 1.0, &
                    unemp_rate_u, real(converged,8)
            END DO
        END DO
        WRITE(*,'(A,F6.3)') 'welfare: A = ', Ak
    END DO
    CLOSE(u)
END SUBROUTINE RunWelfare

SUBROUTINE ResetForNewSolve
    !reset per-solve state (parameters persist across solves)
    gg_prev = 0.0
END SUBROUTINE ResetForNewSolve

END PROGRAM Main
