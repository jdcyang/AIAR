SUBROUTINE SavePolicyCSV(fname)
!Policy functions, savings drift, value function and density per state.
USE Parameters
USE Globals
IMPLICIT NONE
character(len=*), intent(in) :: fname
integer :: u, i

OPEN(NEWUNIT=u, FILE=fname, STATUS='replace', ACTION='write')
WRITE(u,'(A)') 'a,c_uu,c_ue,c_su,c_se,sav_uu,sav_ue,sav_su,sav_se,' // &
               'V_uu,V_ue,V_su,V_se,g_uu,g_ue,g_su,g_se'
DO i = 1, Igrid
    WRITE(u,'(*(ES14.6,:,","))') a(i), ccons(i,1), ccons(i,2), ccons(i,3), ccons(i,4), &
        savings(i,1), savings(i,2), savings(i,3), savings(i,4), &
        V(i,1), V(i,2), V(i,3), V(i,4), &
        gmat(i,1), gmat(i,2), gmat(i,3), gmat(i,4)
END DO
CLOSE(u)

END SUBROUTINE SavePolicyCSV


SUBROUTINE PrintSummary(tag)
!One-line equilibrium summary to stdout.
USE Parameters
USE Globals
IMPLICIT NONE
character(len=*), intent(in) :: tag
real(8) :: gap

gap = Y - token_leak - C_agg - I_inv
WRITE(*,'(A)') '--- ' // trim(tag) // ' ---'
WRITE(*,'(A,F9.4,A,F9.5,A,F9.4)') '  K = ', K, '  r = ', r, '  Y = ', Y
WRITE(*,'(A,F9.4,A,F9.4,A,F9.5)') '  w_u = ', w_u, '  w_s = ', w_s, '  tau = ', tau
WRITE(*,'(A,F9.4,A,F9.4,A,F9.4)') '  unemp_u = ', unemp_rate_u, '  unemp_s = ', unemp_rate_s, &
    '  I* = ', I_star
WRITE(*,'(A,F9.4,A,F9.4,A,F9.4)') '  gini_w = ', wealth_gini, '  gini_i = ', income_gini, &
    '  gini_c = ', cons_gini
WRITE(*,'(A,ES10.2,A,I2,A,ES10.2)') '  best_err = ', best_err, '  converged = ', converged, &
    '  resource_gap = ', gap

END SUBROUTINE PrintSummary
