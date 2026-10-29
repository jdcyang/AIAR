SUBROUTINE IterateBellman(first)

USE Parameters
USE Globals
IMPLICIT NONE
logical, intent(in) :: first
integer :: ii, s, i
real(8) :: Vdiff, Vold(Igrid,nstates), c0(Igrid,nstates)

IF (first) THEN
    c0 = r * aa + y_base
    DO s = 1, nstates
        DO i = 1, Igrid
            V(i,s) = max(c0(i,s), 1e-8)**(1.0-ggamma) / (1.0-ggamma) / rrho
        END DO
    END DO
END IF

Vdiff = 1.0
ii = 1
DO WHILE (ii <= maxit_hjb .and. Vdiff > crit_hjb)
    Vold = V
    CALL HJBUpdate          ! overwrites V with Vnew
    Vdiff = maxval(abs(V - Vold))
    ii = ii + 1
END DO

END SUBROUTINE IterateBellman
