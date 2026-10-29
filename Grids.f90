SUBROUTINE Grids
USE Parameters
USE Globals
IMPLICIT NONE
integer :: i, s

da = (amax - amin) / real(Igrid - 1, 8)
DO i = 1, Igrid
    a(i) = amin + real(i-1, 8) * da
    DO s = 1, nstates
        aa(i, s) = a(i)
    END DO
END DO

END SUBROUTINE Grids
