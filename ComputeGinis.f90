SUBROUTINE ComputeGinis
!Wealth, income and consumption Gini coefficients from the stationary
!joint density (mirrors compute_all_ginis in the MATLAB solver).
USE Parameters
USE Globals
IMPLICIT NONE
real(8) :: keys(ntot), vals(ntot), cm(ntot), cw(ntot)
real(8) :: gmarg(Igrid), totw, area
integer :: i, s, j, n

!wealth Gini: marginal density over the (already sorted) asset grid
DO i = 1, Igrid
    gmarg(i) = SUM(gmat(i,:))
END DO
cm(1) = gmarg(1) * da
cw(1) = gmarg(1) * a(1) * da
DO i = 2, Igrid
    cm(i) = cm(i-1) + gmarg(i) * da
    cw(i) = cw(i-1) + gmarg(i) * a(i) * da
END DO
totw = cw(Igrid)
area = 0.0
IF (totw > 0.0) THEN
    DO i = 1, Igrid-1
        area = area + (cm(i+1)-cm(i)) * (cw(i+1)/totw + cw(i)/totw) / 2.0
    END DO
    wealth_gini = 1.0 - 2.0*area
ELSE
    wealth_gini = 0.0
END IF

!income Gini: sort (y, g) by y
n = 0
DO s = 1, nstates
    DO i = 1, Igrid
        n = n + 1
        keys(n) = y_base(i,s)
        vals(n) = gmat(i,s)
    END DO
END DO
CALL SortKeys(keys, vals, n)
CALL LorenzGini(keys, vals, n, cm, cw, income_gini)

!consumption Gini: sort (c, g) by c
n = 0
DO s = 1, nstates
    DO i = 1, Igrid
        n = n + 1
        keys(n) = ccons(i,s)
        vals(n) = gmat(i,s)
    END DO
END DO
CALL SortKeys(keys, vals, n)
CALL LorenzGini(keys, vals, n, cm, cw, cons_gini)

END SUBROUTINE ComputeGinis


SUBROUTINE LorenzGini(keys, vals, n, cm, cw, gini)
!Lorenz-curve Gini from sorted (key, density) pairs; keys are the values.
USE Globals, only: da
IMPLICIT NONE
integer, intent(in) :: n
real(8), intent(in) :: keys(n), vals(n)
real(8), intent(out) :: cm(n), cw(n)
real(8), intent(out) :: gini
real(8) :: tot, area
integer :: i

cm(1) = vals(1) * da
cw(1) = vals(1) * keys(1) * da
DO i = 2, n
    cm(i) = cm(i-1) + vals(i) * da
    cw(i) = cw(i-1) + vals(i) * keys(i) * da
END DO
tot = cw(n)
IF (tot > 0.0) THEN
    area = 0.0
    DO i = 1, n-1
        area = area + (cm(i+1)-cm(i)) * (cw(i+1)/tot + cw(i)/tot) / 2.0
    END DO
    gini = 1.0 - 2.0*area
ELSE
    gini = 0.0
END IF

END SUBROUTINE LorenzGini


SUBROUTINE SortKeys(keys, vals, n)
!Heapsort keys ascending, carrying vals along.
IMPLICIT NONE
integer, intent(in) :: n
real(8), intent(inout) :: keys(n), vals(n)
real(8) :: tk, tv
integer :: i, j, k

DO i = n/2, 1, -1
    k = i
    tk = keys(k); tv = vals(k)
    DO
        j = 2*k
        IF (j > n) EXIT
        IF (j < n) THEN
            IF (keys(j) < keys(j+1)) j = j + 1
        END IF
        IF (tk >= keys(j)) EXIT
        keys(k) = keys(j); vals(k) = vals(j)
        k = j
    END DO
    keys(k) = tk; vals(k) = tv
END DO
DO i = n, 2, -1
    tk = keys(1); tv = vals(1)
    keys(1) = keys(i); vals(1) = vals(i)
    keys(i) = tk; vals(i) = tv
    k = 1
    tk = keys(1); tv = vals(1)
    DO
        j = 2*k
        IF (j > i-1) EXIT
        IF (j < i-1) THEN
            IF (keys(j) < keys(j+1)) j = j + 1
        END IF
        IF (tk >= keys(j)) EXIT
        keys(k) = keys(j); vals(k) = vals(j)
        k = j
    END DO
    keys(k) = tk; vals(k) = tv
END DO

END SUBROUTINE SortKeys
