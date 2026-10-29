SUBROUTINE StationaryDistribution(first)

USE Parameters
USE Globals
IMPLICIT NONE

logical, intent(in) :: first
integer, parameter :: KL = nstates, KU = nstates, LDAB = 2*KL+KU+1
real(8) :: AB(LDAB,ntot), gk(ntot), gnew(ntot)
real(8) :: dVf(Igrid,nstates), dVb(Igrid,nstates), c0(Igrid,nstates)
real(8) :: ssf(Igrid,nstates), ssb(Igrid,nstates)
integer :: ipiv(ntot), info, i, s, sp, j, kk
real(8) :: gsum, d
logical :: posf, posb

!rebuild upwind drift branches from current V (same as HJBUpdate)
c0 = r * aa + y_base
DO s = 1, nstates
    DO i = 1, Igrid-1
        dVf(i,s) = (V(i+1,s) - V(i,s)) / da
        dVb(i+1,s) = dVf(i,s)
    END DO
    dVf(Igrid,s) = max(c0(Igrid,s), 1e-8)**(-ggamma)
    dVb(1,s)     = max(c0(1,s), 1e-8)**(-ggamma)
    DO i = 1, Igrid
        posf = dVf(i,s) > 0.0
        posb = dVb(i,s) > 0.0
        IF (posf) THEN
            ssf(i,s) = c0(i,s) - dVf(i,s)**(-1.0/ggamma)
        ELSE
            ssf(i,s) = -1.0e30
        END IF
        IF (posb) THEN
            ssb(i,s) = c0(i,s) - dVb(i,s)**(-1.0/ggamma)
        ELSE
            ssb(i,s) = 1.0e30
        END IF
    END DO
END DO

!assemble Bf = (1/Delta_f) I - A_op^T in banded storage
AB = 0.0
DO s = 1, nstates
    DO i = 1, Igrid
        j = (i-1)*nstates + s
        AB(KL+KU+1, j) = (1.0/Delta_fpk) + max(ssf(i,s),0.0)/da &
                         - min(ssb(i,s),0.0)/da - Q(s,s)
        !Bf(j, j-ns) = -A_op(j-ns, j) = -max(ssf(i-1,s),0)/da
        IF (i > 1) AB(KL+KU+1+j-(j-nstates), j-nstates) = -max(ssf(i-1,s),0.0)/da
        !Bf(j, j+ns) = -A_op(j+ns, j) = min(ssb(i+1,s),0)/da
        IF (i < Igrid) AB(KL+KU+1+j-(j+nstates), j+nstates) = min(ssb(i+1,s),0.0)/da
        !CTMC: Bf(j, j+sp-s) = -Q(sp, s)
        DO sp = 1, nstates
            IF (sp /= s) AB(KL+KU+1+j-(j+sp-s), j+sp-s) = -Q(sp,s)
        END DO
    END DO
END DO

!initial density: chain mass spread flat over assets, or warm start
IF (first .or. SUM(abs(gg_prev)) < 1e-14) THEN
    DO s = 1, nstates
        DO i = 1, Igrid
            j = (i-1)*nstates + s
            gk(j) = pi_ss(s) / (Igrid * da)
        END DO
    END DO
ELSE
    gk = gg_prev
END IF
gsum = SUM(gk) * da
gk = gk / gsum

!factor once (dgbsv would overwrite AB on the first call)
CALL dgbtrf(ntot, ntot, KL, KU, AB, LDAB, ipiv, info)
IF (info /= 0) THEN
    WRITE(*,*) 'dgbtrf failed in StationaryDistribution, info = ', info
    STOP
END IF

d = 1.0
kk = 1
DO WHILE (kk <= maxit_fpk .and. d > tol_fpk)
    gnew = gk / Delta_fpk
    CALL dgbtrs('N', ntot, KL, KU, 1, AB, LDAB, ipiv, gnew, ntot, info)
    IF (info /= 0) THEN
        WRITE(*,*) 'dgbtrs failed in StationaryDistribution, info = ', info
        STOP
    END IF
    gnew = max(gnew, 0.0)
    gsum = SUM(gnew) * da
    gnew = gnew / gsum
    d = maxval(abs(gnew - gk))
    gk = gnew
    kk = kk + 1
END DO

gg_prev = gk
DO s = 1, nstates
    DO i = 1, Igrid
        j = (i-1)*nstates + s
        gmat(i,s) = gk(j)
    END DO
END DO

END SUBROUTINE StationaryDistribution
