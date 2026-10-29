SUBROUTINE HJBUpdate

USE Parameters
USE Globals
IMPLICIT NONE

integer, parameter :: KL = nstates, KU = nstates, LDAB = 2*KL+KU+1
real(8) :: dVf(Igrid,nstates), dVb(Igrid,nstates), c0(Igrid,nstates)
real(8) :: cf(Igrid,nstates), cb(Igrid,nstates)
real(8) :: ssf(Igrid,nstates), ssb(Igrid,nstates)
real(8) :: dV0(Igrid,nstates), dVU(Igrid,nstates), u(Igrid,nstates)
real(8) :: AB(LDAB,ntot), rhs(ntot), Vv(Igrid,nstates)
integer :: ipiv(ntot), info, i, s, sp, j
logical :: posf, posb

Vv = V
c0 = r * aa + y_base

!derivatives and boundary conditions
DO s = 1, nstates
    DO i = 1, Igrid-1
        dVf(i,s) = (Vv(i+1,s) - Vv(i,s)) / da
        dVb(i+1,s) = dVf(i,s)
    END DO
    dVf(Igrid,s) = max(c0(Igrid,s), 1e-8)**(-ggamma)
    dVb(1,s)     = max(c0(1,s), 1e-8)**(-ggamma)
END DO

!upwind branches (positive derivatives only)
DO s = 1, nstates
    DO i = 1, Igrid
        posf = dVf(i,s) > 0.0
        posb = dVb(i,s) > 0.0
        IF (posf) THEN
            cf(i,s) = dVf(i,s)**(-1.0/ggamma)
            ssf(i,s) = c0(i,s) - cf(i,s)
        ELSE
            cf(i,s) = 0.0
            ssf(i,s) = -1.0e30
        END IF
        IF (posb) THEN
            cb(i,s) = dVb(i,s)**(-1.0/ggamma)
            ssb(i,s) = c0(i,s) - cb(i,s)
        ELSE
            cb(i,s) = 0.0
            ssb(i,s) = 1.0e30
        END IF
        dV0(i,s) = max(c0(i,s), 1e-8)**(-ggamma)
    END DO
END DO

!upwind value derivative selection
DO s = 1, nstates
    DO i = 1, Igrid
        IF (ssf(i,s) > 0.0) THEN
            dVU(i,s) = dVf(i,s)
        ELSE IF (ssb(i,s) < 0.0) THEN
            dVU(i,s) = dVb(i,s)
        ELSE
            dVU(i,s) = dV0(i,s)
        END IF
    END DO
END DO
ccons = dVU**(-1.0/ggamma)
u = ccons**(1.0-ggamma) / (1.0-ggamma)
savings = c0 - ccons

!assemble B = (1/Delta + rho) I - A_op in banded storage
AB = 0.0
DO s = 1, nstates
    DO i = 1, Igrid
        j = (i-1)*nstates + s
        !drift: diag = max(ssf,0)/da - min(ssb,0)/da moved with minus sign into B
        AB(KL+KU+1, j) = (1.0/Delta_hjb + rrho) &
            + max(ssf(i,s),0.0)/da - min(ssb(i,s),0.0)/da - Q(s,s)
        !drift lower: A_op(j, j-ns) = -min(ssb,0)/da -> B entry = +min(ssb,0)/da
        IF (i > 1) AB(KL+KU+1+j-(j-nstates), j-nstates) = min(ssb(i,s),0.0)/da
        !drift upper: A_op(j, j+ns) = max(ssf,0)/da -> B entry = -max(ssf,0)/da
        IF (i < Igrid) AB(KL+KU+1+j-(j+nstates), j+nstates) = -max(ssf(i,s),0.0)/da
        !CTMC coupling: B(j, j') for state s' at same asset
        DO sp = 1, nstates
            IF (sp /= s) AB(KL+KU+1+j-(j+sp-s), j+sp-s) = -Q(s,sp)
        END DO
    END DO
END DO

!rhs = u + V/Delta
DO s = 1, nstates
    DO i = 1, Igrid
        j = (i-1)*nstates + s
        rhs(j) = u(i,s) + Vv(i,s)/Delta_hjb
    END DO
END DO

CALL dgbsv(ntot, KL, KU, 1, AB, LDAB, ipiv, rhs, ntot, info)
IF (info /= 0) THEN
    WRITE(*,*) 'dgbsv failed in HJBUpdate, info = ', info
    STOP
END IF

!unpack Vnew into V
DO s = 1, nstates
    DO i = 1, Igrid
        j = (i-1)*nstates + s
        V(i,s) = rhs(j)
    END DO
END DO

END SUBROUTINE HJBUpdate
