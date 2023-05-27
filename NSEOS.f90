!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads neturon star EOS table extracted form CompOSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READEOSTABLE
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i

! Dummy variables !
REAL (DP) :: dummy

! Read the number of lines in the file !
nlines = 0 
OPEN (999, file = 'eos.table') 
DO 
    READ (999,*, END=10) 
    nlines = nlines + 1 
END DO 
10 CLOSE (999) 

! Allocate arrays !
ALLOCATE(nbtable(0:nlines-1))
ALLOCATE(epstable(0:nlines-1))
ALLOCATE(ptable(0:nlines-1))

! Read !
OPEN(UNIT=999, FILE = 'eos.table', ACTION='READ')
DO i = 1, nlines
	READ(999,*) dummy, nbtable(i-1), dummy, ptable(i-1), epstable(i-1), dummy
	!IF(ptable(i-1) < 0.0D0) WRITE(*,*) ptable(i-1), i-1
ENDDO
CLOSE(999)

! Subtract the rest mass contribution from epsilon !
DO i = 1, nlines
	epstable(i-1) = epstable(i-1) - mnuev
ENDDO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine destroy arrays alllocated for storing neutron star EOS table 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DESTROYNS
USE DEFINITION
IMPLICIT NONE

! Deallocate arrays !
DEALLOCATE(nbtable)
DEALLOCATE(epstable)
DEALLOCATE(ptable)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine convert density to pressure by using neutron star EOS tables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE NSRtoP(rho_in, p_out)
USE DEFINITION
IMPLICIT NONE

! Input and output !
REAL (DP), INTENT(IN) :: rho_in
REAL (DP), INTENT(OUT) :: p_out

! Integer !
INTEGER :: left, right, m

! Target density !
REAL (DP) :: rhotarget

! Convert to number density !
rhotarget = ((rho_in/density)/cm3tofm3/mnu)

! Case by case !
IF(rhotarget == nbtable(0)) THEN

	! Table minimum !
	p_out = ptable(0)

ELSE

	! Binary search !
	left = 0
	right = nlines - 1
	DO
		IF(left > right) THEN
			EXIT
		END IF
		m = int((right + left)/2)
		IF(nbtable(m) < rhotarget) THEN
			left = m + 1	
		ELSEIF(nbtable(m) > rhotarget) THEN
			right = m - 1
		END IF
	END DO
	m = m - 1

	! Interpolate !
	IF(m >= 2) THEN
		CALL AKIMA(nbtable(m-2), nbtable(m-1), nbtable(m), nbtable(m+1), nbtable(m+2), nbtable(m+3), &
		ptable(m-2), ptable(m-1), ptable(m), ptable(m+1), ptable(m+2), ptable(m+3), rhotarget, p_out)
	ELSE
		CALL LINEAR(nbtable(m), nbtable(m+1), ptable(m), ptable(m+1), rhotarget, p_out)
	END IF

END IF

! Convert to CGS !
p_out = (p_out*mev2erg*cm3tofm3)*(pressure)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine convert density to internal energy by using neutron star EOS tables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE NSRtoE(rho_in, eps_out)
USE DEFINITION
IMPLICIT NONE

! Input and output !
REAL (DP), INTENT(IN) :: rho_in
REAL (DP), INTENT(OUT) :: eps_out

! Integer !
INTEGER :: left, right, m

! Target density !
REAL (DP) :: rhotarget

! Convert to number density !
rhotarget = ((rho_in/density)/mnu/cm3tofm3)

! Case by case !
IF(rhotarget == nbtable(0)) THEN

	! Table minimum !
	eps_out = epstable(0)

ELSE

	! Binary search !
	left = 0
	right = nlines - 1
	DO
		IF(left > right) THEN
			EXIT
		END IF
		m = int((right + left)/2)
		IF(nbtable(m) < rhotarget) THEN
			left = m + 1	
		ELSEIF(nbtable(m) > rhotarget) THEN
			right = m - 1
		END IF
	END DO
	m = m - 1

	! Interpolate !
	IF(m >= 2) THEN
		CALL AKIMA(nbtable(m-2), nbtable(m-1), nbtable(m), nbtable(m+1), nbtable(m+2), nbtable(m+3), &
		epstable(m-2), epstable(m-1), epstable(m), epstable(m+1), epstable(m+2), epstable(m+3), rhotarget, eps_out)
	ELSE
		CALL LINEAR(nbtable(m), nbtable(m+1), epstable(m), epstable(m+1), rhotarget, eps_out)
	END IF
	
END IF

! Convert to CGS !
eps_out = (rhotarget*(eps_out*mev2erg/mnu))*(pressure)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine convert pressure to density by using neutron star EOS tables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE NSPtoR(p_in, rho_out)
USE DEFINITION
IMPLICIT NONE

! Input and output !
REAL (DP), INTENT(IN) :: p_in
REAL (DP), INTENT(OUT) :: rho_out

! Integer !
INTEGER :: left, right, m

! Target density !
REAL (DP) :: ptarget

! Convert to number density !
ptarget = ((p_in/pressure)/cm3tofm3/mev2erg)

! Case by case !
IF(ptarget == ptable(0)) THEN

	! Table minimum !
	rho_out = nbtable(0)
	
ELSE

	! Binary search !
	left = 0
	right = nlines - 1
	DO
		IF(left > right) THEN
			EXIT
		END IF
		m = int((right + left)/2)
		IF(ptable(m) < ptarget) THEN
			left = m + 1	
		ELSEIF(ptable(m) > ptarget) THEN
			right = m - 1
		END IF
	END DO
	m = m - 1

	! Interpolate !
	IF(m >= 2) THEN
		CALL AKIMA(ptable(m-2), ptable(m-1), ptable(m), ptable(m+1), ptable(m+2), ptable(m+3), &
		nbtable(m-2), nbtable(m-1), nbtable(m), nbtable(m+1), nbtable(m+2), nbtable(m+3), ptarget, rho_out)
	ELSE
		CALL LINEAR(ptable(m), ptable(m+1), nbtable(m), nbtable(m+1), ptarget, rho_out)
	END IF
	
END IF

! Convert to CGS !
rho_out = (rho_out*mnu*cm3tofm3)*(density)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine convert pressure to density by using neutron star EOS tables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE NSPtoE(p_in, rho_in, eps_out)
USE DEFINITION
IMPLICIT NONE

! Input and output !
REAL (DP), INTENT(IN) :: p_in, rho_in
REAL (DP), INTENT(OUT) :: eps_out

! Integer !
INTEGER :: left, right, m

! Target density !
REAL (DP) :: ptarget

! Convert to number density !
ptarget = ((p_in/pressure)/cm3tofm3/mev2erg)

! Case by case !
IF(ptarget == ptable(0)) THEN

	! Table minimum !
	eps_out = epstable(0)

ELSE

	! Binary search !
	left = 0
	right = nlines - 1
	DO
		IF(left > right) THEN
			EXIT
		END IF
		m = int((right + left)/2)
		IF(ptable(m) < ptarget) THEN
			left = m + 1	
		ELSEIF(ptable(m) > ptarget) THEN
			right = m - 1
		END IF
	END DO
	m = m - 1

	! Interpolate !
	IF(m >= 2) THEN
		CALL AKIMA(ptable(m-2), ptable(m-1), ptable(m), ptable(m+1), ptable(m+2), ptable(m+3), &
		epstable(m-2), epstable(m-1), epstable(m), epstable(m+1), epstable(m+2), epstable(m+3), ptarget, eps_out)
	ELSE
		CALL LINEAR(ptable(m), ptable(m+1), epstable(m), epstable(m+1), ptarget, eps_out)
	END IF

END IF

! Convert to CGS !
eps_out = ((rho_in/density)*(eps_out*mev2erg/mnu))*(pressure)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine do 2-point interpolation for a given point of existing 
!datum and output the quantity that the user wished to interpolate 	 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LINEAR(x0, x1, y0, y1, x_in, y_out)
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! Input parameter !
REAL (DP), INTENT(IN) :: x0, x1, y0, y1, x_in

! Input parameter !
REAL (DP), INTENT(OUT) :: y_out

! Intermediate polynominal !
REAL (DP) :: nu0, nu1
REAL (DP) :: de0, de1
REAL (DP) :: l0, l1

! Assign numerator !
nu0 = (x_in - x1)
nu1 = (x_in - x0)

! Assign denominator !
de0 = (x0 - x1)
de1 = (x1 - x0)

! Assign polynominal coefficient !
l0 = nu0/de0
l1 = nu1/de1

! Compute the output !
y_out = l0*y0 + l1*y1

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Akima spline interpolation. See Hiroshi Akima 1970 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AKIMA(xm2, xm1, x0, xp1, xp2, xp3, ym2, ym1, y0, yp1, yp2, yp3, x_in, y_out)
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! Input parameter !
REAL (DP), INTENT(IN) :: xm2, xm1, x0, xp1, xp2, xp3, ym2, ym1, y0, yp1, yp2, yp3, x_in

! Input parameter !
REAL (DP), INTENT(OUT) :: y_out

! Intermediate polynominal !
REAL (DP) :: dm2, dm1, d0, dp1, dp2
REAL (DP) :: s0, s1

! Weights !
REAL (DP) :: w1, w2, w3, w4

! Coefficient of polynominal !
REAL (DP) :: p0, p1, p2, p3

! Temporal arrays !
REAL (DP) :: diff, temp

! Assign slopes !
dm2 = (ym1 - ym2)/(xm1 - xm2)
dm1 = (y0 - ym1)/(x0 - xm1)
d0 = (yp1 - y0)/(xp1 - x0)
dp1 = (yp2 - yp1)/(xp2 - xp1)
dp2 = (yp3 - yp2)/(xp3 - xp2)

! Assign weights !
w1 = abs(dp1 - d0)
w2 = abs(dm1 - dm2)
w3 = abs(dp2 - dp1)
w4 = abs(d0 - dm1)

! assign slopes !
IF(w1 == 0.0D0 .AND. w2 == 0.0D0) THEN
	s0 = 0.5D0*(dm1 + d0)
ELSE
	s0 = (w1*dm1 + w2*d0)/(w1 + w2)
END IF
IF(w3 == 0.0D0 .AND. w4 == 0.0D0) THEN
	s1 = 0.5D0*(d0 + dm1)
ELSE
	s1 = (w3*d0 + w4*dp1)/(w3 + w4)
END IF

! Assign temp !
diff = xp1 - x0
temp = x_in - x0

! assign coefficients !
p0 = y0
p1 = s0
p2 = (3.0D0*d0 - 2.0D0*s0 - s1)/diff
p3 = (s0 + s1 - 2.0D0*d0)/diff**2

! Output the interpolation !
y_out = p0 + p1*temp + p2*temp**2 + p3*temp**3

END SUBROUTINE