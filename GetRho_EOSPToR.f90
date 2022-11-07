!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine aims at finding the density corresponds !
! to a given pressure in the case of ideal degenerate     !
! fermi gas equation of state by interpolation            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETRHO_EOSPTOR (ini_rho, ini_p, eosline, type)
USE DEFINITION
IMPLICIT NONE

! the input pressure !
REAL (DP), INTENT (IN) :: ini_p

! Integer parameter !
INTEGER, INTENT (IN) :: type

! the outputed density !
REAL (DP), INTENT (OUT) :: ini_rho

! distinguish between NM and DM !
INTEGER, INTENT (INOUT) :: eosline

! Other essential variables !
INTEGER :: i
REAL (DP) :: x1, x2, x3, x4, y1, y2, y3, y4

! We find the density according to the type, 1 is DM, 2 is NM !
IF (type == 1) THEN

	! We locate the first where the pressure located !
	i = eosline	

	DO
		IF (ini_p > eostable1 (i, 1)) THEN
			i = i + 1
		ELSE
			IF (ini_p > eostable1 (i - 1, 1)) THEN
				x1 = eostable1 (i - 2, 1)
				x2 = eostable1 (i - 1, 1)
				x3 = eostable1 (i, 1)
				x4 = eostable1 (i + 1, 1)
				y1 = eostable1 (i - 2, 2)
				y2 = eostable1 (i - 1, 2)
				y3 = eostable1 (i, 2)
				y4 = eostable1 (i + 1, 2)
				eosline = i
				EXIT
			ELSE
				i = i - 1
			END IF
		END IF
	END DO

	! We do the interpolation !
	ini_rho = ((ini_p - x2) * (ini_p - x3) * (ini_p - x4)) / ((x1 - x2) * (x1 - x3) * (x1 - x4)) * y1 &
		+ ((ini_p - x1) * (ini_p - x3) * (ini_p - x4)) / ((x2 - x1) * (x2 - x3) * (x2 - x4)) * y2 &
		+ ((ini_p - x1) * (ini_p - x2) * (ini_p - x4)) / ((x3 - x1) * (x3 - x2) * (x3 - x4)) * y3 &
		+ ((ini_p - x1) * (ini_p - x2) * (ini_p - x3)) / ((x4 - x1) * (x4 - x2) * (x4 - x3)) * y4
	IF(x1 == x2 .OR. x2 == x3 .OR. x3 == x4) THEN
		STOP 'Error Table DM' 
	END IF

ELSEIF (type == 2) THEN
	
	! We locate the first where the pressure located !
	i = eosline

	DO
		IF (ini_p > eostable2 (i, 1)) THEN
			i = i + 1
		ELSE
			IF (ini_p > eostable2 (i - 1, 1)) THEN
				x1 = eostable2 (i - 2, 1)
				x2 = eostable2 (i - 1, 1)
				x3 = eostable2 (i, 1)
				x4 = eostable2 (i + 1, 1)
				y1 = eostable2 (i - 2, 2)
				y2 = eostable2 (i - 1, 2)
				y3 = eostable2 (i, 2)
				y4 = eostable2 (i + 1, 2)
				eosline = i
				EXIT
			ELSE
				i = i - 1
			END IF
		END IF
	END DO

	! We do the interpolation !
	ini_rho = ((ini_p - x2) * (ini_p - x3) * (ini_p - x4)) / ((x1 - x2) * (x1 - x3) * (x1 - x4)) * y1 &
		+ ((ini_p - x1) * (ini_p - x3) * (ini_p - x4)) / ((x2 - x1) * (x2 - x3) * (x2 - x4)) * y2 &
		+ ((ini_p - x1) * (ini_p - x2) * (ini_p - x4)) / ((x3 - x1) * (x3 - x2) * (x3 - x4)) * y3 &
		+ ((ini_p - x1) * (ini_p - x2) * (ini_p - x3)) / ((x4 - x1) * (x4 - x2) * (x4 - x3)) * y4
	IF(x1 == x2 .OR. x2 == x3 .OR. x3 == x4) THEN
		STOP 'Error Table NM'
	END IF
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!