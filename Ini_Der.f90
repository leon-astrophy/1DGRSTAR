!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine is the differential equation governing the hydrostatic star !
! The differential equation represent hydro static equilibrium of two fluid   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INI_DER_1F (der, x, y, no_of_eq_ini)
USE DEFINITION
IMPLICIT NONE

! Integer parameter governing the no of equation !
INTEGER, INTENT (IN) :: no_of_eq_ini

! Real parameter !
REAL (DP), INTENT (IN) :: x
REAL (DP), INTENT (IN), DIMENSION (1 : no_of_eq_ini) :: y
REAL (DP), INTENT (OUT), DIMENSION (1 : no_of_eq_ini) :: der

! The differential equation of two fluid hydrostatic equilibrium!
der (1) = 4.0E0_DP * pi_old * x ** 2 * y (3)
der (2) = - (y(2) + y(3))*(y(1) + 4.0D0*pi_old*x**3*y(2))/(x*(x - 2.0D0*y(1)))
der (3) = 0.0E0_DP

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine is the differential equation governing the hydrostatic star !
! The differential equation represent hydro static equilibrium of two fluid   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INI_DER_2F (der, x, y, no_of_eq_ini)
USE DEFINITION
IMPLICIT NONE

! Integer parameter governing the no of equation !
INTEGER, INTENT (IN) :: no_of_eq_ini

! Real parameter !
REAL (DP), INTENT (IN) :: x
REAL (DP), INTENT (IN), DIMENSION (1 : no_of_eq_ini) :: y
REAL (DP), INTENT (OUT), DIMENSION (1 : no_of_eq_ini) :: der

! The differential equation of two fluid hydrostatic equilibrium!
der (1) = 4.0E0_DP * pi_old * x ** 2 * y (3)
der (2) = - (y(2) + y(3))*((y(1) + y(4)) + 4.0D0*pi_old*x**3*(y(2) + y(5)))/(x*(x - 2.0D0*(y(1) + y(4))))
der (3) = 0.0E0_DP
der (4) = 4.0E0_DP * pi_old * x ** 2 * y (6)
der (5) = - (y(5) + y(6))*((y(1) + y(4)) + 4.0D0*pi_old*x**3*(y(2) + y(5)))/(x*(x - 2.0D0*(y(1) + y(4))))
der (6) = 0.0E0_DP

END SUBROUTINE