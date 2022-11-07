SUBROUTINE ANALYTIC_TABLE
USE DEFINITION, only : density, pressure, eoslineno2
IMPLICIT NONE

! starting and ending density !
INTEGER, PARAMETER :: start = 1
INTEGER, PARAMETER :: end = 16
INTEGER, PARAMETER :: width = 5000

! fitting coefficient !
REAL*8, DIMENSION(1:12) :: c_fit 
REAL*8, DIMENSION(1:11) :: a_fit

! dummy integer !
INTEGER :: i, j

! fitting variables !
REAL*8 :: logrho, logp
REAL*8 :: cat, cat_low, cat_high, dummy, drho, x_dummy1, x_dummy2, x_dummy3, x_dummy4

! Initialize !
c_fit = (/10.6557D0, 3.7863D0, 0.8124D0, 0.6823D0, 3.5279D0, 11.8100D0, 12.0584d0, 1.4663D0, 3.4952D0, 11.8007d0, 14.4114d0, 14.4081d0/)
a_fit = (/4.3290D0, 4.3622D0, 9.1131D0, -0.4751D0, 3.4614D0, 14.8800D0, 21.3141D0, 0.1023D0, 0.0495D0, 4.9401D0, 10.2957D0/)

! We open the EOS table for input !
OPEN (UNIT = 100, FILE = 'EOS_Table2.eos', STATUS = 'REPLACE')

! Loop !
DO i = start, end
	DO j = 0, width - 1
        logrho = (10.0D0 ** (DBLE(i+1)) - 10.0D0 ** (DBLE(i))) * DBLE(j) / DBLE(width) + 10.0D0 ** (DBLE(i))
        logrho = log10(logrho)

        x_dummy1 = c_fit(5)*(logrho - c_fit(6))
        x_dummy2 = c_fit(9)*(c_fit(10) - logrho)
        x_dummy3 = (c_fit(1) + c_fit(2)*(logrho - c_fit(3))**(c_fit(4)))
        x_dummy4 = (c_fit(7) + c_fit(8)*logrho)
        cat_low = x_dummy3*fitting(x_dummy1) + x_dummy4*fitting(x_dummy2)
    
        x_dummy1 = a_fit(5)*(a_fit(6) - logrho)
        x_dummy2 = a_fit(10)*(a_fit(11) - logrho)
        x_dummy3 = (a_fit(3) + a_fit(4)*logrho)
        x_dummy4 = (a_fit(7) + a_fit(8)*logrho + a_fit(9)*logrho**2)
        cat_high = x_dummy3*fitting(x_dummy1) + x_dummy4*fitting(x_dummy2)

        x_dummy1 = a_fit(1)*(logrho - c_fit(11))
        x_dummy2 = a_fit(2)*(c_fit(12) - logrho)
        cat = cat_low*fitting(x_dummy1) + cat_high*fitting(x_dummy2)

        logp = cat

        ! print !
        WRITE (100,*) 10.0D0**logp*pressure, 10.0D0**logrho*density
    END DO
END DO

! We assign the eosline number !
eoslineno2 = width * (end - start + 1)

! We close the file !
close(100)

contains 

	real*8 function fitting(x_in)
	implicit none
	real*8 :: x_in
	fitting = 1.0D0/(1.0D0 + exp(x_in))
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ANALYTIC_RTOP (ini_p, ini_rho)
USE DEFINITION, only : density, pressure
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: ini_rho
REAL*8, INTENT(OUT) :: ini_p

! fitting coefficient !
REAL*8, DIMENSION(1:12) :: c_fit 
REAL*8, DIMENSION(1:11) :: a_fit

! dummy integer !
INTEGER :: i, j

! fitting variables !
REAL*8 :: logrho, logp
REAL*8 :: cat, cat_low, cat_high, dummy, drho, x_dummy1, x_dummy2, x_dummy3, x_dummy4

! Initialize !
c_fit = (/10.6557D0, 3.7863D0, 0.8124D0, 0.6823D0, 3.5279D0, 11.8100D0, 12.0584d0, 1.4663D0, 3.4952D0, 11.8007d0, 14.4114d0, 14.4081d0/)
a_fit = (/4.3290D0, 4.3622D0, 9.1131D0, -0.4751D0, 3.4614D0, 14.8800D0, 21.3141D0, 0.1023D0, 0.0495D0, 4.9401D0, 10.2957D0/)

logrho = log10(ini_rho/density)

x_dummy1 = c_fit(5)*(logrho - c_fit(6))
x_dummy2 = c_fit(9)*(c_fit(10) - logrho)
x_dummy3 = (c_fit(1) + c_fit(2)*(logrho - c_fit(3))**(c_fit(4)))
x_dummy4 = (c_fit(7) + c_fit(8)*logrho)
cat_low = x_dummy3*fitting(x_dummy1) + x_dummy4*fitting(x_dummy2)
    
x_dummy1 = a_fit(5)*(a_fit(6) - logrho)
x_dummy2 = a_fit(10)*(a_fit(11) - logrho)
x_dummy3 = (a_fit(3) + a_fit(4)*logrho)
x_dummy4 = (a_fit(7) + a_fit(8)*logrho + a_fit(9)*logrho**2)
cat_high = x_dummy3*fitting(x_dummy1) + x_dummy4*fitting(x_dummy2)

x_dummy1 = a_fit(1)*(logrho - c_fit(11))
x_dummy2 = a_fit(2)*(c_fit(12) - logrho)
cat = cat_low*fitting(x_dummy1) + cat_high*fitting(x_dummy2)

logp = cat
ini_p = 10.0D0**(cat)*pressure

contains 

	real*8 function fitting(x_in)
	implicit none
	real*8 :: x_in
	fitting = 1.0D0/(1.0D0 + exp(x_in))
	end function

END SUBROUTINE