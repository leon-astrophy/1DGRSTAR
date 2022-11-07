!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve TOV equations for hydrostatic stars !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM MASSRADIUS
USE DEFINITION
IMPLICIT NONE

! parameter !
INTEGER, parameter :: n_p = 100
INTEGER, parameter :: n_f = 10

! Integer !
INTEGER :: i, j, n, m, q

! Real !
REAL (DP) :: check_m_last, check_m, drho1c

! Array !
REAL (DP), DIMENSION (1:n_p) :: factor 
REAL (DP), DIMENSION (1:n_f) :: m_target 

! for file input output !
character(LEN=10) :: p_file
character(LEN=10) :: f_file
real*8:: p_mass, f_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize !
IF(ns_analytic == 1) THEN
	CALL ANALYTIC_TABLE
ELSE
	CALL EOSTABLE_NM
END IF
CALL BUILDHYDRO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For Neutron Star EOS !
IF(ns_flag == 1) THEN
	CALL READEOSTABLE
ELSE
	! For NM !
	OPEN (UNIT = 100, FILE = 'EOS_Table2.eos', STATUS = 'OLD')
	DO i = 1, eoslineno2
		READ (100, *) eostable2 (i, 1), eostable2 (i, 2)
	END DO
	eosline2 = 1
	CLOSE (100)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! assign !
DO i = 1, n_p
	factor(i) = 0.1D0 + 0.01D0*DBLE(i)
END DO
DO i = 1, n_f
	m_target(i) = 0.01D0 + 0.001D0*DBLE(i-1)
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now do for DM !
IF(DM_flag == 1) THEN

	! Print !
	WRITE (*,*) 'Solving For 2F TOV ...'
	WRITE (*,*)	

	! Loop over DM particle mass !
	DO m = 1, n_p

		! get particle mass !
		mb1 = mb1_def*factor(m)
		me1 = me1_def*factor(m)
		a_max1 = (me1**4)/(2.4D1*pi_old**2*h_bar**3)

		! We read the EOS table for DM !
		CALL EOSTABLE_DM
		OPEN (UNIT = 99, FILE = 'EOS_Table1.eos', STATUS = 'OLD')
		DO i = 1, eoslineno1
			READ (99, *) eostable1 (i, 1), eostable1 (i, 2)
		END DO
		eosline1 = 1
		CLOSE (99)

		! Loop over DM fluid mass !
		DO q = 1, n_f

			p_mass = factor(m)
			f_mass = m_target(q)

			write(p_file,'(f0.3)') factor(m)
			write(f_file,'(f0.3)') m_target(q)

			! Openfile !
			OPEN (UNIT = 101, FILE = './Outfile/MR-Relation-Particle-'// trim(adjustl(p_file)) //'-Fluid-'// trim(adjustl(f_file)) //'.dat', STATUS = 'REPLACE')

			! Write header !
			WRITE (101,701) p_mass, f_mass

			! Compute M-Rho relations !
			DO j = 1, n_rho
		
				! Assign NM maximum density !
				rho2_c = rhostart + (DBLE(j) - 1.0D0)*drho
				rho2_c = (1.0D1**(rho2_c)*density)

				! Guess an initial density !
				log10rho1_c = log10(rho2_c)

				! Set step size !
				drho1c = 0.1D0
		
				! Bisection method !
				DO i = 0, n_max

					! Assign central density !
					rho1_c = 10.0D0**(log10rho1_c)			

					! Solve the equilbrium structure !
					CALL GETRHO2F
	
					! Check if the energy is balanced
					check_m = mass1 - m_target(q)

					! Make sure you go to the right direction of dtemp
					if(i == 0) then
						if(check_m > 0.0D0) THEN
							drho1c = -drho1c
						end if
					endif
		
					! Use shooting method    
					if(check_m_last * check_m < 0.0D0 .and. i /= 0) then
						drho1c = -drho1c * 0.5D0
					end if

					! Update !
					log10rho1_c = log10rho1_c + drho1c
					check_m_last = check_m
			
					WRITE (*,*) i, j, mass1, m_target(q)

					! Exit condition !
					IF(ABS(check_m/m_target(q)) < tor) THEN
						EXIT
					END IF

				END DO

				! Print results !
				WRITE (101,701) log10(rho2_c/density), mass1, mass2, rad1, rad2

			END DO

			! Close file !
			CLOSE (101)

		END DO
	END DO
ELSE

	! Print !
	WRITE (*,*) 'Solving For 1F TOV ...'
	WRITE (*,*)	

	DO j = 1, n_rho
		WRITE (*,*) j
		! Assign NM maximum density !
		rho2_c = rhostart + (DBLE(j) - 1.0D0)*drho
		rho2_c = (1.0D1**(rho2_c)*density)	

		! Solve the equilbrium structure !
		CALL GETRHO1F

	END DO	
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Destroy arrays !
CALL DESTROYHYDRO

! Neutron Star EOS !
IF(ns_flag == 1) THEN
	CALL DESTROYNS
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Print !
WRITE (*,*) 'DONE!'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Format !
701 FORMAT (10ES33.15)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM
