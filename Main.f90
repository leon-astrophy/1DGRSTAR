!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve TOV equations for hydrostatic stars !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM MASSRADIUS
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, n, m, q

! Real !
REAL (DP) :: check_m_last, check_m, drho1c

! for file input output !
character(LEN=10) :: mb_file
character(LEN=10) :: mdm_file
real*8:: mb_mass, mdm_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! section for initializing and building arrays !

! Initialize !
IF(ns_analytic == 1) THEN
	CALL ANALYTIC_TABLE
ELSE
	CALL EOSTABLE_NM
END IF

! Build hydro variables !
CALL BUILDHYDRO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! section for constructing EOS table !

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
! Section for the main functions !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now do for the DM case !
IF(DM_flag == 1) THEN

	! Print !
	WRITE (*,*) 'Solving For 2F TOV ...'
	WRITE (*,*)	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Allocate arrays !
	ALLOCATE(y_rk4(1 : 3, -4 : length_step + 5))

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Loop over DM particle mass !
	DO m = 1, n_p

		! get particle mass !
		mb1 = mbstart + (DBLE(m) - 1.0D0)*dnb
		mb_mass = mb1/GeV/mass
		me1 = mb1

		! Scaling factor for 
		a_max1 = (me1**4)/(2.4D1*pi_old**2*h_bar**3)

		! We read the EOS table for DM !
		CALL EOSTABLE_DM
		OPEN (UNIT = 99, FILE = 'EOS_Table1.eos', STATUS = 'OLD')
		DO i = 1, eoslineno1
			READ (99, *) eostable1 (i, 1), eostable1 (i, 2)
		END DO
		eosline1 = 1
		CLOSE (99)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Loop over DM fluid mass !
		DO q = 1, n_dm

			! DM fluid mass !
			mdm_mass = mdmstart + (DBLE(q) - 1.0D0)*dndm

			write(mb_file,'(f0.4)') mb_mass
			write(mdm_file,'(f0.4)') mdm_mass

			! the leading zero is missing, add it back !
			if(mb_file(1:1) == '.') mb_file = '0' // mb_file
			if(mb_file(1:2) == '-.') mb_file = '-0.' // mb_file(3:)

			if(mdm_file(1:1) == '.') mdm_file = '0' // mdm_file
			if(mdm_file(1:2) == '-.') mdm_file = '-0.' // mdm_file(3:)

			! Openfile !
			OPEN (UNIT = 101, FILE = './outfile/MR-Relation-Particle-'// trim(adjustl(mb_file)) //'-Fluid-'// trim(adjustl(mdm_file)) //'.dat', STATUS = 'REPLACE')

			! Write header !
			WRITE (101,701) mb_mass, mdm_mass

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
	
					! Check if the deviation from the target DM mass !
					check_m = mass1 - mdm_mass

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

					! Exit condition !
					IF(ABS(check_m/mdm_mass) < tor) THEN
						EXIT
					END IF

				END DO

				! Print results !
				WRITE (101,701) log10(rho2_c/density), mass1, mass2, rad1, rad2

			END DO
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! Close file !
			CLOSE (101)

		END DO
	END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now do for the NM case !
ELSE

	! Print !
	WRITE (*,*) 'Solving For 1F TOV ...'
	WRITE (*,*)	

	! Allocate arrays !
	ALLOCATE(y_rk4(1 : 6, -4 : length_step + 5))

	! Openfile !
	OPEN (UNIT = 101, FILE = './outfile/MR-Relation-NM.dat', STATUS = 'REPLACE')

	! Loop over the densit yrange !
	DO j = 1, n_rho

		! Assign NM maximum density !
		rho2_c = rhostart + (DBLE(j) - 1.0D0)*drho
		rho2_c = (1.0D1**(rho2_c)*density)	

		! Solve the equilbrium structure !
		CALL GETRHO1F

	END DO	

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reaching the end of the calcualtions, deallocate all arrays !
 
! Destroy arrays !
CALL DESTROYHYDRO

! Neutron Star EOS !
IF(ns_flag == 1) THEN
	CALL DESTROYNS
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final section !

! Print !
WRITE (*,*) 'DONE!'

! Format !
701 FORMAT (10ES33.15)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM
