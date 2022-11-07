!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine allocate arrays related to hydrodynamics variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDHYDRO
USE DEFINITION
IMPLICIT NONE

! DM hydrodynamic variables !
ALLOCATE(rho1(-4 : length_step + 5))
ALLOCATE(epsilon1(-4 : length_step + 5))
ALLOCATE(eostable1(eoslineno1,2))

! NM hydrodynamic variables !
ALLOCATE(rho2(-4 : length_step + 5))
ALLOCATE(epsilon2(-4 : length_step + 5))
ALLOCATE(eostable2(eoslineno2,2))

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocate arrays related to hydro variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DESTROYHYDRO
USE DEFINITION
IMPLICIT NONE

! DM hydrodynamic variables !
DEALLOCATE(rho1)
DEALLOCATE(epsilon1)
DEALLOCATE(eostable1)

! NM hydrodynamic variables !
DEALLOCATE(rho2)
DEALLOCATE(epsilon2)
DEALLOCATE(eostable2)

END SUBROUTINE