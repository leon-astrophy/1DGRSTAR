!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the parameter files containing constants necessary for computing MR relations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define double precision !
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND (15, 307)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Conversion between units !

! Mathematical constants and physical constants !
REAL (DP), PARAMETER :: pi_old = 3.1415926535897932384626433832795E0_DP

! Physical constants to be as one !
REAL (DP), PARAMETER :: gconst = 6.67430D-8
REAL (DP), PARAMETER :: clight = 2.99792458D10
REAL (DP), PARAMETER :: solar = 1.98847D33
REAL (DP), PARAMETER :: rsolar = 6.96342D10

! Conversion between units !
REAL (DP), PARAMETER :: length = (clight**2)/(solar*gconst)
REAL (DP), PARAMETER :: mass = (1.0D0/solar)
REAL (DP), PARAMETER :: time = (clight**3)/(solar*gconst)

! Derived conversion !
REAL (DP), PARAMETER :: density = (mass/length**3)
REAL (DP), PARAMETER :: epsilon = (1.0D0/clight**2)
REAL (DP), PARAMETER :: h_bar = (1.054571817D-27)*(length**2*mass/time)
REAL (DP), PARAMETER :: pressure = density*epsilon
REAL (DP), PARAMETER :: qdot = pressure/time

! 1 GeV !
REAL (DP), PARAMETER :: GeV = 1.78266191D-24

! We use GK as default temperature unit !
REAL (DP), PARAMETER :: temperature = 1.0D-9

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section is related to the basic parameters governing the simulation !

! Physical length (dimensionless) of the simulation box !
REAL (DP), PARAMETER :: total_length = 2.0E4_DP

! Value of spatial grid size dx !
REAL (DP), PARAMETER :: dx = 1.0E-1_DP

! The total number of array stored by each variables !
INTEGER, PARAMETER :: length_step = INT (total_length / dx)      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! This section governs the physics of the EOS of NM or DM    !
! CAUTION : We assumed an ideal completely degenerate fermi  !
! gas EOS. To change the EOS, you need to input the required !
! parameters by yourself. For example, temperature           !     

! Baryonic mass for normal matter !			
REAL (DP), PARAMETER :: mb2 = 1.66053906660D-24*mass

! Fermionic mass (electrons) for normal matter !			
REAL (DP), PARAMETER :: me2 = 9.1093837015D-28*mass

! Electron fraction for normal matter !
REAL (DP), PARAMETER :: ye2 = 5.0E-1_DP

! Multiplicity factor for normal matter				
REAL (DP), PARAMETER :: gs2 = 2.0E0_DP

! For fermi gas EOS !
REAL (DP), PARAMETER :: a_max2 = (me2**4)/(2.4D1*pi_old**2*h_bar**3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for neutron Star EOS !

! Neutron mass !
REAL (DP), PARAMETER :: mnu = 1.6749274980495D-24
REAL (DP), PARAMETER :: mnuev = 939.56533000000002D0

! Conversion to cubic meter !
REAL (DP), PARAMETER :: cm3tofm3 = 1.0D39
REAL (DP), PARAMETER :: mev2erg = 1.60217662D-6

! Want Neutron Star ? !
INTEGER, PARAMETER :: ns_flag = 0
INTEGER, PARAMETER :: ns_analytic = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for NM central density

! Atmospeheric density factor !
REAL (DP), PARAMETER :: rhofac = 1.0D-10

! Number of models per maximum density !
INTEGER, PARAMETER :: n_rho = 101

! Starting log maximum density !
REAL (DP), PARAMETER :: rhostart = 9.0D0

! Ending log maximum density !
REAL (DP), PARAMETER :: rhoend = 15.0D0

! Step size !
REAL (DP), PARAMETER :: drho = (rhoend - rhostart)/(DBLE(n_rho) - 1.0D0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for DM parameters

! Want DM? !
INTEGER, PARAMETER :: dm_flag = 1

! Maximum bisection step !
INTEGER, PARAMETER :: n_max = 1000

! Tolerance !
REAL (DP), PARAMETER :: tor = 1.0D-6

! Electron fraction for normal matter !
REAL (DP), PARAMETER :: ye1 = 1.0E0_DP

! Multiplicity factor for normal matter				
REAL (DP), PARAMETER :: gs1 = 2.0E0_DP	

! Starting DM particle mass for parameter search 
REAL (DP), PARAMETER :: mbstart = 1.78266191D-25*mass

! Ending DM particle mass for parameter search 
REAL (DP), PARAMETER :: mbend = 2.0d0*1.78266191D-25*mass

! Number of points for the range of DM particle mass !
INTEGER, parameter :: n_p = 2

! Step size !
REAL (DP), PARAMETER :: dnb = (rhoend - rhostart)/(DBLE(n_p) - 1.0D0)

! Starting DM fluid mass for parameter search 
REAL (DP), PARAMETER :: mdmstart = 0.05D0

! Ending DM fluid mass for parameter search 
REAL (DP), PARAMETER :: mdmend = 0.1D0

! Number of points for the range of DM fluid mass !
INTEGER, parameter :: n_dm = 2

! Step size !
REAL (DP), PARAMETER :: dndm = (mdmend - mdmstart)/(DBLE(n_dm) - 1.0D0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!