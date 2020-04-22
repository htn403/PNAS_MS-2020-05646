      SUBROUTINE VUMAT( 
! READ ONLY -
     $     NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL,
     $     STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH,
     $     PROPS, DENSITY, STRAININC, RELSPININC,
     $     TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,
     $     STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,
     $     TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW,
! WRITE ONLY -
     $     STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW )

      INCLUDE 'VABA_PARAM.INC'

      DIMENSION PROPS(NPROPS), DENSITY(NBLOCK),
     $COORDMP(NBLOCK,*),
     $CHARLENGTH(NBLOCK), STRAININC(NBLOCK,NDIR+NSHR),
     $RELSPININC(NBLOCK,NSHR), TEMPOLD(NBLOCK),
     $STRETCHOLD(NBLOCK,NDIR+NSHR), DEFGRADOLD(NBLOCK,NDIR+NSHR+NSHR),
     $FIELDOLD(NBLOCK,NFIELDV), STRESSOLD(NBLOCK,NDIR+NSHR),
     $STATEOLD(NBLOCK,NSTATEV), ENERINTERNOLD(NBLOCK),
     $ENERINELASOLD(NBLOCK), TEMPNEW(NBLOCK),
     $STRETCHNEW(NBLOCK,NDIR+NSHR), DEFGRADNEW(NBLOCK,NDIR+NSHR+NSHR), 
     $FIELDNEW(NBLOCK,NFIELDV),
     $STRESSNEW(NBLOCK,NDIR+NSHR), STATENEW(NBLOCK,NSTATEV),
     $ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK)
C  !
      CHARACTER*80 CMNAME
C  !

C      WRITE(*,'(A," HAS BEEN CALLED.")') "FERHUN: VUMAT"

      IF(CMNAME(1:13).EQ."M7FMATERIAL") THEN
C      WRITE(6,*) "ABOUT TO ENTER M7FMATERIAL."

         CALL M7FMATERIAL( 
     $ NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL,
     $ STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH,
     $ PROPS, DENSITY, STRAININC, RELSPININC,
     $ TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,
     $ STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,
     $ TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW,
     $ STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW )
      END IF
      RETURN
      END SUBROUTINE VUMAT
	  
C This is the model working for concrete casted at NU in July 2018, age 1 year
c Calibrated tests include compression (cylinders, prisms), Brazillian, and 3-point-bending

C +--------------------------------------------------------------------+
C |                 SUBROUTINE C_NORM_ELASTIC                          |
C +--------------------------------------------------------------------+
      SUBROUTINE C_NORM_ELASTIC(EF,DEF,SN0,EPS_N0_POS,EPS_N0_NEG,SV0,
     $SNF_NO_SPLIT,E_N) 
      INCLUDE 'VABA_PARAM.INC'
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7FBOUNDS_M7FIO/K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9,
     $     SV0_P, C_1, C_2, C_3, C_4, C_5,
     $     C_6, C_7, C_8, C_9, C_10, C_11,C_12,C_13,C_14,C_15,C_16,
     $     C_17,C_18,C_19,C_20,C_21,C_R0,C_R2
      REAL*8:: K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9, SV0_P
      REAL*8 :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8,
     $     C_9, C_10
      REAL*8 :: C_11,C_12,C_13,C_14,C_15,C_16,C_17,C_18,
     $     C_19,C_20,C_21,C_R0,C_R2
      COMMON /KM7FBOUNDS_M7FIO_1/ C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1,
     $     C_6_2, C_6_M4, C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4,
     $     C_9_0, C_9_1, C_9_2, C_9_M4
      REAL*8 :: C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1, C_6_2, C_6_M4
      REAL*8 :: C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4
      REAL*8 :: C_9_0, C_9_1, C_9_2, C_9_M4
      COMMON /KM7FBOUNDS_M7FIO_2/SIG_FIB_0, P_1,P_2, P_5, P_6,XP_1,XP_2,
     $     XP_3, V_F, CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3, CF_4,
     $     CF_5, CF_4_4, AP1, BP1, CP2, DP2, EP2, D_1, D_2, D_3, D_4,
     $     D_5, D_6
      REAL*8 :: SIG_FIB_0, P_1, P_2, P_5, P_6
      REAL*8 :: XP_1, XP_2, XP_3, V_F
      REAL*8 :: CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3
      REAL*8 :: CF_4, CF_5, CF_4_4
      REAL*8 :: AP1, BP1, CP2, DP2, EP2
      REAL*8 :: D_1, D_2, D_3, D_4, D_5, D_6
      SAVE :: /KCOMMON_VARS/, /KM7FBOUNDS_M7FIO/, /KM7FBOUNDS_M7FIO_1/,
     $     /KM7FBOUNDS_M7FIO_2/
      REAL*8, DIMENSION(1:NP) :: EF, DEF, SN0
      REAL*8, DIMENSION(1:NP) :: EPS_N0_POS, EPS_N0_NEG
      REAL*8 :: SV0
      REAL*8, DIMENSION(1:NP) :: SNF_NO_SPLIT 
      REAL*8, DIMENSION(1:NP) :: E_N
C      REAL*8, DIMENSION() :: E_N
C      REAL*8 :: C_N0, RATE_FACTOR
      REAL*8 :: E_N_0! , C_1_N , ALPHA , XC
      INTEGER :: ISIZE, I !, ALLOCSTAT
      E_N_0= YOUNG / (1.0D0 - 2.0D0*POISSON)
C THE LARGER THE STIFFNESS E_N, THE SHARPER IS THE REDUCTION IN CONTRACTION. 
C UNLOADING TO ORIGIN IN ALL MICROPLANES USING SECANT STIFFNESS
      ISIZE=SIZE(EF)
C      E_N = C_1_N * YOUNG / (1.0D0 - 2.0D0*POISSON)
      DO I=1,ISIZE
         IF ((SN0(I) > 0.D0)) THEN ! LOADING AND UNLOADING UNDER TENSION, 09.08.2011
            E_N(I) = E_N_0*EXP(-C_19*EPS_N0_POS(I))
            IF (SN0(I) > E_N_0*EF(I) .AND. SN0(I)*DEF(I)<0.0D0) THEN !20.09.2011
               E_N(I)=E_N_0                !20.09.2011
            END IF                         !20.09.2011
         ELSE ! LOADING AND UNLOADING UNDER COMPRESSION, 09.08.2011
            E_N(I) = E_N_0*(EXP(-C_20*ABS(EPS_N0_NEG(I))/
     $           (1.D0+C_18*MAX(-SV0,0.D0)/E_N_0))+
     $           C_21*MAX(-SV0,0.D0)/E_N_0)
         END IF
      END DO
C BEGIN TEST 05.08.2011
C      E_N = YOUNG / (1.0D0 - 2.0D0*POISSON)
C END TEST 05.08.2011
      SNF_NO_SPLIT = SN0*C0 + E_N*DEF
C      SNF_NO_SPLIT = SN0 + DEF*YOUNG/(1.-2.*POISSON)
      RETURN
      END SUBROUTINE C_NORM_ELASTIC

C +--------------------------------------------------------------------+
C |                         SUBROUTINE C_NORM                          |
C +--------------------------------------------------------------------+
      SUBROUTINE C_NORM(EF,SV0,R_N,SNB) 
      INCLUDE 'VABA_PARAM.INC'
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7FBOUNDS_M7FIO/K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9,
     $     SV0_P, C_1, C_2, C_3, C_4, C_5,
     $     C_6, C_7, C_8, C_9, C_10, C_11,C_12,C_13,C_14,C_15,C_16,
     $     C_17,C_18,C_19,C_20,C_21,C_R0,C_R2
      REAL*8 :: K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9, SV0_P
      REAL*8 :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8,
     $     C_9, C_10
      REAL*8 :: C_11,C_12,C_13,C_14,C_15,C_16,C_17,C_18,
     $     C_19,C_20,C_21,C_R0,C_R2
      COMMON /KM7FBOUNDS_M7FIO_1/ C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1,
     $     C_6_2, C_6_M4, C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4,
     $     C_9_0, C_9_1, C_9_2, C_9_M4
      REAL*8 :: C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1, C_6_2, C_6_M4
      REAL*8 :: C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4
      REAL*8 :: C_9_0, C_9_1, C_9_2, C_9_M4
      COMMON /KM7FBOUNDS_M7FIO_2/SIG_FIB_0, P_1, P_2, P_5,P_6,XP_1,XP_2,
     $     XP_3, V_F, CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3, CF_4,
     $     CF_5, CF_4_4, AP1, BP1, CP2, DP2, EP2, D_1, D_2, D_3, D_4,
     $     D_5, D_6
      REAL*8 :: SIG_FIB_0, P_1, P_2, P_5, P_6
      REAL*8 :: XP_1, XP_2, XP_3, V_F
      REAL*8 :: CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3
      REAL*8 :: CF_4, CF_5, CF_4_4
      REAL*8 :: AP1, BP1, CP2, DP2, EP2
      REAL*8 :: D_1, D_2, D_3, D_4, D_5, D_6
      SAVE :: /KCOMMON_VARS/, /KM7FBOUNDS_M7FIO/, /KM7FBOUNDS_M7FIO_1/,
     $     /KM7FBOUNDS_M7FIO_2/
      REAL*8, DIMENSION(1:NP) :: EF 
      REAL*8, DIMENSION(1:NP) :: SNB 
      REAL*8 :: SV0, R_N 
      REAL*8, DIMENSION(1:NP)  :: AUX1
      REAL*8 :: C_N0, RATE_FACTOR, FSTAR, BETA_N
      REAL*8 :: E_V, EB_N0, EB_NK, EB_N
      INTEGER :: ALLOCSTAT, ISIZE
      
C     D_1=2.94D-9
C     D_2=-0.384D0
C     D_3=-3.51D-5
C     D_4=-0.0203D0
C     D_5=0.29

C      C_1 = D_1*EXP(D_2*V_F / (D_3 + V_F*D_4)) + D_5

C      C_1 = D_1*EXP(D_2*V_F / (1.0D0 + V_F*D_3)) + D_4 !02.08.2010
C      C_1 = D_1*TANH(D_2*V_F - D_3) + D_4 !03.08.2010
C      WRITE(*,*) 'C_1 = ', C_1

C      C_1 = D_5

C      WRITE(*,'("C_1=",F10.5)') C_1
C      READ(*,*)

      ISIZE = SIZE(EF)
C      ALLOCATE(AUX1(1:ISIZE), STAT=ALLOCSTAT)
C      IF(ALLOCSTAT.NE.0) THEN
C        WRITE(*,*) 'CAN NOT ALLOCATE AUX1 (M7FBOUNDS.F90)'
C        STOP
C      ENDIF

      E_V    = YOUNG/(1.0D0-2.0D0*POISSON)
      EB_N0  = C_3*K_1 
      EB_NK  = C_4 !*K_1 

C IT MAY BE THAT THE FACTOR OF D_4 IN WHAT FOLLOWS MUST INCLUDE THE NEXT THERM TOO.02.08.2011
      C_1 = (D_1*TANH(D_2*V_F - D_3)+
     $  D_4*EXP(-MAX(-SV0-D_6,0.D0)/E_V*D_5))/2.0D0 !27.07.2011

C THE IF BLOCK HERE, COUPLED WITH SVF = MIN(SVF, SVFB) IN M7F.F90, WHICH EFFECTS THE EVOLUTION OF SV0 AS A
C HISTORY VARIABLE, RESULTS IN A SMOOTHLY SOFTENING UNIAXIAL COMPRESSION POST PEAK. HOWEVER, THIS IS ALSO THE
C REASON FOR THE TAIL NOT APPROACHING ZERO FAST ENOUGH. 22.07.2011
      IF (SV0.LT.0.0D0) THEN 
         EB_N = EB_N0 - EB_NK / E_V * SV0  
      ELSE 
         EB_N = EB_N0 
      END IF       
C BEGIN TEST 22.07.2011
C      EB_N = EB_N0 
C END TEST 22.07.2011

      FSTAR  = K_1*YOUNG*C_1 
      BETA_N = C_2*C_1*K_1 
C RATE EFFECT PARAMETER FOR THE NORMAL BOUNDARY      
      C_N0=C_R2 

      AUX1 = MAX(EF-BETA_N,0.0D0)/EB_N/1.D0
      SNB  = FSTAR*EXP(-AUX1)

C THE RATE FACTOR FOR VOLUMETRIC BOUNDARY BOUNDARY      
      RATE_FACTOR = C_N0*R_N
      SNB = SNB * (1.D0 + RATE_FACTOR)

C      DEALLOCATE(AUX1)
      RETURN 
      END SUBROUTINE C_NORM

! +--------------------------------------------------------------------+
! |                       SUBROUTINE C_NORM_FIB                        |
! +--------------------------------------------------------------------+
!      SUBROUTINE C_NORM_FIB(EF,SV0,R_N,SNB_FIB) 
      SUBROUTINE C_NORM_FIB(EF,R_N,SNB_FIB) 
      INCLUDE 'VABA_PARAM.INC'
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7FBOUNDS_M7FIO/K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9,
     $     SV0_P, C_1, C_2, C_3, C_4, C_5,
     $     C_6, C_7, C_8, C_9, C_10, C_11,C_12,C_13,C_14,C_15,C_16,
     $     C_17,C_18,C_19,C_20,C_21,C_R0,C_R2
      REAL*8 :: K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9, SV0_P
      REAL*8 :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8,
     $     C_9, C_10
      REAL*8 :: C_11,C_12,C_13,C_14,C_15,C_16,C_17,C_18,
     $     C_19,C_20,C_21,C_R0,C_R2
      COMMON /KM7FBOUNDS_M7FIO_1/ C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1,
     $     C_6_2, C_6_M4, C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4,
     $     C_9_0, C_9_1, C_9_2, C_9_M4
      REAL*8 :: C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1, C_6_2, C_6_M4
      REAL*8 :: C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4
      REAL*8 :: C_9_0, C_9_1, C_9_2, C_9_M4
      COMMON /KM7FBOUNDS_M7FIO_2/SIG_FIB_0, P_1, P_2, P_5,P_6,XP_1,XP_2,
     $     XP_3, V_F, CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3, CF_4,
     $     CF_5, CF_4_4, AP1, BP1, CP2, DP2, EP2, D_1, D_2, D_3, D_4,
     $     D_5, D_6
      REAL*8 :: SIG_FIB_0, P_1, P_2, P_5, P_6
      REAL*8 :: XP_1, XP_2, XP_3, V_F
      REAL*8 :: CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3
      REAL*8 :: CF_4, CF_5, CF_4_4
      REAL*8 :: AP1, BP1, CP2, DP2, EP2
      REAL*8 :: D_1, D_2, D_3, D_4, D_5, D_6
      SAVE :: /KCOMMON_VARS/, /KM7FBOUNDS_M7FIO/, /KM7FBOUNDS_M7FIO_1/,
     $     /KM7FBOUNDS_M7FIO_2/
      REAL*8, DIMENSION(1:NP) :: EF 
      REAL*8, DIMENSION(1:NP) :: SNB_FIB 
C      REAL*8 :: SV0 
      REAL*8 :: R_N 
C      REAL*8, DIMENSION()  :: AUX1
      REAL*8 :: C_N0, RATE_FACTOR !, FSTAR, BETA_N
C      REAL*8 :: E_V, EB_N0, EB_NK, EB_N
      REAL*8 :: E_N_MAX, E_N_0, E_N_HAT, EPS_SHIFT, POW !, SIG_N_MAX
      INTEGER :: ISIZE, I !, ALLOCSTAT

C =========================================
C K_1 (RADIAL SCALING) IN FIBER LAW DONE.
C E   (VERTICAL SCALING) IN FIBER LAW DONE.
C =========================================
      
      ISIZE = SIZE(EF)

      EPS_SHIFT = XP_2
      POW = XP_3
C      ALLOCATE(AUX1(1:ISIZE), STAT=ALLOCSTAT)
C      IF(ALLOCSTAT.NE.0) THEN
C        WRITE(*,*) 'CAN NOT ALLOCATE AUX1 (M7FBOUNDS.F90)'
C        STOP
C      ENDIF

C      E_V    = YOUNG/(1.0D0-2.0D0*POISSON)
C      EB_N0  = C_3*K_1 
C      EB_NK  = C_4 

C      IF (SV0.LT.0.0D0) THEN 
C         EB_N = EB_N0 - EB_NK / E_V * SV0  
C      ELSE 
C         EB_N = EB_N0 
C      END IF       

C      FSTAR  = K_1*YOUNG*C_1 
C      BETA_N = C_2*C_1*K_1 

C ==============================================================================
C E_N_MAX IS THE NORMAL STRAIN THAT CORRESPONDS TO STRENGTH AS IN THE FIBER LAW.
C EPS_SHIFT IS THE SHIFT OF THE FIBER LAW ALONG THE NORMAL STRAIN AXIS. 
C ==============================================================================
C      E_N_MAX = 1.0D0/P_2 + EPS_SHIFT
      E_N_MAX = POW/P_2 + EPS_SHIFT
C ==============================================================================
C E_N_0 IS THE NORMAL STRAIN AT WHICH THE FIBER LAW BEGINS TO DECREASE.
C E_N_0 >= E_N_MAX IN ALL CASES. IF E_N_0=E_N_MAX, CONVENTIONAL FIBER LAW IS
C RECOVERED. IF E_N_0 > E_N_MAX, AT THE STRENGTH VALUE, THERE IS A HORIZONTAL
C PLATEAU OF LENGTH (E_N_0 - E_N_MAX). THE PLATEAU IS JUSTIFIED USING THE 
C MICROMECHANICS OF UNIFORMLY SPACED FIBERS BRIDGING A CRACK. A PLATEAU OF 
C NON-ZERO LENGTH IS REQUIRED TO FIT THE DATA AND AT THE END OF THE DATA FOR 
C THE STRESSES TO GO TO ZERO. IF THE PLATEAU IS KEPT AS ZERO LENGTH AND INSTEAD 
C THE CONVENTIONAL FIBER LAW IS SCALED TO FIT THE DATA, GOOD FITS CAN BE 
C OBTAINED, BUT IMMEDIATELY BEYOND THE AVAILABLE DATA THE STRESS WILL NOT DROP 
C TO ZERO BUT WILL CONTINUE TO DROP TO ZERO GRADUALLY WITH INCREASING STRAINS. 
C THE LATTER CASE CORRESPONDS TO INCORRECT MATERIAL BEHAVIOR. 
C ==============================================================================
C      ! E_N_0 = 20000*E_N_MAX
C      E_N_0 = E_N_MAX*(1.0D0 + XP_1) ! 0.037*3
C E_N_0 SHOULD BE A FUNCTION OF V_F?
      E_N_0 = XP_1  ! >= XP_3/P_2 + XP_2 MUST BE SATISFIED, 07/29/2010
C SIG_FIB_0 IS A MISLEADING NAME. ACTUALLY IT IS THE NORMAL MICROPLANE STRAIN 
C BEYOND WHICH THE BOUNDARY IS ZERO.
      E_N_HAT = SIG_FIB_0 ! NOT USED ANY MORE... 01.08.2010

C      SIG_N_MAX=P_1/P_2*EXP(-1.0D0 + P_2/P_1*SIG_FIB_0)


C RATE EFFECT PARAMETER FOR THE NORMAL BOUNDARY      
      C_N0=C_R2 

      DO I=1, ISIZE
         IF ( EF(I)/K_1 < E_N_MAX ) THEN
C            SNB_FIB(I) = P_1*MAX(EF(I) - EPS_SHIFT,0.0D0)*EXP(-P_2*(EF(I) - EPS_SHIFT)) 
            SNB_FIB(I) = YOUNG*P_1*K_1*
     $           MAX(EF(I)/K_1 - EPS_SHIFT,0.0D0)**POW*
     $           EXP(-P_2*MAX(EF(I)/K_1 - EPS_SHIFT,0.0D0)) 
         ELSE IF ((EF(I)/K_1 >= E_N_MAX) .AND. (EF(I)/K_1 < E_N_0)) THEN
C            SNB_FIB(I) = P_1/P_2*EXP(-1.0D0)
            SNB_FIB(I) = YOUNG*P_1*K_1*(POW/P_2)**POW*EXP(-1.0D0*POW)
C1         ELSE IF ((EF(I) >= E_N_0) .AND. (EF(I) < E_N_HAT)) THEN ! THE FOLLOWING IS THE TYPICAL EXPONENTIAL SOFTENING PART OF THE FIBER PULLOUT LAW
         ELSE IF (EF(I)/K_1 >= E_N_0) THEN ! THE FOLLOWING IS THE TYPICAL EXPONENTIAL SOFTENING PART OF THE FIBER PULLOUT LAW
            SNB_FIB(I) = YOUNG*P_1*K_1*
     $           (EF(I)/K_1 - E_N_0  + E_N_MAX - EPS_SHIFT)**POW*
     $           EXP(-P_2*(EF(I)/K_1 - E_N_0 + E_N_MAX - EPS_SHIFT))  
C1            SNB_FIB(I) = P_1*(POW/P_2)**POW*EXP(-1.0D0*POW)/(E_N_HAT - EF(I))
C1         ELSE IF (EF(I) >= E_N_HAT) THEN
C1            SNB_FIB(I)=0.0D0
         END IF
      END DO
      
C THE RATE FACTOR FOR VOLUMETRIC BOUNDARY BOUNDARY      
      RATE_FACTOR = C_N0*R_N
      SNB_FIB = SNB_FIB * (1.D0 + RATE_FACTOR)
C	  PRINT*,"SNB_FIB",SNB_FIB
      RETURN
      END SUBROUTINE C_NORM_FIB


! +--------------------------------------------------------------------+
! |                      SUBROUTINE C_DEV                              |
! +--------------------------------------------------------------------+
      SUBROUTINE C_DEV(DED, ED0, DEV, EV0, SV0, SD0, SDF, C_D, R_D,
     $     SDNEG, SDPOS) 
      INCLUDE 'VABA_PARAM.INC'
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7FBOUNDS_M7FIO/K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9,
     $     SV0_P, C_1, C_2, C_3, C_4, C_5,
     $     C_6, C_7, C_8, C_9, C_10, C_11,C_12,C_13,C_14,C_15,C_16,
     $     C_17,C_18,C_19,C_20,C_21,C_R0,C_R2
      REAL*8 :: K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9, SV0_P
      REAL*8 :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8,
     $     C_9, C_10
      REAL*8 :: C_11,C_12,C_13,C_14,C_15,C_16,C_17,C_18,
     $     C_19,C_20,C_21,C_R0,C_R2
      COMMON /KM7FBOUNDS_M7FIO_1/ C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1,
     $     C_6_2, C_6_M4, C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4,
     $     C_9_0, C_9_1, C_9_2, C_9_M4
      REAL*8 :: C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1, C_6_2, C_6_M4
      REAL*8 :: C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4
      REAL*8 :: C_9_0, C_9_1, C_9_2, C_9_M4
      COMMON /KM7FBOUNDS_M7FIO_2/SIG_FIB_0, P_1, P_2,P_5,P_6,XP_1,XP_2,
     $     XP_3, V_F, CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3, CF_4,
     $     CF_5, CF_4_4, AP1, BP1, CP2, DP2, EP2, D_1, D_2, D_3, D_4,
     $     D_5, D_6
      REAL*8 :: SIG_FIB_0, P_1, P_2, P_5, P_6
      REAL*8 :: XP_1, XP_2, XP_3, V_F
      REAL*8 :: CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3
      REAL*8 :: CF_4, CF_5, CF_4_4
      REAL*8 :: AP1, BP1, CP2, DP2, EP2
      REAL*8 :: D_1, D_2, D_3, D_4, D_5, D_6
      SAVE :: /KCOMMON_VARS/, /KM7FBOUNDS_M7FIO/, /KM7FBOUNDS_M7FIO_1/,
     $     /KM7FBOUNDS_M7FIO_2/
      REAL*8, DIMENSION(1:NP) :: DED, ED0, SD0
      REAL*8, DIMENSION(1:NP) :: SDF, C_D
C##############
      REAL*8, DIMENSION(1:NP) :: SDNEG, SDPOS
C##############
      REAL*8 :: DEV, EV0, SV0, R_D
      REAL*8 :: USI_C, USI_T, PAR5, PAR6, E_S, CD0, EV
      REAL*8, DIMENSION(1:NP) :: EDF, SDE, SDLOWER, SDUPPER
      REAL*8, DIMENSION(1:NP) :: SDLOWER_FIB, SDUPPER_FIB
      REAL*8 :: RATE_FACTOR_C, RATE_FACTOR_T, C_DC0, C_DT0, EVF, E_V
C      REAL*8 :: CC_N, CC_I, CC_J, PP_N, PPSI, FM 
      REAL*8 :: F_C0, E_0, F_CP, C_40, BETA_15,BETA_16,BETA_17, BETA_18,
     $     BETA_19 
      LOGICAL :: L0,L1,L2,L3,L4
      INTEGER :: J, ISIZE, ALLOCSTAT
      
      ISIZE = SIZE(DED)

      EVF = EV0 + DEV

      EV=YOUNG/(1.0D0 - 2.0D0*POISSON)

      F_C0=15.08D0
      E_0=20000.D0
      C_40=1.0D+0
C      !ACI FORMULA IN MPA. AT LEAST WHEN F_C' IS NOT KNOWN.
C      F_CP=(YOUNG/5150.226D0)**2
      F_CP=50.95
      BETA_15=C_5_1*EXP(-C_40*(F_CP/YOUNG-F_C0/E_0))
      BETA_16=C_8_1*EXP(-C_40*(F_CP/YOUNG-F_C0/E_0))
      BETA_17=C_7_1*EXP(-C_40*(F_CP/YOUNG-F_C0/E_0))
      BETA_18=C_6_0*EXP(-C_40*(F_CP/YOUNG-F_C0/E_0))
      BETA_19=C_9_0*EXP(-C_40*(F_CP/YOUNG-F_C0/E_0))
C      C_5 = (C_5_1*TANH(C_5_0*MAX(-EVF,0.0D0)))*(1 + 25.D0*V_F) + C_5_M4
C      C_5 = (C_5_1*TANH(C_5_0*MAX(-EV0,0.0D0)))*(1.0D0 + CF_1*V_F) + C_5_M4
C      C_5 = (C_5_1*TANH(C_5_0*MAX(-EVF,0.0D0)/K_1))*(1.0D0 + CF_1*TANH(CF_1_1*V_F)) + C_5_M4
      C_5 = (BETA_15*TANH(C_5_0*MAX(-EVF,0.0D0)/K_1))*
     $     (1.0D0 + CF_1*TANH(CF_1_1*V_F))+CF_1*TANH(CF_1_1*V_F)+C_5_M4
C      C_5 = (C_5_1*TANH(C_5_0*MAX(-EVF,0.0D0)/K_1)) + CF_1*TANH(CF_1_1*V_F) + C_5_M4
C      C_5 = (C_5_1*TANH(C_5_0*MAX(-EVF,0.0D0)/K_1)) + CF_1*V_F + C_5_M4
C      C_5 = (C_5_1*TANH(C_5_0*MAX(-SV0/EV,0.0D0)))*(1.0D0 + CF_1*V_F) + C_5_M4
C      C_5 = (1.0D0 + CF_1*V_F)*C_5_M4
C      C_5 = C_5_M4
C      C_8 = (C_8_1*TANH(C_8_0*MAX(-EVF,0.0D0)))*(1 + CF_3*V_F) + C_8_M4
C      C_8 = (C_8_1*TANH(C_8_0*MAX(-EVF,0.0D0)/K_1))*(1.0D0 + CF_3*TANH(CF_3_3*V_F)) + C_8_M4
      C_8 = (BETA_16*TANH(C_8_0*MAX(-EVF,0.0D0)/K_1))*
     $     (1.0D0 + CF_3*TANH(CF_3_3*V_F))+CF_3*TANH(CF_3_3*V_F)+C_8_M4
C      C_8 = (C_8_1*TANH(C_8_0*MAX(-EVF,0.0D0)/K_1)) + CF_3*TANH(CF_3_3*V_F) + C_8_M4
C      C_8 = (C_8_1*TANH(C_8_0*MAX(-EVF,0.0D0)/K_1)) + CF_3*V_F + C_8_M4
C      C_8 = C_8_M4

C      C_7 = (C_7_1*TANH(C_7_0*MAX(-EVF,0.0D0)))*(1 + 340.D0*V_F) + C_7_M4
C      C_7 = C_7_1*TANH(C_7_0*MAX(-EVF,0.0D0)) + C_7_M4 + 24000.D0*V_F
C      C_7 = C_7_1*TANH(C_7_0*MAX(-EV0,0.0D0)) + C_7_M4 + CF_2*V_F
      C_7 = BETA_17*TANH(C_7_0*MAX(-EVF,0.0D0)/K_1)+
     $     CF_2*TANH(CF_2_2*V_F) + C_7_M4 
C      C_7 = C_7_1*TANH(C_7_0*MAX(-EVF-MINVAL(EDF),0.0D0)) + C_7_M4 + CF_2*V_F
C      C_7 = C_7_1*TANH(C_7_0*MAX(-MINVAL(EDF),0.0D0)) + C_7_M4 + CF_2*V_F
C      C_7 = C_7_1*TANH(C_7_0*MAX(-SV0/EV,0.0D0)) + C_7_M4 + CF_2*V_F
C      C_7 = C_7_M4 

C      C_6 = C_6*MIN(EXP(2.D-7*MAX(-EVF-5.D-1*K_1,0.0D0)/K_1),C_6*3.D0)
C      C_9 = C_9*MIN(EXP(2.D-7*MAX(-EVF-5.D-1*K_1,0.0D0)/K_1),C_9*3.D0)

C*      C_6 = 1.3*MIN(EXP(4.D2*MAX(-EVF/K_1-4.D1,0.0D0)),1.3*1.D1)
C*      C_9 = 1.3*MIN(EXP(4.D2*MAX(-EVF/K_1-4.D1,0.0D0)),1.3*1.D1)
C NUMERICAL EXCEPTION IN THE FOLLOWING TWO STATEMENTS WHEN EVF<0. REWRITTEN. 02.07.2011
C      C_6 = C_6_M4*MIN(EXP(BETA_18*MAX(-EVF/K_1-C_6_1,0.0D0)),C_6_2)
C      C_9 = C_9_M4*MIN(EXP(BETA_19*MAX(-EVF/K_1-C_9_1,0.0D0)),C_9_2)
      IF (BETA_18*MAX(-EVF/K_1-C_6_1,0.0D0)>=LOG(C_6_2)) THEN
         C_6 = C_6_M4*C_6_2
      ELSE
         C_6 = C_6_M4*EXP(BETA_18*MAX(-EVF/K_1-C_6_1,0.0D0))
      END IF
      IF (BETA_19*MAX(-EVF/K_1-C_9_1,0.0D0)>=LOG(C_9_2)) THEN
         C_9 = C_9_M4*C_9_2
      ELSE
         C_9 = C_9_M4*EXP(BETA_19*MAX(-EVF/K_1-C_9_1,0.0D0))
      END IF
      
C      C_6 = 1.3*(2.0D0 + 1.0D0*TANH(1.D2*MAX(-EVF/K_1,0.0D0) - 1.D0))
C      C_9 = 1.3*(2.0D0 + 1.0D0*TANH(1.D2*MAX(-EVF/K_1,0.0D0) - 1.D0))
C      C_6 = 1.3*(2.0D0 + 1.0D0*TANH(1.D2*MAX(-SV0/E_V/K_1,0.0D0) - 1.D0))
C      C_9 = 1.3*(2.0D0 + 1.0D0*TANH(1.D2*MAX(-SV0/E_V/K_1,0.0D0) - 1.D0))
C      WRITE(*,*) 'C_6=',C_6

C*      C_20 = 0.232D0*EXP(-440.0D0*V_F*MAX(SIGN(1.0D0,-SV0),0D0)) + 0.05D0 !0.0679D0
C      C_20 = C_20_1*EXP(-C_20_0*V_F*MAX(SIGN(1.0D0,-SV0),0D0)) + C_20_2 !0.0679D0

C C_20 IS NOW USED AS THE LOADING/UNLOADING SLOPE PARAMETER.09.08.2011
C*      C_20 = CF_4*EXP(-CF_4_4*V_F*MAX(SIGN(1.0D0,-SV0),0D0)) + CF_5 !0.0679D0

      PAR5 = K_1*YOUNG*C_8
      PAR6 = YOUNG*C_5*K_1
      CD0=YOUNG/(1.D0+POISSON)*(1.D0-4.D0*POISSON)/(1.D0-2.D0*POISSON)
      E_V  = YOUNG/(1.D0-2*POISSON)
C$      USI_C = C_19
C$      USI_T = C_21 !THE UNLOADING SLOPE INTERPOLATOR IN TENSILE DEV. B.

C COMPRESSIVE DEVIATORIC BOUNDARY RATE PARAMETERS       
      C_DC0  = C_R2

C TENSILE DEVIATORIC BOUNDARY RATE PARAMETERS      
      C_DT0  = C_R2

C#######################
      C_D=CD0
C#######################
      EDF = ED0+DED 
C C0 IS THE CREEP COEFFICIENT, FOR NO CREEP C0=1
C$      SDE = SD0*C0+C_D*DED 

C EQUATION (12D)      
C      CALL DEV_COMP(EDF, EVF, SV0, SDLOWER) 
      CALL DEV_COMP(EDF, SDLOWER) 
C      CALL DEV_COMP_FIB(EDF, SDLOWER_FIB) 
C EQUATION (12E)      
C      CALL DEV_TENS(EDF, EVF, SV0, SDUPPER)
      CALL DEV_TENS(EDF, SDUPPER)
      CALL DEV_TENS_FIB(EDF, SDUPPER_FIB)

C! DEFINE THE MODIFYING FUNCTION HERE AND REDUCE THE BOUNDARIES. 14.09.2010
C      CC_N = 0.0D0  ! OR ?
C      CC_I = 5.0D-6 ! OR 1.0D-3
C      CC_J = 3.0D+3 ! OR 1.0D+3
C      PP_N = 1.0D0 ! OR 2.0D0
C      PPSI = MAXVAL(ED0+DED+EV0+DEV,1) - MINVAL(ED0+DED+EV0+DEV,1) ! \EPS_I - \EPS_III
C      FM = CC_N + (1.0D0 - CC_N)/(1.0D0 + MAX((PPSI-CC_I)/K_1/CC_J,0.0D0)**PP_N)
C! CHOOSE AMONG THE FOLLOWING POSSIBILITIES (MAY BE ALL OF THEM), 14.09.2010:      
C      SDLOWER = SDLOWER*FM
C      SDUPPER = SDUPPER*FM
C      SDUPPER_FIB = SDUPPER_FIB*FM

C THE RATE FACTOR FOR COMPRESSIVE DEVIATORIC BOUNDARY      
      RATE_FACTOR_C = C_DC0*R_D
      RATE_FACTOR_T = C_DT0*R_D
C THIS ONE IS AFTER MARCH 21:
      SDLOWER = SDLOWER * (1.D0 + RATE_FACTOR_C)      
C THIS ONE IS AFTER MARCH 21:
      SDUPPER = SDUPPER * (1.D0 + RATE_FACTOR_T)
C################
      SDNEG=SDLOWER
      SDPOS=SDUPPER
C################
      
C EQUATION (15A)      
C      SDF = MIN(MAX(SDE,SDLOWER+SDLOWER_FIB),SDUPPER+SDUPPER_FIB)
C$      SDF = MIN(MAX(SDE,SDLOWER),SDUPPER+SDUPPER_FIB)
      RETURN
      END SUBROUTINE C_DEV

! +--------------------------------------------------------------------+
! |                        SUBROUTINE DEV_COMP                         |
! +--------------------------------------------------------------------+
!      SUBROUTINE DEV_COMP(EDF, EVF, SV0, SDLOWER) 
      SUBROUTINE DEV_COMP(EDF, SDLOWER) 
      INCLUDE 'VABA_PARAM.INC'
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7FBOUNDS_M7FIO/K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9,
     $     SV0_P, C_1, C_2, C_3, C_4, C_5,
     $     C_6, C_7, C_8, C_9, C_10, C_11,C_12,C_13,C_14,C_15,C_16,
     $     C_17,C_18,C_19,C_20,C_21,C_R0,C_R2
      REAL*8 :: K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9, SV0_P
      REAL*8 :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8,
     $     C_9, C_10
      REAL*8 :: C_11,C_12,C_13,C_14,C_15,C_16,C_17,C_18,
     $     C_19,C_20,C_21,C_R0,C_R2
      COMMON /KM7FBOUNDS_M7FIO_1/ C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1,
     $     C_6_2, C_6_M4, C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4,
     $     C_9_0, C_9_1, C_9_2, C_9_M4
      REAL*8 :: C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1, C_6_2, C_6_M4
      REAL*8 :: C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4
      REAL*8 :: C_9_0, C_9_1, C_9_2, C_9_M4
      COMMON /KM7FBOUNDS_M7FIO_2/SIG_FIB_0, P_1, P_2, P_5,P_6,XP_1,XP_2,
     $     XP_3, V_F, CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3, CF_4,
     $     CF_5, CF_4_4, AP1, BP1, CP2, DP2, EP2, D_1, D_2, D_3, D_4,
     $     D_5, D_6
      REAL*8 :: SIG_FIB_0, P_1, P_2, P_5, P_6
      REAL*8 :: XP_1, XP_2, XP_3, V_F
      REAL*8 :: CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3
      REAL*8 :: CF_4, CF_5, CF_4_4
      REAL*8 :: AP1, BP1, CP2, DP2, EP2
      REAL*8 :: D_1, D_2, D_3, D_4, D_5, D_6
      SAVE :: /KCOMMON_VARS/, /KM7FBOUNDS_M7FIO/, /KM7FBOUNDS_M7FIO_1/,
     $     /KM7FBOUNDS_M7FIO_2/
      REAL*8, DIMENSION(1:NP) :: EDF
C      REAL*8 :: EVF, SV0
      REAL*8, DIMENSION(1:NP) :: SDLOWER
      REAL*8 :: PAR1, PAR2, PAR5

      PAR1 = C_9*C_8*K_1
      PAR2 = C_7*K_1 
      PAR5 = K_1*YOUNG*C_8

      SDLOWER=-PAR5/(1.D0+(MAX(-EDF-PAR1,0.0D0)/PAR2)**2.0D0) 

      RETURN 
      END SUBROUTINE DEV_COMP 

C! +--------------------------------------------------------------------+
C! |                        SUBROUTINE DEV_COMP_FIB                     |
C! +--------------------------------------------------------------------+
C      SUBROUTINE DEV_COMP_FIB(EDF,SDLOWER_FIB) 
C      INCLUDE 'VABA_PARAM.INC'
C      REAL*8 :: EDF
C      REAL*8 :: SDLOWER_FIB
C      REAL*8 :: PAR1, PAR2, PAR5
C
C      PAR1 = C_9*C_8*K_1
C      PAR2 = C_7*K_1 
C      PAR5 = K_1*YOUNG*C_8
C
C      SDLOWER_FIB=-PAR5/(1.D0+(MAX(-EDF+PAR1,0.0D0)/PAR2)**2.0D0) 
C
C      RETURN 
C      END SUBROUTINE DEV_COMP_FIB

C +--------------------------------------------------------------------+
C |                         SUBROUTINE DEV_TENS                        |
C +--------------------------------------------------------------------+
C      SUBROUTINE DEV_TENS(EDF, EVF, SV0, SDUPPER) 
      SUBROUTINE DEV_TENS(EDF, SDUPPER) 
      INCLUDE 'VABA_PARAM.INC'
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7FBOUNDS_M7FIO/K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9,
     $     SV0_P, C_1, C_2, C_3, C_4, C_5,
     $     C_6, C_7, C_8, C_9, C_10, C_11,C_12,C_13,C_14,C_15,C_16,
     $     C_17,C_18,C_19,C_20,C_21,C_R0,C_R2
      REAL*8 :: K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9, SV0_P
      REAL*8 :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8,
     $     C_9, C_10
      REAL*8 :: C_11,C_12,C_13,C_14,C_15,C_16,C_17,C_18,
     $     C_19,C_20,C_21,C_R0,C_R2
      COMMON /KM7FBOUNDS_M7FIO_1/ C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1,
     $     C_6_2, C_6_M4, C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4,
     $     C_9_0, C_9_1, C_9_2, C_9_M4
      REAL*8 :: C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1, C_6_2, C_6_M4
      REAL*8 :: C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4
      REAL*8 :: C_9_0, C_9_1, C_9_2, C_9_M4
      COMMON /KM7FBOUNDS_M7FIO_2/SIG_FIB_0, P_1, P_2, P_5,P_6,XP_1,XP_2,
     $     XP_3, V_F, CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3, CF_4,
     $     CF_5, CF_4_4, AP1, BP1, CP2, DP2, EP2, D_1, D_2, D_3, D_4,
     $     D_5, D_6
      REAL*8 :: SIG_FIB_0, P_1, P_2, P_5, P_6
      REAL*8 :: XP_1, XP_2, XP_3, V_F
      REAL*8 :: CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3
      REAL*8 :: CF_4, CF_5, CF_4_4
      REAL*8 :: AP1, BP1, CP2, DP2, EP2
      REAL*8 :: D_1, D_2, D_3, D_4, D_5, D_6
      SAVE :: /KCOMMON_VARS/, /KM7FBOUNDS_M7FIO/, /KM7FBOUNDS_M7FIO_1/,
     $     /KM7FBOUNDS_M7FIO_2/
      REAL*8, DIMENSION(1:NP) :: EDF
C      REAL*8 :: EVF, SV0
      REAL*8, DIMENSION(1:NP) :: SDUPPER
      REAL*8 :: PAR2, PAR3, PAR4, PAR6

      PAR2 = C_7*K_1 
      PAR3 = C_20
      PAR4 = C_6*C_5*K_1
      PAR6 = YOUNG*C_5*K_1

      SDUPPER=PAR6/(1.0D0+(MAX(EDF-PAR4,0.0D0)/PAR2/PAR3)**2.0D0)

      RETURN 
      END SUBROUTINE DEV_TENS 

! +--------------------------------------------------------------------+
! |                         SUBROUTINE DEV_TENS_FIB                    |
! +--------------------------------------------------------------------+
      SUBROUTINE DEV_TENS_FIB(EDF, SDUPPER_FIB) 
      INCLUDE 'VABA_PARAM.INC'
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7FBOUNDS_M7FIO/K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9,
     $     SV0_P, C_1, C_2, C_3, C_4, C_5,
     $     C_6, C_7, C_8, C_9, C_10, C_11,C_12,C_13,C_14,C_15,C_16,
     $     C_17,C_18,C_19,C_20,C_21,C_R0,C_R2
      REAL*8 :: K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9, SV0_P
      REAL*8 :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8,
     $     C_9, C_10
      REAL*8 :: C_11,C_12,C_13,C_14,C_15,C_16,C_17,C_18,
     $     C_19,C_20,C_21,C_R0,C_R2
      COMMON /KM7FBOUNDS_M7FIO_1/ C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1,
     $     C_6_2, C_6_M4, C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4,
     $     C_9_0, C_9_1, C_9_2, C_9_M4
      REAL*8 :: C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1, C_6_2, C_6_M4
      REAL*8 :: C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4
      REAL*8 :: C_9_0, C_9_1, C_9_2, C_9_M4
      COMMON /KM7FBOUNDS_M7FIO_2/SIG_FIB_0, P_1, P_2, P_5,P_6,XP_1,XP_2,
     $     XP_3, V_F, CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3, CF_4,
     $     CF_5, CF_4_4, AP1, BP1, CP2, DP2, EP2, D_1, D_2, D_3, D_4,
     $     D_5, D_6
      REAL*8 :: SIG_FIB_0, P_1, P_2, P_5, P_6
      REAL*8 :: XP_1, XP_2, XP_3, V_F
      REAL*8 :: CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3
      REAL*8 :: CF_4, CF_5, CF_4_4
      REAL*8 :: AP1, BP1, CP2, DP2, EP2
      REAL*8 :: D_1, D_2, D_3, D_4, D_5, D_6
      SAVE :: /KCOMMON_VARS/, /KM7FBOUNDS_M7FIO/, /KM7FBOUNDS_M7FIO_1/,
     $     /KM7FBOUNDS_M7FIO_2/
      REAL*8, DIMENSION(1:NP) :: EDF
      REAL*8, DIMENSION(1:NP) :: SDUPPER_FIB
C !     REAL*8 :: PAR2, PAR3, PAR4, PAR6

C =========================================
C K_1 (RADIAL SCALING) IN FIBER LAW DONE.
C E   (VERTICAL SCALING) IN FIBER LAW DONE.
C =========================================

      SDUPPER_FIB=YOUNG*P_5*K_1*MAX(EDF/K_1,0.0D0)*
     $     EXP(-P_6*MAX(EDF/K_1,0.0D0))

      RETURN 
      END SUBROUTINE DEV_TENS_FIB 

      SUBROUTINE C_SHEAR_TENS(SIG_O, S_0, FSP_0)
      INCLUDE 'VABA_PARAM.INC'
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7FBOUNDS_M7FIO/K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9,
     $     SV0_P, C_1, C_2, C_3, C_4, C_5,
     $     C_6, C_7, C_8, C_9, C_10, C_11,C_12,C_13,C_14,C_15,C_16,
     $     C_17,C_18,C_19,C_20,C_21,C_R0,C_R2
      REAL*8 :: K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9, SV0_P
      REAL*8 :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8,
     $     C_9, C_10
      REAL*8 :: C_11,C_12,C_13,C_14,C_15,C_16,C_17,C_18,
     $     C_19,C_20,C_21,C_R0,C_R2
      COMMON /KM7FBOUNDS_M7FIO_1/ C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1,
     $     C_6_2, C_6_M4, C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4,
     $     C_9_0, C_9_1, C_9_2, C_9_M4
      REAL*8 :: C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1, C_6_2, C_6_M4
      REAL*8 :: C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4
      REAL*8 :: C_9_0, C_9_1, C_9_2, C_9_M4
      COMMON /KM7FBOUNDS_M7FIO_2/SIG_FIB_0, P_1, P_2, P_5,P_6,XP_1,XP_2,
     $     XP_3, V_F, CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3, CF_4,
     $     CF_5, CF_4_4, AP1, BP1, CP2, DP2, EP2, D_1, D_2, D_3, D_4,
     $     D_5, D_6
      REAL*8 :: SIG_FIB_0, P_1, P_2, P_5, P_6
      REAL*8 :: XP_1, XP_2, XP_3, V_F
      REAL*8 :: CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3
      REAL*8 :: CF_4, CF_5, CF_4_4
      REAL*8 :: AP1, BP1, CP2, DP2, EP2
      REAL*8 :: D_1, D_2, D_3, D_4, D_5, D_6
      SAVE :: /KCOMMON_VARS/, /KM7FBOUNDS_M7FIO/, /KM7FBOUNDS_M7FIO_1/,
     $     /KM7FBOUNDS_M7FIO_2/
      REAL*8 :: SIG_O, S_0
      REAL*8 :: FSP_0
      REAL*8 :: C_6_P

      C_6_P = C_10
      FSP_0 = C_6_P * MAX(SIG_O, 0.0D0) / (1.D0 + C_6_P / S_0 *
     $     MAX(SIG_O, 0.0D0)) 
      RETURN
      END SUBROUTINE C_SHEAR_TENS

! +--------------------------------------------------------------------+
! |                          SUBROUTINE C_SHEAR2                       |
! +--------------------------------------------------------------------+
      SUBROUTINE C_SHEAR2(EPS_L, EPS_M, C_D, SNF, R_S, DEL, DEM, SL0,
     $     SM0, SLF, SMF, EPS_V) 
      INCLUDE 'VABA_PARAM.INC'
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7FBOUNDS_M7FIO/K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9,
     $     SV0_P, C_1, C_2, C_3, C_4, C_5,
     $     C_6, C_7, C_8, C_9, C_10, C_11,C_12,C_13,C_14,C_15,C_16,
     $     C_17,C_18,C_19,C_20,C_21,C_R0,C_R2
      REAL*8 :: K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9, SV0_P
      REAL*8 :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8,
     $     C_9, C_10
      REAL*8 :: C_11,C_12,C_13,C_14,C_15,C_16,C_17,C_18,
     $     C_19,C_20,C_21,C_R0,C_R2
      COMMON /KM7FBOUNDS_M7FIO_1/ C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1,
     $     C_6_2, C_6_M4, C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4,
     $     C_9_0, C_9_1, C_9_2, C_9_M4
      REAL*8 :: C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1, C_6_2, C_6_M4
      REAL*8 :: C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4
      REAL*8 :: C_9_0, C_9_1, C_9_2, C_9_M4
      COMMON /KM7FBOUNDS_M7FIO_2/SIG_FIB_0, P_1, P_2, P_5,P_6,XP_1,XP_2,
     $     XP_3, V_F, CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3, CF_4,
     $     CF_5, CF_4_4, AP1, BP1, CP2, DP2, EP2, D_1, D_2, D_3, D_4,
     $     D_5, D_6
      REAL*8 :: SIG_FIB_0, P_1, P_2, P_5, P_6
      REAL*8 :: XP_1, XP_2, XP_3, V_F
      REAL*8 :: CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3
      REAL*8 :: CF_4, CF_5, CF_4_4
      REAL*8 :: AP1, BP1, CP2, DP2, EP2
      REAL*8 :: D_1, D_2, D_3, D_4, D_5, D_6
      SAVE :: /KCOMMON_VARS/, /KM7FBOUNDS_M7FIO/, /KM7FBOUNDS_M7FIO_1/,
     $     /KM7FBOUNDS_M7FIO_2/
      REAL*8, DIMENSION(1:NP)  :: EPS_L, EPS_M, C_D, SNF, DEL,
     $     DEM, SL0, SM0 
      REAL*8  :: R_S, EPS_V
      REAL*8, DIMENSION(1:NP) :: SLF, SMF  
      REAL*8 :: E_T, C_12_P, C_13_P, C_6_P, S_0
      REAL*8, DIMENSION(1:NP) :: C_T, FSP
      REAL*8, DIMENSION(1:NP) :: STF, FSP_0, EPSF, DEPSF
      REAL*8 :: RATE_FACTOR_1, C_S0, SIG_O
      INTEGER :: ISIZE, ALLOCSTAT, J
      
      ISIZE = SIZE(SNF)
C     ----------------------
C ... PARAMETER ASSIGNEMENTS
C     ----------------------

      E_T=YOUNG/(1.D0+POISSON)*(1.D0-4.D0*POISSON)/(1.D0-2.D0*POISSON) 
      C_6_P  = C_10
      C_12_P = C_11 
      C_13_P = C_12*3.371D-4
C      C_13_P = C_12*1.D1
      S_0    = K_1*K_2*E_T
       
C RATE EFFECT PARAMETER      
      C_S0 = C_R2

C THE RATE FACTOR FOR FRICTION BOUNDARY      
      RATE_FACTOR_1 = C_S0*R_S  ! RATE EFFECT ON THE ASYMPTOTE

C THE RATE EFFECT !NEGATIVE AND POSITIVE BOUNDARIES ARE SYMMETRIC!      
C NEW RATE EFFECT ON THE SHEAR BOUNDARY 23.08.2011
C      S_0 = S_0 * (1.D0 + RATE_FACTOR_1)  ! THE NEW FORMULA: ASYMPOTE IS 
                                          ! RATE DEPENDENT.
C END NEW RATE EFFECT ON THE SHEAR BOUNDARY 23.08.2011
C THE FRICTION LAW (NOTE: EXP(-C_13*ABS(EPS_V)) IS APPROX. EQUAL 
C TO (1+C_13*ABS(EPS_V))**(-1)
C      SIG_O = E_T * K_1 * C_12_P / (1.0D0 + C_13_P * MAX(EPS_V, 0.0D0))
C 21/04/2010:
C      SIG_O = MAX(E_T * K_1 * C_12_P - C_13_P*MAX(EPS_V, 0.0D0)/K_1,0.0D0)
C 25/07/2010:
      SIG_O = MAX(E_T * K_1 *
     $     (C_12_P - C_13_P*MAX(EPS_V, 0.0D0)/K_1),0.0D0)
C      SIG_O = MAX(E_T * K_1 * C_12_P - C_13_P*MAX(EPS_V, 0.0D0),0.0D0)
C NEW RATE EFFECT ON THE SHEAR BOUNDARY 23.08.2011
      FSP = C_6_P * MAX(-SNF + SIG_O, 0.0D0) / (1.D0 + C_6_P / S_0 *
     $     MAX(-SNF + SIG_O, 0.0D0)) * (1.D0 + RATE_FACTOR_1)
C END NEW RATE EFFECT ON THE SHEAR BOUNDARY 23.08.2011
C      FSP_0 = C_6_P * MAX(SIG_O, 0.0D0) / (1.D0 + C_6_P / S_0 * &
C           &MAX(SIG_O, 0.0D0))
C!      FSP_0 = 1.D-1*C_6_P * E_T * K_1 * C_12_P / (1.D0 + C_6_P / S_0 * &
C!           &E_T * K_1 * C_12_P)
      EPSF = SQRT(EPS_L*EPS_L + EPS_M*EPS_M)
      DEPSF = SQRT(DEL*DEL + DEM*DEM)
      DO J=1, ISIZE
C         CALL C_SHEAR_TENS(EPSF(J), DEPSF(J), SIG_O, S_0, FSP_0(J))
         CALL C_SHEAR_TENS(SIG_O, S_0, FSP_0(J))
C NEW RATE EFFECT ON THE SHEAR BOUNDARY 23.08.2011
         FSP_0(J)= FSP_0(J) * (1.D0 + RATE_FACTOR_1)
C END NEW RATE EFFECT ON THE SHEAR BOUNDARY 23.08.2011
      END DO

C UNLOADING IN FRICTION BOUNDARY
      DO J=1,ISIZE
        IF(SL0(J)*DEL(J) < 0.0D0 .OR. SM0(J)*DEM(J) < 0.0D0) THEN
          C_T(J) = C_D(J)
        ELSE                 !ALL OTHER CASES
          C_T(J) = E_T
        END IF  
      END DO

C C0 IS THE CREEP COEFFICIENT, FOR NO CREEP C0=1
      SLF = SL0*C0 + C_T*DEL
      SMF = SM0*C0 + C_T*DEM
      STF = SQRT((SLF*SLF + SMF*SMF))
      DO J = 1, ISIZE
C         IF (SNF(J) < 0.0) THEN ! FRICTION BOUNDARY IS ENFORCED ONLY WHEN NORMAL STRESS IS COMPRESSIVE
            IF (STF(J) .GT.1.0D-10) THEN ! PREVENT DIVISION BY ZERO.
               SLF(J) = SLF(J) / STF(J)
               SMF(J) = SMF(J) / STF(J)
            END IF
C         ELSE
C            STF(J) = FSP_0
C            SLF(J) = SLF(J) / STF(J)
C            SMF(J) = SMF(J) / STF(J)           
C         END IF
      END DO

C EQUATIONS (18A,B), (19A,B) 
      DO J=1, ISIZE
         IF (SNF(J) < 0.0 ) THEN ! FRICTION BOUNDARY IS ENFORCED ONLY WHEN NORMAL STRESS IS COMPRESSIVE
            IF (STF(J) > FSP(J)) THEN
               STF(J) = FSP(J)
            END IF
         ELSE
            IF (STF(J) > FSP_0(J)) THEN
               STF(J) = FSP_0(J)
            END IF
C            IF (STF(J) > FSP(J)) THEN
C               STF(J) = FSP(J)
C            END IF
         END IF
         SLF(J) = STF(J) * SLF(J)
         SMF(J) = STF(J) * SMF(J)
      END DO
C      SLF = STF * SLF
C      SMF = STF * SMF

      RETURN 
      END SUBROUTINE C_SHEAR2 

C +--------------------------------------------------------------------+
C |                        SUBROUTINE C_VOL                            |
C +--------------------------------------------------------------------+
      SUBROUTINE C_VOL(DEV, EV0, SV0, DEPS_N, EPS_N, R_N, SVNEG) 
      INCLUDE 'VABA_PARAM.INC'
      EXTERNAL FB
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7FBOUNDS_M7FIO/K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9,
     $     SV0_P, C_1, C_2, C_3, C_4, C_5,
     $     C_6, C_7, C_8, C_9, C_10, C_11,C_12,C_13,C_14,C_15,C_16,
     $     C_17,C_18,C_19,C_20,C_21,C_R0,C_R2
      REAL*8 :: K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9, SV0_P
      REAL*8 :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8,
     $     C_9, C_10
      REAL*8 :: C_11,C_12,C_13,C_14,C_15,C_16,C_17,C_18,
     $     C_19,C_20,C_21,C_R0,C_R2
      COMMON /KM7FBOUNDS_M7FIO_1/ C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1,
     $     C_6_2, C_6_M4, C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4,
     $     C_9_0, C_9_1, C_9_2, C_9_M4
      REAL*8 :: C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1, C_6_2, C_6_M4
      REAL*8 :: C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4
      REAL*8 :: C_9_0, C_9_1, C_9_2, C_9_M4
      COMMON /KM7FBOUNDS_M7FIO_2/SIG_FIB_0, P_1, P_2, P_5,P_6,XP_1,XP_2,
     $     XP_3, V_F, CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3, CF_4,
     $     CF_5, CF_4_4, AP1, BP1, CP2, DP2, EP2, D_1, D_2, D_3, D_4,
     $     D_5, D_6
      REAL*8 :: SIG_FIB_0, P_1, P_2, P_5, P_6
      REAL*8 :: XP_1, XP_2, XP_3, V_F
      REAL*8 :: CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3
      REAL*8 :: CF_4, CF_5, CF_4_4
      REAL*8 :: AP1, BP1, CP2, DP2, EP2
      REAL*8 :: D_1, D_2, D_3, D_4, D_5, D_6
      SAVE :: /KCOMMON_VARS/, /KM7FBOUNDS_M7FIO/, /KM7FBOUNDS_M7FIO_1/,
     $     /KM7FBOUNDS_M7FIO_2/
      REAL*8 :: DEV, EV0, SV0
      REAL*8, DIMENSION(1:NP) :: EPS_N, DEPS_N
C      REAL*8 :: DFB, CV 
      REAL*8 :: SVNEG 
      REAL*8 :: R_N 
      REAL*8 :: CV0, EVF, SVB_NEG !, SVE, C_10_P, C_11_P 
      REAL*8 :: XK0, E0_V
C      REAL*8 :: DFB, CV
C      LOGICAL :: L1, L2
      REAL*8 :: PRSTRAINDIFF, XK4 !, K_ST , K_8, K_9, SV0_P
      
      XK0  = K_3*K_1*YOUNG 
      E0_V = K_4*K_1 
      CV0  = YOUNG / (1.0D0 - 2.0D0 * POISSON) !=E_V

      PRSTRAINDIFF = MAXVAL(EPS_N+DEPS_N,1) - MINVAL(EPS_N+DEPS_N,1)

      XK4 = (K_6*(PRSTRAINDIFF/K_1)**K_7)
     $     /(1.0D0 + MIN(MAX(-SV0,0.D0),SV0_P)/CV0) + K_4

C      XK4 = K_4
      E0_V = XK4*K_1

      EVF = EV0+DEV 

C$      IF ( MAX(-EVF/K_1,0.0D0) < K_5) THEN ! DON'T USE K_5, USE SV0_P
C$         CV = FB(1,EV0,XK0,E0_V) 
C$      END IF
C$      CV = FB(1,EV0,XK0,E0_V) 
C$      SVB_NEG = DFB + CV/(1.0D0+K_6*PRSTRAINDIFF**K_7)*DEV
C$      SVB_NEG = DFB + CV*DEV
C$      DFB=SVB_NEG

      SVB_NEG = FB(0,EVF,XK0,E0_V) 
      SVNEG=SVB_NEG*(1.0D0 + C_R2*R_N)

      END SUBROUTINE C_VOL
 
C +--------------------------------------------------------------------+
C |                            FUNCTION FB                             |    
C +--------------------------------------------------------------------+
      FUNCTION FB(I,EF,XK0,E0_V) 
      INCLUDE 'VABA_PARAM.INC'
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7FBOUNDS_M7FIO/K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9,
     $     SV0_P, C_1, C_2, C_3, C_4, C_5,
     $     C_6, C_7, C_8, C_9, C_10, C_11,C_12,C_13,C_14,C_15,C_16,
     $     C_17,C_18,C_19,C_20,C_21,C_R0,C_R2
      REAL*8 :: K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9, SV0_P
      REAL*8 :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8,
     $     C_9, C_10
      REAL*8 :: C_11,C_12,C_13,C_14,C_15,C_16,C_17,C_18,
     $     C_19,C_20,C_21,C_R0,C_R2
      COMMON /KM7FBOUNDS_M7FIO_1/ C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1,
     $     C_6_2, C_6_M4, C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4,
     $     C_9_0, C_9_1, C_9_2, C_9_M4
      REAL*8 :: C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1, C_6_2, C_6_M4
      REAL*8 :: C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4
      REAL*8 :: C_9_0, C_9_1, C_9_2, C_9_M4
      COMMON /KM7FBOUNDS_M7FIO_2/SIG_FIB_0, P_1, P_2, P_5,P_6,XP_1,XP_2,
     $     XP_3, V_F, CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3, CF_4,
     $     CF_5, CF_4_4, AP1, BP1, CP2, DP2, EP2, D_1, D_2, D_3, D_4,
     $     D_5, D_6
      REAL*8 :: SIG_FIB_0, P_1, P_2, P_5, P_6
      REAL*8 :: XP_1, XP_2, XP_3, V_F
      REAL*8 :: CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3
      REAL*8 :: CF_4, CF_5, CF_4_4
      REAL*8 :: AP1, BP1, CP2, DP2, EP2
      REAL*8 :: D_1, D_2, D_3, D_4, D_5, D_6
      SAVE :: /KCOMMON_VARS/, /KM7FBOUNDS_M7FIO/, /KM7FBOUNDS_M7FIO_1/,
     $     /KM7FBOUNDS_M7FIO_2/
C THIS FUNCTION IS FOR THE VOLUMETRIC COMPRESSIVE BOUNDARY.
      INTEGER ::  I
      REAL*8 :: EF,XK0,E0_V
      REAL*8 :: FB, AUX, AUX2
      
      AUX=-XK0*EXP(-EF/E0_V)/0.85D0
      IF (I.EQ.0) FB = AUX 
C SLOPE OF THE BOUNDARY      
      IF (I.EQ.1) FB = -AUX/E0_V 
      RETURN
      END FUNCTION FB

C **********************************************************************
C *** SUBROUTINE SETSYSTEM *********************************************
C **********************************************************************
 
      SUBROUTINE SETSYSTEM() 
      INCLUDE 'VABA_PARAM.INC'
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7F_MICROPLANE_1/ QN, QL, QM, W 
      COMMON /KM7F_MICROPLANE_2/ DEPS_N, DEPS_L, DEPS_M, EPS_N,  EPS_L,
     $     EPS_M,  SIG_N
      COMMON /KM7F_MICROPLANE_3/ SD0, SN0, DED, ED0, EDF, ENF, SNF, SNB,
     $     SDF
      COMMON /KM7F_MICROPLANE_4/ SL0, SM0, SLF, SMF, DEL, DEM, VH_INI,
     $     VH_FIN, C_D, E_N
      COMMON /KM7F_MICROPLANE_5/ SDNEG, SDPOS, EPS_N0_NEG, EPS_N0_POS,
     $     SNB_FIB, CV
C      INTEGER      :: NVHI
C THESE ARE THE DYNAMIC SYSTEM VARIABLES. THEY ARE ALLOCATED AND INITIALIZED 
C IN SUBROUTINE ALLOC. THEY MUST BE DEALLOCATED AFTER THE LAST CALL. 
      REAL*8, DIMENSION(1:6,1:NP) :: QN, QL, QM
      REAL*8, DIMENSION(1:NP)     :: W
      REAL*8, DIMENSION(1:NP)     :: DEPS_N,
     $     DEPS_L, DEPS_M, EPS_N,  EPS_L,  EPS_M,  SIG_N
      REAL*8, DIMENSION(1:NP)     :: SD0, SN0,
     $     DED, ED0, EDF
      REAL*8, DIMENSION(1:NP)     :: ENF, SNF,
     $     SNB, SDF, SDNEG, SDPOS, SNB_FIB
      REAL*8, DIMENSION(1:NP)     :: SL0, SM0, SLF,
     $     SMF, DEL, DEM, EPS_N0_NEG, EPS_N0_POS
C      REAL*8, DIMENSION(1:NVHI)    :: VH_INI, VH_FIN 
      REAL*8, DIMENSION(1:NVHI+2)    :: VH_INI, VH_FIN ! ADDITIONAL HIST. VAR. FOR ELEMENT DELETION 
      REAL*8, DIMENSION(1:NP)     :: C_D, E_N
      REAL*8 :: CV
      SAVE    :: /KCOMMON_VARS/
      SAVE    :: /KM7F_MICROPLANE_1/,/KM7F_MICROPLANE_2/,
     $     /KM7F_MICROPLANE_3/,/KM7F_MICROPLANE_4/,/KM7F_MICROPLANE_5/
      INTEGER :: JP, IJ(1:2,1:6), I, J, K, ALLOCSTAT
      INTEGER :: RS_SIZE
      REAL*8, DIMENSION(1:4,1:NP) :: TE
      REAL*8, DIMENSION(1:3) :: XN, XM, XL, RAND_VEC 
      REAL*8 :: LENGTHN, LENGTHM, LENGTHL

       COMMON /EIG_VAL_VEC/PR_SIG_N1,PR_SIG_N2,PR_SIG_N3,PR_SIG_N1_T,
     $   PR_SIG_N2_T,PR_SIG_N3_T,STRESS_PR,STRESS_PR1,SIG_NEW

C      NP=37
C NVHM IS THE NUMBER OF HISTORY VARIABLES PER MICROPLANE TO STORE 
C NVHF IS THE HISTORY VARIABLES COMMON TO ALL MICROPLANES
C      NVHM = 5; NVHF=2
C      NVHI = NVHM*NP+NVHF
C
C HISTORY VARIABLES:
C 
C VH_INI(1)=VOLUMETRIC STRESS
C VH_INI(2:NVHI:NVHM)=MICROPLANE NORMAL STRESS
C VH_INI(3:NVHI:NVHM)=MICROPLANE L-DIR SHEAR STRESS
C VH_INI(4:NVHI:NVHM)=MICROPLANE M-DIR SHEAR STRESS
C VH_INI(5:NVHI:NVHM)=MAX STRAIN THAT CAN CHARACTERIZE TENSION DAMAGE
C VH_INI(6:NVHI:NVHM)=MAX STRAIN THAT CAN CHARACTERIZE COMPRESSION DAMAGE

C +---------------------------------------------------------------------+
C |ALLOCATE THE VECTORS AND TENSORS. DO NOT FORGET DEALLOCATING THEM IN |
C |THE FE DRIVER BY CALLING SUBROUTINE DONE() WHEN ANALYSIS ENDS.       |
C +---------------------------------------------------------------------+
C         ALLOCATE(VH_INI(1:NVHI), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C            WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VH_INI (MICROPLANE.F90)'
C            STOP
C         END IF 
         VH_INI = 0.0D0  
C INITIALIZED IN THE SUBROUTINE M7FMATERIAL.
C         VH_INI(2)=1.0D0
C         ALLOCATE(VH_FIN(1:NVHI), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C            WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VH_FIN (MICROPLANE.F90)'
C            STOP
C         END IF 
         VH_FIN = 0.0D0  
C         ALLOCATE(QN(1:6,1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE TENSOR QN (MICROPLANE.F90)'
C           STOP
C         END IF
         QN = 0.0D0
C         ALLOCATE(QL(1:6,1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE TENSOR QL (MICROPLANE.F90)'
C           STOP
C         END IF
         QL = 0.0D0
C         ALLOCATE(QM(1:6,1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE TENSOR QM (MICROPLANE.F90)'
C           STOP
C         END IF
         QM = 0.0D0
C         ALLOCATE(W(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE WEIGHTS W (MICROPLANE.F90)'
C           STOP
C         END IF
         W = 0.0D0
C         ALLOCATE(DEPS_N(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR DEPS_N (MICROPLANE.F90)'
C           STOP
C         END IF
         DEPS_N = 0.0D0
C         ALLOCATE(DEPS_M(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR DEPS_M (MICROPLANE.F90)'
C           STOP
C         END IF
         DEPS_M = 0.0D0
C         ALLOCATE(DEPS_L(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR DEPS_L (MICROPLANE.F90)'
C           STOP
C         END IF
         DEPS_L = 0.0D0
C         ALLOCATE(EPS_N(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR EPS_N (MICROPLANE.F90)'
C           STOP
C         END IF
         EPS_N = 0.0D0
C         ALLOCATE(EPS_M(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR EPS_M (MICROPLANE.F90)'
C           STOP
C         END IF
         EPS_M = 0.0D0
C         ALLOCATE(EPS_L(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR EPS_L (MICROPLANE.F90)'
C           STOP
C         END IF
         EPS_L = 0.0D0
C         ALLOCATE(SIG_N(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR SIG_N (MICROPLANE.F90)'
C           STOP
C         END IF
         SIG_N = 0.0D0
C         ALLOCATE(SD0(1:NP), SN0(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR SD0 OR SN0 (MICROPLANE.F90)'
C           STOP
C         END IF
         SD0 = 0.0D0
         SN0 = 0.0D0
C         ALLOCATE(DED(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR DED (MICROPLANE.F90)'
C           STOP
C         END IF
         DED = 0.0D0
C         ALLOCATE(ED0(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR ED0 (MICROPLANE.F90)'
C           STOP
C         END IF
         ED0 = 0.0D0
C         ALLOCATE(EDF(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR EDF (MICROPLANE.F90)'
C           STOP
C         END IF
         EDF = 0.0D0
C         ALLOCATE(ENF(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR ENF (MICROPLANE.F90)'
C           STOP
C         END IF
         ENF = 0.0D0
C         ALLOCATE(SNF(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR "SNF" '//&
C                &'(MICROPLANE.F90)'
C           STOP
C         END IF
         SNF = 0.0D0
C         ALLOCATE(SNB(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR SNB (MICROPLANE.F90)'
C           STOP
C         END IF
         SNB = 0.0D0
C         ALLOCATE(SNB_FIB(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR SNB_FIB (MICROPLANE.F90)'
C           STOP
C         END IF
         SNB_FIB = 0.0D0
C         ALLOCATE(SDF(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR SDF (MICROPLANE.F90)'
C           STOP
C         END IF
         SDF = 0.0D0
C##############
C         ALLOCATE(SDPOS(1:NP), SDNEG(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR SDPOS OR SDNEG (MICROPLANE.F90)'
C           STOP
C         END IF
         SDPOS = 0.0D0 
         SDNEG=0.0D0
C#############
C         ALLOCATE(EPS_N0_POS(1:NP), EPS_N0_NEG(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR EPS_N0_POS OR EPS_N0_NEG (MICROPLANE.F90)'
C           STOP
C         END IF
         EPS_N0_POS = 0.0D0
         EPS_N0_NEG=0.0D0
C         ALLOCATE(SL0(1:NP), SM0(1:NP), SLF(1:NP), SMF(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR SL0, SM0, SLF OR SMF (MICROPLANE.F90)'
C           STOP
C         END IF
         SL0 = 0.0D0
         SM0=0.0D0
         SLF=0.0D0
         SMF=0.0D0
C         ALLOCATE(DEL(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR DEL (MICROPLANE.F90)'
C           STOP
C         END IF
         DEL = 0.0D0
C         ALLOCATE(DEM(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE VECTOR DEM (MICROPLANE.F90)'
C           STOP
C         END IF
         DEM = 0.0D0
C         ALLOCATE(C_D(1:NP), E_N(1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) ' CAN NOT ALLOCATE C_D OR E_N (MICROPLANE.F90).'
C           STOP
C         END IF   
         C_D=0.0D0
         E_N=0.0D0
         
       PR_SIG_N1 = 0.0D0
       PR_SIG_N2 = 0.0D0
       PR_SIG_N3 = 0.0D0
       PR_SIG_N1_T = 0.0D0
       PR_SIG_N2_T = 0.0D0
       PR_SIG_N3_T = 0.0D0
       STRESS_PR = 0.0D0
       STRESS_PR1 = 0.0D0
       SIG_NEW = 0.0D0


C ALLOCATE TE
C      IF(.NOT.(ALLOCATED(TE))) THEN
C         ALLOCATE(TE(1:4,1:NP), STAT=ALLOCSTAT)
C         IF(ALLOCSTAT.NE.0) THEN
C           WRITE(*,*) 'ERROR: CAN NOT ALLOCATE MATRIX TE (MICROPLANE.F90).'
C           STOP
C         END IF  
C      END IF   

      TE=0.0D0

      IJ=RESHAPE((/1,1,2,2,3,3,1,2,2,3,3,1/),(/2,6/))

      IF(NP.EQ.37) THEN
C                37 PT. FORMULA, DEGREE 9, NO ORTHOGONAL SYMMETRIES
C                     N(3)         N(2)          N(1)            W
C                ------------  ------------  ------------   -------------       
CTE = RESHAPE((/9.822469460D-01,0.000000000D+00,1.875924740D-01,1.984126980D-02,&
C              &3.035309990D-01,-5.257311120D-01,7.946544720D-01,1.984126980D-02,&
C              &3.035309990D-01,5.257311120D-01,7.946544720D-01,1.984126980D-02,&
C              &-4.911234730D-01,-8.506508080D-01,1.875924740D-01,1.984126980D-02,&
C              &-6.070619980D-01,0.000000000D+00,7.946544720D-01,1.984126980D-02,&
C              &-4.911234730D-01,8.506508080D-01,1.875924740D-01,1.984126980D-02,&
C              &7.557613140D-01,-3.090169940D-01,5.773502690D-01,2.539682540D-02,&
C              &7.557613140D-01,3.090169940D-01,5.773502690D-01,2.539682540D-02,&
C              &3.568220900D-01,0.000000000D+00,9.341723590D-01,2.539682540D-02,&
C              &-1.102640900D-01,-8.090169940D-01,5.773502690D-01,2.539682540D-02,&
C              &-1.784110450D-01,-3.090169940D-01,9.341723590D-01,2.539682540D-02,&
C              &-1.784110450D-01,3.090169940D-01,9.341723590D-01,2.539682540D-02,&
C              &-1.102640900D-01,8.090169940D-01,5.773502690D-01,2.539682540D-02,&
C              &-6.454972240D-01,-5.000000000D-01,5.773502690D-01,2.539682540D-02,&
C              &-6.454972240D-01,5.000000000D-01,5.773502690D-01,2.539682540D-02,&
C              &4.670861790D-01,-8.090169940D-01,3.568220900D-01,2.539682540D-02,&
C              &-9.341723590D-01,0.000000000D+00,3.568220900D-01,2.539682540D-02,&
C              &4.670861790D-01,8.090169940D-01,3.568220900D-01,2.539682540D-02,&
C              &8.660254040D-01,-5.000000000D-01,0.000000000D+00,2.539682540D-02,&
C              &-8.660254040D-01,-5.000000000D-01,0.000000000D+00,2.539682540D-02,&
C              &0.000000000D+00,1.000000000D+00,0.000000000D+00,2.539682540D-02/)&      
C              &,(/4,NP/))     
C                21 PT. FORMULA, DEGREE 9, ORTHOGONAL SYMMETRIES
C                     N(3)         N(2)          N(1)            W
C                ------------  ------------  ------------   -------------       
      TE = RESHAPE(
     $(/0.000000000000D+00,0.000000000000D+00,1.000000000000D+00,
     $        1.072388573030D-02,
     $0.000000000000D+00,1.000000000000D+00,0.000000000000D+00,
     $        1.072388573030D-02,
     $1.000000000000D+00,0.000000000000D+00,0.000000000000D+00,
     $        1.072388573030D-02,
     $0.000000000000D+00,7.071067811870D-01,7.071067811870D-01,
     $        2.114160951980D-02,
     $0.000000000000D+00,-7.071067811870D-01,7.071067811870D-01,
     $        2.114160951980D-02,
     $7.071067811870D-01,0.000000000000D+00,7.071067811870D-01,
     $        2.114160951980D-02,
     $-7.071067811870D-01,0.000000000000D+00,7.071067811870D-01,
     $        2.114160951980D-02,
     $7.071067811870D-01,7.071067811870D-01,0.000000000000D+00,
     $        2.114160951980D-02,
     $-7.071067811870D-01,7.071067811870D-01,0.000000000000D+00,
     $        2.114160951980D-02,
     $0.000000000000D+00,3.089512677750D-01,9.510778696510D-01,
     $        5.355055908370D-03,
     $0.000000000000D+00,-3.089512677750D-01,9.510778696510D-01,
     $        5.355055908370D-03,
     $0.000000000000D+00,9.510778696510D-01,3.089512677750D-01,
     $        5.355055908370D-03,
     $0.000000000000D+00,-9.510778696510D-01,3.089512677750D-01,
     $        5.355055908370D-03,
     $3.089512677750D-01,0.000000000000D+00,9.510778696510D-01,
     $        5.355055908370D-03,
     $-3.089512677750D-01,0.000000000000D+00,9.510778696510D-01,
     $        5.355055908370D-03,
     $9.510778696510D-01,0.000000000000D+00,3.089512677750D-01,
     $        5.355055908370D-03,
     $-9.510778696510D-01,0.000000000000D+00,3.089512677750D-01,
     $        5.355055908370D-03,
     $3.089512677750D-01,9.510778696510D-01,0.000000000000D+00,
     $        5.355055908370D-03,
     $-3.089512677750D-01,9.510778696510D-01,0.000000000000D+00,
     $        5.355055908370D-03,
     $9.510778696510D-01,3.089512677750D-01,0.000000000000D+00,
     $        5.355055908370D-03,
     $-9.510778696510D-01,3.089512677750D-01,0.000000000000D+00,
     $        5.355055908370D-03,
     $8.805355183100D-01,3.351545919390D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $-8.805355183100D-01,3.351545919390D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $8.805355183100D-01,-3.351545919390D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $-8.805355183100D-01,-3.351545919390D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $3.351545919390D-01,8.805355183100D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $-3.351545919390D-01,8.805355183100D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $3.351545919390D-01,-8.805355183100D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $-3.351545919390D-01,-8.805355183100D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $3.351545919390D-01,3.351545919390D-01,8.805355183100D-01,
     $        1.677709091560D-02,
     $-3.351545919390D-01,3.351545919390D-01,8.805355183100D-01,
     $        1.677709091560D-02,
     $3.351545919390D-01,-3.351545919390D-01,8.805355183100D-01,
     $        1.677709091560D-02,
     $-3.351545919390D-01,-3.351545919390D-01,8.805355183100D-01,
     $        1.677709091560D-02,
     $5.773502691900D-01,5.773502691900D-01,5.773502691900D-01,
     $        1.884823095080D-02,
     $-5.773502691900D-01,5.773502691900D-01,5.773502691900D-01,
     $        1.884823095080D-02,
     $5.773502691900D-01,-5.773502691900D-01,5.773502691900D-01,
     $        1.884823095080D-02,
     $-5.773502691900D-01,-5.773502691900D-01,5.773502691900D-01,
     $        1.884823095080D-02/)
     $,(/4,NP/))

      END IF   

C     ---------------------------------------
C ... ASSEMBLE TENSORS FROM DIRECTION COSINES
C     ---------------------------------------
      QN = 0.0D0
      QL = 0.0D0
      QM = 0.0D0
       W = 0.0D0
      CALL RANDOM_SEED(SIZE=RS_SIZE)
      DO JP=1,NP 
        W(JP) = TE(4,JP)*6.0D0 
        XN(1) = TE(3,JP) 
        XN(2) = TE(2,JP) 
        XN(3) = TE(1,JP) 
C ROTATE THE GLOBAL AXES HERE: 
C        CALL ROTATE(XN) ! ROTATE ABOUT X AND Y AXES.

C RANDOMLY CHOOSE VECTOR XM. IF XN=XM HAPPENS, CHOOSE XM AGAIN:
        LENGTHM = 0.0D0
        DO WHILE (LENGTHM .LT. EPSILON(LENGTHM))
           CALL RANDOM_NUMBER(RAND_VEC)
           XM = RAND_VEC - DOT_PRODUCT(XN,RAND_VEC)*XN
           LENGTHM = SQRT(DOT_PRODUCT(XM,XM))
        END DO
        XM = XM/LENGTHM

C CALCULATE VECTOR XL = XN .CROSS. XM
        XL(1) = XN(2)*XM(3)-XN(3)*XM(2) 
        XL(2) = XN(3)*XM(1)-XN(1)*XM(3) 
        XL(3) = XN(1)*XM(2)-XN(2)*XM(1) 
        LENGTHL = SQRT(DOT_PRODUCT(XL,XL))
        XL=XL/LENGTHL
        LENGTHN = SQRT(DOT_PRODUCT(XN,XN))
        LENGTHM = SQRT(DOT_PRODUCT(XM,XM))
        LENGTHL = SQRT(DOT_PRODUCT(XL,XL))
        DO K=1,6 
          I=IJ(1,K) 
          J=IJ(2,K) 
          QN(K,JP) = XN(I)*XN(J) 
          QM(K,JP) = 0.5D0*(XN(I)*XM(J)+XN(J)*XM(I)) 
          QL(K,JP) = 0.5D0*(XN(I)*XL(J)+XN(J)*XL(I)) 
        END DO 
      END DO 
      RETURN
      END SUBROUTINE SETSYSTEM 

C *******************************************************************
C *** SUBROUTINE M7FMATERIAL *****************************************
C *******************************************************************
C      SUBROUTINE M7FMATERIAL(CONVERGED, EPS, DEPS, SIG, SIG_OLD) 
      SUBROUTINE M7FMATERIAL(
     $     NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL,
     $     STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH,
     $     PROPS, DENSITY, STRAININC, RELSPININC,
     $     TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,
     $     STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,
     $     TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW,
     $     STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW )

C  !
C  ! EXPLANATIONS OF THE ARGUMENTS OF SUBROUTINE VUMAT
C  !
C  ! VARIABLES TO BE DEFINED
C  !
C  ! STRESSNEW(:) : STRESSES AT THE END OF THE INCREMENT,                   (OUTPUT)
C  ! STATENEW(:)  : STATE VARIABLES AT THE END OF THE INCREMENT,            (OUTPUT)
C  !
C  ! VARIABLES THAT CAN BE UPDATED
C  !
C  ! ENERINTERNNEW(:): INTERNAL ENERGY PER UNIT MASS AT EACH MAT. POINT AT
C  !                THE END OF THE INCREMENT                                (OUTPUT)
C  ! ENERINELASNEW(:): DISSIPATED INELASTIC ENERGY PER UNIT MASS AT EACH
C  !                MAT. POINT AT THE END OF THE INCREMENT                  (OUTPUT)
C  !
C  ! VARIABLES FOR INFORMATION (OR INPUT ONLY)
C  !
C  ! NBLOCK       : NO. OF MAT. PTS TO BE PROCESSED IN THIS CALL TO VUMAT    (INPUT)
C  ! NDIR         : NO. OF DIRECT COMPONENTS IN A SYMMETRIC TENSOR           (INPUT)
C  ! NSHR         : NO. OF INDIRECT COMPONENTS IN A SYMMETRIC TENSOR         (INPUT)
C  ! NSTATEV      : NO. OF USER-DEFINED STATE VARS FOR THIS MATERIAL         (INPUT)
C  ! NFIELDV      : NO. OF USER-DEFINED EXTERNAL FIELD VARS                  (INPUT)
C  ! NPORPS       : NO. OF USER-DEFINED MAT. PROPERTIES                      (INPUT)
C  ! LANNEAL      : =1 INDICATES THAT VUMAT IS CALLED DURING ANNEALING       (INPUT)   
C  ! STEPTIME     : TIME SINCE THE BEGINNING OF STEP                         (INPUT)
C  ! TOTALTIME    : TOTAL TIME SINCE THE BEGINNING OF FIRST STEP             (INPUT)
C  ! DT           : TIME INCREMENT SIZE                                      (INPUT)
C  ! CMNAME       : USER SPECIFIED MATERIAL NAME, LEFT JUSTIFIED.   
C  ! COORDMP(:,:) : MAT. PT. COORDS.                                         (INPUT)
C  ! CHARLENGTH(:): CHARACTERISTIC ELEMENT LENGTH                            (INPUT)
C  ! PROPS(:)     : USER-SUPPLIED MATERIAL PROPS                             (INPUT)
C  ! DENSITY(:)   : CURRENT DENSITY AT MAT. PTS. IN THE MIDSTEP CONF.        (INPUT)
C  ! STRAININC(:,:) : STRAIN INCREMENT TENSOR AT EACH MATERIAL POINT         (INPUT)
C  ! RELSPININC(:,:): INCREMENTAL RELATIVE ROTATION VECTOR AT EACH MAT. PT.  (INPUT)
C  ! TEMPOLD(:,:)   : TEMPERATURES AT THE MAT. PTS. AT THE BEGINNING OF STEP (INPUT)
C  ! STRETCHOLD(:,:): STRETCH TENSOR U, AT EACH MAT. PT. AT THE BEG. OF STEP (INPUT)
C  ! DEFGRADOLD(:,:): DEF. GRAD. TENSOR F, AT EACH MAT.PT. AT BEG. OF STEP   (INPUT)
C  ! FIELDOLD(:,:)  : VALUES OF THE USER DEFINED FIELD VARS AT EACH MAT. PT.
C  !                AT THE BEGINNING OF THE STEP                             (INPUT)
C  ! STRESSOLD(:,:) : STRESS TENSOR AT EACH MAT. PT. AT THE BEG. OF THE STEP (INPUT)
C  ! STATEOLD(:,:)  : STATE VARIABLES AT EACH MAT. PT. AT THE BEG. OF STEP   (INPUT)
C  ! ENERINTERNOLD(:) : INTERNAL ENERGY PER UNIT MASS AT EACH MAT. PT. AT  
C  !                    THE BEG. OF THE INCR.                                (INPUT)
C  ! ENERINELASOLD(:) : DISSIPATED INELASTIC ENERGY PER UNIT MASS AT EACH  
C  !                    MAT. PT. AT THE BEG. OF THE INCR.                    (INPUT)
C  ! TEMPNEW(:)   : TEMPERATURES AT EACH MAT. PT. AT THE END OF INCREMENT    (INPUT)
C  ! STRETCHNEW(:,:): STRETCH TENSOR, U, AT EACH MAT. PT. AT THE END OF INCR.(INPUT)
C  ! DEFGRADNEW(:,:):  OF HOW TO PROVIDE SOME OF THE MATERIAL SPECIFIC
C  ! FIELDNEW(:,:)  : VALUES OF USER-DEFINED FIELD VARS AT EACH MAT.PT. AT THE
C  !                END OF INCREMENT
C  !
C  !
C  !
C  ! OTHER USEFUL INFORMATION 
C  ! (1) SHEAR STRAINS ARE USED INSTEAD OF SHEAR ANGLES AS OPPOSED TO UMAT.  
C  ! (2) STRAIN ORDER IS \EPS_1, \EPS_2, \EPS_3, \EPS_{12}, \EPS_{13}, \EPS_{23}.
C  ! (3) WHEN THE LOCAL SYSTEM IS DEFINED OVER THE ELEMENTS FOR WHICH THIS ROUTINE
C  !     IS CALLED AT EVERY GAUSS POINT, THE ROTATION MATRIX IS NOT NEEDED, BECAUSE
C  !     STRESSES AND STRAINS IN THIS SUBROUTINE ARE DEFINED WITH RESPECT TO THAT
C  !     LOCAL COORDINATE SYSTEM.
C  !
C  !
      INCLUDE 'VABA_PARAM.INC'

      EXTERNAL RATEFUNC
C
      DIMENSION PROPS(NPROPS), DENSITY(NBLOCK),
     $COORDMP(NBLOCK,*),
     $CHARLENGTH(NBLOCK), STRAININC(NBLOCK,NDIR+NSHR),
     $RELSPININC(NBLOCK,NSHR), TEMPOLD(NBLOCK),
     $STRETCHOLD(NBLOCK,NDIR+NSHR), DEFGRADOLD(NBLOCK,NDIR+NSHR+NSHR),
     $FIELDOLD(NBLOCK,NFIELDV), STRESSOLD(NBLOCK,NDIR+NSHR),
     $STATEOLD(NBLOCK,NSTATEV), ENERINTERNOLD(NBLOCK),
     $ENERINELASOLD(NBLOCK), TEMPNEW(NBLOCK),
     $STRETCHNEW(NBLOCK,NDIR+NSHR), DEFGRADNEW(NBLOCK,NDIR+NSHR+NSHR), 
     $FIELDNEW(NBLOCK,NFIELDV),
     $STRESSNEW(NBLOCK,NDIR+NSHR), STATENEW(NBLOCK,NSTATEV),
     $ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK)
C  ! ALL ARRAYS DIMENSIONED BY (*) ARE NOT USED IN THIS ALGORITHM
C      DIMENSION PROPS(NPROPS), DENSITY(NBLOCK),
C     $     COORDMP(NBLOCK,*),
C     $     CHARLENGTH(*), STRAININC(NBLOCK,NDIR+NSHR),
C     $     RELSPININC(*), TEMPOLD(*),
C     $     STRETCHOLD(*), DEFGRADOLD(*),
C     $     FIELDOLD(*), STRESSOLD(NBLOCK,NDIR+NSHR),
C     $     STATEOLD(NBLOCK,NSTATEV), ENERINTERNOLD(NBLOCK),
C     $     ENERINELASOLD(NBLOCK), TEMPNEW(*),
C     $     STRETCHNEW(*), DEFGRADNEW(*), FIELDNEW(*),
C     $     STRESSNEW(NBLOCK,NDIR+NSHR), STATENEW(NBLOCK,NSTATEV),
C     $     ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK)

C  !      REAL*8 :: STEPTIME, TOTALTIME, DT
C  !
      CHARACTER*80 CMNAME

C  !      INTEGER :: NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL, NTENSC
C
      REAL*8 :: EQUIVSTRESS, FRACTUREWORKINC, SMEAN, STRESSPOWER, TWO
      ! DIMENSION STRESS_PR(NBLOCK,NDIR+NSHR)
C  !
C  ! ADD USER DEFINED LOCAL VARIABLES HERE
C  !
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KIFED_M7F/ STRAINRATETENSORINV2, DEVSTRAINRATETENSORINV2,
     $     OLD_DTIME, LRUNCREEP, LFIRSTPASS 
      REAL*8 :: STRAINRATETENSORINV2,
     $     DEVSTRAINRATETENSORINV2, OLD_DTIME
      COMMON /KM7F_M7FIO/ ON_RATE_DEP_FRAC, ON_CREEP
      COMMON /KM7F_M7FIO_1/ C_R1, Q_1
      COMMON /KM7F_M7FIO_2/ T_0, T_CH, F_C, A_TO_C, G_C, W_TO_C, E28
      INTEGER      :: ON_RATE_DEP_FRAC, ON_CREEP
      REAL*8 :: C_R1, Q_1  
      REAL*8 :: T_0, T_CH, F_C, A_TO_C, G_C, W_TO_C, E28
      COMMON /KM7F_MICROPLANE_1/ QN, QL, QM, W 
      COMMON /KM7F_MICROPLANE_2/ DEPS_N, DEPS_L, DEPS_M, EPS_N,  EPS_L,
     $     EPS_M,  SIG_N
      COMMON /KM7F_MICROPLANE_3/ SD0, SN0, DED, ED0, EDF, ENF, SNF, SNB,
     $     SDF
      COMMON /KM7F_MICROPLANE_4/ SL0, SM0, SLF, SMF, DEL, DEM, VH_INI,
     $     VH_FIN, C_D, E_N
      COMMON /KM7F_MICROPLANE_5/ SDNEG, SDPOS, EPS_N0_NEG, EPS_N0_POS,
     $     SNB_FIB, CV
      REAL*8, DIMENSION(1:6,1:NP) :: QN, QL, QM
      REAL*8, DIMENSION(1:NP)     :: W
      REAL*8, DIMENSION(1:NP)     :: DEPS_N,
     $     DEPS_L, DEPS_M, EPS_N,  EPS_L,  EPS_M,  SIG_N
      REAL*8, DIMENSION(1:NP)     :: SD0, SN0,
     $     DED, ED0, EDF
      REAL*8, DIMENSION(1:NP)     :: ENF, SNF,
     $     SNB, SDF, SDNEG, SDPOS, SNB_FIB
      REAL*8, DIMENSION(1:NP)     :: SL0, SM0, SLF,
     $     SMF, DEL, DEM, EPS_N0_NEG, EPS_N0_POS
C      REAL*8, DIMENSION(1:NVHI)    :: VH_INI, VH_FIN 
      REAL*8, DIMENSION(1:NVHI+2)    :: VH_INI, VH_FIN ! ADDITIONAL HIST. VAR. FOR ELEMENT DELETION 
      REAL*8, DIMENSION(1:NP)     :: C_D, E_N
      REAL*8 :: CV
      SAVE    :: /KCOMMON_VARS/, /KIFED_M7F/
      SAVE    :: /KM7F_M7FIO/, /KM7F_M7FIO_1/, /KM7F_M7FIO_2/
      SAVE    :: /KM7F_MICROPLANE_1/, /KM7F_MICROPLANE_2/,
     $     /KM7F_MICROPLANE_3/, /KM7F_MICROPLANE_4/, /KM7F_MICROPLANE_5/

      REAL*8, DIMENSION(1:6) :: EPS, DEPS, SIG_OLD 
      REAL*8, DIMENSION(1:6) :: SIG
      REAL*8, DIMENSION(1:6) :: DEPS_S
      REAL*8 :: R_N, DFB
      REAL*8 :: SV0, EV0, DEV, EVF, SUM_SNF
      REAL*8 :: SVFB, SVF, PHI0, PHI, SIG_I_OLD, EPS_I
      INTEGER ::  JP, I , JM
      REAL*8 :: SVNEG, SUM1, STRAIN_RATE, EK(1:6,1:6)
      
      COMMON /EIG_VAL_VEC/PR_SIG_N1,PR_SIG_N2,PR_SIG_N3,PR_SIG_N1_T,
     $   PR_SIG_N2_T,PR_SIG_N3_T,STRESS_PR,STRESS_PR1,SIG_NEW

	  REAL*8,DIMENSION(3) :: PR_SIG_N1,PR_SIG_N2,PR_SIG_N3
	  REAL*8,DIMENSION(3) :: STRESS_PR,STRESS_PR1
	  REAL*8,DIMENSION(3,3) :: PR_SIG_N1_T,PR_SIG_N2_T,PR_SIG_N3_T
	  REAL*8,DIMENSION(6) :: SIG_NEW 
	  REAL*8,DIMENSION(nblock,ndir) :: eigVal
	  REAL*8,DIMENSION(nblock,3,3) :: eigVec
!	  REAL*8,DIMENSION(NBLOCK,ndir+nshr)::SIG_PR
!	  REAL*8,DIMENSION(3) :: eigVal
!	  REAL*8,DIMENSION(3,3) :: eigVec
	  REAL*8,DIMENSION(NBLOCK,ndir+nshr)::SIG_PR
	  
      NTENS = 6
      
C     -------------------------------------------------------------------
C ... ALLOCATE AND ASSIGN DATA TO TENSORS AND GAUSSIAN QUADRATURE WEIGHTS;
C     ALLOCATE THE VECTORS TO USE: SET THE MICROPLANE SYSTEM
C     -------------------------------------------------------------------
      
C      IF (LFIRSTPASS) THEN 
C        CALL SETSYSTEM() 
C        LFIRSTPASS=.FALSE.
C      END IF
C      PRINT *, "ENTERED M7FMATERIAL."
      IF (TOTALTIME<=DT) THEN
         CALL INPUTPARAMS()
         CALL SETSYSTEM() 
      END IF
C                   !$OMPXPARALLEL DO
      DO I = 1, NBLOCK

C HISTORY VARIABLES
         VH_INI = STATEOLD(I,:) 
C INITIALIZE THE HISTORY VARIABLES THAT ARE NONZERO AT THE BEGINNING OF THE SIMULATION.
         IF (TOTALTIME<=DT) THEN
            VH_INI(2)=1.0D0
         END IF

C     ! INITIALIZE THE STRESS TENSOR; CALCULATE THE STRAIN INCREMENT TENSOR  
         DO IC = 1, NTENS
            SIG(IC)=STRESSOLD(I,IC)
            SIG_OLD(IC)=STRESSOLD(I,IC)
C ---------------------------------------------------------
C WHAT THE PARAMETER UNIT_CONV DOES IS, IT MAKES SURE THAT
C ON ENTRY TO M7FMATERIAL, THE STRESSES ARE IN MPA AND ON
C EXIT FROM M7FMATERIAL, THE STRESSES AND FORCES ARE IN THE 
C UNITS DICTATED BY THE MESH UNITS.
C IN M7F, THE STRESSES ARE IN MPA. THIS REQUIRES THE MESH
C TO BE IN MM, FORCES IN N AND MODULI IN MPA. IN THAT CASE,
C UNIT_CONV=1.
C HOWEVER, IF THE MESH IS IN M, FORCES IN N AND MODULI ARE 
C IN PA, UNIT_CONV SHOULD BE SET TO 1.D-6 IN THE SUBROUTINE 
C INPUTPARAMS(). 
C OTHER UNIT CONVERSIONS CAN ALSO BE HANDLED THROUGH SETTING
C UNIT_CONV TO THE RIGHT PARAMETER. FOR EXAMPLE, IF THE MESH
C IS IN INCHES, AND THE FORCES IN POUNDS, UNITCONV MUST BE
C SET TO 6.89475729D-3.
C ----------------------------------------------------------
            SIG(IC)=SIG(IC)*UNIT_CONV
            SIG_OLD(IC)=SIG_OLD(IC)*UNIT_CONV
            IF (IC > 3) THEN    ! SHEAR STRAINS ARE SHEAR ANGLES IN MICROPLANE MODEL
               DEPS(IC)=2.0*STRAININC(I,IC)
            ELSE
               DEPS(IC)=STRAININC(I,IC)
            END IF
         END DO
C      ! CALCULATE THE TOTAL STRAIN TENSOR NEEDED IN MICROPLANE MODEL
         EPS(1)=(DEFGRADOLD(I,1)*DEFGRADOLD(I,1) + 
     $           DEFGRADOLD(I,7)*DEFGRADOLD(I,7) +
     $           DEFGRADOLD(I,6)*DEFGRADOLD(I,6) - 1.0)/2.0
         EPS(2)=(DEFGRADOLD(I,4)*DEFGRADOLD(I,4) + 
     $           DEFGRADOLD(I,2)*DEFGRADOLD(I,2) +
     $           DEFGRADOLD(I,8)*DEFGRADOLD(I,8) - 1.0)/2.0
         EPS(3)=(DEFGRADOLD(I,9)*DEFGRADOLD(I,9) + 
     $           DEFGRADOLD(I,5)*DEFGRADOLD(I,5) +
     $           DEFGRADOLD(I,3)*DEFGRADOLD(I,3) - 1.0)/2.0
         EPS(4)=(DEFGRADOLD(I,1)*DEFGRADOLD(I,4) +
     $           DEFGRADOLD(I,2)*DEFGRADOLD(I,7) +
     $           DEFGRADOLD(I,6)*DEFGRADOLD(I,8))/2.0
         EPS(5)=(DEFGRADOLD(I,9)*DEFGRADOLD(I,4) + 
     $           DEFGRADOLD(I,2)*DEFGRADOLD(I,5) +
     $           DEFGRADOLD(I,3)*DEFGRADOLD(I,8))/2.0
         EPS(6)=(DEFGRADOLD(I,1)*DEFGRADOLD(I,9) + 
     $           DEFGRADOLD(I,7)*DEFGRADOLD(I,5) +
     $           DEFGRADOLD(I,3)*DEFGRADOLD(I,6))/2.0
         DO IC = 1, NTENS
            IF (IC > 3) THEN    ! SHEAR STRAINS ARE SHEAR ANGLES IN MICROPLANE MODEL
               EPS(IC)=2.0*EPS(IC)
            END IF
         END DO

C         EPS = EPS + DEPS !11.09.2011: IS THIS REALLY NECESSARY? IT SHOULD NOT BE.

C      WRITE(6,*), "EPS=",EPS
C      WRITE(6,*), "SIG=",SIG
C      IF (CONVERGED) THEN
C NEW RATE FORMULATION 23.08.2011
         DEPS_S(1)=DEPS(1)*DEPS(1)+DEPS(4)*DEPS(4)+DEPS(6)*DEPS(6)
         DEPS_S(2)=DEPS(4)*DEPS(4)+DEPS(2)*DEPS(2)+DEPS(5)*DEPS(5)
         DEPS_S(3)=DEPS(6)*DEPS(6)+DEPS(5)*DEPS(5)+DEPS(3)*DEPS(3)
         DEPS_S(4)=DEPS(1)*DEPS(4)+DEPS(2)*DEPS(4)+DEPS(5)*DEPS(6)
         DEPS_S(5)=DEPS(4)*DEPS(6)+DEPS(2)*DEPS(5)+DEPS(5)*DEPS(3)
         DEPS_S(6)=DEPS(1)*DEPS(6)+DEPS(4)*DEPS(5)+DEPS(6)*DEPS(3)
         SUM1=0.D0
         DO JP = 1, NP
            SUM1=SUM1+(QN(1,JP)*DEPS_S(1)+QN(4,JP)*DEPS_S(4)+
     $           QN(6,JP)*DEPS_S(6)+QN(4,JP)*DEPS_S(4)+
     $           QN(2,JP)*DEPS_S(2)+QN(5,JP)*DEPS_S(5)+
     $           QN(6,JP)*DEPS_S(6)+QN(5,JP)*DEPS_S(5)+
     $           QN(3,JP)*DEPS_S(3))*W(JP)
         END DO
         STRAIN_RATE=SQRT(SUM1/2.D0)/DT
C END NEW RATE FORMULATION 23.08.2011
C     ------------------------ 
C ... UPDATE HISTORY VARIABLES
C     ------------------------ 
C        VH_INI = VH_FIN
C        CONVERGED = .FALSE.
C      ENDIF
 
C     -------------------------------------------------------------------
C ... CALL CREEP ROUTINE ONLY WHEN DT IS DIFFERENT THAN PREVIOUS STEP:
C     -------------------------------------------------------------------
      
C      IF(DT.EQ.OLD_DT) THEN
C         LRUNCREEP=.FALSE.
C         OLD_DT = DT
C         WRITE(*,*) 'TIME ICREMENT IN THIS STEP = ', DT
C      END IF

C     -------------------------------------------------------------------
C     THE CREEP SUBROUTINE. WHEN CALLED, IT MODIFIES THE ELASTIC MODULUS
C ... AND STORES THE NEW C0 AS C0 FOR USE IN THE MATERIAL SUBROUTINES.
C     ON_CREEP IS THE ON/OFF SWITCH FOR THE CREEP ROUTINE.
C     -------------------------------------------------------------------
      
C      IF((ON_CREEP.NE.0).AND.LRUNCREEP) THEN
C         CALL M7F_CREEP(DT)
C         LRUNCREEP=.FALSE.
C         WRITE(*,*) 'E IN THIS LOADING STEP IS ', YOUNG
C         WRITE(*,*) 'C0 IN THIS LOADING STEP IS', C0
C      END IF

C     -----------------
C ... HISTORY VARIABLES
C     -----------------
 
C     -------------------------------------------------------------------
C     COMPUTE MATERIAL RESPONSE :
C     BOUNDING CURVE FORMULATION FOR VOLUMETRIC,
C ... DEVIATORIC AND SHEAR VARIABLES.
C     SHEAR STRAINS ARE SHEAR ANGLES.
C     -------------------------------------------------------------------

C     ---------------
C ... INITIALIZATIONS
C     ---------------

      JP = SIZE(W) ! NO. OF MICROPLANES      
      
      SIG = 0.0D0
       
C     -------------------------------------------------------------------
C ... CALCULATE THE RATE COEFFICIENTS USING THE RATE FUNCTION
C     -------------------------------------------------------------------
      
      IF(ON_RATE_DEP_FRAC.EQ.0) THEN !RATE DEPENDENT FRACTURING EFFECTS ARE OFF
        R_N=0.0D0
      ELSE  
C NEW RATE FORMULATION 23.08.2011
        R_N = RATEFUNC(STRAIN_RATE/C_R1) 
C END NEW RATE FORMULATION 23.08.2011
      END IF  

C     --------------------------------------------
C ... EVALUATE LOCAL STRAINS (EPS_N, EPS_L, EPS_M)
C     --------------------------------------------
  
       EPS_N =  MATMUL( EPS,QN) 
       EPS_L =  MATMUL( EPS,QL)
       EPS_M =  MATMUL( EPS,QM)
      DEPS_N =  MATMUL(DEPS,QN)
      DEPS_L =  MATMUL(DEPS,QL)
      DEPS_M =  MATMUL(DEPS,QM)

C MATERIAL POINT DELETION FOR ELEMENTS DISTORTED TOO MUCH
      IF (MAXVAL(EPS_N)>TH_DEL) THEN
         VH_FIN(188)=0.
      ELSE
         VH_FIN(188)=1.
      END IF

C     -----------------------------------------
C ... VOLUMETRIC LAW (SAME FOR ALL MICROPLANES)
C     -----------------------------------------
 
      SV0 = VH_INI(1)
      PHI0 = VH_INI(2)

C HISTORY VARIABLES (NVHI = NVHM*NP+1)
C WHERE NVHM=5 (5 VARIABLES PER MICROPLANE ARE STORED) AND
C VH_INI(1)=VOLUMETRIC STRESS
C VH_INI(2:NVHI:NVHM)=MICROPLANE NORMAL STRESS
C VH_INI(3:NVHI:NVHM)=MICROPLANE L-DIR SHEAR STRESS
C VH_INI(4:NVHI:NVHM)=MICROPLANE M-DIR SHEAR STRESS
C VH_INI(5:NVHI:NVHM)=MAX STRAIN THAT CAN CHARACTERIZE TENSION DAMAGE
C VH_INI(6:NVHI:NVHM)=MAX STRAIN THAT CAN CHARACTERIZE COMPRESSION DAMAGE

C     -------------------------------------------------------------------
C ... THE FOLLOWING ARE VALID ONLY FOR SMALL STRAINS 
C     -------------------------------------------------------------------
      
      EV0 =  (EPS(1)+ EPS(2)+ EPS(3))/3.0D0
      DEV = (DEPS(1)+DEPS(2)+DEPS(3))/3.0D0
      EVF = EV0+DEV 

      CALL C_VOL(DEV, EV0, SV0, DEPS_N, EPS_N, R_N, SVNEG)

C     ---------------------
C ... NORMAL DEVIATORIC LAW
C     ---------------------

      SN0 = VH_INI(3:NVHI:NVHM)

      DED = DEPS_N - DEV
      ED0 = EPS_N - EV0     
      EDF = ED0 + DED

      CALL C_DEV(DED, ED0, DEV, EV0, SV0, SD0, SDF, C_D, R_N,
     $     SDNEG, SDPOS)

      ENF = EVF + EDF 

      EPS_N0_POS=VH_INI(6:NVHI:NVHM)
      EPS_N0_NEG=VH_INI(7:NVHI:NVHM)
      CALL C_NORM_ELASTIC(EPS_N, DEPS_N, SN0, EPS_N0_POS,
     $     EPS_N0_NEG, SV0, SNF, E_N)
C THE SUBROUTINE C_NORM_ELASTIC DOES NOT HAVE A BOUNDARY.

C     ---------------------
C ... TENSILE NORMAL LAW
C     ---------------------
      CALL C_NORM(ENF, SV0, R_N, SNB)
      CALL C_NORM_FIB(ENF, R_N, SNB_FIB)

C     ------------------------------------------------------------------------
C ... WIPE OUT THE NORMAL BOUNDARY FOR LARGE MAX. PRINCIPAL STRAIN. 28.08.2011
C     ------------------------------------------------------------------------
      PHI0 = VH_INI(2)
!      PHI = 1.0D0
!      SIG_I_OLD =MAXVAL(MATMUL( SIG_OLD,QN)) ! APPROXIMATE VALUE OF MAX. PRINCIPAL STRESS
!      IF (SIG_I_OLD > 1.0D-3) THEN ! WE ARE IN TENSION.
!         EPS_I = MAXVAL(EPS_N)
!         PHI = EXP(-K_5*MAX(EPS_I/K_1-C_16,0.0D0))
!      END IF
!      IF (PHI < PHI0) THEN
!         SNB = SNB*PHI
!         SNB_FIB = SNB_FIB*PHI
!         VH_FIN(2)=PHI
!      ELSE
!         SNB = SNB*PHI0
!         SNB_FIB = SNB_FIB*PHI0
!         VH_FIN(2)=PHI0
!      END IF

C IMPORTANT: SNB_FIB SHOULD BE ZERO WHEN V_F=0% (I.E. FOR PLAIN CONCRETE). 
      SIG_N = MAX(MIN(SNF,SNB+SNB_FIB),SVNEG+SDNEG) 
C      SIG_N = MAX(MIN(SNF,SNB+SNB_FIB),SVNEG+SDNEG) 
C ALSO, NOTE THAT SNF AND SIG_N BOTH ARE NEEDED TO DETECT THE EPS_N0_POS AND EPS_N0_NEG.
      DO IC = 1, JP
         IF (SNF(IC) > SNB(IC)+SNB_FIB(IC)) THEN ! THE TENSILE NORMAL B. IS EXCEEDED
C            IF(TOTALTIME<1.D0) THEN
               EPS_N0_POS(IC) = ENF(IC)
C            END IF
         ELSE IF (SNF(IC) < SVNEG+SDNEG(IC)) THEN ! THE COMPRESSIVE NORMAL B. IS EXCEEDED
            EPS_N0_NEG(IC) = ENF(IC)
         END IF
      END DO

C     ------------------------------------------------------
C ... ADJUST VOLUMETRIC STRESS (ACCORDING TO SIG_N <= SIG_B)
C     SUM_SNF IS $3/(2*PI)\INT_\OMEGA \SIGMA_N D\OMEGA$, BECAUSE THE 
C     WEIGHTS W(JP) ARE SCALED ACCORDINGLY.      
C     ------------------------------------------------------
      
      SUM_SNF = DOT_PRODUCT(SIG_N,W)              
      SVFB = SUM_SNF/3.0D0 !SVFB IS THE CORRECT VALUE OF THE AVERAGE OF \SIGMA_N OVER THE SPHERE. 
C      SVF = MIN(SVF,SVFB)
C NOTE THE FOLLOWING: SVF IS NEEDED FOR CORRECT BEHAVIOR UNDER UNIAXIAL COMPRESSION. OTHERWISE THE 
C UNIAXIAL COMPRESSION IS SIMILAR TO UNIAXIAL TENSION WITH A SHARP POSTPEAK THAT CANNOT BE CONTROLLED. 
      SVF = SVFB

      SL0 = VH_INI(4:NVHI:NVHM) 
      SM0 = VH_INI(5:NVHI:NVHM) 
 
      SNF = SIG_N

      ENF = EPS_N + DEPS_N
      DEL = DEPS_L
      DEM = DEPS_M

C     ---------
C ... SHEAR LAW
C     ---------

      CALL C_SHEAR2(EPS_L, EPS_M, C_D, SNF, R_N, DEL, DEM, SL0, SM0,
     $     SLF, SMF, EVF)

C     -------------------------------------------------
C     FINAL STRESS VECTOR
C     -------------------------------------------------
      
      SIG = MATMUL(QN,SNF*W) + MATMUL(QM,SMF*W) + MATMUL(QL,SLF*W)

C STRESSES ARE IN PA IN ABAQUS, AND IN MPA IN M7FMATERIAL
C *******************************************************
      SIG = SIG/UNIT_CONV

C      EEE=25000.
C      XNU=0.18
C! 
CC      WRITE(*,'("I AM AT 2360. YOUNG =",6G15.5)') PI, YOUNG, POISSON, K_1, C0, DTIME
C! 
C      EK = 0.0D0
C! 
C       EK(1,1)=EEE*(1.0D0-XNU)/((1.0D0+XNU)*(1.0D0-2.0D0*XNU))
C       EK(1,2)=EEE*XNU/((1.0D0+XNU)*(1.0D0-2.0D0*XNU))
C       EK(1,3)=EK(1,2)
C       EK(2,1)=EK(1,2)
C       EK(2,2)=EK(1,1)
C       EK(2,3)=EK(1,2)
C       EK(3,1)=EK(1,2)
C       EK(3,2)=EK(1,2)
C       EK(3,3)=EK(1,1)
C       EK(4,4)=EEE/(2.0D0+2.0D0*XNU)
C       EK(5,5)=EK(4,4)
C       EK(6,6)=EK(4,4)
C
C       SIG=SIG+MATMUL(EK,DEPS)
        
C     -------------------------------------------
C ... UPDATE MICROPLANE NORMAL AND SHEAR STRESSES
C     -------------------------------------------
      
      VH_FIN(1) = SVF
C      VH_FIN(2) = UPDATED IN THE CODE.
      VH_FIN(3:NVHI:NVHM) = SNF 
      VH_FIN(4:NVHI:NVHM) = SLF 
      VH_FIN(5:NVHI:NVHM) = SMF 
      VH_FIN(6:NVHI:NVHM) = EPS_N0_POS
      VH_FIN(7:NVHI:NVHM) = EPS_N0_NEG 
C      PRINT*, "EPS=",EPS
c      PRINT*, "SIG=",SIG
C----- END OF ORIGINAL SUBROUTINE
      STATENEW(I,:)=VH_FIN
      STRESSNEW(I,:)=SIG


	  
C     ! UPDATE THE SPECIFIC INTERNAL ENERGY -
      STRESSPOWER = (( STRESSOLD(I,1)+STRESSNEW(I,1) )*
     1      STRAININC(I,1) + ( STRESSOLD(I,2)+STRESSNEW(I,2) )*
     2        STRAININC(I,2) + ( STRESSOLD(I,3)+STRESSNEW(I,3) )*
     3        STRAININC(I,3) + 2.0*( STRESSOLD(I,4)+STRESSNEW(I,4) )*
     4        STRAININC(I,4)+2.0*( STRESSOLD(I,5)+STRESSNEW(I,5) )*
     5	      STRAININC(I,5)+2.0*( STRESSOLD(I,6)+STRESSNEW(I,6) )*
     6        STRAININC(I,6))/2.0
C     !
      ENERINTERNNEW(I) = ENERINTERNOLD(I) + STRESSPOWER / DENSITY(I)
C     !
C     ! UPDATE THE DISSIPATED INELASTIC SPECIFIC ENERGY -
      SMEAN = ( STRESSNEW(I,1) + STRESSNEW(I,2) +
     $        STRESSNEW(I,3) )/3.0
      EQUIVSTRESS = SQRT( 3.0/2.0 * ( (STRESSNEW(I,1)-SMEAN)**2 +
     $        (STRESSNEW(I,2)-SMEAN)**2 + (STRESSNEW(I,3)-SMEAN)**2 +
     $        2. * STRESSNEW(I,4)**2+2. * STRESSNEW(I,5)**2+
     $        2. * STRESSNEW(I,6)**2) )

      FRACTUREWORKINC = STRESSNEW(I,1)*DEPS(1) +
     $        STRESSNEW(I,2)*DEPS(2) + STRESSNEW(I,3)*DEPS(3) +
     $        STRESSNEW(I,4)*DEPS(4) + STRESSNEW(I,5)*DEPS(5) +
     $        STRESSNEW(I,6)*DEPS(6) 
      ENERINELASNEW(I) = ENERINELASOLD(I)+FRACTUREWORKINC/DENSITY(I)

      END DO
C                                  !$OMPXPARALLEL END DO
      RETURN
      END SUBROUTINE M7FMATERIAL


C***********************************************************************
C**** THIS FUNCTION CHARACTERIZES THE THEORETICAL RATE EFFECT      *****
C***********************************************************************

      FUNCTION RATEFUNC(ARG)
      INCLUDE 'VABA_PARAM.INC'
      REAL*8 :: RATEFUNC
      REAL*8 :: ARG

C THE FUNCTION ASINH(X) IS DEFINED BELOW FOR X IN [0, INFINITY).
      RATEFUNC = LOG(ARG + SQRT(ARG * ARG + 1.0D0))
      RETURN
      END FUNCTION RATEFUNC

C **********************************************************************
C *** SUBROUTINE CREEP *************************************************
C **********************************************************************
      SUBROUTINE M7F_CREEP(DT)
C+---------------------------------------------------------------------+
C|  REPLACES THE ELASTIC MODULUS E BY THE NEW ONE EP WHICH ACCOUNTS FOR|
C|  CREEP EFFECTS AND STORES THE COEFFICIENT C0                        |
C+---------------------------------------------------------------------+
      INCLUDE 'VABA_PARAM.INC'
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER ::  NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7F_M7FIO/ ON_RATE_DEP_FRAC, ON_CREEP
      COMMON /KM7F_M7FIO_1/ C_R1, Q_1
      COMMON /KM7F_M7FIO_2/ T_0, T_CH, F_C, A_TO_C, G_C, W_TO_C, E28
      INTEGER      :: ON_RATE_DEP_FRAC, ON_CREEP
      REAL*8 :: C_R1, Q_1  
      REAL*8 :: T_0, T_CH, F_C, A_TO_C, G_C, W_TO_C, E28
      REAL*8 :: DT
      REAL*8, PARAMETER  :: DAYSPERSEC = 1.1574074074074074D-5 
      REAL*8 :: T_AGE, T_CH_P, DENSITY, Q1_OPT
      REAL*8 :: Q1, Q2, Q3, Q4, N, M, T, J, R, Z, QF, Q, J_DOT
      REAL*8 :: EP, ETA, TAU1, C1
      
C      DAYSPERSEC=1.0D0/864.0D+2
      T_AGE=T_0*DAYSPERSEC
      T_CH_P =T_CH*DAYSPERSEC/2.0D0 ! FACTOR 2 IS FROM T_CH=T_LOAD/2
      DENSITY=G_C
      Q1_OPT=Q_1

C      WRITE(*,*)'T_AGE = ', T_AGE, ' [DAYS]'
C      WRITE(*,*)'T_CH  = ', T_CH_P, ' [DAYS]'
C      WRITE(*,*)'A/C   = ', A_TO_C
C      WRITE(*,*)'\RHO  = ', DENSITY, ' [KG/M^3]'

C THE UNIT IS MPA. THE FOLLOWING IS BASED ON ACI FORMULA:
      IF(F_C.EQ.0.0D0) THEN
        F_C=E28/4734.0D0*E28/4734.0D0 ! E28 IS ALWAYS THE INITIAL YOUNG.
      ELSE
        E28=4734.0D0*SQRT(F_C)
      END IF
C      WRITE(*,*)'F_C   = ', F_C, ' MPA'
      Q2  = 185.4D0*SQRT(DENSITY)*F_C**(-0.9D0)
      Q4  = 20.3D0*A_TO_C**(-0.7D0)
      IF(W_TO_C.EQ.0.0D0) W_TO_C=1.0D0/(F_C/22.8D0+0.535D0)
C      WRITE(*,*)'W/C   = ', W_TO_C
      Q1 = 0.60D+6/E28
      Q3 = 0.29D0*W_TO_C**4*Q2
      N=0.1D0
      M=0.5D0
      T=T_AGE+T_CH_P
C      WRITE(*,*) 'Q1     = ', Q1
C      WRITE(*,*) 'Q1_OPT = ', Q1_OPT
      R=1.7D0*T_AGE**0.12D0+8.0D0
      Z=T_AGE**(-M)*LOG(1.0D0+(T-T_AGE)**N)
      QF=1.0D0/(0.086D0*T_AGE**(2.D0/9.D0)+1.21D0*T_AGE**(4.0D0/9.0D0))
      Q=QF*Z/(Z**R+QF**R)**(1.0D0/R)
      J=Q1+Q2*Q+Q3*LOG(1.0D0+(T-T_AGE)**N)+Q4*LOG(T/T_AGE)
      J=Q1_OPT*(J/Q1)
      
C EP = 1.0D+6/(J-T_CH*J_DOT)
      J_DOT=N*(Q2*T**(-M)+Q3)/(T-T_AGE+(T-T_AGE)**(1.0D0-N))+Q4/T
      J_DOT=Q1_OPT*(J_DOT/Q1)
      EP=1.0D+6/(J-T_CH_P*J_DOT)

C ETA = 1.0D+6/J_DOT
      ETA=1.0D+6/J_DOT
      TAU1=ETA/EP

C THIS PART IS TO FACILITATE COMPUTATIONS
      IF(DT/TAU1.GE.1.0D+2) THEN
        C0=0.0D0
      ELSE
        C0=EXP(-DT/TAU1*DAYSPERSEC)
      END IF
      IF(DT/TAU1.LE.1.0D-3) THEN
        C1=1.0D0 - DT/TAU1/2.0D0*DAYSPERSEC
      ELSE
        C1=(1.0D0-C0)/(DT/TAU1*DAYSPERSEC)
      END IF
      EP=C1*EP
      YOUNG = EP
C      WRITE(*,*) 'J_DOT = ', J_DOT, ' 1/[MPA DAYS]'
C      WRITE(*,*) 'E_28 = E_STAT = ', E28, ' [MPA]'
C      WRITE(*,*) 'EP = ', EP, ' [MPA]'
C      WRITE(*,*) 'TAU1=', TAU1, ' [DAYS]'
C      WRITE(*,*) 'C1 = ', C1
C      WRITE(*,*) 'ETA = ', ETA, ' [MPA DAYS]'

      RETURN
      END SUBROUTINE M7F_CREEP

      SUBROUTINE INPUTPARAMS() 
      INCLUDE 'VABA_PARAM.INC'
      REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795D0
      INTEGER, PARAMETER :: NP=37, NVHM=5, NVHF=2, NVHI=NP*NVHM+NVHF
      COMMON /KCOMMON_VARS/YOUNG, POISSON,K_1, C0,DTIME,TH_DEL,UNIT_CONV  ! TH_DEL: FOR ELEMENT DELETION
      REAL*8 :: YOUNG, POISSON, K_1, C0, DTIME, TH_DEL, UNIT_CONV
      COMMON /KM7F_M7FIO/ ON_RATE_DEP_FRAC, ON_CREEP
      COMMON /KM7F_M7FIO_1/ C_R1, Q_1
      COMMON /KM7F_M7FIO_2/ T_0, T_CH, F_C, A_TO_C, G_C, W_TO_C, E28
      INTEGER      :: ON_RATE_DEP_FRAC, ON_CREEP
      REAL*8 :: C_R1, Q_1  
      REAL*8 :: T_0, T_CH, F_C, A_TO_C, G_C, W_TO_C, E28
      COMMON /KM7FBOUNDS_M7FIO/K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9,
     $     SV0_P, C_1, C_2, C_3, C_4, C_5,
     $     C_6, C_7, C_8, C_9, C_10, C_11,C_12,C_13,C_14,C_15,C_16,
     $     C_17,C_18,C_19,C_20,C_21,C_R0,C_R2
      REAL*8 :: K_2, K_3, K_4, K_5, K_6, K_7, K_8, K_9, SV0_P
      REAL*8 :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8,
     $     C_9, C_10
      REAL*8 :: C_11,C_12,C_13,C_14,C_15,C_16,C_17,C_18,
     $     C_19,C_20,C_21,C_R0,C_R2
      COMMON /KM7FBOUNDS_M7FIO_1/ C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1,
     $     C_6_2, C_6_M4, C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4,
     $     C_9_0, C_9_1, C_9_2, C_9_M4
      REAL*8 :: C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1, C_6_2, C_6_M4
      REAL*8 :: C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4
      REAL*8 :: C_9_0, C_9_1, C_9_2, C_9_M4
      COMMON /KM7FBOUNDS_M7FIO_2/SIG_FIB_0, P_1, P_2, P_5,P_6,XP_1,XP_2,
     $     XP_3, V_F, CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3, CF_4,
     $     CF_5, CF_4_4, AP1, BP1, CP2, DP2, EP2, D_1, D_2, D_3, D_4,
     $     D_5, D_6
      REAL*8 :: SIG_FIB_0, P_1, P_2, P_5, P_6
      REAL*8 :: XP_1, XP_2, XP_3, V_F
      REAL*8 :: CF_1, CF_2, CF_3, CF_1_1, CF_2_2, CF_3_3
      REAL*8 :: CF_4, CF_5, CF_4_4
      REAL*8 :: AP1, BP1, CP2, DP2, EP2
      REAL*8 :: D_1, D_2, D_3, D_4, D_5, D_6
      SAVE :: /KCOMMON_VARS/,/KM7F_M7FIO/,/KM7F_M7FIO_1/,/KM7F_M7FIO_2/,
     $     /KM7FBOUNDS_M7FIO/,/KM7FBOUNDS_M7FIO_1/,/KM7FBOUNDS_M7FIO_2/

C      CHARACTER(LEN=100) :: INP_PATH
C      CHARACTER(LEN=5)   :: PROCEED     
C!      INTEGER            :: I, I1, J1, JLOC
C      INTEGER            :: IOERR, OPENERR
C      NAMELIST/MASTER_INP_FILE/MASTER
C      NAMELIST/ADJUSTABLE_PARAMS/NP,YOUNG,POISSON,K_1,K_2,K_3,K_4,&
C           &C_5_0, C_5_1, C_5_M4, C_6_0, C_6_1, C_6_2, C_6_M4, &
C           &C_7_0, C_7_1, C_7_M4, C_8_0, C_8_1, C_8_M4, &
C           &C_9_0, C_9_1, C_9_2, C_9_M4, D_1, D_2, D_3, D_4, D_5, D_6 
C      NAMELIST/ADJUSTABLE_COUPLING_PARAMS/K_5,K_6,K_7,K_8,K_9,SV0_P
C      NAMELIST/FIXED_PARAMS/C_2,C_3,C_4,C_10,&
C           &C_11,C_12,C_16,C_17,C_18,C_19,C_20,C_21
C      NAMELIST/ADJUSTABLE_FIBER_PARAMS/SIG_FIB_0,P_1,P_2,P_5,P_6,&
C           &XP_1,XP_2,XP_3,V_F,CF_1,CF_2,CF_3,CF_1_1,CF_2_2,CF_3_3,&
C           &CF_4,CF_4_4,CF_5
C!      NAMELIST/ROTATION_PARAMS/R_1,R_2 
C      NAMELIST/RATE_DEP_FRAC_PARAMS/ON_RATE_DEP_FRAC,C_R1,C_R2 
C      NAMELIST/CREEP_PARAMS/ON_CREEP,T_0,T_CH,F_C,A_TO_C,G_C,W_TO_C
C      NAMELIST/OUTFILES/OUTFILE, FAILFILE, SYSFILE
C      NAMELIST/LOADPARAMS/NSTEP,MITER,TOLER 

C THIS IS A MASTER INPUT FILE. IT ONLY HAS INFORMATION ABOUT THE MATERIAL
C PARAMETERS. IT DOES NOT HAVE ANY INFORMATION ABOUT THE RUN PARAMETERS 
C SUCH AS BOUNDARY CONDITIONS, NUMBER OF LOADING STEPS, OUTPUT FILES, ETC.
C THE RUN PARAMETERS ARE SPECIFIED FOR EACH SUBFIGURE IN A DIFFERENT INP FILE.
C      WRITE(7,*) "INSIDE INPUTPARAMS...."

      TH_DEL = 0.005D10 !THRESHOLD TENSILE STRAIN FOR ELEMENT DELETION
C ---------------------------------------------------------
C WHAT THE PARAMETER UNIT_CONV DOES IS, IT MAKES SURE THAT
C ON ENTRY TO M7FMATERIAL, THE STRESSES ARE IN MPA AND ON
C EXIT FROM M7FMATERIAL, THE STRESSES AND FORCES ARE IN THE 
C UNITS DICTATED BY THE MESH UNITS.
C IN M7F, THE STRESSES ARE IN MPA. THIS REQUIRES THE MESH
C TO BE IN MM, FORCES IN N AND MODULI IN MPA. IN THAT CASE,
C UNIT_CONV=1.
C HOWEVER, IF THE MESH IS IN M, FORCES IN N AND MODULI ARE 
C IN PA, UNIT_CONV SHOULD BE SET TO 1.D-6 IN THE SUBROUTINE 
C INPUTPARAMS(). 
C OTHER UNIT CONVERSIONS CAN ALSO BE HANDLED THROUGH SETTING
C UNIT_CONV TO THE RIGHT PARAMETER. FOR EXAMPLE, IF THE MESH
C IS IN INCHES, AND THE FORCES IN POUNDS, UNITCONV MUST BE
C SET TO 6.89475729D-3.
C ----------------------------------------------------------
C COMMENT OUT ALL UNIT_CONV LINES EXCEPT THE ONE THAT APPLIES.
C      UNIT_CONV=1.D-6             !MULTIPLY WITH 1.D-6 TO CONVERT FROM MPA TO PA
      UNIT_CONV=1.                !MULTIPLY WITH 1. TO CONVERT FROM MPA TO MPA
C      UNIT_CONV=6.89475729D-3     !MULTIPLY WITH 6.89475729D-3 TO CONVERT FROM PSI TO MPA

      YOUNG  = 24541.D0 !25000.0D0
      POISSON = 0.23000D0
      K_1    = 171.80D-6 !110.00D-6!
      K_2    = 40.000D0
      K_3    = 20.3D0 ! 30.0000D0
      K_4    = 36.6D0 ! 65.000D0
      C_5_0 = 1.3D-2! 1.2D-2 !1.2D-2
      C_5_1 = 4.0D0 ! 4.0D0 !#4.0D0 !7.0D0 !11.0D0
      C_7_0  = 1.2D-2 !1.3D-2
      C_7_1 = 35.0D2 ! 35.0D2 !#35.0D2
      C_8_0  = 1.2D-2 !1.2D-2
      C_8_1 = 20.0D0 ! 20.0D0 !#20.0D0 !7.0D0!11.0D0
      C_6_0 = 4.0D2 ! 4.0D2 !#4.0D2 
      C_6_1  = 4.0D1 ! 4.0D1
      C_6_2  = 13.0D0
      C_9_0 = 4.0D2 ! 4.0D2 !#4.0D2 
      C_9_1  = 4.0D1 ! 4.0D1
      C_9_2  = 13.0D0
      C_5_M4 = 3.0D0    ! 2.55D0 !#2.55D0 !3.05D0
      C_6_M4 = 1.30D0 !1.40D0
      C_7_M4 = 10.D-1   !# 10.D0 !200D0 ! 90.0D0 !200.D0!#200.0D0 !70.0D0
      C_8_M4 = 8.0D0  !8.0D0 !2.55D0 !3.05D0
      C_9_M4 = 0.D-1 !#1.30D0!#1.40D0
C     ! FIBER EFFECT ON THE NORMAL BOUNDARY
      D_1=0.095 !2.94D-9
      D_2=35.   !10940.171D0 ! -0.384D0
      D_3=1.7   !578.348D0   ! -3.51D-5
      D_4=1.7D0 !#0.55  !0.29D0      ! -0.0203D0
      D_5 = 1.D3
      D_6 = 25.D0
     
      SV0_P=250.D0
      K_8=3.D0
      K_9=5.D-1
      K_5=1.0D0
      K_6=1.0D-4 !1.0D-4 ! 5.D-3 ! 1.0D-4
      K_7=1.8D0 !1.7D0 !1.6D0  ! 1.D0  ! 1.6D-0
     
      C_2 = 7.1D0 !2.76D0         !NORMAL BOUND. PARAM. ......... C_2
      C_3 = 4.D0 ! 4.00D0         !NORMAL PLASTICITY  ........... C_3
      C_4 = 250.D0 !20.00D0        !STRAIN RATIO: NORMAL/VOL ..... C_4
      C_10 = 3.3D-1 !3.3D-1 !7.30D-1       !FRIC.B. INITIAL SLOPE ........ C_10 = 7.30D-1
      C_11 = 5.D-1 !5.D-1 !9.00D-1       !FRIC.B.\SIG_N INTER.@\SIG_V=0. C_11 = 2.00D-1
      C_12 = 7.00D3 !#7.00D+3       !FRIC.B.SPEED\SIG_N GOES ZERO.. C_12 = 7.00D+3
      C_16 = 10.0D0 !2.00D-2 ! 19.09.2011       
      C_17 = 1.00D-2       
      C_18 = 4.D3 !6.D3 ! UNLOADING/RELOADING IN COMPRESSION
      C_19 = 4.5D3 !6.3D3 !19.09.2011 ! UNLOADING/RELOADING IN TENSION        
      C_20 = 3.D2 !3.D2 ! UNLOADING/RELOADING IN COMPRESSION
      C_21 = 6.D1  ! UNLOADING/RELOADING IN COMPRESSION      
     
      SIG_FIB_0 = -4.5D0 ! M5F:   2.5D0
      P_1    = 0.0D0     ! M5F: 290.0D0
      P_2    = 300.D-3   ! M5F:  22.0D0
      P_5    = 0.0D0     ! M5F:1500.0D0
      P_6    = 0.0D0     ! M5F: 200.0D0
      XP_1   = 3.34D0    ! 8.4D-4   ! >= XP_3/P_2 + XP_2 MUST BE SATISFIED!
      XP_2   = 0.        ! HORIZONTAL SHIFT OF FIBER LAW
      XP_3   = 1.0D0     ! POWER ON THE FIBER LAW.
      V_F    = 0.0D0     ! FIBER VOLUME FRACTION
      CF_1=104.0D-2
      CF_2=54.0D+1
      CF_3=104.0D-2
      CF_1_1=3.5D+1
      CF_2_2=1.5D+1
      CF_3_3=3.5D+1
      CF_4_4=4.4D2 
      CF_4  = 3.2D-2 ! 23.2D-2
      CF_5  = 5.0D-2 ! 5.0D-2
     
      ON_RATE_DEP_FRAC = 0
      C_R1 = 4.0D-6
      C_R2 = 22.0D-3
     
      ON_CREEP = 0
      T_0  = 241.920D+4
      T_CH = 180.000D+0
      F_C  =  48.000D+0
      A_TO_C = 7.000D+0
      G_C  = 304.340D+0
      W_TO_C = 0.000D+0
     
      Q_1=21.00D0        !Q1 (TO BE OPTIMIZED)(1/MPA)... Q_1
      C0=1.0D0           !WHEN THERE IS NO CREEP........ C0
     
C     ! YOUNG WILL CHANGE IF DT CHANGES (IN SUBROUTINE CREEP),
C     ! BUT E28 BELOW IS FIXED, AND IS COMMUNICATED TO SUBROUTINE
C     ! CREEP. 
      E28 = YOUNG
      IF(T_CH.EQ.0.0D0) THEN
         WRITE(*,*) 'ERROR: DURATION OF LOADING IS NOT ALLOWED TO'//
     $        ' BE ZERO.'
         WRITE(*,*) 'IF YOU WANT TO TURN OFF CREEP EFFECTS, USE THE '//
     $        'PROVIDED SWITCH.'
         STOP
      END IF
      RETURN
      END SUBROUTINE INPUTPARAMS
