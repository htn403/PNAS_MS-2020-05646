C
C ADJUSTING PARAMETERS TO FIT GIVEN MATERIAL TEST DATA
C 
C Using the free parameters to fit given data may sometimes lead to less
C than perfect fits of test data. In that case, it is necessary to
C adjust some of the fixed, or hard-to-change, parameters c_i. This is
C quite challenging. Improving the fit of one type of test by changing
C one parameter, c_j, may spoil the fit of another type of test. It may
C also necessitate changing several parameters c_j, which is a rather
C difficult nonlinear optimization problem. Attempts by several
C researchers to automatically adjust these parameters using
C optimization software has been unsuccessful in the past because of the
C complexity of the response. To succeed requires considerable
C experience. In that case, you can obtain assistance from
C Prof. F.C. Caner ( <ferhun.caner@upc.edu>, Tel. +34659816715, Skype
C ferhun.caner). He has ample experience in the fitting of test data
C with M7 using the finite element method.
C 
C The program is based on the following papers, which can be freely
C downloaded as 527.pdf, 528.pdf, 519.pdf and 547.pdf, respectively,
C from ZP Bazant's website:
C http://www.civil.northwestern.edu/people/bazant/PDFs/Papers/
C 
C [1] Caner, F.C., and Ba\v zant, Z.P. (2013). ``Microplane model M7 for
C     plain concrete: I. formulation." {\em ASCE J. of Engrg. Mechanics}
C     139 (12), Dec., 1714--1723.
C 
C [2] Caner, F.C., and Ba\v zant, Z.P. (2013). ``Microplane model M7 for
C     plain concrete: II. calibration and verification." {\em ASCE J. of
C     Engrg. Mechanics} 139 (12), Dec., 1724--1735.
C 
C [3] Caner, F.C., Ba\v zant, Z.P., and Wendner, R. (2013). ``Microplane
C     model M7f for fiber reinforced concrete." {\em Engrg. Fracture
C     Mechanics} 105, 41--57.
C 
C [4] Kirane, K., and Ba\v zant, Z.P. (2014). ``Microplane damage model
C     for fatigue of quasibrittle materials: Sub-critical crack growth,
C     lifetime and residual strength." {\em International Journal of
C     Fatigue} 70, 93--105.
C      
      subroutine vumat( 
     $     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     $     stepTime, totalTime, dt, cmname, coordMp, charLength,
     $     props, density, strainInc, relSpinInc,
     $     tempOld, stretchOld, defgradOld, fieldOld,
     $     stressOld, stateOld, enerInternOld, enerInelasOld,
     $     tempNew, stretchNew, defgradNew, fieldNew,
     $     stressNew, stateNew, enerInternNew, enerInelasNew )
      include 'vaba_param.inc'
      dimension props(nprops), density(nblock),
     $coordMp(nblock,*),
     $charLength(nblock), strainInc(nblock,ndir+nshr),
     $relSpinInc(nblock,nshr), tempOld(nblock),
     $stretchOld(nblock,ndir+nshr), defgradOld(nblock,ndir+nshr+nshr),
     $fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     $stateOld(nblock,nstatev), enerInternOld(nblock),
     $enerInelasOld(nblock), tempNew(nblock),
     $stretchNew(nblock,ndir+nshr), defgradNew(nblock,ndir+nshr+nshr), 
     $fieldNew(nblock,nfieldv),
     $stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     $enerInternNew(nblock), enerInelasNew(nblock)
      character*80 cmname

      IF(CMNAME(1:11).EQ."M7FMATERIAL") THEN
         CALL M7fMaterial( 
     $ nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     $ stepTime, totalTime, dt, cmname, coordMp, charLength,
     $ props, density, strainInc, relSpinInc,
     $ tempOld, stretchOld, defgradOld, fieldOld,
     $ stressOld, stateOld, enerInternOld, enerInelasOld,
     $ tempNew, stretchNew, defgradNew, fieldNew,
     $ stressNew, stateNew, enerInternNew, enerInelasNew )
      END IF
      return
      END subroutine vumat

C +--------------------------------------------------------------------+
C |                 SUBROUTINE C_NORM_ELASTIC                          |
C +--------------------------------------------------------------------+
      subroutine c_norm_elastic(ef,def,sn0,eps_N0_pos,eps_N0_neg,sv0,
     $snf_no_split,E_N,zeta) 
      include 'vaba_param.inc'
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv, zeta
      common /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     sv0_p, c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21,C_R0,C_R2
      real :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, sv0_p
      real :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      real :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21,C_R0,C_R2
      common /KM7fBounds_M7fIO_1/ c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1,
     $     c_6_2, c_6_M4, c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4,
     $     c_9_0, c_9_1, c_9_2, c_9_M4
      real :: c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_6_M4
      real :: c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4
      real :: c_9_0, c_9_1, c_9_2, c_9_M4
      common /KM7fBounds_M7fIO_2/sig_fib_0, p_1, p_2, p_5, p_6,xp_1,xp_2,
     $     xp_3, V_f, cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3, cf_4,
     $     cf_5, cf_4_4, ap1, bp1, cp2, dp2, ep2, d_1, d_2, d_3, d_4,
     $     d_5, d_6
      real :: sig_fib_0, p_1, p_2, p_5, p_6
      real :: xp_1, xp_2, xp_3, V_f
      real :: cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3
      real :: cf_4, cf_5, cf_4_4
      real :: ap1, bp1, cp2, dp2, ep2
      real :: d_1, d_2, d_3, d_4, d_5, d_6
      save :: /KCommon_Vars/, /KM7fBounds_M7fIO/, /KM7fBounds_M7fIO_1/,
     $     /KM7fBounds_M7fIO_2/
      real, dimension(1:np) :: ef, def, sn0
      real, dimension(1:np) :: eps_N0_pos, eps_N0_neg
      real :: sv0
      real, dimension(1:np) :: snf_no_split 
      real, dimension(1:np) :: E_N

      real :: E_N_0
      real :: weight_factor, t0, t1, t2, t3, t4, tsum, fzeta
      integer :: isize, i 
      E_N_0= young / (1.0d0 - 2.0d0*poisson)

      isize=size(ef)

      do i=1,isize
         if ((sn0(i) > 0.d0)) then 
            weight_factor=4.0d1
            t0=1.d0
            t1=(weight_factor*zeta)**2.d0
            t2=(weight_factor*zeta)**4.d0
            t3=(weight_factor*zeta)**6.d0
            t4=(weight_factor*zeta)**8.d0
            tsum=t0+t1+t2
            fzeta=1.d0/tsum
            E_N(i) = E_N_0*fzeta*exp(-c_19*eps_N0_pos(i))
            if (sn0(i) > E_N_0*ef(i) .and. sn0(i)*def(i)<0.0d0) then 
               E_N(i)=E_N_0                
            end if                         
         else 
            E_N(i) = E_N_0*(exp(-c_20*abs(eps_N0_neg(i))/
     $           (1.d0+c_18*max(-sv0,0.d0)/E_N_0))+
     $           c_21*max(-sv0,0.d0)/E_N_0)
         end if
      end do

      snf_no_split = sn0*C0 + E_N*def

      return
      end subroutine c_norm_elastic

C +--------------------------------------------------------------------+
C |                         SUBROUTINE C_NORM                          |
C +--------------------------------------------------------------------+
      subroutine c_norm(ef,sv0,R_N,snb) 
      include 'vaba_param.inc'
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv
      common /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     sv0_p, c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21,C_R0,C_R2
      real :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, sv0_p
      real :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      real :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21,C_R0,C_R2
      common /KM7fBounds_M7fIO_1/ c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1,
     $     c_6_2, c_6_M4, c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4,
     $     c_9_0, c_9_1, c_9_2, c_9_M4
      real :: c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_6_M4
      real :: c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4
      real :: c_9_0, c_9_1, c_9_2, c_9_M4
      common /KM7fBounds_M7fIO_2/sig_fib_0, p_1, p_2, p_5,p_6,xp_1,xp_2,
     $     xp_3, V_f, cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3, cf_4,
     $     cf_5, cf_4_4, ap1, bp1, cp2, dp2, ep2, d_1, d_2, d_3, d_4,
     $     d_5, d_6
      real :: sig_fib_0, p_1, p_2, p_5, p_6
      real :: xp_1, xp_2, xp_3, V_f
      real :: cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3
      real :: cf_4, cf_5, cf_4_4
      real :: ap1, bp1, cp2, dp2, ep2
      real :: d_1, d_2, d_3, d_4, d_5, d_6
      save :: /KCommon_Vars/, /KM7fBounds_M7fIO/, /KM7fBounds_M7fIO_1/,
     $     /KM7fBounds_M7fIO_2/
      real, dimension(1:np) :: ef 
      real, dimension(1:np) :: snb 
      real :: sv0, R_N 
      real, dimension(1:np)  :: aux1
      real :: C_N0, RATE_FACTOR, fstar, beta_N
      real :: E_V, eb_N0, eb_NK, eb_N
      integer :: allocstat, isize

      isize = size(ef)

      E_V    = young/(1.0d0-2.0d0*poisson)
      eb_N0  = c_3*k_1 
      eb_NK  = c_4 

      c_1 = d_1*tanh(d_2*V_f - d_3)+d_4*exp(-max(-sv0-d_6,0.d0)/E_V*d_5) 

      if (sv0.lt.0.0d0) then 
         eb_N = eb_N0 - eb_NK / E_V * sv0  
      else 
         eb_N = eb_N0 
      end if       

      fstar  = k_1*young*c_1 
      beta_N = c_2*c_1*k_1 

      C_N0=C_R2 

      aux1 = max(ef-beta_N,0.0d0)/eb_N
      snb  = fstar*exp(-aux1)

      RATE_FACTOR = C_N0*R_N
      snb = snb * (1.d0 + RATE_FACTOR)

      return 
      end subroutine c_norm

! +--------------------------------------------------------------------+
! |                       SUBROUTINE C_NORM_FIB                        |
! +--------------------------------------------------------------------+
      subroutine c_norm_fib(ef,R_N,snb_fib) 
      include 'vaba_param.inc'
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv
      common /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     sv0_p, c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21,C_R0,C_R2
      real :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, sv0_p
      real :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      real :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21,C_R0,C_R2
      common /KM7fBounds_M7fIO_1/ c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1,
     $     c_6_2, c_6_M4, c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4,
     $     c_9_0, c_9_1, c_9_2, c_9_M4
      real :: c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_6_M4
      real :: c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4
      real :: c_9_0, c_9_1, c_9_2, c_9_M4
      common /KM7fBounds_M7fIO_2/sig_fib_0, p_1, p_2, p_5,p_6,xp_1,xp_2,
     $     xp_3, V_f, cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3, cf_4,
     $     cf_5, cf_4_4, ap1, bp1, cp2, dp2, ep2, d_1, d_2, d_3, d_4,
     $     d_5, d_6
      real :: sig_fib_0, p_1, p_2, p_5, p_6
      real :: xp_1, xp_2, xp_3, V_f
      real :: cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3
      real :: cf_4, cf_5, cf_4_4
      real :: ap1, bp1, cp2, dp2, ep2
      real :: d_1, d_2, d_3, d_4, d_5, d_6
      save :: /KCommon_Vars/, /KM7fBounds_M7fIO/, /KM7fBounds_M7fIO_1/,
     $     /KM7fBounds_M7fIO_2/
      real, dimension(1:np) :: ef 
      real, dimension(1:np) :: snb_fib 

      real :: R_N 

      real :: C_N0, RATE_FACTOR 

      real :: e_N_max, e_N_0, e_N_hat, eps_shift, pow 
      integer :: isize, i 
      
      isize = size(ef)

      eps_shift = xp_2
      pow = xp_3

      e_N_max = pow/p_2 + eps_shift

      e_N_0 = xp_1  

      e_N_hat = sig_fib_0 

      C_N0=C_R2 

      do i=1, isize
         if ( ef(i)/k_1 < e_N_max ) then
            snb_fib(i) = young*p_1*k_1*
     $           max(ef(i)/k_1 - eps_shift,0.0d0)**pow*
     $           exp(-p_2*max(ef(i)/k_1 - eps_shift,0.0d0)) 
         else if ((ef(i)/k_1 >= e_N_max) .and. (ef(i)/k_1 < e_N_0)) then
            snb_fib(i) = young*p_1*k_1*(pow/p_2)**pow*exp(-1.0d0*pow)
         else if (ef(i)/k_1 >= e_N_0) then 
            snb_fib(i) = young*p_1*k_1*
     $           (ef(i)/k_1 - e_N_0  + e_N_max - eps_shift)**pow*
     $           exp(-p_2*(ef(i)/k_1 - e_N_0 + e_N_max - eps_shift))  

         end if
      end do

      RATE_FACTOR = C_N0*R_N
      snb_fib = snb_fib * (1.d0 + RATE_FACTOR)
      return
      end subroutine c_norm_fib

! +--------------------------------------------------------------------+
! |                      SUBROUTINE C_DEV                              |
! +--------------------------------------------------------------------+
      subroutine c_dev(ded, ed0, dev, ev0, sv0, sd0, sdf, C_d, R_D,
     $     sdneg, sdpos) 
      include 'vaba_param.inc'
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv
      common /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     sv0_p, c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21,C_R0,C_R2
      real :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, sv0_p
      real :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      real :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21,C_R0,C_R2
      common /KM7fBounds_M7fIO_1/ c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1,
     $     c_6_2, c_6_M4, c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4,
     $     c_9_0, c_9_1, c_9_2, c_9_M4
      real :: c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_6_M4
      real :: c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4
      real :: c_9_0, c_9_1, c_9_2, c_9_M4
      common /KM7fBounds_M7fIO_2/sig_fib_0, p_1, p_2,p_5,p_6,xp_1,xp_2,
     $     xp_3, V_f, cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3, cf_4,
     $     cf_5, cf_4_4, ap1, bp1, cp2, dp2, ep2, d_1, d_2, d_3, d_4,
     $     d_5, d_6
      real :: sig_fib_0, p_1, p_2, p_5, p_6
      real :: xp_1, xp_2, xp_3, V_f
      real :: cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3
      real :: cf_4, cf_5, cf_4_4
      real :: ap1, bp1, cp2, dp2, ep2
      real :: d_1, d_2, d_3, d_4, d_5, d_6
      save :: /KCommon_Vars/, /KM7fBounds_M7fIO/, /KM7fBounds_M7fIO_1/,
     $     /KM7fBounds_M7fIO_2/
      real, dimension(1:np) :: ded, ed0, sd0
      real, dimension(1:np) :: sdf, C_d

      real, dimension(1:np) :: sdneg, sdpos

      real :: dev, ev0, sv0, R_D
      real :: usi_C, usi_T, par5, par6, E_s, Cd0, Ev
      real, dimension(1:np) :: edf, sde, sdlower, sdupper
      real, dimension(1:np) :: sdlower_fib, sdupper_fib
      real :: RATE_FACTOR_C, RATE_FACTOR_T, C_DC0, C_DT0, evf, E_v

      real :: f_c0, E_0, f_cp, c_40, beta_15, beta_16, beta_17, beta_18,
     $     beta_19 
      logical :: l0,l1,l2,l3,l4
      integer :: j, isize, allocstat
      
      isize = size(ded)

      evf = ev0 + dev

      Ev=young/(1.0d0 - 2.0d0*poisson)

      f_c0=15.08d0
      E_0=20000.d0
      c_40=1.0d+0

      f_cp=90.3
      beta_15=c_5_1*exp(-c_40*(f_cp/young-f_c0/E_0))
      beta_16=c_8_1*exp(-c_40*(f_cp/young-f_c0/E_0))
      beta_17=c_7_1*exp(-c_40*(f_cp/young-f_c0/E_0))
      beta_18=c_6_0*exp(-c_40*(f_cp/young-f_c0/E_0))
      beta_19=c_9_0*exp(-c_40*(f_cp/young-f_c0/E_0))

      c_5 = (beta_15*tanh(c_5_0*max(-evf,0.0d0)/k_1))*
     $     (1.0d0 + cf_1*tanh(cf_1_1*V_f))+cf_1*tanh(cf_1_1*V_f)+c_5_M4

      c_8 = (beta_16*tanh(c_8_0*max(-evf,0.0d0)/k_1))*
     $     (1.0d0 + cf_3*tanh(cf_3_3*V_f))+cf_3*tanh(cf_3_3*V_f)+c_8_M4

      c_7 = beta_17*tanh(c_7_0*max(-evf,0.0d0)/k_1)+
     $     cf_2*tanh(cf_2_2*V_f) + c_7_M4 

      if (beta_18*max(-evf/k_1-c_6_1,0.0d0)>=log(c_6_2)) then
         c_6 = c_6_M4*c_6_2
      else
         c_6 = c_6_M4*exp(beta_18*max(-evf/k_1-c_6_1,0.0d0))
      end if
      if (beta_19*max(-evf/k_1-c_9_1,0.0d0)>=log(c_9_2)) then
         c_9 = c_9_M4*c_9_2
      else
         c_9 = c_9_M4*exp(beta_19*max(-evf/k_1-c_9_1,0.0d0))
      end if

      par5 = k_1*young*c_8
      par6 = young*c_5*k_1
      Cd0=young/(1.d0+poisson)*(1.d0-4.d0*poisson)/(1.d0-2.d0*poisson)
      E_v  = young/(1.d0-2*poisson)

      C_DC0  = C_R2

      C_DT0  = C_R2

      C_d=Cd0

      edf = ed0+ded 

      call dev_comp(edf, sdlower) 

      call dev_tens(edf, sdupper)
      call dev_tens_fib(edf, sdupper_fib)

      RATE_FACTOR_C = C_DC0*R_D
      RATE_FACTOR_T = C_DT0*R_D

c      sdlower = sdlower * (1.d0 + RATE_FACTOR_C)      

c      sdupper = sdupper * (1.d0 + RATE_FACTOR_T)

      sdneg=sdlower
      sdpos=sdupper

      return
      end subroutine c_dev

! +--------------------------------------------------------------------+
! |                        SUBROUTINE DEV_COMP                         |
! +--------------------------------------------------------------------+
      subroutine dev_comp(edf, sdlower) 
      include 'vaba_param.inc'
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv
      common /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     sv0_p, c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21,C_R0,C_R2
      real :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, sv0_p
      real :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      real :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21,C_R0,C_R2
      common /KM7fBounds_M7fIO_1/ c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1,
     $     c_6_2, c_6_M4, c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4,
     $     c_9_0, c_9_1, c_9_2, c_9_M4
      real :: c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_6_M4
      real :: c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4
      real :: c_9_0, c_9_1, c_9_2, c_9_M4
      common /KM7fBounds_M7fIO_2/sig_fib_0, p_1, p_2, p_5,p_6,xp_1,xp_2,
     $     xp_3, V_f, cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3, cf_4,
     $     cf_5, cf_4_4, ap1, bp1, cp2, dp2, ep2, d_1, d_2, d_3, d_4,
     $     d_5, d_6
      real :: sig_fib_0, p_1, p_2, p_5, p_6
      real :: xp_1, xp_2, xp_3, V_f
      real :: cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3
      real :: cf_4, cf_5, cf_4_4
      real :: ap1, bp1, cp2, dp2, ep2
      real :: d_1, d_2, d_3, d_4, d_5, d_6
      save :: /KCommon_Vars/, /KM7fBounds_M7fIO/, /KM7fBounds_M7fIO_1/,
     $     /KM7fBounds_M7fIO_2/
      real, dimension(1:np) :: edf

      real, dimension(1:np) :: sdlower
      real :: par1, par2, par5

      par1 = c_9*c_8*k_1
      par2 = c_7*k_1 
      par5 = k_1*young*c_8

      sdlower=-par5/(1.d0+(max(-edf-par1,0.0d0)/par2)**2.0d0) 

      return 
      end subroutine dev_comp 
C
C
C
C +--------------------------------------------------------------------+
C |                         SUBROUTINE DEV_TENS                        |
C +--------------------------------------------------------------------+

      subroutine dev_tens(edf, sdupper) 
      include 'vaba_param.inc'
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv
      common /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     sv0_p, c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21,C_R0,C_R2
      real :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, sv0_p
      real :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      real :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21,C_R0,C_R2
      common /KM7fBounds_M7fIO_1/ c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1,
     $     c_6_2, c_6_M4, c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4,
     $     c_9_0, c_9_1, c_9_2, c_9_M4
      real :: c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_6_M4
      real :: c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4
      real :: c_9_0, c_9_1, c_9_2, c_9_M4
      common /KM7fBounds_M7fIO_2/sig_fib_0, p_1, p_2, p_5,p_6,xp_1,xp_2,
     $     xp_3, V_f, cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3, cf_4,
     $     cf_5, cf_4_4, ap1, bp1, cp2, dp2, ep2, d_1, d_2, d_3, d_4,
     $     d_5, d_6
      real :: sig_fib_0, p_1, p_2, p_5, p_6
      real :: xp_1, xp_2, xp_3, V_f
      real :: cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3
      real :: cf_4, cf_5, cf_4_4
      real :: ap1, bp1, cp2, dp2, ep2
      real :: d_1, d_2, d_3, d_4, d_5, d_6
      save :: /KCommon_Vars/, /KM7fBounds_M7fIO/, /KM7fBounds_M7fIO_1/,
     $     /KM7fBounds_M7fIO_2/
      real, dimension(1:np) :: edf

      real, dimension(1:np) :: sdupper
      real :: par2, par3, par4, par6

      par2 = c_7*k_1 
      par3 = c_20
      par4 = c_6*c_5*k_1
      par6 = young*c_5*k_1

      sdupper=par6/(1.0d0+(max(edf-par4,0.0d0)/par2/par3)**2.0d0)

      return 
      end subroutine dev_tens 

! +--------------------------------------------------------------------+
! |                         SUBROUTINE DEV_TENS_FIB                    |
! +--------------------------------------------------------------------+
      subroutine dev_tens_fib(edf, sdupper_fib) 
      include 'vaba_param.inc'
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv
      common /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     sv0_p, c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21,C_R0,C_R2
      real :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, sv0_p
      real :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      real :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21,C_R0,C_R2
      common /KM7fBounds_M7fIO_1/ c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1,
     $     c_6_2, c_6_M4, c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4,
     $     c_9_0, c_9_1, c_9_2, c_9_M4
      real :: c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_6_M4
      real :: c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4
      real :: c_9_0, c_9_1, c_9_2, c_9_M4
      common /KM7fBounds_M7fIO_2/sig_fib_0, p_1, p_2, p_5,p_6,xp_1,xp_2,
     $     xp_3, V_f, cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3, cf_4,
     $     cf_5, cf_4_4, ap1, bp1, cp2, dp2, ep2, d_1, d_2, d_3, d_4,
     $     d_5, d_6
      real :: sig_fib_0, p_1, p_2, p_5, p_6
      real :: xp_1, xp_2, xp_3, V_f
      real :: cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3
      real :: cf_4, cf_5, cf_4_4
      real :: ap1, bp1, cp2, dp2, ep2
      real :: d_1, d_2, d_3, d_4, d_5, d_6
      save :: /KCommon_Vars/, /KM7fBounds_M7fIO/, /KM7fBounds_M7fIO_1/,
     $     /KM7fBounds_M7fIO_2/
      real, dimension(1:np) :: edf
      real, dimension(1:np) :: sdupper_fib

      sdupper_fib=young*p_5*k_1*max(edf/k_1,0.0d0)*
     $     exp(-p_6*max(edf/k_1,0.0d0))

      return 
      end subroutine dev_tens_fib 

      subroutine c_shear_tens(sig_o, s_0, fsp_0)
      include 'vaba_param.inc'
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv
      common /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     sv0_p, c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21,C_R0,C_R2
      real :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, sv0_p
      real :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      real :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21,C_R0,C_R2
      common /KM7fBounds_M7fIO_1/ c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1,
     $     c_6_2, c_6_M4, c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4,
     $     c_9_0, c_9_1, c_9_2, c_9_M4
      real :: c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_6_M4
      real :: c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4
      real :: c_9_0, c_9_1, c_9_2, c_9_M4
      common /KM7fBounds_M7fIO_2/sig_fib_0, p_1, p_2, p_5,p_6,xp_1,xp_2,
     $     xp_3, V_f, cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3, cf_4,
     $     cf_5, cf_4_4, ap1, bp1, cp2, dp2, ep2, d_1, d_2, d_3, d_4,
     $     d_5, d_6
      real :: sig_fib_0, p_1, p_2, p_5, p_6
      real :: xp_1, xp_2, xp_3, V_f
      real :: cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3
      real :: cf_4, cf_5, cf_4_4
      real :: ap1, bp1, cp2, dp2, ep2
      real :: d_1, d_2, d_3, d_4, d_5, d_6
      save :: /KCommon_Vars/, /KM7fBounds_M7fIO/, /KM7fBounds_M7fIO_1/,
     $     /KM7fBounds_M7fIO_2/
      real :: sig_o, s_0
      real :: fsp_0
      real :: c_6_p

      c_6_p = c_10
      fsp_0 = c_6_p * max(sig_o, 0.0d0) / (1.d0 + c_6_p / s_0 *
     $     max(sig_o, 0.0d0)) 
      return
      end subroutine c_shear_tens

! +--------------------------------------------------------------------+
! |                          SUBROUTINE C_SHEAR2                       |
! +--------------------------------------------------------------------+
      subroutine c_shear2(eps_L, eps_M, C_d, snf, R_S, del, dem, sl0,
     $     sm0, slf, smf, eps_V) 
      include 'vaba_param.inc'
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv
      common /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     sv0_p, c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21,C_R0,C_R2
      real :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, sv0_p
      real :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      real :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21,C_R0,C_R2
      common /KM7fBounds_M7fIO_1/ c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1,
     $     c_6_2, c_6_M4, c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4,
     $     c_9_0, c_9_1, c_9_2, c_9_M4
      real :: c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_6_M4
      real :: c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4
      real :: c_9_0, c_9_1, c_9_2, c_9_M4
      common /KM7fBounds_M7fIO_2/sig_fib_0, p_1, p_2, p_5,p_6,xp_1,xp_2,
     $     xp_3, V_f, cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3, cf_4,
     $     cf_5, cf_4_4, ap1, bp1, cp2, dp2, ep2, d_1, d_2, d_3, d_4,
     $     d_5, d_6
      real :: sig_fib_0, p_1, p_2, p_5, p_6
      real :: xp_1, xp_2, xp_3, V_f
      real :: cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3
      real :: cf_4, cf_5, cf_4_4
      real :: ap1, bp1, cp2, dp2, ep2
      real :: d_1, d_2, d_3, d_4, d_5, d_6
      save :: /KCommon_Vars/, /KM7fBounds_M7fIO/, /KM7fBounds_M7fIO_1/,
     $     /KM7fBounds_M7fIO_2/
      real, dimension(1:np)  :: eps_L, eps_M, C_d, snf, del,
     $     dem, sl0, sm0 
      real  :: R_S, eps_V
      real, dimension(1:np) :: slf, smf  
      real :: E_T, c_12_p, c_13_p, c_6_p, s_0
      real, dimension(1:np) :: C_t, fsp
      real, dimension(1:np) :: stf, fsp_0, epsf, depsf
      real :: RATE_FACTOR_1, C_S0, sig_o
      integer :: isize, allocstat, j
      
      isize = size(snf)

      E_T=young/(1.d0+poisson)*(1.d0-4.d0*poisson)/(1.d0-2.d0*poisson) 
      c_6_p  = c_10
      c_12_p = c_11 
      c_13_p = c_12*3.371d-4

      s_0    = k_1*k_2*E_T

      C_S0 = C_R2

      RATE_FACTOR_1 = C_S0*R_S  

      sig_o = max(E_T * k_1 *
     $     (c_12_p - c_13_p*max(eps_V, 0.0d0)/k_1),0.0d0)

      fsp = c_6_p * max(-snf + sig_o, 0.0d0) / (1.d0 + c_6_p / s_0 *
     $     max(-snf + sig_o, 0.0d0))! * (1.d0 + RATE_FACTOR_1)

      epsf = sqrt(eps_L*eps_L + eps_M*eps_M)
      depsf = sqrt(del*del + dem*dem)
      do j=1, isize
         call c_shear_tens(sig_o, s_0, fsp_0(j))
c         fsp_0(j)= fsp_0(j) * (1.d0 + RATE_FACTOR_1)
      end do

      do j=1,isize
        if(sl0(j)*del(j) < 0.0d0 .or. sm0(j)*dem(j) < 0.0d0) then
          C_t(j) = C_d(j)
        else                 
          C_t(j) = E_T
        end if  
      end do

      slf = sl0*C0 + C_t*del
      smf = sm0*C0 + C_t*dem
      stf = sqrt((slf*slf + smf*smf))
      do j = 1, isize
            if (stf(j) /= 0.0d0) then 
               slf(j) = slf(j) / stf(j)
               smf(j) = smf(j) / stf(j)
            end if
      end do

      do j=1, isize
         if (snf(j) < 0.0 ) then 
            if (stf(j) > fsp(j)) then
               stf(j) = fsp(j)
            end if
         else
            if (stf(j) > fsp_0(j)) then
               stf(j) = fsp_0(j)
            end if
         end if
         slf(j) = stf(j) * slf(j)
         smf(j) = stf(j) * smf(j)
      end do

      return 
      end subroutine c_shear2 

C +--------------------------------------------------------------------+
C |                        SUBROUTINE C_VOL                            |
C +--------------------------------------------------------------------+
      subroutine c_vol(dev, ev0, sv0, deps_N, eps_N, R_N, svneg) 
      include 'vaba_param.inc'
      external fb
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv
      common /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     sv0_p, c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21,C_R0,C_R2
      real :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, sv0_p
      real :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      real :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21,C_R0,C_R2
      common /KM7fBounds_M7fIO_1/ c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1,
     $     c_6_2, c_6_M4, c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4,
     $     c_9_0, c_9_1, c_9_2, c_9_M4
      real :: c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_6_M4
      real :: c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4
      real :: c_9_0, c_9_1, c_9_2, c_9_M4
      common /KM7fBounds_M7fIO_2/sig_fib_0, p_1, p_2, p_5,p_6,xp_1,xp_2,
     $     xp_3, V_f, cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3, cf_4,
     $     cf_5, cf_4_4, ap1, bp1, cp2, dp2, ep2, d_1, d_2, d_3, d_4,
     $     d_5, d_6
      real :: sig_fib_0, p_1, p_2, p_5, p_6
      real :: xp_1, xp_2, xp_3, V_f
      real :: cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3
      real :: cf_4, cf_5, cf_4_4
      real :: ap1, bp1, cp2, dp2, ep2
      real :: d_1, d_2, d_3, d_4, d_5, d_6
      save :: /KCommon_Vars/, /KM7fBounds_M7fIO/, /KM7fBounds_M7fIO_1/,
     $     /KM7fBounds_M7fIO_2/
      real :: dev, ev0, sv0
      real, dimension(1:np) :: eps_N, deps_N
      real :: svneg 
      real :: R_N 
      real :: Cv0, evf, svb_neg 
      real :: xk0, e0_V

      real :: prStrainDiff, xk4 
      
      xk0  = k_3*k_1*young 
      e0_V = k_4*k_1 
      Cv0  = young / (1.0d0 - 2.0d0 * poisson) 

      prStrainDiff = maxval(eps_N+deps_N,1) - minval(eps_N+deps_N,1)

      xk4 = (k_6*(prStrainDiff/k_1)**k_7)
     $     /(1.0d0 + min(max(-sv0,0.d0),sv0_p)/Cv0) + k_4

      e0_V = xk4*k_1

      evf = ev0+dev 

      svb_neg = fb(0,evf,xk0,e0_V) 
      svneg=svb_neg*(1.0d0 + C_R2*R_N)

      end subroutine c_vol
 
C +--------------------------------------------------------------------+
C |                            FUNCTION FB                             |    
C +--------------------------------------------------------------------+
      function fb(i,ef,xk0,e0_V) 
      include 'vaba_param.inc'
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv
      common /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     sv0_p, c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21,C_R0,C_R2
      real :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, sv0_p
      real :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      real :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21,C_R0,C_R2
      common /KM7fBounds_M7fIO_1/ c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1,
     $     c_6_2, c_6_M4, c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4,
     $     c_9_0, c_9_1, c_9_2, c_9_M4
      real :: c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_6_M4
      real :: c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4
      real :: c_9_0, c_9_1, c_9_2, c_9_M4
      common /KM7fBounds_M7fIO_2/sig_fib_0, p_1, p_2, p_5,p_6,xp_1,xp_2,
     $     xp_3, V_f, cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3, cf_4,
     $     cf_5, cf_4_4, ap1, bp1, cp2, dp2, ep2, d_1, d_2, d_3, d_4,
     $     d_5, d_6
      real :: sig_fib_0, p_1, p_2, p_5, p_6
      real :: xp_1, xp_2, xp_3, V_f
      real :: cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3
      real :: cf_4, cf_5, cf_4_4
      real :: ap1, bp1, cp2, dp2, ep2
      real :: d_1, d_2, d_3, d_4, d_5, d_6
      save :: /KCommon_Vars/, /KM7fBounds_M7fIO/, /KM7fBounds_M7fIO_1/,
     $     /KM7fBounds_M7fIO_2/

      integer ::  i
      real :: ef,xk0,e0_V
      real :: fb, aux, aux2
      
      aux=-xk0*exp(-ef/(e0_V*1.d0)) 
      if (i.eq.0) fb = aux 

      if (i.eq.1) fb = -aux/e0_V 
      return
      end function fb

C **********************************************************************
C *** SUBROUTINE SETSYSTEM *********************************************
C **********************************************************************
 
      subroutine setsystem() 
      include 'vaba_param.inc'
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv, zeta
      real :: zeta0
      common /KM7f_Microplane_1/ qn, ql, qm, w 
      common /KM7f_Microplane_2/ deps_N, deps_L, deps_M, eps_N,  eps_L,
     $     eps_M,  sig_N
      common /KM7f_Microplane_3/ sd0, sn0, ded, ed0, edf, enf, snf, snb,
     $     sdf
      common /KM7f_Microplane_4/ sl0, sm0, slf, smf, del, dem, vh_ini,
     $     vh_fin, C_d, E_N
      common /KM7f_Microplane_5/ sdneg, sdpos, eps_N0_neg, eps_N0_pos,
     $     snb_fib, Cv, zeta, zeta0
      real, dimension(1:6,1:np) :: qn, ql, qm
      real, dimension(1:np)     :: w
      real, dimension(1:np)     :: deps_N,
     $     deps_L, deps_M, eps_N,  eps_L,  eps_M,  sig_N
      real, dimension(1:np)     :: sd0, sn0,
     $     ded, ed0, edf
      real, dimension(1:np)     :: enf, snf,
     $     snb, sdf, sdneg, sdpos, snb_fib
      real, dimension(1:np)     :: sl0, sm0, slf,
     $     smf, del, dem, eps_N0_neg, eps_N0_pos
      real, dimension(1:nvhi+2)    :: vh_ini, vh_fin 
      real, dimension(1:np)     :: C_d, E_N
      real :: Cv
      save    :: /KCommon_Vars/
      save    :: /KM7f_Microplane_1/,/KM7f_Microplane_2/,
     $     /KM7f_Microplane_3/,/KM7f_Microplane_4/,/KM7f_Microplane_5/
      integer :: jp, ij(1:2,1:6), i, j, k, allocstat
      integer :: rs_size
      real, dimension(1:4,1:np) :: te
      real, dimension(1:3) :: xn, xm, xl, rand_vec 
      real :: lengthn, lengthm, lengthl
C

         vh_ini = 0.0d0  

         vh_fin = 0.0d0  

         qn = 0.0d0

         ql = 0.0d0

         qm = 0.0d0

         w = 0.0d0

         deps_N = 0.0d0

         deps_M = 0.0d0

         deps_L = 0.0d0

         eps_N = 0.0d0

         eps_M = 0.0d0

         eps_L = 0.0d0

         sig_N = 0.0d0

         sd0 = 0.0d0
         sn0 = 0.0d0

         ded = 0.0d0

         ed0 = 0.0d0

         edf = 0.0d0

         enf = 0.0d0

         snf = 0.0d0

         snb = 0.0d0

         snb_fib = 0.0d0

         sdf = 0.0d0

         sdpos = 0.0d0 
         sdneg=0.0d0

         eps_N0_pos = 0.0d0
         zeta = 0.d0
         zeta0 = 0.d0
         eps_N0_neg=0.0d0

         sl0 = 0.0d0
         sm0=0.0d0
         slf=0.0d0
         smf=0.0d0

         del = 0.0d0

         dem = 0.0d0

         C_d=0.0d0
         E_N=0.0d0

      te=0.0d0

      ij=reshape((/1,1,2,2,3,3,1,2,2,3,3,1/),(/2,6/))

      if(np.eq.21) then

      te = reshape(
     $(/0.000000000000D+00,0.000000000000D+00,1.000000000000D+00,
     $        2.652142440930D-02,
     $0.000000000000D+00,1.000000000000D+00,0.000000000000D+00,
     $        2.652142440930D-02,
     $1.000000000000D+00,0.000000000000D+00,0.000000000000D+00,
     $        2.652142440930D-02,
     $0.000000000000D+00,7.071067811870D-01,7.071067811870D-01,
     $        1.993014763120D-02,
     $0.000000000000D+00,-7.071067811870D-01,7.071067811870D-01,
     $        1.993014763120D-02,
     $7.071067811870D-01,0.000000000000D+00,7.071067811870D-01,
     $        1.993014763120D-02,
     $-7.071067811870D-01,0.000000000000D+00,7.071067811870D-01,
     $        1.993014763120D-02,
     $7.071067811870D-01,7.071067811870D-01,0.000000000000D+00,
     $        1.993014763120D-02,
     $-7.071067811870D-01,7.071067811870D-01,0.000000000000D+00,
     $        1.993014763120D-02,
     $8.360955967490D-01,3.879073040670D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $-8.360955967490D-01,3.879073040670D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $8.360955967490D-01,-3.879073040670D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $-8.360955967490D-01,-3.879073040670D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $3.879073040670D-01,8.360955967490D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $-3.879073040670D-01,8.360955967490D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $3.879073040670D-01,-8.360955967490D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $-3.879073040670D-01,-8.360955967490D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $3.879073040670D-01,3.879073040670D-01,8.360955967490D-01,
     $        2.507123674870D-02,
     $-3.879073040670D-01,3.879073040670D-01,8.360955967490D-01,
     $        2.507123674870D-02,
     $3.879073040670D-01,-3.879073040670D-01,8.360955967490D-01,
     $        2.507123674870D-02,
     $-3.879073040670D-01,-3.879073040670D-01,8.360955967490D-01,
     $        2.507123674870D-02/)
     $,(/4,np/))
      else if(np.eq.28) then

      te = reshape(
     $(/.577350259D0, .577350259D0, .577350259D0, .0160714276D0,
     $.577350259D0, .577350259D0,-.577350259D0, .0160714276D0,
     $.577350259D0,-.577350259D0, .577350259D0, .0160714276D0,
     $.577350259D0,-.577350259D0,-.577350259D0, .0160714276D0,
     $.935113132D0, .250562787D0, .250562787D0, .0204744730D0,
     $.935113132D0, .250562787D0,-.250562787D0, .0204744730D0,
     $.935113132D0,-.250562787D0, .250562787D0, .0204744730D0,
     $.935113132D0,-.250562787D0,-.250562787D0, .0204744730D0,
     $.250562787D0, .935113132D0, .250562787D0, .0204744730D0,
     $.250562787D0, .935113132D0,-.250562787D0, .0204744730D0,
     $.250562787D0,-.935113132D0, .250562787D0, .0204744730D0,
     $.250562787D0,-.935113132D0,-.250562787D0, .0204744730D0,
     $.250562787D0, .250562787D0, .935113132D0, .0204744730D0,
     $.250562787D0, .250562787D0,-.935113132D0, .0204744730D0,
     $.250562787D0,-.250562787D0, .935113132D0, .0204744730D0,
     $.250562787D0,-.250562787D0,-.935113132D0, .0204744730D0,
     $.186156720D0, .694746614D0, .694746614D0, .0158350505D0,
     $.186156720D0, .694746614D0,-.694746614D0, .0158350505D0,
     $.186156720D0,-.694746614D0, .694746614D0, .0158350505D0,
     $.186156720D0,-.694746614D0,-.694746614D0, .0158350505D0,
     $.694746614D0, .186156720D0, .694746614D0, .0158350505D0,
     $.694746614D0, .186156720D0,-.694746614D0, .0158350505D0,
     $.694746614D0,-.186156720D0, .694746614D0, .0158350505D0,
     $.694746614D0,-.186156720D0,-.694746614D0, .0158350505D0,
     $.694746614D0, .694746614D0, .186156720D0, .0158350505D0,
     $.694746614D0, .694746614D0,-.186156720D0, .0158350505D0,
     $.694746614D0,-.694746614D0, .186156720D0, .0158350505D0,
     $.694746614D0,-.694746614D0,-.186156720D0, .0158350505D0/)
     $,(/4,np/))     
      else if(np.eq.37) then

      te = reshape(
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
     $,(/4,np/))

      end if   

C     ---------------------------------------
C ... Assemble tensors from direction cosines
C     ---------------------------------------
      qn = 0.0d0
      ql = 0.0d0
      qm = 0.0d0
       w = 0.0d0
      call random_seed(size=rs_size)
      do jp=1,np 
        w(jp) = te(4,jp)*6.0D0 
        xn(1) = te(3,jp) 
        xn(2) = te(2,jp) 
        xn(3) = te(1,jp) 

        lengthm = 0.0d0
        do while (lengthm .lt. epsilon(lengthm))
           call random_number(rand_vec)
           xm = rand_vec - dot_product(xn,rand_vec)*xn
           lengthm = sqrt(dot_product(xm,xm))
        end do
        xm = xm/lengthm

        xl(1) = xn(2)*xm(3)-xn(3)*xm(2) 
        xl(2) = xn(3)*xm(1)-xn(1)*xm(3) 
        xl(3) = xn(1)*xm(2)-xn(2)*xm(1) 
        lengthl = sqrt(dot_product(xl,xl))
        xl=xl/lengthl
        lengthn = sqrt(dot_product(xn,xn))
        lengthm = sqrt(dot_product(xm,xm))
        lengthl = sqrt(dot_product(xl,xl))
        do k=1,6 
          i=ij(1,k) 
          j=ij(2,k) 
          qn(k,jp) = xn(i)*xn(j) 
          qm(k,jp) = 0.5D0*(xn(i)*xm(j)+xn(j)*xm(i)) 
          ql(k,jp) = 0.5D0*(xn(i)*xl(j)+xn(j)*xl(i)) 
        end do 
      end do 
      return
      end subroutine setsystem 

C *******************************************************************
C *** SUBROUTINE M7fMATERIAL *****************************************
C *******************************************************************

      subroutine M7fMaterial(
     $     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     $     stepTime, totalTime, dt, cmname, coordMp, charLength,
     $     props, density, strainInc, relSpinInc,
     $     tempOld, stretchOld, defgradOld, fieldOld,
     $     stressOld, stateOld, enerInternOld, enerInelasOld,
     $     tempNew, stretchNew, defgradNew, fieldNew,
     $     stressNew, stateNew, enerInternNew, enerInelasNew )

      include 'vaba_param.inc'

      external RateFunc
C
      dimension props(nprops), density(nblock),
     $coordMp(nblock,*),
     $charLength(nblock), strainInc(nblock,ndir+nshr),
     $relSpinInc(nblock,nshr), tempOld(nblock),
     $stretchOld(nblock,ndir+nshr), defgradOld(nblock,ndir+nshr+nshr),
     $fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     $stateOld(nblock,nstatev), enerInternOld(nblock),
     $enerInelasOld(nblock), tempNew(nblock),
     $stretchNew(nblock,ndir+nshr), defgradNew(nblock,ndir+nshr+nshr), 
     $fieldNew(nblock,nfieldv),
     $stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     $enerInternNew(nblock), enerInelasNew(nblock)

      character*80 cmname

C
      real :: equivStress, fractureWorkInc, smean, stressPower, two

      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv, zeta
      real :: sum_bulk, bulk, zeta0
      common /KIFED_M7f/ StrainRateTensorInv2, DevStrainRateTensorInv2,
     $     old_dtime, lruncreep, lfirstpass 
      real :: StrainRateTensorInv2,
     $     DevStrainRateTensorInv2, old_dtime
      common /KM7f_M7fIO/ on_rate_dep_frac, on_creep
      common /KM7f_M7fIO_1/ C_R1, q_1
      common /KM7f_M7fIO_2/ t_0, t_ch, f_c, a_to_c, g_c, w_to_c, E28
      integer      :: on_rate_dep_frac, on_creep
      real :: C_R1, q_1  
      real :: t_0, t_ch, f_c, a_to_c, g_c, w_to_c, E28
      common /KM7f_Microplane_1/ qn, ql, qm, w 
      common /KM7f_Microplane_2/ deps_N, deps_L, deps_M, eps_N,  eps_L,
     $     eps_M,  sig_N
      common /KM7f_Microplane_3/ sd0, sn0, ded, ed0, edf, enf, snf, snb,
     $     sdf
      common /KM7f_Microplane_4/ sl0, sm0, slf, smf, del, dem, vh_ini,
     $     vh_fin, C_d, E_N
      common /KM7f_Microplane_5/ sdneg, sdpos, eps_N0_neg, eps_N0_pos,
     $     snb_fib, Cv, zeta, zeta0
      real, dimension(1:6,1:np) :: qn, ql, qm
      real, dimension(1:np)     :: w
      real, dimension(1:np)     :: deps_N,
     $     deps_L, deps_M, eps_N,  eps_L,  eps_M,  sig_N
      real, dimension(1:np)     :: sd0, sn0,
     $     ded, ed0, edf
      real, dimension(1:np)     :: enf, snf,
     $     snb, sdf, sdneg, sdpos, snb_fib
      real, dimension(1:np)     :: sl0, sm0, slf,
     $     smf, del, dem, eps_N0_neg, eps_N0_pos

      real, dimension(1:nvhi+2)    :: vh_ini, vh_fin 
      real, dimension(1:np)     :: C_d, E_N
      real :: Cv
      save    :: /KCommon_Vars/, /KIFED_M7f/
      save    :: /KM7f_M7fIO/, /KM7f_M7fIO_1/, /KM7f_M7fIO_2/
      save    :: /KM7f_Microplane_1/, /KM7f_Microplane_2/,
     $     /KM7f_Microplane_3/, /KM7f_Microplane_4/, /KM7f_Microplane_5/

      real, dimension(1:6) :: eps, deps, sig_old 
      real, dimension(1:6) :: sig
      real, dimension(1:6) :: deps_s
      real :: R_N, dfb
      real :: sv0, ev0, dev, evf, sum_snf, devv
      real :: svfB, svf, phi0, phi, sig_I_old, eps_I
      integer ::  jp, i 
      real :: svneg, sum1, strain_rate, ek(1:6,1:6)

      NTENS = 6

      if (totalTime<=dt) then
         call inputparams()
         call setsystem() 
      end if

      do i = 1, nblock

         vh_ini = stateOld(i,:) 

         if (totalTime<=dt) then
            vh_ini(2)=1.0d0
         end if

         DO IC = 1, NTENS
            sig(IC)=stressOld(i,IC)
            sig_old(IC)=stressOld(i,IC)
            sig(IC)=sig(IC)*unit_conv
            sig_old(IC)=sig_old(IC)*unit_conv
            IF (IC > 3) THEN    
               deps(IC)=2.0*strainInc(i,IC)
            ELSE
               deps(IC)=strainInc(i,IC)
            END IF
         END DO

         eps(1)=(defgradOld(i,1)*defgradOld(i,1) + 
     $           defgradOld(i,7)*defgradOld(i,7) +
     $           defgradOld(i,6)*defgradOld(i,6) - 1.0)/2.0
         eps(2)=(defgradOld(i,4)*defgradOld(i,4) + 
     $           defgradOld(i,2)*defgradOld(i,2) +
     $           defgradOld(i,8)*defgradOld(i,8) - 1.0)/2.0
         eps(3)=(defgradOld(i,9)*defgradOld(i,9) + 
     $           defgradOld(i,5)*defgradOld(i,5) +
     $           defgradOld(i,3)*defgradOld(i,3) - 1.0)/2.0
         eps(4)=(defgradOld(i,1)*defgradOld(i,4) +
     $           defgradOld(i,2)*defgradOld(i,7) +
     $           defgradOld(i,6)*defgradOld(i,8))/2.0
         eps(5)=(defgradOld(i,9)*defgradOld(i,4) + 
     $           defgradOld(i,2)*defgradOld(i,5) +
     $           defgradOld(i,3)*defgradOld(i,8))/2.0
         eps(6)=(defgradOld(i,1)*defgradOld(i,9) + 
     $           defgradOld(i,7)*defgradOld(i,5) +
     $           defgradOld(i,3)*defgradOld(i,6))/2.0
         DO IC = 1, NTENS
            IF (IC > 3) THEN    ! SHEAR STRAINS ARE SHEAR ANGLES IN MICROPLANE MODEL
               eps(IC)=2.0*eps(IC)
            END IF
         END DO

         deps_s(1)=deps(1)*deps(1)+deps(4)*deps(4)+deps(6)*deps(6)
         deps_s(2)=deps(4)*deps(4)+deps(2)*deps(2)+deps(5)*deps(5)
         deps_s(3)=deps(6)*deps(6)+deps(5)*deps(5)+deps(3)*deps(3)
         deps_s(4)=deps(1)*deps(4)+deps(2)*deps(4)+deps(5)*deps(6)
         deps_s(5)=deps(4)*deps(6)+deps(2)*deps(5)+deps(5)*deps(3)
         deps_s(6)=deps(1)*deps(6)+deps(4)*deps(5)+deps(6)*deps(3)
         sum1=0.d0
         do jp = 1, np
            sum1=sum1+(qn(1,jp)*deps_s(1)+qn(4,jp)*deps_s(4)+
     $           qn(6,jp)*deps_s(6)+qn(4,jp)*deps_s(4)+
     $           qn(2,jp)*deps_s(2)+qn(5,jp)*deps_s(5)+
     $           qn(6,jp)*deps_s(6)+qn(5,jp)*deps_s(5)+
     $           qn(3,jp)*deps_s(3))*w(jp)
         end do
         strain_rate=sqrt(sum1/2.d0)/dt
 
C     -------------------------------------------------------------------
C     Compute material response :
C     Bounding curve formulation for volumetric,
C ... deviatoric and shear variables.
C     Shear strains are shear angles.
C     -------------------------------------------------------------------

C     ---------------
C ... Initializations
C     ---------------

      jp = size(w) 
      
      sig = 0.0d0
       
C     -------------------------------------------------------------------
C ... Calculate the rate coefficients using the rate function
C     -------------------------------------------------------------------
      
      if(on_rate_dep_frac.eq.0) then 
        R_N=0.0d0
      else  

        R_N = RateFunc(strain_rate/C_R1) 

      end if  

C     --------------------------------------------
C ... Evaluate local strains (eps_N, eps_L, eps_M)
C     --------------------------------------------
  
       eps_N =  matmul( eps,qn) 
       eps_L =  matmul( eps,ql)
       eps_M =  matmul( eps,qm)
      deps_N =  matmul(deps,qn)
      deps_L =  matmul(deps,ql)
      deps_M =  matmul(deps,qm)

      if (maxval(eps_N)>th_del) then
         vh_fin(188)=0.
      else
         vh_fin(188)=1.
      end if

C     -----------------------------------------
C ... Volumetric law (same for all microplanes)
C     -----------------------------------------
 
      sv0 = vh_ini(1)
      phi0 = vh_ini(2)
      
      ev0 =  (eps(1)+ eps(2)+ eps(3))/3.0d0
      dev = (deps(1)+deps(2)+deps(3))/3.0d0
      evf = ev0+dev 
      zeta0 = vh_ini(189)
      zeta = zeta0

      call c_vol(dev, ev0, sv0, deps_N, eps_N, R_N, svneg)

C     ---------------------
C ... Normal deviatoric law
C     ---------------------

      sn0 = vh_ini(3:nvhi:nvhm)

      ded = deps_N - dev
      ed0 = eps_N - ev0     
      edf = ed0 + ded

      call c_dev(ded, ed0, dev, ev0, sv0, sd0, sdf, C_d, R_N,
     $     sdneg, sdpos)

      enf = evf + edf 

      eps_N0_pos=vh_ini(6:nvhi:nvhm)
      eps_N0_neg=vh_ini(7:nvhi:nvhm)
      call c_norm_elastic(eps_N, deps_N, sn0, eps_N0_pos,
     $     eps_N0_neg, sv0, snf, E_N, zeta)

      call c_norm(enf, sv0, R_N, snb)
      call c_norm_fib(enf, R_N, snb_fib)

      phi0 = vh_ini(2)

      sig_N = max(min(snf,snb+snb_fib),svneg+sdneg) 

      do ic = 1, jp
         if (snf(ic) > snb(ic)+snb_fib(ic)) then 
               eps_N0_pos(ic) = enf(ic)
         else if (snf(ic) < svneg+sdneg(ic)) then 
            eps_N0_neg(ic) = enf(ic)
         end if
      end do
      
      sum_snf = dot_product(sig_N,w)              
      svfB = sum_snf/3.0d0 

      svf = svfB

      sum_bulk = dot_product(E_N,w)   
      bulk = sum_bulk/3.0d0

      if(sv0 > 0.d0) then
         if(svf > 0.d0) then
            devv = abs(dev - ((svf-sv0)/bulk))
            zeta = zeta0 + abs(devv)

         end if
      end if
C     -----------------
C ... History variables
C     -----------------

      sl0 = vh_ini(4:nvhi:nvhm) 
      sm0 = vh_ini(5:nvhi:nvhm) 
 
      snf = sig_N

      enf = eps_N + deps_N
      del = deps_L
      dem = deps_M

C     ---------
C ... Shear law
C     ---------

      call c_shear2(eps_L, eps_M, C_d, snf, R_N, del, dem, sl0, sm0,
     $     slf, smf, evf)

C     -------------------------------------------------
C     Final Stress Vector
C     -------------------------------------------------
      
      sig = matmul(qn,snf*w) + matmul(qm,smf*w) + matmul(ql,slf*w)

      sig = sig/unit_conv
C
C     -------------------------------------------
C ... Update microplane normal and shear stresses
C     -------------------------------------------
      
      vh_fin(1) = svf

      vh_fin(3:nvhi:nvhm) = snf 
      vh_fin(4:nvhi:nvhm) = slf 
      vh_fin(5:nvhi:nvhm) = smf 
      vh_fin(6:nvhi:nvhm) = eps_N0_pos
      vh_fin(7:nvhi:nvhm) = eps_N0_neg 
      vh_fin(189)=zeta

      stateNew(i,:)=vh_fin
      stressNew(i,:)=sig

      stressPower = (( stressOld(i,1)+stressNew(i,1) )*
     $        strainInc(i,1) + ( stressOld(i,2)+stressNew(i,2) )*
     $        strainInc(i,2) + ( stressOld(i,3)+stressNew(i,3) )*
     $        strainInc(i,3) + 2.0*( stressOld(i,4)+stressNew(i,4) )*
     $        strainInc(i,4) )/2.0

      enerInternNew(i) = enerInternOld(i) + stressPower / density(i)

      smean = ( stressNew(i,1) + stressNew(i,2) +
     $        stressNew(i,3) )/3.0
      equivStress = sqrt( 3.0/2.0 * ( (stressNew(i,1)-smean)**2 +
     $        (stressNew(i,2)-smean)**2 + (stressNew(i,3)-smean)**2 +
     $        two * stressNew(i,4)**2 ) )

      fractureWorkInc = stressNew(i,1)*deps(1) +
     $        stressNew(i,2)*deps(2) + stressNew(i,3)*deps(3) +
     $        stressNew(i,4)*deps(4) + stressNew(i,5)*deps(5) +
     $        stressNew(i,6)*deps(6) 
      enerInelasNew(i) = enerInelasOld(i)+fractureWorkInc/density(i)

      end do

      return
      end subroutine M7fMaterial

C***********************************************************************
C**** This function characterizes the theoretical rate effect      *****
C***********************************************************************
      function RateFunc(arg)
      include 'vaba_param.inc'
      real :: RateFunc
      real :: arg

      RateFunc = log(arg + sqrt(arg * arg + 1.0d0))
      return
      end function RateFunc

C **********************************************************************
C *** SUBROUTINE CREEP *************************************************
C **********************************************************************
      subroutine M7f_Creep(dt)
C+---------------------------------------------------------------------+
C|  Replaces the Elastic Modulus E by the new one Ep which accounts for|
C|  creep effects and stores the coefficient C0                        |
C+---------------------------------------------------------------------+
      include 'vaba_param.inc'
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv
      common /KM7f_M7fIO/ on_rate_dep_frac, on_creep
      common /KM7f_M7fIO_1/ C_R1, q_1
      common /KM7f_M7fIO_2/ t_0, t_ch, f_c, a_to_c, g_c, w_to_c, E28
      integer      :: on_rate_dep_frac, on_creep
      real :: C_R1, q_1  
      real :: t_0, t_ch, f_c, a_to_c, g_c, w_to_c, E28
      real :: dt
      real, parameter  :: DaysPerSec = 1.1574074074074074d-5 
      real :: t_age, t_ch_p, density, q1_opt
      real :: q1, q2, q3, q4, n, m, t, J, r, Z, Qf, Q, J_dot
      real :: Ep, eta, tau1, C1

      t_age=t_0*DaysPerSec
      t_ch_p =t_ch*DaysPerSec/2.0d0 
      density=g_c
      q1_opt=q_1

      if(f_c.eq.0.0d0) then
        f_c=E28/4734.0d0*E28/4734.0d0 
      else
        E28=4734.0d0*sqrt(f_c)
      end if

      q2  = 185.4d0*sqrt(density)*f_c**(-0.9d0)
      q4  = 20.3d0*a_to_c**(-0.7d0)
      if(w_to_c.eq.0.0d0) w_to_c=1.0d0/(f_c/22.8d0+0.535d0)

      q1 = 0.60d+6/E28
      q3 = 0.29d0*w_to_c**4*q2
      n=0.1d0
      m=0.5d0
      t=t_age+t_ch_p

      r=1.7d0*t_age**0.12d0+8.0d0
      Z=t_age**(-m)*log(1.0d0+(t-t_age)**n)
      Qf=1.0d0/(0.086d0*t_age**(2.d0/9.d0)+1.21d0*t_age**(4.0d0/9.0d0))
      Q=Qf*Z/(Z**r+Qf**r)**(1.0d0/r)
      J=q1+q2*Q+q3*log(1.0d0+(t-t_age)**n)+q4*log(t/t_age)
      J=q1_opt*(J/q1)

      J_dot=n*(q2*t**(-m)+q3)/(t-t_age+(t-t_age)**(1.0d0-n))+q4/t
      J_dot=q1_opt*(J_dot/q1)
      Ep=1.0d+6/(J-t_ch_p*J_dot)

      eta=1.0d+6/J_dot
      tau1=eta/Ep

      if(dt/tau1.ge.1.0d+2) then
        C0=0.0d0
      else
        C0=exp(-dt/tau1*DaysPerSec)
      end if
      if(dt/tau1.le.1.0d-3) then
        C1=1.0d0 - dt/tau1/2.0d0*DaysPerSec
      else
        C1=(1.0d0-C0)/(dt/tau1*DaysPerSec)
      end if
      Ep=C1*Ep
      young = Ep

      return
      end subroutine M7f_Creep

C **********************************************************************
C *** SUBROUTINE INPUTPARAMS *******************************************
C **********************************************************************
      subroutine inputparams() 
      include 'vaba_param.inc'
      real, parameter :: PI=3.1415926535897932384626433832795d0
      integer, parameter :: np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      common /KCommon_Vars/young, poisson,k_1, C0,dtime,th_del,unit_conv  
      real :: young, poisson, k_1, C0, dtime, th_del, unit_conv
      common /KM7f_M7fIO/ on_rate_dep_frac, on_creep
      common /KM7f_M7fIO_1/ C_R1, q_1
      common /KM7f_M7fIO_2/ t_0, t_ch, f_c, a_to_c, g_c, w_to_c, E28
      integer      :: on_rate_dep_frac, on_creep
      real :: C_R1, q_1  
      real :: t_0, t_ch, f_c, a_to_c, g_c, w_to_c, E28
      common /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     sv0_p, c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21,C_R0,C_R2
      real :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, sv0_p
      real :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      real :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21,C_R0,C_R2
      common /KM7fBounds_M7fIO_1/ c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1,
     $     c_6_2, c_6_M4, c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4,
     $     c_9_0, c_9_1, c_9_2, c_9_M4
      real :: c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_6_M4
      real :: c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4
      real :: c_9_0, c_9_1, c_9_2, c_9_M4
      common /KM7fBounds_M7fIO_2/sig_fib_0, p_1, p_2, p_5,p_6,xp_1,xp_2,
     $     xp_3, V_f, cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3, cf_4,
     $     cf_5, cf_4_4, ap1, bp1, cp2, dp2, ep2, d_1, d_2, d_3, d_4,
     $     d_5, d_6
      real :: sig_fib_0, p_1, p_2, p_5, p_6
      real :: xp_1, xp_2, xp_3, V_f
      real :: cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3
      real :: cf_4, cf_5, cf_4_4
      real :: ap1, bp1, cp2, dp2, ep2
      real :: d_1, d_2, d_3, d_4, d_5, d_6
      save :: /KCommon_Vars/,/KM7f_M7fIO/,/KM7f_M7fIO_1/,/KM7f_M7fIO_2/,
     $     /KM7fBounds_M7fIO/,/KM7fBounds_M7fIO_1/,/KM7fBounds_M7fIO_2/

      th_del = 0.005d10 

      unit_conv=1.                

      young  = 24541.0d0
      poisson= 0.18000d0
      k_1    = 127.80d-6 
      k_2    = 40.000d0
      k_3    = 20.d0 
      k_4    = 36.d0 
      c_5_0 = 1.3d-2
      c_5_1 = 4.0d0 
      c_7_0  = 1.2d-2 
      c_7_1 = 35.0d2 
      c_8_0  = 1.2d-2 
      c_8_1 = 20.0d0 
      c_6_0 = 4.0d2 
      c_6_1  = 4.0d1 
      c_6_2  = 13.0d0
      c_9_0 = 4.0d2 
      c_9_1  = 4.0d1 
      c_9_2  = 13.0d0
      c_5_M4 = 3.0d0    
      c_6_M4 = 1.30d0 
      c_7_M4 = 10.d-1   
      c_8_M4 = 8.0d0  
      c_9_M4 = 0.d-1 

      d_1=0.095d0 
      d_2=35.d0   
      d_3=1.7d0   
      d_4=1.7d0 
      d_5 = 1.d3
      d_6 = 25.d0
     
      sv0_p=250.d0
      k_8=3.d0
      k_9=5.d-1
      k_5=1.0d0
      k_6=1.0d-4 
      k_7=1.8d0 
     
      c_2 = 2.1d-0 
      c_3 = 4.d0 
      c_4 = 60.d0 
      c_10 = 3.3d-1 
      c_11 = 5.d-1 
      c_12 = 7.00d3 
      c_16 = 10.0d0 
      c_17 = 1.00d-2       
      c_18 = 4.d3 
      c_19 = 4.5d3 
      c_20 = 3.d2 
      c_21 = 6.d1  
     
      sig_fib_0 = -4.5d0 
      p_1    = 0.0d0     
      p_2    = 300.d-3   
      p_5    = 0.0d0     
      p_6    = 0.0d0     
      xp_1   = 3.34d0    
      xp_2   = 0.        
      xp_3   = 1.0d0     
      V_f    = 0.0d0     
      cf_1=104.0d-2
      cf_2=54.0d+1
      cf_3=104.0d-2
      cf_1_1=3.5d+1
      cf_2_2=1.5d+1
      cf_3_3=3.5d+1
      cf_4_4=4.4d2 
      cf_4  = 3.2d-2 
      cf_5  = 5.0d-2 
     
      on_rate_dep_frac = 0
      C_R1 = 4.0d-6
      C_R2 = 22.0d-3
     
      on_creep = 0
      t_0  = 241.920d+4
      t_ch = 180.000d+0
      f_c  =  48.000d+0
      a_to_c = 7.000d+0
      g_c  = 304.340d+0
      w_to_c = 0.000d+0
     
      q_1=21.00d0        
      C0=1.0d0           

      E28 = young
      if(t_ch.eq.0.0d0) then
         write(*,*) 'Error: Duration of loading is not allowed to'//
     $        ' be zero.'
         write(*,*) 'If you want to turn off creep effects, use the '//
     $        'provided switch.'
         stop
      end if
      return
      end subroutine inputparams
