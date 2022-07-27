c...  ------------------------------------------------------------------
      subroutine sdvini(statev,coords,nstatv,ncrds,noel,npt,layer,kspt)
c...  ------------------------------------------------------------------
      include 'aba_param.inc'


      dimension statev(nstatv)

      statev(1)=1.0d0
      statev(2)=1.0d0
      statev(3)=1.0d0
      statev(4)=1.0d0
      statev(5)=1.0d0
      statev(6)=1.0d0
      statev(7)=1.0d0
      statev(8)=1.0d0
      statev(9)=1.0d0
      statev(10)=1.0d0
      statev(11)=0.0d0

      return
      end

c...  ------------------------------------------------------------------
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     #rpl,ddsddt,drplde,drpldt,
     #stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     #ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     #celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
c...  ------------------------------------------------------------------
      include 'aba_param.inc'

      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     #ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     #stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     #props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      call umat_H(stress,statev,ddsdde,sse,
     #                       time,dtime,coords,props,dfgrd1,
     #                       ntens,ndi,nshr,nstatv,nprops,
     #                       noel,npt,kstep,kinc)

      return
      end
c...  
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     #                             NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     #                             JMAC,JMATYP,MATLAYO,LACCFLA) 
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C Stress tensor:
      CALL GETVRM('SDV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     #                              MATLAYO,LACCFLA)
      UVAR(1) = ARRAY(1)
      UVAR(2) = ARRAY(2)
      UVAR(3) = ARRAY(3)
      UVAR(4) = ARRAY(4)
      UVAR(5) = ARRAY(5)
      UVAR(6) = ARRAY(6)
      UVAR(7) = ARRAY(7)
      UVAR(8) = ARRAY(8)
      UVAR(9) = ARRAY(9)
      UVAR(10) = ARRAY(10)
      UVAR(11) = ARRAY(11)

      RETURN
      END

c...  ------------------------------------------------------------------
      subroutine umat_H(stress,statev,ddsdde,sse,
     #                             time,dtime,coords,props,dfgrd1,
     #                             ntens,ndi,nshr,nstatv,nprops,
     #                             noel,npt,kstep,kinc)


c...  ------------------------------------------------------------------
 
c...  ------------------------------------------------------------------

      implicit none

c...  variables to be defined
      real*8  stress(ntens), ddsdde(ntens,ntens), statev(nstatv), sse

c...  variables passed in for information
      real*8  time(2), dtime, coords(3), props(nprops), dfgrd1(3,3)
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt, kstep, kinc

c...  local variables (mostly mechanics part)
      real*8 a0(3), da(6), da2D(3)
      real*8 detf, detF2d, finv(3,3)
      real*8 c(6), cmat(3,3), I1b, bvoigt(6), bmat(3,3), biso(6)
      real*8 Fa0(3), I1bar
      real*8 kronmat(3,3), kronmat2D(2,2), kron(6), kron2D(3)
      real*8 p, sigma2D(3),  sigma(6), sigmamat(3,3), sigmamat2D(2,2)
      real*8 sigmaiso(6), sigmaiso2D(3), sigmabar(6), dpdJ, Psivol, tr_sigmabar 
      real*8 PKfiber, E_damaged, E_undamaged, Edam_Eundam, sigmaf(6) 
      real*8 CCfiber1, CCfiber2, ccfiber
      real*8 IIII, III
      real*8 cc(4,4), ccf(4,4), cciso(4,4), ccvol(4,4) 
      

c...  some auxiliar variables, tensors
      integer i, j, l, m, n, q, r, s, t, v, nitl, II, JJ, Itoi(6), Itoj(6)
      integer Itoi2D(3), Itoj2D(3)

c...  material properties
      real*8 mu0, mu1, kappa, mutheta, I0kappa, gama, beta, c1, c2, lam_m, K_vol

c...  variables and constants for Weibull Distribution
      real*8 pi, theta, delta, lam, VM, delta0, dtheta
      real*8 a, b, p_lams_a, p_lams_b, f_a, f_b
      real*8 f_Edam_a, f_Edam_b, tot_sum_Edam, f_Edam, integral_Edam, E_fiber
      real*8 p0_lams_a, p0_lams_b, f_Eundam_a, f_Eundam_b, tot_sum_Eundam, f_Eundam, integral_Eundam, Eun_fiber
      real*8 p_lams, lam_bar, f, h, tot_sum, integral
      real*8 f_a1, f_b1, f_a2, f_b2, f1, f2, tot_sum1, tot_sum2, integral1, integral2
      real*8 lam_m_n(10)

c      print all props
c      print *, 'all props'
c      do i=1,10
c        print *, props(i)
c      end do
c...  initialize material parameters
      mu0   = props(1)
      mu1 = props(2)
      kappa = props(3)
      mutheta = props(4)
      I0kappa = props(5)
      gama = props(6)
      beta = props(7)
      c1 = props(8)
      c2 = props(9)
      K_vol = props(10)
c      delta = 0.5

      lam_m_n(1) = statev(1)
      lam_m_n(2) = statev(2)
      lam_m_n(3) = statev(3)
      lam_m_n(4) = statev(4)
      lam_m_n(5) = statev(5)
      lam_m_n(6) = statev(6)
      lam_m_n(7) = statev(7)
      lam_m_n(8) = statev(8)
      lam_m_n(9) = statev(9)
      lam_m_n(10) = statev(10)
      Edam_Eundam = statev(11)

c      print *, lam_m_n(10)
c      print *, dfgrd1(1,1),dfgrd1(1,2),dfgrd1(1,3)
c      print *, dfgrd1(2,1),dfgrd1(2,2),dfgrd1(2,3)
c      print *, dfgrd1(3,1),dfgrd1(3,2),dfgrd1(3,3)

c..................................................................
c... For plane strain f13=f23=f31=f32=0 and f33=1
      detF2d = +dfgrd1(1,1)*dfgrd1(2,2) - dfgrd1(1,2)*dfgrd1(2,1)
      dfgrd1(1,3) = 0
      dfgrd1(2,3) = 0
      dfgrd1(3,1) = 0
      dfgrd1(3,2) = 0
      dfgrd1(3,3) = 1
c...  fnew = [[f11, f12,0],[f21, f22, 0],[0,0,1]]
c..................................................................
c...  calculate determinant of deformation gradient
      detf = +dfgrd1(1,1)*(dfgrd1(2,2)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,2))
     #       -dfgrd1(1,2)*(dfgrd1(2,1)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,1))
     #       +dfgrd1(1,3)*(dfgrd1(2,1)*dfgrd1(3,2)-dfgrd1(2,2)*dfgrd1(3,1))
c      print *, 'detF'
c      print *, detf


c...  calculate inverse of F
      finv(1,1) = (+dfgrd1(2,2)*dfgrd1(3,3) - dfgrd1(2,3)*dfgrd1(3,2))/detf
      finv(1,2) = (-dfgrd1(1,2)*dfgrd1(3,3) + dfgrd1(1,3)*dfgrd1(3,2))/detf
      finv(1,3) = (+dfgrd1(1,2)*dfgrd1(2,3) - dfgrd1(1,3)*dfgrd1(2,2))/detf
      finv(2,1) = (-dfgrd1(2,1)*dfgrd1(3,3) + dfgrd1(2,3)*dfgrd1(3,1))/detf
      finv(2,2) = (+dfgrd1(1,1)*dfgrd1(3,3) - dfgrd1(1,3)*dfgrd1(3,1))/detf
      finv(2,3) = (-dfgrd1(1,1)*dfgrd1(2,3) + dfgrd1(1,3)*dfgrd1(2,1))/detf
      finv(3,1) = (+dfgrd1(2,1)*dfgrd1(3,2) - dfgrd1(2,2)*dfgrd1(3,1))/detf
      finv(3,2) = (-dfgrd1(1,1)*dfgrd1(3,2) + dfgrd1(1,2)*dfgrd1(3,1))/detf
      finv(3,3) = (+dfgrd1(1,1)*dfgrd1(2,2) - dfgrd1(1,2)*dfgrd1(2,1))/detf

c...      C = F^T*F
c...      b = F*F^T -> full notation [[b11, b12, b13],[b12,b22,b23],[b13,b23,b33]]
c...      b = [b11,b22,b33,b12,b13,b23] -> voigt notation

c...  calculate right cauchy-green deformation tensor c = f^t*f
      c(1) = dfgrd1(1,1)*dfgrd1(1,1) + dfgrd1(2,1)*dfgrd1(2,1) + dfgrd1(3,1)*dfgrd1(3,1)
      c(2) = dfgrd1(1,2)*dfgrd1(1,2) + dfgrd1(2,2)*dfgrd1(2,2) + dfgrd1(3,2)*dfgrd1(3,2)
      c(3) = dfgrd1(1,3)*dfgrd1(1,3) + dfgrd1(2,3)*dfgrd1(2,3) + dfgrd1(3,3)*dfgrd1(3,3)
      c(4) = dfgrd1(1,1)*dfgrd1(1,2) + dfgrd1(2,1)*dfgrd1(2,2) + dfgrd1(3,1)*dfgrd1(3,2)
      c(5) = dfgrd1(1,1)*dfgrd1(1,3) + dfgrd1(2,1)*dfgrd1(2,3) + dfgrd1(3,1)*dfgrd1(3,3)
      c(6) = dfgrd1(1,2)*dfgrd1(1,3) + dfgrd1(2,2)*dfgrd1(2,3) + dfgrd1(3,2)*dfgrd1(3,3)
c...      print *, 'c'
c...      print *, c(1), c(2), c(3), c(4), c(5), c(6)

c..............................................................
c...  calculate new right cauchy-green deformation tensor c = fnew^T * fnew
c...  commenting because it should be taken care of based on new F
c      c(1) = dfgrd1(1,1)*dfgrd1(1,1) + dfgrd1(2,1)*dfgrd1(2,1)
c      c(2) = dfgrd1(1,2)*dfgrd1(1,2) + dfgrd1(2,2)*dfgrd1(2,2)
c      c(3) = 1
c      c(4) = dfgrd1(1,1)*dfgrd1(1,2) + dfgrd1(2,1)*dfgrd1(2,2)
c      c(5) = 0
c      c(6) = 0
c.............................................................. 

c... c full notation
      cmat(1,1) = c(1)
      cmat(2,2) = c(2)
      cmat(3,3) = c(3)
      cmat(1,2) = c(4)
      cmat(2,1) = c(4)
      cmat(1,3) = c(5)
      cmat(3,1) = c(5)
      cmat(2,3) = c(6)
      cmat(3,2) = c(6)

c... Setting the PKfiber with zeros
      sigmaf(1) = 0.0
      sigmaf(2) = 0.0
      sigmaf(3) = 0.0
      sigmaf(4) = 0.0
      sigmaf(5) = 0.0
      sigmaf(6) = 0.0

      E_damaged = 0.0                                                  ! *Added for Damaged Energy Calculation
      E_undamaged = 0.0                                                ! *Added for Undamaged Energy Calculation

c... Fiber part
      pi = 3.14159265
      do i=1,10
        theta = (i*pi/10.0)
        dtheta = 0.3141592                                 ! Added for integral Von Mises 1/14/21
        a0(1) = cos(theta)
        a0(2) = sin(theta)
        a0(3) = 0.0
c...    Deformed fiber a = F*a0
        Fa0(1) = dfgrd1(1,1)*a0(1) + dfgrd1(1,2)*a0(2) + dfgrd1(1,3)*a0(3)
        Fa0(2) = dfgrd1(2,1)*a0(1) + dfgrd1(2,2)*a0(2) + dfgrd1(2,3)*a0(3)
        Fa0(3) = dfgrd1(3,1)*a0(1) + dfgrd1(3,2)*a0(2) + dfgrd1(3,3)*a0(3)
        da(1) = Fa0(1)*Fa0(1)
        da(2) = Fa0(2)*Fa0(2)
        da(3) = Fa0(3)*Fa0(3)
        da(4) = Fa0(1)*Fa0(2)
        da(5) = Fa0(1)*Fa0(3)
        da(6) = Fa0(2)*Fa0(3)

        lam = sqrt(a0(1)*(cmat(1,1)*a0(1)+cmat(1,2)*a0(2)+cmat(1,3)*a0(3)) 
     #            +a0(2)*(cmat(1,2)*a0(1)+cmat(2,2)*a0(2)+cmat(2,3)*a0(3)) 
     #            +a0(3)*(cmat(1,3)*a0(1)+cmat(2,3)*a0(2)+cmat(3,3)*a0(3))  )
        
        lam_m = lam_m_n(i) 
c        print *, lam_m
        IF (lam > lam_m) THEN
            statev(i) = lam
	  END IF

        PKfiber = 0.0

c...    Integration over the Weibull

        n = 100		! number of divisions for trapezoidal method
        a = gama	! lower limit, gamma
        b = lam		! upper limit, lamba, unless lam is less than gama 
        if (lam<gama) then
          b = gama+1e-5
        end if

        delta = (c2*lam_m)+c1
        p_lams_a = (beta/delta)*((a-gama)/delta)**(beta - 1)*exp(-((a-gama)/delta)**beta)
        p_lams_b = (beta/delta)*((b-gama)/delta)**(beta - 1)*exp(-((b-gama)/delta)**beta)
        f_a = p_lams_a*(((lam/a)**2) - (1.0/(lam/a)))
        f_b = p_lams_b*(((lam/b)**2) - (1.0/(lam/b)))

	  f_Edam_a = p_lams_a*(((lam/a)**2) + 2*(1.0/(lam/a)) - 3.0)     ! *Added for Damaged Energy Calculation
	  f_Edam_b = p_lams_b*(((lam/b)**2) + 2*(1.0/(lam/b)) - 3.0)     ! *Added for Damaged Energy Calculation

        h = (b-a)/n
        tot_sum = 0.5*(f_a + f_b)
        tot_sum_Edam = 0.5*(f_Edam_a + f_Edam_b)                       ! *Added for Damaged Energy Calculation

        do m = 1, n-1
          p_lams = (beta/delta)*(((a+m*h)-gama)/delta)**(beta-1)*exp(-(((a+m*h)-gama)/delta)**beta)
          lam_bar = lam/(a+m*h)
          f = p_lams*((lam_bar**2)-(1.0/lam_bar))         
          tot_sum = tot_sum + f

          f_Edam = p_lams*((lam_bar**2) + 2*(1.0/lam_bar) - 3.0)       ! *Added for Damaged Energy Calculation
          tot_sum_Edam = tot_sum_Edam + f_Edam                         ! *Added for Damaged Energy Calculation

        end do
    
        integral = h*tot_sum
        PKfiber = mu1*(1.0/(lam**2))*integral

        integral_Edam = h*tot_sum_Edam                                 ! *Added for Damaged Energy Calculation
        E_fiber = mu1*integral_Edam                                    ! *Added for Damaged Energy Calculation

c...    Calculating Undamaged Energy

        v = 100                                                           ! *Added for Undamaged Energy Calculation
        delta0 = c1 + c2                                                  ! *Added for Undamaged Energy Calculation
        p0_lams_a = (beta/delta0)*((a-gama)/delta0)**(beta - 1)*exp(-((a-gama)/delta0)**beta)     ! *Added for Undamaged Energy Calculation
        p0_lams_b = (beta/delta0)*((b-gama)/delta0)**(beta - 1)*exp(-((b-gama)/delta0)**beta)     ! *Added for Undamaged Energy Calculation

	  f_Eundam_a = p0_lams_a*(((lam/a)**2) + 2*(1.0/(lam/a)) - 3.0)     ! *Added for Undamaged Energy Calculation
	  f_Eundam_b = p0_lams_b*(((lam/b)**2) + 2*(1.0/(lam/b)) - 3.0)     ! *Added for Undamaged Energy Calculation

        tot_sum_Eundam = 0.5*(f_Eundam_a + f_Eundam_b)                    ! *Added for Undamaged Energy Calculation

        do m = 1, v-1
          p_lams = (beta/delta0)*(((a+m*h)-gama)/delta0)**(beta-1)*exp(-(((a+m*h)-gama)/delta0)**beta)
          lam_bar = lam/(a+m*h)
          f_Eundam = p_lams*((lam_bar**2) + 2*(1.0/lam_bar) - 3.0)        ! *Added for Undamaged Energy Calculation
          tot_sum_Eundam = tot_sum_Eundam + f_Eundam                      ! *Added for Undamaged Energy Calculation
        end do

        integral_Eundam = h*tot_sum_Eundam                                ! *Added for Undamaged Energy Calculation
        Eun_fiber = mu1*integral_Eundam                                   ! *Added for Undamaged Energy Calculation

c...    End Calculation Undamaged Energy ..........................................................................

c...    Von mises according to wikipedia
c       VM = exp(kappa*cos(2*(theta-mutheta)))/(pi*I0kappa)
        VM = (exp(kappa*cos(2*(theta-mutheta))))/(2*pi*I0kappa)

c...    ABT Jan 9, 2021, debugging 
c...    The push forward of the fiber stress is actually just the deformed fiber 
c...    Note the 1/1 is from 1/J, but J=1 because we can impose incompressible exactly
c...    JT Jan 14, 2021 Added dtheta missing.
c...    VDS Jul 25, 2022 Replacing 1/1 with 1/J for this case since incompressibilty is enfoced with penalty
        sigmaf(1) = sigmaf(1) + (1.0/detf)*VM*PKfiber*da(1)*dtheta 
        sigmaf(2) = sigmaf(2) + (1.0/detf)*VM*PKfiber*da(2)*dtheta 
        sigmaf(3) = sigmaf(3) + (1.0/detf)*VM*PKfiber*da(3)*dtheta  
        sigmaf(4) = sigmaf(4) + (1.0/detf)*VM*PKfiber*da(4)*dtheta  
        sigmaf(5) = sigmaf(5) + (1.0/detf)*VM*PKfiber*da(5)*dtheta  
        sigmaf(6) = sigmaf(6) + (1.0/detf)*VM*PKfiber*da(6)*dtheta 

        E_damaged = E_damaged + (1.0/1.0)*E_fiber*VM*dtheta 
        E_undamaged = E_undamaged + (1.0/1.0)*Eun_fiber*VM*dtheta

      end do

c...      Edam_Eundam = 1.0 - (E_damaged/E_undamaged)
      Edam_Eundam = E_undamaged - E_damaged                ! Added for integral Von Mises 1/14/21 
      statev(11) = Edam_Eundam

c...  ABT edit nov 27, 2020
c...  Need to calculate b = F*F^T 
      bvoigt(1) = dfgrd1(1,1)*dfgrd1(1,1) + dfgrd1(1,2)*dfgrd1(1,2) + dfgrd1(1,3)*dfgrd1(1,3)
      bvoigt(2) = dfgrd1(2,1)*dfgrd1(2,1) + dfgrd1(2,2)*dfgrd1(2,2) + dfgrd1(2,3)*dfgrd1(2,3)
      bvoigt(3) = dfgrd1(3,1)*dfgrd1(3,1) + dfgrd1(3,2)*dfgrd1(3,2) + dfgrd1(3,3)*dfgrd1(3,3)
      bvoigt(4) = dfgrd1(1,1)*dfgrd1(2,1) + dfgrd1(1,2)*dfgrd1(2,2) + dfgrd1(1,3)*dfgrd1(2,3)
      bvoigt(5) = dfgrd1(1,1)*dfgrd1(3,1) + dfgrd1(1,2)*dfgrd1(3,2) + dfgrd1(1,3)*dfgrd1(3,3)
      bvoigt(6) = dfgrd1(2,1)*dfgrd1(3,1) + dfgrd1(2,2)*dfgrd1(3,2) + dfgrd1(2,3)*dfgrd1(3,3)
c...  Calculating some more tensors to check stress calculated in both reference and deformed
c...  See if they match, they do! Turned off the fiber part, so it is only the Neo-Hookean
      bmat(1,1) = bvoigt(1)
      bmat(2,2) = bvoigt(2)
      bmat(3,3) = bvoigt(3)
      bmat(1,2) = bvoigt(4)
      bmat(2,1) = bvoigt(4)
      bmat(1,3) = bvoigt(5)
      bmat(3,1) = bvoigt(5)
      bmat(2,3) = bvoigt(6)
      bmat(3,2) = bvoigt(6)
c...  pressure
c...  Note, in plane strain case going back to the nearly incompressible case
      biso(1) = detf**(-2./3.)*bvoigt(1)
      biso(2) = detf**(-2./3.)*bvoigt(2)
      biso(3) = detf**(-2./3.)*bvoigt(3)
      biso(4) = detf**(-2./3.)*bvoigt(4)
      biso(5) = detf**(-2./3.)*bvoigt(5)
      biso(6) = detf**(-2./3.)*bvoigt(6)
      I1bar = biso(1)+biso(2)+biso(3)
      sigmabar(1) = (1.0/detf)*(mu0*biso(1) )
      sigmabar(2) = (1.0/detf)*(mu0*biso(2) )
      sigmabar(3) = (1.0/detf)*(mu0*biso(3) )
      sigmabar(4) = (1.0/detf)*(mu0*biso(4) )
      sigmabar(5) = (1.0/detf)*(mu0*biso(5) )
      sigmabar(6) = (1.0/detf)*(mu0*biso(6) )
c...  sigmaiso = Projection::sigmabar, super simple in eulerian 
      tr_sigmabar = sigmabar(1)+sigmabar(2)+sigmabar(3)
      sigmaiso(1) = sigmabar(1) -(1./3.)*tr_sigmabar
      sigmaiso(2) = sigmabar(2) -(1./3.)*tr_sigmabar
      sigmaiso(3) = sigmabar(3) -(1./3.)*tr_sigmabar
      sigmaiso(4) = sigmabar(4)
      sigmaiso(5) = sigmabar(5)
      sigmaiso(6) = sigmabar(6)
c...  volumetric part (we can change this later)
c...  sigmavol = J*p*Identity, with p = dPsivol/dJ
      Psivol = K_vol*((detf-1)**2)
      p = 2*K_vol*(detf-1)
      dpdJ = 2*K_vol
c... full stress tensor 
      sigma(1) = sigmaiso(1) + sigmaf(1) + p 
      sigma(2) = sigmaiso(2) + sigmaf(2) + p
      sigma(3) = sigmaiso(3) + sigmaf(3) + p
      sigma(4) = sigmaiso(4) + sigmaf(4)
      sigma(5) = sigmaiso(5) + sigmaf(5)
      sigma(6) = sigmaiso(6) + sigmaf(6)
c...      print *, 'Sigma in 3D' 
c...      print *, sigma(1), sigma(2), sigma(3), sigma(4), sigma(5), sigma(6)
c     FULL stress, because this is plane strain the sigma_33 might be non-zero
      sigmamat(1,1) = sigma(1)
      sigmamat(2,2) = sigma(2)
      sigmamat(3,3) = sigma(3)
      sigmamat(1,2) = sigma(4)
      sigmamat(2,1) = sigma(4)
      sigmamat(1,3) = sigma(5)
      sigmamat(3,1) = sigma(5)
      sigmamat(2,3) = sigma(6)
      sigmamat(3,2) = sigma(6)

c...  2D stress
      sigma2D(1) = sigma(1)
      sigma2D(2) = sigma(2)
      sigma2D(3) = sigma(4)
      sigmaiso2D(1) = sigmaiso(1)
      sigmaiso2D(2) = sigmaiso(2)
      sigmaiso2D(3) = sigmaiso(4)
      sigmamat2D(1,1) = sigma(1)
      sigmamat2D(2,2) = sigma(2)
      sigmamat2D(1,2) = sigma(4)
      sigmamat2D(2,1) = sigma(4)
c      print *, 'Sigma 2D' 
c      print *, sigma2D(1), sigma2D(2), sigma2D(3)

c...  FIRST PART OF TANGENT
c...  Need the kron delta in voigt
      kron(1) = 1.0
      kron(2) = 1.0
      kron(3) = 1.0
      kron(4) = 0.0
      kron(5) = 0.0
      kron(6) = 0.0
      kronmat(1,1) = 1.0
      kronmat(2,2) = 1.0
      kronmat(3,3) = 1.0
      kronmat(1,2) = 0.0
      kronmat(2,1) = 0.0
      kronmat(1,3) = 0.0
      kronmat(3,1) = 0.0
      kronmat(2,3) = 0.0
      kronmat(3,2) = 0.0
      Itoi(1) = 1
      Itoi(2) = 2
      Itoi(3) = 3
      Itoi(4) = 1
      Itoi(5) = 1
      Itoi(6) = 2
      Itoj(1) = 1
      Itoj(2) = 2
      Itoj(3) = 3
      Itoj(4) = 2
      Itoj(5) = 3
      Itoj(6) = 3
c...  2D versions 
      kron2D(1) = 1.0
      kron2D(2) = 1.0
      kron2D(3) = 0.0
      kronmat2D(1,1) = 1.0
      kronmat2D(2,2) = 1.0
      kronmat2D(1,2) = 0.0
      kronmat2D(2,1) = 0.0
      Itoi2D(1) = 1
      Itoi2D(2) = 2
      Itoi2D(3) = 1
      Itoj2D(1) = 1
      Itoj2D(2) = 2
      Itoj2D(3) = 2

c...  Fill part of the tangent in voigt notation
c...  Reduced tangent, non fiber part
c...  1->11, 2->22, 3->12
c...  Itoi2D = [1,2,1] (definition above)
c...  Itoj2D = [1,2,2] (definition above)
      do II=1,4
        do JJ=1,4
c...   The tangent should be the same as the 3D case but just take the corresponding
c...   components, because anything to do with the z coordinate is simply fixed, so just 
c...   eliminating rows and columns, there is no coupling like in plane stress
c...   VDS edit july 25, 2022,  Notice all tangents are 4x4 (11,22,33,12)
          q = Itoi(II)
          r = Itoj(II)
          s = Itoi(JJ)
          t = Itoj(JJ)
       
c...    ABT edit july 20, 2022
c...    Changing to spatial fourth order tensor, 
          IIII = 0.5*(kronmat(q,s)*kronmat(r,t)+kronmat(q,t)*kronmat(r,s))
          III = kronmat(q,r)*kronmat(s,t)
          cciso(II,JJ) = (2.0/3.0)*tr_sigmabar*(IIII+(1.0/3.0)*kron(II)*kron(JJ))
     #        -(2.0/3.0)*(kron(II)*sigmaiso(JJ)+sigmaiso(II)*kron(JJ))
c...        print *, 'CCiso from directly deformed'
c...        print *, cciso(II,JJ)
          ccvol(II,JJ) = (p + detf*dpdJ)*kron(II)*kron(JJ) - 2.0*p*IIII
c...        print *, 'CCvol directly deformed'
c...        print *, ccvol(II,JJ)
        end do
      end do
   
c...  SECOND PART OF THE TANGENT
      do II=1,4
        do JJ=1,4
          ccf(II,JJ) = 0.0
        end do
      end do

      do i=1,10
        theta = (i*pi/10.0)
        dtheta =  0.3141592
        a0(1) = cos(theta)
        a0(2) = sin(theta)
        a0(3) = 0.0
c...    Deformed fiber a = F*a0
        Fa0(1) = dfgrd1(1,1)*a0(1) + dfgrd1(1,2)*a0(2) + dfgrd1(1,3)*a0(3)
        Fa0(2) = dfgrd1(2,1)*a0(1) + dfgrd1(2,2)*a0(2) + dfgrd1(2,3)*a0(3)
        Fa0(3) = dfgrd1(3,1)*a0(1) + dfgrd1(3,2)*a0(2) + dfgrd1(3,3)*a0(3)
        da(1) = Fa0(1)*Fa0(1)
        da(2) = Fa0(2)*Fa0(2)
        da(3) = Fa0(3)*Fa0(3)
        da(4) = Fa0(1)*Fa0(2)
        da(5) = Fa0(1)*Fa0(3)
        da(6) = Fa0(2)*Fa0(3)
        da2D(1) = da(1)
        da2D(2) = da(2)
        da2D(3) = da(4)

        lam = sqrt(a0(1)*(cmat(1,1)*a0(1)+cmat(1,2)*a0(2)+cmat(1,3)*a0(3)) 
     #            +a0(2)*(cmat(1,2)*a0(1)+cmat(2,2)*a0(2)+cmat(2,3)*a0(3)) 
     #            +a0(3)*(cmat(1,3)*a0(1)+cmat(2,3)*a0(2)+cmat(3,3)*a0(3))  )
    
        lam_m = lam_m_n(i)

        CCfiber1 = 0.0
	    CCfiber2 = 0.0

c...    Integration over the Weibull

        n = 100		! number of divisions for trapezoidal method
        a = gama	! lower limit, gamma
        b = lam		! upper limit, lambda, unless lam < gamma in which case basically integrate zero 
        if (lam<gama) then
          b = gama + 1e-5
        end if

        delta = c1+(c2*lam_m)
        p_lams_a = (beta/delta)*((a-gama)/delta)**(beta-1)*exp(-((a-gama)/delta)**beta)
        p_lams_b = (beta/delta)*((b-gama)/delta)**(beta-1)*exp(-((b-gama)/delta)**beta)
    	
        f_a1 = p_lams_a*(((lam/a)**2) - (1.0/(lam/a)))
        f_b1 = p_lams_b*(((lam/b)**2) - (1.0/(lam/b)))
        h = (b-a)/n
        tot_sum1 = 0.5*(f_a1 + f_b1)

        do l = 1, n-1
          p_lams = (beta/delta)*(((a+l*h)-gama)/delta)**(beta-1)*exp(-(((a+l*h)-gama)/delta)**beta)
          lam_bar = lam/(a+l*h)
          f1 = p_lams*((lam_bar**2)-(1.0/lam_bar))         
          tot_sum1 = tot_sum1 + f1
        end do
    
        integral1 = h*tot_sum1

        f_a2 = p_lams_a*((2*(lam/a)**2) + (1.0/(lam/a)))
        f_b2 = p_lams_b*((2*(lam/b)**2) + (1.0/(lam/b)))
        tot_sum2 = 0.5*(f_a2 + f_b2)

        do m = 1, n-1
          p_lams = (beta/delta)*(((a+m*h)-gama)/delta)**(beta-1)*exp(-(((a+m*h)-gama)/delta)**beta)
          lam_bar = lam/(a+m*h)
          f2 = p_lams*((2*lam_bar**2) + (1.0/lam_bar))         
          tot_sum2 = tot_sum2 + f2
        end do
    
        integral2 = h*tot_sum2

        CCfiber1 = -2*mu1*(1.0/(lam**4))*integral1
        CCfiber2 =  1*mu1*(1.0/(lam**4))*integral2

	
        ccfiber = CCfiber1 + CCfiber2

c...    Von mises according to wikipedia
c...    https://en.wikipedia.org/wiki/Von_Mises_distribution
c        VM = exp(kappa*cos(2*(theta-mutheta)))/(pi*I0kappa)
        VM = (exp(kappa*cos(2*(theta-mutheta))))/(2*pi*I0kappa)

	
      do II=1,4
          do JJ=1,4

            q = Itoi(II)
            r = Itoj(II)
            s = Itoi(JJ)
            t = Itoj(JJ)

c...   JT edit April 28. I will extract the 11,22 and 12 terms after the loop into the new ccf matrix
c...   I'm changing to 2D even from now. Added dtheta missing edit January 14th, 2021
            ccf(II,JJ) = ccf(II,JJ) + ccfiber*VM*Fa0(q)*Fa0(r)*Fa0(s)*Fa0(t)*dtheta
          end do
        end do
      end do


c...  ABT April 2021, just arranging everything to 2D
      do II=1,4
         STRESS(II) = sigma(II)
        do JJ=II,4
          q = Itoi(II)
          r = Itoj(II)
          s = Itoi(JJ)
          t = Itoj(JJ)

c...  Abaqus corrections
c...  ABT aprl 2021, also changing the different cciso, ccvol, etc to just cc 
          ddsdde(II,JJ) = cciso(II,JJ)+ccvol(II,JJ)+ccf(II,JJ)+0.5*(kronmat(q,s)*sigmamat(r,t)
     #                    +kronmat(q,t)*sigmamat(r,s)+kronmat(r,s)*sigmamat(q,t)
     #                    +kronmat(r,t)*sigmamat(q,s))

          if (JJ>II) then
            ddsdde(JJ,II) = ddsdde(II,JJ)
          end if
        end do
      end do
c...      print *, 'STRESS'
c...      print *, stress(1), stress(2), stress(3), stress(4), stress(5), stress(6)
c...      print *, 'DDSDDE'
c...      print *, ddsdde(1,1), ddsdde(1,2), ddsdde(1,3), ddsdde(1,4), ddsdde(1,5), ddsdde(1,6)
c...      print *, ddsdde(2,1), ddsdde(2,2), ddsdde(2,3), ddsdde(2,4), ddsdde(2,5), ddsdde(2,6)
c...      print *, ddsdde(3,1), ddsdde(3,2), ddsdde(3,3), ddsdde(3,4), ddsdde(3,5), ddsdde(3,6)
c...      print *, ddsdde(4,1), ddsdde(4,2), ddsdde(4,3), ddsdde(4,4), ddsdde(4,5), ddsdde(4,6)
c...      print *, ddsdde(5,1), ddsdde(5,2), ddsdde(5,3), ddsdde(5,4), ddsdde(5,5), ddsdde(5,6)
c...      print *, ddsdde(6,1), ddsdde(6,2), ddsdde(6,3), ddsdde(6,4), ddsdde(6,5), ddsdde(6,6)
      return
      end

c...  ------------------------------------------------------------------

      subroutine cross(aa, bb,cc1)
      implicit none

      real*8 :: cc1(3)
      real*8 :: aa(3), bb(3)

      cc1(1) = aa(2) * bb(3) - aa(3) * bb(2)
      cc1(2) = aa(3) * bb(1) - aa(1) * bb(3)
      cc1(3) = aa(1) * bb(2) - aa(2) * bb(1)

      return
      end
c...  ------------------------------------------------------------------

c...  ------------------------------------------------------------------


c...  ------------------------------------------------------------------
      end
c...  ------------------------------------------------------------------
