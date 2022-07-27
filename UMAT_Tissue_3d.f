c... Copyright(c) 2022, John Toaquiza Tubon, Omar Moreno Flores and Adrian Buganza Tepole
c... All rights reserved.
c...
c... Redistribution and use in source and binary forms, with or without
c... modification, are permitted provided that the following conditions are met:
c... 
c... 1. Redistributions of source code must retain the above copyright notice, this
c... list of conditions and the following disclaimer. 
c... 2. Redistributions in binary form must reproduce the above copyright notice,
c...    this list of conditions and the following disclaimer in the documentation
c...    and/or other materials provided with the distribution.

c... THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
c... ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
c... WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
c... DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
c... ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
c... (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
c... LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
c... ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
c... (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c... SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
      statev(11)=1.0d0

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
c... --------------------------------------------------------------------- 
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

      implicit none

c...  variables to be defined
      real*8  stress(ntens), ddsdde(ntens,ntens), statev(nstatv), sse

c...  variables passed in for information
      real*8  time(2), dtime, coords(3), props(nprops), dfgrd1(3,3)
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt, kstep, kinc

c...  local variables (mostly mechanics part)
      real*8 finv(3,3), detf, c(6), lnJ
      real*8 cmat(3,3), kronmat(3,3)
      real*8 I1b, p, Psivol, dpdJ, a0(3), da(6), I1bar		
      real*8 PKfiber						
      real*8 CCfiber1, CCfiber2, ccfiber
      real*8 tr_sigmabar, Jccbar(6,6), kron(6), IIII
      real*8 cciso(6,6), ccvol(6,6), ccf(6,6), cc(6,6)
      real*8 E_damaged, E_undamaged, Edam_Eundam
      real*8 bvoigt(6), bmat(3,3), Fa0(3), biso(6), sigmaiso(6), sigma(6), sigmamat(3,3)
      real*8 sigmabar(6), sigmaf(6) 

c...  some auxiliar variables, tensors
      integer i, j, l, m, n, q, r, s, t, v, w, nitl, II, JJ, Itoi(6), Itoj(6)

c...  material properties
      real*8 mu0, mu1, kappa, mutheta, I0kappa, gama, beta, c1, c2, lam_m, K_vol

c...  variables and constants for Weibull Distribution
      real*8 pi, theta, delta, lam, VM, delta0, deltad, dtheta
      real*8 a, b, p_lams_a, p_lams_b, f_a, f_b
      real*8 p0_lams_a, p0_lams_b, f_Edam_a, f_Edam_b, tot_sum_Edam, f_Edam, integral_Edam, E_fiber
      real*8 pd_lams_a, pd_lams_b, f_Eundam_a, f_Eundam_b, tot_sum_Eundam, f_Eundam, integral_Eundam, Eun_fiber    
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
c      deltad = delta

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

c      print *, 'deformation gradient'
c      print *, dfgrd1(1,1),dfgrd1(1,2),dfgrd1(1,3)
c      print *, dfgrd1(2,1),dfgrd1(2,2),dfgrd1(2,3)
c      print *, dfgrd1(3,1),dfgrd1(3,2),dfgrd1(3,3)
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
      pi = 3.1415926
      do i=1,10
        theta = (i*pi/10.0)
        dtheta = 0.3141592 
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

c       ABT jan 9, 2021, trying to debug 
        lam = sqrt(a0(1)*(cmat(1,1)*a0(1)+cmat(1,2)*a0(2)+cmat(1,3)*a0(3)) 
     #            +a0(2)*(cmat(1,2)*a0(1)+cmat(2,2)*a0(2)+cmat(2,3)*a0(3)) 
     #            +a0(3)*(cmat(1,3)*a0(1)+cmat(2,3)*a0(2)+cmat(3,3)*a0(3))  )
        
        lam_m = lam_m_n(i)
        IF (lam .GT. lam_m) THEN
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

        delta = c1 + c2*lam_m
        deltad = c1 + c2*lam_m
        p_lams_a = (beta/delta)*((a-gama)/delta)**(beta-1)*exp(-((a-gama)/delta)**beta)
        p_lams_b = (beta/delta)*((b-gama)/delta)**(beta-1)*exp(-((b-gama)/delta)**beta)
      	f_a = p_lams_a*(((lam/a)**2) - (1.0/(lam/a)))
      	f_b = p_lams_b*(((lam/b)**2) - (1.0/(lam/b)))

        h = (b-a)/n
        tot_sum = 0.5*(f_a + f_b)

        do m = 1, n-1
          p_lams = (beta/delta)*(((a+m*h)-gama)/delta)**(beta-1)*exp(-(((a+m*h)-gama)/delta)**beta)
          lam_bar = lam/(a+m*h)
          f = p_lams*((lam_bar**2)-(1.0/lam_bar))         
          tot_sum = tot_sum + f
        end do
    
        integral = h*tot_sum
      	PKfiber = 2*mu1*(1.0/(lam**2))*integral

c...    Calculating Damaged Energy
        w = 100
        pd_lams_a = (beta/deltad)*((a-gama)/deltad)**(beta - 1)*exp(-((a-gama)/deltad)**beta)     ! *Added for Damaged Energy Calculation
        pd_lams_b = (beta/deltad)*((b-gama)/deltad)**(beta - 1)*exp(-((b-gama)/deltad)**beta)     ! *Added for Damaged Energy Calculation
	f_Edam_a = pd_lams_a*(((lam/a)**2) + 2*(1.0/(lam/a)) - 3.0)    ! *Added for Damaged Energy Calculation
	f_Edam_b = pd_lams_b*(((lam/b)**2) + 2*(1.0/(lam/b)) - 3.0)    ! *Added for Damaged Energy Calculation

        tot_sum_Edam = 0.5*(f_Edam_a + f_Edam_b)                       ! *Added for Damaged Energy Calculation

        do m = 1, w-1
          p_lams = (beta/deltad)*(((a+m*h)-gama)/deltad)**(beta-1)*exp(-(((a+m*h)-gama)/deltad)**beta)
          lam_bar = lam/(a+m*h)

          f_Edam = p_lams*((lam_bar**2) + 2*(1.0/lam_bar) - 3.0)       ! *Added for Damaged Energy Calculation
          tot_sum_Edam = tot_sum_Edam + f_Edam                         ! *Added for Damaged Energy Calculation

        end do

        integral_Edam = h*tot_sum_Edam                                 ! *Added for Damaged Energy Calculation
        E_fiber = mu1*integral_Edam                                    ! *Added for Damaged Energy Calculation

c...    End Calculating Damaged Energy


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
c...    https://en.wikipedia.org/wiki/Von_Mises_distribution
        VM = exp(kappa*cos(2*(theta-mutheta)))/(2*pi*I0kappa)

c...    ABT Jan 9, 2021, debugging 
c...    The push forward of the fiber stress is actually just the deformed fiber 
        sigmaf(1) = sigmaf(1) + (1.0/detf)*VM*PKfiber*da(1)*dtheta
        sigmaf(2) = sigmaf(2) + (1.0/detf)*VM*PKfiber*da(2)*dtheta
        sigmaf(3) = sigmaf(3) + (1.0/detf)*VM*PKfiber*da(3)*dtheta 
        sigmaf(4) = sigmaf(4) + (1.0/detf)*VM*PKfiber*da(4)*dtheta 
        sigmaf(5) = sigmaf(5) + (1.0/detf)*VM*PKfiber*da(5)*dtheta 
        sigmaf(6) = sigmaf(6) + (1.0/detf)*VM*PKfiber*da(6)*dtheta
 
        E_damaged = E_damaged + E_fiber*VM*dtheta
        E_undamaged = E_undamaged + Eun_fiber*VM*dtheta

      end do
c...      Edam_Eundam = 1.0 - (E_damaged/E_undamaged)
      Edam_Eundam = E_undamaged-E_damaged
      statev(11) = Edam_Eundam

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

c...  ABT edit nov 27, 2020
c...  Need to calculate b = F*F^T for the eulerian tangent 
      bvoigt(1) = dfgrd1(1,1)*dfgrd1(1,1) + dfgrd1(1,2)*dfgrd1(1,2) + dfgrd1(1,3)*dfgrd1(1,3)
      bvoigt(2) = dfgrd1(2,1)*dfgrd1(2,1) + dfgrd1(2,2)*dfgrd1(2,2) + dfgrd1(2,3)*dfgrd1(2,3)
      bvoigt(3) = dfgrd1(3,1)*dfgrd1(3,1) + dfgrd1(3,2)*dfgrd1(3,2) + dfgrd1(3,3)*dfgrd1(3,3)
      bvoigt(4) = dfgrd1(1,1)*dfgrd1(2,1) + dfgrd1(1,2)*dfgrd1(2,2) + dfgrd1(1,3)*dfgrd1(2,3)
      bvoigt(5) = dfgrd1(1,1)*dfgrd1(3,1) + dfgrd1(1,2)*dfgrd1(3,2) + dfgrd1(1,3)*dfgrd1(3,3)
      bvoigt(6) = dfgrd1(2,1)*dfgrd1(3,1) + dfgrd1(2,2)*dfgrd1(3,2) + dfgrd1(2,3)*dfgrd1(3,3)
c...  Calculating some more tensors to check stress calculated in both reference and deformed
c...  See if they match, they do! Turned off the fiber part, so it is only the Neo-Hookean
c...  part that I will be compering
      biso(1) = detf**(-2./3.)*bvoigt(1)
      biso(2) = detf**(-2./3.)*bvoigt(2)
      biso(3) = detf**(-2./3.)*bvoigt(3)
      biso(4) = detf**(-2./3.)*bvoigt(4)
      biso(5) = detf**(-2./3.)*bvoigt(5)
      biso(6) = detf**(-2./3.)*bvoigt(6)
      bmat(1,1) = bvoigt(1)
      bmat(2,2) = bvoigt(2)
      bmat(3,3) = bvoigt(3)
      bmat(1,2) = bvoigt(4)
      bmat(2,1) = bvoigt(4)
      bmat(1,3) = bvoigt(5)
      bmat(3,1) = bvoigt(5)
      bmat(2,3) = bvoigt(6)
      bmat(3,2) = bvoigt(6)
c...  Will calculate the stress for neo hookean directly in the deformed and see if the two match 
c      I1bar = detf**(-2./3.)*I1c
      I1bar = biso(1)+biso(2)+biso(3)
      sigmabar(1) = (1.0/detf)*(2*mu0*biso(1) )
      sigmabar(2) = (1.0/detf)*(2*mu0*biso(2) )
      sigmabar(3) = (1.0/detf)*(2*mu0*biso(3) )
      sigmabar(4) = (1.0/detf)*(2*mu0*biso(4) )
      sigmabar(5) = (1.0/detf)*(2*mu0*biso(5) )
      sigmabar(6) = (1.0/detf)*(2*mu0*biso(6) )
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
      Psivol = K_vol*(detf-1)**2
      p = 2*K_vol*(detf-1)
      dpdJ = 2*K_vol
      sigma(1) = sigmaiso(1) + sigmaf(1) + p 
      sigma(2) = sigmaiso(2) + sigmaf(2) + p
      sigma(3) = sigmaiso(3) + sigmaf(3) + p
      sigma(4) = sigmaiso(4) + sigmaf(4)
      sigma(5) = sigmaiso(5) + sigmaf(5)
      sigma(6) = sigmaiso(6) + sigmaf(6)
c      print *, 'Sigma directly in deformed' 
c      print *, sigma(1), sigma(2), sigma(3), sigma(4), sigma(5), sigma(6)
      sigmamat(1,1) = sigma(1)
      sigmamat(2,2) = sigma(2)
      sigmamat(3,3) = sigma(3)
      sigmamat(1,2) = sigma(4)
      sigmamat(2,1) = sigma(4)
      sigmamat(1,3) = sigma(5)
      sigmamat(3,1) = sigma(5)
      sigmamat(2,3) = sigma(6)
      sigmamat(3,2) = sigma(6)

c...  Fill part of the tangent in voigt notation
c...  1->11, 2->22, 3->33, 4->12, 5->13, 6->23
c...  Itoi = [1,2,3,1,1,2] (definition above)
c...  Itoj = [1,2,3,2,3,3] (definition above)
      do II=1,6
        do JJ=1,6
c...   There is ont term in the tangent which requires special tensor product
          q = Itoi(II)
          r = Itoj(II)
          s = Itoi(JJ)
          t = Itoj(JJ)
       
c...    ABT edit nov 27 
c...    Changing to spatial fourth order tensor, 
          IIII = 0.5*(kronmat(q,s)*kronmat(r,t)+kronmat(q,t)*kronmat(r,s))
          cciso(II,JJ) = (2.0/3.0)*tr_sigmabar*(IIII-(1.0/3.0)*kron(II)*kron(JJ))
     #       -(2.0/3.0)*(kron(II)*sigmaiso(JJ)+sigmaiso(II)*kron(JJ))
c...        print *, 'CCiso from directly deformed'
c...        print *, cciso(II,JJ)
          ccvol(II,JJ) = (p + detf*dpdJ)*kron(II)*kron(JJ) - 2.0*p*IIII
c...        print *, 'CCvol directly deformed'
c...        print *, ccvol(II,JJ)
c...    There seemse to be a factor of 2 missing in ccvol
        end do
      end do
      
      do II=1,6
        do JJ=1,6
          ccf(II,JJ) = 0.0
        end do
      end do

      do i=1,10
        theta = (i*pi/10.0)
        dtheta = 0.3141592
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

        CCfiber1 = 0.0
	CCfiber2 = 0.0

c...    Integration over the Weibull

        n = 100		! number of divisions for trapezoidal method
        a = gama	! lower limit, gamma
        b = lam		! upper limit, lambda, unless lam < gamma in which case basically integrate zero 
        if (lam<gama) then
          b = gama + 1e-5
        end if

        delta = c1 + c2*lam_m
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

       	CCfiber1 = -4*mu1*(1.0/(lam**4))*integral1
       	CCfiber2 =  2*mu1*(1.0/(lam**4))*integral2

	
        ccfiber = CCfiber1 + CCfiber2

c...    Von mises according to wikipedia
c...    https://en.wikipedia.org/wiki/Von_Mises_distribution
        VM = exp(kappa*cos(2*(theta-mutheta)))/(2*pi*I0kappa)

	
      	do II=1,6
          do JJ=1,6

            q = Itoi(II)
            r = Itoj(II)
            s = Itoi(JJ)
            t = Itoj(JJ)

c...   ABT edit Nov 27 
c...   Commenting this line
c            ccf(II,JJ) = ccf(II,JJ) + ccfiber*VM*a0(q)*a0(r)*a0(s)*a0(t)
c...   elasticity tensor pushed to the deformed configuration is almost the same but with deformed fiber
            ccf(II,JJ) = ccf(II,JJ) + (1.0/detf)*ccfiber*VM*Fa0(q)*Fa0(r)*Fa0(s)*Fa0(t)*dtheta

          end do
        end do
      end do
  
      do II=1,6
c        STRESS(II)=CS(II)
         STRESS(II) = sigma(II)
        do JJ=II,6
          q = Itoi(II)
          r = Itoj(II)
          s = Itoi(JJ)
          t = Itoj(JJ)

          cc(II,JJ) = cciso(II,JJ) + ccvol(II,JJ) + ccf(II,JJ)

c...  Abaqus corrections
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
c...  calculate strain energy
      sse = Psivol + mu0*(I1bar-3)

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
