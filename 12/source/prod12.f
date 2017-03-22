c234567

!     =================================================================
!     Productivity1: Photosynthesis, Plant Respiration and NPP
!     Carbon2: Microbial Respiration

!     Code written by David Lapola and Helena Alves do Prado
!     revised and modified by jpdarela
!     LAST UPDATE: Ter 17 Jan 2017 00:37:49 BRST JP 
!     =================================================================

      subroutine productivity1 (pft,ocp_pft,ligth_limit,temp,p0,w,
     &    wmax,ca,ipar,rh,cl1,ca1,cf1,beta_leaf,beta_awood,      
     &    beta_froot,emax,ph,ar,nppa,laia,f5,f1,vpd,rm,rml, 
     &    rmf,rms,rg,rgl,rgf,rgs,rc)

      implicit none

!     Variables
!     =========
!     Input
!     -----     
      integer pft
      real ocp_pft              !PFT area occupation (%)
      real temp                 !Mean monthly temperature (oC)
      real p0                   !Mean surface pressure (hPa)
      real wa,w,wmax            !Soil moisture (dimensionless)
      real ca                   !Atmospheric CO2 concentration (Pa)
      real ipar                 !Incident photosynthetic active radiation (w/m2)'
      real rh                   !Relative humidity
      real cl1, cf1, ca1        !Carbon in plant tissues (kg/m2)
      real beta_leaf            !npp allocation to carbon pools (kg/m2/day)
      real beta_awood
      real beta_froot
      logical ligth_limit       !True for no ligth limitation
      
!     Output
!     ------
      real ph                   !Canopy gross photosynthesis (kgC/m2/yr)
      real rc                   !Stomatal resistence (not scaled to canopy!) (s/m)
      real laia                 !Autotrophic respiration (kgC/m2/yr)
      real ar                   !Leaf area index (m2 leaf/m2 area)
      real nppa, vpd            !Net primary productivity (kgC/m2/yr)
      real f5                   !Water stress response modifier (unitless)
      real*16 :: f5_64          !f5 auxiliar (more precision) 
      real rm, rml, rmf, rms    !autothrophic respiration (kgC/m2/day)
      real rg, rgl, rgf, rgs 
      
!     Internal
!     --------
      real es, es2              !Saturation partial pressure (hPa) == mbar
      real aux_ipar
      double precision vm       !Rubisco maximum carboxylaton rate (molCO2/m2/s)
      double precision mgama    !Photo-respiration compensation point (Pa)
      double precision rmax     !Saturated mixing ratio (kg/kg)
      double precision r        !Moisture deficit at leaf level (kg/kg)
      double precision ci       !Internal leaf CO2 partial pressure (Pa)
      double precision a,b,c,a2,b2,c2 !Auxiliars
      double precision delta,delta2 !Auxiliars
      double precision sunlai,shadelai !Sunlai/Shadelai
      double precision vpd64,d, vpd_real !Vapor pressure deficit (kPa)
      double precision laia64, ph64, ar64, rm64, rms64, rml64
      double precision rmf64, rg64, rgf64, rgs64, rgl64
      DOUBLE PRECISION  leaf_t_months,leaf_turnover,leaf_t_coeff
            
      
!     Rubisco, light and transport limited photosynthesis rate
!     --------------------------------------------------------
      double precision jc,jl,je,jp,jp1,jp2,j1,j2 !Auxiliars
      
!     Functions
!     ---------
      double precision f1       !Leaf level gross photosynthesis (molCO2/m2/s)
      double precision f1a      !auxiliar_f1
      double precision f2       !Michaelis-Menton CO2 constant (Pa)
      double precision f3       !Michaelis-Menton O2 constant (Pa)
      double precision f4,f4sun,f4shade !Scaling-up LAI to canopy level (dimensionless)
      
      real tleaf(12)             !leaf turnover time (yr)
      real p21(12)
      real emax                 !potential evapotranspiration (mm/dia)
      double precision sla      !specific leaf area (m2/kg)
      double precision csa      !sapwood compartment´s carbon content (5% of woody tissues) (kgC/m2)
      double precision ncl      !leaf N:C ratio (gN/gC)
      double precision ncf      !fine roots N:C ratio (gN/gC)
      double precision ncs      !sapwood N:C ratio(gN/gC)
      double precision csai 
      double precision pt       !taxa potencial de fornecimento para transpiração (mm/dia)
      double precision csru     !Specific root water uptake (0.5 mm/gC/dia; based in jedi
      double precision alfm     !maximum Priestley-Taylor coefficient (based in Gerten et al. 2004
      double precision gm       !scaling conductance (dia/mm)(based in Gerten et al. 2004;
      double precision gc       !Canopy conductance (dia/mm)(based in Gerten et al. 2004;

!     MODEL PARAMETERS
      DOUBLE PRECISION p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14
     &    ,p15,p19,p20,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31
     
      a   = 0.8300              !Photosynthesis co-limitation coefficient
      a2  = 0.930               !Photosynthesis co-limitation coefficient
      p3  = 21200.0             !Atmospheric oxygen concentration (Pa)
      p4  = 0.080               !Quantum efficiency (mol electrons/Ein)
      p5  = 0.150               !Light scattering rate
      p6  = 2.0                 !Parameter for jl
      p7  = 0.50                !Ratio of light limited photosynthesis to Rubisco carboxylation
      p8  = 5200.0              !Photo-respiration compensation point
      p9  = 0.570               !Photosynthesis co-limitation coefficient
      p10 = 0.100               !Q10 function
      p11 = 25.0                !Q10 function reference temperature (oC)
      p12 = 30.0                !Michaelis-Menten constant for CO2 (Pa)
      p13 = 2.100               !Michaelis-Menten constant for CO2
      p14 = 30000.0             !Michaelis-Menten constant for O2 (Pa)
      p15 = 1.20                !Michaelis-Menten constant for O2
      p19 = 0.90                !Maximum ratio of internal to external CO2
      p20 = 0.10                !Critical humidity deficit (kg/kg)
      p22 = 2.0                 !Rubisco carboxylation rate
      p23 = 0.30                !Rubisco carboxylation rate
      p24 = 36.0                !Rubisco carboxylation rate (oC)
      p25 = 0.000008            !Maximum gross photosynthesis rate (molCO2/m2/s)
      p26 = 0.50                !light extinction coefficient for IPAR/sun (0.5/sen90)
      p27 = 1.50                !light extinction coefficient for IPAR/shade (0.5/sen20)
      p28 = 0.500               !Soil moisture at field capacity
      p29 = 0.205               !Soil moisture at wilting point
      p30 = 0.015               !Ratio of respiration to Rubisco carboxylation rates
      p31 = 3.850               !Whole plant to leaf respiration ratio


!     getting pft parameters
      call pft_par(2, p21)
      call pft_par(6, tleaf)

!     ==============
!     Photosynthesis 
!     ==============
      
!     Rubisco maximum carboxylaton rate (molCO2/m2/s)
!     -----------------------------------------------
      vm = (p21(pft))
c      vm = (p21(pft)*(p22**(p10*(temp-p11))))/ !Free-range parameter --> 0.0358>vm>840 (micromol)
c     &    (1.0+exp(p23*(temp-p24)))

      
!     Photo-respiration compensation point (Pa)
!     -----------------------------------------
      mgama = p3/(p8*(p9**(p10*(temp-p11))))

      
!     Michaelis-Menton CO2 constant (Pa)
!     ----------------------------------     
      f2 = p12*(p13**(p10*(temp-p11)))


!     Michaelis-Menton O2 constant (Pa)
!     ---------------------------------
      f3 = p14*(p15**(p10*(temp-p11)))

      
!     Saturation partial pressure of water vapour (Pa)
      call tetens(temp,es)

!     Saturated mixing ratio (kg/kg)
!     -----------------------------     
      rmax = 0.622*(es/(p0-es))
     
!     Moisture deficit at leaf level (kg/kg)
!     --------------------------------------    
      r = -0.315*rmax

!     Internal leaf CO2 partial pressure (Pa)
!     ---------------------------------------     
      ci = p19* (1.-(r/p20)) * (ca-mgama) + mgama

!     Rubisco carboxilation limited photosynthesis rate (molCO2/m2/s)
!     ---------------------------------------------------------------     
      jc = vm*((ci-mgama)/(ci+(f2*(1.+(p3/f3)))))

!     Light limited photosynthesis rate (molCO2/m2/s)
!     -----------------------------------------------
      if (ligth_limit) then
         aux_ipar= ipar   
      else
         aux_ipar = ipar - (ipar * 0.20)
      endif
      jl = p4*(1.0-p5)*aux_ipar*((ci-mgama)/(ci+(p6*mgama)))
     
!     Transport limited photosynthesis rate (molCO2/m2/s)
!     ---------------------------------------------------
      je = p7*vm

!     Jp (minimum between jc and jl)
!     ------------------------------   
      a = 0.83
      b = (-1.)*(jc+jl)
      c = jc*jl
      delta = (b**2)-4.0*a*c
    
      jp1=(-b-(sqrt(delta)))/(2.0*a)
      jp2=(-b+(sqrt(delta)))/(2.0*a)
      jp= amin1(jp1,jp2)
         
!     Leaf level gross photosynthesis (minimum between jc, jl and je)
!     ---------------------------------------------------------------
      a2 = 0.93
      b2 = (-1.)*(jp+je)
      c2 = jp*je
      delta2 = (b2**2)-4.0*a2*c2

      j1=(-b2-(sqrt(delta2)))/(2.0*a2)
      j2=(-b2+(sqrt(delta2)))/(2.0*a2)
      f1a = amin1(j1,j2)

!     VPD
!     ===
!     buck equation...references:
!     http://www.hygrometers.com/wp-content/uploads/CR-1A-users-manual-2009-12.pdf
!     Hartmann 1994 - Global Physical Climatology p.351
!     https://en.wikipedia.org/wiki/Arden_Buck_equation#CITEREFBuck1996
      
!     ES2 = VPD-POTENTIAL - Saturation Vapor Pressure (millibar)
      call tetens(temp, es2)
      
!     VPD-REAL = Actual vapor pressure
      vpd_real = es2 * rh       ! RESULTS are IN hPa == mbar! we want kPa (DIVIDE by 10.)

!     Vapor Pressure Deficit
      vpd64 = (es2 - VPD_REAL) / 10. 
      vpd = real(vpd64, 4)
      
!     Soil water
!     ==========      
      wa = (w/wmax)
      
!     Stomatal resistence
!     ===================
      call canopy_resistence(pft , vpd64, f1a, rc)
      
!     Water stress response modifier (dimensionless)
!     ----------------------------------------------

c     Water stress response modifier (dimensionless)
c     [f5 ; Eq. 21]
c      ! vamos deixar o F5 como antigamente ate este problema ser resolvido
c      if (wa.gt.0.5) then
c         f5 = 1.0               !Not too lower in e.g. Amazonian dry season
c      else if ((wa.ge.0.205).and.(wa.le.0.5)) then
c         f5 = (wa-0.205)/(0.5-0.205)
c      else if (wa.lt.0.205) then
c         f5 = wa !Below wilting point f5 accompains wa (then Sahara is well represented)
c      endif

      
      csru = 0.5 
      pt = csru*(cf1*1000.)*wa  !(based in Pavlick et al. 2013; *1000. converts kgC/m2 to gC/m2)
      alfm = 1.391
      gm = 3.26 * 86400.           !(*86400 transform s/mm to dia/mm)
c      
      if(rc .gt. 0.001) then
         gc = rc * 1.15741e-08 ! transfor s/m  to dia/mm)  !testamos, nao muda nada! Bia vai rever
         gc = (1./gc)  ! molCO2/mm2/dia
      else
         gc =  1.0/0.001 ! BIANCA E HELENA - Mudei este esquema..   
      endif                     ! tentem entender o algoritmo
                                ! e tenham certeza que faz sentido ecologico
      d =(emax*alfm)/(1. + gm/gc) !(based in Gerten et al. 2004)
!     BIanca- Eu alterei a estrutura desta equacao para evitar erros
!     Isso faz a mesma coisa que o calculo que vc implementou - jp
      if(d .gt. 0.0) then
         f5_64 = pt/d
         f5_64 = exp(-1 * f5_64)
         f5_64 = 1.0 - f5_64
      else
         f5_64 = wa  ! eu acrescentei esta parte caso d seja igual a zero
      endif          ! nao sei se faz sentido!

      f5 = real(f5_64,4) !esta funcao transforma o f5 (double precision)
!     em real (32 bits) 
      
!     Photosysthesis minimum and maximum temperature
!     ----------------------------------------------
      
      if ((temp.ge.-10.0).and.(temp.le.50.0)) then
c         f1 = f1a*real(f5_64,8) !f5:water stress factor-- Notem que aqui a tranformacao eh de 128 pra 64 bits
          f1 = f1a * f5
      else
         f1 = 0.0               !Temperature above/below photosynthesis windown
      endif

!     Leaf area index (m2/m2)
      leaf_t_months = tleaf(pft)*12. ! turnover time in months
      leaf_t_coeff = leaf_t_months/100. !1 - 100 months == ~ 1/12 to 8.3 years (TRY-kattge et al. 2011; Jedi-Pavlick 2012)
c      if(leaf_t_months .gt. 0) print*, leaf_t_months, leaf_t_coeff, pft
      if (leaf_t_coeff .gt. 1.) leaf_t_coeff = 1. 
      leaf_turnover =  (365.0/12.0) * (10. **(2.0*leaf_t_coeff))
      sla = (3e-2 * (365.0/leaf_turnover)**(-0.46))
      
      

      
      laia64 = ((cl1*1000.)  * sla) ! * 1000 transform kg to g - laia64 in m2 m-2
c      if(laia64 .gt. 0.0) print*, laia64, 'jp'
!     ---------------------------------
!     Specifc Leaf Area----------------
!      sla=((0.0300*1000.)*((365./(((tleaf(pft))/365.)/12.))**(-0.46)))    
!      laia64 = (cl1 * 365.0 * sla)


c      sla=(0.03*(365/(tleaf(pft)/365))**(-0.46))
c      laia64 = (cl1 * 1000 * 365.0 * sla)
c      if(laia64 .gt. 0.0) print*, laia64, 'bia'
!     LAI
!     ------
      sunlai = (1.0-(exp(-p26*laia64)))/p26
!     --------
      shadelai = laia64 - sunlai

!     Scaling-up to canopy level (dimensionless)
!     ------------------------------------------
      f4 = (1.0-(exp(-p26*laia64)))/p26 !Sun 90 degrees in the whole canopy, to be used for respiration
      
!     Sun/Shade approach to canopy scaling !Based in de Pury & Farquhar (1997)
!     ------------------------------------------------------------------------
      f4sun = ((1.0-(exp(-p26*sunlai)))/p26) !sun 90 degrees
      f4shade = ((1.0-(exp(-p27*shadelai)))/p27) !sun ~20 degrees

      laia  = real((f4sun + f4shade), 4) !real(laia64,4) ! pra mim faz sentido que a laia final seja
                                         ! a soma da lai em nivel de dossel (sun + shade) - jp
!     Canopy gross photosynthesis (kgC/m2/yr)
!     =======================================
!     (0.012 converts molCO2 to kgC)
!     (31557600 converts seconds to year [with 365.25 days])
      ph64 = 0.012*31557600.0*f1*f4sun*f4shade
      ph = real(ph64, 4) * ocp_pft       ! kg m-2 year-1
c      write(*,*) 'ph', ph
c     ============================================================
c     Autothrophic respiration
!     ========================
!     Maintenance respiration (kgC/m2/yr) (based in Ryan 1991)
      csa= 0.05 * ca1           !sapwood carbon content (kgC/m2). 5% of woody tissues (Pavlick, 2013)
            
      ncl = 1./29.               !(gN/gC) from lpj3 
      ncf = 1./29.               !(gN/gC)
      ncs = 1./330.              !(gN/gC)
 
      rml64 = (ncl * cl1) * 27. * exp(0.035*temp)
      rml =  real(rml64,4)

      rmf64 = (ncf * cf1) * 27. * exp(0.035*temp)
      rmf =  real(rmf64,4)

      rms64 = (ncs * csa) * 27. * exp(0.035*temp)
      rms = real(rms64,4) 
       
      rm64 = (rml64 + rmf64 + rms64)
      rm = real(rm64, 4)

c     Growth respiration (KgC/m2/yr)(based in Ryan 1991; Sitch et al.
c     2003; Levis et al. 2004)         
       
      csai =  (beta_awood * 0.05)

      rgl64 = 0.25 * beta_leaf * 365.
      rgl = real(rgl64,4) 

      rgf64 =  0.25* beta_froot * 365.
      rgf = real(rgf64,4) 

      rgs64 = (0.25 * csai * 365.)
      rgs = real(rgs64,4) 
     
      rg64 = (rgl64 + rgf64 + rgs64)
      rg = real(rg64,4) 
     
      if (rg.lt.0) then
         rg = 0.0
      endif
      
!     c Autotrophic (plant) respiration -ar- (kgC/m2/yr)
!     Respiration minimum and maximum temperature
!     -------------------------------------------     
      if ((temp.ge.-10.0).and.(temp.le.50.0)) then
         ar64 = rm64 + rg64
         ar = real(ar64,4) * ocp_pft
         if(ar .lt. 0.00001) ar = 0.0
      else
         ar = 0.0               !Temperature above/below respiration windown
      endif
      
c     -----------------------------------------------------------------
!     NPP
!     ============
!     Productivity
!     ============
!     Net primary productivity(kgC/m2/yr)
!     ====================================
      nppa = ph - ar
      if(nppa .lt. 0.0) nppa = 0.0 
!     No futuro proximo poderiamos usar uma npp negativa para indicar uma perda de carbono
!     dos tecidos vegetais... e uma possivel extincao do pft/pls da celula de grid.
      return
      end subroutine productivity1



      
!     ==================================================================
      subroutine canopy_resistence(m,vpd_in,f1_in,rc2_in)
      implicit none
!     Variables
!     =========
      
!     Inputs
!     ------      
      integer m
      double precision f1_in    !Photosynthesis (molCO2/m2/s)
      double precision vpd_in   !kPa
      
!     Outputs
!     -------
      real rc2_in               !Canopy resistence (s/m)
      
!     Internal
!     --------
      real g1(12)
      real rcmax, rcmin
      double precision f1b      !Photosynthesis (micromolCO2/m2/s)
      double precision gs2      !Canopy conductance (m/s)
      double precision gs       !Canopy conductance (molCO2/m2/s)
      double precision g0       !Residual stomatance conductance
      double precision D1       !kPA
      double precision aa

      call pft_par(1, g1)
      
      f1b = (f1_in*10e5)        ! Helena - Mudei algumas coisas aqui
      aa = (f1b/363.)           ! Entenda o algoritmo e tenha certeza de que  
      g0 = 0.01                 ! condiz com a realidade esperada =)
      rcmax = 553.000
      rcmin = 100.000

      if(f1_in .le. 0.0) then 
         rc2_in = rcmax
         goto 110
      endif
      if(vpd_in .gt. 0.1) then
         goto 10
      else
         rc2_in = rcmin
         goto 110
      endif
 10   continue
      if (vpd_in .gt. 0.95) then
         rc2_in = rcmax
         goto 110
      else
         D1 = sqrt(vpd_in)
         gs = g0 + 1.6 * (1. + (g1(m)/D1)) * (aa) !Based on Medlyn et al. 2011
         if(gs .le. 0.0) then
            rc2_in = rcmax
            goto 110
         endif
      endif
      
      gs2 = (gs/41.)
      if(gs2 .le. 0.0) rc2_in = rcmax
      if(gs2 .gt. 0.0) rc2_in = real((gs2**(-1)),4)
 110  continue
      
      return
      end subroutine canopy_resistence
      
!     =====================================
!     Microbial (heterotrophic) respiration
!     =====================================
      
      subroutine carbon2 (tsoil,f5c,evap,laia,d_litter,     !Inputs
     &    cl,cs,hr)                                         !Outputs

      implicit none

!     Variables
!     =========
!     
!     Inputs
!     ------
      
      real tsoil                !Mean monthly soil temperature (oC)
      real f5c                  !Stress response to soil moisture (dimensionless)
      real evap                 !Actual evapotranspiration (mm/day)
      real laia
      real d_litter

!     xOutputs 
!     -------
!     
      real cl                   !Litter carbon (kgC/m2)
      real cs                   !Soil carbon (kgC/m2)
      real hr                   !Heterotrophic (microbial) respiration (kgC/m2)

!     Internal
!     --------
!     
      real lf                   !Litterfall (kgC/m2)
      real f6                   !Litter decayment function
      real f7                   !Soil carbon storage function
      real p10                  !Q10 function
      real p11                  !Q10 function reference temperature (oC)
      real p32                  !Q10 parameter of soil respiration sensibility to temperature
      real p33                  !Litterfall (kg/C/yr)
      real p34                  !Average fraction of litter carbon lost to atmosphere
      real p35                  !Carbon soil turnover (1/20yr)
      real p36                  !Specific rate heterotrophic respiration
!     
!     Initialize
!     ----------
!     
      lf  = 0.0
      f6  = 0.0
      f7  = 0.0
      cl  = 0.0
      cs  = 0.0
!     
      p10 = 0.10
      p11 = 25.0
      p32 = 2.00
      p33 = 0.10
      p34 = 0.30
      p35 = 0.05
      p36 = 0.25
!     
!     Litter decayment function                                             !Controlled by annual evapotranspiration
!     -------------------------
!     
      f6 = 1.16*10.**(-1.4553+0.0014175*(evap*365.0))
C      call critical_value(f6)
!     
!     Soil carbon storage function                                          !Controlled by temperature
!     ----------------------------
!     
      f7 = p32**(p10*(tsoil-p11))
C      call critical_value(f7)
!     Litterfall (kgC/m2)
!     ------------------
!     
      lf = p33 * (laia + d_litter)
C      call critical_value(lf)
!     
!     Litter carbon (kgC/m2)
!     ----------------------
!     
      cl = lf/f6
C      call critical_value(cl)
!     
!     Soil carbon(kgC/m2)
!     -------------------
!     
      cs = ((p34*cl)/(p35*f7))*f5c
C      call critical_value(cs)
!     
!     Respiration minimum and maximum temperature
!     -------------------------------------------
!     
      if ((tsoil.ge.-10.0).and.(tsoil.le.50.0)) then
         hr = p36*(cl*(f6**2)+(cs*f5c*evap*(f7**2))) !Litter and Soil respectively
C         call critical_value(hr)
      else
         hr = 0.0               !Temperature above/below respiration windown
      endif
!     
      return
      end
!     ===================================================

!     Estas 3 subrotinas critical_value eu utilizei durante o debug
!     Ainda tem algumas espalhadas pelo codigo, pretendo
!     ir retirando aos poucos agora, mas no futuro eu vou utiliza-las
      
      subroutine critical_value(var)
      implicit none
      real var
      
      if(abs(var) .lt. 0.000001) var = 0.0

      return
      end subroutine critical_value


      
      subroutine critical_value2(var1)
      implicit none
      double precision var1
      
      if(abs(var1) .lt. 0.00000001) var1 = 0.0

      return
      end subroutine critical_value2


      
      subroutine critical_value3(var2)
      implicit none
      real*16 var2
      
      if(abs(var2) .lt. 0.000000000001) var2 = 0.0

      return
      end subroutine critical_value3
!     =================================================================

c23456
      SUBROUTINE PFT_AREA_FRAC(CLEAF, CFROOT, CAWOOD, OCP_COEFFS,
     &    OCP_WOOD)

      IMPLICIT NONE
 
      INTEGER, PARAMETER :: NPFT = 12
      INTEGER :: P, MAX_INDEX(1), I
      REAL :: CLEAF(NPFT), CFROOT(NPFT), CAWOOD(NPFT)
      REAL :: TOTAL_BIOMASS_PFT(NPFT),OCP_COEFFS(NPFT)
      REAL :: TOTAL_BIOMASS, TOTAL_WOOD, TOTAL_W_PFT(NPFT)
      LOGICAL :: OCP_WOOD(NPFT)
      
      TOTAL_BIOMASS = 0.0
      TOTAL_WOOD = 0.0
      do p = 1,npft
         TOTAL_W_PFT(P) = 0.0
         TOTAL_BIOMASS_PFT(P) = 0.0
         OCP_COEFFS(P) = 0.0
         OCP_WOOD(P) = .FALSE.
      enddo
      

      DO P = 1,NPFT
         TOTAL_BIOMASS_PFT(P) = CLEAF(P) + CFROOT(P) + CAWOOD(P)
         TOTAL_BIOMASS = TOTAL_BIOMASS + TOTAL_BIOMASS_PFT(P)
         TOTAL_WOOD = TOTAL_WOOD + CAWOOD(P)
         TOTAL_W_PFT(P) = CAWOOD(P)
      ENDDO

!     GRID CELL OCCUPATION COEFFICIENTS
      IF(TOTAL_BIOMASS .GT. 0.0) THEN
         DO P = 1,NPFT   
            OCP_COEFFS(P) = TOTAL_BIOMASS_PFT(P) / TOTAL_BIOMASS
            IF(OCP_COEFFS(P) .LT. 0.0) OCP_COEFFS(P) = 0.0
         ENDDO
      ELSE
         DO P = 1,NPFT
            OCP_COEFFS(P) = 0.0
         ENDDO
      ENDIF

!     GRIDCELL PFT LIGTH LIMITATION BY WOOD CONTENT 
      IF(TOTAL_WOOD .GT. 0.0) THEN
         MAX_INDEX = MAXLOC(TOTAL_W_PFT)
         I = MAX_INDEX(1)
         OCP_WOOD(I) = .TRUE.
      ENDIF

      RETURN
      END SUBROUTINE PFT_AREA_FRAC
!     =========================================================
     
      subroutine penman (spre,temp,ur,rn,rc2,evap)
      implicit none
!     Inputs
!     ------
!     
      real spre                 !Surface pressure (mb)
      real temp                 !Temperature (oC)
      real ur                   !Relative humidity (0-1,dimensionless)
      real rn                   !Radiation balance (W/m2)
      real rc2                  !Canopy resistence (s/m)
!     
!     
!     Output
!     ------
!     
      real evap                 !Evapotranspiration (mm/day)
!     
!     Parameters
!     ----------
      real ra, h5, t1, t2, es, es1, es2, delta_e, delta
      real gama, gama2
      ra = 100.                  !s/m
      h5 = 0.0275               !mb-1
!     
!     Delta
!     -----
!     
      t1 = temp + 1.
      t2 = temp - 1.
      call tetens(t1,es1)       !Saturation partial pressure of water vapour at temperature T
      call tetens(t2,es2)
      delta = (es1-es2)/(t1-t2) !mb/oC
!     
!     Delta_e
!     -------
!     
      call tetens (temp,es)
      delta_e = es*(1. - ur)    !mb
!     
      if ((delta_e.ge.(1./h5)-0.5).or.(rc2.ge.550.0)) evap = 0.
      if ((delta_e.lt.(1./h5)-0.5).or.(rc2.lt.550.0)) then
!     
!     Gama and gama2
!     --------------
         gama  = spre*(1004.)/(2.45e6*0.622)
         gama2 = gama*(ra + rc2)/ra

!     Real evapotranspiration
!     -----------------------     
         evap = (delta* rn + (1.20*1004./ra)*delta_e)/(delta+gama2) !W/m2
         evap = evap*(86400./2.45e6) !mm/day
         evap = amax1(evap,0.)  !Eliminates condensation
      endif

      return
      end subroutine penman     
!     ============================================
!     
      subroutine evpot2 (spre,temp,ur,rn,evap) 
      implicit none
!     Inputs
!     ------
!     
      real spre                 !Surface pressure (mb)
      real temp                 !Temperature (oC)
      real ur                   !Relative humidity (0-1,dimensionless)
      real rn                   !Irradiation balance (W/m2)
!     
!     Output
!     ------
!     
      real evap                 !Potencial evapotranspiration without stress (mm/day)
!     
!     Parameters
!     ----------
      real ra, rcmin, t1, t2, es, es1, es2, delta_e, delta
      real gama, gama2, rc
!     
      ra      = 100.            !s/m
      rcmin   = 100.            !s/m
!     
!     Delta
!     -----
!     
      t1 = temp + 1.
      t2 = temp - 1.
      call tetens(t1,es1)
      call tetens(t2,es2)
      delta = (es1-es2)/(t1-t2) !mb/oC
!     
!     Delta_e
!     -------
!     
      call tetens (temp,es)
      delta_e = es*(1. - ur)    !mb
!     
!     Stomatal Conductance
!     --------------------
!     
      rc = rcmin
!     
!     Gama and gama2
!     --------------
!     
      gama  = spre*(1004.)/(2.45e6*0.622)
      gama2 = gama*(ra + rc)/ra
!     
!     Potencial evapotranspiration (without stress)
!     ---------------------------------------------
!     
      evap = (delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) !W/m2
      evap = evap*(86400./2.45e6) !mm/day
      evap = amax1(evap,0.)     !Eliminates condensation
!     
      return
      end
!     
!     =================================================================
!     ===
!     
      subroutine runoff (wa,roff)
      implicit none
      real*8 :: wa
      real :: roff
      real*16 :: roff64
      roff64 = 38.*(wa**11.)
c      roff64 = 11.5*(wa**6.6) * 1000. !From NCEP-NCAR Reanalysis data
      call critical_value3(roff64)
      roff = real(roff64, 4)
      return
      end
!     
!     =================================================================
!     ====

     
      subroutine tetens (t,es)  !Saturation Vapor Pressure (hPa)!
      implicit none
      real t,es
      if (t.ge.0.) then
         es = 6.1121 * exp((18.678-(t/234.5))*(t/(257.14+t))) ! mbar == hPa
      else
         es = 6.1115 * exp((23.036-(t/333.7))*(t/(279.82+t))) ! mbar == hPa
      endif
c      es = es/10. ! converting to kPa
!     
      return
      end        
!     =============================================================
c=====================================================================
c     
c     subroutine allocation calculates the daily carbon content of each
c     compartment
c     
c     code written by Bianca Rius & David Lapola (27.Ago.2015)
c     
c=====================================================================
      
      subroutine allocation (pft, npp ,scl1,sca1,scf1,
     &    scl2,sca2,scf2,bio_litter)          !output
      implicit none
c     
!     variables
      integer, parameter :: npfts = 12
      integer pft   
      real npp                  !potential npp (KgC/m2/yr)
      real*16 npp_aux           !auxiliary variable to calculate potential npp in KgC/m2/day
      real scl1                  !previous day carbon content on leaf compartment (KgC/m2)
      real scl2                  !final carbon content on leaf compartment (KgC/m2)
      real sca1                  !previous day carbon content on aboveground woody biomass compartment(KgC/m2)
      real sca2                  !final carbon content on aboveground woody biomass compartment (KgC/m2)
      real scf1                  !previous day carbon content on fine roots compartment (KgC/m2)
      real scf2                  !final carbon content on fine roots compartment (KgC/m2)
      real bio_litter            !carbon that goes to litter when the carbon is less than 1e-12 (KgC/m2)     
      real*16 scf2_128, sca2_128, scl2_128
      real aleaf(npfts)             !npp percentage allocated compartment
      real aawood(npfts)
      real afroot(npfts)
      real tleaf(npfts)             !turnover time (yr)
      real tawood(npfts)
      real tfroot(npfts)            

      call pft_par(3, aleaf)
      call pft_par(4, aawood)
      call pft_par(5, afroot)
      call pft_par(6, tleaf)
      call pft_par(7, tawood)
      call pft_par(8, tfroot)
    
!     Carbon content of each compartment(KgC/m2)
c     
c     
c     initialization
      if((scl1 .lt. 1e-12) .or. (scf1 .lt. 1e-12)) then
         bio_litter = scl1 + scf1 + sca1
            scl2 = 0.0
            scf2 = 0.0
            sca2 = 0.0 
            goto 10
         
      endif   
      npp_aux = npp/365.0       !transform (KgC/m2/yr) in (KgC/m2/day)
      scl2_128 = scl1 + (aleaf(pft) * npp_aux) -(scl1 /(tleaf(pft)
     &    *365.0))
      scf2_128 = scf1 +(afroot(pft) * npp_aux)-(scf1 /(tfroot(pft)
     &    *365.0))
      if(aawood(pft) .gt. 0.0) then
         sca2_128 = sca1 +(aawood(pft)*npp_aux)-(sca1/(tawood(pft)
     &       *365.0))
         sca2 = real(sca2_128,4)
      else
         sca2 = 0.0
      endif

      scf2 = real(scf2_128,4)
      scl2 = real(scl2_128,4)


      if(scl2 .lt. 1e-12) scl2 = 0.0
      if(scf2 .lt. 1e-12) scf2 = 0.0
      if(sca2 .lt. 1e-12) sca2 = 0.0
      
 10   continue
      return
      end
c
