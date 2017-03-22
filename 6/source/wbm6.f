C23456
      subroutine wbm (prec,temp,lsmk,p0,ca,par,rhs,cleaf_ini,cawood_ini
     &    ,cfroot_ini,emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft
     &    ,clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,evapm_pft
     &    ,wsoil_pft,rml_pft,rmf_pft,rms_pft,rm_pft,rgl_pft,rgf_pft
     &    ,rgs_pft,rg_pft,cleaf_pft,cawood_pft,cfroot_pft,grid_area
     &    ,betal,betaw,betaf)
      

      implicit none


!     =================================================================
!     Water balance model (WBM).
!     From monthly climatologies of precipitation and surface
!     temperature, 
!     the WBM calculates the environmental variables.
!     
!     05 Jul 2005, MDO: ae, rh & runoff are changed
!     11 Jul 2005, MDO: wsoil2 is written (for testing purpose)
!     31 Ago 2006, DML: carbon cycle is included
!     31 dez 2016, DML, HAP, BR and JP- several changes        
!     =================================================================

!     Variables
!     =========
!     
      integer ,parameter :: nx=120,ny=160,q=6
      real ,parameter :: no_data = -9999.0
      
c     --------------------------I N P U T S----------------------------
      real ca                   !CO2 atmospheric concentration (ppmv)
      real p0(nx,ny,12)         !Atmospheric pressure (mb)
      real lsmk(nx,ny)          !Land=1/Ocean=0
      real prec(nx,ny,12)       !Precipitation (mm/month)
      real temp(nx,ny,12)       !Temperature (oC)
      real par(nx,ny,12)        !IPAR (Ein/m2/s)
      real rhs(nx,ny,12)        !Relative humidity

      real  cleaf_ini(nx,ny,q)  ! Initial carbon content in leaves (kg m-2)
      real  cawood_ini(nx,ny,q) ! Initial carbon content in aboveground wood (kg m-2)
      real  cfroot_ini(nx,ny,q) ! Initial carbon content in fineroots (kg m-2)
C     -----------------------------E N D-------------------------------
      
c     -------------------------O U T P U T S---------------------------
       
      real tsoil (nx,ny,12)     !soil temperature
      
      real photo_pft(nx,ny,12,q) !Monthly photosynthesis   (kgC m-2)
      real aresp_pft(nx,ny,12,q) !Monthly autotrophic res  (kgC m-2)
      real npp_pft  (nx,ny,12,q) !Monthly net primary produ (kgC m-2)
      
      real lai_pft  (nx,ny,12,q) !Monthly leaf area index
      real clit_pft (nx,ny,12,q) !Monthly litter carbon
      real csoil_pft(nx,ny,12,q) !Monthly soil carbon
      real hresp_pft(nx,ny,12,q) !Monthly het resp  (kgC/m2)
      real rcm_pft  (nx,ny,12,q) 

      real emaxm    (nx,ny,12) !Max.evapotranspiration (kg m-2 day-1)
      real runom_pft(nx,ny,12,q) !Runoff 
      real evapm_pft(nx,ny,12,q) !Actual evapotranspiration        
      real wsoil_pft(nx,ny,12,q) !Soil moisture (mm)
      
c     NOVOS OUTPUTS DA BIANCA
      real rml_pft  (nx,ny,12,q) ! Maintenance respiration 
      real rmf_pft  (nx,ny,12,q)
      real rms_pft  (nx,ny,12,q)
      real rm_pft   (nx,ny,12,q)
      
      real rgl_pft  (nx,ny,12,q) ! Growth respiration
      real rgf_pft  (nx,ny,12,q)
      real rgs_pft  (nx,ny,12,q)
      real rg_pft   (nx,ny,12,q)
      
      real cleaf_pft (nx,ny,q) ! leaf biomass (KgC/m2)
      real cawood_pft(nx,ny,q) ! aboveground wood biomass (KgC/m2)
      real cfroot_pft(nx,ny,q)  ! fine root biomass
      
      real grid_area(nx,ny,q)   ! gridcell area fraction of pfts!
      real betal(nx,ny,12,q)
      real betaw(nx,ny,12,q)
      real betaf(nx,ny,12,q)

c     --------------------------------E N D----------------------------

c     ------------------------- internal variables---------------------

      integer i, j, k, kk, n, mes, nerro, p
      
      
      real wg0  (nx,ny,12)      !Previous year soil moisture
      real wsoilt(nx,ny,12)     !soil water (check wbm equilibrium)
      real gsoilt(nx,ny,12)     !soil ice   (check wbm equilibrium)
      real gsoil(nx,ny,12,q)    !Soil ice 
      real ssoil(nx,ny,12,q)    !Soil snow
      real snowm(nx,ny,12,q)    !Snowmelt

      real sini(q), sfim(q)
      real wini(q), wfim(q)
      real gini(q), gfim(q)

      real cleaf1_pft (q)
      real cawood1_pft(q)
      real cfroot1_pft(q)
      real wood(q)

!     outputs for budget
      real epmes ! equal for all pfts - Maximum evapotranspiration (mm day-1)
      real rmes(q),phmes(q),smes(q),rcmes(q),hrmes(q)
      real nppmes(q),laimes(q), armes(q),clmes(q),csmes(q)
      real emes(q),rmlmes(q),rmfmes(q),rmsmes(q)
      real rmmes(q),rglmes(q),rgfmes(q),rgsmes(q),rgmes(q)
      real cleafmes(q),cawoodmes(q),cfrootmes(q), gridocpmes(q)
      real betalmes(q), betawmes(q), betafmes(q)
      
      real t0,t1                !Soil temperature aux variables 
      real wsaux1,dwww,wmax     !auxiliar to equilibrium check
      real ae                   !Available energy     
      real pr,spre,ta,td,ipar,ru
      
      real leaf0(q), froot0(q),awood0(q)
      real biomass, biomass0
      real sensi
      logical check

      real, parameter :: H = 1.0 !Soil layer(m) 
      real, parameter :: diffu = 4.e-7*(30.*86400.0) !Soil thermal diffusivity (m2/month)
      real, parameter :: tau = (H**2)/(2.0*diffu) !E-folding time (months)
      real, parameter :: auxs = -100.0 !Auxiliar for calculation of Snpp

      
            

!     ================      
!     Soil temperature
!     ================

!     For all grid points
!     -------------------
      do i=1,nx
         do j=1,ny

!     Initialize soil temperature
!     ---------------------------
            do k=1,12
               tsoil(i,j,k) = no_data
            enddo

!     Only for land grid points
!     -------------------------!     
            if (nint(lsmk(i,j)).ne.0) then
               t0 = 0.          !Initialization
               do n = 1,1200    !1200 months run to attain equilibrium
                  k = mod(n,12)
                  if (k.eq.0) k = 12
                  t1 = t0*exp(-1.0/tau)+(1.0-exp(-1.0/tau))*temp(i,j,k) 
                  tsoil(i,j,k) = (t0 + t1)/2.0
                  t0 = t1
               enddo
            endif     
         enddo
      enddo

!     ============
!     Water budget
!     ============
!     
!     For all grid points
!     -------------------     
      do i=1,nx
         do j=1,ny
            do p = 1,q   
               cleaf1_pft (p) =  cleaf_ini(i,j,p)
               cawood1_pft(p) = cawood_ini(i,j,p)
               cfroot1_pft(p) = cfroot_ini(i,j,p)
               cleaf_pft(i,j,p)  = no_data ! leaf biomass (KgC/m2)
               cawood_pft(i,j,p) = no_data ! aboveground biomass (KgC/m2)
               cfroot_pft(i,j,p) = no_data ! fine root biomass (KgC/m2)
               grid_area(i,j,p)  = no_data ! gridcell area fraction of pfts(%)
            enddo

c     Write to track program execution     
            if ((mod(j,ny).eq.0).and.(mod(i,10).eq.0))
     &          write(*,*) 'working: ', (real(i)/real(nx))*100.0, '%' 
          
!     Initialize variables
            do k=1,12
               wg0  (i,j,k) = no_data !Soil water content in preceeding year integration
               emaxm(i,j,k) = no_data !Maximum evapotranspiration
               do p=1,q
                  photo_pft(i,j,k,p)  = no_data !Monthly photosynthesis (kgC/m2)
                  aresp_pft(i,j,k,p)  = no_data !Monthly autotrophic respiration (kgC/m2)
                  npp_pft(i,j,k,p)    = no_data !Monthly net primary productivity (average between PFTs) (kgC/m2)
                  lai_pft(i,j,k,p)    = no_data !Monthly leaf area index
                  clit_pft(i,j,k,p)   = no_data !Monthly litter carbon
                  csoil_pft(i,j,k,p)  = no_data !Monthly soil carbon
                  hresp_pft(i,j,k,p)  = no_data !Monthly heterotrophic respiration (kgC/m2)
                  rcm_pft(i,j,k,p)    = no_data
                  gsoil(i,j,k,p)      = no_data !Soil ice
                  ssoil(i,j,k,p)      = no_data !Soil snow
                  runom_pft(i,j,k,p)  = no_data !Runoff
                  evapm_pft(i,j,k,p)  = no_data !Actual evapotranspiration        
                  wsoil_pft(i,j,k,p)  = no_data !Soil moisture (mm)
                  rml_pft(i,j,k,p)    = no_data !MAINTENANCE ARESP
                  rmf_pft(i,j,k,p)    = no_data
                  rms_pft(i,j,k,p)    = no_data
                  rm_pft(i,j,k,p)     = no_data
                  rgl_pft(i,j,k,p)    = no_data !GROWTH RESPIRATION 
                  rgf_pft(i,j,k,p)    = no_data
                  rgs_pft(i,j,k,p)    = no_data
                  rg_pft(i,j,k,p)     = no_data
                  betal(i,j,k,p) = no_data ! to store accumulated beta_leaf(kgC/m2)
                  betaw(i,j,k,p) = no_data ! 
                  betaf(i,j,k,p) = no_data

               enddo              
            enddo
     
!     Only for land grid points
!     -------------------------
            if (nint(lsmk(i,j)).ne.0) then
               do k=1,12
                  wg0(i,j,k) = -1.0 
               enddo

!     Set some variables
!     ------------------
               do p = 1,q
                  wini(p)  = 0.01  !Soil moisture_initial condition (mm)
                  gini(p)  = 0.0   !Soil ice_initial condition (mm)
                  sini(p)  = 0.0   !Overland snow_initial condition (mm)
               enddo

!     =================
!     START INTEGRATION
!     =================
               n = 0

 10            continue

               n = n + 1

               k = mod(n,12)
               if (k.eq.0) k = 12
               mes = k
               spre = p0(i,j,k)
               td = tsoil(i,j,k)
               ta = temp(i,j,k)
               pr = prec(i,j,k)
               ipar = par(i,j,k)
               ru = rhs(i,j,k)
               ae = 2.895*ta+52.326 !Available energy (W/m2) - From NCEP-NCAR reanalysis data
!     
!     Monthly water budget
!     ====================

               call budget (mes,wini,gini,sini,td,ta,pr,spre,ae
     &             ,ca,ipar,ru,cleaf1_pft,cawood1_pft,cfroot1_pft
     &             ,wfim,gfim, sfim,smes,rmes,emes,epmes,phmes,armes
     &             ,nppmes,laimes,clmes,csmes,hrmes,rcmes,rmlmes,rmfmes
     &             ,rmsmes,rmmes,rglmes,rgfmes,rgsmes,rgmes,cleafmes
     &             ,cawoodmes,cfrootmes, gridocpmes,betalmes, betawmes
     &             ,betafmes)

               do p=1,q
                  if(p .eq. 1) emaxm(i,j,k) = epmes
                  gsoil    (i,j,k,p) = gfim(p)
                  ssoil    (i,j,k,p) = sfim(p)
                  wsoil_pft(i,j,k,p) = wfim(p)
                  snowm    (i,j,k,p) = smes(p)
                  runom_pft(i,j,k,p) = rmes(p)
                  evapm_pft(i,j,k,p) = emes(p)
                  rcm_pft  (i,j,k,p) = rcmes(p)
                  lai_pft  (i,j,k,p) = laimes(p)
                  photo_pft(i,j,k,p) = phmes(p)
                  aresp_pft(i,j,k,p) = armes(p)
                  npp_pft  (i,j,k,p) = nppmes(p)
                  clit_pft (i,j,k,p) = clmes(p)
                  csoil_pft(i,j,k,p) = csmes(p)
                  hresp_pft(i,j,k,p) = hrmes(p)
                  rml_pft  (i,j,k,p) = rmlmes(p)
                  rmf_pft  (i,j,k,p) = rmfmes(p)
                  rms_pft  (i,j,k,p) = rmsmes(p)
                  rm_pft   (i,j,k,p) = rmmes(p)
                  rgl_pft   (i,j,k,p) = rglmes(p)
                  rgf_pft   (i,j,k,p) = rgfmes(p)
                  rgs_pft   (i,j,k,p) = rgsmes(p)
                  rg_pft    (i,j,k,p) = rgmes(p)
                  
                  wini(p) = wfim(p)
                  gini(p) = gfim(p)
                  sini(p) = sfim(p)
                  
                  cleaf1_pft(p)  = cleafmes(p) 
                  cawood1_pft(p) = cawoodmes(p)
                  cfroot1_pft(p) = cfrootmes(p)

                  betal(i,j,k,p) = betalmes(p)
                  betaw(i,j,k,p) = betawmes(p)
                  betaf(i,j,k,p) = betafmes(p)
                  
                  if(k .eq. 12) then
                      cleaf_pft (i,j,p) = cleafmes(p)
                      cawood_pft(i,j,p) = cawoodmes(p)
                      cfroot_pft(i,j,p) = cfrootmes(p)
                      grid_area(i,j,p) = gridocpmes(p)
                  endif
               enddo
               
!     Check if equilibrium is attained
!     --------------------------------
               if (k.eq.12) then
                  wsoilt(i,j,k) = 0.0
                  gsoilt(i,j,k) = 0.0

                  do p = 1,q
                     wsoilt(i,j,k) = wsoilt(i,j,k) + wsoil_pft(i,j,k,p)
                     gsoilt(i,j,k) = gsoilt(i,j,k) + gsoil(i,j,k,p)
                  enddo

                  wmax = 500.
                  nerro = 0

                  do kk=1,12
                     wsaux1 = wsoilt(i,j,kk) + gsoilt(i,j,kk)                     
                     dwww = (wsaux1 - wg0(i,j,kk)) / wmax
                     if (abs(dwww).gt.0.01) nerro = nerro + 1
                  enddo

                  if (nerro.ne.0) then

                     do kk=1,12
                        wg0(i,j,kk) = wsoilt(i,j,kk) + gsoilt(i,j,kk)
                     enddo
                  else
                     goto 100
                  endif
               endif
               goto 10               
 100           continue
 

!     PFTs equilibrium check
!     ==================================================
!     tentativa 1 - usando variacao no pool de C vegetal
!     --------------------------------------------------
               if (k.eq.12) then
                  nerro = 0
                  biomass = 0.0
                  biomass0 = 0.0
                  check = .false.
                  sensi = 0.5 ! (kg/m2/y) if biomas change .le. sensi: equilibrium
                  call pft_par(4,wood)

                  do p = 1,q
                     if(cleaf_pft(i,j,p) .gt. 0.0 .and.
     &                   cfroot_pft(i,j,p).gt. 0.0) then
                        ! GRASS
                        if(wood(p) .le. 0.0) then
                           check = .true.
                           biomass = biomass + cleaf1_pft(p) +
     &                         cfroot1_pft(p)  
                           biomass0 = biomass0 + leaf0(p) +
     &                         froot0(p)
                        ! WOOD
                        else
                           check = .true.
                           biomass = biomass + cleaf1_pft(p) +
     &                         cfroot1_pft(p) + cawood1_pft(p)
                           biomass0 = biomass0 + leaf0(p) +
     &                         froot0(p) + awood0(p)
                        endif
                     endif
                  enddo
                  if(check) then
                     if (abs(biomass0-biomass) .gt. sensi) then
                        do p=1,q
                           leaf0(p) = cleaf1_pft(p)
                           froot0(p) = cfroot1_pft(p)
                           awood0(p) = cawood1_pft(p)
                           nerro = nerro + 1
                        enddo
                     endif
                  endif
                  if(nerro .gt. 0) then
!                     print*, abs(biomass-biomass0), n
                     goto 10
                  else
 !                    print*, 'eq attained',n
                     continue
                  endif
               endif
            endif               ! endif lsmk
c     finalize ny loop
         enddo                  ! j
c    finalize nx loop
      enddo                     ! I
      return
            
      end subroutine wbm

      

      subroutine budget (month,w1,g1,s1,ts,temp,prec,p0,ae,ca
     $    ,ipar,rh,cl1_pft,ca1_pft,cf1_pft,w2,g2,s2,smavg,ruavg,evavg
     &    ,epavg,phavg,aravg,nppavg,laiavg,clavg,csavg,hravg,rcavg
     &    ,rmlavg,rmfavg,rmsavg,rmavg,rglavg,rgfavg,rgsavg,rgavg
     &    ,cleafavg_pft,cawoodavg_pft,cfrootavg_pft,ocpavg,betalavg
     &    ,betawavg,betafavg)
      implicit none

      integer, parameter :: npft = 6 

!     ----------------------------INPUTS-------------------------------
!        
      integer month             !Actual month (1-12)
      real w1(npft)             !Initial (previous month last day) soil moisture storage (mm)
      real g1(npft)             !Initial soil ice storage (mm)
      real s1(npft)             !Initial overland snow storage (mm)
      real cl1_pft(npft)        ! initial BIOMASS cleaf compartment
      real cf1_pft(npft)        !                 froot
      real ca1_pft(npft)        !                 cawood
      real ts                   !Soil temperature (oC)
      real temp                 !Surface air temperature (oC)
      real prec                 !Precipitation (mm/day)
      real p0                   !Surface pressure (mb)
      real ae                   !Available energy (W/m2)
      real ca                   !Atmospheric carbon
      real ipar                 !Incident photosynthetic active radiation
      real rh                   !Relative humidity

!     ----------------------------OUTPUTS------------------------------
      real w2(npft)             !Final (last day) soil moisture storage (mm)
      real g2(npft)             !Final soil ice storage (mm)
      real s2(npft)             !Final overland snow storage (mm)
      real smavg(npft)          !Snowmelt monthly average (mm/day)
      real ruavg(npft)          !Runoff monthly average (mm/day)
      real evavg(npft)          !Actual evapotranspiration monthly average (mm/day)
      real epavg                !Maximum evapotranspiration monthly average (mm/day)
      real phavg(npft)          !Monthly photosynthesis
      real aravg(npft)          !Monthly autotrophic respiration
      real nppavg(npft)         !Monthly NPP (average between PFTs)
      real laiavg(npft)         !Monthly leaf area Index
      real clavg(npft)          !Monthly carbon litter
      real csavg(npft)          !Monthly carbon soil
      real hravg(npft)          !Monthly heterotrophic respiration
      real rcavg(npft)          !Monthly canopy resistence
      real rmlavg(npft),rmfavg(npft) ! maintenance/growth respiration
      real rmsavg(npft),rmavg(npft)  
      real rglavg(npft),rgfavg(npft)
      real rgsavg(npft),rgavg(npft)
      real cleafavg_pft(npft) ! Carbon in plant tissues
      real cawoodavg_pft(npft) 
      real cfrootavg_pft(npft)
      real ocpavg(npft), betalavg(npft), betawavg(npft), betafavg(npft)
      
!     -----------------------Internal Variables------------------------
      integer p,i
      
      real alfa_leaf(npft), alfa_awood(npft), alfa_froot(npft)
      real beta_leaf(npft), beta_awood(npft), beta_froot(npft)

!     RELATED WITH GRIDCELL OCUPATION
      
      REAL OCP_COEFFS(NPFT), ocp_mm(npft)
      LOGICAL OCP_WOOD(NPFT)
      
!     WBM COMMUNICATION (water balance)
      real wmax                 !Soil moisture availability (mm)
      real tsnow                !Temperature threshold for snowfall (oC)
      real tice                 !Temperature threshold for soil freezing (oC)
      real psnow                !Snowfall (mm/day)
      real prain                !Rainfall (mm/day)
      real rimelt(npft)               !Runoff due to soil ice melting
      real smelt(npft)                !Snowmelt (mm/day)
      real w   (npft)           !Daily soil moisture storage (mm)
      real g   (npft)           !Daily soil ice storage (mm)
      real s   (npft)           !Daily overland snow storage (mm)
      real ds  (npft)
      real dw  (npft)
      real roff(npft)           !Total runoff
      real evap(npft)           !Actual evapotranspiration (mm/day)
      real emax
      real*8 wapft                     !Maximum evapotranspiration
      
      integer ndmonth(12)       !Number of months
      data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/ !Number of days for each month
      
c     Carbon Cycle
      real ph  (npft)           !Canopy gross photosynthesis (kgC/m2/yr)
      real ar  (npft)           !Autotrophic respiration (kgC/m2/yr)
      real nppa(npft)           !Net primary productivity / auxiliar
      real laia(npft)           !Leaf area index (m2 leaf/m2 area)
      real cl  (npft)           !Litter carbon (kgC/m2)
      real cs  (npft)           !Soil carbon (kgC/m2) 
      real hr  (npft)           !Heterotrophic (microbial) respiration (kgC/m2/yr)
      real rc2 (npft)           !Canopy resistence (s/m)
      real f1  (npft)           !
      real f5  (npft)           !Photosynthesis (mol/m2/s)
c      real f1b (npft)           !Photosynthesis (micromol/m2/s)
      real vpd (npft)           !Vapor Pressure deficit
      real rm(npft),rml(npft),rmf(npft),rms(npft) ! maintenance & growth a.resp
      real rg(npft),rgl(npft),rgf(npft),rgs(npft)
      real cl1(npft),cf1(npft),ca1(npft) ! carbon pre-allocation 
      real cl2(npft),cf2(npft),ca2(npft) ! carbon pos-allocation
      

!     Initialize Canopy Resistence Parameters
!     ---------------------------------------
      do p = 1,npft
         rc2(p) = 0.0
         f1(p)  = 0.0
c         f1b(p) = 0.0
         f5(p) = 0.0
      enddo
      
!     Parameters
!     ----------
c      rh    = 0.685
      wmax  = 500.0
      tsnow = -1.0
      tice  = -2.5

!     Precipitation
!     =============     
      psnow = 0.0
      prain = 0.0
      if (temp.lt.tsnow) then
         psnow = prec/real(ndmonth(month))
      else
         prain = prec/real(ndmonth(month))
      endif
      
!     Initialization
!     --------------
      epavg = 0.0
      do p = 1,npft
         w(p)       = w1(p)     ! hidrological pools state vars  
         g(p)       = g1(p)
         s(p)       = s1(p)
         smavg(p)   = 0.0       ! accumulators for i days of month m for pft p
         ruavg(p)   = 0.0
         evavg(p)   = 0.0
         rcavg(p)   = 0.0
         laiavg(p)  = 0.0
         phavg(p)   = 0.0
         aravg(p)   = 0.0
         nppavg(p)  = 0.0
         clavg(p)   = 0.0
         csavg(p)   = 0.0
         hravg(p)   = 0.0
         rmlavg(p)  = 0.0
         rmfavg(p)  = 0.0
         rmsavg(p)  = 0.0
         rmavg(p)   = 0.0
         rglavg(p)  = 0.0
         rgfavg(p)  = 0.0
         rgsavg(p)  = 0.0
         rgavg(p)   = 0.0
         ocpavg(p)  = 0.0
         ocp_mm(p)  = 0.0
         betalavg(p) = 0.0
         betawavg(p) = 0.0
         betafavg(p) = 0.0
         alfa_leaf(p) = 0.0
         alfa_awood(p) = 0.0
         alfa_froot(p) = 0.0

         
      enddo
      
!     Numerical integration
!     ---------------------
      do i=1,ndmonth(month)
         emax  = 0.0
         do p=1,npft
            cl1(p) = cl1_pft(p)
            ca1(p) = ca1_pft(p)
            cf1(p) = cf1_pft(p)
            
            beta_leaf(p) = alfa_leaf(p)
            beta_awood(p) = alfa_awood(p)
            beta_froot(p) = alfa_froot(p)

            nppa(p)  = 0.0
            ph(p)    = 0.0
            ar(p)    = 0.0
            laia(p)  = 0.0
            f5(p)    = 0.0
            f1(p)    = 0.0
            vpd(p)   = 0.0
            rc2(p)   = 0.0
            rm(p)    = 0.0
            rml(p)   = 0.0
            rmf(p)   = 0.0
            rms(p)   = 0.0
            rg(p)    = 0.0
            rgl(p)   = 0.0
            rgf(p)   = 0.0
            rgs(p)   = 0.0

            if ((i.eq.1).and.(month.eq.1)) then    
               beta_leaf(p) = 0.00000001
               beta_awood(p) = 0.00000001
               beta_froot(p)= 0.00000001
            endif
         enddo
         
!     Grid cell area fraction (%) ocp_coeffs(pft(1), pft(2), ...,pft(p))
!     =================================================================     
         CALL PFT_AREA_FRAC(CL1, CF1, CA1, OCP_COEFFS, OCP_WOOD) ! def in productivity1.f
            
!     Maximum evapotranspiration   (emax)
!     ================================= 
         call evpot2 (p0,temp,rh,ae,emax)
             
!     Productivity (ph, aresp, vpd, rc2 & etc.) for each PFT
!     ================================= 
         do p = 1,npft
            
            call productivity1 (p,ocp_coeffs(p),OCP_WOOD(P),temp,p0,w(p)
     &          ,wmax,ca,ipar,rh,cl1(p),ca1(p),cf1(p),beta_leaf(p)
     &          ,beta_awood(p),beta_froot(p),emax,ph(p),ar(p)
     &          ,nppa(p),laia(p),f5(p),f1(p),vpd(p),rm(p),rml(p)
     &          ,rmf(p),rms(p),rg(p),rgl(p),rgf(p),rgs(p),rc2(p))

            
c     Carbon allocation (carbon content on each compartment)
!     =====================================================
            call allocation (p, nppa(p), cl1(p), ca1(p), !output !input
     &          cf1(p),cl2(p), ca2(p), cf2(p)) 

            
            alfa_leaf(p)  = amax1((cl2(p) - cl1(p)), 0.0) 
            alfa_awood(p) = amax1((ca2(p) - ca1(p)), 0.0)
            alfa_froot(p) = amax1((cf2(p) - cf1(p)), 0.0)
            
!     Snow budget
!     ===========     
            smelt(p) = 2.63 + 2.55*temp + 0.0912*temp*prain !Snowmelt (mm/day)
            smelt(p) = amax1(smelt(p),0.)
            smelt(p) = amin1(smelt(p),s(p)+psnow)
            ds(p) = psnow - smelt(p)
            s(p) = s(p) + ds(p)
            
!     Water budget
!     ============
            if (ts.le.tice) then !Frozen soil
               g(p) = g(p) + w(p) !Soil moisture freezes
               w(p) = 0.0
               roff(p) = smelt(p) + prain !mm/day
               evap(p) = 0.0
               ph(p) = 0.0
               ar(p) = 0.0
               nppa(p) = 0.0
               laia(p) = 0.0
               cl(p) = 0.0
               cs(p) = 0.0
               hr(p) = 0.0
            else                !Non-frozen soil
               w(p) = w(p) + g(p)
               g(p) = 0.0
               rimelt(p) = 0.0
               if (w(p).gt.wmax) then
                  rimelt(p) = w(p) - wmax !Runoff due to soil ice melting
                  w(p) = wmax
               endif
     
               wapft = (w(p)/wmax)
               call runoff (wapft,roff(p))       !Soil moisture runoff (roff, mm/day)

               call penman (p0,temp,rh,ae,rc2(p),evap(p)) !Actual evapotranspiration (evap, mm/day)
               dw(p) = prain + smelt(p) - evap(p) - roff(p)
               w(p) = w(p) + dw(p)
               if (w(p).gt.wmax) then
                  roff(p) = roff(p) + (w(p) - wmax)
                  w(p) = wmax
               endif
               if (w(p).lt.0.) w(p) = 0.
               roff(p) = roff(p) + rimelt(p) !Total runoff

               
!     Carbon cycle (Microbial respiration, litter and soil carbon)
!     ============================================================     
               call carbon2 (ts,f5(p),evap(p),laia(p)
     &             ,cl(p),cs(p),hr(p))   
            endif

!     Accumulate daily budgets weighted by occupation coefficients

            ocp_mm(p) = ocp_mm(p) + ocp_coeffs(p)
            
            if(p .eq. 1) epavg = epavg + emax !mm/day
            smavg(p) = smavg(p) + smelt(p)
            ruavg(p) = ruavg(p) + roff(p)  ! mm day-1
            evavg(p) = evavg(p) + evap(p) * ocp_coeffs(p)  ! mm day-1
            rcavg(p) = rcavg(p) + rc2(p) * ocp_coeffs(p) ! s m -1
            
            phavg(p) = phavg(p) +   ph(p) !kgC/m2/day
            aravg(p) = aravg(p) +   ar(p)  !kgC/m2/year
            nppavg(p) = nppavg(p) + nppa(p) !kgC/m2/day
            
            laiavg(p) = laiavg(p) + laia(p)
            clavg(p) = clavg(p) + cl(p) * ocp_coeffs(p) !kgC/m2/day
            csavg(p) = csavg(p) + cs(p) * ocp_coeffs(p) !kgC/m2/day
            hravg(p) = hravg(p) + hr(p) * ocp_coeffs(p) !kgC/m2/day
            rmlavg(p) = rmlavg(p) + rml(p) * ocp_coeffs(p)
            rmfavg(p) = rmfavg(p) + rmf(p) * ocp_coeffs(p)
            rmsavg(p) = rmsavg(p) + rms(p) * ocp_coeffs(p)
            rmavg(p) = rmavg(p) + rm(p) * ocp_coeffs(p)
            rglavg(p) = rglavg(p) + rgl(p) * ocp_coeffs(p)
            rgfavg(p) = rgfavg(p) + rgf(p) * ocp_coeffs(p)
            rgsavg(p) = rgsavg(p) + rgs(p) * ocp_coeffs(p)
            rgavg(p) = rgavg(p) + rg(p) * ocp_coeffs(p)
            cleafavg_pft(p)  =  cl2(p)
            cawoodavg_pft(p) =  ca2(p)
            cfrootavg_pft(p) =  cf2(p)
            betalavg(p) = betalavg(p) + alfa_leaf(p)  
            betawavg(p) = betawavg(p) + alfa_awood(p)
            betafavg(p) = betafavg(p) + alfa_froot(p)
            cl1_pft(p) = cl2(p)     ! Adicionado ------ para fazer transforcaoes diarias
            ca1_pft(p) = ca2(p)     ! Adicionado ------ para fazer transforcaoes diarias
            cf1_pft(p) = cf2(p)     ! Adicionado ------ para fazer transforcaoes diarias
         enddo
      enddo                     ! end ndmonth loop
      
!     Final calculations
!     ------------------
!     monthly values
      do p=1,NPFT
         if (p .eq. 1) epavg = epavg/real(ndmonth(month))
         w2(p) = w(p) * (ocp_mm(p)/real(ndmonth(month)))
         g2(p) = g(p) * (ocp_mm(p)/real(ndmonth(month)))
         s2(p) = s(p) * (ocp_mm(p)/real(ndmonth(month)))
         smavg(p) = smavg(p)/real(ndmonth(month))
         ruavg(p) = ruavg(p)/real(ndmonth(month))
         evavg(p) = evavg(p)/real(ndmonth(month))
         rcavg(p) = rcavg(p)/real(ndmonth(month))
         phavg(p) = (phavg(p)/365.0) * 12.0 !kgC/m2/yr
         aravg(p) = (aravg(p)/365.0) * 12.0   !kgC/m2/yr
         nppavg(p) = (nppavg(p)/365.0) * 12.0 !kgC/m2/yr
         laiavg(p) = (laiavg(p)/365.0) * 12.0
         clavg(p) = (clavg(p)/365.0) * 12.0 !kgC/m2
         csavg(p) = (csavg(p)/365.0) * 12.0 !kgC/m2
         hravg(p) = (hravg(p)/365.0) * 12.0 !kgC/m2/yr
         rmlavg(p) = (rmlavg(p)/365.0) * 12.0 
         rmfavg(p) = (rmfavg(p)/365.0) * 12.0
         rmsavg(p) = (rmsavg(p)/365.0) * 12.0
         rmavg(p) = (rmavg(p)/365.0) * 12.0 
         rglavg(p) = (rglavg(p)/365.0) * 12.0 
         rgfavg(p) = (rgfavg(p)/365.0) * 12.0 
         rgsavg(p) = (rgsavg(p)/365.0) * 12.0 
         rgavg(p) = (rgavg(p)/365.0) * 12.0
         betalavg(p) = (betalavg(p)/365.0) * 12.0   
         betawavg(p) = (betawavg(p)/365.0) * 12.0 
         betafavg(p) = (betafavg(p)/365.0) * 12.0 
         ocpavg(p) = ocp_coeffs(p) * 100.
            
      enddo
      return
      end subroutine budget
