      program pls_test

      real data_pls(580)

      call pft_par(7, data_pls)

      print*, data_pls

      stop
      
      end program pls_test

      subroutine pft_par(par, dt) 
      implicit none
      integer, parameter :: vars = 580 
      integer :: par            ! parameter number 
      real, dimension(vars) :: dt
      
!     dt1 = aleaf
!     dt2 = aawood
!     dt3 = afroot
!     dt4 = tleaf
!     dt5 = tawood
!     dt6 = tfroot
!     dt7 = g1
!     dt8 = p21
!     DT9 = JMAX

      open(23,file='../inputs/pls.bin',status='old',
     &    form='unformatted',access='direct',recl=4*580)

      if(par .gt. 0 .and. par .lt. 10) then
         read(23,rec=par) dt
      else
         print*, 'search failed'
      endif
      return
      end subroutine pft_par

