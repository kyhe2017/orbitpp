module class_momenta
   use orbitpp_private
   implicit none

   private

   integer :: itoggle = 0

   type, public :: momenta
      integer :: nrun
      real*8, dimension(:), allocatable :: rmpol, vol, sur, qfact, bb
      real*8, dimension(:), allocatable :: power, powermurf, powerloss
      real*8, dimension(:,:,:,:), allocatable :: pmom
   contains
      procedure, pass :: init
      procedure, pass :: ld_data
      procedure, pass :: wrt_power
      procedure, pass :: wrt2d
      procedure, pass :: wrt2d_moment
      final           :: free
      
      generic :: wrt => wrt2d, wrt2d_moment
   end type momenta

contains

   subroutine init(self, nrun)
      class(momenta), intent(inout) :: self
      integer, intent(in) :: nrun

      self%nrun = nrun
      allocate(self%rmpol(lsp))
      allocate(self%vol(lsp))
      allocate(self%sur(lsp))
      allocate(self%qfact(lsp))
      allocate(self%bb(lsp))
      allocate(self%power(0:navr-1))
      allocate(self%powermurf(0:navr-1))
      allocate(self%powerloss(0:navr-1))
      allocate(self%pmom(iks,imo,lsp,0:navr))
      return
   end subroutine init
   
   subroutine ld_data(self, fid)
      class(momenta), intent(inout) :: self
      integer, intent(in) :: fid
      integer :: l, k, n
      character(len=128) :: dum
! --- radial profiles
      read(fid, *) dum
      print*, dum
      read(fid, *) self%rmpol, self%vol, self%sur, self%qfact, self%bb

      associate (pmom => self%pmom)
! --- particle density
        read(fid, *) dum
        print*, dum
        read(fid, *) (((pmom(l,iden,k,n), k = 1, lsp), n = 0, navr), l = 1, iks)
! --- energy density
        read(fid, *) dum
        print*, dum
        read(fid, *) (((pmom(l,ieng,k,n), k = 1, lsp), n = 0, navr), l = 1, iks)
! --- power density
        read(fid, *) dum
        print*, dum
        read(fid, *) (((pmom(l,ipwr,k,n), k = 1, lsp), n = 0, navr), l = 1, iks)
! --- kick density
        read(fid, *) dum
        print*, dum
        read(fid, *) (((pmom(l,ipwrnrf,k,n), k = 1, lsp), n = 0, navr), l = 1, iks)
! --- lost power density
        read(fid, *) dum
        print*, dum
        read(fid, *) (((pmom(l,ipwrloss,k,n), k = 1, lsp), n = 0, navr), l = 1, iks)
! --- drag density
        read(fid, *) dum
        print*, dum
        read(fid, *) (((pmom(l,idrg,k,n), k = 1, lsp), n = 0, navr), l = 1, iks)
! --- radial flux
        read(fid, *) dum
        print*, dum
        read(fid, *) (((pmom(l,iflx,k,n), k = 1, lsp), n = 0, navr), l = 1, iks)
! --- density source
        read(fid, *) dum
        print*, dum
        read(fid, *) (((pmom(l,isrc,k,n), k = 1, lsp), n = 0, navr), l = 1, iks)
! --- energy source
        read(fid, *) dum
        print*, dum
        read(fid, *) (((pmom(l,isre,k,n), k = 1, lsp), n = 0, navr), l = 1, iks)
! --- parallel energy        
        read(fid, *) dum
        print*, dum
        read(fid, *) (((pmom(l,ipar,k,n), k = 1, lsp), n = 0, navr), l = 1, iks)
! --- perpendicular energy        
        read(fid, *) dum
        print*, dum
        read(fid, *) (((pmom(l,iper,k,n), k = 1, lsp), n = 0, navr), l = 1, iks)
      end associate

      call ps_data(self)
      return
   end subroutine ld_data

   subroutine ps_data(self)
      class(momenta), intent(inout) :: self
      integer :: m, n, k, l, i
      real*8, dimension(iks,lsp,0:navr) :: tmp, tqm, tqf, tqs
      real*8 :: dp, scalen, dNdt, dta

      dp = pw/(lsp-1)
      scalen = 1.d0
         
      associate (pmom => self%pmom)
! --- power deposition of fast particles

! --- time integration of EP power
        tmp = 0.d0
        do n = 0, navr-1
           do m = 0, n
              tmp(:,:,n) = tmp(:,:,n) + pmom(:,ipwr,:,m)
           end do
        end do
! --- EP power deposition rate
        do n = 0, navr-1
           dta = (n+1)*tavr
           pmom(:,ipwr,:,n) = tmp(:,:,n)/dta
        end do
! --- EP power deposition rate for all species and positions
        do n = 0, navr-1
           self%power(n) = 0.d0
           do k = 1, lsp
              do l = 1, iks
                 self%power(n) = self%power(n) + pmom(l,ipwr,k,n)
              end do
           end do
        end do
! --- power deposition of RF wave

! --- time integration of RF power
        tmp = 0.d0
        do n = 0, navr-1
           do m = 0, n
              tmp(:,:,n) = tmp(:,:,n) + pmom(:,ipwrnrf,:,m)
           end do
        end do
! --- RF power deposition rate
        do n = 0, navr-1
           dta = (n+1)*tavr
           pmom(:,ipwrnrf,:,n) = tmp(:,:,n)/dta
        end do
! --- RF power deposition rate for all species and positions
        do n = 0, navr-1
           self%powermurf(n) = 0.d0
           do k = 1, lsp
              do l = 1, iks
                 self%powermurf(n) = self%powermurf(n) + pmom(l,ipwrnrf,k,n)
              end do
           end do
        end do
! --- power losss of EPs

! --- time integration of lost EP power
        tmp = 0.d0
        do n = 0, navr-1
           do m = 0, n
              tmp(:,:,n) = tmp(:,:,n) + pmom(:,ipwrloss,:,m)
           end do
        end do
! --- EP power loss rate
        do n = 0, navr-1
           dta = (n+1)*tavr
           pmom(:,ipwrloss,:,n) = tmp(:,:,n)/dta
        end do
! --- EP power loss rate for all species and positions(only edge makes sense)
        do n = 0, navr-1
           self%powerloss(n) = 0.d0
           do k = 1, lsp
              do l = 1, iks
                 self%powerloss(n) = self%powerloss(n) + pmom(l,ipwrloss,k,n)
              end do
           end do
        end do

! --- magnetic and friction torques
        tqm = 0.d0
        tqf = 0.d0
        tqs = 0.d0
        do l = 1, iks
           do n = 0, navr 
              do k = 1, lsp 
                 do m = 1, k
                    if (itoggle == 0) then
                       tqm(l,k,n) = tqm(l,k,n) + pmom(l,iden,m,n)
                       tqs(l,k,n) = tqs(l,k,n) + pmom(l,isrc,m,n)
                    else
                       tqm(l,k,n) = tqm(l,k,n) + pmom(l,iden,m,n) - &
                            & pmom(l,iden,m,0)
                    end if
                    tqf(l,k,n) = tqf(l,k,n) + pmom(l,idrg,m,n)
                 end do
              end do
           end do
        end do

! --- just for output
        pmom(:,isrc,:,:) = pmom(:,isrc,:,:)/tavr
        pmom(:,isre,:,:) = pmom(:,isre,:,:)/tavr
        
! --- volume-integrated jxB torque
        pmom(:,iflx,:,:) = 0.
        do n = 0, navr-1 
           do k = 2, lsp 
              do l = 1, iks
                 if (itoggle == 0) then
                    dNdt = (tqm(l,k,n+1) - tqm(l,k,n))/tavr - tqs(l,k,n)/tavr
                    pmom(l,iflx,k,n) = pmom(l,iflx,k-1,n) + dNdt*dp
                    pmom(l,idrg,k,n) = tqf(l,k,n)/tavr
                 else
                    pmom(l,iflx,k,n) = pmom(l,iflx,k-1,n) + tqm(l,k,n)*dp
                 end if
              end do
           end do
        end do

! --- time-integrated collisional torque
        if (itoggle /= 0) then
           pmom(:,idrg,:,n) = 0.d0
           do n = 0, navr-1 
              do k = 1, lsp 
                 do l = 1, iks
                    do m = 1, n
                       pmom(l,idrg,k,n) = pmom(l,idrg,k,n) + tqf(l,k,m)
                    end do
                 end do
              end do
           end do
        end if

! --- normalized moments
        do n = 0, navr-1 
           do k = 1, lsp 
              do i = 1, imo
                 do l = 1, iks 
                    if (i == iflx .or. i == idrg) then
                       pmom(l,i,k,n) = scalen*pmom(l,i,k,n)
                       if(itoggle /= 0) pmom(l,i,k,n)=pmom(l,i,k,n)/(navr*tavr)
                    else
                       pmom(l,i,k,n) = scalen*(pmom(l,i,k,n)/self%vol(k))
                    end if
                 end do
              end do
           end do
        end do

! --- convert particle energy density to beta (added 07-02-02)
        do l = 1, iks
           do n = 0, navr
              do k = 1, lsp
                 pmom(l,ibet,k,n) = pmom(l,ieng,k,n)/self%bb(k)
              end do
           end do
        end do
      end associate
      return
   end subroutine ps_data

   subroutine wrt_power(self)
      class(momenta), intent(in) :: self
      integer :: k, fid
      character(len=20) :: filename
      character(len=10) :: numstr

      fid = 21

      write(numstr, '(I2)') self%nrun
      numstr = adjustl(numstr)
      filename = 'mm_POWR_'//trim(numstr)//'.dat'

      open(fid, file=filename, status='unknown', action='write')

      write(fid, '(I15,3E15.6)') (k, self%powermurf(k), self%power(k), &
           & self%powerloss(k), k = 0, navr-1)

      close(fid)
      return
   end subroutine wrt_power
   
   subroutine wrt2d(self, str)
      class(momenta), intent(in) :: self
      character(len=*), intent(in) :: str
      integer :: j, k, l, ll, fid
      character(len=20) :: filename
      character(len=10) :: numstr

      fid = 21
      ll  = idex(str)
      
      do l = 1, iks
         if (self%nrun <= 9) then
            write(numstr, '(I1,A1,I1)') self%nrun, '_', l
            filename = 'mm_'//trim(str)//'_'//numstr(1:3)//'.dat'
         else if (self%nrun <= 99) then
            write(numstr, '(I2,A1,I1)') self%nrun, '_', l
            filename = 'mm_'//trim(str)//'_'//numstr(1:4)//'.dat'
         end if

         open(fid, file=filename, status='unknown', action='write')

         do j = 0, navr-1
            write(fid, '(I15,2E15.6)') (j, self%rmpol(k), self%pmom(l,ll,k,j), &
                 & k = 1, lsp)
            write(fid, '(a)')
         end do

         close(fid)
      end do
      return
   end subroutine wrt2d
   
   subroutine wrt2d_moment(self, str, moment)
      class(momenta), intent(in) :: self
      character(len=*), intent(in) :: str
      integer, intent(in) :: moment
      integer :: k, l, ll, fid
      character(len=20) :: filename
      character(len=10) :: numstr1, numstr2, numstr3

      fid = 21
      ll  = idex(str)
      if (moment < 0 .or. moment > navr-1) stop 'moment beyond range !!!'
      
      do l = 1, iks
         write(numstr1, '(I2)') self%nrun
         write(numstr2, '(A1,I1,A1)') '_', l, '_'
         write(numstr3, '(I2)') moment
         numstr1 = adjustl(numstr1)
         numstr3 = adjustl(numstr3)
         filename = 'mm_'//trim(str)//'_'//trim(numstr1)//trim(numstr2) &
              & //trim(numstr3)//'.dat'

         open(fid, file=filename, status='unknown', action='write')

         write(fid, '(2E15.6)') (self%rmpol(k), self%pmom(l,ll,k,moment), &
                 & k = 1, lsp)

         close(fid)
      end do
      return
   end subroutine wrt2d_moment

   function idex(str)
      integer :: idex
      character(len=*), intent(in) :: str

      select case(trim(str))
      case('den')
         idex = iden
      case('eng')
         idex = ieng
      case('pwr')
         idex = ipwr
      case('drg')
         idex = idrg
      case('flx')
         idex = iflx
      case('src')
         idex = isrc
      case('sre')
         idex = isre
      case('bet')
         idex = ibet
      case('pwrnrf')
         idex = ipwrnrf
      case('pwrloss')
         idex = ipwrloss
      case('par')
         idex = ipar
      case('per')
         idex = iper
      case default
         stop 'undefined momenta attribution !!!'
      end select
      return
   end function idex

   subroutine free(self)
      type(momenta), intent(inout) :: self

      deallocate(self%rmpol)
      deallocate(self%vol)
      deallocate(self%sur)
      deallocate(self%qfact)
      deallocate(self%bb)
      deallocate(self%power)
      deallocate(self%powermurf)
      deallocate(self%powerloss)
      deallocate(self%pmom)
      return
   end subroutine free
      
end module class_momenta
