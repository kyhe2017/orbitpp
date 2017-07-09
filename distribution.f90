module class_distribution
   implicit none

   private

   type, public :: distribution
      integer :: time, ltag, nlive
      integer, dimension(:), allocatable :: nout, opt, reg
      real*8, dimension(:), allocatable :: &
           & pol, thet, zet, xx, zz, &
           & vpar, vper, vr, &
           & en, rmu, ptch, wegtF, &
           & enloss, ptchloss, zloss, &
           & philoss, xloss, timeloss
   contains
      procedure, pass :: init
      procedure, pass :: ld_data
      procedure, pass :: wrt1d
      procedure, pass :: wrt2d
      procedure, pass :: slctvec
      final           :: free

      generic :: wrt => wrt1d, wrt2d
   end type distribution

contains

   subroutine init(self, time, ltag, nlive)
      class(distribution), intent(inout) :: self
      integer, intent(in) :: time, ltag, nlive

      self%time  = time
      self%ltag  = ltag
      self%nlive = nlive
      allocate(self%nout(nlive))
      allocate(self%opt(nlive))
      allocate(self%reg(nlive))
      allocate(self%pol(nlive))
      allocate(self%thet(nlive))
      allocate(self%zet(nlive))
      allocate(self%xx(nlive))
      allocate(self%zz(nlive))
      allocate(self%vpar(nlive))
      allocate(self%vper(nlive))
      allocate(self%vr(nlive))
      allocate(self%en(nlive))
      allocate(self%rmu(nlive))
      allocate(self%ptch(nlive))
      allocate(self%wegtF(nlive))
      allocate(self%enloss(nlive))
      allocate(self%ptchloss(nlive))
      allocate(self%zloss(nlive))
      allocate(self%philoss(nlive))
      allocate(self%xloss(nlive))
      allocate(self%timeloss(nlive))
      return
   end subroutine init

   subroutine ld_data(self, fid)
      class(distribution), intent(inout) :: self
      integer, intent(in) :: fid
      integer :: k

      do k = 1, self%nlive
         read(fid, '(3I15,18E15.6)') &
              & self%nout(k), self%opt(k), self%reg(k), &
              & self%pol(k), self%thet(k), self%zet(k), &
              & self%xx(k), self%zz(k),&
              & self%vpar(k), self%vper(k), self%vr(k), &
              & self%en(k), self%rmu(k), self%ptch(k), self%wegtF(k), &
              & self%enloss(k), self%ptchloss(k), self%zloss(k), &
              & self%philoss(k), self%xloss(k), self%timeloss(k)
      end do
      return
   end subroutine ld_data

   subroutine wrt1d(self, str)
      class(distribution), intent(in) :: self
      character(len=*), intent(in) :: str
      integer :: fid, k
      real*8, dimension(self%nlive) :: vec
      character(len=20) :: filename
      character(len=8) :: numstr

      fid = 21

      if (self%time <= 9) then
         write(numstr, '(A1,I1,A1,I1)') '_', self%ltag, '_', self%time
         filename = 'ds1_'//trim(str)//numstr(1:4)//'.dat'
      else if (self%time <= 99) then
         write(numstr, '(A1,I1,A1,I2)') '_', self%ltag, '_', self%time
         filename = 'ds1_'//trim(str)//numstr(1:5)//'.dat'
      end if

      open(fid, file=filename, status='unknown', action='write')

      call self%slctvec(str, vec)

      write(fid, '(1E15.6)') (vec(k), k = 1, self%nlive)

      close(fid)
      return
   end subroutine wrt1d

   subroutine wrt2d(self, str1, str2)
      class(distribution), intent(in) :: self
      character(len=*), intent(in) :: str1, str2
      integer :: fid, k
      real*8, dimension(self%nlive) :: vec1, vec2
      character(len=25) :: filename
      character(len=8) :: numstr

      fid = 21

      if (self%time <= 9) then
         write(numstr, '(A1,I1,A1,I1)') '_', self%ltag, '_', self%time
         filename = 'ds2_'//trim(str1)//'_'//trim(str2)//numstr(1:4)//'.dat'
      else if (self%time <= 99) then
         write(numstr, '(A1,I1,A1,I2)') '_', self%ltag, '_', self%time
         filename = 'ds2_'//trim(str1)//'_'//trim(str2)//numstr(1:5)//'.dat'
      end if

      open(fid, file=filename, status='unknown', action='write')

      call self%slctvec(str1, vec1)
      call self%slctvec(str2, vec2)

      write(fid, '(2E15.6)') (vec1(k), vec2(k), k = 1, self%nlive)

      close(fid)
      return
   end subroutine wrt2d

   subroutine free(self)
      type(distribution), intent(inout) :: self

      deallocate(self%nout)
      deallocate(self%opt)
      deallocate(self%reg)
      deallocate(self%pol)
      deallocate(self%thet)
      deallocate(self%zet)
      deallocate(self%xx)
      deallocate(self%zz)
      deallocate(self%vpar)
      deallocate(self%vper)
      deallocate(self%vr)
      deallocate(self%en)
      deallocate(self%rmu)
      deallocate(self%ptch)
      deallocate(self%wegtF)
      deallocate(self%enloss)
      deallocate(self%ptchloss)
      deallocate(self%zloss)
      deallocate(self%philoss)
      deallocate(self%xloss)
      deallocate(self%timeloss)
      return
   end subroutine free

   subroutine slctvec(self, str, vec)
      class(distribution), intent(in) :: self
      character(len=*), intent(in) :: str
      real*8, dimension(:), intent(out) :: vec

      select case(trim(str))
      case('nout')
         vec = self%nout
      case('opt')
         vec = self%opt
      case('reg')
         vec = self%reg
      case('pol')
         vec = self%pol
      case('thet')
         vec = self%thet
      case('zet')
         vec = self%zet
      case('xx')
         vec = self%xx
      case('zz')
         vec = self%zz
      case('vpar')
         vec = self%vpar
      case('vper')
         vec = self%vper
      case('vr')
         vec = self%vr
      case('en')
         vec = self%en
      case('rmu')
         vec = self%rmu
      case('ptch')
         vec = self%ptch
      case('wegtF')
         vec = self%wegtF
      case('enloss')
         vec = self%enloss
      case('ptchloss')
         vec = self%ptchloss
      case('zloss')
         vec = self%zloss
      case('philoss')
         vec = self%philoss
      case('xloss')
         vec = self%xloss
      case('timeloss')
         vec = self%timeloss
      case default
         stop 'undefined distribution variables !!!'
      end select
      return
   end subroutine slctvec

end module class_distribution
