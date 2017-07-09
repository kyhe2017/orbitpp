module class_statistics

   use class_distribution
   
   implicit none

   private

   type, public :: statistics
      integer :: nd1, nd2
      integer :: orbit_types = 9
      integer :: orbit_regimes = 14
      class(distribution), pointer :: dist => null()
   contains
      procedure, pass :: init
      procedure, pass :: wrt1d
      procedure, pass :: wrt2d
      procedure, pass :: count_opt
      procedure, pass :: count_reg
      final           :: free

      generic :: wrt => wrt1d, wrt2d
   end type statistics

contains

   subroutine init(self, dist, nd1, nd2)
      class(statistics), intent(inout) :: self
      class(distribution), target, intent(in) :: dist
      integer, intent(in) :: nd1, nd2

      self%dist => dist
      self%nd1 = nd1
      self%nd2 = nd2
      self%orbit_types   = 9
      self%orbit_regimes = 14
      return
   end subroutine init

   subroutine wrt1d(self, str)
      class(statistics), intent(inout) :: self
      character(len=*), intent(in) :: str
      integer :: fid, k, xd
      real*8, dimension(self%dist%nlive) :: vec, wegtF
      real*8, dimension(self%nd1) :: x, f1d
      real*8 :: x1, x2, dx, dx1, dx2, tmp
      character(len=20) :: filename
      character(len=8) :: numstr

      fid = 21

      associate (time => self%dist%time, ltag => self%dist%ltag)
        if (time <= 9) then
           write(numstr, '(A1,I1,A1,I1)') '_', ltag, '_', time
           filename = 'st1_'//trim(str)//numstr(1:4)//'.dat'
        else if (time <= 99) then
           write(numstr, '(A1,I1,A1,I2)') '_', ltag, '_', time
           filename = 'st1_'//trim(str)//numstr(1:5)//'.dat'
        end if
      end associate

      open(fid, file=filename, status='unknown', action='write')

      call self%dist%slctvec(str, vec)
      call self%dist%slctvec('wegtF', wegtF)

      x1 = minval(vec)
      x2 = maxval(vec)
      dx = (x2-x1)/(self%nd1 - 1)
      do k = 1, self%nd1
         x(k) = (k-1)*dx + x1
      end do
      f1d = 0.d0

      do k = 1, self%dist%nlive
         tmp = vec(k) - x1
         xd  = tmp/dx + 1
         xd  = max(xd, 1)
         xd  = min(xd, self%nd1 - 1)
         dx2 = tmp/dx - xd + 1.d0
         dx1 = 1.d0 - dx2

         f1d(xd)   = f1d(xd) + dx1*wegtF(k)
         f1d(xd+1) = f1d(xd+1) + dx2*wegtF(k)
      end do

      write(fid, '(2E15.6)') (x(k), f1d(k)/dx, k = 1, self%nd1)

      close(fid)
      return
   end subroutine wrt1d

   subroutine wrt2d(self, str1, str2)
      class(statistics), intent(inout) :: self
      character(len=*), intent(in) :: str1, str2
      integer :: fid, k, l, xd, yd
      real*8, dimension(self%dist%nlive) :: vec1, vec2, wegtF
      real*8, dimension(self%nd1) :: x
      real*8, dimension(self%nd2) :: y
      real*8, dimension(self%nd1, self%nd2) :: f2d
      real*8 :: x1, x2, dx, dx1, dx2, y1, y2, dy, dy1, dy2, tmp
      character(len=25) :: filename
      character(len=8) :: numstr

      fid = 21

      associate (time => self%dist%time, ltag => self%dist%ltag)
        if (time <= 9) then
           write(numstr, '(A1,I1,A1,I1)') '_', ltag, '_', time
           filename = 'st2_'//trim(str1)//'_'//trim(str2)//numstr(1:4)//'.dat'
        else if (time <= 99) then
           write(numstr, '(A1,I1,A1,I2)') '_', ltag, '_', time
           filename = 'st2_'//trim(str1)//'_'//trim(str2)//numstr(1:5)//'.dat'
        end if
      end associate

      open(fid, file=filename, status='unknown', action='write')

      call self%dist%slctvec(str1, vec1)
      call self%dist%slctvec(str2, vec2)
      call self%dist%slctvec('wegtF', wegtF)

      x1 = minval(vec1)
      x2 = maxval(vec1)
      dx = (x2-x1)/(self%nd1 -1)
      y1 = minval(vec2)
      y2 = maxval(vec2)
      dy = (y2-y1)/(self%nd2 -1)
      do k = 1, self%nd1
         x(k) = (k-1)*dx + x1
      end do
      do k = 1, self%nd2
         y(k) = (k-1)*dy + y1
      end do
      f2d = 0.d0

      do k = 1, self%dist%nlive
         tmp = vec1(k) - x1
         xd  = tmp/dx + 1
         xd  = max(xd, 1)
         xd  = min(xd, self%nd1 - 1)
         dx2 = tmp/dx - xd + 1.d0
         dx1 = 1.d0 - dx2

         tmp = vec2(k) - y1
         yd  = tmp/dy + 1
         yd  = max(yd, 1)
         yd  = min(yd, self%nd2 - 1)
         dy2 = tmp/dy - yd + 1.d0
         dy1 = 1.d0 - dy2

         f2d(xd,yd)     = f2d(xd,yd) + dx1*dy1*wegtF(k)
         f2d(xd,yd+1)   = f2d(xd,yd+1) + dx1*dy2*wegtF(k)
         f2d(xd+1,yd)   = f2d(xd+1,yd) + dx2*dy1*wegtF(k)
         f2d(xd+1,yd+1) = f2d(xd+1,yd+1) + dx2*dy2*wegtF(k)
      end do

      do l = 1, self%nd2
         write(fid, '(3E15.6)') (x(k), y(l), f2d(k,l)/(dx*dy), k = 1, self%nd1)
         write(fid, '(a)')
      end do

      close(fid)
      return
   end subroutine wrt2d

   subroutine count_opt(self)
      class(statistics), intent(in) :: self
      integer :: fid, k
      real*8, dimension(self%orbit_types) :: y
      real*8, dimension(self%dist%nlive) :: opt, wegtF
      character(len=25) :: filename
      character(len=8) :: numstr

      fid = 21

      associate (time => self%dist%time, ltag => self%dist%ltag)
        if (time <= 9) then
           write(numstr, '(A1,I1,A1,I1)') '_', ltag, '_', time
           filename = 'st1_orbits'//numstr(1:4)//'.dat'
        else if (time <= 99) then
           write(numstr, '(A1,I1,A1,I2)') '_', ltag, '_', time
           filename = 'st1_orbits'//numstr(1:5)//'.dat'
        end if
      end associate

      open(fid, file=filename, status='unknown', action='write')
      
      call self%dist%slctvec('opt', opt)
      call self%dist%slctvec('wegtF', wegtF)

      y = 0.d0
      do k = 1, self%dist%nlive
         if (opt(k) /= 0) then
            y(opt(k)) = y(opt(k)) + wegtF(k)
         end if
      end do

      write(fid, '(I15, E15.6)') (k, y(k), k = 1, self%orbit_types)

      close(fid)
      return
   end subroutine count_opt
   
   subroutine count_reg(self)
      class(statistics), intent(in) :: self
      integer :: fid, k
      real*8, dimension(self%orbit_regimes) :: y
      real*8, dimension(self%dist%nlive) :: reg, wegtF
      character(len=25) :: filename
      character(len=8) :: numstr

      fid = 21

      associate (time => self%dist%time, ltag => self%dist%ltag)
        if (time <= 9) then
           write(numstr, '(A1,I1,A1,I1)') '_', ltag, '_', time
           filename = 'st1_regimes'//numstr(1:4)//'.dat'
        else if (time <= 99) then
           write(numstr, '(A1,I1,A1,I2)') '_', ltag, '_', time
           filename = 'st1_regimes'//numstr(1:5)//'.dat'
        end if
      end associate

      open(fid, file=filename, status='unknown', action='write')
      
      call self%dist%slctvec('reg', reg)
      call self%dist%slctvec('wegtF', wegtF)

      y = 0.d0
      do k = 1, self%dist%nlive
         if (reg(k) /= 0) then
            y(reg(k)) = y(reg(k)) + wegtF(k)
         end if
      end do

      write(fid, '(I15, E15.6)') (k, y(k), k = 1, self%orbit_regimes)

      close(fid)
      return
   end subroutine count_reg

   subroutine free(self)
      type(statistics), intent(inout) :: self

      nullify(self%dist)
      return
   end subroutine free

end module class_statistics
