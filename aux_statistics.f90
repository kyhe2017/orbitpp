! --- The derived class introduced to improve statistics of distribution, by
! --- use of global exponential weight assignment.
module class_aux_statistics
   use class_statistics
   implicit none

   private

   type, extends(statistics), public :: aux_statistics
   contains
      procedure, pass :: wrt1d
      procedure, pass :: wrt2d
   end type aux_statistics

contains

   subroutine wrt1d(self, str)
      class(aux_statistics), intent(inout) :: self
      character(len=*), intent(in) :: str
      integer :: fid, k, n
      real*8, dimension(self%dist%nlive) :: vec, wegtF
      real*8, dimension(self%nd1) :: x, pow, f1d
      real*8 :: x1, x2, dx, tolpow
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

      do n = 1, self%dist%nlive
         tolpow = 0.d0
         do k = 1, self%nd1
            pow(k) = dexp(-abs(x(k)-vec(n))/dx)
            tolpow = tolpow + pow(k)
         end do
         pow = pow/tolpow

         do k = 1, self%nd1
            f1d(k) = f1d(k) + pow(k)*wegtF(n)
         end do
      end do

      write(fid, '(2E15.6)') (x(k), f1d(k), k = 1, self%nd1)

      close(fid)
      return
   end subroutine wrt1d

   subroutine wrt2d(self, str1, str2)
      class(aux_statistics), intent(inout) :: self
      character(len=*), intent(in) :: str1, str2
      integer :: fid, k, l, n
      real*8, dimension(self%dist%nlive) :: vec1, vec2, wegtF
      real*8, dimension(self%nd1,self%nd2) :: x, y, pow, f2d
      real*8 :: x1, x2, dx, y1, y2, dy, tolpow
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
         do l = 1, self%nd2
            x(k,l) = (k-1)*dx + x1
            y(k,l) = (l-1)*dy + y1
         end do
      end do
      f2d = 0.d0

      do n = 1, self%dist%nlive
         tolpow = 0.d0
         do k = 1, self%nd1
            do l = 1,self%nd2
               pow(k,l) = dexp(-abs(x(k,l)-vec1(n))/dx)*dexp(-abs(y(k,l)-vec2(n))/dy)
               tolpow = tolpow + pow(k,l)
            end do
         end do
         pow = pow/tolpow

         do k = 1, self%nd1
            do l = 1, self%nd2
               f2d(k,l) = f2d(k,l) + pow(k,l)*wegtF(n)
            end do
         end do
      end do

      do l = 1, self%nd2
         write(fid, '(3E15.6)') (x(k,l), y(k,l), f2d(k,l), k = 1, self%nd1)
         write(fid, '(a)')
      end do

      close(fid)
      return
   end subroutine wrt2d

end module class_aux_statistics
