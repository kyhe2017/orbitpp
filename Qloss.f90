module class_Qloss

  use class_diagloss

  implicit none

  private
  public :: Qloss

  type :: Qloss
     integer :: nd1, nd2
     integer :: orbit_types = 9
     class(diagloss), pointer :: dist => null()
   contains
     procedure, pass :: init
     procedure, pass :: wrt1d
     procedure, pass :: wrt2d
     procedure, pass :: wrt3d
     procedure, pass :: count_opt
     final           :: free

     generic :: wrt => wrt1d, wrt2d, wrt3d
  end type Qloss

contains

  subroutine init(self, dist, nd1, nd2)
    class(Qloss), intent(inout) :: self
    class(diagloss), target :: dist
    integer, intent(in) :: nd1, nd2

    self%dist => dist
    self%nd1 = nd1
    self%nd2 = nd2
    self%orbit_types = 9
  end subroutine init

  subroutine wrt1d(self, str)
    class(Qloss), intent(inout) :: self
    character(len=*), intent(in) :: str
    integer :: fid, k, xd
    real*8, dimension(self%dist%nlive) :: vec, wegtF
    real*8, dimension(self%nd1) :: x, f1d
    real*8 :: x1, x2, dx, dx1, dx2, tmp
    character(len=20) :: filename
    character(len=8) :: numstr

    fid = 21

    write(numstr, '(I2)') self%dist%nrun
    numstr = adjustl(numstr)
    filename = 'loss_1d_'//trim(str)//'_'//trim(numstr)//'.dat'

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

    write(fid, '(2E15.6)') (x(k), f1d(k), k = 1, self%nd1)

    close(fid)

  end subroutine wrt1d

  subroutine wrt2d(self, str1, str2)
    class(Qloss), intent(inout) :: self
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

    write(numstr, '(I2)') self%dist%nrun
    numstr = adjustl(numstr)
    filename = 'loss_2d_'//trim(str1)//'_'//trim(str2)//'_'// &
         & trim(numstr)//'.dat'

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
       write(fid, '(3E15.6)') (x(k), y(l), f2d(k,l), k = 1, self%nd1)
       write(fid, '(a)')
    end do

    close(fid)

  end subroutine wrt2d

  subroutine wrt3d(self, str1, str2, str3)
    class(Qloss), intent(inout) :: self
    character(len=*), intent(in) :: str1, str2, str3
    integer :: fid, k, l, xd, yd
    real*8, dimension(self%dist%nlive) :: vec1, vec2, vec3, wegtF
    real*8, dimension(self%nd1) :: x
    real*8, dimension(self%nd2) :: y
    real*8, dimension(self%nd1, self%nd2) :: f2d
    real*8 :: x1, x2, dx, dx1, dx2, y1, y2, dy, dy1, dy2, tmp
    character(len=25) :: filename
    character(len=8) :: numstr

    fid = 21

    write(numstr, '(I2)') self%dist%nrun
    numstr = adjustl(numstr)
    filename = 'loss_2d_'//trim(str1)//'_'//trim(str2)//'_'//trim(str3)// &
         & '_'//trim(numstr)//'.dat'

    open(fid, file=filename, status='unknown', action='write')

    call self%dist%slctvec(str1, vec1)
    call self%dist%slctvec(str2, vec2)
    call self%dist%slctvec(str3, vec3)
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

       f2d(xd,yd)     = f2d(xd,yd) + dx1*dy1*wegtF(k)*vec3(k)
       f2d(xd,yd+1)   = f2d(xd,yd+1) + dx1*dy2*wegtF(k)*vec3(k)
       f2d(xd+1,yd)   = f2d(xd+1,yd) + dx2*dy1*wegtF(k)*vec3(k)
       f2d(xd+1,yd+1) = f2d(xd+1,yd+1) + dx2*dy2*wegtF(k)*vec3(k)
    end do

    do l = 1, self%nd2
       write(fid, '(3E15.6)') (x(k), y(l), f2d(k,l), k = 1, self%nd1)
       write(fid, '(a)')
    end do

    close(fid)

  end subroutine wrt3d

  subroutine count_opt(self)
    class(Qloss), intent(in) :: self
    integer :: fid, k
    real*8, dimension(self%orbit_types) :: y
    real*8, dimension(self%dist%nlive) :: opt, wegtF
    character(len=25) :: filename
    character(len=8) :: numstr
    character(len=20), dimension(9) :: str_of_orbit

    fid = 21

    write(numstr, '(I2)') self%dist%nrun
    numstr = adjustl(numstr)
    filename = 'loss_orbits_'//trim(numstr)//'.dat'

    open(fid, file=filename, status='unknown', action='write')

    call self%dist%slctvec('opt', opt)
    call self%dist%slctvec('wegtF', wegtF)

    str_of_orbit = (/'"co-current conf"', '"co-current loss"', &
         & '"ctr-current conf"', '"ctr-current loss"', &
         & '"trapped conf"', '"trapped loss"', &
         & '"stagnation"', '"patato conf"', '"patato loss"'/)

    y = 0.d0
    do k = 1, self%dist%nlive
       if (opt(k) /= 0) then
          y(opt(k)) = y(opt(k)) + wegtF(k)
       end if
    end do

    write(fid, '(I15, E15.6, A21)') (k, y(k), str_of_orbit(k), k = 1, self%orbit_types)

    close(fid)

  end subroutine count_opt

  subroutine free(self)
    type(Qloss), intent(inout) :: self

    nullify(self%dist)

  end subroutine free

end module class_Qloss
