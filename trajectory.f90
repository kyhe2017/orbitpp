module class_trajectory

  implicit none

  private

  type, public :: trajectory
     integer :: nrun, itag, nlast
     real*8, dimension(:), allocatable :: &
          & tp, xp, zp, plp, thetp, rmup, rhop, torp, &
          & enp, bp, vrp, vthetp, orbp, ptchp, zetap, &
          & taurfintp, qprop, alp, brmp, bphimp, bzmp
   contains
     procedure, pass :: init
     procedure, pass :: ld_data
     procedure, pass :: wrt1d
     procedure, pass :: wrt2d
     procedure, pass :: wrt3d
     procedure, pass :: slctvec
     final           :: free

     generic :: wrt => wrt1d, wrt2d, wrt3d
  end type trajectory

contains

  subroutine init(self, nrun, itag, nlast)
    class(trajectory), intent(inout) :: self
    integer, intent(in) :: nrun, itag, nlast

    self%nrun = nrun
    self%itag = itag
    self%nlast = nlast
    allocate(self%tp(nlast))
    allocate(self%xp(nlast))
    allocate(self%zp(nlast))
    allocate(self%plp(nlast))
    allocate(self%thetp(nlast))
    allocate(self%rmup(nlast))
    allocate(self%rhop(nlast))
    allocate(self%torp(nlast))
    allocate(self%enp(nlast))
    allocate(self%bp(nlast))
    allocate(self%vrp(nlast))
    allocate(self%vthetp(nlast))
    allocate(self%orbp(nlast))
    allocate(self%ptchp(nlast))
    allocate(self%zetap(nlast))
    allocate(self%taurfintp(nlast))
    allocate(self%qprop(nlast))
    allocate(self%alp(nlast))
    allocate(self%brmp(nlast))
    allocate(self%bphimp(nlast))
    allocate(self%bzmp(nlast))
    return
  end subroutine init

  subroutine ld_data(self, fid)
    class(trajectory), intent(inout) :: self
    integer, intent(in) :: fid
    integer :: i

    do i = 1, self%nlast
       read(fid, '(21E15.6)') &
            & self%tp(i), self%xp(i), self%zp(i), self%plp(i), &
            & self%thetp(i), self%rmup(i), self%rhop(i), &
            & self%torp(i), self%enp(i), self%bp(i), self%vrp(i), &
            & self%vthetp(i), self%orbp(i), self%ptchp(i), &
            & self%zetap(i), self%taurfintp(i), self%qprop(i), &
            & self%alp(i), self%brmp(i), self%bphimp(i), &
            & self%bzmp(i)
    end do
    return
  end subroutine ld_data

  subroutine wrt1d(self, str)
    class(trajectory), intent(in) :: self
    character(len=*), intent(in) :: str
    integer :: fid, k
    real*8, dimension(self%nlast) :: vec
    character(len=15) :: filename
    character(len=6)  :: numstr

    fid = 21

    if (self%itag <= 9) then
       write(numstr, '(I1,A1,I1)') self%nrun, '_', self%itag
       filename = numstr(1:3)//trim(str)//'.dat'
    else if (self%itag <= 99) then
       write(numstr, '(I1,A1,I2)') self%nrun, '_', self%itag
       filename = numstr(1:4)//trim(str)//'.dat'
    else if (self%itag <= 999) then
       write(numstr, '(I1,A1,I3)') self%nrun, '_', self%itag
       filename = numstr(1:5)//trim(str)//'.dat'
    end if

    open(fid, file=filename, status='unknown', action='write')

    call self%slctvec(str, vec)

    write(fid, '(1E15.6)') (vec(k), k = 1, self%nlast)

    close(fid)
    return
  end subroutine wrt1d

  subroutine wrt2d(self, str1, str2)
    class(trajectory), intent(in) :: self
    character(len=*), intent(in) :: str1, str2
    integer :: fid, k
    real*8, dimension(self%nlast) :: vec1, vec2
    character(len=15) :: filename
    character(len=6)  :: numstr

    fid = 21

    if (self%itag <= 9) then
       write(numstr, '(I1,A1,I1)') self%nrun, '_', self%itag
       filename = numstr(1:3)//trim(str1)//'_'//trim(str2)//'.dat'
    else if (self%itag <= 99) then
       write(numstr, '(I1,A1,I2)') self%nrun, '_', self%itag
       filename = numstr(1:4)//trim(str1)//'_'//trim(str2)//'.dat'
    else if (self%itag <= 999) then
       write(numstr, '(I1,A1,I3)') self%nrun, '_', self%itag
       filename = numstr(1:5)//trim(str1)//'_'//trim(str2)//'.dat'
    end if

    open(fid, file=filename, status='unknown', action='write')

    call self%slctvec(str1, vec1)
    call self%slctvec(str2, vec2)

    write(fid, '(2E15.6)') (vec1(k), vec2(k), k = 1, self%nlast)

    close(fid)
    return
  end subroutine wrt2d

! --- definition if necessary
  subroutine wrt3d(self, str1, str2, str3)
    class(trajectory), intent(in) :: self
    character(len=*), intent(in) :: str1, str2, str3
    integer :: fid, k
    real*8, dimension(self%nlast) :: vec1, vec2, vec3
    character(len=15) :: filename
    character(len=6)  :: numstr

    return
  end subroutine wrt3d

  subroutine free(self)
    type(trajectory), intent(inout) :: self

    deallocate(self%tp)
    deallocate(self%xp)
    deallocate(self%zp)
    deallocate(self%plp)
    deallocate(self%thetp)
    deallocate(self%rmup)
    deallocate(self%rhop)
    deallocate(self%torp)
    deallocate(self%enp)
    deallocate(self%bp)
    deallocate(self%vrp)
    deallocate(self%vthetp)
    deallocate(self%orbp)
    deallocate(self%ptchp)
    deallocate(self%zetap)
    deallocate(self%taurfintp)
    deallocate(self%qprop)
    deallocate(self%alp)
    deallocate(self%brmp)
    deallocate(self%bphimp)
    deallocate(self%bzmp)
    return
  end subroutine free

  subroutine slctvec(self, str, vec)
    class(trajectory), intent(in) :: self
    character(len=*), intent(in) :: str
    real*8, dimension(:), intent(out) :: vec

    select case(trim(str))
    case('t')
       vec = self%tp
    case('x')
       vec = self%xp
    case('z')
       vec = self%zp
    case('pol')
       vec = self%plp
    case('thet')
       vec = self%thetp
    case('mu')
       vec = self%rmup
    case('rho')
       vec = self%rhop
    case('tor')
       vec = self%torp
    case('en')
       vec = self%enp
    case('b')
       vec = self%bp
    case('vr')
       vec = self%vrp
    case('vthet')
       vec = self%vthetp
    case('ptch')
       vec = self%ptchp
    case('zeta')
       vec = self%zetap
    case('q')
       vec = self%qprop
    case default
       stop 'undefined trajectory variables !!!'
    end select
    return
  end subroutine slctvec

end module class_trajectory
