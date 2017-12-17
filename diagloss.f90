module class_diagloss
! --- only for particles tracing modeling, not for steady state!!!
  implicit none

  private

  type, public :: diagloss
! --- nlive donating number of lost markers, for compatibility with distribution.o
     integer :: nrun, nlive
     integer, dimension(:), allocatable :: opt0, opt, reg0, reg
     real*8,  dimension(:), allocatable :: &
          & pol0, thet0, zet0, xx0, zz0, &
          & vpar0, vper0, en0, rmu0, tm0, ptch0, &
          & pol, thet, zet, xx, zz, &
          & vpar, vper, en, rmu, tm, ptch, &
          & wegtF, time
   contains
     procedure, pass :: init
     procedure, pass :: ld_data
     procedure, pass :: wrt1d
     procedure, pass :: wrt2d
     procedure, pass :: wrt3d
     procedure, pass :: slctvec
     final           :: free

     generic :: wrt => wrt1d, wrt2d, wrt3d
  end type diagloss

contains

  subroutine init(self, nrun, nlive)
    class(diagloss), intent(inout) :: self
    integer, intent(in) :: nrun, nlive

    self%nrun  = nrun
    self%nlive = nlive
    allocate(self%opt0(nlive), self%opt(nlive))
    allocate(self%reg0(nlive), self%reg(nlive))
    allocate(self%pol0(nlive), self%pol(nlive))
    allocate(self%thet0(nlive), self%thet(nlive))
    allocate(self%zet0(nlive), self%zet(nlive))
    allocate(self%xx0(nlive), self%xx(nlive))
    allocate(self%zz0(nlive), self%zz(nlive))
    allocate(self%vpar0(nlive), self%vpar(nlive))
    allocate(self%vper0(nlive), self%vper(nlive))
    allocate(self%en0(nlive), self%en(nlive))
    allocate(self%rmu0(nlive), self%rmu(nlive))
    allocate(self%tm0(nlive), self%tm(nlive))
    allocate(self%ptch0(nlive), self%ptch(nlive))
    allocate(self%wegtF(nlive), self%time(nlive))
    return
  end subroutine init

  subroutine ld_data(self, fid)
    class(diagloss), intent(inout) :: self
    integer, intent(in) :: fid
    integer :: k

    do k = 1, self%nlive
       read(fid, '(4I15,24E15.6)') &
            & self%opt0(k), self%opt(k), &
            & self%reg0(k), self%reg(k), &
            & self%pol0(k), self%pol(k), &
            & self%thet0(k), self%thet(k), &
            & self%zet0(k), self%zet(k), &
            & self%xx0(k), self%xx(k), &
            & self%zz0(k), self%zz(k), &
            & self%vpar0(k), self%vpar(k), &
            & self%vper0(k), self%vper(k), &
            & self%en0(k), self%en(k), &
            & self%rmu0(k), self%rmu(k), &
            & self%tm0(k), self%tm(k), &
            & self%ptch0(k), self%ptch(k), &
            & self%wegtF(k), self%time(k)
    end do
    return
  end subroutine ld_data

  subroutine wrt1d(self, str)
    class(diagloss), intent(in) :: self
    character(len=*), intent(in) :: str
    integer :: fid, k
    real*8, dimension(self%nlive) :: vec
    character(len=20) :: filename
    character(len=8) :: numstr

    fid = 21

    write(numstr, '(I2)') self%nrun
    filename = 'loss_'//trim(str)//'_'//trim(numstr)//'.dat'

    open(fid, file=filename, status='unknown', action='write')

    call self%slctvec(str, vec)

    write(fid, '(1E15.6)') (vec(k), k = 1, self%nlive)

    close(fid)
    return
  end subroutine wrt1d

  subroutine wrt2d(self, str1, str2)
    class(diagloss), intent(in) :: self
    character(len=*), intent(in) :: str1, str2
    integer :: fid, k
    real*8, dimension(self%nlive) :: vec1, vec2
    character(len=25) :: filename
    character(len=8) :: numstr

    fid = 21

    write(numstr, '(I2)') self%nrun
    filename = 'loss_'//trim(str1)//'_'//trim(str2)//'_'// &
         & trim(numstr)//'.dat'

    open(fid, file=filename, status='unknown', action='write')

    call self%slctvec(str1, vec1)
    call self%slctvec(str2, vec2)

    write(fid, '(2E15.6)') (vec1(k), vec2(k), k = 1, self%nlive)

    close(fid)
    return
  end subroutine wrt2d

  subroutine wrt3d(self, str1, str2, str3)
    class(diagloss), intent(in) :: self
    character(len=*), intent(in) :: str1, str2, str3
    integer :: fid, k
    real*8, dimension(self%nlive) :: vec1, vec2, vec3
    character(len=25) :: filename
    character(len=8) :: numstr

    fid = 21

    write(numstr, '(I2)') self%nrun
    filename = 'loss_'//trim(str1)//'_'//trim(str2)//'_'//trim(str3)// &
         & '_'//trim(numstr)//'.dat'

    open(fid, file=filename, status='unknown', action='write')

    call self%slctvec(str1, vec1)
    call self%slctvec(str2, vec2)
    call self%slctvec(str3, vec3)

    write(fid, '(3E15.6)') (vec1(k), vec2(k), vec3(k), k = 1, self%nlive)

    close(fid)
    return
  end subroutine wrt3d

  subroutine free(self)
    type(diagloss), intent(inout) :: self

    deallocate(self%opt0, self%opt)
    deallocate(self%reg0, self%reg)
    deallocate(self%pol0, self%pol)
    deallocate(self%thet0, self%thet)
    deallocate(self%zet0, self%zet)
    deallocate(self%xx0, self%xx)
    deallocate(self%zz0, self%zz)
    deallocate(self%vpar0, self%vpar)
    deallocate(self%vper0, self%vper)
    deallocate(self%en0, self%en)
    deallocate(self%rmu0, self%rmu)
    deallocate(self%tm0, self%tm)
    deallocate(self%ptch0, self%ptch)
    deallocate(self%wegtF, self%time)
    return
  end subroutine free

  subroutine slctvec(self, str, vec)
    class(diagloss), intent(in) :: self
    character(len=*), intent(in) :: str
    real*8, dimension(:), intent(out) :: vec

    select case(trim(str))
    case('opt0')
       vec = self%opt0
    case('opt')
       vec = self%opt
    case('reg0')
       vec = self%reg0
    case('reg')
       vec = self%reg
    case('pol0')
       vec = self%pol0
    case('pol')
       vec = self%pol
    case('thet0')
       vec = self%thet0
    case('thet')
       vec = self%thet
    case('zet0')
       vec = self%zet0
    case('zet')
       vec = self%zet
    case('xx0')
       vec = self%xx0
    case('xx')
       vec = self%xx
    case('zz0')
       vec = self%zz0
    case('zz')
       vec = self%zz
    case('vpar0')
       vec = self%vpar0
    case('vpar')
       vec = self%vpar
    case('vper0')
       vec = self%vper0
    case('vper')
       vec = self%vper
    case('en0')
       vec = self%en0
    case('en')
       vec = self%en
    case('rmu0')
       vec = self%rmu0
    case('rmu')
       vec = self%rmu
    case('tm0')
       vec = self%tm0
    case('tm')
       vec = self%tm
    case('ptch0')
       vec = self%ptch0
    case('ptch')
       vec = self%ptch
    case('wegtF')
       vec = self%wegtF
    case('time')
       vec = self%time
    case default
       stop 'undefined distribution variables !!!'
    end select
    return
  end subroutine slctvec

end module class_diagloss
