module class_geometry

  implicit none

  private

  type, public :: geometry
     integer :: nsp, nst
     real*8, dimension(:),   allocatable :: rm, psim, qm, gm, rim
     real*8, dimension(:,:), allocatable :: xm, zm, pom, bm, bripm, gim
   contains
     procedure, pass :: init
     procedure, pass :: ld_data
     procedure, pass :: wrt1d
     procedure, pass :: wrt2d
     procedure, pass :: slct1d
     procedure, pass :: slct2d
     final           :: free

     generic :: wrt => wrt1d, wrt2d
  end type geometry

contains

  subroutine init(self, nsp, nst)
    class(geometry), intent(inout) :: self
    integer, intent(in) :: nsp, nst

    self%nsp = nsp
    self%nst = nst
    allocate(self%rm(nsp))
    allocate(self%psim(nsp))
    allocate(self%qm(nsp))
    allocate(self%gm(nsp))
    allocate(self%rim(nsp))
    allocate(self%xm(nsp,nst))
    allocate(self%zm(nsp,nst))
    allocate(self%pom(nsp,nst))
    allocate(self%bm(nsp,nst))
    allocate(self%bripm(nsp,nst))
    allocate(self%gim(nsp,nst))
    return
  end subroutine init

  subroutine ld_data(self, fid)
    class(geometry), intent(inout) :: self
    integer, intent(in) :: fid
    integer :: j, k

    do k = 1, self%nsp
       read(fid, '(5E15.6)') &
            & self%rm(k), self%psim(k), self%qm(k), &
            & self%gm(k), self%rim(k) 
    end do

    do j = 1, self%nst
       do k = 1, self%nsp
          read(fid, '(6E15.6)') &
               & self%xm(k,j), self%zm(k,j), &
               & self%pom(k,j), self%bm(k,j), &
               & self%bripm(k,j), self%gim(k,j) 
       end do
    end do
    return
  end subroutine ld_data

  subroutine wrt1d(self, str1, str2)
    class(geometry), intent(in) :: self
    character(len=*), intent(in) :: str1, str2
    integer :: fid, k
    character(len=15) :: filename
    real*8, dimension(self%nsp) :: vec1, vec2

    fid = 21     
    filename = 'gm1_'//trim(str1)//'_'//trim(str2)//'.dat'

    open(fid, file=filename, status='unknown', action='write')

    call self%slct1d(str1, vec1)
    call self%slct1d(str2, vec2)

    write(fid, '(2E15.6)') (vec1(k), vec2(k), k = 1, self%nsp)

    close(fid)
    return
  end subroutine wrt1d

  subroutine wrt2d(self, str)
    class(geometry), intent(in) :: self
    character(len=*), intent(in) :: str
    integer :: fid, k, j
    character(len=15) :: filename
    real*8, dimension(self%nsp, self%nst) :: vec

    fid = 21
    filename = 'gm2_'//trim(str)//'.dat'

    open(fid, file=filename, status='unknown', action='write')

    call self%slct2d(str, vec)

    do j = 1, self%nst
       write(fid, '(3E15.6)') (self%xm(k,j), self%zm(k,j), &
            & vec(k,j), k = 1, self%nsp)
       write(fid, '(a)')
    end do

    close(fid)
    return
  end subroutine wrt2d

  subroutine free(self)
    type(geometry), intent(inout) :: self

    deallocate(self%rm)
    deallocate(self%psim)
    deallocate(self%qm)
    deallocate(self%gm)
    deallocate(self%rim)
    deallocate(self%xm)
    deallocate(self%zm)
    deallocate(self%pom)
    deallocate(self%bm)
    deallocate(self%bripm)
    deallocate(self%gim)
    return
  end subroutine free

  subroutine slct1d(self, str, vec)
    class(geometry), intent(in) :: self
    character(len=*), intent(in) :: str
    real*8, dimension(:), intent(out) :: vec

    select case(trim(str))
    case('r')
       vec = self%rm
    case('psip')
       vec = self%psim
    case('q')
       vec = self%qm
    case('g')
       vec = self%gm
    case('I')
       vec = self%rim
    case default
       stop 'undefined 1d array !!!'
    end select
    return
  end subroutine slct1d

  subroutine slct2d(self, str, vec)
    class(geometry), intent(in) :: self
    character(len=*), intent(in) :: str
    real*8, dimension(:,:), intent(out) :: vec

    select case(trim(str))
    case('psip')
       vec = self%pom
    case('b')
       vec = self%bm
    case('rip')
       vec = self%bripm
    case('giac')
       vec = self%gim
    case default
       stop 'undefined 2d array !!!'
    end select
    return
  end subroutine slct2d

end module class_geometry
