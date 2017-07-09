module class_csdata
   use orbitpp_private

   implicit none

   private

   type, public :: csdata
      integer :: nrun
      real*8, dimension(:,:), allocatable :: xdata, zdata
      real*8, dimension(:,:,:,:), allocatable :: power2d
   contains
      procedure, pass :: init
      procedure, pass :: ld_data
      final           :: free
   end type csdata

contains

   subroutine init(self, nrun)
      class(csdata), intent(inout) :: self
      integer, intent(in) :: nrun

      self%nrun = nrun
      allocate(self%xdata(lsp,lst))
      allocate(self%zdata(lsp,lst))
      allocate(self%power2d(iks,lsp,lst,0:navr))
      return
   end subroutine init

   subroutine ld_data(self, fid)
      class(csdata), intent(inout) :: self
      integer, intent(in) :: fid
      return
   end subroutine ld_data

   subroutine free(self)
      type(csdata), intent(inout) :: self

      deallocate(self%xdata)
      deallocate(self%zdata)
      deallocate(self%power2d)
      return
   end subroutine free

end module class_csdata
