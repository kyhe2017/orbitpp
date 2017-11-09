! --- The derived class is introduced to collect particle distribution 
! --- of all spieces.
! ---  self%time with same meanning as base class; self%ltag=0 donating 
! ---  sum of 1, 2, 3, 4; self%nlive as offset for stacking.
module class_aux_distribution
  
  use orbitpp_private, only : nprt0
  
  use class_distribution
  
  implicit none

  private

  type, extends(distribution), public :: aux_distribution
   contains
     procedure, pass :: init
     procedure, pass :: stacks
  end type aux_distribution
   
contains

  subroutine init(self, time, ltag, nlive)
    class(aux_distribution) ,intent(inout) :: self
    integer, intent(in) :: time, ltag, nlive

    self%time  = time
    self%ltag  = 0
    self%nlive = 0
    allocate(self%nout(nprt0))
    allocate(self%opt(nprt0))
    allocate(self%reg(nprt0))
    allocate(self%pol(nprt0))
    allocate(self%thet(nprt0))
    allocate(self%zet(nprt0))
    allocate(self%xx(nprt0))
    allocate(self%zz(nprt0))
    allocate(self%vpar(nprt0))
    allocate(self%vper(nprt0))
    allocate(self%vr(nprt0))
    allocate(self%en(nprt0))
    allocate(self%rmu(nprt0))
    allocate(self%tm(nprt0))
    allocate(self%ptch(nprt0))
    allocate(self%wegtF(nprt0))
    allocate(self%enloss(nprt0))
    allocate(self%ptchloss(nprt0))
    allocate(self%zloss(nprt0))
    allocate(self%philoss(nprt0))
    allocate(self%xloss(nprt0))
    allocate(self%thetloss(nprt0))
    allocate(self%timeloss(nprt0))
    return
  end subroutine init

  subroutine stacks(self, other)
    class(aux_distribution) ,intent(inout) :: self
    class(distribution), intent(in) :: other
    integer :: i, j

    i = self%nlive + 1
    j = self%nlive + other%nlive

    self%nout(i:j) = other%nout
    self%opt(i:j)  = other%opt
    self%reg(i:j)  = other%reg
    self%pol(i:j)  = other%pol
    self%thet(i:j) = other%thet
    self%zet(i:j)  = other% zet
    self%xx(i:j)   = other%xx
    self%zz(i:j)   = other%zz
    self%vpar(i:j) = other%vpar
    self%vper(i:j) = other%vper
    self%vr(i:j)   = other%vr
    self%en(i:j)   = other%en
    self%rmu(i:j)  = other%rmu
    self%tm(i:j)   = other%tm
    self%ptch(i:j) = other%ptch
    self%wegtF(i:j)    = other%wegtF
    self%enloss(i:j)   = other%enloss
    self%ptchloss(i:j) = other%ptchloss
    self%zloss(i:j)    = other%zloss
    self%philoss(i:j)  = other%philoss
    self%xloss(i:j)    = other%xloss
    self%thetloss(i:j) = other%thetloss
    self%timeloss(i:j) = other%timeloss

    self%nlive = self%nlive + other%nlive
    return
  end subroutine stacks

end module class_aux_distribution
