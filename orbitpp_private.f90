module orbitpp_private
  implicit none

  integer :: nshot

  real*8 :: pi   = 4.d0*datan(1.d0)
  real*8 :: pi2  = 8.d0*datan(1.d0)

! --- variables below in synchronization with globals.f90
  integer :: iks = 4
  integer :: imo = 12

  integer :: iden = 1
  integer :: ieng = 2
  integer :: ipwr = 3
  integer :: idrg = 4
  integer :: iflx = 5
  integer :: isrc = 6
  integer :: isre = 7
  integer :: ibet = 8
  integer :: ipwrnrf  = 9
  integer :: ipwrloss = 10
  integer :: ipar  = 11
  integer :: iper = 12

! --- variables for unit conversion loaded for every run
  integer :: &
       & nprt0, ndist, &
       & protma, zprtma, ptest, ztest, &
       & ntor, navr, &
       & lsp, lst, &
       & icurrorbit, ibtororbit

  real*8  :: &
       & trun, tavr, tran, &
       & pw, rmaj, xc, bkg, bax, &
       & omeg0, engn, ekev, &
       & dene0, Te0, Ti0, &
       & xright, xleft, ztop, zbot

contains

  subroutine load_units(fid)
    integer, intent(in) :: fid

    read(fid, *) nprt0, ndist
    read(fid, *) protma, zprtma, ptest, ztest
    read(fid, *) ntor, navr, trun, tavr, tran
    read(fid, *) lsp, lst, pw
    read(fid, *) rmaj, xc, bkg, bax, omeg0
    read(fid, *) engn, ekev
    read(fid, *) dene0, Te0, Ti0
    read(fid, *) xright, xleft, ztop, zbot
    read(fid, *) icurrorbit, ibtororbit
    return
  end subroutine load_units

end module orbitpp_private
