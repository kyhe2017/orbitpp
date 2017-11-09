! --- Unit conversion from C.U. to P.U. & Preparation for output
! ---  orbitpp serves as an output processing program for ORBIT-RF,
! ---  and the old version is revised for more convenient usage.
! ---  Moreover, the new version is object oriented by use of 
! ---  Fortran standard 2003.
! ---  K.Y. He, May 19th, 2017
program orbitpp
  use orbitpp_private
  implicit none

! --- number of shots
  open(7, file='timer', status='old', action='read')

  read(7, *) nshot

  close(7)

  print*, 'nshot = ', nshot
  nshot = nshot + 1 ! removed when run

  call print_icfg

! --- magnetic geometry
  call wrt_geometry

! --- particle distributions
  call dist

! --- particle trajactories
  call traj

! --- particle momenta and power
  call pmomenta

contains

  subroutine print_icfg
    implicit none

    open(14, file='cu2pu', status='old', action='read')

    call load_units(14)

    close(14)

! --- print initial configurations
    write(*, *) 'Initial configurations :'
    write(*, *) 'nprt0 = ', nprt0
    write(*, *) 'navr = ', navr
    write(*, *) 'lsp = ', lsp
    write(*, *) 'lst = ', lst
    write(*, *) 'pw = ', pw
    write(*, *) 'prot_maj=', protma
    write(*, *) 'ptest = ', ptest
    write(*, *) 'tran = ', tran
    write(*, *) 'trun = ', trun
    write(*, *) 'tavr = ', tavr
    write(*, *) 'omeg0 = ', omeg0
    write(*, *) 'xc = ', xc
    write(*, *) 'bax = ', bax
    write(*, *) 'bkg (kG) = ', bkg
    write(*, *) 'engn = ', engn
    write(*, *) 'ekev (keV) = ', ekev
    write(*, *) 'dene0 (* 1.0e13 cm-3) = ', dene0/1.0e13
    write(*, *) 'Te0 (keV) = ', Te0
    write(*, *) 'Ti0 (keV) = ', Ti0
    write(*, *) 'ndist = ', ndist
    write(*, *) 'xright (cm) = ', xright*rmaj/xc
    write(*, *) 'xleft (cm) = ', xleft*rmaj/xc
    write(*, *) 'ztop (cm) = ', ztop*rmaj/xc
    write(*, *) 'zbot (cm) = ', zbot*rmaj/xc
    return
  end subroutine print_icfg

end program orbitpp



