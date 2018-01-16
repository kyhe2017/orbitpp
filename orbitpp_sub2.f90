subroutine wrt_geometry
  use orbitpp_private
  use mod_interfaces
  implicit none
  integer :: nsp, nst, fid=7
  class(geometry), pointer :: gg

  print*, '#######################################'
  print*, 'Starting magnetic geometry processing...'
  print*, ''

! --- load units conversion(only onece)
  open(14, file='cu2pu', status='old', action='read')

  call load_units(14)

! --- magnetic geometry
  open(fid, file='gdata', status='old', action='read')

  read(fid, '(2I5)') nsp, nst

  allocate(gg)

  call gg % init(nsp, nst)
  call gg % ld_data(fid)

! --- C.U. to P.U.
  gg % psim  = gg % psim/pw
  gg % xm    = gg % xm*rmaj/xc
  gg % zm    = gg % zm*rmaj/xc
  gg % pom   = gg % pom/pw
  gg % bm    = gg % bm*(bkg/bax)
  gg % bripm = gg % bripm/bax

! --- 1D and 2D output
  call gg % wrt('psip', 'q')
  call gg % wrt('psip')
  call gg % wrt('b')

  deallocate(gg)
  
  close(14)
  close(fid)

  print*, ''
  print*, 'Magnetic geometry processed...'
  print*, ''
  return
end subroutine wrt_geometry

subroutine traj
  use orbitpp_private
  use mod_interfaces
  implicit none
  integer :: n, k, itag, nlast, ntrj, fid=7
  class(trajectory), pointer :: tt

  print*, '#######################################'
  print*, 'Starting particle trajectories processing...'
  print*, ''

  open(14, file='cu2pu', status='old', action='read')
  open(fid, file='tdata', status='old', action='read')

  do n = 1, nshot
! --- load units for every run
     call load_units(14)

     read(fid, *) ntrj
     if (ntrj > 999) stop 'ntrj>999 --- stop'

     print*, 'Processing ', ntrj, ' trajectories for run', n

     do k = 1, ntrj
        read(fid, *) itag, nlast
        print*,'particle=', itag,' last step=', nlast

        allocate(tt)

        call tt % init(n, itag, nlast)         

        call tt % ld_data(fid)

! --- C.U. to P.U.
        tt % tp = tt % tp*(bax/omeg0)
        tt % xp = tt % xp*(rmaj/xc)
        tt % zp = tt % zp*(rmaj/xc)
        tt % rmup = tt % rmup*(ekev/engn)/(bkg/bax)
        tt % enp = tt % enp*(ekev/engn)
        tt % bp = tt % bp*(bkg/bax)
        tt % vthetp = tt % vthetp*(omeg0/bax)
        tt % plp = tt % plp/pw
        tt % torp = tt % torp/pw*icurrorbit

! --- write data file you need         
        call tt % wrt('pol')
        call tt % wrt('x', 'z')

        deallocate(tt)
     end do
  end do

  close(14)
  close(fid)

  print*, ''
  print*, 'particle trajectories processed...'
  print*, ''
  return
end subroutine traj
