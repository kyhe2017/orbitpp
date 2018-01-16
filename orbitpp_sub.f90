
subroutine loss
  use orbitpp_private
  use mod_interfaces
  implicit none
  integer :: m, nloss, fid=7
  class(diagloss), pointer :: dd
  class(Qloss), pointer :: ss

  print*, '#######################################'
  print*, 'Starting particle loss processing...'
  print*, ''

  open(14, file='cu2pu', status='old', action='read')
  open(fid, file='lossdata', status='old', action='read')

  do m = 1, nshot
! --- load units for every run
     call load_units(14)
     
     read(fid, *) nloss

     if (nloss > 0) then
        allocate(dd)
        allocate(ss)
     
        call dd % init(m, nloss)
     
        call dd % ld_data(fid)

! ---
        dd%time = dd%time*bax/omeg0*1.d3
        
        call dd % wrt('opt0')
        call dd % wrt('opt')
        call dd % wrt('xx0', 'zz0')
        call dd % wrt('xx', 'zz')
        
        call ss%init(dd, 50, 50)
        
        call ss%wrt('time')
        call ss%count_opt

        deallocate(dd)
        deallocate(ss)
     end if
  end do
  
  close(14)
  close(fid)

  print*, ''
  print*, 'loss distribution processed...'
  print*, ''
  return
end subroutine loss

subroutine dist
  use orbitpp_private
  use mod_interfaces
  implicit none
  integer :: l, n, nlive, fid=7
  real*8  :: t
  class(distribution),     pointer :: dd
  class(aux_distribution), pointer :: aux_dd
  class(statistics),       pointer :: ss

  print*, '#######################################'
  print*, 'Starting particle distribution processing...'
  print*, ''

  open(14, file='cu2pu', status='old', action='read')
  open(fid, file='pdata', status='old', action='read')

  do n = 1, nshot + 1
! --- load units for every run
     if (n /= 2) call load_units(14)

     allocate(aux_dd)

     call aux_dd % init(n-1, 0, 0)

     do l = 1, iks
        read(fid, *) t, nlive
        print*, ' t=', t, ' ltag=', l,' nlive=', nlive

        if (nlive > 0) then
           allocate(dd)
           allocate(ss)
           !allocate(aux_statistics :: ss)

           call dd % init(n-1, l, nlive)

           call dd % ld_data(fid)

! --- C.U. to P.U.            
           dd % pol = dd % pol/pw
           dd % xx = dd % xx*(rmaj/xc)
           dd % zz = dd % zz*(rmaj/xc)
           dd % vpar = dd % vpar*(rmaj/xc)*(omeg0/bax)
           dd % vper = dd % vper*(rmaj/xc)*(omeg0/bax)
           dd % vr = dd % vr*(rmaj/xc)*(omeg0/bax)
           !dd % en = dd % en*(ekev/engn)
           !dd % rmu = dd % rmu*(ekev/engn)/(bkg/bax)
           dd % rmu = dd % rmu*bax/dd % en
           dd % tm = dd % tm/pw*icurrorbit

! --- Data stack for all spiecies            
           call aux_dd % stacks(dd)

! --- write output file            
           call dd % wrt('zz')
           call dd % wrt('xx', 'zz')           
           call dd % wrt('xx', 'zz', 'en')           
           call dd % wrt('tm', 'rmu')

           call ss % init(dd, 50, 50)

           call ss % wrt('zz')
           call ss % wrt('xx', 'zz')
           call ss % wrt('xx', 'zz', 'en')
           call ss % count_opt      
           call ss % count_reg

           deallocate(dd)
           deallocate(ss)
        end if
     end do

! --- write output for all spiecies      
     call aux_dd % wrt('xx', 'zz')

     deallocate(aux_dd)
  end do

  close(14)
  close(fid)

  print*, ''
  print*, 'particle distribution processed...'
  print*, ''
  return
end subroutine dist

subroutine pmomenta
  use orbitpp_private
  use mod_interfaces
  implicit none
  integer :: m, fid=7
  real*8  :: pwscl, totvol, den0, ekev2erg, beta0
  class(momenta), pointer :: mm

  print*, '#######################################'
  print*, 'Starting to process momenta history...'
  print*, ''

  open(14, file='cu2pu', status='old', action='read')
  open(fid, file='wdata', status='old', action='read')

  do m = 1, nshot
     print*, 'momenta history for run', m
! --- load units for every run
     call load_units(14)

     pwscl = (1.6021892d-22*ekev/engn)/(bax/omeg0)  ! Unit in MW
     den0  = 1.d0/(rmaj/xc)**3                      ! Unit in cm^-3
     ekev2erg = 1.6021892d-16*ekev*1e7
     beta0 = 4.0*pi2*den0*ekev2erg/engn/(bkg*1e3)**2

     allocate(mm)

     call mm % init(m)

     call mm % ld_data(fid)

! --- C.U. to P.U.
     mm % powermurf = mm % powermurf*pwscl
     mm % power     = mm % power*pwscl
     mm % powerloss = mm % powerloss*pwscl

     mm % pmom(:,iden,:,:) = mm % pmom(:,iden,:,:)*den0               ! cm^-3
     mm % pmom(:,ieng,:,:) = mm % pmom(:,ieng,:,:)*(ekev/engn)*den0   ! keV/cm^3
     mm % pmom(:,ipar,:,:) = mm % pmom(:,ipar,:,:)*(ekev/engn)*den0   ! keV/cm^3
     mm % pmom(:,iper,:,:) = mm % pmom(:,iper,:,:)*(ekev/engn)*den0   ! keV/cm^3
     mm % pmom(:,ipwr,:,:) = mm % pmom(:,ipwr,:,:)*den0*pwscl         ! MW/cm^3
     mm % pmom(:,ipwrnrf,:,:) = mm % pmom(:,ipwrnrf,:,:)*den0*pwscl   ! MW/cm^3
     mm % pmom(:,ipwrloss,:,:) = mm % pmom(:,ipwrloss,:,:)*den0*pwscl ! MW/cm^3
     mm % pmom(:,isrc,:,:) = mm % pmom(:,isrc,:,:)*den0/(bax/omeg0)   ! sec^-1*cm^-3
     mm % pmom(:,isre,:,:) = mm % pmom(:,isre,:,:)*(ekev/engn)*den0/(bax/omeg0) ! keV*sec^-1*cm^-3
     mm % pmom(:,ibet,:,:) = mm % pmom(:,ibet,:,:)*beta0

! --- write output files you need
     call mm % wrt('den')

     call mm % wrt('eng', navr-1)
     call mm % wrt('par', navr-1)
     call mm % wrt('per', navr-1)

     call mm % wrt_power

! --- print deposited power
     print*,'RF power (MW) :',        mm % powermurf(navr-1)
     print*,'deposited power (MW) :', mm % power(navr-1)
     print*,'lost power (MW) :',      mm % powerloss(navr-1)

     deallocate(mm)

  end do

  close(fid)
  close(14)

  print*, ''
  print*, 'Particle momenta processed...'
  print*, ''
  return
end subroutine pmomenta



