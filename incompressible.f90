 program incompressible
!**********************************************
!this is initial-version of the generalized solver for the simulation of incompressible-flows.
!this is a code written in Fortran-90 by Haroon Ahmad (2015-AMZ-8217)
!for submission to Dr. Balaji Srinivasan in APL-720.
!this program uses the correction-pressure based modified SMAC scheme (refer Hasan & Sanghi,JHT,126,963-984) alongwith QUICK-upwind to deal convective-terms.
!this program uses Rhie-Chow momentum interpolation to deal with the grid-scale oscillations due to CPPE.
!**********************************************
  implicit none
!**********************************************
 integer, parameter :: maxstep = 1000000
 double precision, parameter :: dt=1.0d-03
 integer, parameter :: maxiter = 100000
 integer, parameter :: snapshot = 1000
 double precision, parameter :: errtol=1.0d-03
 integer :: Nx,Ny,nstep 
 integer :: i,j,jcnt
 double precision, parameter :: Re = 100.0d0
 double precision :: xcell_Pe,ycell_Pe
 double precision, allocatable, dimension(:) :: x,y,de,dw,dee,dww,dn,ds,dnn,dss
 double precision, allocatable, dimension(:) :: etaEE_2cen,etaE_2cen,etaW_2cen,etaWW_2cen
 double precision, allocatable, dimension(:) :: etaNN_2cen,etaN_2cen,etaS_2cen,etaSS_2cen
 double precision, allocatable, dimension(:) :: etaEE_1cen,etaE_1cen,etaW_1cen,etaWW_1cen
 double precision, allocatable, dimension(:) :: etaNN_1cen,etaN_1cen,etaS_1cen,etaSS_1cen
 double precision, allocatable, dimension(:) :: etaE_gradP,etaW_gradP,etaN_gradP,etaS_gradP
 double precision, allocatable, dimension(:) :: etaEE_upw_pos,etaE_upw_pos,etaW_upw_pos,etaWW_upw_pos  
 double precision, allocatable, dimension(:) :: etaNN_upw_pos,etaN_upw_pos,etaS_upw_pos,etaSS_upw_pos  
 double precision, allocatable, dimension(:) :: etaEE_upw_neg,etaE_upw_neg,etaW_upw_neg,etaWW_upw_neg  
 double precision, allocatable, dimension(:) :: etaNN_upw_neg,etaN_upw_neg,etaS_upw_neg,etaSS_upw_neg 
 double precision, allocatable, dimension(:,:) :: Avp,Ave,Avw,Avn,Avs,Avee,Avww,Avnn,Avss 
 double precision :: de_bnd,dw_bnd,dee_bnd,dww_bnd,dn_bnd,ds_bnd,dnn_bnd,dss_bnd
 double precision :: etaE_bnd2,etaEE_bnd2
 double precision :: etaW_bnd2,etaWW_bnd2
 double precision :: etaN_bnd2,etaNN_bnd2
 double precision :: etaS_bnd2,etaSS_bnd2
 double precision :: etaE,etaEE
 double precision :: etaW,etaWW
 double precision :: etaN,etaNN
 double precision :: etaS,etaSS,ae,aw,an,as,ap
 double precision, allocatable, dimension(:,:) :: App,Ape,Apw,Apn,Aps 
 double precision, allocatable, dimension(:,:,:) :: Ut,forc,Q,conv
 double precision, allocatable, dimension(:,:) :: press,pcorr,Qp
 double precision, dimension(2) :: gradP
 double precision :: dx1,dx2,dy1,dy2 
 double precision :: time,aa
 character(len=30) :: filname
 integer :: waste
!**********************************************
   aa = dt/Re
   print*, aa
!*****read the grid-data
  open(unit=1,file='mesh.xy',status='unknown')

  read(1,*) Nx,Ny

  print*, Nx,Ny
!**************allocate the flow variables***************
     
   allocate( x(Nx) )
   allocate( y(Ny) )
   allocate( de(Nx) )
   allocate( dw(Nx) )
   allocate( dee(Nx) )
   allocate( dww(Nx) )
   allocate( dn(Ny) )
   allocate( ds(Ny) )
   allocate( dnn(Ny) )
   allocate( dss(Ny) )
   allocate( Ut(Nx,Ny,4) )
   allocate( press(Nx,Ny) )
   allocate( pcorr(Nx,Ny) )
   allocate( forc(Nx,Ny,2) )
   allocate( Q(Nx,Ny,6) )
   allocate( conv(Nx,Ny,2) )
   allocate( Qp(Nx,Ny) )
   allocate( etaEE_2cen(Nx) )
   allocate( etaE_2cen(Nx) )
   allocate( etaW_2cen(Nx) )
   allocate( etaWW_2cen(Nx) )
   allocate( etaNN_2cen(Ny) )
   allocate( etaN_2cen(Ny) )
   allocate( etaS_2cen(Ny) )
   allocate( etaSS_2cen(Ny) )
   allocate( etaEE_1cen(Nx) )
   allocate( etaE_1cen(Nx) )
   allocate( etaW_1cen(Nx) )
   allocate( etaWW_1cen(Nx) )
   allocate( etaNN_1cen(Ny) )
   allocate( etaN_1cen(Ny) )
   allocate( etaS_1cen(Ny) )
   allocate( etaSS_1cen(Ny) )
   allocate( etaEE_upw_pos(Nx) )
   allocate( etaE_upw_pos(Nx) )
   allocate( etaW_upw_pos(Nx) )
   allocate( etaWW_upw_pos(Nx) )
   allocate( etaNN_upw_pos(Ny) )
   allocate( etaN_upw_pos(Ny) )
   allocate( etaS_upw_pos(Ny) )
   allocate( etaSS_upw_pos(Ny) )
   allocate( etaEE_upw_neg(Nx) )
   allocate( etaE_upw_neg(Nx) )
   allocate( etaW_upw_neg(Nx) )
   allocate( etaWW_upw_neg(Nx) )
   allocate( etaNN_upw_neg(Ny) )
   allocate( etaN_upw_neg(Ny) )
   allocate( etaS_upw_neg(Ny) )
   allocate( etaSS_upw_neg(Ny) )
   allocate( etaE_gradP(Nx) )
   allocate( etaW_gradP(Nx) )
   allocate( etaN_gradP(Ny) )
   allocate( etaS_gradP(Ny) )
   allocate( App(Nx,Ny) )
   allocate( Ape(Nx,Ny) )
   allocate( Apw(Nx,Ny) )
   allocate( Apn(Nx,Ny) )
   allocate( Aps(Nx,Ny) )
   allocate( Avp(Nx,Ny) )
   allocate( Ave(Nx,Ny) )
   allocate( Avw(Nx,Ny) )
   allocate( Avn(Nx,Ny) )
   allocate( Avs(Nx,Ny) )
   allocate( Avee(Nx,Ny) )
   allocate( Avww(Nx,Ny) )
   allocate( Avnn(Nx,Ny) )
   allocate( Avss(Nx,Ny) )

!********************************************************
     do i=1,Nx
   read(1,11) waste,x(i)
     end do
 
     do j=1,Ny
   read(1,11) waste,y(j)
     end do

  close(1)
11 format(1X,I4,1X,F23.15)

 print*, "grid data read" 
!********************************************************
!***debug-statement
!  print*, 'first',Nx,Ny,Re
!  print*, 'second',dt,maxiter
!   print*, 'third',errtol
!pause
!*****calculate the step-sizes

       do i=2,Nx-1
   de(i) = x(i+1)-x(i)
   dw(i) = x(i)-x(i-1)
  if(i.ge.3.and.i.le.(Nx-2))then
  dee(i) = x(i+2)-x(i)
  dww(i) = x(i)-x(i-2)
  else
  dee(i) = 0.0d0
  dww(i) = 0.0d0
  end if
       end do
!***debug-statement
!  print*, 'step-size',de(2),dw(2)
!  print*, 'step-size',dee(2),dww(2)
!  print*, 'step-size',de(5),dw(5)
!  print*, 'step-size',dee(5),dww(5)
!  print*, 'step-size',de(81),dw(81)
!  print*, 'step-size',dee(81),dww(81)
!  print*, 'step-size',de(98),dw(98)
! do i=1,Nx
!  if(de(i).eq.0.0d0) print*, 'de',i
!  if(dw(i).eq.0.0d0) print*, 'dw',i
!  if(dee(i).eq.0.0d0) print*, 'dee',i
!  if(dww(i).eq.0.0d0) print*, 'dww',i
! end do
!pause

       do j=2,Ny-1
   dn(j) = y(j+1)-y(j)
   ds(j) = y(j)-y(j-1)
  if(j.ge.3.and.j.le.(Ny-2))then
  dnn(j) = y(j+2)-y(j)
  dss(j) = y(j)-y(j-2)
  else
  dnn(j) = 0.0d0
  dss(j) = 0.0d0
  end if
       end do

!  print*, 'step-size',dn(5),ds(5)
!  print*, 'step-size',dnn(5),dss(5)
!  print*, 'step-size',dn(81),ds(81)
!  print*, 'step-size',dnn(81),dss(81)
! do j=1,Ny
!  if(dn(j).eq.0.0d0) print*, 'dn',j
!  if(ds(j).eq.0.0d0) print*, 'ds',j
!  if(dnn(j).eq.0.0d0) print*, 'dnn',j
!  if(dss(j).eq.0.0d0) print*, 'dss',j
! end do
!pause
!********************************************************
!*****calculate the weights for different discretization schemes
!********************************************************
!
!*****for one-sided approximations of the first and second-order derivative at boundaries
!****************************************
         de_bnd = x(2)-x(1)
        dee_bnd = x(3)-x(1)
!****************************************
       etaE = dee_bnd/(de_bnd*(dee_bnd-de_bnd))
      etaEE = (-1.0d0*de_bnd)/(dee_bnd*(dee_bnd-de_bnd))

       etaE_bnd2 = -2.0d0/(de_bnd*(dee_bnd-de_bnd))
      etaEE_bnd2 = 2.0d0/(dee_bnd*(dee_bnd-de_bnd))
!****************************************
!  print*, 'east',de_bnd,dee_bnd,etaE,etaEE,etaE_bnd2,etaEE_bnd2
! pause
!****************************************
         dw_bnd = x(Nx)-x(Nx-1)
        dww_bnd = x(Nx)-x(Nx-2) 
!****************************************
       etaW = (-1.0d0*dww_bnd)/(dw_bnd*(dww_bnd-dw_bnd))
      etaWW = dw_bnd/(dww_bnd*(dww_bnd-dw_bnd))

       etaW_bnd2 = -2.0d0/(dw_bnd*(dww_bnd-dw_bnd))
      etaWW_bnd2 = 2.0d0/(dww_bnd*(dww_bnd-dw_bnd))
!****************************************
!  print*, 'west',dw_bnd,dww_bnd,etaW,etaWW,etaW_bnd2,etaWW_bnd2
!pause
!****************************************
        dn_bnd = y(2)-y(1)
       dnn_bnd = y(3)-y(1)
!****************************************
       etaN = dnn_bnd/(dn_bnd*(dnn_bnd-dn_bnd))
      etaNN = (-1.0d0*dn_bnd)/(dnn_bnd*(dnn_bnd-dn_bnd))

       etaN_bnd2 = -2.0d0/(dn_bnd*(dnn_bnd-dn_bnd))
      etaNN_bnd2 = 2.0d0/(dnn_bnd*(dnn_bnd-dn_bnd))
!****************************************
!  print*, 'north',dn_bnd,dnn_bnd,etaN,etaNN,etaN_bnd2,etaNN_bnd2
!pause
!****************************************
        ds_bnd = y(Ny)-y(Ny-1)
       dss_bnd = y(Ny)-y(Ny-2)
!****************************************
       etaS = (-1.0d0*dss_bnd)/(ds_bnd*(dss_bnd-ds_bnd))
      etaSS = ds_bnd/(dss_bnd*(dss_bnd-ds_bnd))

       etaS_bnd2 = -2.0d0/(ds_bnd*(dss_bnd-ds_bnd))
      etaSS_bnd2 = 2.0d0/(dss_bnd*(dss_bnd-ds_bnd))
!****************************************
!  print*, 'south',ds_bnd,dss_bnd,etaS,etaSS,etaS_bnd2,etaSS_bnd2
!pause
!****************************************
!
!*****for Second derivative central-approximations 
!***close all 4th order discrete-molecules for time-being.
       do i=2,Nx-1
!  if(i.ge.4.and.i.le.(Nx-3))then
!   etaWW_2cen(i) = -2.0d0*( (dee(i)*dw(i))-(dee(i)*de(i))+(de(i)*dw(i)) )/( dww(i)*(dww(i)-dw(i))*(dww(i)+dee(i))*(dww(i)+de(i)) )  
!    etaW_2cen(i) = 2.0d0*( (dee(i)*dww(i))-(dee(i)*de(i))+(de(i)*dww(i)) )/( dw(i)*(dww(i)-dw(i))*(dw(i)+dee(i))*(dw(i)+de(i)) )
!    etaE_2cen(i) = 2.0d0*( (dee(i)*dw(i))-(dww(i)*dw(i))+(dee(i)*dww(i)) )/( de(i)*(dee(i)-de(i))*(dww(i)+de(i))*(dw(i)+de(i)) )
!   etaEE_2cen(i) = -2.0d0*( (de(i)*dw(i))-(dww(i)*dw(i))+(de(i)*dww(i)) )/( dee(i)*(dee(i)-de(i))*(dw(i)+dee(i))*(dww(i)+dee(i)) )
!  else
!   etaEE_2cen(i) = 0.0d0
    etaE_2cen(i) = 2.0d0/( de(i)*( de(i)+dw(i) ) )
    etaW_2cen(i) = 2.0d0/( dw(i)*( de(i)+dw(i) ) )
!   etaWW_2cen(i) = 0.0d0
!  end if
       end do
!***debug-statement
!  print*, 'weight-E',etaE_2cen(5),etaEE_2cen(5)
!  print*, 'weight-W',etaW_2cen(5),etaWW_2cen(5)
! do i=1,Nx
!  if(etaEE_2cen(i).eq.0.0d0) print*, 'EE',i
!  if(etaE_2cen(i).eq.0.0d0) print*, 'E',i
!  if(etaW_2cen(i).eq.0.0d0) print*, 'W',i
!  if(etaWW_2cen(i).eq.0.0d0) print*, 'WW',i
! end do
!pause
!**********************
       do j=2,Ny-1
!***close all 4th order discrete-molecules for time-being.
!  if(j.ge.4.and.j.le.(Ny-3))then
!   etaSS_2cen(j) = -2.0d0*( (dnn(j)*ds(j))-(dnn(j)*dn(j))+(dn(j)*ds(j)) )/( dss(j)*(dss(j)-ds(j))*(dss(j)+dnn(j))*(dss(j)+dn(j)) )  
!    etaS_2cen(j) = 2.0d0*( (dnn(j)*dss(j))-(dnn(j)*dn(j))+(dn(j)*dss(j)) )/( ds(j)*(dss(j)-ds(j))*(ds(j)+dnn(j))*(ds(j)+dn(j)) )
!    etaN_2cen(j) = 2.0d0*( (dnn(j)*ds(j))-(dss(j)*ds(j))+(dnn(j)*dss(j)) )/( dn(j)*(dnn(j)-dn(j))*(dss(j)+dn(j))*(ds(j)+dn(j)) )
!   etaNN_2cen(j) = -2.0d0*( (dn(j)*ds(j))-(dss(j)*ds(j))+(dn(j)*dss(j)) )/( dnn(j)*(dnn(j)-dn(j))*(ds(j)+dnn(j))*(dss(j)+dnn(j)) )
!  else
!   etaNN_2cen(j) = 0.0d0
    etaN_2cen(j) = 2.0d0/( dn(j)*( dn(j)+ds(j) ) )
    etaS_2cen(j) = 2.0d0/( ds(j)*( dn(j)+ds(j) ) )
!   etaSS_2cen(j) = 0.0d0
!  end if
       end do
!***debug-statement
!  print*, 'weight-E',etaE_2cen(2),etaEE_2cen(2)
!  print*, 'weight-W',etaW_2cen(2),etaWW_2cen(2)
!  print*, 'weight-N',etaN_2cen(2),etaNN_2cen(2)
!  print*, 'weight-S',etaS_2cen(2),etaSS_2cen(2)
!  print*, 'weight-N',etaN_2cen(5),etaNN_2cen(5)
!  print*, 'weight-S',etaS_2cen(5),etaSS_2cen(5)
! do j=1,Ny
!  if(etaNN_2cen(j).eq.0.0d0) print*, 'NN',j
!  if(etaN_2cen(j).eq.0.0d0) print*, 'N',j
!  if(etaS_2cen(j).eq.0.0d0) print*, 'S',j
!  if(etaSS_2cen(j).eq.0.0d0) print*, 'SS',j
! end do
!pause
!**********************

!*****for First derivative central-approximations
       do i=2,Nx-1
  if(i.ge.3.and.i.le.(Nx-2))then
   etaWW_1cen(i) = ( dee(i)*dw(i)*de(i) )/( dww(i)*(dww(i)-dw(i))*(dww(i)+dee(i))*(dww(i)+de(i)) )  
    etaW_1cen(i) = ( -dee(i)*dww(i)*de(i) )/( dw(i)*(dww(i)-dw(i))*(dw(i)+dee(i))*(dw(i)+de(i)) )
    etaE_1cen(i) = ( dee(i)*dw(i)*dww(i) )/( de(i)*(dee(i)-de(i))*(dww(i)+de(i))*(dw(i)+de(i)) )
   etaEE_1cen(i) = ( -de(i)*dw(i)*dww(i) )/( dee(i)*(dee(i)-de(i))*(dw(i)+dee(i))*(dww(i)+dee(i)) )
  end if
       end do
!       do i=2,Nx-1
!    etaE_1Lcen(i) = dw(i)/( de(i)*( de(i)+dw(i) ) )
!    etaW_1Lcen(i) = -de(i)/( dw(i)*( de(i)+dw(i) ) )
!       end do
!***debug-statement
!  print*, 'weight-E',etaE_1cen(5),etaEE_1cen(5)
!  print*, 'weight-W',etaW_1cen(5),etaWW_1cen(5)
! do i=1,Nx
!  if(etaEE_1cen(i).eq.0.0d0) print*, 'EE_1',i
!  if(etaE_1cen(i).eq.0.0d0) print*, 'E_1',i
!  if(etaW_1cen(i).eq.0.0d0) print*, 'W_1',i
!  if(etaWW_1cen(i).eq.0.0d0) print*, 'WW_1',i
! end do
!pause
!**********************
       do j=2,Ny-1
  if(j.ge.3.and.j.le.(Ny-2))then
   etaSS_1cen(j) = ( dnn(j)*ds(j)*dn(j) )/( dss(j)*(dss(j)-ds(j))*(dss(j)+dnn(j))*(dss(j)+dn(j)) )  
    etaS_1cen(j) = ( -dnn(j)*dss(j)*dn(j) )/( ds(j)*(dss(j)-ds(j))*(ds(j)+dnn(j))*(ds(j)+dn(j)) )
    etaN_1cen(j) = ( dnn(j)*ds(j)*dss(j) )/( dn(j)*(dnn(j)-dn(j))*(dss(j)+dn(j))*(ds(j)+dn(j)) )
   etaNN_1cen(j) = ( -dn(j)*ds(j)*dss(j) )/( dnn(j)*(dnn(j)-dn(j))*(ds(j)+dnn(j))*(dss(j)+dnn(j)) )
  end if
       end do
!       do j=2,Ny-1
!    etaN_1Lcen(j) = ds(j)/( dn(j)*( dn(j)+ds(j) ) )
!    etaS_1Lcen(j) = -dn(j)/( ds(j)*( dn(j)+ds(j) ) )
!  end if
!       end do
!***debug-statement
!  print*, 'weight-E',etaE_1cen(2),etaEE_1cen(2)
!  print*, 'weight-W',etaW_1cen(2),etaWW_1cen(2)
!  print*, 'weight-N',etaN_1cen(2),etaNN_1cen(2)
!  print*, 'weight-S',etaS_1cen(2),etaSS_1cen(2)
!  print*, 'weight-N',etaN_1cen(5),etaNN_1cen(5)
!  print*, 'weight-S',etaS_1cen(5),etaSS_1cen(5)
! do j=1,Ny
!  if(etaNN_1cen(j).eq.0.0d0) print*, 'NN_1',j
!  if(etaN_1cen(j).eq.0.0d0) print*, 'N_1',j
!  if(etaS_1cen(j).eq.0.0d0) print*, 'S_1',j
!  if(etaSS_1cen(j).eq.0.0d0) print*, 'SS_1',j
! end do
!pause
!**********************
       do i=2,Nx-1
    etaE_gradP(i) = dw(i)/( de(i)*( de(i)+dw(i) ) )
    etaW_gradP(i) = -de(i)/( dw(i)*( de(i)+dw(i) ) )
       end do
!**********************
       do j=2,Ny-1
    etaN_gradP(j) = ds(j)/( dn(j)*( dn(j)+ds(j) ) )
    etaS_gradP(j) = -dn(j)/( ds(j)*( dn(j)+ds(j) ) )
       end do
!***debug-statement
!  print*, 'weight-E_P',etaE_gradP(2)
!  print*, 'weight-W_P',etaW_gradP(2)
!  print*, 'weight-N_P',etaN_gradP(2)
!  print*, 'weight-S_P',etaS_gradP(2)
!pause
!*****for First derivative upwind-approximations
!***close FOU 
!***for up>0.0
       do i=2,Nx-1
  if(i.ge.3.and.i.le.(Nx-2))then
   etaWW_upw_pos(i) = ( dw(i)*de(i) )/( dww(i)*(dww(i)-dw(i))*(dww(i)+de(i)) )  
    etaW_upw_pos(i) = ( -dww(i)*de(i) )/( dw(i)*(dww(i)-dw(i))*(dw(i)+de(i)) )
    etaE_upw_pos(i) = ( dw(i)*dww(i) )/( de(i)*(dww(i)+de(i))*(dw(i)+de(i)) )
   etaEE_upw_pos(i) = 0.0d0
!  else
!   etaWW_upw_pos(i) = 0.0d0
!    etaW_upw_pos(i) = -1.0d0/dw(i)
!    etaE_upw_pos(i) = 0.0d0
!   etaEE_upw_pos(i) = 0.0d0
  end if
       end do

!***debug-statement
!  print*, 'weight-upwind-E',etaE_upw_pos(5),etaEE_upw_pos(5)
!  print*, 'weight-upwind-W',etaW_upw_pos(5),etaWW_upw_pos(5)
! do i=1,Nx
!  if(etaEE_upw_pos(i).ne.0.0d0) print*, 'EE_U',i
!  if(etaE_upw_pos(i).eq.0.0d0) print*, 'E_U',i
!  if(etaW_upw_pos(i).eq.0.0d0) print*, 'W_U',i
!  if(etaWW_upw_pos(i).eq.0.0d0) print*, 'WW_U',i
! end do
!pause
!***for vp>0.0
!***close FOU 
       do j=2,Ny-1
  if(j.ge.3.and.j.le.(Ny-2))then
   etaSS_upw_pos(j) = ( ds(j)*dn(j) )/( dss(j)*(dss(j)-ds(j))*(dss(j)+dn(j)) )  
    etaS_upw_pos(j) = ( -dss(j)*dn(j) )/( ds(j)*(dss(j)-ds(j))*(ds(j)+dn(j)) )
    etaN_upw_pos(j) = ( ds(j)*dss(j) )/( dn(j)*(dss(j)+dn(j))*(ds(j)+dn(j)) )
   etaNN_upw_pos(j) = 0.0d0
!  else
!   etaSS_upw_pos(j) = 0.0d0
!    etaS_upw_pos(j) = -1.0d0/ds(j)
!    etaN_upw_pos(j) = 0.0d0
!   etaNN_upw_pos(j) = 0.0d0
  end if
       end do
!***debug-statement
!  print*, 'weight-upwind-E',etaE_upw_pos(2),etaEE_upw_pos(2)
!  print*, 'weight-upwind-W',etaW_upw_pos(2),etaWW_upw_pos(2)
!  print*, 'weight-upwind-N',etaN_upw_pos(2),etaNN_upw_pos(2)
!  print*, 'weight-upwind-S',etaS_upw_pos(2),etaSS_upw_pos(2)
!  print*, 'weight-upwind-N',etaN_upw_pos(5),etaNN_upw_pos(5)
!  print*, 'weight-upwind-S',etaS_upw_pos(5),etaSS_upw_pos(5)
! do j=1,Ny
!  if(etaNN_upw_pos(j).ne.0.0d0) print*, 'NN_U',j
!  if(etaN_upw_pos(j).eq.0.0d0) print*, 'N_U',j
!  if(etaS_upw_pos(j).eq.0.0d0) print*, 'S_U',j
!  if(etaSS_upw_pos(j).eq.0.0d0) print*, 'SS_U',j
! end do
!pause
!***for up<0.0
!***close FOU 
       do i=2,Nx-1
  if(i.ge.3.and.i.le.(Nx-2))then
   etaWW_upw_neg(i) = 0.0d0      
    etaW_upw_neg(i) = ( -dee(i)*de(i) )/( dw(i)*(dw(i)+dee(i))*(dw(i)+de(i)) )
    etaE_upw_neg(i) = ( dee(i)*dw(i) )/( de(i)*(dee(i)-de(i))*(dw(i)+de(i)) )
   etaEE_upw_neg(i) = ( -de(i)*dw(i) )/( dee(i)*(dee(i)-de(i))*(dw(i)+dee(i)) )
!  else
!   etaWW_upw_neg(i) = 0.0d0
!    etaW_upw_neg(i) = 0.0d0        
!    etaE_upw_neg(i) = 1.0d0/de(i)
!   etaEE_upw_neg(i) = 0.0d0
  end if
       end do

!***debug-statement
!  print*, 'weight-upwind-E',etaE_upw_neg(5),etaEE_upw_neg(5)
!  print*, 'weight-upwind-W',etaW_upw_neg(5),etaWW_upw_neg(5)
! do i=1,Nx
!  if(etaEE_upw_neg(i).eq.0.0d0) print*, 'EEU-',i
!  if(etaE_upw_neg(i).eq.0.0d0) print*, 'EU-',i
!  if(etaW_upw_neg(i).eq.0.0d0) print*, 'WU-',i
!  if(etaWW_upw_neg(i).ne.0.0d0) print*, 'WWU-',i
! end do
!pause
!***for vp<0.0
!***close FOU 
       do j=2,Ny-1
  if(j.ge.3.and.j.le.(Ny-2))then
   etaSS_upw_neg(j) = 0.0d0      
    etaS_upw_neg(j) = ( -dnn(j)*dn(j) )/( ds(j)*(ds(j)+dnn(j))*(ds(j)+dn(j)) )
    etaN_upw_neg(j) = ( dnn(j)*ds(j) )/( dn(j)*(dnn(j)-dn(j))*(ds(j)+dn(j)) )
   etaNN_upw_neg(j) = ( -dn(j)*ds(j) )/( dnn(j)*(dnn(j)-dn(j))*(ds(j)+dnn(j)) )
!  else
!   etaSS_upw_neg(j) = 0.0d0
!    etaS_upw_neg(j) = 0.0d0        
!    etaN_upw_neg(j) = 1.0d0/dn(j)
!   etaNN_upw_neg(j) = 0.0d0
  end if
       end do

!***debug-statement
!  print*, 'weight-upwind-E',etaE_upw_neg(2),etaEE_upw_neg(2)
!  print*, 'weight-upwind-W',etaW_upw_neg(2),etaWW_upw_neg(2)
!  print*, 'weight-upwind-N',etaN_upw_neg(2),etaNN_upw_neg(2)
!  print*, 'weight-upwind-S',etaS_upw_neg(2),etaSS_upw_neg(2)
!  print*, 'weight-upwind-N',etaN_upw_neg(5),etaNN_upw_neg(5)
!  print*, 'weight-upwind-S',etaS_upw_neg(5),etaSS_upw_neg(5)
! do j=1,Ny
!  if(etaNN_upw_neg(j).eq.0.0d0) print*, 'NNU-',j
!  if(etaN_upw_neg(j).eq.0.0d0) print*, 'NU-',j
!  if(etaS_upw_neg(j).eq.0.0d0) print*, 'SU-',j
!  if(etaSS_upw_neg(j).ne.0.0d0) print*, 'SSU-',j
! end do
!pause
!********************************************************
!*****calculate the Stencil-parameters for different elliptic-equations.
!********************************************************
!*****for the transport-equation of velocity
!***use second order central discretization
      do i=2,Nx-1
      do j=2,Ny-1
!****************************** 
!   Avp(i,j) = -1.0d0*( etaE_2cen(i)+etaW_2cen(i)+etaN_2cen(j)+etaS_2cen(j) ) 
!   Ave(i,j) = etaE_2cen(i)
!   Avw(i,j) = etaW_2cen(i)
!   Avn(i,j) = etaN_2cen(j)
!   Avs(i,j) = etaS_2cen(j)
!   Avp(i,j) = 1.0d0+( (dt/Re)*(etaEE_2cen(i)+etaE_2cen(i)+etaW_2cen(i)+etaWW_2cen(i)+etaNN_2cen(j)+etaN_2cen(j)+etaS_2cen(j)+etaSS_2cen(j)) )
   Avp(i,j) = 1.0d0+( aa*( etaE_2cen(i)+etaW_2cen(i)+etaN_2cen(j)+etaS_2cen(j) ) )
   Ave(i,j) = (-aa)*etaE_2cen(i)
   Avw(i,j) = (-aa)*etaW_2cen(i)
   Avn(i,j) = (-aa)*etaN_2cen(j)
   Avs(i,j) = (-aa)*etaS_2cen(j)
!  Avee(i,j) = (-dt/Re)*etaEE_2cen(i)
!  Avww(i,j) = (-dt/Re)*etaWW_2cen(i)
!  Avnn(i,j) = (-dt/Re)*etaNN_2cen(j)
!  Avss(i,j) = (-dt/Re)*etaSS_2cen(j)
!****************************** 
       end do
       end do

!***debug-statement
!  print*, 'diffusion-Stencil-E',Ave(2,2)
!  print*, 'diffusion-Stencil-W',Avw(2,2)
!  print*, 'diffusion-Stencil-N',Avn(2,2)
!  print*, 'diffusion-Stencil-S',Avs(2,2)
!  print*, 'diffusion-Stencil-P',Avp(2,2)
!pause
!*****for the poisson-equation of correction-pressure
      do i=2,Nx-1
      do j=2,Ny-1
!****************************** 
     Ape(i,j) = 2.0d0/(de(i)*(de(i)+dw(i)))
     Apw(i,j) = 2.0d0/(dw(i)*(de(i)+dw(i)))
     Apn(i,j) = 2.0d0/(dn(j)*(dn(j)+ds(j)))
     Aps(i,j) = 2.0d0/(ds(j)*(dn(j)+ds(j)))
!     App(i,j) = -1.0d0*(Ape(i,j)+Apw(i,j)+Apn(i,j)+Aps(i,j))
     App(i,j) = (-2.0d0/(de(i)*dw(i)))+(-2.0d0/(dn(j)*ds(j))) 
!   Ape(i,j) = etaE_2cen(i)
!   Apw(i,j) = etaW_2cen(i)
!   Apn(i,j) = etaN_2cen(j)
!   Aps(i,j) = etaS_2cen(j)
!   App(i,j) = -1.0d0*( etaE_2cen(i)+etaW_2cen(i)+etaN_2cen(j)+etaS_2cen(j) ) 
!****************************** 
       end do
       end do
!***debug-statement
!     ae = 2.0d0/(de(2)*(de(2)+dw(2)))
!     aw = 2.0d0/(dw(2)*(de(2)+dw(2)))
!     an = 2.0d0/(dn(2)*(dn(2)+ds(2)))
!     as = 2.0d0/(ds(2)*(dn(2)+ds(2)))
!     ap = -1.0d0*(ae+aw+an+as)
!  print*, 'diffusion-Stencil-E',Ape(2,2),ae
!  print*, 'diffusion-Stencil-W',Apw(2,2),aw
!  print*, 'diffusion-Stencil-N',Apn(2,2),an
!  print*, 'diffusion-Stencil-S',Aps(2,2),as
!  print*, 'diffusion-Stencil-P',App(2,2),ap
!pause
!****************************** 
        call initialize()
!        call savesnapshot()
!***debug statement
! do i=1,Nx
! do j=1,Ny

! if(Ut(i,j,1).ne.0.0d0) print*, 'U',i,j,Ut(i,j,1)
! if(Ut(i,j,2).ne.0.0d0) print*, 'V',i,j,Ut(i,j,2)
! if(pcorr(i,j).ne.0.0d0) print*, 'pcorr',i,j,pcorr(i,j)
! if(press(i,j).ne.1.0d0) print*, 'press',i,j,press(i,j)
 
! end do
! end do

!pause
!**********
   time=0.0d0
  nstep=0

  do while(nstep.lt.maxstep)

    nstep=nstep+1
    time=time+dt

        call convective()
!***debug statement
! do i=1,Nx
! do j=1,Ny
! if(xcell_Pe.ne.0.0d0) print*, i,j,xcell_Pe
! if(ycell_Pe.ne.0.0d0) print*, i,j,ycell_Pe
! end do
! end do
! print*, 'conv'
!pause
!**********

        call solve_transport()

!***debug statement
! print*, 'transport'
!    if(mod(nstep,snapshot).eq.0) call savesnapshot()
! pause
!**********

        call rhei_chow()

!***debug statement
! print*, 'rhie-chow'
! pause
!**********

        call pressure_poisson()

!***debug statement
! print*, 'pressure-poisson'
! pause
!**********

        call update_primitives()

!***debug statement
! print*, 'primitive'
!    if(mod(nstep,snapshot).eq.0) call savesnapshot()
! pause
!**********

        call boundary_condition()

        call divergence()

    if(mod(nstep,snapshot)==0)then
    call savesnapshot()
    end if

  end do
        
!********************************************************
!*****finally the clean-up act

   deallocate(x)
   deallocate(y)
   deallocate(Ut)
   deallocate(press)
   deallocate(pcorr)
   deallocate(Q)
   deallocate(conv)
   deallocate(Qp)

!********************************************************
   contains
!********************************************************
!
!
!
!********************************************************
    subroutine initialize()
!************************************
    implicit none 
!************************************
         do i=1,Nx
         do j=1,Ny

   Ut(i,j,1) = dsin(x(i))*dcos(y(j))!*dexp(-2.0d0*dt/Re)
    Ut(i,j,2) = -1.0d0*dsin(y(j))*dcos(x(i))!*dexp(-2.0d0*dt/Re)
    Ut(i,j,3) = dsin(x(i))*dcos(y(j))!*dexp(-2.0d0*dt/Re)
    Ut(i,j,4) = -1.0d0*dsin(y(j))*dcos(x(i))!*dexp(-2.0d0*dt/Re)
   press(i,j) = 0.0d0 !(dcos(2.0d0*x(i))+dcos(2.0d0*y(j)))/4.0d0
   pcorr(i,j) = 0.0d0 !(dcos(2.0d0*x(i))+dcos(2.0d0*y(j)))/4.0d0
  forc(i,j,1) = 0.0d0
  forc(i,j,2) = 0.0d0
 
         end do
         end do
!********************************
!***at the top-layer


!************************************
    end subroutine initialize
!********************************************************
!
!
!
!********************************************************
     subroutine convective()
!************************************
    implicit none 
!************************************
 double precision :: dx_int,dy_int,x_conv,y_conv,xintPLUS,xintMINUS,yintPLUS,yintMINUS
!************************************

       do i=2,Nx-1
       do j=2,Ny-1
       do jcnt=1,2

!************************************   
   xintPLUS = ( x(i+1)+x(i) )/2.0d0
  xintMINUS = ( x(i)+x(i-1) )/2.0d0

   yintPLUS = ( y(j+1)+y(j) )/2.0d0
  yintMINUS = ( y(j)+y(j-1) )/2.0d0
!************************************
     dx_int = xintPLUS-xintMINUS
     dy_int = yintPLUS-yintMINUS
!************************************
    xcell_Pe = dabs(Re*Ut(i,j,1)*dx_int)
    ycell_Pe = dabs(Re*Ut(i,j,2)*dy_int)
!*****************************
!***along x-direction
!*****************************
!   if(i.ge.4.and.i.le.(Nx-3))then       !*****###*****!
!*****************************
   if(i.le.2.or.i.ge.(Nx-1))then       !*****###*****!
 x_conv=Ut(i,j,1)*( ( etaE_gradP(i)*(Ut(i+1,j,jcnt)-Ut(i,j,jcnt)) )+( etaW_gradP(i)*(Ut(i-1,j,jcnt)-Ut(i,j,jcnt)) ) )
   elseif(xcell_Pe.ge.2.0d0.and.Ut(i,j,1).gt.0.0d0)then
 x_conv = Ut(i,j,1)*( ( etaEE_upw_pos(i)*(Ut(i+2,j,jcnt)-Ut(i,j,jcnt)) )+( etaE_upw_pos(i)*(Ut(i+1,j,jcnt)-Ut(i,j,jcnt)) )+( etaW_upw_pos(i)*(Ut(i-1,j,jcnt)-Ut(i,j,jcnt)) )+( etaWW_upw_pos(i)*(Ut(i-2,j,jcnt)-Ut(i,j,jcnt)) ) )
   elseif(xcell_Pe.ge.2.0d0.and.Ut(i,j,1).lt.0.0d0)then
   x_conv = Ut(i,j,1)*( ( etaEE_upw_neg(i)*(Ut(i+2,j,jcnt)-Ut(i,j,jcnt)) )+( etaE_upw_neg(i)*(Ut(i+1,j,jcnt)-Ut(i,j,jcnt)) )+( etaW_upw_neg(i)*(Ut(i-1,j,jcnt)-Ut(i,j,jcnt)) )+( etaWW_upw_neg(i)*(Ut(i-2,j,jcnt)-Ut(i,j,jcnt)) ) )
   else
   x_conv = Ut(i,j,1)*( ( etaEE_1cen(i)*(Ut(i+2,j,jcnt)-Ut(i,j,jcnt)) )+( etaE_1cen(i)*(Ut(i+1,j,jcnt)-Ut(i,j,jcnt)) )+( etaW_1cen(i)*(Ut(i-1,j,jcnt)-Ut(i,j,jcnt)) )+( etaWW_1cen(i)*(Ut(i-2,j,jcnt)-Ut(i,j,jcnt)) ) )
   end if
!*****************************
!   else                                 !*****###*****!
!*****************************
!   if(xcell_Pe.ge.2.0d0.and.Ut(i,j,jcnt).gt.0.0d0)then
!   x_conv = Ut(i,j,1)*( etaE_upw_pos(i)*(Ut(i+1,j,jcnt)-Ut(i,j,jcnt))+etaW_upw_pos(i)*(Ut(i-1,j,jcnt)-Ut(i,j,jcnt)) )
!   elseif(xcell_Pe.ge.2.0d0.and.Ut(i,j,jcnt).lt.0.0d0)then
!   x_conv = Ut(i,j,1)*( etaE_upw_neg(i)*(Ut(i+1,j,jcnt)-Ut(i,j,jcnt))+etaW_upw_neg(i)*(Ut(i-1,j,jcnt)-Ut(i,j,jcnt)) )
!   else
!   x_conv = Ut(i,j,1)*( etaE_1cen(i)*(Ut(i+1,j,jcnt)-Ut(i,j,jcnt))+etaW_1cen(i)*(Ut(i-1,j,jcnt)-Ut(i,j,jcnt)) )
!   end if
!*****************************
!   end if                               !*****###*****!
!*****************************
!***along y-direction
!*****************************
!   if(j.ge.4.and.j.le.(Ny-3))then       !*****###*****!
!*****************************
   if(j.le.2.or.j.ge.(Ny-1))then       !*****###*****!
!   y_conv = ( ( etaN_1cen(j)*(Ut(i,j+1,jcnt)-Ut(i,j,jcnt)) )+( etaS_1cen(j)*(Ut(i,j-1,jcnt)-Ut(i,j,jcnt)) ) )
   y_conv = Ut(i,j,2)*( ( etaN_gradP(j)*(Ut(i,j+1,jcnt)-Ut(i,j,jcnt)) )+( etaS_gradP(j)*(Ut(i,j-1,jcnt)-Ut(i,j,jcnt)) ) )
   elseif(ycell_Pe.ge.2.0d0.and.Ut(i,j,2).gt.0.0d0)then
   y_conv = Ut(i,j,2)*( ( etaNN_upw_pos(j)*(Ut(i,j+2,jcnt)-Ut(i,j,jcnt)) )+( etaN_upw_pos(j)*(Ut(i,j+1,jcnt)-Ut(i,j,jcnt)) )+( etaS_upw_pos(j)*(Ut(i,j-1,jcnt)-Ut(i,j,jcnt)) )+( etaSS_upw_pos(j)*(Ut(i,j-2,jcnt)-Ut(i,j,jcnt)) ) )
   elseif(ycell_Pe.ge.2.0d0.and.Ut(i,j,2).lt.0.0d0)then
   y_conv = Ut(i,j,2)*( ( etaNN_upw_neg(j)*(Ut(i,j+2,jcnt)-Ut(i,j,jcnt)) )+( etaN_upw_neg(j)*(Ut(i,j+1,jcnt)-Ut(i,j,jcnt)) )+( etaS_upw_neg(j)*(Ut(i,j-1,jcnt)-Ut(i,j,jcnt)) )+( etaSS_upw_neg(j)*(Ut(i,j-2,jcnt)-Ut(i,j,jcnt)) ) )
   else
   y_conv = Ut(i,j,2)*( ( etaNN_1cen(j)*(Ut(i,j+2,jcnt)-Ut(i,j,jcnt)) )+( etaN_1cen(j)*(Ut(i,j+1,jcnt)-Ut(i,j,jcnt)) )+( etaS_1cen(j)*(Ut(i,j-1,jcnt)-Ut(i,j,jcnt)) )+( etaSS_1cen(j)*(Ut(i,j-2,jcnt)-Ut(i,j,jcnt)) ) )
   end if
!*****************************
!   else                                 !*****###*****!
!*****************************
!   if(ycell_Pe.ge.2.0d0.and.Ut(i,j,jcnt).gt.0.0d0)then
!   y_conv = Ut(i,j,2)*( etaN_upw_pos(j)*(Ut(i,j+1,jcnt)-Ut(i,j,jcnt))+etaS_upw_pos(j)*(Ut(i,j-1,jcnt)-Ut(i,j,jcnt)) )
!   elseif(ycell_Pe.ge.2.0d0.and.Ut(i,j,jcnt).lt.0.0d0)then
!   y_conv = Ut(i,j,2)*( etaN_upw_neg(j)*(Ut(i,j+1,jcnt)-Ut(i,j,jcnt))+etaS_upw_neg(j)*(Ut(i,j-1,jcnt)-Ut(i,j,jcnt)) )
!   else
!   y_conv = Ut(i,j,2)*( etaN_1cen(j)*(Ut(i,j+1,jcnt)-Ut(i,j,jcnt))+etaS_1cen(j)*(Ut(i,j-1,jcnt)-Ut(i,j,jcnt)) )
!   end if
!*****************************
!   end if                               !*****###*****!
!*****************************
!***debug statement
! if(xcell_Pe.ne.0.0d0) print*, i,j,xcell_Pe
! if(ycell_Pe.ne.0.0d0) print*, i,j,ycell_Pe
!  if(x_conv.ne.0.0d0) print*, i,j,jcnt,x_conv
!  if(y_conv.ne.0.0d0) print*, i,j,jcnt,y_conv,etaN_1cen(j)
!**********
!***obtain the final-convective term

   conv(i,j,jcnt) = x_conv+y_conv 

!***debug statement
!  if(conv(i,j,jcnt).ne.0.0d0) print*, i,j,jcnt,conv(i,j,jcnt)
!*****************************
       end do
       end do
       end do
!*****************************

       do i=2,Nx-1
       do j=2,Ny-1
!*****************************

    gradP(1) = (etaE_gradP(i)*(press(i+1,j)-press(i,j)))+(etaW_gradP(i)*(press(i-1,j)-press(i,j))) 
    gradP(2) = (etaN_gradP(j)*(press(i,j+1)-press(i,j)))+(etaS_gradP(j)*(press(i,j-1)-press(i,j))) 
!***debug statement
!  if(gradP(1).ne.0.0d0) print*, i,j,gradP(1)
!  if(gradP(2).ne.0.0d0) print*, i,j,gradP(2)

!   Ut(i,j,3) = Ut(i,j,1)+dt*(-conv(i,j,1)+( ( (Ave(i,j)*Ut(i+1,j,1))+(Avw(i,j)*Ut(i-1,j,1))+(Avn(i,j)*Ut(i,j+1,1))+(Avs(i,j)*Ut(i,j-1,1))+(Avp(i,j)*Ut(i,j,1)) )/Re ) ) 
!   Ut(i,j,1) = Ut(i,j,3)+(dt*-gradP(1))

!   Ut(i,j,4) = Ut(i,j,2)+dt*(-conv(i,j,2)+( ( (Ave(i,j)*Ut(i+1,j,2))+(Avw(i,j)*Ut(i-1,j,2))+(Avn(i,j)*Ut(i,j+1,2))+(Avs(i,j)*Ut(i,j-1,2))+(Avp(i,j)*Ut(i,j,2)) )/Re ) ) 
!   Ut(i,j,2) = Ut(i,j,4)+(dt*-gradP(2))

    Q(i,j,3) = Ut(i,j,1)+dt*(-conv(i,j,1)) 
    Q(i,j,1) = Q(i,j,3)+(dt*-gradP(1))

    Q(i,j,4) = Ut(i,j,2)+dt*(-conv(i,j,2)) 
    Q(i,j,2) = Q(i,j,4)+(dt*-gradP(2))
!***debug statement
!  if(Q(i,j,1).ne.0.0d0) print*, i,j,1,Q(i,j,1)
!  if(Q(i,j,2).ne.0.0d0) print*, i,j,2,Q(i,j,2)
!  if(Q(i,j,3).ne.0.0d0) print*, i,j,3,Q(i,j,3)
!  if(Q(i,j,4).ne.0.0d0) print*, i,j,4,Q(i,j,4)
!*****************************
       end do
       end do
!************************************
    end subroutine convective
!********************************************************
!
!
!
!********************************************************
     subroutine solve_transport()
!************************************
    implicit none 
!************************************
 integer :: iter
 double precision :: sumav,sumnor,tmp,res
!************************************
!    iter=0
    sumav=0.0d0
!***debug-statement
!  print*, 'diffusion-Stencil-E',Ave(2,2)
!  print*, 'diffusion-Stencil-W',Avw(2,2)
!  print*, 'diffusion-Stencil-N',Avn(2,2)
!  print*, 'diffusion-Stencil-S',Avs(2,2)
!  print*, 'diffusion-Stencil-P',Avp(2,2)
!************************************
  do jcnt=1,4
!************************************
!   do while(sumav.gt.errtol)
    do iter=1,maxiter
!    do iter=1,9
!************************************
!    iter=iter+1
    tmp=0.0d0

   do i=2,Nx-1
   do j=2,Ny-1
!*****************************
!   if(i.ge.4.and.j.ge.4.and.i.le.(Nx-3).and.j.le.(Ny-3))then       !*****###*****!
!*****************************
!   Ut(i,j,jcnt)=(Q(i,j,jcnt)-( (Avee(i+2,j)*Ut(i+2,j,jcnt))+(Avww(i,j)*Ut(i-2,j,jcnt))+(Avnn(i,j)*Ut(i,j+2,jcnt))+(Avss(i,j)*Ut(i,j-2,jcnt))+(Ave(i,j)*Ut(i+1,j,jcnt))+(Avw(i,j)*Ut(i-1,j,jcnt))+(Avn(i,j)*Ut(i,j-1,jcnt))+(Avs(i,j)*Ut(i,j+1,jcnt)) ))/Avp(i,j)
!*****************************
!   else                                 !*****###*****!
!*****************************
   Ut(i,j,jcnt)=( Q(i,j,jcnt)-( (Ave(i,j)*Ut(i+1,j,jcnt))+(Avw(i,j)*Ut(i-1,j,jcnt))+(Avn(i,j)*Ut(i,j+1,jcnt))+(Avs(i,j)*Ut(i,j-1,jcnt)) ) )/Avp(i,j)
!*****************************
!***debug statement
!   if(Ut(50,j,1).ne.0.0d0) print*, 'U',iter,j,jcnt,Ut(50,j,1) 
!   if(Ut(i,j,2).ne.0.0d0) print*, 'V',i,j,jcnt,Ut(i,j,2) 
!   if(Ut(i,j,3).ne.0.0d0) print*, 'Up',i,j,jcnt,Ut(i,j,3) 
!   if(Ut(i,j,4).ne.0.0d0) print*, 'Vp',i,j,jcnt,Ut(i,j,4) 
!*****************************
!   end if                               !*****###*****!
!*****************************
   end do
   end do

   do i=2,Nx-1
   do j=2,Ny-1
!*****************************
!   if(i.ge.4.and.j.ge.4.and.i.le.(Nx-3).and.j.le.(Ny-3))then       !*****###*****!
!*****************************
!   res = Ut(i,j,jcnt)-(Q(i,j,jcnt)-( (Avee(i+2,j)*Ut(i+2,j,jcnt))+(Avww(i,j)*Ut(i-2,j,jcnt))+(Avnn(i,j)*Ut(i,j+2,jcnt))+(Avss(i,j)*Ut(i,j-2,jcnt))+(Ave(i,j)*Ut(i+1,j,jcnt))+(Avw(i,j)*Ut(i-1,j,jcnt))+(Avn(i,j)*Ut(i,j-1,jcnt))+(Avs(i,j)*Ut(i,j+1,jcnt)) ))/Avp(i,j)
!*****************************
!   else                                 !*****###*****!
!*****************************
   res = Ut(i,j,jcnt)-((Q(i,j,jcnt)-( (Ave(i,j)*Ut(i+1,j,jcnt))+(Avw(i,j)*Ut(i-1,j,jcnt))+(Avn(i,j)*Ut(i,j+1,jcnt))+(Avs(i,j)*Ut(i,j-1,jcnt)) ))/Avp(i,j))
!*****************************
!   end if                               !*****###*****!
!*****************************
    tmp = tmp+dabs(res)
! print*, iter,sumav
   end do
   end do

  if(tmp.eq.0) goto 111
  if(iter.eq.1) sumnor=tmp
!  if(sumnor.ne.0.0d0) sumav=tmp/sumnor
  sumav=tmp/sumnor

  if(sumav.le.errtol)goto 111
!************************************
   end do
!************************************
      print *,sumav,iter,jcnt
111    continue
!      print *,sumav,iter,jcnt
  end do 
!************************************
    end subroutine solve_transport
!********************************************************
!
!
!
!********************************************************
     subroutine rhei_chow()
!************************************
    implicit none 
!************************************
    double precision, dimension(2) :: divVELstar
    double precision :: Ustar_xint_for,Ustar_xint_bck,Vstar_yint_for,Vstar_yint_bck
!************************************
!*****now we have got up,vp,wp and provisional u,v,w********************                   
!****************************************************************************************
!now applying Rhie-Chow Momentum Interpolation to avoid spurious oscillations in pressure
!****************************************************************************************
!*****boundary conditions for the w/o pressure field velocity (up,vp)*****************
         Ut(1:Nx,1,3) = Ut(1:Nx,Ny-1,3)    ! predicted up is zero at the bottom wall. 
         Ut(1:Nx,1,4) = Ut(1:Nx,Ny-1,4)    ! predicted up is zero at the bottom wall. 

        Ut(1:Nx,Ny,3) = Ut(1:Nx,2,3)    ! predicted up is unity at the top wall. 
        Ut(1:Nx,Ny,4) = Ut(1:Nx,2,4)    ! predicted up is unity at the top wall. 

        Ut(1,2:Ny-1,3)= Ut(Nx-1,2:Ny-1,3)    ! predicted up is zero at the back wall. 
        Ut(1,2:Ny-1,4)= Ut(Nx-1,2:Ny-1,4)    ! predicted up is zero at the back wall. 

       Ut(Nx,2:Ny-1,3)= Ut(2,2:Ny-1,3)    ! predicted up is zero at the front wall. 
       Ut(Nx,2:Ny-1,4)= Ut(2,2:Ny-1,4)    ! predicted up is zero at the front wall. 
!****************************************************************************************    
!***debug statement
!   print*, de(98),dw(98)
!***********************
         do i=2,Nx-1
         do j=2,Ny-1 


   Ustar_xint_for = ( 0.5d0*(Ut(i+1,j,3)+Ut(i,j,3)) )-( dt*(press(i+1,j)-press(i,j))/de(i) ) 
   Ustar_xint_bck = ( 0.5d0*(Ut(i,j,3)+Ut(i-1,j,3)) )-( dt*(press(i,j)-press(i-1,j))/dw(i) ) 


   Vstar_yint_for = ( 0.5d0*(Ut(i,j+1,4)+Ut(i,j,4)) )-( dt*(press(i,j+1)-press(i,j))/dn(j) ) 
   Vstar_yint_bck = ( 0.5d0*(Ut(i,j,4)+Ut(i,j-1,4)) )-( dt*(press(i,j)-press(i,j-1))/ds(j) ) 

    divVELstar(1) =  2.0d0*(Ustar_xint_for-Ustar_xint_bck)/(de(i)+dw(i)) 
    divVELstar(2) =  2.0d0*(Vstar_yint_for-Vstar_yint_bck)/(dn(j)+ds(j)) 

          Qp(i,j) = ( divVELstar(1)+divVELstar(2) )/dt
!***debug-statement
!  if(i.eq.5) print*, de(5),dw(5)
!  if(j.eq.5) print*, dn(5),ds(5)
!   if(Ustar_xint_for.ne.0.0d0) print*, 'Ui+1/2',i,j,Ustar_xint_for  
!   if(Ustar_xint_bck.ne.0.0d0) print*, 'Ui-1/2',i,j,Ustar_xint_bck  
!   if(Vstar_yint_for.ne.0.0d0) print*, 'Vj+1/2',i,j,Vstar_yint_for  
!   if(Vstar_yint_bck.ne.0.0d0) print*, 'Vj-1/2',i,j,Vstar_yint_bck  
!   if(Qp(i,j).ne.0.0d0) print*, 'Qp',i,j,Qp(i,j)
!    if(i.eq.98.and.j.eq.98)then
!    print*, Ustar_xint_for,Ustar_xint_bck,Vstar_yint_for,Vstar_yint_bck 
!    print*, divVELstar(1),divVELstar(2),Qp(i,j)
!    end if
!****************************************************************************************    
         end do
         end do
!************************************
    end subroutine rhei_chow
!********************************************************
!
!
!
!********************************************************
     subroutine pressure_poisson()
!************************************
    implicit none 
!************************************
    double precision :: sumav,tmp,sumnor,res
    integer :: iter
!************************************
!*****initial condition for correction pressure*****************
         pcorr(1:Nx,1) = 0.0d0    ! pcorr is zero at the bottom wall. 
         pcorr(1:Nx,1) = 0.0d0    ! pcorr is zero at the bottom wall. 

        pcorr(1:Nx,Ny) = 0.0d0    ! pcorr is zero at the top wall. 
        pcorr(1:Nx,Ny) = 0.0d0    ! pcorr is zero at the top wall. 

        pcorr(1,2:Ny-1)= 0.0d0    ! pcorr is zero at the back wall. 
        pcorr(1,2:Ny-1)= 0.0d0    ! pcorr is zero at the back wall. 

       pcorr(Nx,2:Ny-1)= 0.0d0    ! pcorr is zero at the front wall. 
       pcorr(Nx,2:Ny-1)= 0.0d0    ! pcorr is zero at the front wall. 
!************************************
   iter=0
   sumav=0.0d0
!************************************

!   do while(sumav.gt.errtol)
!    iter=iter+1
   do iter = 1,maxiter
    tmp=0.0d0

   do i=2,Nx-1
   do j=2,Ny-1
    res = pcorr(i,j)-( (Qp(i,j)-( (Ape(i,j)*pcorr(i+1,j))+(Apw(i,j)*pcorr(i-1,j))+(Apn(i,j)*pcorr(i,j+1))+(Aps(i,j)*pcorr(i,j-1)) ))/App(i,j) )
    tmp = tmp+dabs(res)
   end do
   end do

  if(tmp.eq.0)then
  goto 191
  end if
  if(iter.eq.1)then
  sumnor=tmp
  end if
!  if(sumnor.ne.0.0d0) sumav=tmp/sumnor
  sumav=tmp/sumnor

   do i=2,Nx-1
   do j=2,Ny-1
   pcorr(i,j)=(Qp(i,j)-( (Ape(i,j)*pcorr(i+1,j))+(Apw(i,j)*pcorr(i-1,j))+(Apn(i,j)*pcorr(i,j+1))+(Aps(i,j)*pcorr(i,j-1)) ))/App(i,j)
   end do
   end do

  if(sumav.lt.errtol)then
  goto 191
  end if

       end do

 print*, 'P-P',iter,sumav
191    continue
! print*, 'P-P',iter,sumav
! print*, pcorr(2,5),pcorr(3,5)
! print*, pcorr(2,87),pcorr(3,87)
! print*, pcorr(99,5),pcorr(98,5)
! print*, pcorr(99,87),pcorr(98,87)
!****************************************

!************************************
    end subroutine pressure_poisson
!************************************
!
!
!
!********************************************************
     subroutine update_primitives()
!************************************
    implicit none 
!************************************
    double precision :: ucorr,vcorr
    double precision, dimension(2) :: gradPcorr
!************************************
       do i=2,Nx-1
       do j=2,Ny-1
!*****************************
      press(i,j) = press(i,j)+pcorr(i,j)

    gradPcorr(1) = (etaE_gradP(i)*(pcorr(i+1,j)-pcorr(i,j)))+(etaW_gradP(i)*(pcorr(i-1,j)-pcorr(i,j))) 
    gradPcorr(2) = (etaN_gradP(j)*(pcorr(i,j+1)-pcorr(i,j)))+(etaS_gradP(j)*(pcorr(i,j-1)-pcorr(i,j)))

           ucorr = -dt*gradPcorr(1)
           vcorr = -dt*gradPcorr(2)

       Ut(i,j,1) = Ut(i,j,1)+ucorr
       Ut(i,j,2) = Ut(i,j,2)+vcorr
!*****************************
         end do
         end do
!************************************
    end subroutine update_primitives
!************************************
!
!
!
!************************************
     subroutine boundary_condition()
!************************************
    implicit none 
!************************************
    double precision :: del2U_delx2,del2U_dely2
    double precision :: del2V_delx2,del2V_dely2
    double precision :: sumtmp
    double precision :: conv_v_1,conv_v_2
!************************************
          Ut(1:Nx,1,1) = Ut(1:Nx,Ny-1,1)    ! predicted up is zero at the bottom wall. 
         Ut(1:Nx,1,2) = Ut(1:Nx,Ny-1,2)    ! predicted up is zero at the bottom wall. 

        Ut(1:Nx,Ny,1) = Ut(1:Nx,2,1)    ! predicted up is unity at the top wall. 
        Ut(1:Nx,Ny,2) = Ut(1:Nx,2,2)    ! predicted up is unity at the top wall. 

        Ut(1,2:Ny-1,1)= Ut(Nx-1,2:Ny-1,1)    ! predicted up is zero at the back wall. 
        Ut(1,2:Ny-1,2)= Ut(Nx-1,2:Ny-1,2)    ! predicted up is zero at the back wall. 

       Ut(Nx,2:Ny-1,1)= Ut(2,2:Ny-1,1)    ! predicted up is zero at the front wall. 
       Ut(Nx,2:Ny-1,2)= Ut(2,2:Ny-1,2)    ! predicted up is zero at the front wall. 

       press(1,2:Ny-1) = press(Nx-1,2:Ny-1)    ! predicted up is zero at the back wall. 
       press(Nx,2:Ny-1) = press(2,2:Ny-1)    ! predicted up is zero at the front wall. 
       press(1:Nx,1) = press(1:Nx,Ny-1)    ! predicted up is zero at the bottom wall. 
       press(1:Nx,Ny) = press(1:Nx,2)    ! predicted up is unity at the top wall. 
!************************************
    end subroutine boundary_condition
!************************************
!
!
!
!************************************
     subroutine divergence()
!************************************
    implicit none 
!************************************
    double precision :: div_x,div_y,summ,div
!************************************

       summ = 0.0d0

       do i=2,Nx-1
       do j=2,Ny-1
!*****************************
     div_x = (etaE_gradP(i)*(Ut(i+1,j,1)-Ut(i,j,1)))+(etaW_gradP(i)*(Ut(i-1,j,1)-Ut(i,j,1))) 
     div_y = (etaN_gradP(j)*(Ut(i,j+1,2)-Ut(i,j,2)))+(etaS_gradP(j)*(Ut(i,j-1,2)-Ut(i,j,2)))

     summ = summ+div_x+div_y
!*****************************
         end do
         end do

      div = dsqrt( dabs(summ)/dble((Nx-2)*(Ny-2)) )

      print*, nstep,time,div
!************************************
    end subroutine divergence
!************************************
!
!
!
!************************************
     subroutine savesnapshot()
!************************************
    implicit none 
!************************************
!*****write the flow-field data in ASCII-Tecplot format directly*****
 write(filname,'(a,i7.7,a)')'2D',nstep/snapshot,'.tp' 
 open(unit=7,file=filname,status='unknown')
            
 write(7,13) 'TITLE ="',time,nstep,'"'
 write(7,*) 'VARIABLES = "x","y","U","V","Up","Vp","press","pcorr"'
 write(7,14) 'ZONE j=',Nx,',i=',Ny,',DATAPACKING ="POINT"'
! write(7,14) 'ZONE i=',Ny,',j=',Nx,',DATAPACKING ="POINT"'

    do i=1,Nx
    write(7,*)
           do j=1,Ny
!***debug statements
!    do j=1,Ny
!    write(7,*)
!           do i=1,Nx

!  write(7,59) i,j,Ut(i,j,1),Ut(i,j,2),Ut(i,j,3),Ut(i,j,4),press(i,j),pcorr(i,j)
!*******************
  write(7,59) x(i),y(j),Ut(i,j,1),Ut(i,j,2),Ut(i,j,3),Ut(i,j,4),press(i,j),pcorr(i,j)
           end do
    end do

   close(7)

13 format(A,F10.6,I8,A)
14 format(A,I3,A,I3,A)
59 format(2(1X,F10.7),6(1X,E17.7)) 
!59 format(2(1X,I4),6(1X,E15.7)) 
!************************************
    end subroutine savesnapshot
!************************************
!
!
!
!********************************************************
  end program incompressible
!********end of program**********************************
