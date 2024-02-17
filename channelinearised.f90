program channelinearised
! Sat May 11 10:31:03 IST 2013
! This is the main Program 
! based on original source-code developed and written by Fuaad P A
! Under the Guidance of Prof Mirza Faisal Baig  
! Aligarh Muslim University
! 
! for integer parametric array zv the elemets change with value of delta
!
implicit none
!
! number types
!
INTEGER, PARAMETER :: nstart = 50
! parameters
! data read in from files
!
!
LOGICAL ::  binary       =.false.           ! logical operator for backup file.
LOGICAL ::  control      =.false.          ! logical operator to switch the control on & off.
! if FALSE then control is off.
INTEGER ::  nstep                           ! current number of time integrations 
INTEGER ::  maxstep      =20000             ! maximum number of time integrations 
INTEGER ::  maxiter      =1000               ! maximum number of iterations
DOUBLE PRECISION :: errorstop =1.0d-8      ! Convergence Criteria of Transport equations
INTEGER ::  disp         =1                 ! display frequency 0f divergence 
INTEGER ::  bckp         =500              ! backup file display frequency.
INTEGER ::  snapshot     =1000              ! display frequency 
INTEGER ::  icase        =1                 ! an integer to select type of file for restart.
DOUBLE PRECISION :: time                          ! current non-dimensional time
DOUBLE PRECISION , PARAMETER :: dt=1.0d-4        ! timestep
DOUBLE PRECISION :: dXi,dx,dy          
DOUBLE PRECISION :: KEmax,Pi

INTEGER :: i,j,k,jcnt,it                            ! Local loop variables
INTEGER :: im1,im2,im3,ip1,ip2,ip3,jm1,jm2,jm3,jp1,jp2,jp3  ! loop variables
! Real Thing the Flow field variables
!
DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: Ut,Q,conv,linforc,Um,gradp       ! 3d Flow Field (t:i:j:k)
DOUBLE PRECISION, dimension(:,:,:), allocatable   :: Qp,Pcorr,Press                   ! 3d terms for pressure
DOUBLE PRECISION, dimension(:,:,:,:), allocatable   :: Utnorm, Utsq              ! 3d terms for the measure and quantification of the response.
!
! gmres internal Parameters
!
DOUBLE PRECISION,allocatable,dimension(:,:,:,:):: v
DOUBLE PRECISION,allocatable,dimension(:,:,:)  :: r,rhat
DOUBLE PRECISION,allocatable,dimension(:,:)    :: h
DOUBLE PRECISION,allocatable,dimension(:)      :: g,c,s,y,mag
DOUBLE PRECISION                               :: rho,hik,hipk,nu,gk,gkp,w,errtol
DOUBLE PRECISION                               :: tempgmres
INTEGER                            :: kmax,kit,mmax,m,sgn
!
!real variables used in Gauss-Seidal
!
DOUBLE PRECISION :: res                                 !Hold the Residual rho
DOUBLE PRECISION :: residual, tmpsum ,t,sumnor,sumav    !Variables for Gauss-siedel
!
!
DOUBLE PRECISION, parameter     :: Re  =180.0d0    ! Friction-Reynolds Number 
!
! Mesh Parameters
!
INTEGER, parameter, dimension(3) :: n3d=(/131,67,131/)                                  ! grid points in 3 dimensions
DOUBLE PRECISION, dimension(3) :: Umax ,ke                                              ! Maximum Velocity in 3 dimensions
DOUBLE PRECISION, dimension(:,:,:), allocatable :: Anw,Ann,Ane,Ans,Anp,Ant,Anb,rAnp    ! Stencil Parameters for Pressure
DOUBLE PRECISION, dimension(:,:,:), allocatable :: Avw,Avn,Ave,Avs,Avp,Avt,Avb         ! Stencil Parameters for transport
DOUBLE PRECISION, dimension(:),allocatable :: x1,x2,x3,Jac,Jac2,Xi                     ! Physical Geometery (Cartesian)
!
!
DOUBLE PRECISION, parameter ::   delta =5.7d0            !   stretching parameter
DOUBLE PRECISION ::   beta,b,Gmag,a0,b0                   ! parametric values to be used in nonlinear source-terms
!
!
!
DOUBLE PRECISION :: Pe,term1,term2
DOUBLE PRECISION, dimension(3) :: divcon          ! variable used for computing components of continuity equation.    
DOUBLE PRECISION, dimension(3) :: divuu,Linear    ! variable used for computation of convective terms
DOUBLE PRECISION, dimension(3) :: divuStar        ! variable used for computation of gradient Ustar using Rhie-Chow momentum interpolation 
DOUBLE PRECISION, dimension(3) :: divgrad  ! variable used for computation of diffusion terms in control subroutine and final update of pressre.
DOUBLE PRECISION, dimension(3) :: gradpcorr   ! variable used for gradient of correction pressure.
DOUBLE PRECISION :: fder,Flux,KE1
INTEGER :: npln,kpln,zpln
DOUBLE PRECISION :: tmp,diver  
!
DOUBLE PRECISION, dimension(3,3) :: Uat12                             ! Grad of Ustar components
DOUBLE PRECISION,parameter,dimension(3) :: zp=(/6.0,12.0,18.0/)      ! zplus 
INTEGER, parameter, dimension(3) :: zv=(/27,32,35/)! zplus for planes.
!variables used in control-subroutine**************************************
DOUBLE PRECISION, parameter :: W0 =0.35
DOUBLE PRECISION, parameter :: Omega =30.0d0
DOUBLE PRECISION, parameter :: kappa =0.0d0
!************************************************************************************
DOUBLE PRECISION :: h1
   CHARACTER(LEN=25) :: filname
!*********************************Decleration of all CONSTANTS***********************  
               Pi = 4.0d0*datan(1.0d0)
               beta = 1.0*Re
               b =0.00325/Re
               Gmag =1.0d0
               a0 = 110.0d0
               b0 = 500.0d0


        errtol=1.0d-4
        kmax=50
        mmax=1000
!!        cnt=0
!!        h1=0.0d0
        w = 1.835
 
!*************************************************************************************
!**********************************Allocation of arrays*******************************
     allocate( Ut( 6,n3d(1), n3d(2), n3d(3)) )
     allocate( Um( 3,n3d(1), n3d(2), n3d(3)) )
     allocate( gradp( 3,n3d(1),n3d(2),n3d(3)) )
     allocate( linforc( 3,n3d(1), n3d(2), n3d(3)) )
     allocate( Press(n3d(1), n3d(2), n3d(3)) )

     allocate( x1(n3d(1)))
     allocate( x2(n3d(2))) 
     allocate( x3(n3d(3)))
     allocate( Xi(n3d(3)))

     allocate(Jac(n3d(3)) )
     allocate(Jac2(n3d(3)) )

     allocate( Conv( 3,n3d(1), n3d(2), n3d(3)) )
     allocate( Q ( 6,n3d(1), n3d(2), n3d(3)) )
     allocate( Pcorr(n3d(1), n3d(2), n3d(3)) )
     allocate(    Qp(n3d(1), n3d(2), n3d(3)) )

     allocate( Utnorm(3, n3d(1), n3d(2), n3d(3)) )
     allocate( Utsq(3, n3d(1), n3d(2), n3d(3)) )
    
     allocate(Anw(n3d(1), n3d(2), n3d(3)) )
     allocate(Ann(n3d(1), n3d(2), n3d(3)) )
     allocate(Ane(n3d(1), n3d(2), n3d(3)) )
     allocate(Ans(n3d(1), n3d(2), n3d(3)) )
     allocate(Anb(n3d(1), n3d(2), n3d(3)) )
     allocate(Ant(n3d(1), n3d(2), n3d(3)) )
     allocate(Anp(n3d(1), n3d(2), n3d(3)) )
     allocate( rAnp(n3d(1),n3d(2),n3d(3)) )

     allocate(Avw(n3d(1), n3d(2), n3d(3)) )
     allocate(Avn(n3d(1), n3d(2), n3d(3)) )
     allocate(Ave(n3d(1), n3d(2), n3d(3)) )
     allocate(Avs(n3d(1), n3d(2), n3d(3)) )
     allocate(Avb(n3d(1), n3d(2), n3d(3)) )
     allocate(Avt(n3d(1), n3d(2), n3d(3)) )
     allocate(Avp(n3d(1), n3d(2), n3d(3)) )
!These variables are internal to GMRES 

        allocate(v(nstart+1,n3d(1),n3d(2),n3d(3)))
        allocate(r(n3d(1),n3d(2),n3d(3)))
        allocate(rhat(n3d(1),n3d(2),n3d(3)))
        allocate(h(nstart+1,nstart+1))
        allocate(g(nstart+1))
        allocate(c(nstart+1))
        allocate(s(nstart+1))
        allocate(y(nstart+1))
        allocate(mag(nstart+1))

!Allocation ends here*****************************************************************************
!*************************************************************************************************
!  Initial guess used in GMRES.  
        rhat = 0.0d0
        h = 0.0d0
        v= 0.0d0              
        c= 0.0d0
        s= 0.0d0
        g = 0.0d0
        y = 0.0d0
!********************mesh-generation takes place here*********************************************
                x1(n3d(1)) = 4.0*Pi
                x2(n3d(2)) = 2.0*Pi
                x3(n3d(3)) = 2.0d0


      dx = x1(n3d(1))/dble(n3d(1)-1)
      dy = x2(n3d(2))/dble(n3d(2)-1)
      dXi= x3(n3d(3))/dble(n3d(3)-1)


    x1(1) =0.0
    x2(1) =0.0
    Xi(1) =0.0

 do i=1,n3d(1)-1
      x1(i+1) =x1(i) + dx    
 enddo

 do j=1,n3d(2)-1
      x2(j+1) =x2(j) + dy    
 enddo

 do k=1,n3d(3)-1
      Xi(k+1) =Xi(k) + dXi    
 enddo
!**************************************************************************************************
!Mapping of Non-Uniform Points from Uniform Points
!Give Symmetric spacing with respect to 1

        do k=1,n3d(3)
     x3(k) =1.0 + (dTanh(delta*0.50*(Xi(k)-1.0))/(dTanh(0.50*delta))) 
! Find the Local Derivative 
     Jac(k) = (0.5*delta)/(dTanh(0.5 *delta)*((dCosh(0.5 *delta *(Xi(k)-1.0)))**2))
! Find the Second Local Derivative 
     Jac2(k) = -(0.5*(delta**2)*dTanh(0.5*delta*(Xi(k)-1.0)))/((dTanh(0.5*delta)*(dCosh(0.5*delta*(Xi(k)-1.0)))**2))
       enddo
!********************saving the non-uniform points in a file****************************************
    open (unit=5,file='zpoints.dat',status='unknown')
    do k=1,n3d(3)
    write(5,23) k,x3(k),x3(k)*Re
    end do

    close(5)
23  format(1X,I4,2(1X,F20.10))
!****************************************************************************************************
! **************************initialization of the flow-field*****************************************
                 Um( 1:3,1:n3d(1),1:n3d(2),1:n3d(3) ) = 0.0d0
              
       open(unit=4,file='U_inner.dat')
                 do k=1,n3d(3)/2+1
       read(4,101) Um(1,1,1,k)
                 enddo
       close(4)

101    format(1x,E18.8)

! mirror the array 
      do k=n3d(3),n3d(3)/2+2,-1
          Um(1,1,1,k) =Um(1,1,1,n3d(3)-k+1) 
      enddo
! Copy the array to all points 

      do i=1,n3d(1)
              do j=1,n3d(2)
                      do k=1,n3d(3)
 Um(1,i,j,k) = Um(1,1,1,k) 
                      enddo
              enddo
      enddo

 ! Set up equations
        do k=1,n3d(3)
                do j=1,n3d(2)
                        do i=1,n3d(1)
      Anw(i,j,k)= 1.0d0/(dx*dx)
      Ann(i,j,k)= 1.0d0/(dy*dy)
      Ane(i,j,k)= 1.0d0/(dx*dx)
      Ans(i,j,k)= 1.0d0/(dy*dy)
      Anb(i,j,k)=(1.0d0/(dXi*dXi*Jac(k)*Jac(k)))+(Jac2(k)/(2.0*dXi*Jac(k)*Jac(k)*Jac(k)))
      Ant(i,j,k)=(1.0d0/(dXi*dXi*Jac(k)*Jac(k)))-(Jac2(k)/(2.0*dXi*Jac(k)*Jac(k)*Jac(k)))
      Anp(i,j,k)=-2.0d0*((1.0d0/(dx*dx)) +(1.0d0/(dy*dy))+(1.0/(dXi*dXi*Jac(k)*Jac(k))))
     rAnp(i,j,k)= 1.0d0/Anp(i,j,k)

                        enddo
                enddo
        enddo

    ! Set up equations
        do k=1,n3d(3)
                do j=1,n3d(2)
                        do i=1,n3d(1)
      Avw(i,j,k)=1.0d0/(dx*dx)
      Avn(i,j,k)=1.0d0/(dy*dy)
      Ave(i,j,k)=1.0d0/(dx*dx)
      Avs(i,j,k)=1.0d0/(dy*dy)
      Avb(i,j,k)=(1.0d0/(dXi*dXi*Jac(k)*Jac(k)))+(Jac2(k)/(2.0*dXi*Jac(k)*Jac(k)*Jac(k)))
      Avt(i,j,k)=(1.0d0/(dXi*dXi*Jac(k)*Jac(k)))-(Jac2(k)/(2.0*dXi*Jac(k)*Jac(k)*Jac(k)))
      Avp(i,j,k)=-2.0d0*((1.0d0/(dx*dx)) +(1.0d0/(dy*dy))+(1.0/(dXi*dXi*Jac(k)*Jac(k))))

                        enddo
                enddo
        enddo    

!************allocation of initial-values to the field variables**************************************
     if (binary.eqv..true.) then
!**************************************************
    print *,'binary is true'
! initialization of the arrays.
                 time=0.0d0
                 nstep=0
                 Ut(1,:,:,:) =0.0d0 
                 Ut(2,:,:,:) =0.0d0 
                 Ut(3,:,:,:) =0.0d0 
             Utnorm(:,:,:,:) =0.0d0
               Utsq(:,:,:,:) =0.0d0

  do i=1,n3d(1)
              do j=1,n3d(2)
                      do k=1,n3d(3)
! Initialize all the mean values to be zero 
   Press(i,j,k) = 10.0d0
                      enddo
              enddo
      enddo

!Apply the nonlinear forcing to the bottom plate of the channel.

        do k=1,n3d(3)
        do j=1,n3d(2)
        do i=1,n3d(1)
                 linforc(1,i,j,k)= 0.0d0
                 linforc(2,i,j,k)= 0.0d0
                 linforc(3,i,j,k)= Gmag*(x3(k)*x3(k)*cos(beta*x2(j))*exp((-a0*(x1(i)-1.014)**2)+(-b0*(x3(k)-0.025)**2))*exp(-1.0*b*x3(k)*x3(k)))
!    linforc(3,i,j,k)=Gmag*(cos(beta*x2(j))*exp((-a0*(x1(i)-1.014)**2)+(-b0*(x3(k)-0.025)**2)))                
        end do
        end do
        end do
!***************************************************
  else
!***************************************************
   if(icase.eq.1)then 
  print *,'binary is false and you might have faced a power failure.'       
     open (unit=10,FORM='unformatted',file='backup.xy',status='unknown')
      read(10) time,nstep
      read(10) Ut,Press
    close(10)
 else
 open (unit=12,file='snapat00010.tp',status='unknown')
 read (12,07) time,nstep
      read (12,*)
      read (12,*)
      do i=1,n3d(1)
        read(12,*)
              do j=1,n3d(2)
                read(12,*)
                      do k=1,n3d(3)
    read(12,87)x1(i),x2(j),x3(k),Ut(1,i,j,k),Ut(2,i,j,k),Ut(3,i,j,k),Press(i,j,k),Pcorr(i,j,k),gradp(3,i,j,k) 
                      end do 
              end do
      end do
      close(12)
!*********************************
07    format(A,F20.10,I9,A)
87    format(3(1X,F12.6),6(1X,E22.15))
   end if
!**********************************************************************
   end if
!**************************************************************************************
!************now starting the while-loop time-integrations***************************** 
do while( nstep < maxstep)

 time = time+dt

 if (time.gt.(1.0*dt)) then          ! forcing time is allocated here 
  linforc=0.0d0
  endif 
 
 nstep = nstep +1
!*********subroutine for implementation of control of streamwise travelling waves******** 
if(control.eq..true.)then
!**************************************************************************************** 
!at the Lower wall 
do i=1,n3d(1)
  Um(3,i,1:n3d(2),1) =W0*dcos(kappa*x1(i)-Omega*time)
enddo

!at the Upper wall 
  Um(3,1:n3d(1),1:n3d(2),n3d(3)) =0.0


            j=1

     do k=2,n3d(3)-1
       do i=2,n3d(1)-1

!***********************************************************************
                im1 =i-1
                im2 =i-2
                ip1= i+1
                ip2 =i+2
                im3 =i-3
                ip3 =i+3
!***********************************************************************

              Pe = abs(Re*Um(1,i,j,k)*(dx))   !taking only absolute quantity
       if(i.le.4.or.i.ge.(n3d(1)-3)) then
        divuu(1) = Um(1,i,j,k)*(Um(3,ip1,j,k)-Um(3,im1,j,k))/(2.0d0*dx)
       elseif(Pe.lt.2.0) then
 divuu(1)=Um(1,i,j,k)*((-Um(3,im3,j,k)+9.0*Um(3,im2,j,k)-45.0*Um(3,im1,j,k)+45.0*Um(3,ip1,j,k)-9.0*Um(3,ip2,j,k)+Um(3,ip3,j,k))/(60.0d0*dx))
         else
 term1 = (Abs(Um(1,i,j,k)))*(Um(3,ip2,j,k)-(4.0*(Um(3,ip1,j,k)))+(6.0*(Um(3,i,j,k)))-(4.0*(Um(3,im1,j,k)))+Um(3,im2,j,k))/(4.0d0*dx)    
 term2 =     (Um(1,i,j,k))*(-Um(3,ip2,j,k)+(8.0*(Um(3,ip1,j,k)))-(8.0*(Um(3,im1,j,k)))+Um(3,im2,j,k))/(12.0d0*dx)    
  divuu(1) = term1+term2
         endif

           Pe= abs((Re*Um(3,i,j,k)*dXi)/Jac(k))  !taking only absolute quantity
          if(k.le.4.or.k.ge.(n3d(3)-3)) then
  divuu(3) = Um(3,i,j,k)*(Um(3,i,j,k+1)-Um(3,i,j,k-1))/(2.0*dXi)
       elseif(Pe.lt.2.0) then
  divuu(3)=Um(3,i,j,k)*((-Um(3,i,j,k-3)+9*Um(3,i,j,k-2)-45*Um(3,i,j,k-1)+45*Um(3,i,j,k+1)-9*Um(3,i,j,k+2)+Um(3,i,j,k+3))/(60.0d0*dXi))
         else
term1 = (Abs(Um(3,i,j,k)))*(Um(3,i,j,k+2)-(4*(Um(3,i,j,k+1)))+(6*(Um(3,i,j,k)))-(4*(Um(3,i,j,k-1)))+Um(3,i,j,k-2))/(4.0d0*dXi)    
term2 =     (Um(3,i,j,k))*(-Um(3,i,j,k+2)+(8*(Um(3,i,j,k+1)))-(8*(Um(3,i,j,k-1)))+Um(3,i,j,k-2))/(12.0d0*dXi)    
  divuu(3) = term1+term2
         endif
      
  divgrad(1) =(Um(3,i+1,j,k)+Um(3,i-1,j,k)-2.0d0*Um(3,i,j,k))/(dx*dx)

  fder =(Um(3,i,j,k+1)-Um(3,i,j,k-1))/(2.0d0*dXi)

  divgrad(3) =(Um(3,i,j,k)+Um(3,i,j,k-1)-2.0d0*Um(3,i,j,k))/(dXi*dXi)

  divgrad(3)  = (divgrad(3)/(Jac(k)*Jac(k))) - ((Jac2(k)*fder)/(Jac(k)*Jac(k)*Jac(k)))

! Update Um(3) explicitly
! 
  Um(3,i,j,k) =Um(3,i,j,k) + dt*(-divuu(1)-(divuu(3)/Jac(k)) + (divgrad(1)+divgrad(3))/Re) 
!
       end do
       end do  

!Inflow boundary condition 
Um(3,1,j,1:n3d(3)) = 0.0d0
!Outflow boundary condition 
Um(3,n3d(1),j,1:n3d(3)) =(5.0*Um(3,n3d(1)-1,j,1:n3d(3))-4.0*Um(3,n3d(1)-2,j,1:n3d(3))+Um(3,n3d(1)-3,j,1:n3d(3)))/2.0d0    ! second derivative of Wmean zero at exit

       do j=2,n3d(2)
Um(3,1:n3d(1),j,1:n3d(3)) =Um(3,1:n3d(1),1,1:n3d(3))
       end do
!************************************************
end if
!****************************************************************************************                     
! Calculation of Convective Terms and the Linear terms
!****************************************************************************************
     do k=2,n3d(3)-1
      do j=2,n3d(2)-1
       do i=2,n3d(1)-1
        do jcnt=1,3
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
                im2 =im1-1
!              if (i.eq.3) im2 =n3d(1)-1
                im3 =im2-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
                jm2 =jm1-1
              if (j.eq.3) jm2 =n3d(2)-1
                jm3 =jm2-1
              if (j.eq.4) jm3 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
                ip2 =ip1+1
!              if (i.eq.n3d(1)-2) ip2 =2
                ip3 =ip2+1
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
                jp2 =jp1+1
              if (j.eq.n3d(2)-2) jp2 =2
                jp3 =jp2+1
              if (j.eq.n3d(2)-3) jp3 =2
!***********************************************************************
!      Find  Cell Peclet Number         
              Pe = abs(Re*Um(1,i,j,k)*(dx))      !taking only absolute quantity
          if(i.le.4.or.i.ge.(n3d(1)-3)) then
        divuu(1) = Um(1,i,j,k)*(Ut(jcnt,ip1,j,k)-Ut(jcnt,im1,j,k))/(2.0d0*dx)
       elseif(Pe.lt.2.0) then
!   divuu(1)= Um(1,i,j,k)*( (Ut(jcnt,im2,j,k)-8.0*Ut(jcnt,im1,j,k)+8.0*Ut(jcnt,ip1,j,k)-Ut(jcnt,ip2,j,k))/(12.0d0*dx) )
   divuu(1)=Um(1,i,j,k)*((-Ut(jcnt,im3,j,k)+9.0*Ut(jcnt,im2,j,k)-45.0*Ut(jcnt,im1,j,k)+45.0*Ut(jcnt,ip1,j,k)-9.0*Ut(jcnt,ip2,j,k)+Ut(jcnt,ip3,j,k))/(60.0d0*dx))
!   divuu(1)=Ut(1,i,j,k)*(Ut(jcnt,i-2,j,k)-8*Ut(jcnt,i-1,j,k)+8*Ut(jcnt,i+1,j,k)-Ut(jcnt,i+2,j,k))/(12.0d0*dx)
         else
term1 = (Abs(Um(1,i,j,k)))*(Ut(jcnt,ip2,j,k)-(4.0*(Ut(jcnt,ip1,j,k)))+(6.0*(Ut(jcnt,i,j,k)))-(4.0*(Ut(jcnt,im1,j,k)))+Ut(jcnt,im2,j,k))/(4.0d0*dx)    
term2 =     (Um(1,i,j,k))*(-Ut(jcnt,ip2,j,k)+(8.0*(Ut(jcnt,ip1,j,k)))-(8.0*(Ut(jcnt,im1,j,k)))+Ut(jcnt,im2,j,k))/(12.0d0*dx)    
  divuu(1) = term1+term2
        endif        
!***********************************************************************
             Pe= abs(Re*Um(2,i,j,k)*(dy))        !taking only absolute quantity
         if(j.le.4.or.j.ge.(n3d(2)-3)) then
        divuu(2) = Um(2,i,j,k)*(Ut(jcnt,i,jp1,k)-Ut(jcnt,i,jm1,k))/(2.0d0*dy)
        elseif(Pe.lt.2.0d0) then
!   divuu(2)= Um(2,i,j,k)*( (Ut(jcnt,i,jm2,k)-8.0*Ut(jcnt,i,jm1,k)+8.0*Ut(jcnt,i,jp1,k)-Ut(jcnt,i,jp2,k))/(12.0d0*dy) )
  divuu(2)=Um(2,i,j,k)*((-Ut(jcnt,i,jm3,k)+9.0*Ut(jcnt,i,jm2,k)-45.0*Ut(jcnt,i,jm1,k)+45.0*Ut(jcnt,i,jp1,k)-9.0*Ut(jcnt,i,jp2,k)+Ut(jcnt,i,jp3,k))/(60.0d0*dy))
! divuu(2)=Um(2,i,j,k)*((-Ut[i,j+2,k]+8*Ut[i,j+1,k]-8*Ut[i,j-1,k]+Ut[i,j-2,k])/(12.0d0*dy))
         else
term1 = (Abs(Um(2,i,j,k)))*(Ut(jcnt,i,jp2,k)-(4.0*(Ut(jcnt,i,jp1,k)))+(6.0*(Ut(jcnt,i,j,k)))-(4.0*(Ut(jcnt,i,jm1,k)))+Ut(jcnt,i,jm2,k))/(4.0d0*dy)    
term2 =     (Um(2,i,j,k))*(-Ut(jcnt,i,jp2,k)+(8.0*(Ut(jcnt,i,jp1,k)))-(8.0*(Ut(jcnt,i,jm1,k)))+Ut(jcnt,i,jm2,k))/(12.0d0*dy)    
  divuu(2) = term1+term2
        endif        
!***********************************************************************
           Pe= abs((Re*Um(3,i,j,k)*dXi)/Jac(k))   !taking only absolute quantity
       if(k.le.4.or.k.ge.(n3d(3)-3)) then
        divuu(3) = Um(3,i,j,k)*(Ut(jcnt,i,j,k+1)-Ut(jcnt,i,j,k-1))/(2.0d0*dXi)
       elseif(Pe.lt.2.0) then
!    divuu(3)= Um(3,i,j,k)*( (Ut(jcnt,i,j,k-2)-8.0*Ut(jcnt,i,j,k-1)+8.0*Ut(jcnt,i,j,k+1)-Ut(jcnt,i,j,k+2))/(12.0d0*dXi) )
  divuu(3)=Um(3,i,j,k)*((-Ut(jcnt,i,j,k-3)+9.0*Ut(jcnt,i,j,k-2)-45.0*Ut(jcnt,i,j,k-1)+45.0*Ut(jcnt,i,j,k+1)-9.0*Ut(jcnt,i,j,k+2)+Ut(jcnt,i,j,k+3))/(60.0d0*dXi))
!divuu(3)=Um(3,i,j,k)*((-Ut(jcnt,i,j,k+2)+8*Ut(jcnt,i,j,k+1)-8*Ut(jcnt,i,j,k-1)+Ut(jcnt,i,j,k-2))/(12*dz))
       else
term1 = (Abs(Um(3,i,j,k)))*(Ut(jcnt,i,j,k+2)-(4.0*(Ut(jcnt,i,j,k+1)))+(6.0*(Ut(jcnt,i,j,k)))-(4.0*(Ut(jcnt,i,j,k-1)))+Ut(jcnt,i,j,k-2))/(4.0*dXi)    
term2 =     (Um(3,i,j,k))*(-Ut(jcnt,i,j,k+2)+(8.0*(Ut(jcnt,i,j,k+1)))-(8.0*(Ut(jcnt,i,j,k-1)))+Ut(jcnt,i,j,k-2))/(12.0*dXi)    
  divuu(3) = (term1+term2)
       endif

if(jcnt.eq.1) then
!**************************************************************************** 
   Linear(1)=0.0d0
   Linear(2)=0.0d0
  if(k.le.4.or.k.ge.(n3d(3)-3)) then
  Linear(3)=Ut(3,i,j,k)*(Um(1,i,j,k+1)-Um(1,i,j,k-1))/(2.0d0*dXi)
  else
 Linear(3)=Ut(3,i,j,k)*((-Um(1,i,j,k-3)+9.0*Um(1,i,j,k-2)-45.0*Um(1,i,j,k-1)+45.0*Um(1,i,j,k+1)-9.0*Um(1,i,j,k+2)+Um(1,i,j,k+3))/(60.0d0*dXi))
! Linear(3)= Ut(3,i,j,k)*( (Um(1,i,j,k-2)-8*Um(1,i,j,k-1)+8*Um(1,i,j,k+1)-Um(1,i,j,k+2))/(12.0d0*dXi) )
  endif
!*****************************************************************************
endif

if(jcnt.eq.2) then 
!*****************************************************************************
 if(i.le.4.or.i.ge.(n3d(1)-3)) then
 Linear(1)=Ut(1,i,j,k)*(Um(2,ip1,j,k)-Um(2,im1,j,k))/(2.0d0*dx)
 else
 Linear(1)=Ut(1,i,j,k)*((-Um(2,im3,j,k)+9.0*Um(2,im2,j,k)-45.0*Um(2,im1,j,k)+45.0*Um(2,ip1,j,k)-9.0*Um(2,ip2,j,k)+Um(2,ip3,j,k))/(60.0d0*dx))
 endif

   Linear(2)=0.0d0

  if(k.le.4.or.k.ge.(n3d(3)-3)) then
  Linear(3)=Ut(3,i,j,k)*(Um(2,i,j,k+1)-Um(2,i,j,k-1))/(2.0d0*dXi)
  else
  Linear(3)=Ut(3,i,j,k)*((-Um(2,i,j,k-3)+9.0*Um(2,i,j,k-2)-45.0*Um(2,i,j,k-1)+45.0*Um(2,i,j,k+1)-9.0*Um(2,i,j,k+2)+Um(2,i,j,k+3))/(60.0d0*dXi))
  endif
!******************************************************************************
endif

if(jcnt.eq.3) then 
!******************************************************************************
  if(i.le.4.or.i.ge.(n3d(1)-3)) then
  Linear(1)=Ut(1,i,j,k)*(Um(3,ip1,j,k)-Um(3,im1,j,k))/(2.0d0*dx)
  else
  Linear(1)=Ut(1,i,j,k)*((-Um(3,im3,j,k)+9.0*Um(3,im2,j,k)-45.0*Um(3,im1,j,k)+45.0*Um(3,ip1,j,k)-9.0*Um(3,ip2,j,k)+Um(3,ip3,j,k))/(60.0d0*dx))
  endif

  Linear(2)=0.0d0

  if(k.le.4.or.k.ge.(n3d(3)-3)) then
  Linear(3)=Ut(3,i,j,k)*(Um(3,i,j,k+1)-Um(3,i,j,k-1))/(2.0d0*dXi)
  else
  Linear(3)=Ut(3,i,j,k)*((-Um(3,i,j,k-3)+9.0*Um(3,i,j,k-2)-45.0*Um(3,i,j,k-1)+45.0*Um(3,i,j,k+1)-9.0*Um(3,i,j,k+2)+Um(3,i,j,k+3))/(60.0d0*dXi))
  endif
!*******************************************************************************
endif
                                
 conv(jcnt,i,j,k) = divuu(1)+ Linear(1) +divuu(2)+Linear(2) + ((divuu(3)+ Linear(3))/Jac(k))

        end do
      end do  
    end do
 end do
    
     do k=2, n3d(3)-1
       do j=2,n3d(2)-1
        do i=2,n3d(1)-1
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
!***********************************************************************

      gradp(1,i,j,k)= (Press(ip1,j,k)-Press(im1,j,k))/(2.0d0*dx)

      gradp(2,i,j,k)= (Press(i,jp1,k)-Press(i,jm1,k))/(2.0d0*dy)

      gradp(3,i,j,k)= (Press(i,j,k+1)-Press(i,j,k-1))/(2.0d0*dXi*Jac(k))

      Q(4,i,j,k)=Ut(1,i,j,k)+(dt*(linforc(1,i,j,k)-conv(1,i,j,k)))
      Q(1,i,j,k)= Q(4,i,j,k)-(dt*gradp(1,i,j,k))
      Q(5,i,j,k)=Ut(2,i,j,k)+(dt*(linforc(2,i,j,k)-conv(2,i,j,k)))
      Q(2,i,j,k)= Q(5,i,j,k)-(dt*gradp(2,i,j,k))
      Q(6,i,j,k)=Ut(3,i,j,k)+(dt*(linforc(3,i,j,k)-conv(3,i,j,k)))
      Q(3,i,j,k)= Q(6,i,j,k)-(dt*gradp(3,i,j,k))

       end do
       end do  
       end do
!****************Source terms of RHS in the transport equations provisionally calculated****************************************
!****************Now itratively-solve for the provisional velocities and w/o pressure vlocities*********************************
    do jcnt =1,6 
             t=-dt/(Re)
!Gauss siedel iteration

      do it=1,maxiter

       tmpsum=0.0d0
      do k=2,n3d(3)-1
             do j=2,n3d(2)-1
                     do i=2,n3d(1)-1
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
!***********************************************************************
        res=Q(jcnt,i,j,k)-                       &
     &(1.0d0+ (t*Avp(i,j,k)))*Ut(jcnt,i,j,k)-    &
     &(t*Ave(i,j,k))*Ut(jcnt,ip1,j,k)-           &
     &(t*Avw(i,j,k))*Ut(jcnt,im1,j,k)-           &
     &(t*Avn(i,j,k))*Ut(jcnt,i,jp1,k)-           &
     &(t*Avs(i,j,k))*Ut(jcnt,i,jm1,k)-           &
     &(t*Avt(i,j,k))*Ut(jcnt,i,j,k+1)-           &
     &(t*Avb(i,j,k))*Ut(jcnt,i,j,k-1)

       tmpsum=tmpsum+abs(res)

                     end do
             end do
      end do
 
       if(tmpsum.eq.0)goto 191
       if(it.eq.1)sumnor=tmpsum 
                  sumav=tmpsum/sumnor 
 
                
      do k=2,n3d(3)-1
             do j=2,n3d(2)-1
                     do i=2,n3d(1)-1
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
!***********************************************************************
      Ut(jcnt,i,j,k)=(Q(jcnt,i,j,k)-                &
     &           (t*Ave(i,j,k))*Ut(jcnt,ip1,j,k)-   &
     &           (t*Avw(i,j,k))*Ut(jcnt,im1,j,k)-   &
     &           (t*Avn(i,j,k))*Ut(jcnt,i,jp1,k)-   &
     &           (t*Avs(i,j,k))*Ut(jcnt,i,jm1,k)-   &
     &           (t*Avt(i,j,k))*Ut(jcnt,i,j,k+1)-   &
     &           (t*Avb(i,j,k))*Ut(jcnt,i,j,k-1))/  &
     &      ( 1.0d0 +(t*Avp(i,j,k) ))
 
                     end do
             end do
      end do
      if(sumav.lt.errorstop)goto 191
      end do 
      print *,sumav,it,jcnt
191    CONTINUE
      END DO
! **********************now we have got up,vp,wp and provisional u,v,w********************                   
!****************************************************************************************
! now applying Momentum Interpolation to avoid spurious oscillations in pressure
!****************************************************************************************
!******** boundary conditions for the w/o pressure field velocity (up,vp,wp) ******************************************

          Ut(4,1,2:n3d(2)-1,1:n3d(3)) = 0.0d0      ! predicted up,vp,wp are zero at the entrance of the channel. 
          Ut(5,1,2:n3d(2)-1,1:n3d(3)) = 0.0d0      ! 
          Ut(6,1,2:n3d(2)-1,1:n3d(3)) = 0.0d0      !   

 Ut(4,n3d(1),2:n3d(2)-1,1:n3d(3)) =  (5.0*Ut(4,n3d(1)-1,2:n3d(2)-1,1:n3d(3))-4.0*Ut(4,n3d(1)-2,2:n3d(2)-1,1:n3d(3))+    &
                &         Ut(4,n3d(1)-3,2:n3d(2)-1,1:n3d(3)))/2.0d0 !second derivative of predicted up zero at exit!

 Ut(5,n3d(1),2:n3d(2)-1,1:n3d(3)) =  (5.0*Ut(5,n3d(1)-1,2:n3d(2)-1,1:n3d(3))-4.0*Ut(5,n3d(1)-2,2:n3d(2)-1,1:n3d(3))+    &
                &         Ut(5,n3d(1)-3,2:n3d(2)-1,1:n3d(3)))/2.0d0 !second derivative of predicted vp zero at exit!

 Ut(6,n3d(1),2:n3d(2)-1,1:n3d(3)) = (5.0*Ut(6,n3d(1)-1,2:n3d(2)-1,1:n3d(3))-4.0*Ut(6,n3d(1)-2,2:n3d(2)-1,1:n3d(3))+   &
                &         Ut(6,n3d(1)-3,2:n3d(2)-1,1:n3d(3)))/2.0d0 !second derivative of predicted wp zero at exit!      
!***********************************************************************************************************
                       Ut(4,1:n3d(1),2:n3d(2)-1,1) = 0.0    ! predicted up,vp,wp are zero at the bottom plate. 
                       Ut(5,1:n3d(1),2:n3d(2)-1,1) = 0.0    !    
                       Ut(6,1:n3d(1),2:n3d(2)-1,1) = 0.0    !

                      
                       Ut(4,1:n3d(1),2:n3d(2)-1,n3d(3)) = 0.0    ! predicted up,vp,wp are zero at the top plate. 
                       Ut(5,1:n3d(1),2:n3d(2)-1,n3d(3)) = 0.0    !  
                       Ut(6,1:n3d(1),2:n3d(2)-1,n3d(3)) = 0.0    !
!**********************************************************************************************************        
         ! ensure periodicity by mapping in spanwise direction for the up,vp,wp.

              Ut(4,1:n3d(1),1,1:n3d(3))= Ut(4,1:n3d(1),n3d(2)-1,1:n3d(3))                               
              Ut(4,1:n3d(1),n3d(2),1:n3d(3))= Ut(4,1:n3d(1),2,1:n3d(3)) 

              Ut(5,1:n3d(1),1,1:n3d(3))= Ut(5,1:n3d(1),n3d(2)-1,1:n3d(3))                               
              Ut(5,1:n3d(1),n3d(2),1:n3d(3))= Ut(5,1:n3d(1),2,1:n3d(3)) 

              Ut(6,1:n3d(1),1,1:n3d(3))= Ut(6,1:n3d(1),n3d(2)-1,1:n3d(3))                               
              Ut(6,1:n3d(1),n3d(2),1:n3d(3))= Ut(6,1:n3d(1),2,1:n3d(3)) 
!***********************************************************************************************************

        do k=2,n3d(3)-1
        do j=2,n3d(2)-1
        do i=2,n3d(1)-1
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
!***********************************************************************

               divuStar(1)= (((Ut(4,ip1,j,k) + Ut(4,i,j,k))/2.0d0) &
      &                                      -                     &
      &                (dt*((Press(ip1,j,k) - Press(i,j,k))/dx))&
      &                                      -                     & 
      &          (((Ut(4,i,j,k) + Ut(4,im1,j,k))/2.0d0)            &
      &                                      -                     &
      &          (dt*((Press(i,j,k) - Press(im1,j,k))/dx))))/dx
 
 
               divuStar(2)=(((Ut(5,i,jp1,k) + Ut(5,i,j,k))/2.0d0)   &
      &                                      -                      &
      &                 (dt*((Press(i,jp1,k) - Press(i,j,k))/dy))&
      &                                      -                      &
      &          (((Ut(5,i,j,k) + Ut(5,i,jm1,k))/2.0d0)             &
      &                                      -                      &
      &          (dt*((Press(i,j,k) - Press(i,jm1,k))/dy))))/dy
 
 
               divuStar(3)=2.0*(((Ut(6,i,j,k+1) + Ut(6,i,j,k))/2.0d0)   &
     &                                      -                       &
     &                 (dt*((Press(i,j,k+1) - Press(i,j,k))/(x3(k+1)-x3(k)))) &
     &                                      -                       &
     &          (((Ut(6,i,j,k) + Ut(6,i,j,k-1))/2.0d0)              &
     &                                      -                       &
     &           (dt*((Press(i,j,k) - Press(i,j,k-1))/(x3(k)-x3(k-1))))))/(x3(k+1)-x3(k-1))
         
         Qp(i,j,k)=(divuStar(1)+divuStar(2)+divuStar(3))/dt
           
          END DO
        END DO
      END DO
!*******initializing Pcorr before solution of pressure-poisson equation by GMRES*************************************
!on inflow and outflow faces in streamwise directions
                   Pcorr ( 1,2:n3d(2)-1,1:n3d(3) ) =0.0d0
                   Pcorr ( n3d(1),2:n3d(2)-1,1:n3d(3) ) =0.0d0
!on both the faces along spanwise directions   
!                   Pcorr ( 1:n3d(1),1,1:n3d(3) ) =0.0d0
!                   Pcorr ( 1:n3d(1),n3d(2),1:n3d(3) ) =0.0d0
!on both the walls of the channel 
                   Pcorr ( 1:n3d(1),2:n3d(2)-1,1) =0.0d0
                   Pcorr ( 1:n3d(1),2:n3d(2)-1,n3d(3)) =0.0d0
!*********************************************************************************************
!*******************Start GMRES for solution of Pressure-Poisson equation*********************
        sumav=1.0d0
        m = 0
!  Begin restart loop.
        do while ((sumav>errtol).and.(m<mmax))
        m = m+1              
        h = 0.0d0
        v= 0.0d0
        c= 0.0d0
        s= 0.0d0
        g = 0.0d0
        y = 0.0d0
! Matrix vector product for the initial residual.      

           do k = 2,n3d(3)-1
           do j = 2,n3d(2)-1
           do i = 2,n3d(1)-1
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
!***********************************************************************
    r(i,j,k)=Qp(i,j,k)-( Anw(i,j,k)*PCorr(im1,j  ,k)  + Ane(i,j,k)*PCorr(ip1,j  ,k)+     &
                        Ans(i,j,k)*PCorr(i  ,jm1,k)  + Ann(i,j,k)*PCorr(i  ,jp1,k)+     &
                        Anb(i,j,k)*PCorr(i  ,j  ,k-1)+ Ant(i,j,k)*PCorr(i  ,j  ,k+1)+   &
                        Anp(i,j,k)*PCorr(i  ,j  ,k))        
           enddo
           enddo
           enddo

!       This preconditioner changes with the number of processors!
        do k = 2,n3d(3)-1
          do j = 2,n3d(2)-1
           do i = 2,n3d(1)-1
            rhat(i,j,k) = -w*(r(i,j,k)+Anw(i,j,k)*rhat(i-1,j,k)+Ans(i,j,k)*rhat(i,j-1,k)+Anb(i,j,k)*rhat(i,j,k-1))*rAnp(i,j,k)
           end do
          end do
        end do

        rhat(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1) =((2-w)/w)*Anp(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)*rhat(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)

        do k = n3d(3)-1,2,-1
          do j = n3d(2)-1,2,-1
             do i = n3d(1)-1,2,-1
                rhat(i,j,k) =-w*(rhat(i,j,k)+Ane(i,j,k)*rhat(i+1,j,k)+Ann(i,j,k)*rhat(i,j+1,k)+Ant(i,j,k)*rhat(i,j,k+1))*rAnp(i,j,k)
             end do
          end do
        end do

        r(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1) = rhat(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)                                  

        rho=(sum(r(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)*r(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)))

!        call mpi_allreduce(loc_rho,rho,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

        rho = sqrt(rho)
        if(m.eq.1)sumnor=rho
                  sumav =rho/sumnor 

        g(1) =rho
        v(1,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)=r(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)/rho
        kit=0
! Begin gmres loop.
        do while((sumav > errtol).and.(kit < kmax))          
                   kit=kit+1                  
! Matrix vector product.      


           do k = 2,n3d(3)-1
           do j = 2,n3d(2)-1
           do i = 2,n3d(1)-1
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
!***********************************************************************

    v(kit+1,i,j,k)=  Anw(i,j,k)*v(kit,im1,j  ,k) +Ane(i,j,k)*v(kit,ip1,j  ,k)&
                    +Ans(i,j,k)*v(kit,i  ,jm1,k) +Ann(i,j,k)*v(kit,i  ,jp1,k)&
                    +Anb(i,j,k)*v(kit,i  ,j  ,k-1) +Ant(i,j,k)*v(kit,i  ,j  ,k+1)&
                    +Anp(i,j,k)*v(kit,i  ,j  ,k)        
           enddo
           enddo
           enddo

!       This is a preconditioner!
           do k = 2,n3d(3)-1
           do j = 2,n3d(2)-1
           do i = 2,n3d(1)-1
            rhat(i,j,k) =-w*(v(kit+1,i,j,k) + Anw(i,j,k)*rhat(i-1,j,k) + Ans(i,j,k)*rhat(i,j-1,k) + Anb(i,j,k)*rhat(i,j,k-1))*rAnp(i,j,k)
           end do
           end do
           end do

        rhat(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1) =  ((2-w)/w)*Anp(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)*rhat(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)

          do k= n3d(3)-1,2,-1
          do j = n3d(2)-1,2,-1
          do i = n3d(1)-1,2,-1
                rhat(i,j,k) =-w*(rhat(i,j,k)+Ane(i,j,k)*rhat(i+1,j,k) + Ann(i,j,k)*rhat(i,j+1,k)+ Ant(i,j,k)*rhat(i,j,k+1))*rAnp(i,j,k)
          end do
          end do
          end do
        v(kit+1,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1) = rhat(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)                                                 
! Begin modified GS. May need to reorthogonalize.   
                do k=1,kit                               

!       tempgmres=sum( v(k,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)*v(kit+1,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1) )

       h(k,kit)=sum( v(k,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)*v(kit+1,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1) )

!        call mpi_allreduce(tempgmres,h(k,kit),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
        v(kit+1,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)=v(kit+1,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)-h(k,kit)*v(k,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)

!      v(kit+1,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)=v(kit+1,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)-tempgmres*v(k,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)
        
                end do

!                      tempgmres=(sum(v(2:n3d(1)-1,2:n3d(2)-1,bn:en,kit+1)*v(2:n3d(1)-1,2:n3d(2)-1,bn:en,kit+1)))

!                    call mpi_allreduce(tempgmres,h(kit+1,kit),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

            h(kit+1,kit) =(sum(v(kit+1,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)*v(kit+1,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)))

                       h(kit+1,kit) = sqrt(h(kit+1,kit))

                       if (h(kit+1,kit).gt.0.0.or.h(kit+1,kit).lt.0.0) then

                           v(kit+1,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)=v(kit+1,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)/h(kit+1,kit)

                       end if

                       if (kit>1)  then                           
! Apply old Givens rotations to h(1:kit,kit).
                               do i=1,kit-1
                                 hik    =c(i)*h(i,kit)-s(i)*h(i+1,kit)
                                 hipk   =s(i)*h(i,kit)+c(i)*h(i+1,kit)
                                 h(i,kit) = hik
                                 h(i+1,kit) = hipk
                              end do
                    end if
                    nu=sqrt(h(kit,kit)**2 + h(kit+1,kit)**2)                
! May need better Givens implementation.
! Define and Apply new Givens rotations to h(kit:kit+1,kit).  
                    if (nu.gt.0.0) then
                        c(kit)=h(kit,kit)/nu
                        s(kit)=-h(kit+1,kit)/nu
                        h(kit,kit)=c(kit)*h(kit,kit)-s(kit)*h(kit+1,kit)
                        h(kit+1,kit)=0
                        gk    =c(kit)*g(kit) -s(kit)*g(kit+1)
                        gkp  =s(kit)*g(kit) +c(kit)*g(kit+1)
                        g(kit) = gk
                        g(kit+1) = gkp
                end if
                rho=abs(g(kit+1))
                  sumav =rho/sumnor 
                    mag(kit) = rho
! End of gmres loop.

     print *, sumav,kit

        end do
! h(1:kit,1:kit) is upper triangular matrix in QR.                                        
        y(kit) = g(kit)/h(kit,kit)
        do i = kit-1,1,-1
                y(i) = g(i)
                do k = i+1,kit
                        y(i) = y(i) -h(i,k)*y(k)
                end do
                y(i) = y(i)/h(i,i)
        end do
! Form linear combination.
         do i = 1,kit                                   
           PCorr(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1) = PCorr(2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1) + v(i,2:n3d(1)-1,2:n3d(2)-1,2:n3d(3)-1)*y(i)
         end do
! End restart loop.         
        end do

!*********************************************************************************************
! update of the correction pressure at the top and bottom walls.
                   Pcorr (1:n3d(1),2:n3d(2)-1,1)      =Pcorr (1:n3d(1),2:n3d(2)-1,2)
                   Pcorr (1:n3d(1),2:n3d(2)-1,n3d(3)) =Pcorr (1:n3d(1),2:n3d(2)-1,n3d(3)-1)

! ensure periodicity by mapping in spanwise direction for the correction pressure (Pcorr).
                   Pcorr (1:n3d(1),1,1:n3d(3)) =Pcorr(1:n3d(1),n3d(2)-1,1:n3d(3))
                   Pcorr (1:n3d(1),n3d(2),1:n3d(3)) =Pcorr(1:n3d(1),2,1:n3d(3)) 

!*********************************************************************************************
! Update  all Flow variables
!*********************************************************************************************
 Umax(:)=0.0
 KEmax=0.0
        do k=2,n3d(3)-1
        do j=2,n3d(2)-1
        do i=2,n3d(1)-1
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
!***********************************************************************
!************************************************************************
!* Find grad(P')        
!************************************************************************
!*                  P^n+1   =P^n      + P'
                Press(i,j,k)=Press(i,j,k)+Pcorr(i,j,k)

      gradpcorr(1)= (Pcorr(ip1,j,k)-Pcorr(im1,j,k))/(2.0d0*dx)

      gradpcorr(2)= (Pcorr(i,jp1,k)-Pcorr(i,jm1,k))/(2.0d0*dy)

      gradpcorr(3)= (Pcorr(i,j,k+1)-Pcorr(i,j,k-1))/(2.0d0*dXi*Jac(k))
 
!************************************************************************
!*                  V^n+1   =V^n      - dt*grad( P')
                 Ut(1,i,j,k)= Ut(1,i,j,k)-(dt*gradpcorr(1))
                 Ut(2,i,j,k)= Ut(2,i,j,k)-(dt*gradpcorr(2))
                 Ut(3,i,j,k)= Ut(3,i,j,k)-(dt*gradpcorr(3))

 KE1 = (Ut(1,i,j,k)* Ut(1,i,j,k)) +  (Ut(2,i,j,k)* Ut(2,i,j,k)) + (Ut(3,i,j,k)* Ut(3,i,j,k))

if (abs(Ut(1,i,j,k)).gt.Umax(1)) Umax(1)=abs(Ut(1,i,j,k)) 
if (abs(Ut(2,i,j,k)).gt.Umax(2)) Umax(2)=abs(Ut(2,i,j,k)) 
if (abs(Ut(3,i,j,k)).gt.Umax(3)) Umax(3)=abs(Ut(3,i,j,k)) 

     if (KE1.gt.KEmax) KEmax=KE1


          END DO
        END DO
      END DO
!applying the no-slip conditions
                       Ut(:,1:n3d(1),2:n3d(2)-1,1) =0.0
                       Ut(:,1:n3d(1),2:n3d(2)-1,n3d(3)) =0.0

! Update Pressure at boundaries using momentum equation

 !!! Bottom Boundary 
      
          do j=2,n3d(2)-1
          do i=2,n3d(1)-1

    fder =(-11.0*Ut(3,i,j,1) + 18.0*Ut(3,i,j,2) - 9.0*Ut(3,i,j,3) + 2.0*Ut(3,i,j,4))/(6.0*dXi) 

    divuu(3)=Um(3,i,j,1)*(fder/Jac(1)) 
       

divgrad(1)=(Ut(3,i-1,j,1)+Ut(3,i+1,j,1)-2*Ut(3,i,j,1))/(dx*dx)
divgrad(2)=(Ut(3,i,j-1,1)+Ut(3,i,j+1,1)-2*Ut(3,i,j,1))/(dy*dy)
divgrad(3)=(2*Ut(3,i,j,1) - 5*Ut(3,i,j,2) + 4*Ut(3,i,j,3) - Ut(3,i,j,4))/(dXi*dXi)

divgrad(3)  = (divgrad(3)/(Jac(1)*Jac(1))) - ((Jac2(1)*fder)/(Jac(1)*Jac(1)*Jac(1)))

! Finding The H and Pressure at the Top Wall          

Flux= -divuu(3) +((divgrad(1)+divgrad(2)+divgrad(3))/Re) +  linforc(3,i,j,1) 

Press(i,j,1)=((-6.0*dXi*Flux*Jac(1))+ 18.0*Press(i,j,2)-9.0*Press(i,j,3)+2.d0*Press(i,j,4))/11.0d0
  
            enddo
          enddo


 !!! Top Boundary 

           npln =n3d(3)   

          do j=2,n3d(2)-1
          do i=2,n3d(1)-1

fder =(11.0*Ut(3,i,j,npln) - 18.0*Ut(3,i,j,npln-1) + 9.0*Ut(3,i,j,npln-2) - 2.0*Ut(3,i,j,npln-3))/(6.0*dXi) 

divuu(3)=Um(3,i,j,npln)*(fder/Jac(npln)) 
       

divgrad(1)=(Ut(3,i-1,j,npln)+Ut(3,i+1,j,npln)-2*Ut(3,i,j,npln))/(dx*dx)
divgrad(2)=(Ut(3,i,j-1,npln)+Ut(3,i,j+1,npln)-2*Ut(3,i,j,npln))/(dy*dy)

divgrad(3) = (2*Ut(3,i,j,npln) - 5*Ut(3,i,j,npln-1) + 4*Ut(3,i,j,npln-2) - Ut(3,i,j,npln-3))/(dXi*dXi)

divgrad(3)  = (divgrad(3)/(Jac(npln)*Jac(npln))) - ((Jac2(npln)*fder)/(Jac(npln)*Jac(npln)*Jac(npln)))

! Finding The H and Pressure at the Top Wall          
Flux=-(divuu(3)) +((divgrad(1)+divgrad(2)+divgrad(3))/Re)  + linforc(3,i,j,npln)  

Press(i,j,npln)=((6.0*dXi*Flux*Jac(npln))+ 18.0*Press(i,j,npln-1)-9.0*Press(i,j,npln-2)+2*Press(i,j,npln-3))/11.0d0
  
            enddo
          enddo
! Outflow BCs at the streamwise 
Ut(:,n3d(1),2:n3d(2)-1,1:n3d(3)) =(5*Ut(:,n3d(1)-1,2:n3d(2)-1,1:n3d(3))-4.0*Ut(:,n3d(1)-2,2:n3d(2)-1,1:n3d(3))      &
               &           +Ut(:,n3d(1)-3,2:n3d(2)-1,1:n3d(3)))/2.0! second derivative of u' zero at front

Press(n3d(1),2:n3d(2)-1,1:n3d(3)) = (4.0*Press(n3d(1)-1,2:n3d(2)-1,1:n3d(3))-Press(n3d(1)-2,2:n3d(2)-1,1:n3d(3)))/3.0d0 ! first derivative  of p zero  at left

! Inflow BCs at the streamwise inflow
Ut(:,1,2:n3d(2)-1,1:n3d(3)) =0.0d0

do k=1,n3d(3)
do j=2,n3d(2)-1
        divuu(1) = Um(1,i,j,k)*((4*Ut(1,2,j,k)-Ut(1,3,j,k)-3*Ut(1,1,j,k)))/(2.0d0*dx)
       divgrad(1) =(Ut(1,1,j,k)+Ut(1,3,j,k)-2*Ut(1,2,j,k))/(dx*dx)
Press(1,j,k) =(2*dx*(divuu(1)-(divgrad(1)/Re)-linforc(1,1,j,k))-Press(3,j,k)+(4*Press(2,j,k))) /3.0d0
enddo
enddo
! Periodic BCs at the spanwise boundaries

                Ut(1,:,1,:)= Ut(1,:,n3d(2)-1,:)              ! ensure periodicity of u-prime by mapping in spanwise direction
                Ut(1,:,n3d(2),:)= Ut(1,:,2,:)                ! ensure periodicity of u-prime by mapping in spanwise direction

                Ut(2,:,1,:)= Ut(2,:,n3d(2)-1,:)              ! ensure periodicity  of v-prime by mapping in spanwise direction
                Ut(2,:,n3d(2),:)= Ut(2,:,2,:)                ! ensure periodicity of v-prime by mapping in spanwise direction

                Ut(3,:,1,:)= Ut(3,:,n3d(2)-1,:)              ! ensure periodicity of w-prime by mapping in spanwise direction
                Ut(3,:,n3d(2),:)= Ut(3,:,2,:)                ! ensure periodicity of w-prime by mapping in spanwise direction

                Press(:,1,:)= Press(:,n3d(2)-1,:)            ! ensure periodicity by mapping in spanwise direction
                Press(:,n3d(2),:)= Press(:,2,:)              ! ensure periodicity by mapping in spanwise direction

!*********************************************************************************************                      
      ! compute Divergence
  if(mod(nstep,1*disp)==0) then
!*********************************************************************************************
 tmp =0.0d0
!************************************************************************
        do k=2,n3d(3)-1
        do j=2,n3d(2)-1
        do i=2,n3d(1)-1
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
                im2 =im1-1
!              if (i.eq.3) im2 =n3d(1)-1
                im3 =im2-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
                jm2 =jm1-1
              if (j.eq.3) jm2 =n3d(2)-1
                jm3 =jm2-1
              if (j.eq.4) jm3 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
                ip2 =ip1+1
!              if (i.eq.n3d(1)-2) ip2 =2
                ip3 =ip2+1
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
                jp2 =jp1+1
              if (j.eq.n3d(2)-2) jp2 =2
                jp3 =jp2+1
              if (j.eq.n3d(2)-3) jp3 =2
!***********************************************************************

if(i.le.4.or.i.ge.(n3d(1)-3)) then
  divcon(1)=(Ut(1,ip1,j,k)-Ut(1,im1,j,k))/(2.0*dx)
else
  divcon(1)=(-Ut(1,im3,j,k)+9*Ut(1,im2,j,k)-45*Ut(1,im1,j,k)+45*Ut(1,ip1,j,k)-9*Ut(1,ip2,j,k)+Ut(1,ip3,j,k))/(60.0d0*dx)
endif

if(j.le.4.or.j.ge.(n3d(2)-3)) then
  divcon(2)=(Ut(2,i,jp1,k)-Ut(2,i,jm1,k))/(2.0*dy)
else
  divcon(2)=(-Ut(2,i,jm3,k)+9*Ut(2,i,jm2,k)-45*Ut(2,i,jm1,k)+45*Ut(2,i,jp1,k)-9*Ut(2,i,jp2,k)+Ut(2,i,jp3,k))/(60.0d0*dy)
endif
 
if(k.le.4.or.k.ge.(n3d(3)-3)) then
  divcon(3)=(Ut(3,i,j,k+1)-Ut(3,i,j,k-1))/(2.0*dXi*Jac(k))
else
  divcon(3)=(-Ut(3,i,j,k-3)+9*Ut(3,i,j,k-2)-45*Ut(3,i,j,k-1)+45*Ut(3,i,j,k+1)-9*Ut(3,i,j,k+2)+Ut(3,i,j,k+3))/(60.0d0*dXi*Jac(k))
endif
          tmp =tmp + (divcon(1)+divcon(2)+divcon(3))
           
          END DO
        END DO
      END DO

            print *,nstep,time,dsqrt(abs(tmp)/((n3d(2)-2)*(n3d(1)-2)*(n3d(3)-2)))
 end if 
!*********************************************************************************************
 Uat12(:,:)=0.0

!find k corresponding to z+ 
do kpln=1,3
        do k=2,n3d(3)
if (k.eq.zv(kpln)) then 
    

        do j=2,n3d(2)-1
        do i=2,n3d(1)-1


if (abs(Ut(1,i,j,k)).gt.Uat12(1,kpln)) Uat12(1,kpln)=abs(Ut(1,i,j,k)) 
if (abs(Ut(2,i,j,k)).gt.Uat12(2,kpln)) Uat12(2,kpln)=abs(Ut(2,i,j,k)) 
if (abs(Ut(3,i,j,k)).gt.Uat12(3,kpln)) Uat12(3,kpln)=abs(Ut(3,i,j,k)) 

           
          end do
        end do
!***************
 end if
!***************
end do
enddo

  open (unit=12,file='timeseriesatzp.xy',status='unknown',access='append')
  write(12,99)time*Re,Umax(1),Umax(2),Umax(3),Uat12(1,1),Uat12(2,1),Uat12(3,1),Uat12(1,2),Uat12(2,2),Uat12(3,2),Uat12(1,3),Uat12(2,3),Uat12(3,3),KEmax
  close(12)

99    format(15(1X,E20.10))
!*********************************************************************************************
   ke(:)=0.0

   do zpln=1,3
     do k=2,n3d(3)-1
    if(k.eq.zv(zpln)) then

      do i=1,n3d(1)
      do j=1,n3d(2)

     Utsq(zpln,i,j,k) =(Ut(1,i,j,k)* Ut(1,i,j,k))      ! calculation of u^2.
     Utnorm(zpln,i,j,k) =(Utsq(zpln,i,j,k)*dx*dy)      ! calculation of u^2*dx*dy.
     ke(zpln) = ke(zpln) + Utnorm(zpln,i,j,k)          ! calculation of the k.e. of the disturbance on a particular plane.

     end do
     end do
!****   
     end if
!****
     end do
     end do
!****

  open (unit=15,file='measureZplus.xy',status='unknown',access='append')

  
 
  write(15,77)time*Re,ke(1),ke(2),ke(3) 
 

  close(15)
77    format(4(1X,E20.10))
!*********************************************************************************************
! Writing Tecplot ASCII format Directly 

  !
  ! Save the snapshot of the flowfield 
  !
!***********************************************************************
  if(mod(nstep,snapshot)==0) then
!      saving snapshots
      write(filname,'(a,i5.5,a)')'snapat',nstep/snapshot,'.tp'
!      writing tecplot ascii format directly
      open (unit=12,file=filname,status='unknown')
!***********************************************************************

!      open (unit=12,file='3D.tp',status='unknown')
      write (12,01)'TITLE ="',time,nstep,'"'
      write (12,*)'Variables = "x","y","z","u","v","w","press","Pcorr","gradPz"'
      write (12,02)'ZONE k=',n3d(1),',j=',n3d(2),',i=',n3d(3),',DATAPACKING="POINT"'
      do i=1,n3d(1)
        write(12,*)
              do j=1,n3d(2)
                write(12,*)
                      do k=1,n3d(3)
    write(12,97)x1(i),x2(j),x3(k),Ut(1,i,j,k),Ut(2,i,j,k),Ut(3,i,j,k),Press(i,j,k),Pcorr(i,j,k),gradp(3,i,j,k) 
                      end do 
              end do
      end do
      close(12)
!***********************************************************************
01    format(A,F20.10,I9,A)
02    format(A,I3,A,I3,A,I3,A)
97    format(3(1X,F12.6),6(1X,E22.15))  
!*********************************************************************************************
  end if 
!************************************************
 if(mod(nstep,snapshot)==0) then 
 !
  ! Save the snapshot of the wmean-flowfield 
  !
!***********************************************************************
!      saving snapshots
      write(filname,'(a,i5.5,a)')'wmeanat',nstep/snapshot,'.tp'
      open (unit=12,file=filname ,status='unknown')
      write (12,03)'title ="',time,nstep,'"'
      write (12,*)'Variables = "x","z","wmean"'
      write (12,04)'ZONE k=',n3d(3),',i=',n3d(1),',DATAPACKING="POINT"'
      do k=1,n3d(3)
        write(12,*)
      
      do i=1,n3d(1)
    write(12,98)x1(i),x3(k),Um(3,i,34,k)
              end do
        
       end do
      close(12)
03    format(A,F20.10,I9,A)
04    format(A,I3,A,I3,A)
98    format(2(1X,F10.6),1(1X,E22.15))
!**********************************
end if
!*********************************************************************************************
 if(mod(nstep,bckp)==0) then
   open (unit=10,form='unformatted',file='backup.xy',status='unknown')
      write(10) time,nstep
      write(10) Ut,Press
      close(10)
 end if
!*********************************************************************************************
            
 end do ! while loop
!*********************************************************************************************                            
! Clean Up act
!*********************************************************************************************
     deallocate(Anw)
     deallocate(Ann)
     deallocate(Ane)
     deallocate(Ans)
     deallocate(Anb)
     deallocate(Ant)
     deallocate(Anp)
     deallocate(Ut)
     deallocate(gradp)
     deallocate(linforc)
     deallocate(Press)
     deallocate(Pcorr)
     deallocate(Q)
     deallocate(Qp)
     deallocate(x1)
     deallocate(x2)
     deallocate(x3)
     deallocate(Avw)
     deallocate(Avn)
     deallocate(Ave)
     deallocate(Avs)
     deallocate(Avb)
     deallocate(Avt)
     deallocate(Avp)
     deallocate(rAnp)
     deallocate(Conv)
     deallocate(Jac)
     deallocate( Xi)
     deallocate(Jac2)
     deallocate(Utnorm)
     deallocate(Utsq)
!*********************************************************************************************
end program channelinearised

