MODULE mathmod
  use inMod
  implicit none

  real,parameter :: pi = acos(-1.0)
  contains

  subroutine Init_GRID(X,Y,run_p) 
    real,dimension(:,:),intent(inout) :: X
    real,dimension(:,:),intent(inout) :: Y
    type(RunParms),intent(IN)   :: run_p

    integer :: I,J

    real    :: dx
    real    :: R ! the distance of the stright section from the centerline
    real    :: Xmax ! right most point
    real    :: Louter ! the length of all the outer section
    real    :: ds     ! the distance between points on the outer section
    integer :: i0     ! the index of the last stright point
    real    :: theta0 ! angle to the last point the stright section
    real    :: theta_front ! the angle of the rounded arc
    real    :: Lfront ! the length of the rounded section
    integer :: Nfront ! number  of points in the rounded section
    real    :: dT     ! the angle partition in the rounded section
    real    :: theta  ! the angle used for the calculation

    X = 0.
    Y = 0.

    dx = 1./(run_p%i_LE - run_p%i_tel)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !                              j=1 Wall
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=run_p%i_TEL,run_p%i_LE
      X(i,1) = 1 - cos(pi*(run_p%I_LE - I)*dx/2.)
      Y(i,1) = -naca_profile(x(i,1),run_p%t_profile)
    end do

    do i= (run_p%i_LE + 1),run_p%I_TEU
      X(i,1) = 1 - cos(pi*(run_p%I_LE - I)*dx/2)
      Y(i,1) = naca_profile(x(i,1),run_p%t_profile)
    end do

    ! create the wake upper part
    do i=(run_p%I_TEU + 1),run_p%I_MAX
      X(i,1) = X(i-1,1) + ( X(i-1,1) - X(i-2,1))*run_p%XSF
      Y(i,1) = 0.
    end do
   
    ! create the wake lower part
    do i=1,(run_p%I_TEL - 1)
      X(i,1) = X(run_p%I_MAX - i + 1,1)
      Y(i,1) = 0.
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !                              i_max Wall
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Y(run_p%I_MAX,2) = run_p%dy
    do j=3,run_p%J_MAX
      Y(run_p%I_MAX,j) = Y(run_p%I_MAX,j-1) + & 
                        & (Y(run_p%I_MAX,j-1) - Y(run_p%I_MAX,j-2))*run_p%YSF
    end do

    do j=2,run_p%J_MAX
      X(run_p%I_MAX,j) = X(run_p%I_MAX,1) 
      Y(1,j) = -Y(run_p%I_MAX,j) 
      X(1,j) =  X(1,1)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !                              Outer wall
    !                              j = j_max
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! stright part is the same as the wake, only at differnt Y value
    do i=1,run_p%I_LE
      X(i,run_p%J_MAX) = X(i,1)
      Y(i,run_p%J_MAX) = Y(1,run_p%J_MAX)
    end do

    R    = ABS(Y(1,run_p%j_max))
    Xmax = ABS(X(1,1))

    Louter = pi*R + 2. * Xmax
    ds     = Louter / (run_p%I_MAX - 1)

    i0 = floor(Xmax/ds) + 1

    do i=2,i0
      Y(i,run_p%j_max) = -R
      X(i,run_p%j_max) = X(1,run_p%j_max) - ds*(i-1)
    end do

    theta0 = (pi/2. + atan(Y(i0,run_p%j_max)/X(i0,run_p%j_max)))
    theta_front = pi+2.*theta0

    Lfront = R*pi + 2.*abs(X(i0,run_p%j_max))
    Nfront = NINT(Lfront/ds)
    dT     = theta_front/real(Nfront)

    do i=(i0+1),(run_p%I_MAX - i0)
      theta = -pi/2 - dT*(i-i0) + theta0
      Y(i,run_p%j_max) = R*sin(theta)
      X(i,run_p%j_max) = R*cos(theta)
    end do

    do i=run_p%i_max , (run_p%i_max-i0+1), -1
      Y(i,run_p%j_max) = R
      X(i,run_p%j_max) = X(run_p%i_max,run_p%j_max) - ds*(run_p%i_max-i)
    end do

    ! Interpolate between the two walls!
    call interp2Dmat_eta(X)
    call interp2Dmat_eta(Y)
  end subroutine

  subroutine calc_BC(X,Y,x_xi,x_eta,y_xi,y_eta,case_p,PHI,PSI)
    real,dimension(:,:),intent(IN)    :: x_xi,x_eta
    real,dimension(:,:),intent(IN)    :: y_xi,y_eta
    type(CaseParms),intent(in)         :: case_p

    real,dimension(:,:),intent(inout)       :: X,Y
    real,dimension(:,:),intent(inout)       :: PHI,PSI

    integer i,j,N,M

    N = size(X,1)
    M = size(X,2)

    PHI = 0.
    PSI = 0.

    if (case_p%isPSI) THEN
      ! for xi = 1,N
      do i=1,N,(n-1)
        do j=2,M-1
          ! BC on xi = 1
          if ( abs(y_eta(i,j)) >= abs(x_eta(i,j)) ) THEN
            PSI(i,j) = - calc_2nd_diff(Y,i,j,2)/calc_diff(Y,i,j,2)
          else
            PSI(i,j) = - calc_2nd_diff(X,i,j,2)/calc_diff(X,i,j,2)
          end if
        end do
      end do
    call interp2Dmat_xi(PSI)
    end if

    if (case_p%isPHY) THEN
      do j=1,M,(M-1)
        do i=2,N-1
          ! BC on eta = 1
          if (abs(x_xi(i,j)) >= abs(y_xi(i,j) ) ) THEN
            PHI(i,j) = - calc_2nd_diff(X,i,j,1) / calc_diff(X,i,j,1)
          else
            PHI(i,j) = - calc_2nd_diff(Y,i,j,1) / calc_diff(Y,i,j,1)
          end if
        end do
      end do
      call interp2Dmat_eta(PHI)
    end if
  end subroutine

  subroutine calc_coef(x_xi,x_eta,y_xi,y_eta,alpha,beta,gama)
    ! returns alpha beta gamma
    real,dimension(:,:),intent(IN)    :: x_xi,x_eta
    real,dimension(:,:),intent(IN)    :: y_xi,y_eta
    real,dimension(:,:),intent(INOUT) :: alpha
    real,dimension(:,:),intent(INOUT) :: beta
    real,dimension(:,:),intent(INOUT) :: gama

    alpha = x_eta**2 + y_eta**2
    beta  = x_xi * x_eta + y_xi * y_eta
    gama  = x_xi**2 + y_xi**2
  end subroutine

  subroutine calc_metrics(X,Y,X_xi,X_eta,Y_xi,Y_eta)
    ! this function calcs the partial derivatives of X,Y 
    ! Inputs - X,Y pointers
    ! Outputs - X_xi, X_eta, Y_xi, Y_xi - externally allocated matrecies
    real,dimension(:,:),intent(inout) :: X
    real,dimension(:,:),intent(inout) :: Y

    real,dimension(:,:),intent(INOUT) :: X_xi,X_eta,Y_xi,Y_eta

    integer :: N,M
    integer :: I,J

    N = size(X,1)
    M = size(X,2)

    X_xi = 0.
    X_eta = 0.
    Y_xi = 0.
    Y_eta = 0.

    ! calculate the derivatives in the xi direction
    do j = 1,M
      do i = 2,N-1
        X_xi(i,j) = calc_diff(X,i,j,1)
        Y_xi(i,j) = calc_diff(Y,i,j,1)
      end do
    end do

    ! calculate the derivatives in the eta direction
    do i = 1,N
      do j = 2,M-1
        X_eta(i,j) = calc_diff(X,i,j,2)
        Y_eta(i,j) = calc_diff(Y,i,j,2)
      end do
    end do
  end subroutine

  real function calc_diff(X,i,j,dir)
    real,dimension(:,:),intent(inout) :: X
    integer,intent(IN) :: i,j
    integer,intent(IN) :: dir !may be 1 or 2

    select case (dir)
      case (1)
      calc_diff = ( X(i+1,j) - X(i-1,j) )/2.
      case (2)
      calc_diff = ( X(i,j+1) - X(i,j-1) )/2.
    end select
  end function

  real function calc_2nd_diff(X,i,j,dir)
    real,dimension(:,:),intent(inout) :: X
    integer,intent(IN) :: i,j
    integer,intent(IN) :: dir !may be 1 or 2

    select case (dir)
      case (1)
        calc_2nd_diff =  X(i+1,j) -2.*X(i,j) + X(i-1,j) 
      case (2)
        calc_2nd_diff =  X(i,j+1) -2.*X(i,j) + X(i,j-1) 
    end select
  end function

  real function naca_profile(x,t)
    real,intent(in) :: x
    real,intent(in) :: t

    real      :: x_int = 1.008930411365

    naca_profile = 5.*t * (0.2969*sqrt(x_int*x) - 0.126*(x_int*x) &
                  &-0.3516*((x_int*x)**2) + 0.2843*((x_int*x)**3) &
                  &-0.1015*((x_int*x)**4) )
  end function

  subroutine interp2Dmat_eta(X)
  ! Subroutine for interpulation in the eta/j direction
  ! INPUT/OUTPUT: X - the interpulated matrix
    real,dimension(:,:),intent(inout) :: X
    integer   :: N 
    integer   :: M 

    integer   :: i,j
    real      :: X0,Xmax
    
    N = SIZE(X,1)
    M = SIZE(X,2)

    do i = 2,N-1
      X0 = X(i,1)
      Xmax = X(i,M)
      do j=2,M-1
        X(i,j) = X0*(M-j)/M + Xmax*j/M
      end do
    end do
  end subroutine

  subroutine interp2Dmat_xi(X)
  ! Subroutine for interpulation in the xi/i direction
  ! INPUT/OUTPUT: X - the interpulated matrix
    real,dimension(:,:),intent(inout) :: X
    integer   :: N 
    integer   :: M 

    integer   :: i,j
    real      :: X0,Xmax

    N = SIZE(X,1)
    M = SIZE(X,2)

    do j = 2,M-1
      X0 = X(1,j)
      Xmax = X(M,j)
      do i=2,N-1
        X(i,j) = X0*(N-i)/N + Xmax*i/N
      end do
    end do
  end subroutine



end module
