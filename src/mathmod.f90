MODULE mathmod
  use inMod
  use outMod
  implicit none


  ! Use external TRIDIAG subroutine
  interface
    subroutine TRIDIAG(A,B,C,D,U,N,IS,IE)
      REAL ::  A(0:N),B(0:N),C(0:N),D(0:N),U(0:N)
      INTEGER :: N,IS,IE
    end subroutine
  end interface

  real,parameter :: pi = acos(-1.0)
  contains

    subroutine step(case_p,run_p)
    ! master subroutine of the solution
    ! INPUT:
    !     run_p - run parameters structure
    !     case_p- case parameters structure
    ! OUTPUT (Exported to file)
    !     X - the grid x coordinates
    !     Y - the grid y coordinates
    !     Err_x Err_y - the error 
      type(CaseParms),intent(in)          :: case_p
      type(RunParms),intent(in)           :: run_p

      real,dimension(:,:),allocatable :: X,Y
      real,dimension(:,:),allocatable :: X_xi,X_eta,Y_xi,Y_eta ! the derivs
      real,dimension(:,:),allocatable :: PHI,PSI ! the control functions
      real,dimension(:,:),allocatable :: alpha,beta,gama ! factor matrecies
      real,dimension(:,:),allocatable :: Cx_n,Fx_n ! the solutions of the two eq in x
      real,dimension(:,:),allocatable :: Cy_n,Fy_n ! the solutions of the two eq in y
      real,dimension(:,:),allocatable :: Lx,Ly ! the RHS operator Lx and Ly

      real    :: error1_x = BAD_REAL ! the error of the first iteration in x
      real    :: error1_y = BAD_REAL ! the error of the first iteration in y
      real    :: error_x,error_y ! the current error in x and y
      real    :: partial_error_x,partial_error_y ! partial error 
      logical :: cont_run   ! flag which determins if the run converged
      integer :: N,M        ! help variables for the matrix size
      integer :: i_loop     ! current loop index
      integer :: outUnit_XY = 20 ! XY output file unit number
      integer :: outUnit_Err = 25! Err output file unit number
      integer :: outMod = 1 ! determins how often the error is prined.

      N = run_p%i_max
      M = run_p%j_max

      i_loop = 1 ! initiate the loop index 
      cont_run = .TRUE. ! initiate the run flag

      ! allocate the matrecies
      allocate(X(N,M),Y(N,M))
      allocate(Lx(N,M),Ly(N,M))
      allocate(X_xi(N,M),Y_xi(N,M),X_eta(N,M),Y_eta(N,M))
      allocate(PHI(N,M),PSI(N,M),Cx_n(N,M),Cy_n(N,M),Fx_n(N,M),Fy_n(N,M))
      allocate(alpha(N,M),beta(N,M),gama(N,M))

      ! open both output files
      open(unit=outUnit_XY,file=case_p%path_XY,&
           &STATUS='REPLACE',FORM='FORMATTED')
      open(unit=outUnit_Err,file=case_p%path_Err,&
           &STATUS='REPLACE',FORM='FORMATTED')
      
      ! first Init the grid
      call Init_GRID(X,Y,run_p)

      ! run while cont_run is TRUE and under 10,000 loops
      do while ( cont_run .AND. (i_loop < 10000))
        call calc_metrics(X,Y,X_xi,X_eta,Y_xi,Y_eta) ! calculate the X,Y Drivs
        call calc_control_func(X,Y,x_xi,x_eta,y_xi,y_eta,case_p,PHI,PSI) ! calc phi and psi control function
        call calc_coef(x_xi,x_eta,y_xi,y_eta,alpha,beta,gama) ! calc alpha,beta,gama

        ! solve the equation in the xi direction for X and Y
        call solve_Xi_eq(alpha,beta,gama,PHI,PSI,case_p%r,case_p%w,X,Fx_n,error_X)
        call solve_Xi_eq(alpha,beta,gama,PHI,PSI,case_p%r,case_p%w,Y,Fy_n,error_Y)

        ! write the error file header and save the first errors
        if (i_loop == 1) THEN
          error1_x = error_X
          error1_y = error_Y

          WRITE(outUnit_Err,'(A)'),"N  Err_x  Err_y"
        end if

        ! solve the equation in the eta direction for X and Y
        call solve_eta_eq(gama,case_p%r,Fx_n,Cx_n)
        call solve_eta_eq(gama,case_p%r,Fy_n,Cy_n)

        ! update X and Y and add the differences
        X = (X + Cx_n)
        Y = (Y + Cy_n)

        ! calculate the current errors and check if converged
        partial_error_x = log10(error_X)-log10(error1_X)
        partial_error_y = log10(error_Y)-log10(error1_Y)

        if ((partial_error_x < case_p%eps) .AND. (partial_error_y < case_p%eps)) THEN
          cont_run = .FALSE.
        end if

        ! write the current error to the file and screen (for live conv results)
        if (mod(i_loop,outMod) == 0) THEN
          write(outUnit_Err,'(1X,I6,2(1X,E15.7))'),i_loop,partial_error_x,partial_error_y
          write(*,'(A,I6,2(1X,E15.7))'),"N =",i_loop, partial_error_x,partial_error_y
        end if

        i_loop = i_loop + 1 ! update the loop index
      end do

      ! export the calculated X,Y matrecies
      call write_XY(outUnit_XY,X,Y,run_p) 
      ! close the file units
      close(unit=outUnit_XY) 
      close(unit=outUnit_Err)

      ! deallocate all the matrecies which we allocated before
      deallocate(X,Y)
      deallocate(Lx,Ly)
      deallocate(X_xi,Y_xi,X_eta,Y_eta)
      deallocate(PHI,PSI,Cx_n,Cy_n,Fx_n,Fy_n)
      deallocate(alpha,beta,gama)
    end subroutine

    subroutine solve_Xi_eq(alpha,beta,gama,phi,psi,r,w,X,Fx_n,error)
      ! sub for solving in xi direction either X or Y equation. Gets alpha,r,w,X
      ! and returns the Fx_n matrix and the error as defined
      ! Log10(Max(Abs(Lx)))
      real,dimension(:,:),intent(in)   :: alpha ! alpha factor matrix
      real,dimension(:,:),intent(in)   :: beta  ! beta  factor matrix
      real,dimension(:,:),intent(in)   :: gama  ! gama  factor matrix
      real,dimension(:,:),intent(in)   :: phi,psi ! phi and psi control func
      real,dimension(:,:),intent(in)   :: X     ! the current solution matrix,used for Lx calculation
      real,dimension(:,:),intent(inout)  :: Fx_n ! the soluton matrix
      real,intent(in)   :: r,w ! relaxation parameters 
      real,intent(out)    :: error ! the error Lx of the current iteration

      real,dimension(:),allocatable :: A_xi,B_xi,C_xi,D_xi ! A,B,C vectors of sol
      real,dimension(:),allocatable :: Fnj ! the solution for the j index
      real,dimension(:,:),allocatable :: Lx ! the RHS Lx
      integer :: J,N,M ! help indexes for index and matrix size

      ! get the matrix size
      N = size(X,1)
      M = size(X,2)

      ! init Fx_n
      Fx_n = BAD_REAL

      ! allocate the matrecies
      allocate(A_xi(N),B_xi(N),C_xi(N),D_xi(N))
      allocate(Fnj(N))
      allocate(Lx(N,M))

      ! init Fnj, Lx and D
      Fnj = BAD_REAL
      Lx = BAD_REAL
      D_xi = 0.

      call calc_Lx(X,alpha,beta,gama,phi,psi,Lx) ! calc Lx for RHS
      do j = 2,M-1
        call calc_A_xi(alpha,j,A_xi) ! calculate the A vector
        call calc_B_xi(alpha,r,j,B_xi) ! calculate the B vector
        call calc_C_xi(alpha,j,C_xi) ! calculate the C vector
        D_xi(2:N-1) = Lx(2:N-1,j)*r*w ! calculate the D vector using Lx, D in i=1,imax = 0

        ! get solution for this j index  using tri-diagonal 
        call TRIDIAG(A_xi,B_xi,C_xi,D_xi,Fnj,N,1,N)
        Fx_n(:,j) = Fnj ! update the output solution matrix for the current j index
      end do
      error = calc_error(Lx) ! calculate the error of all the Lx matrix

      ! deallocate the previously allocated matrecies
      deallocate(A_xi)
      deallocate(B_xi)
      deallocate(C_xi)
      deallocate(D_xi)
      deallocate(Fnj)
      deallocate(Lx)

    end subroutine

    subroutine solve_eta_eq(gama,r,Fx_n,Cx_n)
      ! sub for solving in eta direction either X or Y equation. Gets alpha,r,w,X
      ! and returns the Cx_n matrix 
      real,dimension(:,:),intent(in)   :: gama ! gama factor matrix
      real,intent(in)   :: r ! relaxation factor 
      real,dimension(:,:),intent(in)  :: Fx_n ! RHS of the equation
      real,dimension(:,:),intent(inout)  :: Cx_n ! solution matrix

      real,dimension(:),allocatable :: A_eta,B_eta,C_eta,D_eta,Cni ! the A,B,C,D
      ! Vectors and the solution of the ith index
      integer :: I,N,M ! index I and N,M which define matrix size

      ! Get the amtrix size
      N = size(gama,1)
      M = size(gama,2)

      Cx_n = 0. ! initiate the output matrix

      ! allocate the "working" vectors
      allocate(A_eta(M),B_eta(M),C_eta(M),D_eta(M),Cni(M))

      ! init D and the result vector Cni
      D_eta = 0.
      Cni   = 0.

      ! run for all the i indexes except the walls
      do i = 2,N-1
        call calc_A_eta(gama,i,A_eta)   ! calculate A vector
        call calc_B_eta(gama,r,i,B_eta) ! calculate B vector
        call calc_C_eta(gama,i,C_eta)   ! calculate C vector
        D_eta(2:M-1) = Fx_n(i,2:M-1)    ! calculate D vector, in i=1,i_max D=0

        ! solve for this i index using tri-diagonal solver
        call TRIDIAG(A_eta,B_eta,C_eta,D_eta,Cni,M,1,M)
        ! put the current row result into the output matrix Cx_n
        Cx_n(i,:) = Cni
      end do

      ! deallocate the matrecies 
      deallocate(A_eta)
      deallocate(B_eta)
      deallocate(C_eta)
      deallocate(D_eta)
      deallocate(Cni)
    end subroutine

    real function calc_error(Lx)
      ! this function calculates the error in Lx matrix
      ! Namely it finds the largest value |Lx| and returns it
      ! Input: Lx matrix
      ! Output: calc_error
      real,dimension(:,:),intent(IN) :: Lx
      integer             :: I,J,N,M ! indexes for loops and matrix size
      real                :: current_error ! the current error

      ! get matrix size
      N = size(Lx,1)
      M = size(Lx,2)

      calc_error = -1. ! init current error 
      ! run on all the values except for walls
      do i=2,N-1
        do j=2,M-1
          current_error = abs(Lx(i,j)) ! get current error
          if (calc_error < current_error) calc_error = current_error 
        end do
      end do
    end function

    subroutine calc_A_xi(alpha,j,A_xi)
      ! subroutine for calculating the A vector in xi eq
      ! INPUT:
      !       - alpha - factor matrix
      !       - j     - the col index 
      ! OUTPUT:
      !       - A_xi  -the output vector
      real,dimension(:,:),intent(in)  ::  alpha
      real,dimension(:),intent(inout) ::  A_xi
      integer,intent(in)              ::  j ! the eta location

      integer I,N ! I help index, N - vector length

      N = size(A_xi)
      A_xi = 0.
      
      ! calculate A(i) = -alpha(i,j)
      do i=2,N-1
        A_xi(i) = -alpha(i,j)
      end do
    end subroutine

    subroutine calc_B_xi(alpha,r,j,B_xi)
      ! subroutine for calculating the B vector  in xi eq
      ! INPUT:
      !       - alpha - factor matrix
      !       - r     - relaxation parameter
      !       - j     - the col index 
      ! OUTPUT:
      !       - B_xi  -the output vector
      real,dimension(:,:),intent(in)  ::  alpha
      real,dimension(:),intent(inout) ::  B_xi
      real,intent(in)                 ::  r ! relexation parameter
      integer,intent(in)              ::  j ! the eta location

      integer I,N ! I help index, N - vector length

      N = size(B_xi)
      B_xi = 1.
      
      ! calculate B(i) = r + 2alpha(i,j)
      do i=2,N-1
        B_xi(i) = (r + 2.*alpha(i,j) )
      end do
    end subroutine

    subroutine calc_C_xi(alpha,j,C_xi)
      ! subroutine for calculating the C vector  in xi eq
      ! INPUT:
      !       - alpha - factor matrix
      !       - j     - the col index 
      ! OUTPUT:
      !       - C_xi  -the output vector
      real,dimension(:,:),intent(in)  ::  alpha
      real,dimension(:),intent(inout) ::  C_xi
      integer,intent(in)              ::  j ! the eta location

      integer I,N ! I help index, N - vector length

      N = size(C_xi)
      C_xi = 0.
      
      ! calculate C(i) = -alpha(i,j)
      do i=2,N-1
        C_xi(i) = -alpha(i,j)
      end do
    end subroutine

    subroutine calc_A_eta(gama,i,A_eta)
      ! subroutine for calculating the A vector in eta eq
      ! INPUT:
      !       - alpha - factor matrix
      !       - j     - the col index 
      ! OUTPUT:
      !       - A_eta  -the output vector
      real,dimension(:,:),intent(in)  ::  gama
      real,dimension(:),intent(inout) ::  A_eta
      integer,intent(in)              ::  i ! the eta location

      integer J,M

      M = size(A_eta)
      A_eta = 0.
      
      ! A(j) = -gama(i,j)
      do j=2,M-1
        A_eta(j) = -gama(i,j)
      end do
    end subroutine

    subroutine calc_B_eta(gama,r,i,B_eta)
      ! subroutine for calculating the C vector  in eta eq
      ! INPUT:
      !       - gama - factor matrix
      !       - j     - the col index 
      ! OUTPUT:
      !       - B_eta  -the output vector
      real,dimension(:,:),intent(in)  ::  gama
      real,dimension(:),intent(inout) ::  B_eta
      real,intent(in)                 ::  r ! relexation parameter
      integer,intent(in)              ::  i ! the eta location

      integer J,M

      M = size(B_eta)
      B_eta = 1.
      
      ! B(j) = (R + 2gama(i,j) )
      do j=2,M-1
        B_eta(j) = (r + 2.*gama(i,j) )
      end do
    end subroutine

    subroutine calc_C_eta(gama,i,C_eta)
      ! subroutine for calculating the C vector  in eta eq
      ! INPUT:
      !       - alpha - factor matrix
      !       - j     - the col index 
      ! OUTPUT:
      !       - C_eta  -the output vector
      real,dimension(:,:),intent(in)  ::  gama
      real,dimension(:),intent(inout) ::  C_eta
      integer,intent(in)              ::  i ! the eta location

      integer J,M

      M = size(C_eta)
      C_eta = 0.
      
      ! C(j) = -gama(i,j)
      do j=2,M-1
        C_eta(j) = -gama(i,j)
      end do
    end subroutine

    subroutine calc_Lx(X,alpha,beta,gama,phi,psi,Lx)
      ! This function calculates the Lx
      ! 
      ! Input: 
      !       - X - the matrix at time point n
      !       - alpha,beta,gama - the coefficents of the equation
      !       - phi,psi         - the control functions
      ! Output: Lx - the calculated output matrix, needs to be allocated
      ! outside of the subroutine
      ! 
      real,dimension(:,:),intent(in)  :: X
      real,dimension(:,:),intent(in)  :: alpha,beta,gama
      real,dimension(:,:),intent(in)  :: phi,psi

      real,dimension(:,:),intent(inout) :: Lx


      integer     :: I,J ! help indecies for row and col
      integer     :: N,M ! size of the matrix (N,M)
      ! I assume that Lx was externally allocated

      ! N,M is the size of the matrix
      N = size(X,1)
      M = size(X,2)

      ! init the Lx matrix
      Lx = 0.

      ! for all the indexes except the walls
      do i = 2,N-1
        do j = 2,M-1
          Lx(i,j) = alpha(i,j) * (X(i+1,j) -2.0*X(i,j) + X(i-1,j) + phi(i,j) *&
                    &(X(i+1,j)-X(i-1,j))/2. ) &
                    &- beta(i,j)/2. * (X(i+1,j+1) - X(i+1,j-1) &
                    &- X(i-1,j+1) + X(i-1,j-1) ) &
                    &+ gama(i,j) * (X(i,j+1) - 2.*X(i,j) + X(i,j-1) + psi(i,j) *&
                    &(X(i,j+1) - X(i,j-1))/2. )
        end do
      end do

    end subroutine

  subroutine Init_GRID(X,Y,run_p) 
    ! Subroutine for calculating the initial Grid 
    ! INPUT:
    !       - run_p - a structure which holds all the grid parametrs like
    !                 i_max,j_max etc (see inmod.f90)
    ! OUTPUT:
    !       - X,Y   - the matrecies 
    real,dimension(:,:),intent(inout) :: X
    real,dimension(:,:),intent(inout) :: Y
    type(RunParms),intent(IN)   :: run_p

    integer :: I,J ! help indecies for row and col

    real    :: dx ! the differences between two points in xi direction 
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

    ! init X,Y data
    X = 0.
    Y = 0.

    dx = 1./(run_p%i_LE - run_p%i_tel)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !                              j=1 Wall
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! calc points on profile between the lower TE to the LE
    do i=run_p%i_TEL,run_p%i_LE
      X(i,1) = 1 - cos(pi*(run_p%I_LE - I)*dx/2.)
      Y(i,1) = -naca_profile(x(i,1),run_p%t_profile)
    end do

    ! calc points on profile between the LE and the upper TE
    do i= (run_p%i_LE + 1),run_p%I_TEU
      X(i,1) = 1 - cos(pi*(run_p%I_LE - I)*dx/2)
      Y(i,1) = naca_profile(x(i,1),run_p%t_profile)
    end do

    ! calc wake points in the upper part
    do i=(run_p%I_TEU + 1),run_p%I_MAX
      X(i,1) = X(i-1,1) + ( X(i-1,1) - X(i-2,1))*run_p%XSF
      Y(i,1) = 0.
    end do
   
    ! calc wake points in the lower part
    do i=1,(run_p%I_TEL - 1)
      X(i,1) = X(run_p%I_MAX - i + 1,1)
      Y(i,1) = 0.
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !                              i_max Wall
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Y(run_p%I_MAX,2) = run_p%dy ! get dy from input
    ! calcualte the upper back wall
    do j=3,run_p%J_MAX
      Y(run_p%I_MAX,j) = Y(run_p%I_MAX,j-1) + & 
                        & (Y(run_p%I_MAX,j-1) - Y(run_p%I_MAX,j-2))*run_p%YSF
    end do

    ! calculate the X points and lower part
    do j=2,run_p%J_MAX
      X(run_p%I_MAX,j) = X(run_p%I_MAX,1) 
      Y(1,j) = -Y(run_p%I_MAX,j) 
      X(1,j) =  X(1,1)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !                              Outer wall
    !                              j = j_max
    !
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

    ! lower part
    do i=2,i0
      Y(i,run_p%j_max) = -R
      X(i,run_p%j_max) = X(1,run_p%j_max) - ds*(i-1)
    end do

    theta0 = (pi/2. + atan(Y(i0,run_p%j_max)/X(i0,run_p%j_max)))
    theta_front = pi+2.*theta0

    Lfront = R*pi + 2.*abs(X(i0,run_p%j_max))
    Nfront = NINT(Lfront/ds)
    dT     = theta_front/real(Nfront)

    ! calculate round part 
    do i=(i0+1),(run_p%I_MAX - i0)
      theta = -pi/2 - dT*(i-i0) + theta0
      Y(i,run_p%j_max) = R*sin(theta)
      X(i,run_p%j_max) = R*cos(theta)
    end do

    ! calculate the upper part
    do i=run_p%i_max , (run_p%i_max-i0+1), -1
      Y(i,run_p%j_max) = R
      X(i,run_p%j_max) = X(run_p%i_max,run_p%j_max) - ds*(run_p%i_max-i)
    end do

    ! Interpolate between the two walls!
    call interp2Dmat_eta(X)
    call interp2Dmat_eta(Y)
  end subroutine

  subroutine calc_control_func(X,Y,x_xi,x_eta,y_xi,y_eta,case_p,PHI,PSI)
    ! calculates the Control funcctions PHI and PSI
    ! INPUT: 
    !       x_xi,x_eta,y_xi,y_eta - the metric derivetives
    !       case_p - case parameters structure
    !       X,Y    - the grid at the current time point n
    ! OUTPUT:
    !       PHI,PSI - the output control functions
    real,dimension(:,:),intent(IN)    :: x_xi,x_eta
    real,dimension(:,:),intent(IN)    :: y_xi,y_eta
    type(CaseParms),intent(in)         :: case_p

    real,dimension(:,:),intent(inout)       :: X,Y
    real,dimension(:,:),intent(inout)       :: PHI,PSI

    integer i,j,N,M

    ! get matrix size
    N = size(X,1)
    M = size(X,2)

    ! init the control functions
    PHI = 0.
    PSI = 0.

    if (case_p%isPSI) THEN
      ! for xi = 1,N
      do i=1,N,(N-1)
        do j=2,M-1
          ! BC on xi = 1
          if ( abs(y_eta(i,j)) >= abs(x_eta(i,j)) ) THEN
            PSI(i,j) = - (calc_2nd_diff(Y,i,j,2))/(calc_diff(Y,i,j,2))
          else
            PSI(i,j) = - (calc_2nd_diff(X,i,j,2))/(calc_diff(X,i,j,2))
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
            PHI(i,j) = - (calc_2nd_diff(X,i,j,1)) / (calc_diff(X,i,j,1))
          else
            PHI(i,j) = - (calc_2nd_diff(Y,i,j,1)) / (calc_diff(Y,i,j,1))
          end if
        end do
      end do
      call interp2Dmat_eta(PHI)
    end if
  end subroutine

  subroutine calc_coef(x_xi,x_eta,y_xi,y_eta,alpha,beta,gama)
    ! subroutine for calculating the alpha,beta,gama factor metrecies
    !
    ! INPUT:
    !     - X_xi,X_eta,Y_xi,Y_eta - the derivatives of X,Y
    ! OUTPUT:
    !     - alpha,beta,gama - the factor matrecies
    real,dimension(:,:),intent(IN)    :: x_xi,x_eta
    real,dimension(:,:),intent(IN)    :: y_xi,y_eta
    real,dimension(:,:),intent(INOUT) :: alpha
    real,dimension(:,:),intent(INOUT) :: beta
    real,dimension(:,:),intent(INOUT) :: gama

    integer :: I,J,N,M ! help indecies and matrix size (N,M)

    N = size(x_xi,1)
    M = size(x_xi,2)

    do i=1,N
      do j=1,M
        alpha(i,j) = x_eta(i,j)*x_eta(i,j) + y_eta(i,j)*y_eta(i,j)
        beta(i,j)  = x_xi(i,j)* x_eta(i,j) + y_xi(i,j)* y_eta(i,j)
        gama(i,j)  = x_xi(i,j)* x_xi(i,j)  + y_xi(i,j)* y_xi(i,j)
      end do
    end do
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

    ! init the metrecies
    X_xi = 0.
    X_eta = 0.
    Y_xi = 0.
    Y_eta = 0.

    ! calculate the derivatives in the xi direction all j indecies and i=2..N-1
    do j = 1,M
      do i = 2,N-1
        X_xi(i,j) = calc_diff(X,i,j,1)
        Y_xi(i,j) = calc_diff(Y,i,j,1)
      end do
    end do

    ! calculate the derivatives in the eta direction all i indecies and j=2..M-1
    do i = 1,N
      do j = 2,M-1
        X_eta(i,j) = calc_diff(X,i,j,2)
        Y_eta(i,j) = calc_diff(Y,i,j,2)
      end do
    end do
  end subroutine

  real function calc_diff(X,i,j,dir)
    ! calcualtes the first derivative in the 'dir' direction of 2d matrix X at
    ! location i,j please note that the function DOES NOT devide by dx or dy for
    ! simplicity. if needed, please devide outside the function
    real,dimension(:,:),intent(in) :: X
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
    ! calculates the second derivative in the 'dir' direction of 2d matrix X at
    ! location i,j. Please note that the function DOES NOT devide by dx**2 or
    ! dy**2 for simplicity. if needed, please devide outside the function
    real,dimension(:,:),intent(in) :: X
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
    ! function for calculating the naca_profile point at x for profiel with 
    ! width t
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
      Xmax = X(N,j)
      do i=2,N-1
        X(i,j) = X0*(N-i)/N + Xmax*i/N
      end do
    end do
  end subroutine

end module
