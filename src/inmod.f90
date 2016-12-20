MODULE inmod
  implicit none

  real,parameter  :: BAD_REAL   = -999.0
  integer,parameter  :: BAD_INT = -999

  type RunParms
    real      :: t_profile     !profile width
    integer   :: i_max          ! max index in xi direction
    integer   :: j_max          ! max index in eta direction
    integer   :: i_TEL          ! index of the lower trailing edge
    integer   :: i_LE           ! index of the leading edge
    integer   :: i_TEU          ! index of the Upper trailing edge
    real      :: dy            ! the y difference near the wall
    real      :: XSF           ! streach in X direction parameter
    real      :: YSF           ! streach in Y direction parameter

  end type 

  type CaseParms
    logical   :: isPHY        ! determins if phy control function is used
    logical   :: isPSI        ! determins if PSI control function is used

    real      :: w            ! relaxation parameter
    real      :: r            ! relaxation parameter

    real   :: eps          ! portion between the last and first errors 
    character(len=300) :: path_XY
    character(len=300) :: path_Err
  end type


contains

  !  Subroutine for reading the basic run parameters which do not change between
  !  cases. 
  !
  !  Input: iUnit - unit number for reading the file
  !  Output: parms structure which holds the various parameters as described
  !  below
  SUBROUTINE readRunParms(iUnit,parms)
    
    integer,intent(INOUT):: iUnit ! the input unit
    type(RunParms),intent(out) :: parms

    integer   :: ISTAT           ! status for if the read was wrong
    logical   :: badInp = .FALSE.! status for bad input

    real      :: t_profile  = BAD_REAL     !profile width
    integer   :: i_max      = BAD_INT      ! max index in xi direction
    integer   :: j_max      = BAD_INT      ! max index in eta direction
    integer   :: i_TEL      = BAD_INT      ! index of the lower trailing edge
    integer   :: i_LE       = BAD_INT      ! index of the leading edge
    integer   :: i_TEU      = BAD_INT      ! index of the Upper trailing edge
    real      :: dy         = BAD_REAL     ! the y difference near the wall
    real      :: XSF        = BAD_REAL     ! streach in X direction parameter
    real      :: YSF        = BAD_REAL     ! streach in Y direction parameter

    namelist /rundata/ t_profile,i_max,j_max,i_TEL,i_LE,i_TEU,dy,XSF,YSF


    read(unit=iUnit,nml=rundata,IOSTAT=ISTAT) ! read the namelist

    IF (ISTAT > 0) THEN
      WRITE(*,*),"Could not read %rundata namelist, please check input!"
      close(unit=iUnit)
      call EXIT()
    end if

    if (t_profile == BAD_REAL ) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %rundata namelist 't_profile' variable, please check input!"
    end if

    if (i_max == BAD_INT ) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %rundata namelist 'i_max' variable, please check input!"
    end if

    if (j_max == BAD_INT) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %rundata namelist 'j_max' variable, please check input!"
    end if

    if (i_TEL == BAD_INT ) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %rundata namelist 'i_TEL' variable, please check input!"
    end if

    if (i_LE == BAD_INT ) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %rundata namelist 'i_LE' variable, please check input!"
    end if

    if (i_TEU == BAD_INT ) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %rundata namelist 'i_TEU' variable, please check input!"
    end if

    if (dy == BAD_REAL ) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %rundata namelist 'dy' variable, please check input!"
    end if

    if (XSF == BAD_REAL ) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %rundata namelist 'XSF' variable, please check input!"
    end if

    if (YSF == BAD_REAL ) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %rundata namelist 'YSF' variable, please check input!"
    end if



    if (badInp) THEN
      close(unit=iUnit)
      CALL EXIT()
    end if

    parms%t_profile = t_profile
    parms%i_max     = i_max    
    parms%j_max     = j_max    
    parms%i_TEL     = i_TEL    
    parms%i_LE      = i_LE     
    parms%i_TEU     = i_TEU    
    parms%dy        = dy       
    parms%XSF       = XSF      
    parms%YSF       = YSF      

  END SUBROUTINE readRunParms
  !  Subroutine for reading the basic case parameters   !

  !  Input: iUnit - unit number for reading the file
  !  Output: parms structure which holds the various parameters as described
  !  below
  SUBROUTINE readCaseParms(iUnit,parms)
    
    integer,intent(INOUT):: iUnit ! the input unit
    type(CaseParms),intent(out) :: parms

    integer   :: ISTAT           ! status for if the read was wrong
    logical   :: badInp = .FALSE.! status for bad input

    logical   :: isPHY = .FALSE.       ! determins if phy control function is used
    logical   :: isPSI = .FALSE.       ! determins if PSI control function is used

    real      :: w    = BAD_REAL        ! relaxation parameter
    real      :: r    = BAD_REAL        ! relaxation parameter
    real      :: eps  = BAD_REAL        ! relaxation parameter
    character(len=300) :: path_XY = ''     ! the path to export the XY data
    character(len=300) :: path_Err = ''    ! the path to export the error

    namelist /casedata/isPHY,isPSI,w,r,eps,path_XY,path_Err

    read(unit=iUnit,nml=casedata,IOSTAT=ISTAT) ! read the namelist

    IF (ISTAT > 0) THEN
      WRITE(*,*),"Could not read %casedata namelist, please check input!"
      close(unit=iUnit)
      call EXIT()
    end if


    if (w == BAD_REAL ) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %casedata namelist 'w' variable, please check input!"
      WRITE(*,*),"w read = ",w
    end if

    if (r == BAD_REAL ) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %casedata namelist 'r' variable, please check input!"
      WRITE(*,*),"r read = ",r
    end if

    if (eps == BAD_REAL ) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %casedata namelist 'eps' variable, please check input!"
      WRITE(*,*),"eps read = ",eps
    end if


    if (path_XY == '' ) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %rundata namelist 'path_XY' variable, please check input!"
    end if

    if (path_Err == '' ) THEN
      badInp = .TRUE.
      WRITE(*,*),"Could not read %rundata namelist 'path_Err' variable, please check input!"
    end if


    if (badInp) then
      close(unit=iUnit)
      CALL EXIT()
    end if

  parms%isPHY       = isPHY
  parms%isPSI       = isPSI
  parms%r           = r
  parms%w           = w
  parms%eps         = eps
  parms%path_XY     = path_XY      
  parms%path_Err    = path_Err      

  END SUBROUTINE readCaseParms

  logical function get_next_filename(funit_main,nextFile)
    integer,intent(inout) :: funit_main
    character(len=300),intent(out) :: nextFile
    character(len=1)  :: firstChar
      
    integer   :: ISTAT = 0

    get_next_filename = .False.
    

    do while ( (ISTAT == 0) .AND. (.NOT. get_next_filename) ) 
      read(funit_main,'(A)',IOSTAT=ISTAT),nextfile
      if (ISTAT == 0) THEN
        nextfile = TRIM(ADJUSTL(nextfile))
        firstChar = nextfile(1:1)
        if (firstChar /= "!" .AND. firstChar /= "") THEN
          get_next_filename = .TRUE.
        end if
      end if
    end do


  end function

END MODULE inmod
