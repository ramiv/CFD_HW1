! CFD HW1 - elliptic problem by Rami Veiberman ID 308498831

program HW1
  use inMod
  use mathMod

  implicit none
  
  type(RunParms) :: run_p
  type(CaseParms) :: case_p


  character(len=300) :: path_in_file ! path to input file
  integer           :: inUnit_main = 10
  integer           :: inUnit= 15
  integer           :: ISTAT

  logical           :: continue_run = .TRUE.
  character(len=300):: sec_input ! path to secondary input

  !integer ::  i,j

  CALL getarg(1,path_in_file) ! get input from user

  open(inUnit_main,file=path_in_file,FORM='FORMATTED',IOSTAT=ISTAT)

  if (ISTAT > 0) THEN
    WRITE(*,*), "Couldn't read" ,path_in_file
    CALL EXIT()
  end if

  continue_run = get_next_filename(inUnit_main,sec_input)
  do while (continue_run)
    open(inUnit,file=sec_input,FORM='FORMATTED',IOSTAT=ISTAT)

    if (ISTAT > 0) THEN
      WRITE(*,*), "Couldn't read" ,sec_input
      close(unit=inUnit_main)
      CALL EXIT()
    end if


    CALL readRunParms(inUnit,run_p)
    CALL readCaseParms(inUnit,case_p)

    CALL step(case_p,run_p)
    continue_run = get_next_filename(inUnit_main,sec_input)
  end do
  close(unit=inUnit_main)
end
