! CFD HW1 - elliptic problem by Rami Veiberman ID 308498831
! This is the main file for running the program
! This program gets a main input file as an argument.
! 
! The main input file points to the run parameters of each case

program HW1
  use inMod
  use mathMod

  implicit none
  
  type(RunParms) :: run_p ! run parameters of each case
  type(CaseParms) :: case_p ! case parameters of each case


  character(len=300) :: path_in_file ! path to main input file
  integer           :: inUnit_main = 10  ! input unit of the main input file
  integer           :: inUnit_case= 15 ! input unit of the case input file
  integer           :: ISTAT

  logical           :: continue_run = .TRUE. ! while there are more cases to be
  
  character(len=300):: sec_input ! path to secondary input


  CALL getarg(1,path_in_file) ! get input from user

  ! open the main input file
  open(inUnit_main,file=path_in_file,FORM='FORMATTED',IOSTAT=ISTAT)

  ! if unsuccessful in opening exit run
  if (ISTAT > 0) THEN
    WRITE(*,*), "Couldn't read" ,path_in_file
    CALL EXIT()
  end if

  ! get the first filename
  continue_run = get_next_filename(inUnit_main,sec_input)
  do while (continue_run) ! while there are more files to be read
    ! open case file
    open(inUnit_case,file=sec_input,FORM='FORMATTED',IOSTAT=ISTAT)

    ! if unsuccessfull exit run
    if (ISTAT > 0) THEN
      WRITE(*,*), "Couldn't read" ,sec_input
      close(unit=inUnit_main)
      CALL EXIT()
    end if


    ! read run and case parameters from the case input file
    CALL readRunParms(inUnit_case,run_p)
    CALL readCaseParms(inUnit_case,case_p)

    ! solve for this case and export the ouptuts to the files
    CALL step(case_p,run_p)
    ! read next file
    continue_run = get_next_filename(inUnit_main,sec_input)
    close(unit=inUnit_case)
  end do
  close(unit=inUnit_main)
end
