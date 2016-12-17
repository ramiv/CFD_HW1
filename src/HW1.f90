! CFD HW1 - elliptic problem by Rami Veiberman ID 308498831

program HW1
  use inMod
  use mathMod

  implicit none
  
  type(RunParms) :: run_p
  type(CaseParms) :: case_p


  character(len=300) :: path_in_file ! path to input file
  integer           :: inUnit = 10
  integer           :: ISTAT

  !integer ::  i,j

  CALL getarg(1,path_in_file) ! get input from user

  open(inUnit,file=path_in_file,FORM='FORMATTED',IOSTAT=ISTAT)

  if (ISTAT > 0) THEN
    WRITE(*,*), "Couldn't read" ,path_in_file
    CALL EXIT()
  end if
  rewind(unit=inUnit)

  CALL readRunParms(inUnit,run_p)
  CALL readCaseParms(inUnit,case_p)
  close(unit=inUnit)

  CALL step(case_p,run_p)
end
