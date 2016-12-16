! CFD HW1 - elliptic problem by Rami Veiberman ID 308498831

program HW1
  use inMod
  use outMod
  use mathMod

  implicit none
  !namelist /rundata/ t_profile, i_max, j_max
  

  type(RunParms) :: run_p
  type(CaseParms) :: case_p


  character(len=300) :: path_in_file ! path to input file
  character(len=300) :: path_out_file ! path to input file
  integer           :: inUnit = 10
  integer           :: outUnit= 12
  integer           :: ISTAT

  !integer ::  i,j

  real,dimension(:,:),allocatable :: X ! x output
  real,dimension(:,:),allocatable :: Y ! y output
  
  path_out_file = 'outTest.txt'

  CALL getarg(1,path_in_file) ! get input from user

  open(inUnit,file=path_in_file,FORM='FORMATTED',IOSTAT=ISTAT)
  open(outUnit,file=path_out_file,FORM='FORMATTED',STATUS='REPLACE')

  if (ISTAT > 0) THEN
    WRITE(*,*), "Couldn't read" ,path_in_file
  end if
  rewind(unit=inUnit)

  CALL readRunParms(inUnit,run_p)
  CALL readCaseParms(inUnit,case_p)

  allocate(X(run_p%i_max,run_p%j_max),Y(run_p%i_max,run_p%j_max))
  X = BAD_REAL
  Y = BAD_REAL


  call Init_GRID(X,Y,run_p)
  call write_XY(outUnit,X,Y,run_p)


  close(unit=inUnit)
  close(unit=outUnit)
  if (allocated(X)) deallocate(X)
  if (allocated(Y)) deallocate(Y)
end
