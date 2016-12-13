MODULE outmod
  use inMod
  implicit none

  contains
    subroutine write_XY(outUnit,X,Y,rundata)
      integer,intent(inout)                 ::  outUnit  !number of the output unit
      real,dimension(:,:),intent(inout)     ::  X        !X result
      real,dimension(:,:),intent(inout)     ::  Y        !Y result
      type(RunParms),intent(IN)             ::  rundata  ! data about the run

      character(len=50)   ::  format_str !help string for creating formats

      integer :: i

      write(format_str,'(A,I2,A)'),"(",rundata%j_max,"(1X,E15.7))"
      
      write(outUnit,'(A,I3)'),"N = ",rundata%i_max
      write(outUnit,'(A,I3)'),"M = ",rundata%j_max

      write(outUnit,'(A)'),"X"
      do i=1,rundata%i_max
        write(outUnit,format_str),X(i,:)
      end do
      write(outUnit,'(A)'),"Y"
      do i=1,rundata%i_max
        write(outUnit,format_str),Y(i,:)
      end do
    end subroutine
end module
