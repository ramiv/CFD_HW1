MODULE outmod
  use inMod
  implicit none

  contains
    subroutine write_XY(outUnit,X,Y,rundata)
      ! subroutine for wrting out XY  grid
      ! INPUT:
      !   - outUnit - externally open output unit outUnit
      !   - X,Y     - the calculated grid
      !   - rundata - type which holds the matrix size  
      integer,intent(inout)                 ::  outUnit  !number of the output unit
      real,dimension(:,:),intent(inout)     ::  X        !X result
      real,dimension(:,:),intent(inout)     ::  Y        !Y result
      type(RunParms),intent(IN)             ::  rundata  ! data about the run

      character(len=50)   ::  format_str !help string for creating formats

      integer :: i

      rewind(unit=outUnit) ! start writing in the beggining of the file
      ! create appropriate format for writing matrix (with space)
      write(format_str,'(A,I2,A)'),"(",rundata%j_max,"(1X,E15.7))"
      
      ! write the matrix size
      write(outUnit,'(A,I3)'),"N = ",rundata%i_max
      write(outUnit,'(A,I3)'),"M = ",rundata%j_max

      ! write header and then X matrix, Row wise starting from i=1
      write(outUnit,'(A)'),"X"
      do i=1,rundata%i_max
        write(outUnit,format_str),X(i,:)
      end do
      ! write header and then Y matrix, Row wise starting from i=1
      write(outUnit,'(A)'),"Y"
      do i=1,rundata%i_max
        write(outUnit,format_str),Y(i,:)
      end do
    end subroutine
end module
