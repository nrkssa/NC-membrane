 Interface
  Subroutine write_vtkformat(samp_no,ch) BIND(C,name="write_vtkformat_")
  Use ISO_C_BINDING, ONLY:C_INT,C_DOUBLE,C_CHAR
  Implicit None
  Integer(C_INT) :: samp_no
  Character(C_CHAR) :: ch
  End Subroutine write_vtkformat
 End Interface			
