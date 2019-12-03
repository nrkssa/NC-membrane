  Interface
  Function getx(i,ch) Result(xcoord) BIND(C,name="get_xcoordinate_")
  Use ISO_C_BINDING, ONLY:C_INT,C_CHAR,C_DOUBLE
  Character(C_CHAR) ::ch
  Integer (C_INT) :: i
  REAL(C_DOUBLE) :: xcoord
  End Function getx

  Function gety(i,ch) Result(ycoord) BIND(C,name="get_ycoordinate_")
  Use ISO_C_BINDING, ONLY:C_INT,C_CHAR,C_DOUBLE
  Character(C_CHAR) ::ch
  Integer (C_INT) :: i
  REAL(C_DOUBLE) :: ycoord
  End Function gety

  Function getz(i,ch) Result(zcoord) BIND(C,name="get_zcoordinate_")
  Use ISO_C_BINDING, ONLY:C_INT,C_CHAR,C_DOUBLE
  Character(C_CHAR) ::ch
  Integer (C_INT) :: i
  REAL(C_DOUBLE) :: zcoord
  End Function getz
  End Interface