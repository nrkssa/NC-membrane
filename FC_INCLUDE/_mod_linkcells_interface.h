   Interface
   Function get_headlist(ch) Result(head) BIND(C,name="get_headlist_")
   Use ISO_C_BINDING, ONLY:C_INT,C_CHAR,C_PTR
   TYPE(C_PTR) :: head
   Character(C_CHAR) ::ch
   End Function get_headlist

   Function get_lscllist(ch) Result(lscl) BIND(C,name="get_lscllist_")
   Use ISO_C_BINDING, ONLY:C_INT,C_CHAR,C_PTR
   Type(C_PTR) :: lscl
   Character(C_CHAR) ::ch
   End Function get_lscllist

   Function getcells(rsize,i,ch)Result(ring) BIND(C,name="getcells_")
   Use ISO_C_BINDING, ONLY:C_INT,C_CHAR,C_PTR,C_INT
   Type(C_PTR) :: ring
   Character(C_CHAR) ::ch
   Integer(C_INT) :: i,rsize
   End Function getcells

  Subroutine print_celldata(vert,ch) BIND(C,name="print_celldata_")
  USE ISO_C_BINDING, ONLY:C_CHAR,C_INT
  Implicit None
  Integer(C_int) :: vert
  Character(C_CHAR) :: ch
  End Subroutine print_celldata

  End Interface
