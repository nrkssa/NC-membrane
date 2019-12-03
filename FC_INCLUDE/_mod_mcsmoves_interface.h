 Interface
  Function antigenflexE(ant) Result(flexE) BIND(C,name="antigenflexure_energy_")
  Use ISO_C_BINDING, ONLY:C_INT,C_DOUBLE
  Implicit None
  Integer(C_int)::ant
  Real(C_double)::flexE
  End Function antigenflexE
  
  Function bondedstate(antigenno) Result(bondstate) BIND(C,name="is_antigen_bonded_")
  Use ISO_C_BINDING, ONLY:C_INT,C_BOOL
  Implicit None
  Integer(C_int):: antigenno
  Logical(C_bool):: bondstate
  End function bondedstate

  Function does_bond_breaks(antigenno) Result(bond_break) BIND(C,name="does_bond_breaks_")
  USE ISO_C_BINDING, ONLY: C_INT,C_BOOL
  Implicit None
  Integer(C_int) :: antigenno
  Logical(C_bool):: bond_break
  End Function does_bond_breaks
  
  Function antigen_reaction_E_change(antigenno,tipold,tipnew) Result(dre) BIND(C,name="antigen_reaction_energy_change_")
  USE ISO_C_BINDING, ONLY:C_double,C_INT
  Implicit None
  Integer(C_int) :: antigenno
  Real(C_double) :: dre,tipold,tipnew
  End Function antigen_reaction_E_change

  Function antigen_reaction_energy(antigenno,tippos) Result(dre) BIND(C,name="antigen_reaction_energy_")
  USE ISO_C_BINDING, ONLY:C_double,C_INT
  Implicit None
  Integer(C_int) :: antigenno
  Real(C_double) :: dre,tippos
  End Function antigen_reaction_energy
 
  Subroutine Update_antigenbondedstate(i) BIND(C,name="update_antigenbondedstate_")
  USE ISO_C_BINDING, ONLY:C_INT
  Implicit None
  Integer(C_int) :: i
  End Subroutine Update_antigenbondedstate

  Subroutine Update_linkcells(vert,xold,yold,zold,ch) BIND(C,name="update_linkcells_")
  USE ISO_C_BINDING, ONLY:C_CHAR,C_INT,C_DOUBLE
  Implicit None
  Integer(C_int) :: vert
  Real(C_Double) :: xold,yold,zold
  Character(C_CHAR) :: ch
  End Subroutine Update_linkcells
  
  Subroutine print_cellnumber(vert,ch) BIND(C,name="print_cellnumber_")
  USE ISO_C_BINDING, ONLY:C_CHAR,C_INT
  Implicit None
  Integer(C_int) :: vert
  Character(C_CHAR) :: ch
  End Subroutine print_cellnumber

 End Interface
