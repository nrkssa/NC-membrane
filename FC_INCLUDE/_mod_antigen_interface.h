Interface

  Function fget_number_antigens() Result(ntotal) BIND(C,name="get_number_antigens_")
  Use ISO_C_BINDING, ONLY : C_INT
  Implicit None
  Integer(C_int):: ntotal
  End Function fget_number_antigens

  Function fget_antigen_types() Result(ntype) BIND(C,name="get_antigen_types_")
  Use ISO_C_BINDING, ONLY : C_INT
  Implicit None
  Integer(C_int):: ntype
  End Function fget_antigen_types
  
  Function fget_num_type_antigens() Result(numtype) BIND(C,name="get_num_type_antigens_")  
  Use ISO_C_BINDING, ONLY : C_PTR
  Implicit None
  Type(C_PTR) :: numtype
  End Function fget_num_type_antigens
  
  Function fget_radius_antigens() Result(antradius) BIND(C,name="get_radius_antigens_")  
  Use ISO_C_BINDING, ONLY : C_PTR
  Implicit None
  Type(C_PTR) :: antradius
  End Function fget_radius_antigens

  Function fget_length_antigens() Result(antlength) BIND(C,name="get_length_antigens_")  
  Use ISO_C_BINDING, ONLY : C_PTR
  Implicit None
  Type(C_PTR) :: antlength
  End Function fget_length_antigens
  
  Function fget_flexure_antigens() Result(kflex) BIND(C,name="get_flexure_antigens_")  
  Use ISO_C_BINDING, ONLY : C_PTR
  Implicit None
  Type(C_PTR) :: kflex
  End Function fget_flexure_antigens
  
  !Function fget_pattern_antigens() Result(pattern) BIND(C,name="get_pattern_antigens_")
  !Use ISO_C_BINDING, ONLY: C_STR
  !Implicit None
  !Type(C_STR) :: pattern
  !End Function fget_pattern_antigens
 
End Interface