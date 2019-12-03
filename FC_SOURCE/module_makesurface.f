!===========================================================================================================================
!===========================================================================================================================
!                                           MODULE TO MAKE A TILED SURFACE
!===========================================================================================================================
!===========================================================================================================================
       MODULE module_makesurface
       Integer :: chkflag,chkimg
       contains

!----------------------------------------------------------------------------------------------------------------------------
!                 SUBROUTINE CALLED FROM CPP CODE TO CONSTRUCT THE MEMBRANE SURFACE       
!----------------------------------------------------------------------------------------------------------------------------
       Subroutine makesurface()
       Use Module_datastruct ; Use Module_Curvcalc
       Use module_initialize_curvature_field ; Use Module_WriteData
       Implicit None
       Real(Kind =8 ):: per_czero_vert,shadow_size_fac
       Integer :: i,clock,startstate,num_czero_vert,annulus_start,annulus_end,czero_loci
       Character (Len=100) :: fname,statename,czero_pattern
       
       drzero=reshape((/zero,zero,zero/),(/3,1/))                                           ! Zero displacement vector for use in vertex move 
       Call SYSTEM_CLOCK(COUNT=clock); seed = -int(clock/8)*(processor_number+1)            !Initialize the Fortran seed from System clock
       
       num_bias_moves=100                                                                   ! number of bias moves for Rosenbluth Sampling of bound antigens
       
       If (debug_mode) Then
       	Print*
       	Print*
       	Print*,str_blue,"<START OF FORTRAN INITIALIZATION>",str_black
        Print*
       	Print*,str_red,"<MEMBRANE STRUCTURE>",str_black
       	Print*
	   Endif
	   
       Call Initialize_datastruct()
       
       If (((gsize)**2+4*(gsize)+2) .Gt. verno) Then
		Print*, 'Allocated vertices are smaller than required. Reallocate and start'
		Print*, 'The current setup can only handle a maximum grid size of  300X300'
		Print*, 'Change the value of verno in "module_datastruct.f" to the required number'
		Stop 
       Endif
       
       If ((2*(gsize)**2+4*(gsize)) .Gt. trino) Then
		Print*, 'Allocated Triangles are smaller than required. Reallocate and start'
		Print*, 'The current setup can only handle a maximum grid size of  300X300'
		Print*, 'Change the value of trino in "module_datastruct.f" to the required number'
		Stop 
       Endif
       
       If ((3*(gsize)**2+4*(gsize)-1) .Gt. linno) Then
		Print*, 'Allocated links are smaller than required. Reallocate and start'
		Print*, 'The current setup can only handle a maximum grid size of  300X300'
		Print*, 'Change the value of linno in "module_datastruct.f" to the required number'
		Stop 
       Endif
       
       fname="../PARAMETERS/membrane_parameters.inp"
       Open(10,File=fname)
       Read(10,*) kappa_unscaled                                                                                                    ! Original value of bending rigidity 
       Read(10,*) czero                                                                                                             ! value of Spontaneous curvature 
       Read(10,*) czero_pattern
       Read(10,*) per_czero_vert
       Read(10,*) annulus_start
       Read(10,*) annulus_end                                                                                                       ! pattern, number of vertices with czero, annulus_pattern_start, annulus_pattern_stop 
       Read(10,*) pressure                                                                                                          ! pressure 
       Read(10,*) memb_geom                                                                                                         ! geometry to initialize the membrane 
       Read(10,*) depth                           
       Read(10,*) period
       Read(10,*) blen_scale_factor                                                                                                 ! (lower,upper) bond len cutoff=(init_memb_blen-scale_factor),(\sqrt{3}init_memb_blen/scale_factor)
       Read(10,*) memb_restart_flag,startstate
       Read(10,*) antigen_restart_flag
       Read(10,*) vermov_step_adjust                                                                                                ! Time interval at which the vertex move size is adjusted 
       Read(10,*) fixed_frame_flag
       Read(10,*) shadow_size_fac                                                                                                   ! The shadow size is shadow_size_fac*size of the nanoparticle
       Close(10)
       
       kappa=kappa_unscaled/2.0                                                                                                     ! kappa is kappa_unscaled/2 
       Call Parse_antigen_specs()
       
       If(abs(pressure) .Gt. 0.0) Then
       	fixed_frame_flag = .True.
       	Print*,'Resetting frame flag to True since pressure is set to be  non-zero'
       Endif
       
       If (debug_mode) Then
        Print* 
        Print*,str_red, '<C++ - FORTRAN>',str_black
        Print*
        Print*,tab,tab,tab,'--> Avg. Bond length : ',blen
        Print*,tab,tab,tab,'--> Size of the grid : ',gsize
        Print*,tab,tab,tab,'--> Number of nanocarriers : ',num_nanocarrier
        Print*,tab,tab,tab,'--> Length of Periodic Box : ',periodic_box_length
        Print*,tab,tab,tab,'--> Height of Periodic Box : ',periodic_box_height
        Print*,tab,tab,tab,'--> Inverse Temperature : ', beta
        Print*
       Endif
       
       If(memb_restart_flag)Then
         Write(statename,*) startstate 
         Call Read_Membrane_Dump(startstate)
         num_czero_vert = NINT(per_czero_vert*nver) 
         czero_loci = find_middlevertex()
         Call Initialize_Spontaneous_curvature(czero_loci,num_czero_vert,czero,czero_pattern,annulus_start,annulus_end)
       Else
         If (debug_mode) Then       
          Print*       
          Print*,str_red,"<MEMBRANE PARAMETERS>",str_black
          Print*,tab,tab,tab,"--> Bending Rigidity (k_BT):",kappa_unscaled
          Print*,tab,tab,tab,"--> Grid size (units): ",gsize
          Print*,tab,tab,tab,"--> System Geometry (initial): ",Trim(Adjustl(memb_geom))
          If(Trim(Adjustl(memb_geom)).Eq.'SINUSOIDAL')Then
          Print*,tab,tab,tab,"--> Amplitude of Sinusoidal: ",depth
          Print*,tab,tab,tab,"--> Period of Sinusoidal:",period
          Endif
          Print*,tab,tab,tab,"--> Fortran Seed Initialized:",seed
          Print*
          PRint*,endmarker
         Endif
         
         If(Trim(Adjustl(memb_geom)) .Eq. 'PLANAR')Then
          Call Make_Planar_Membrane()
         Else If (Trim(Adjustl(memb_geom)) .Eq. 'SINUSOIDAL')Then
          Call Make_sinusoidal_membrane()
         Endif
         
         If (debug_mode) Then 
          Print*
          Print*,str_red,"<MESH DETAILS>",str_black
          Print*,"--> Number of vertices (Original,Phantom)",nver,pbnum
          Print*,"--> Number of triangles",ntr
          Print*,"--> Number of links",tlink
          Print*
          Print*,endmarker
         Endif
         
         Allocate(mp%vertex(pbnum)) ; Allocate(mp%triangle(ntr))
         Allocate(mp%link(-tlink:tlink))
         ver(:)%cur1=0.0 ; ver(:)%cur2=0.0 ; ver(:)%mcur=0.0
         antig(:)%disp1=0.0 ; antig(:)%disp2=0.0 ; antig(:)%diffus_tri=0
         Call Initialize_curvcalc() 
         
         vertice_call: Do i=1,nver                                                                                                    !Normal calculation over each vertex 
          Call normalcalc(i)
         Enddo vertice_call
         
         Do i=nver+1,pbnum,1                                                                                                          ! Initialize the periodic images 
          call mapimagetovertex(i)
         Enddo
         num_czero_vert=NINT(per_czero_vert*nver)        	
         
         If ( (num_czero_vert .Gt. 0) .And. (abs(czero) .Gt. 0.0)) Then
          czero_loci = find_middlevertex()
          If (debug_mode) Then
           Print*,str_green
           Print*,'=========================================================='
           Print*,' Spontaneous Curvature Initialization '
           Print*,'=========================================================='
           Print*,czero,czero_pattern
           Print*,"Pattern Initialized : ",Trim(adjustl(czero_pattern))
           Print*,"Magnitude of Spontaneous curvature : ",czero
           Print*,"Pattern centered at vertex : ",czero_loci
           If (Trim(Adjustl(czero_pattern)) .Eq. 'ANNULUS') Then
            Print*,"Initializing an annulus pattern between rings ",annulus_start, " and ",annulus_end
           Else
            Print*,"Number of vertices with czero :",num_czero_vert
           Endif
          Endif
          Call Initialize_Spontaneous_curvature(czero_loci,num_czero_vert,czero,czero_pattern,annulus_start,annulus_end)
         Endif
         
        blLcut=(blen-blen_scale_factor)**2 ; blUcut=3.0*blLcut                                                                       ! Fix the smallest bond length and then fix the largest bond length from there.
        If (debug_mode)  Print*,'Lower and upper bond length',blLcut,blUcut                                                          ! (1.2*blLcut)**0.5=blen, blUcut=sqrt(3.0)*blLcut**0.5 (a_0<1.2a_0<sqrt(3.0)a_0) 
          If(blUcut .Le. (2*blen**2)) Then
           Print*,str_red
           Print*, 'Max bondlength is less than the diagonal size of the  input mesh'
           Print*, 'Exiting: restart with correct values'
           Print*,str_black
           Stop
          Else If(blLcut .Gt. blen**2) Then
          Print*,str_red
          Print*, 'Min bondlength is larger than the bond size in the  input mesh'
          Print*, 'Exiting: restart with correct values'
          Print*,str_black
          Stop
         Endif
        Endif
        
      Call Initialize_Antigens()
      
      
       minshadow_size = 10                               ! Minimum number of vertices that should under the shadow ( we choose this to be one plaquette)
       mcsmovecounter=0
       
       If (num_nanocarrier .Gt. maxnc) Then
        Print*,str_red,' Number of nanocarrier is larger than the  maximum'
        Print*,str_blue,'Reset varible maxnc in module_datastruct.f'
        Print*,'exiting'
        Stop
       Endif
       
       
       Call set_selfavoidance_distances()
       Call Create_periodicbox()
       Call Membrane_Dump1(0)
       
       vermov_step_size=blen/10.0 ; vermov_attempt=0 ; vermov_accepted=0
       If (debug_mode) Then
       	Print*
       	Print*,tab,tab,tab,'--> Initial vertex move step size ',vermov_step_size
       	Print*
       	Print*,endmarker
       	Print*
       	Print*,str_blue,"<END OF FORTRAN INITIALIZATION>",str_black
	   Endif
       		
       
      Call check_triangle_in_pbbox()
      Call Check_vert_in_pbbox()
      Call compute_total_volume()
     End Subroutine makesurface
    
!===========================================================================================================================
!                 Compute the various self avoidance distances
!===========================================================================================================================     
	Subroutine set_selfavoidance_distances()
	   Use module_datastruct
	   Implicit None
	   Integer :: i, j
	   Allocate(sad%antant(num_antigens,num_antigens)); sad%antant = 0.0
	   Allocate(sad%antves(num_antigens,num_nanocarrier)); sad%antves = 0.0
	   Allocate(sad%vesmem(num_nanocarrier)); sad%vesmem = 0.0
	   
	   Do i = 1, num_antigens,1
	    Do j = 1, num_antigens,1
	    sad%antant(i,j) = (antig(i)%radius + antig(j)%radius)**2
	    Enddo
	   Enddo
	   
	      
	   Do i = 1, num_antigens, 1
	    Do j = 1, num_nanocarrier,1
	      sad%antves(i,j) = (nc_f(j)%radius)**2
	    Enddo
	   Enddo
	   
	   Do i=1,num_nanocarrier,1
	    sad%vesmem(i) = (nc_f(i)%soft_radius)**2
	   Enddo
	     
   End Subroutine set_selfavoidance_distances
  
!===========================================================================================================================
!                 Parse the specification for the antigens
!===========================================================================================================================     
       Subroutine Parse_antigen_specs()
       Use module_datastruct; Use, Intrinsic :: ISO_C_BINDING
       Implicit None
       Include '_mod_antigen_interface.h'
       Integer, Pointer  :: tempint(:)
       Real(Kind=8), Pointer  :: tempdouble(:)
       
        antspec%ntotal = fget_number_antigens()
        antspec%ntype = fget_antigen_types()
        Allocate(antspec%num_ntype(antspec%ntype))
        Call C_F_Pointer(fget_num_type_antigens(),tempint,[1])
        antspec%num_ntype(:) = tempint(1:antspec%ntype)
         
        Allocate(antspec%radius(antspec%ntype))
        Call C_F_Pointer(fget_radius_antigens(),tempdouble,[1])
        antspec%radius(:) = tempdouble(1:antspec%ntype)
         
        Allocate(antspec%length(antspec%ntype))
        Call C_F_Pointer(fget_length_antigens(),tempdouble,[1])
        antspec%length(:) = tempdouble(1:antspec%ntype)
         
        Allocate(antspec%kflex(antspec%ntype))
        Call C_F_Pointer(fget_flexure_antigens(),tempdouble,[1])
        antspec%kflex(:) = tempdouble(1:antspec%ntype)
         
        antspec%pattern = 'random'   ! Need to find a way to parse string from C++
         
        If (debug_mode) Then
        Print*,str_red,'<ANTIGEN INITIALIZATION IN FORTRAN>',str_black
        Print*,tab,tab,tab,'--> Number of Antigens = ',antspec%ntotal
        Print*,tab,tab,tab,'--> Number of Antigen Types  = ',antspec%ntype
        Print*,tab,tab,tab,'--> Number of Antigen in each type = ', antspec%num_ntype
        Print*,tab,tab,tab,'--> Radius of Antigen types = ', antspec%radius
        Print*,tab,tab,tab,'--> Length of Antigen types = ', antspec%length
        Print*,tab,tab,tab,'--> Flexure of Antigen types = ', antspec%kflex
        Print*,tab,tab,tab,'--> Antigens are arranged in ',antspec%pattern
        Print*,tab,tab,tab,'--> Initializing Antigens on the surface'
        Print*
        Print*,endmarker
        Endif
       
       End Subroutine Parse_antigen_specs
       
!===========================================================================================================================
!                 Function to Chose the middle vertex in a given membrane data
!===========================================================================================================================

       Function find_middlevertex()Result(midver)
        USE module_datastruct
        IMPLICIT NONE
        Integer :: midver,i1(1)
        Real(Kind=8), Allocatable :: errx(:),erry(:),err(:)
        Real(Kind=8) :: minx,maxx,midx,miny,maxy,midy
        Allocate(errx(pbnum)) ;  Allocate(erry(pbnum)) ;  Allocate(err(pbnum))
        minx=minval(ver(:)%vcoord(1,1)) ; maxx=maxval(ver(:)%vcoord(1,1)) ; midx=(minx+maxx)/2
        miny=minval(ver(:)%vcoord(2,1)) ; maxy=maxval(ver(:)%vcoord(2,1)) ; midy=(miny+maxy)/2
        errx=abs(ver(:)%vcoord(1,1)-midx); erry=abs(ver(:)%vcoord(2,1)-midy) ;	err=errx+erry
        i1=minloc(err)
        midver=i1(1)
       End Function find_middlevertex       

!------------------------------------------------------------------------------------------------------------------------
!                                    SUBROUTINE TO READ THE STATE OF THE SYSTEM
!------------------------------------------------------------------------------------------------------------------------
       Subroutine Read_membrane_dump(time)
       Use module_datastruct ; Use module_curvcalc
       Integer :: time,i,trinum
       Character(50):: filename,tname
       Character(Len=100) :: marker_string
       Write(tname,*) time
       filename='../SYS_STATE/RESTART/membrane_state-'//Trim(Adjustl(tname))//'.dump'
       If (debug_mode) Print*,str_red, '<MEMBRANE DUMP>   ',filename,str_black
       Open(01,File=filename,Form='Formatted')

       tri(:)%nantigen=0
       Read(01,*)blen,blLcut,blUcut                                                                                          	    ! blen, blLcut,blUcut
       Read(01,*)nver,pbnum,ntr,tlink
       If (debug_mode) Then
       Print*
       Print*,str_red,"<MESH DETAILS>",str_black
       Print*,tab,tab,tab,"--> (nver, pbnum): ",nver,pbnum
       Print*,tab,tab,tab,"--> (ntr): ",ntr
       Print*,tab,tab,tab,"--> (tlink): ",tlink
       Print*
       Print*,endmarker
       Endif

        Read(01,*)marker_string
		If(Trim(Adjustl(marker_string)) .Eq. "---BeginVertexData")Then
    	    Do i=1,pbnum,1
        	Read(01,*)ver(i)%vnor
        	Read(01,*)ver(i)%vcoord
        	Read(01,*)ver(i)%vneipt
        	Read(01,*)ver(i)%vneitr
        	Read(01,*)ver(i)%PBCver
        	Read(01,*)ver(i)%L2G(:,1),ver(i)%L2G(:,2),ver(i)%L2G(:,3)
        	Read(01,*)ver(i)%HHM(:,1),ver(i)%HHM(:,2),ver(i)%HHM(:,3)
        	Read(01,*)ver(i)%nonei,ver(i)%boundary,ver(i)%pbimno,ver(i)%pbmap,ver(i)%imver,ver(i)%nover,ver(i)%boxvert
			Read(01,*)ver(i)%mcur,ver(i)%cur1,ver(i)%cur2,ver(i)%totarea
    	    Read(01,*)ver(i)%czero_flag,ver(i)%czero
        	Read(01,*)ver(i)%cellno
        	Enddo
		Endif
		Read(01,*)marker_string
        Allocate(mp%vertex(pbnum)) 

       Read(01,*)marker_string
       If(Trim(Adjustl(marker_string)).Eq. "---BeginAntigenData")Then	
		If (antigen_restart_flag .Eqv. .True.) Then
         If(antspec%ntotal .Gt.0)Then                                                                                 ! Neglect Antigen data if num_antigens in set to zero
 	  		Do i=1,antspec%ntotal,1
	   		Read(01,*)antig(i)%vertex
           	Read(01,*)antig(i)%base_coord,antig(i)%tip_coord
			Read(01,*)antig(i)%diffus_tri,antig(i)%disp1,antig(i)%disp2
	   		trinum=antig(i)%diffus_tri
	   		If(trinum .Gt. 0) Then                                  ! Antigen is on the vertex
			tri(trinum)%nantigen=tri(trinum)%nantigen+1
	    	tri(trinum)%antigen_list(tri(trinum)%nantigen)=i       ! connect the antigen to the triangle
	   		Endif
	   		ver(antig(i)%vertex)%antigen_flag=1
	   		ver(antig(i)%vertex)%nantigen=ver(antig(i)%vertex)%nantigen+1                       ! increment the number of antigen connected to the vertex by 1 
	   		ver(antig(i)%vertex)%antigen_list(ver(antig(i)%vertex)%nantigen)=i                      ! append 'i' to the antigen list of antig(i)%vertex
	   		Read(01,*)antig(i)%theta,antig(i)%phi
	  		Enddo
 	  		Allocate(mp%antig(antspec%ntotal))
 	  		Read(01,*)marker_string
		Else
	  		Do While(Trim(Adjustl(marker_string)).Ne."---EndofAntigenData")
 	   		Read(01,*)marker_string
	  		Enddo	
	 	Endif
	   Else	
	 	Read(01,*)marker_string
	 	Do While(Trim(Adjustl(marker_string)).Ne."---EndofAntigenData")
	  	Read(01,*)marker_string
	 	Enddo	
	  Endif
	 Endif

       Read(01,*)marker_string
       If(Trim(Adjustl(marker_string)).Eq. "---BeginTriangleData")Then
        Do i=1,ntr,1
         Read(01,*)tri(i)%ar,tri(i)%pbflag,tri(i)%boxtriangle,tri(i)%li,tri(i)%vert,tri(i)%fnor
        Enddo 
       Endif
       Allocate(mp%triangle(ntr))
       Read(01,*)marker_string
       Read(01,*)marker_string

       If(Trim(Adjustl(marker_string)).Eq. "---BeginLinkData")Then
        Do i=1,tlink
         Read(01,*)lin(i)%tr,lin(-i)%tr,lin(i)%boundary,lin(i)%pbflag,lin(i)%sep
         lin(-i)%boundary=lin(i)%boundary
         lin(-i)%sep(1)=lin(i)%sep(2)
         lin(-i)%sep(2)=lin(i)%sep(1)
         lin(-i)%pbflag=lin(i)%pbflag
        Enddo
       Endif
       Read(01,*)marker_string
       Close(01)

       Allocate(mp%link(-tlink:tlink))

       Call Initialize_curvcalc()        

       vertice_call: Do i=1,nver                                                                                                    ! Normal calculation over each vertex 
        Call normalcalc(i)
        Enddo vertice_call

	triangle_vol: Do i=1,ntr,1                                                                                           	    ! Initialize the volume associated with each triangle
	Call areacalc(i)
	Call  compute_tri_volume(i)
	Enddo triangle_vol
        
	Do i=nver+1,pbnum,1                                                                                                         ! Initialize the periodic images 
	call mapimagetovertex(i)
	Enddo

	If(antigen_restart_flag .Eqv. .False.) Then
	  If (debug_mode) Then
      	Print*,str_red,'<ANTIGEN INTIALIZATION>',str_black
      	Print*,tab,tab,tab,'--> Initializing Antigens randomly on the surface'
      	PRint*
      	Print*,endmarker
      Endif
	 Call Initialize_Antigens()
	Endif

    End Subroutine Read_membrane_dump      
!---------------------------------------------------------------------------------------------------------------------------
!                            SUBROUTINE TO INITIALIZE THE ANTIGENS
!---------------------------------------------------------------------------------------------------------------------------
    Subroutine Initialize_Antigens()
       Use module_datastruct ; Use module_randomnumber
       Implicit None
       Integer :: i,i1,ranno
       
       If (debug_mode) Then
       Print*
       Print*,str_red,"<RANDOM ANTIGEN INITIALIZATION>",str_black
       Print*
       Endif
        
       Do i=1,pbnum,1
       ver(i)%antigen_list=0
       ver(i)%nantigen=0
       ver(i)%antigen_flag=0                                                                ! Set antigen flag to zero at all vertices to start with
       Enddo

       Do i=1,ntr,1
       tri(i)%nantigen=0
       tri(i)%antigen_list=0
       Enddo

       antig(:)%vertex=0 ; antig(:)%base_coord(1,1)=0.0
       antig(:)%base_coord(2,1)=0.0 ; antig(:)%base_coord(3,1)=0.0
       antig(:)%tip_coord(1,1)=0.0 ; antig(:)%tip_coord(2,1)=0.0
       antig(:)%tip_coord(3,1)=0.0 ; antig(:)%diffus_tri=0
       
       If (.Not. Allocated(mp%antig)) Allocate(mp%antig(antspec%ntotal))                     ! Allocate mp%antig if not allocated already
       
       ! Set the types of antigens
       antig(1:antspec%num_ntype(1))%ant_type = 0                                             ! we start with 0 for c++ (care needed when handling fortran)
       Do i=2,antspec%ntype,1
        antig(Sum(antspec%num_ntype(1:i-1))+1: Sum(antspec%num_ntype(1:i)))%ant_type = i-1
       Enddo
        
       ! Assign the various parameters
       Do i=1,antspec%ntotal,1
        i1 = antig(i)%ant_type + 1
        antig(i)%length = antspec%length(i1)
        antig(i)%radius = antspec%radius(i1)
        antig(i)%kflex = antspec%kflex(i1)
       Enddo
        
    ! If the antigens are to be randomly distributed
      If (Trim(Adjustl(antspec%pattern)) .Eq. 'random') Then
       i=0
       Do While(i .Lt. antspec%ntotal)
        ranno = Nint(nver*ran2(Seed))+1                                                                                               ! Choose a random vertex to place the antigens 
        If((ranno .Lt. nver) .And. (ver(ranno)%nantigen .Eq. 0)) Then                                                                 ! Place one antigen per vertex in the initialization step 
         i=i+1 ; ver(ranno)%nantigen = ver(ranno)%nantigen + 1
         ver(ranno)%antigen_flag = 1
         ver(ranno)%antigen_list(ver(ranno)%nantigen) = i
         antig(i)%vertex = ranno                                                                                                       ! antigen i is associated to vertex rand
         antig(i)%base_coord = ver(ranno)%vcoord                                                                                       ! Shares the base coordinates 
         antig(i)%tip_coord = antig(i)%base_coord + antig(i)%length*ver(ranno)%vnor                                                    ! Coordinates for the tip
        Endif 
       Enddo
      Endif
      
      If (debug_mode) Then
      	Print*,tab,tab,tab,tab,'--> Completed Initialization of Antigens ',antspec%ntotal
      	Print*
      	Print*,endmarker
	  Endif
	  
    End Subroutine Initialize_Antigens
!---------------------------------------------------------------------------------------------------------------------------
!                            SUBROUTINE TO INITIALIZE THE DATASTRUCTURE
!---------------------------------------------------------------------------------------------------------------------------
       Subroutine Initialize_datastruct()
       Use module_datastruct
       Implicit None
       Integer :: i
       Do i=1,pbnum,1
       ver(i)%vcoord=0.0 ; ver(i)%vneipt=0 ;ver(i)%vneitr=0 ;ver(i)%pbmap=0
       ver(i)%vnor=0.0 ; ver(i)%pbimno=0 ; ver(i)%L2G=0.0 ;ver(i)%HHM=0.0
       ver(i)%pbmap=0 ; ver(i)%antigen_list=0 ; ver(i)%nantigen=0
       ver(i)%antigen_flag=0 ; ver(i)%boxvert=0; ver(i)%shadownc=0
       Enddo

       Do i=1,ntr,1
       tri(i)%nantigen=0 ;  tri(i)%antigen_list=0
       tri(i)%vert=0 ; tri(i)%li=0 ; tri(i)%fnor=0.0
       tri(i)%boxtriangle=1 ; tri(i)%vol=0.0
       Enddo

       antig(:)%vertex=0 ; antig(:)%base_coord(1,1)=0
       antig(:)%base_coord(2,1)=0 ; antig(:)%base_coord(3,1)=0
       antig(:)%tip_coord(1,1)=0 ; antig(:)%tip_coord(2,1)=0
       antig(:)%tip_coord(3,1)=0 ; antig(:)%diffus_tri=0
       End Subroutine Initialize_datastruct
      
!---------------------------------------------------------------------------------------------------------------------------
!                            SUBROUTINE TO MAKE A PLANAR SURFACE
!---------------------------------------------------------------------------------------------------------------------------
      Subroutine Make_Planar_Membrane()
      USE module_datastruct
      IMPLICIT NONE
      Integer :: row,col,vno,v,v1,v2,lno,trn,i
      Integer :: nv1,nv2,nntr,tm1,tm2,tmv(3),prchk,nechk,j
      Integer :: tvno,flagallot,linflag,vflag1,vflag2
      Real(KIND=8) :: pbcshift,pbx,pby,pbz
      
      Print*,str_red,"<INITIALIZING A PLANAR MEMBRANE>",str_black
               
      lin(:)%sep(1)=0 ; lin(:)%sep(2)=0 ; lin(:)%pbflag=0                                                                           ! set all link details to zero
      lin(:)%tr=0     ; lin(:)%boundary=0 ; pbnum=0
      ver(:)%nover=0  ; ver(:)%pbimno=0
      
      vno=0 ; tlink=4000 ; pbcshift=gsize*blen                                                                                      ! Initializing the position 
      row_loop : Do row=1,gsize,1 
      col_loop: Do col=1,gsize,1
      vno=vno+1
      ver(vno)%vcoord= reshape((/row*blen,col*blen,membrane_init_z/),(/3,1/))                                               ! Position of vertex vno 
      Enddo col_loop
      Enddo row_loop 
      nver=vno ; pbnum=vno                                                                                                          !Periodic images will have index greater than the number of vertex  

      If (debug_mode) Print*,tab,tab,tab,"--> Membrane Z location",membrane_init_z

        find_neigh_pt: Do  v=1,vno,1                                    
          ver(v)%nonei=0 ; ver(v)%vneipt=0 ; ver(v)%boundary=0
          ver(v)%PBCver=0 
          c1:If(v.LE.gsize)Then                                                                                                     ! upper row 
           upper_row:If(v.EQ.1)Then                                                                                                 ! Left upper corner
              ver(v)%boundary=1
              ver(v)%nonei=6
              ver(v)%PBCver(1:6)=(/0,0,0,1,1,1/)
              ver(v)%vneipt(1:6)=(/gsize+v,gsize+v+1,v+1,gsize*(gsize-1)+1,gsize**2,gsize/)                                         ! Vertex 1 with periodic boundary

                tvno=gsize*(gsize-1)+1
                  pbx=ver(tvno)%vcoord(1,1)-pbcshift
                  pby=ver(tvno)%vcoord(2,1)                                                                                         ! New coordinates of the boundary image 
                  pbz=ver(tvno)%vcoord(3,1)
                  If(ver(tvno)%pbimno.Eq.0)Then
                  flagallot=0
                  Else
                  flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(4)=chkimg
                   Endif
                  Endif 
                  If(flagallot.Eq. 0) Then                                                                                          ! The periodic images are stored as seperate variables      
                  call allotposition(pbx,pby,pbz,tvno,v,4)
                  Endif

               tvno=gsize**2
                 pbx=ver(tvno)%vcoord(1,1)-pbcshift
                 pby=ver(tvno)%vcoord(2,1)-pbcshift                                                                                 ! New coordinates of the boundary image 
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else               
                 flagallot=1
                  call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(5)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then                                                                                           ! The periodic images are stored as seperate variables      
                 call allotposition(pbx,pby,pbz,tvno,v,5)
                 Endif
         
               tvno=gsize
                 pbx=ver(tvno)%vcoord(1,1) 
                 pby=ver(tvno)%vcoord(2,1)-pbcshift                     
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(6)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,6)
                 Endif            
                                    
             Else If(v.EQ.gsize)Then                                                                                                ! Right upper corner 
                ver(v)%boundary=1
                ver(v)%nonei=6 
                ver(v)%PBCver(1:6)=(/0,1,1,1,1,0/)
                ver(v)%vneipt(1:6)=(/v+gsize,v+1,1,gsize**2,gsize**2-1,v-1/)

               tvno=v+1
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)+pbcshift                
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                 flagallot=1
                 call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(2)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,2)
                 Endif            

               tvno=1
                 pbx=ver(tvno)%vcoord(1,1) 
                 pby=ver(tvno)%vcoord(2,1)+pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(3)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,3)
                 Endif            

               tvno=gsize**2
                 pbx=ver(tvno)%vcoord(1,1)-pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                  flagallot=1
                  call ifalreadyalloted(pbx,pby,pbz,tvno)
                  If(chkflag.Eq.0)Then
                  flagallot=0
                  Else
                  ver(v)%vneipt(4)=chkimg
                  Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,4)
                 Endif            

               tvno=gsize**2-1
                 pbx=ver(tvno)%vcoord(1,1)-pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(5)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,5)
                 Endif            
                     
             Else
               ver(v)%boundary=1
               ver(v)%nonei=6  
               ver(v)%PBCver(1:6)=(/0,0,0,1,1,0/)
               ver(v)%vneipt(1:6)=(/v+gsize,v+gsize+1,v+1,gsize*(gsize-1)+v, gsize*(gsize-1)+v-1,v-1/)                              ! The rest of the points on the top row 

               tvno=gsize*(gsize-1)+v
                  
                 pbx=ver(tvno)%vcoord(1,1)-pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(4)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,4)
                 Endif            

               tvno=gsize*(gsize-1)+v-1
                 pbx=ver(tvno)%vcoord(1,1)-pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(5)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,5)
                 Endif            

               Endif upper_row



            Else If(v.GT.(vno-gsize))Then                                                                                           ! Lower row 
 
               lower_row:If(v.EQ.(vno-gsize+1))Then                                                                                 ! Vertex at left bottom corner                             
               ver(v)%boundary=1                                                                                                    ! Vertices with spl=1 are those on the boundary            
               ver(v)%nonei=6
               ver(v)%PBCver(1:6)=(/1,1,0,0,1,1/)
               ver(v)%vneipt(1:6)=(/1,2,v+1,v-gsize,v-1,gsize**2/)

               tvno=1    
                 pbx=ver(tvno)%vcoord(1,1)+pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(1)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,1)
                 Endif            

               tvno=2
                 pbx=ver(tvno)%vcoord(1,1)+pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(2)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,2)
                 Endif        
 
               tvno=v-1
                 pbx=ver(tvno)%vcoord(1,1) 
                 pby=ver(tvno)%vcoord(2,1)-pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1                       
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(5)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 lin(:)%sep(1)=0 ; lin(:)%sep(2)=0 ; lin(:)%pbflag=0                                                                           ! set all link details to zero
                 call allotposition(pbx,pby,pbz,tvno,v,5)
                 Endif         

               tvno=gsize**2
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)-pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1                       
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(6)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,6)
                 Endif         


             Else If(v.EQ.vno)Then                                                                                                  ! Right bottom corner 
               ver(v)%boundary=1
               ver(v)%nonei=6
               ver(v)%PBCver(1:6)=(/1,1,1,0,0,0/)
               ver(v)%vneipt(1:6)=(/gsize,1,v-gsize+1,v-gsize,v-gsize-1,v-1/)

               tvno=gsize
                 pbx=ver(tvno)%vcoord(1,1)+pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1                    
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(1)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,1)
                 Endif         

               tvno=1
                 pbx=ver(tvno)%vcoord(1,1)+pbcshift
                 pby=ver(tvno)%vcoord(2,1)+pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1                     
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(2)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,2)
                 Endif         

               tvno=v-gsize+1
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)+pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1                     
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(3)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,3)
                 Endif         
                 
             Else 
               ver(v)%boundary=1
               ver(v)%nonei=6
               ver(v)%PBCver(1:6)=(/1,1,0,0,0,0/)
               ver(v)%vneipt(1:6)=(/v-gsize*(gsize-1),v+1-gsize*(gsize-1),v+1,v-gsize,v-gsize-1,v-1/)

               tvno=v-gsize*(gsize-1)
                 pbx=ver(tvno)%vcoord(1,1)+pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=ver(tvno)%vcoord(3,1)                       
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(1)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,1)
                 Endif         
             
               tvno=v-gsize*(gsize-1)+1
                 pbx=ver(tvno)%vcoord(1,1)+pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(2)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,2)
                 Endif         

               Endif lower_row

           Else If(mod(v,gsize).EQ.1)Then                                                                                           ! left column
               ver(v)%boundary=1
               ver(v)%nonei=6
               ver(v)%PBCver(1:6)=(/0,0,0,0,1,1/)
               ver(v)%vneipt(1:6)=(/v+gsize,v+gsize+1,v+1,v-gsize,v-1,v-1+gsize/)

               tvno=v-1
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)-pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1                     
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(5)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,5)
                 Endif         

               tvno=v-1+gsize
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)-pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                  flagallot=1                
                  call ifalreadyalloted(pbx,pby,pbz,tvno)
                  If(chkflag.Eq.0)Then
                  flagallot=0
                  Else
                  ver(v)%vneipt(6)=chkimg
                  Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,6)
                 Endif         

              Else If(mod(v,gsize).EQ.0)Then                                                                                         ! Right column 
               ver(v)%boundary=1
               ver(v)%nonei=6
               ver(v)%PBCver(1:6)=(/0,1,1,0,0,0/)
               ver(v)%vneipt(1:6)=(/v+gsize,v+1,v-gsize+1,v-gsize,v-gsize-1,v-1/)

               tvno=v+1
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)+pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1                      
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(2)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,2)
                 Endif         

               tvno=v+1-gsize
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)+pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1     
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(3)=chkimg
                   Endif
                 Endif 
                 If(flagallot.Eq. 0) Then
                 lin(:)%sep(1)=0 ; lin(:)%sep(2)=0 ; lin(:)%pbflag=0                                                                           ! set all link details to zero
                 call allotposition(pbx,pby,pbz,tvno,v,3)
                 Endif         

               Else
                 ver(v)%nonei=6
                 ver(v)%PBCver(1:6)=(/0,0,0,0,0,0/)
                 ver(v)%vneipt(1:6)=(/v+gsize,v+gsize+1,v+1,v-gsize,v-gsize-1,v-1/) 
          Endif c1        
        Enddo find_neigh_pt



! So far we have made a square lattice and had defined its neighbours inclusive of a PBC in x and y directions

! The above part may look longer but it is a simple by hand implementation of the PBC. PBC for connected
! triangulated surface with bond length constraints are a bit complicated and to take care of this, I in
! this code, maintain a seperate PBC coordinate corresponding to each vertex and update it too when the 
! vertex under question is disturbed

! ver(v)%PBCver --> a value of 1 here indicates that this position is occupied by PB vertex. So use 
! the coordinates stored in ver()%vcoord
       
! The links from each vertex is alloted and the boundary condition works
! fine as far as the neighbour list is concerned.

! The periodic vertices are marked by another type structre called
! pbvertex

! Start  from alloting the links and making the triangles

! How to allot the links. Should we allot the links that includes the
! periodic vertices too ?  In princple the links between the boundaary
! nodes are to preserved in neighbourhood. So the details of the links
! that arising from PBC are not essential and so are those for such
! triangles


      lno=0
       Do v=1,pbnum,1
         Do v1=1,ver(v)%nonei,1        
         prchk=0  ; linflag=0                                                                                                       ! To check if this link is already alloted ?                   
         nv1=ver(v)%vneipt(v1) ; i=-lno                                                                                           ! v1 th neighbour of v ; linno --> size of array for link list 
           Do WHILE(i.LE.lno .And. prchk.EQ.0)                                                                                      ! check over all the link list                       
           If(i.NE.0)Then                                                                                                           ! There is no link named zero since every l has a -l 

           If((lin(i)%sep(1).Eq.v) .And. (lin(i)%sep(2).Eq.nv1) .And. (linflag .Eq. lin(i)%pbflag))Then
              prchk=1 ; Endif
           If((lin(i)%sep(2).EQ.v) .And. (lin(i)%sep(1).EQ.nv1) .And. (linflag .Eq. lin(i)%pbflag)) Then
              prchk=1 ; Endif ; Endif
           i=i+1
           Enddo

           If (prchk.EQ.0)Then                                                                                                      ! If no existing link has been found then 
           lno=lno+1
           lin(lno)%sep=(/v,nv1/) 
           lin(-lno)%sep=(/nv1,v/)
           If((ver(v)%boundary.EQ.1).And.(ver(nv1)%boundary.EQ.1))Then                                                              ! Mark all the boundary links 
           lin(lno)%boundary=1 ; lin(-lno)%boundary=1
           Endif  ; Endif
       Enddo ; Enddo
       tlink=lno                                                                                                                    ! Total number of links in the system 

! Link allocation was easy making the triangles seems to be difficult.
! Another idea that can be thought of is to remove the periodic vertex
! structure and make it a part of the vertex structure by giving it a
! seperate flag.

! vertex number > nver is a marker for the periodic boundaries          

        trn=0  
        tri(:)%pbflag=0                                                                                                             ! Set all peridoic boundary flag to zero
        Do v=1,nver,1
             nntr=0 
             Do v1=1,ver(v)%nonei,1                                                                                                 ! Run over the neighbourhood of v                                     
             prchk=0 ;nechk=0                                            ! prchk --> check if already alloted                                  
             nv1=ver(v)%vneipt(v1) ; nv2=ver(v)%vneipt(v1+1)             ! the i th and i+1 th neighbour, be careful Periodic images can occur 
             If(v1.EQ.ver(v)%nonei) nv2=ver(v)%vneipt(1)                 ! For all submarginal vertices there exists a circular boundary       
                                                                         ! check with all neig triangles of nv1 for previous allotment         
             tm1=1       
             Do WHILE(tm1.LE.ver(nv1)%nonei.And.prchk.EQ.0)              ! For each neigh triangle around vertex 'v', ^^^^^^^^^^^ 
             tm2=ver(nv1)%vneitr(tm1)                                    ! Pick up each triangle  as tm2                          
	     If (tm2 .Gt. 0) Then
             tmv=tri(tm2)%vert                                           ! vertices of tm2                                        
             If(v.EQ.tmv(1).OR.v.EQ.tmv(2).OR.v.EQ.tmv(3))Then           ! if vertcies of tm2 == (v,nv1,nv2) or (nv2,v,nv1)       
             If(nv1.EQ.tmv(1).OR.nv1.EQ.tmv(2).OR.nv1.EQ.tmv(3))Then     ! or (nv1,nv2,v)                                         
             If(nv2.EQ.tmv(1).OR.nv2.EQ.tmv(2).OR.nv2.EQ.tmv(3))Then
             nntr=nntr+1                                                 ! Number of nearest triangles list is increased by 1      
             ver(v)%vneitr(nntr)=tm2 ; prchk=1 ;nechk=1                  ! If present then vert --> trian link for vertex v is tm2 
	     Endif ; End If; End If ; Endif ;tm1=tm1+1;Enddo                       ! Exit the loop once presence is already detected .     

                                                                         ! check with all neig triangles of nv2 for previous allotment  
             tm1=1                                                       ! spl vertices have no circular boundaries                     
             Do WHILE(tm1.LE.ver(nv2)%nonei.And.prchk.EQ.0)              ! For each neigh triangle around vertex 'v', ^^^^^^^^^^^       
             tm2=ver(nv2)%vneitr(tm1)                                    ! Pick up each triangle  as tm2                                
	     If (tm2 .Gt. 0) Then
             tmv=tri(tm2)%vert                                           ! vertices of tm2                                              
             If(v.EQ.tmv(1).OR.v.EQ.tmv(2).OR.v.EQ.tmv(3))Then           ! if vertcies of tm2 == (v,nv1,nv2) or (nv2,v,nv1)             
             If(nv1.EQ.tmv(1).OR.nv1.EQ.tmv(2).OR.nv1.EQ.tmv(3))Then     ! or (nv1,nv2,v)                                               
             If(nv2.EQ.tmv(1).OR.nv2.EQ.tmv(2).OR.nv2.EQ.tmv(3))Then
             nntr=nntr+1                                                 ! Number of nearest triangles list is increased by 1      
             ver(v)%vneitr(nntr)=tm2 ; prchk=1 ;nechk=1                  ! If present then vert --> trian link for vertex v is tm2 
	     Endif ; End If; End If ; Endif
             tm1=tm1+1;Enddo   

             If(nechk.EQ.0)Then                                          ! If presence not detected a new triangle is formed  &&&&&&&&&&& 
              nntr=nntr+1 ; trn=trn+1                                    ! Both number of nearest triangle and total no is increased      
              ver(v)%vneitr(nntr)=trn                                   
              tri(trn)%vert=reshape((/v,nv1,nv2/),shape(tri(trn)%vert))  ! Vertices of the new traingle   

              If(v.Gt.nver .Or. nv1.Gt.nver .Or. nv2.Gt.nver)Then        ! $ seperately for excluding in rendering later 
              tri(trn)%pbflag=1
              Endif   

 
              tm1=-tlink;prchk=0                                         ! tm1 is a loop variable to see if a link already exists for the vertices 
              Do WHILE(tm1.LE.tlink .And. prchk.EQ.0)                    ! Check for the existence of the links                                    
              If((tm1).NE.0)Then                                         ! Link number zero is not present in the system                           
              If(lin(tm1)%sep(1).EQ.v.And.lin(tm1)%sep(2).EQ.nv1)Then    ! prchk --> variable to break out from the loop the link exists           
              tri(trn)%li(1)=tm1 ; prchk=1                               ! nechk --> flag to allot a new link if a link does not exists            
              lin(tm1)%tr=trn
              Endif ; End If  ;tm1=tm1+1 
              Enddo

              tm1=-tlink ; prchk=0 
              Do WHILE(tm1.LE.tlink .And. prchk.EQ.0)                    ! Check for the existence of the links now for nv1,nv2 as before
              If((tm1).NE.0)Then
              If(lin(tm1)%sep(1).EQ.nv1.And.lin(tm1)%sep(2).EQ.nv2)Then
              tri(trn)%li(2)=tm1 ; prchk=1  
              lin(tm1)%tr=trn
              Endif ; End If  ;tm1=tm1+1 
              Enddo

              tm1=-tlink ; prchk=0 
              Do WHILE(tm1.LE.tlink .And. prchk.EQ.0)                    ! Check for the existence of the links  between nv2 --> v
              If((tm1).NE.0)Then
              If(lin(tm1)%sep(1).EQ.nv2.And.lin(tm1)%sep(2).EQ.v)Then
              tri(trn)%li(3)=tm1 ; prchk=1 ;nechk=1 
              lin(tm1)%tr=trn
              Endif ; End If  ;tm1=tm1+1 
              Enddo

             Endif                                                       ! &&&&&&&&&&& 
            Enddo                                                        ! ^^^^^^^^^^^ 
        Enddo
        ntr=trn



           Do j=1,ntr,1                                                     ! Allot the neighbouring triangles to the periodic image boundaries
           v=tri(j)%vert(1)
           v1=tri(j)%vert(2)
           v2=tri(j)%vert(3)
  
              If(v>nver)Then
              ver(v)%nonei=ver(v)%nonei+1
              ver(v)%vneitr(ver(v)%nonei)=j
              Endif

              If(v1>nver)Then
              ver(v1)%nonei=ver(v1)%nonei+1
              ver(v1)%vneitr(ver(v1)%nonei)=j
              Endif

              If(v2>nver)Then
              ver(v2)%nonei=ver(v2)%nonei+1
              ver(v2)%vneitr(ver(v2)%nonei)=j
              Endif
           Enddo

           ver(nver+1:pbnum)%nonei=ver(nver+1:pbnum)%nonei+1             ! Since the edges do not have a circular boundary no of vertex=no of triang+1


           Do j=1,ntr,1                                                  ! Allot the neighbouring vertices to the boundaries (careful they do not $
           v=tri(j)%vert(1)                                              !$ follow any sequence as done by the inner nodes
           v1=tri(j)%vert(2)
           v2=tri(j)%vert(3)
  
            If(v>nver)Then
              vflag1=0 ; vflag2=0
              Do i=1,ver(v)%nonei,1
              If(ver(v)%vneipt(i).Eq.v1) vflag1=1
              If(ver(v)%vneipt(i).Eq.v2) vflag2=1
              Enddo

              If(vflag1.Eq.0) Then
              ver(v)%nover=ver(v)%nover+1
              ver(v)%vneipt(ver(v)%nover)=v1
              Endif  

              If(vflag2.Eq.0) Then
              ver(v)%nover=ver(v)%nover+1
              ver(v)%vneipt(ver(v)%nover)=v2
              Endif  
            Endif

            If(v1>nver)Then
              vflag1=0 ; vflag2=0
              Do i=1,ver(v1)%nonei,1
              If(ver(v1)%vneipt(i).Eq.v)  vflag1=1
              If(ver(v1)%vneipt(i).Eq.v2) vflag2=1
              Enddo

              If(vflag1.Eq.0) Then
              ver(v1)%nover=ver(v1)%nover+1
              ver(v1)%vneipt(ver(v1)%nover)=v
              Endif  

              If(vflag2.Eq.0) Then
              ver(v1)%nover=ver(v1)%nover+1
              ver(v1)%vneipt(ver(v1)%nover)=v2
              Endif
            Endif


           If(v2>nver)Then
              vflag1=0 ; vflag2=0
              Do i=1,ver(v2)%nonei,1
              If(ver(v2)%vneipt(i).Eq.v)  vflag1=1
              If(ver(v2)%vneipt(i).Eq.v1) vflag2=1
              Enddo

              If(vflag1.Eq.0) Then
              ver(v2)%nover=ver(v2)%nover+1
              ver(v2)%vneipt(ver(v2)%nover)=v
              Endif  

              If(vflag2.Eq.0) Then
              ver(v2)%nover=ver(v2)%nover+1
              ver(v2)%vneipt(ver(v2)%nover)=v1
              Endif
           Endif
       Enddo


        Do i=1,ntr,1                                                     ! Calculate the area of the triangles 
        call areacalc(i)
        Enddo 

        End Subroutine Make_Planar_Membrane

!-----------------------------------------------------------------------------------------------------------------------------------	
!                                   Subroutine  to make a sinusoidal membrane geometry
!-----------------------------------------------------------------------------------------------------------------------------------	
	Subroutine Make_Sinusoidal_Membrane()
        USE module_datastruct
        IMPLICIT NONE
        Integer :: row,col,vno,v,v1,v2,lno,trn,i,pos(3,1)
        Integer :: nv1,nv2,nntr,tm1,tm2,tmv(3),prchk,nechk,j
        Integer :: tvno,flagallot,linflag,vflag1,vflag2
        Real(KIND=8) ::pbcshift,pbx,pby,pbz,cx,cy,cz

        lin(:)%sep(1)=0 ; lin(:)%sep(2)=0 ; lin(:)%pbflag=0              ! set all link details to zero
        lin(:)%tr=0     ; lin(:)%boundary=0 ; pbnum=0
        ver(:)%nover=0

        vno=0 ; tlink=4000 ; pbcshift=gsize*blen                         ! Initializing the position

        row_loop : Do row=1,gsize,1
        cz=depth*sin(2*pi*row*period/gsize)
        cx=(row)*blen
        col_loop: Do col=1,gsize,1
        vno=vno+1 ; cy=(col)*blen
        ver(vno)%vcoord=reshape((/cx,cy,membrane_init_z+cz/),shape(pos))                 ! Position of vertex vno
        Enddo col_loop
        Enddo row_loop
        nver=vno ; pbnum=vno                                                                                                        !Periodic images will have index greater than the number of vertex

        find_neigh_pt: Do  v=1,vno,1
          ver(v)%nonei=0 ; ver(v)%vneipt=0 ; ver(v)%boundary=0 ; ver(v)%PBCver=0

          c1:If(v.LE.gsize)Then                                          ! upper row
           upper_row:If(v.EQ.1)Then                                      ! Left upper corner
              ver(v)%boundary=1
              ver(v)%nonei=6
              ver(v)%PBCver(1:6)=(/0,0,0,1,1,1/)
              ver(v)%vneipt(1:6)=(/gsize+v,gsize+v+1,v+1,gsize*(gsize-1)+1,gsize**2,gsize/)

                  tvno=gsize*(gsize-1)+1
                  pbx=ver(tvno)%vcoord(1,1)-pbcshift
                  pby=ver(tvno)%vcoord(2,1)                               ! New coordinates of the boundary image
                  pbz=ver(tvno)%vcoord(3,1)
                  If(ver(tvno)%pbimno.Eq.0)Then
                  flagallot=0
                  Else
                  flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(4)=chkimg
                   Endif
                  Endif
                  If(flagallot.Eq. 0) Then                               ! The periodic images are stored as seperate variables
                  call allotposition(pbx,pby,pbz,tvno,v,4)
                  Endif

               tvno=gsize**2
                 pbx=ver(tvno)%vcoord(1,1)-pbcshift
                 pby=ver(tvno)%vcoord(2,1)-pbcshift                      ! New coordinates of the boundary image
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                 flagallot=1
                  call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(5)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then                                ! The periodic images are stored as seperate variables
                 call allotposition(pbx,pby,pbz,tvno,v,5)
                 Endif


               tvno=gsize
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)-pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(6)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,6)
                 Endif



             Else If(v.EQ.gsize)Then                                     ! Right upper corner
                ver(v)%boundary=1
                ver(v)%nonei=6
                ver(v)%PBCver(1:6)=(/0,1,1,1,1,0/)
                ver(v)%vneipt(1:6)=(/v+gsize,v+1,1,gsize**2,gsize**2-1,v-1/)

               tvno=v+1
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)+pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                 flagallot=1
                 call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(2)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,2)
                 Endif

               tvno=1
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)+pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(3)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,3)
                 Endif

               tvno=gsize**2
                 pbx=ver(tvno)%vcoord(1,1)-pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                  flagallot=1
                  call ifalreadyalloted(pbx,pby,pbz,tvno)
                  If(chkflag.Eq.0)Then
                  flagallot=0
                  Else
                  ver(v)%vneipt(4)=chkimg
                  Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,4)
                 Endif

               tvno=gsize**2-1
                 pbx=ver(tvno)%vcoord(1,1)-pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(5)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,5)
                 Endif

             Else
               ver(v)%boundary=1
               ver(v)%nonei=6
               ver(v)%PBCver(1:6)=(/0,0,0,1,1,0/)
               ver(v)%vneipt(1:6)=(/v+gsize,v+gsize+1,v+1,gsize*(gsize-1)+v, gsize*(gsize-1)+v-1,v-1/)     ! The rest of the points on the top row

               tvno=gsize*(gsize-1)+v

                 pbx=ver(tvno)%vcoord(1,1)-pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0

          Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(4)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,4)
                 Endif

               tvno=gsize*(gsize-1)+v-1
                 pbx=ver(tvno)%vcoord(1,1)-pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(5)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,5)
                 Endif

               Endif upper_row



            Else If(v.GT.(vno-gsize))Then                                ! Lower row

               lower_row:If(v.EQ.(vno-gsize+1))Then                      ! Vertex at left bottom corner
               ver(v)%boundary=1                                         ! Vertices with spl=1 are those on the boundary
               ver(v)%nonei=6
               ver(v)%PBCver(1:6)=(/1,1,0,0,1,1/)
               ver(v)%vneipt(1:6)=(/1,2,v+1,v-gsize,v-1,gsize**2/)

               tvno=1
                 pbx=ver(tvno)%vcoord(1,1)+pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=depth*sin(2*pi*(gsize+1)*period/gsize)+membrane_init_z
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(1)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,1)
                 Endif

               tvno=2
                 pbx=ver(tvno)%vcoord(1,1)+pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=depth*sin(2*pi*(gsize+1)*period/gsize)+membrane_init_z
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(2)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,2)
                 Endif

               tvno=v-1
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)-pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(5)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,5)
                 Endif

               tvno=gsize**2
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)-pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(6)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,6)
                 Endif


             Else If(v.EQ.vno)Then                                       ! Right bottom corner
               ver(v)%boundary=1
               ver(v)%nonei=6
               ver(v)%PBCver(1:6)=(/1,1,1,0,0,0/)
               ver(v)%vneipt(1:6)=(/gsize,1,v-gsize+1,v-gsize,v-gsize-1,v-1/)

               tvno=gsize
                 pbx=ver(tvno)%vcoord(1,1)+pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=depth*sin(2*pi*(gsize+1)*period/gsize)+membrane_init_z
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(1)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,1)
                 Endif

               tvno=1
                 pbx=ver(tvno)%vcoord(1,1)+pbcshift
                 pby=ver(tvno)%vcoord(2,1)+pbcshift
                 pbz=depth*sin(2*pi*(gsize+1)*period/gsize)+membrane_init_z
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(2)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,2)
                 Endif

               tvno=v-gsize+1
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)+pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(3)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,3)
                 Endif

             Else
               ver(v)%boundary=1
               ver(v)%nonei=6
               ver(v)%PBCver(1:6)=(/1,1,0,0,0,0/)
               ver(v)%vneipt(1:6)=(/v-gsize*(gsize-1),v+1-gsize*(gsize-1),v+1,v-gsize,v-gsize-1,v-1/)

               tvno=v-gsize*(gsize-1)
                 pbx=ver(tvno)%vcoord(1,1)+pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=depth*sin(2*pi*(gsize+1)*period/gsize)+membrane_init_z
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(1)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,1)
                 Endif

               tvno=v-gsize*(gsize-1)+1
                 pbx=ver(tvno)%vcoord(1,1)+pbcshift
                 pby=ver(tvno)%vcoord(2,1)
                 pbz=depth*sin(2*pi*(gsize+1)*period/gsize)+membrane_init_z
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(2)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,2)
                 Endif

               Endif lower_row

           Else If(mod(v,gsize).EQ.1)Then                                ! left column
               ver(v)%boundary=1
               ver(v)%nonei=6
               ver(v)%PBCver(1:6)=(/0,0,0,0,1,1/)
               ver(v)%vneipt(1:6)=(/v+gsize,v+gsize+1,v+1,v-gsize,v-1,v-1+gsize/)

               tvno=v-1
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)-pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(5)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,5)
                 Endif

               tvno=v-1+gsize
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)-pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                  flagallot=1
                  call ifalreadyalloted(pbx,pby,pbz,tvno)
                  If(chkflag.Eq.0)Then
                  flagallot=0
                  Else
                  ver(v)%vneipt(6)=chkimg
                  Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,6)
                 Endif

              Else If(mod(v,gsize).EQ.0)Then                             ! Right column
               ver(v)%boundary=1
               ver(v)%nonei=6
               ver(v)%PBCver(1:6)=(/0,1,1,0,0,0/)
               ver(v)%vneipt(1:6)=(/v+gsize,v+1,v-gsize+1,v-gsize,v-gsize-1,v-1/)

               tvno=v+1
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)+pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(2)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,2)
                 Endif

               tvno=v+1-gsize
                 pbx=ver(tvno)%vcoord(1,1)
                 pby=ver(tvno)%vcoord(2,1)+pbcshift
                 pbz=ver(tvno)%vcoord(3,1)
                 If(ver(tvno)%pbimno.Eq.0)Then
                 flagallot=0
                 Else
                   flagallot=1
                   call ifalreadyalloted(pbx,pby,pbz,tvno)
                   If(chkflag.Eq.0)Then
                   flagallot=0
                   Else
                   ver(v)%vneipt(3)=chkimg
                   Endif
                 Endif
                 If(flagallot.Eq. 0) Then
                 call allotposition(pbx,pby,pbz,tvno,v,3)
                 Endif

               Else
                 ver(v)%nonei=6
                 ver(v)%PBCver(1:6)=(/0,0,0,0,0,0/)
                 ver(v)%vneipt(1:6)=(/v+gsize,v+gsize+1,v+1,v-gsize,v-gsize-1,v-1/)
          Endif c1
        Enddo find_neigh_pt



! So far we have made a square lattice and had defined its neighbours inclusive of a PBC in x and y directions

! The above part may look longer but it is a simple by hand implementation of the PBC. PBC for connected
! triangulated surface with bond length constraints are a bit complicated and to take care of this, I in
! this code, maintain a seperate PBC coordinate corresponding to each vertex and update it too when the
! vertex under question is disturbed

! ver(v)%PBCver --> a value of 1 here indicates that this position is occupied by PB vertex. So use
! the coordinates stored in ver()%vcoord

! The links from each vertex is alloted and the boundary condition works
! fine as far as the neighbour list is concerned.

! The periodic vertices are marked by another type structre called
! pbvertex

! Start  from alloting the links and making the triangles

! How to allot the links. Should we allot the links that includes the
! periodic vertices too ?  In princple the links between the boundaary
! nodes are to preserved in neighbourhood. So the details of the links
! that arising from PBC are not essential and so are those for such
! triangles

      lno=0
       Do v=1,pbnum,1
         Do v1=1,ver(v)%nonei,1
         prchk=0  ; linflag=0                                            ! To check if this link is already alloted ?
         nv1=ver(v)%vneipt(v1) ; i=-lno                                ! v1 th neighbour of v ; linno --> size of array for link list
           Do WHILE(i.LE.lno.And.prchk.EQ.0)                           ! check over all the link list
           If(i.NE.0)Then                                                ! There is no link named zero since every l has a -l
           If(lin(i)%sep(1).Eq.v .And. lin(i)%sep(2).Eq.nv1 .And.linflag .Eq. lin(i)%pbflag) Then
              prchk=1 ; Endif
           If(lin(i)%sep(2).EQ.v .And. lin(i)%sep(1).EQ.nv1 .And.linflag .Eq. lin(i)%pbflag) Then
              prchk=1 ; Endif
           Endif
           i=i+1
           Enddo

           If (prchk.EQ.0)Then                                           ! If no existing link has been found then
           lno=lno+1
           lin(lno)%sep=(/v,nv1/)
           lin(-lno)%sep=(/nv1,v/)
           If((ver(v)%boundary.EQ.1).And.(ver(nv1)%boundary.EQ.1))Then   ! Mark all the boundary links
           lin(lno)%boundary=1 ; lin(-lno)%boundary=1
           Endif  ; Endif
       Enddo ; Enddo
       tlink=lno                                                         ! Total number of links in the system


! Link allocation was easy making the triangles seems to be difficult.
! Another idea that can be thought of is to remove the periodic vertex
! structure and make it a part of the vertex structure by giving it a
! seperate flag.

! vertex number > nver is a marker for the periodic boundaries

        trn=0
        tri(:)%pbflag=0                                                  ! Set all peridoic boundary flag to zero
        Do v=1,nver,1
             nntr=0
             Do v1=1,ver(v)%nonei,1                                      ! Run over the neighbourhood of v
             prchk=0 ;nechk=0                                            ! prchk --> check if already alloted
             nv1=ver(v)%vneipt(v1) ; nv2=ver(v)%vneipt(v1+1)             ! the i th and i+1 th neighbour, be careful Periodic images can occur
             If(v1.EQ.ver(v)%nonei) nv2=ver(v)%vneipt(1)                 ! For all submarginal vertices there exists a circular boundary
                                                                         ! check with all neig triangles of nv1 for previous allotment
             tm1=1
             Do WHILE(tm1.LE.ver(nv1)%nonei.And.prchk.EQ.0)              ! For each neigh triangle around vertex 'v', ^^^^^^^^^^^
             tm2=ver(nv1)%vneitr(tm1)                                    ! Pick up each triangle  as tm2
	     If (tm2 .Gt. 0) Then
             tmv=tri(tm2)%vert                                           ! vertices of tm2
             If(v.EQ.tmv(1).OR.v.EQ.tmv(2).OR.v.EQ.tmv(3))Then           ! if vertcies of tm2 == (v,nv1,nv2) or (nv2,v,nv1)
             If(nv1.EQ.tmv(1).OR.nv1.EQ.tmv(2).OR.nv1.EQ.tmv(3))Then     ! or (nv1,nv2,v)
             If(nv2.EQ.tmv(1).OR.nv2.EQ.tmv(2).OR.nv2.EQ.tmv(3))Then
             nntr=nntr+1                                                 ! Number of nearest triangles list is increased by 1
             ver(v)%vneitr(nntr)=tm2 ; prchk=1 ;nechk=1                  ! If present then vert --> trian link for vertex v is tm2
	     Endif ; Endif ; End If; End If ;tm1=tm1+1;Enddo                       ! Exit the loop once presence is already detected .

                                                                         ! check with all neig triangles of nv2 for previous allotment
             tm1=1                                                       ! spl vertices have no circular boundaries
             Do WHILE(tm1.LE.ver(nv2)%nonei.And.prchk.EQ.0)              ! For each neigh triangle around vertex 'v', ^^^^^^^^^^^
             tm2=ver(nv2)%vneitr(tm1)                                    ! Pick up each triangle  as tm2
	     If (tm2 .Gt. 0) Then
             tmv=tri(tm2)%vert                                           ! vertices of tm2
             If(v.EQ.tmv(1).OR.v.EQ.tmv(2).OR.v.EQ.tmv(3))Then           ! if vertcies of tm2 == (v,nv1,nv2) or (nv2,v,nv1)
             If(nv1.EQ.tmv(1).OR.nv1.EQ.tmv(2).OR.nv1.EQ.tmv(3))Then     ! or (nv1,nv2,v)
             If(nv2.EQ.tmv(1).OR.nv2.EQ.tmv(2).OR.nv2.EQ.tmv(3))Then
             nntr=nntr+1                                                 ! Number of nearest triangles list is increased by 1
             ver(v)%vneitr(nntr)=tm2 ; prchk=1 ;nechk=1                  ! If present then vert --> trian link for vertex v is tm2
	     Endif; Endif ; End If; End If
             tm1=tm1+1;Enddo

             If(nechk.EQ.0)Then                                          ! If presence not detected a new triangle is formed  &&&&&&&&&&&
              nntr=nntr+1 ; trn=trn+1                                    ! Both number of nearest triangle and total no is increased
              ver(v)%vneitr(nntr)=trn
              tri(trn)%vert=reshape((/v,nv1,nv2/),shape(tri(trn)%vert))  ! Vertices of the new traingle

              If(v.Gt.nver .Or. nv1.Gt.nver .Or. nv2.Gt.nver)Then        ! $ seperately for excluding in rendering later
              tri(trn)%pbflag=1
              Endif


              tm1=-tlink;prchk=0                                         ! tm1 is a loop variable to see if a link already exists for the vertices
              Do WHILE(tm1.LE.tlink .And. prchk.EQ.0)                    ! Check for the existence of the links
              If((tm1).NE.0)Then                                         ! Link number zero is not present in the system
              If(lin(tm1)%sep(1).EQ.v.And.lin(tm1)%sep(2).EQ.nv1)Then    ! prchk --> variable to break out from the loop the link exists
              tri(trn)%li(1)=tm1 ; prchk=1                               ! nechk --> flag to allot a new link if a link does not exists
              lin(tm1)%tr=trn
              Endif ; End If  ;tm1=tm1+1
              Enddo

              tm1=-tlink ; prchk=0
              Do WHILE(tm1.LE.tlink .And. prchk.EQ.0)                    ! Check for the existence of the links now for nv1,nv2 as before
              If((tm1).NE.0)Then
              If(lin(tm1)%sep(1).EQ.nv1.And.lin(tm1)%sep(2).EQ.nv2)Then
              tri(trn)%li(2)=tm1 ; prchk=1
              lin(tm1)%tr=trn
              Endif ; End If  ;tm1=tm1+1
              Enddo

              tm1=-tlink ; prchk=0
              Do WHILE(tm1.LE.tlink .And. prchk.EQ.0)                    ! Check for the existence of the links  between nv2 --> v
              If((tm1).NE.0)Then
              If(lin(tm1)%sep(1).EQ.nv2.And.lin(tm1)%sep(2).EQ.v)Then
              tri(trn)%li(3)=tm1 ; prchk=1 ;nechk=1
              lin(tm1)%tr=trn
              Endif ; End If  ;tm1=tm1+1
              Enddo

             Endif                                                       ! &&&&&&&&&&&
            Enddo                                                        ! ^^^^^^^^^^^
        Enddo
        ntr=trn


           Do j=1,ntr,1                                                     ! Allot the neighbouring triangles to the periodic image boundaries
           v=tri(j)%vert(1)
           v1=tri(j)%vert(2)
           v2=tri(j)%vert(3)

              If(v>nver)Then
              ver(v)%nonei=ver(v)%nonei+1
              ver(v)%vneitr(ver(v)%nonei)=j
              Endif

              If(v1>nver)Then
              ver(v1)%nonei=ver(v1)%nonei+1
              ver(v1)%vneitr(ver(v1)%nonei)=j
              Endif

              If(v2>nver)Then
              ver(v2)%nonei=ver(v2)%nonei+1
              ver(v2)%vneitr(ver(v2)%nonei)=j
              Endif
           Enddo

           ver(nver+1:pbnum)%nonei=ver(nver+1:pbnum)%nonei+1             ! Since the edges do not have a circular boundary no of vertex=no of triang+1


           Do j=1,ntr,1                                                  ! Allot the neighbouring vertices to the boundaries (careful they do not $
           v=tri(j)%vert(1)                                              !$ follow any sequence as done by the inner nodes
           v1=tri(j)%vert(2)
           v2=tri(j)%vert(3)

            If(v>nver)Then
              vflag1=0 ; vflag2=0
              Do i=1,ver(v)%nonei,1
              If(ver(v)%vneipt(i).Eq.v1) vflag1=1
              If(ver(v)%vneipt(i).Eq.v2) vflag2=1
              Enddo

              If(vflag1.Eq.0) Then
              ver(v)%nover=ver(v)%nover+1
              ver(v)%vneipt(ver(v)%nover)=v1
              Endif

              If(vflag2.Eq.0) Then
              ver(v)%nover=ver(v)%nover+1
              ver(v)%vneipt(ver(v)%nover)=v2
              Endif
            Endif

            If(v1>nver)Then
              vflag1=0 ; vflag2=0
              Do i=1,ver(v1)%nonei,1
              If(ver(v1)%vneipt(i).Eq.v)  vflag1=1
              If(ver(v1)%vneipt(i).Eq.v2) vflag2=1
              Enddo

              If(vflag1.Eq.0) Then
              ver(v1)%nover=ver(v1)%nover+1
              ver(v1)%vneipt(ver(v1)%nover)=v
              Endif

              If(vflag2.Eq.0) Then
              ver(v1)%nover=ver(v1)%nover+1
              ver(v1)%vneipt(ver(v1)%nover)=v2
              Endif
            Endif


           If(v2>nver)Then
              vflag1=0 ; vflag2=0
              Do i=1,ver(v2)%nonei,1
              If(ver(v2)%vneipt(i).Eq.v)  vflag1=1
              If(ver(v2)%vneipt(i).Eq.v1) vflag2=1
              Enddo

              If(vflag1.Eq.0) Then
              ver(v2)%nover=ver(v2)%nover+1
              ver(v2)%vneipt(ver(v2)%nover)=v
              Endif

              If(vflag2.Eq.0) Then
              ver(v2)%nover=ver(v2)%nover+1
              ver(v2)%vneipt(ver(v2)%nover)=v1
              Endif
           Endif
       Enddo

        Do i=1,ntr,1                                                     ! Calculate the area of the triangles
        call areacalc(i)
        Enddo

        End Subroutine Make_sinusoidal_Membrane	
	
!---------------------------------------------------------------------------------------------------------------------------
!                 Subroutine to map surface quantifiers at a point to its periodic image
!---------------------------------------------------------------------------------------------------------------------------
! Call  this everytime a vertex with a periodic image is moved or
! involved in a successsful flip operation. When used for a code
! with the surface fields, include the update for the spin terms too.
      Subroutine mapimagetovertex(vimg)
      Use module_datastruct
      Implicit None
      Integer :: v,vimg
         v=ver(vimg)%imver                                                                                                          ! Pick the original vertex corresponding to this image
         ver(vimg)%cur1=ver(v)%cur1 ; ver(vimg)%cur2=ver(v)%cur2 
         ver(vimg)%mcur=ver(v)%mcur ; ver(vimg)%L2G=ver(v)%L2G
         ver(vimg)%t1=ver(v)%t1     ; ver(vimg)%t2=ver(v)%t2
         ver(vimg)%vnor=ver(v)%vnor ; ver(vimg)%totarea=ver(v)%totarea
         ver(vimg)%Antigen_flag=ver(v)%Antigen_flag
	   ver(vimg)%czero=ver(v)%czero; ver(vimg)%czero_flag=ver(v)%czero_flag
      End Subroutine mapimagetovertex

!---------------------------------------------------------------------------------------------------------------------------
!                 Subroutine to map surface quantifiers at a point to all its periodic image
!---------------------------------------------------------------------------------------------------------------------------
! Call  this everytime a vertex with a periodic image is moved or
! involved in a successsful flip operation. When used for a code
! with the surface fields, include the update for the spin terms too.

      Subroutine mapvertextoimages(vert)
      Use module_datastruct
      IMPLICIT NONE
      Integer :: vert,vimg,i
      	Do i=1,ver(vert)%pbimno,1
	   vimg = ver(vert)%pbmap(i)
           ver(vimg)%cur1=ver(vert)%cur1 ; ver(vimg)%cur2=ver(vert)%cur2 
           ver(vimg)%mcur=ver(vert)%mcur ; ver(vimg)%L2G=ver(vert)%L2G
           ver(vimg)%t1=ver(vert)%t1     ; ver(vimg)%t2=ver(vert)%t2
           ver(vimg)%vnor=ver(vert)%vnor ; ver(vimg)%totarea=ver(vert)%totarea
           ver(vimg)%Antigen_flag=ver(vert)%Antigen_flag
	   ver(vimg)%czero=ver(vert)%czero; ver(vimg)%czero_flag=ver(vert)%czero_flag
	Enddo
      End Subroutine mapvertextoimages

!---------------------------------------------------------------------------------------------------------------------------
!                 Subroutine to check if the periodic boundary of a vertice is already named 
!---------------------------------------------------------------------------------------------------------------------------
       Subroutine ifalreadyalloted(pbx,pby,pbz,vno)
       Use module_datastruct
       Implicit None
       Integer :: vno,nimg,i,imgVN
       Real(Kind=8) :: pbx,pby,pbz
       i=0
       nimg=ver(vno)%pbimno ; chkflag=0
       Do while(i.Le.nimg .And. chkflag.Eq.0)
       imgVN=ver(vno)%pbmap(i)                                                                                                      ! The vertex corresponding to image number i 
       If(ver(imgVN)%vcoord(1,1).Eq.pbx .And. ver(imgVN)%vcoord(2,1).Eq.pby .And. ver(imgVN)%vcoord(3,1).Eq.pbz )Then
       chkflag=1 ; chkimg=ver(vno)%pbmap(i)
       Endif
      i=i+1
       Enddo
       End Subroutine ifalreadyalloted
!-------------------------------------------------------------------------------------------------------------------------------
!                Subroutine to allocate the position of the periodic image
!-------------------------------------------------------------------------------------------------------------------------------
       Subroutine allotposition(pbx,pby,pbz,vno,v,p)
       Use module_datastruct
       Implicit None
       Integer :: vno,v,p
       Real(Kind=8) :: pbx,pby,pbz
         pbnum=pbnum+1                                                   ! mainly due to the reason being a connected map            
         ver(vno)%pbimno=ver(vno)%pbimno+1                               ! Each vertex can have one or multiple periodic images      
         ver(vno)%pbmap(ver(vno)%pbimno)=pbnum                           ! mapping from real space to an image in the periodic space 
         ver(pbnum)%vcoord(1,1)=pbx
         ver(pbnum)%vcoord(2,1)=pby                                      ! New coordinates of the boundary image 
         ver(pbnum)%vcoord(3,1)=pbz
         ver(v)%vneipt(p)=pbnum
         ver(pbnum)%imver=vno
         ver(pbnum)%boundary=1
       End Subroutine allotposition

!---------------------------------------------------------------------------------------------------------------------------
!                 Subroutine to calculate the area of the given face and update the total area linked to the vertices
!---------------------------------------------------------------------------------------------------------------------------

       Subroutine areacalc(tr)
       USE module_datastruct 
       Implicit None
       Integer::i,j,k,tr
       Real(KIND=8) ::ax,ay,az,area
       Real(KIND=8),DIMENSION(3,1)::r1,r2,r3,r21,r31
       ax=0 ; ay=0 ; az=0 ; area=0

       i=tri(tr)%vert(1) ; j=tri(tr)%vert(2) ; k=tri(tr)%vert(3)
       r1=ver(i)%vcoord ; r2=ver(j)%vcoord ; r3=ver(k)%vcoord                                                                       ! The coordinates of the 3 vertices      
       r21=r2-r1 ; r31=r3-r1                                                                                                        ! Relative position vectors for area calc

       ver(i)%totarea=ver(i)%totarea-tri(tr)%ar*0.3333333                                                                           ! Contrib of exisiting area to totalarea 
       ver(j)%totarea=ver(j)%totarea-tri(tr)%ar*0.3333333
       ver(k)%totarea=ver(k)%totarea-tri(tr)%ar*0.3333333

        ax=(r21(2,1)*r31(3,1)-r21(3,1)*r31(2,1))*0.5
        ay=(r21(3,1)*r31(1,1)-r21(1,1)*r31(3,1))*0.5
        az=(r21(1,1)*r31(2,1)-r21(2,1)*r31(1,1))*0.5
        area=SQRT(ax**2+ay**2+az**2)
        tri(tr)%ar=area                                       
        ax=ax/area ; ay=ay/area ;  az=az/area                                                                                       ! Normalized area vector 

       tri(tr)%fnor=RESHAPE((/ax,ay,az/),SHAPE(tri(tr)%fnor)) 
       ver(i)%totarea=ver(i)%totarea+tri(tr)%ar*0.3333333                                                                           ! Contrib of new area to the total area of vert
       ver(j)%totarea=ver(j)%totarea+tri(tr)%ar*0.3333333
       ver(k)%totarea=ver(k)%totarea+tri(tr)%ar*0.3333333

       tri(tr)%vol=(r1(1,1)*(r2(2,1)*r3(3,1)-r2(3,1)*r3(2,1))+r1(2,1)*(r2(3,1)*r3(1,1)-r2(1,1)*r3(3,1))+ &                          ! volume of the tetrahedron formed by triangle tr and origin (0,0,0)
                    r1(3,1)*(r2(1,1)*r3(2,1)-r2(2,1)*r3(1,1)))/6.0

       End Subroutine areacalc  

       End MODULE module_makesurface


