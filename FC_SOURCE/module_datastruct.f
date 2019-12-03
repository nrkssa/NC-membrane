!! ---------------- Important (Note by Ramakrishnan Natesan) -----------------------------------
!! C++ and Fortran share the same memory elements. If you add or remove a variable from any of the type structures below, 
!! the same procedure should be done, at exactly the same location, in the corresponding structure defined in 
!! CC_INCLUDE/_fortran_structures.h. If not C++ part cannot read/write to the variable you expect it to !!!

!==========================================================================================================================
!===========================================================================================================================
!                   $moddata                   MODULE CONTAINING THE DATASTRUCTURE
!===========================================================================================================================
!==========================================================================================================================

      MODULE module_datastruct
        IMPLICIT NONE
        Include '_string_formatter.h'
        Integer,Parameter :: verno = 91202, trino = 181200, linno = 271199, maxnc = 100, maxshadowvert = 10000        ! g**2+(4*g+2),2*g**2+4*g,3*g**2+4*g-1          
        Integer,Parameter :: maxantigen_ver = 10
        Integer,Parameter :: maxantigen_tri = 3*maxantigen_ver
        Integer :: processor_number
        Integer:: itime, ftime, noli, mcs, seed, iloop, curr_step                                                      ! mcs-->mcs step num |iloop --> inner loop num  
        Integer:: nver, ntr, tlink, num_antigens, nlinkcells, num_nanocarrier, num_bias_moves, num_antibodies, minshadow_size
        Integer:: gsize, pbnum 
        Integer:: datainterval, mcsmovecounter
        Integer :: vermov_attempt, vermov_accepted, vermov_step_adjust
        Real(Kind=8) :: vermov_step_size, reac_energy_change
        Real(KIND=8) :: kappa_unscaled, kappa, czero, beta, blen, spcu, sppar, fangle, pressure      
        Real(KIND=8) :: scur, pi, depth, period, zero, blen_scale_factor                
        Real(Kind=8) :: periodic_box_length, periodic_box_height, membrane_init_z                ! Dimensions of the periodic box used to contain the Nanocarrier    
        Real(KIND=8) :: blUcut, blLcut, drzero(3,1)                                              ! Magnitude of displacement and maximum displacement in vertex move 
        Real(Kind=8),Parameter :: infinity = 100.00                                                ! just a large number 
        Logical :: memb_restart_flag, antigen_restart_flag, fixed_frame_flag, debug_mode
        PARAMETER (scur = 0.0,zero = 0.0,pi = acos(-1.00000000000),fangle = -0.5) 
        Character*100 memb_geom
        
         Type vertex
          Real(KIND=8),DIMENSION(3,1) :: vnor,vcoord,t1,t2                                        ! Orientation of the nematic in local and global sys 
          Integer,DIMENSION(10) :: vneipt,vneitr,PBCver                                           ! List of all neigh points,triangles,links and PBC   
          Real(KIND=8),DIMENSION(3,3) :: L2G,HHM
          Integer :: nonei,boundary,pbimno,pbmap(3),imver,nover,boxvert                           ! Number of Neighbours a vertex has | number of the image
          Real(KIND=8) :: mcur,cur1,cur2,czero,totarea                                            ! pbimno --> number of vertex in the periodic description
          Integer :: antigen_flag,cellno,neigh                                                    ! Antigen connected to a vertex
          Integer :: antigen_list(10),nantigen,shadownc(maxnc),czero_flag
         End Type vertex                                                                          ! A single vertex at the edge can have 2 Periodic image   

         Type Antigen
         Integer :: vertex                                                             ! membrane vertex an antigen is attached to
         Real(Kind=8) :: base_coord(3,1), tip_coord(3,1), theta, phi                   ! base and tip position of the antigen   
         Real(Kind=8) :: disp1, disp2, length, radius, kflex                           ! Variables to determine the current position of the antigen in the triangle
         Integer :: diffus_tri, ant_type                                                ! The triangle it freely diffuses on 
         End Type Antigen

        Type triangle                                                                  ! Triangle details                                 
         Real(KIND=8) :: ar                                                            ! Area of each face and their corresponding volume 
         Integer :: pbflag,boxtriangle                                                 ! Does this has a vertex with a periodic boundary ?
         Integer :: nantigen                                         ! number of antigens living on this triangle (maximum 3, one corresponding to each vertex) 
         Integer,DIMENSION(3) :: li,vert                             ! Links and vertices that make up the triangle, antigens on a given triangle
         Real(KIND=8),DIMENSION(3,1) :: fnor
         Real(Kind=8) :: HHM(3,3)          
         Integer :: antigen_list(30)                                 ! A triangle can accomodate 3*maxantigen (thrice the max. number of antigens on a vertex)
         Real(Kind=8) :: vol
        End Type triangle 
         
        Type link
         Integer :: tr,boundary,pbflag                                                              ! Triangle associated   
         Integer,DIMENSION(2) :: sep                                                                ! Begin and end vertex  
        End Type link

        Type geom_boundary
         Real(Kind=8) :: coord(3,1)
         Integer :: vert
        End Type geom_boundary

        Type Nanocarrier
         Integer :: nshadow_vert, shadow_vertices(maxshadowvert)       ! maxshadowvert shoule be increased in both here and in C++ if particle size becomes large
         Real(Kind=8) :: kbias, biasref, lambda, shadow_cutoff, shadow_cutoff_sq
         Real(Kind=8) :: radius, soft_radius, multivalency
         Real(Kind=8) :: coord(3,1), meanr(3,1), meancurv, meandist
         Character :: bias_mode
         Integer :: pbcx(maxshadowvert), pbcy(maxshadowvert)
        End Type Nanocarrier

        Type membrane_prop
         Type(vertex),Allocatable,Dimension(:) :: vertex
         Type(triangle),Allocatable,Dimension(:) :: triangle
         Type(link),Allocatable,Dimension(:)  :: link
         Type(antigen),Allocatable,Dimension(:) :: antig
         Real(KIND=8) :: enel,ento,enss,ifaclen                                                                     ! The different energy terms involved               
         Real(KIND=8) :: vol,area,rg,blen,projarea                                                                  ! Volume, area, Radius of Gyration,average blength  
         Type(Nanocarrier) :: nc_f(maxnc)
        End Type membrane_prop
        
        Type spec_antigens
         Integer :: ntotal, ntype
         Integer, Allocatable, Dimension(:) :: num_ntype
         Real(Kind=8), Allocatable, Dimension(:) :: radius, length, kflex
         Character(100) :: pattern
        End Type spec_antigens
        
        !@ Self avoidance distances
        Type SA_distances
         Real(Kind=8), Allocatable :: antant(:,:), antves(:,:), vesmem(:)
        End Type SA_distances

        Type(vertex) :: ver(verno)                                             ! Setting the size for type vertex      
        Type(triangle) :: tri(trino)                                           ! Setting the size for type triangle    
        Type(link) :: lin(-linno:linno)                                        ! Setting the size for type link        
        Type(Antigen) :: antig(verno)
        Type(membrane_prop) :: mp
        Type(geom_boundary) :: geomb(8)                                        ! eight vertices that defines the boundary of the membrane + box under it   
        Type(Nanocarrier) :: nc_f(maxnc) 
        Type(spec_antigens) :: antspec
        Type(SA_distances) :: sad
        
!-----------------------------------------------------------------------------------------------------------------------------------

        Contains

!==================================================================================================================================	
!                                   Interactive shell to print the various values during runtime
!==================================================================================================================================	

	Subroutine Interactive_Shell()
	Use module_string_utility
	Implicit None
	Integer :: objno
	Character(200) :: shell_cmd,object,objvalue

		shell_cmd='You are now in the interactive shell mode'
		Print*,str_green,shell_cmd,str_black
		Print*,'Enter your queries as follows: "shell_cmd","object","objvalue","obj_no"'
		Print*,'----------------------------------------------------------------------------'
		Print*,str_blue,'Choose "shell_cmd " from ("print", "exit")',str_black
		Print*,str_blue,'Choose "object" from ("vertex", "triangle", "link", "antigen")',str_black
		Print*,str_red ,'CHOOSE "objvalue" for object', str_black,str_green,'"vertex" from :',str_black
    	Print*,str_blue,'("coord" "nantigen" "neigh_tri" "neigh_vert" "ant_list")',str_black
		Print*,str_red ,'CHOOSE "obj_value" for object', str_black,str_green,'"triangle" from :',str_black
		Print*,str_blue,'("vert" "nantigen" "ant_list")',str_black
		Print*,str_red ,'CHOOSE "obj_value" for object', str_black,str_green,'"link" from:',str_black
		Print*,str_blue,'("vert" "tri")',str_black
		Print*,str_red ,'Choose "obj_value" for object', str_black,str_green,'"antigen": from',str_black
		Print*,str_blue,'("base" "tip" "vert"	"tri")',str_black
		Print*,str_green,'Choose "obj_no" in the following range:',str_black
		Print*,str_blue,'"vertex":  1-->',pbnum
		Print*,str_blue,'"triangle": 1-->',ntr
		Print*,str_blue,'"link":',-tlink,'--> ',tlink,'(zero is not allowed)'
		Print*,str_blue,'"antigen": 1 --> ',num_antigens,str_black
		Print*,'----------------------------------------------------------------------------'

	Do while (Trim(adjustl(StrLowCase(shell_cmd))) .Ne. 'exit')
		Read(*,*)shell_cmd,object,objvalue,objno
		If(Trim(Adjustl(StrLowCase(object))) .Eq. 'vertex') Then
			If(Trim(Adjustl(objvalue)).Eq. 'coord') Then
				Print*,str_red,'Printing Coordinate of vertex ',objno,' :',str_black
				Print*,ver(objno)%vcoord
				Print*,str_equal

			Else If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'nantigen') Then
				Print*,str_red,'Number of antigens on vertex ',objno,' : ',ver(objno)%nantigen,str_black
				Print*,str_equal

			Else If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'ant_list') Then
				Print*,str_red,'Number of antigens on vertex ',objno,' : ',ver(objno)%nantigen,str_black
				Print*,str_red,'Antigen List ',ver(objno)%antigen_list(1:ver(objno)%nantigen),str_black
				Print*,str_equal

			Else If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'neigh_vert') Then
				Print*,str_red,'Vertex ',objno,' has ',ver(objno)%nonei,' neighbouring vertices',str_black
				Print*,str_red,'Nei Vert List ',ver(objno)%vneipt(1:ver(objno)%nonei),str_black
				Print*,str_equal

			Else If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'neigh_tri') Then
				Print*,str_red,'Vertex ',objno,' has ',ver(objno)%nonei,' neighbouring triangles',str_black
				Print*,str_red,'Nei Tri List ',ver(objno)%vneitr(1:ver(objno)%nonei),str_black
				Print*,str_equal

			Else If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'images') Then
				Print*,str_red,'Vertex ',objno,' has ',ver(objno)%pbimno,'images',str_black
				Print*,str_red,'Image List ',ver(objno)%pbmap(1:ver(objno)%pbimno),str_black
				Print*,str_equal
			Endif
				
		Else If(Trim(Adjustl(StrLowCase(object))) .Eq. 'triangle') Then
			If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'vert') Then
				Print*,str_red,'Coordinate forming triangle ',objno,' :(',tri(objno)%vert,')',str_black
				Print*,str_equal

			Else If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'nantigen') Then
				Print*,str_red,'Triangle ',objno,' has ',tri(objno)%nantigen,' antigens',str_black
				Print*,str_equal

			Else If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'ant_list') Then
				Print*,str_red,'Triangle ',objno,' has ',tri(objno)%nantigen,' antigens',str_black
				Print*,str_red,tri(objno)%antigen_list(1:tri(objno)%nantigen),str_black
				Print*,str_equal

			Endif

		Else If(Trim(Adjustl(StrLowCase(object))) .Eq. 'link') Then
			If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'vert') Then
				Print*,str_red,'Link ',objno,' connect vertex',lin(objno)%sep(1),' to ',lin(objno)%sep(2),str_black
				Print*
			Else If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'tri') Then
				Print*,str_red,'Link ',objno,' contains triangle',lin(objno)%tr,str_black
				Print*
			Endif


		Else If(Trim(Adjustl(StrLowCase(object))) .Eq. 'antigen') Then
			If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'base') Then
				Print*,str_red,'Base coord of antigen ',objno,':',str_black
				Print*,'(',antig(objno)%base_coord,')'
				Print*,str_equal

			Else If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'tip') Then
				Print*,str_red,'Tip coord of antigen ',objno,':',str_black
				Print*,'(',antig(objno)%tip_coord,')'
				Print*,str_equal

			Else If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'vertex') Then
				Print*,str_red,'Antigen ',objno,'is on vertex ',antig(objno)%vertex,str_black
				Print*,str_equal

			Else If(Trim(Adjustl(StrLowCase(objvalue))).Eq. 'triangle') Then
				Print*,str_red,'Antigen ',objno,' diffuses on triangle ',antig(objno)%diffus_tri,str_black
				Print*,str_equal

			Endif
		Endif
		Print*,'Enter the next query or enter "exit"' 
		Print*,'>>>'
	Enddo
	Return
	End Subroutine Interactive_Shell

!----------------------------------------------------------------------------------------------------------------------------
!                                SUBROUTINE TO COMPUTE THE TOTAL ENERGY COST DUE TO MEMBRANE BENDING
!----------------------------------------------------------------------------------------------------------------------------
	Function total_membrane_energy() Result(bend_energy)
	Implicit None
		Integer :: i
		Real(Kind=8) :: bend_energy
		bend_energy=0.0
		Do i=1,nver,1
			bend_energy = bend_energy + ((ver(i)%mcur-ver(i)%czero)**2)*ver(i)%totarea
		Enddo
		bend_energy= bend_energy*kappa*4.1E-21
		Return
	End Function total_membrane_energy

!----------------------------------------------------------------------------------------------------------------------------
!                                SUBROUTINE TO ITERATIVELY SET THE VERTEX MOVE STEP SIZE       
!----------------------------------------------------------------------------------------------------------------------------
    Subroutine Reset_vertexmove_Stepsize()
       Implicit None 

       If(Nint(Real(vermov_accepted)/vermov_attempt).Eq.1) Then                                            	                            ! If more that 50% of moves are accepted then increase step size
	       vermov_step_size=1.05*vermov_step_size
       Else
	       vermov_step_size=vermov_step_size/1.05                                                        	                    ! Else decrease stepsize
       Endif
       vermov_attempt=0 ; vermov_accepted=0
    End Subroutine Reset_vertexmove_Stepsize

!----------------------------------------------------------------------------------------------------------------------------
!                                UPDATE MEAN Z-POSITION OF THE MEMBRANE AFTER A  VERTEX MOVE   
!----------------------------------------------------------------------------------------------------------------------------
    Function Update_mean_membrane_z(oldmeanz,zold,znew,N) Result(meanz)
       Implicit None 
       Real(Kind=8) :: meanz,oldmeanz,zold,znew
       Integer :: N
       meanz=oldmeanz+((znew-zold)/N)                                                                                               ! The new meanz is (N<zold>-zold+znew)/N
       Return
    End Function Update_mean_membrane_z

!---------------------------------------------------------------------------------------------------------------------------       
!                                Function to get the location of a vertex within the shadow vertices of an NC
!---------------------------------------------------------------------------------------------------------------------------       
	Function get_shadow_index(ncnum,vert) Result(sindex)
	Implicit None
	Integer :: ncnum,vert,sindex
	sindex=1 
	Do while(sindex .Le. nc_f(ncnum)%nshadow_vert)
  	 If (nc_f(ncnum)%shadow_vertices(sindex) .Eq. vert)  Then
	  exit
	 Else
	  sindex = sindex+1
	 Endif
	Enddo
	Return
	End Function get_shadow_index
!----------------------------------------------------------------------------------------------------------------------------
!                                COMPUTE THE SHADOW OF A NANOCARRIER and its meanz and meancurv values      
!----------------------------------------------------------------------------------------------------------------------------
	Subroutine Compute_nc_shadow_vertices(ncnum) 
	Implicit None
    Include '_mod_getcoordinates_interface.h'
	Real(Kind=8) :: dx,dy,dz
	Integer :: ncnum,i,pbcx,pbcy

	 nc_f(ncnum)%shadow_vertices(1:pbnum) = 0
	 nc_f(ncnum)%nshadow_vert = 0
	 
     Do i=1,nver,1
	 	If(ver(i)%shadownc(ncnum) .Eq. ncnum) ver(i)%shadownc(ncnum) = 0
         dx = nc_f(ncnum)%coord(1,1) - ver(i)%vcoord(1,1) ; pbcx = Nint(dx/periodic_box_length)
         dx = dx - pbcx*periodic_box_length
         If (abs(dx) .Le. nc_f(ncnum)%shadow_cutoff) Then
         dy = nc_f(ncnum)%coord(2,1) - ver(i)%vcoord(2,1) ; pbcy = Nint(dy/periodic_box_length)
		 dy = dy - pbcy*periodic_box_length
          If (abs(dy) .Le. nc_f(ncnum)%shadow_cutoff) Then
            dz = nc_f(ncnum)%coord(3,1) - ver(i)%vcoord(3,1) 
	   	    If ((dx**2 + dy**2 + dz**2) .Le. nc_f(ncnum)%shadow_cutoff_sq) Then
            nc_f(ncnum)%nshadow_vert = nc_f(ncnum)%nshadow_vert +1 
            If (nc_f(ncnum)%nshadow_vert .Gt. maxshadowvert) Call Throw_error_message(1)
            nc_f(ncnum)%shadow_vertices(nc_f(ncnum)%nshadow_vert) = i
            nc_f(ncnum)%pbcx(nc_f(ncnum)%nshadow_vert) = pbcx                    ! Store the shifting factor for later use
            nc_f(ncnum)%pbcy(nc_f(ncnum)%nshadow_vert) = pbcy                    
            ver(i)%shadownc(ncnum) = ncnum
          Endif
         Endif
         Endif 
        Enddo

	Call ncshadow_r(ncnum)                                                               ! Compute the mean z and mean curvature in the shadow
	Call ncshadow_mcurv(ncnum)             
	End subroutine Compute_nc_shadow_vertices
!----------------------------------------------------------------------------------------------------------------------------
!                                COMPUTE THE SHADOW OF A NANOCARRIER and its meanz and meancurv values      
!----------------------------------------------------------------------------------------------------------------------------
       Subroutine Compute_nc_shadow_vertices_local(ncnum,vert) 
        Implicit None
        Include '_mod_getcoordinates_interface.h'
        Real(Kind=8) :: dx,dy,dz,length
        Integer :: ncnum,pbcx,pbcy,sindex,vert
        Logical :: inshadow

        inshadow = .False.
        If (ver(vert)%shadownc(ncnum) .Eq. ncnum) inshadow = .True.
         dx = nc_f(ncnum)%coord(1,1) - ver(vert)%vcoord(1,1) ; pbcx = Nint(dx/periodic_box_length)
         dy = nc_f(ncnum)%coord(2,1) - ver(vert)%vcoord(2,1) ; pbcy = Nint(dy/periodic_box_length)
         dx = dx - pbcx*periodic_box_length ; dy = dy - pbcy*periodic_box_length 
         dz = nc_f(ncnum)%coord(3,1) - ver(vert)%vcoord(3,1) 
         length = dx**2 + dy**2 + dz**2

	!@ if already present and moved out of shadow
	If ((length .Gt. nc_f(ncnum)%shadow_cutoff_sq) .And. inshadow) Then
		ver(vert)%shadownc(ncnum) = 0
		sindex = get_shadow_index(ncnum,vert)                                                                ! get location of vert in the shadow of ncnum
		nc_f(ncnum)%shadow_vertices(sindex) = nc_f(ncnum)%shadow_vertices(nc_f(ncnum)%nshadow_vert)          ! push last entry to sindex
		nc_f(ncnum)%pbcx(sindex) = pbcx 
		nc_f(ncnum)%pbcy(sindex) = pbcy                     
		nc_f(ncnum)%nshadow_vert = nc_f(ncnum)%nshadow_vert - 1

	!@ if within the cutoff
	Else If (length .Le. nc_f(ncnum)%shadow_cutoff_sq) Then
	  !@ update if already in shadow
	  If (inshadow) Then
	   sindex = get_shadow_index(ncnum,vert)                                                               
	   nc_f(ncnum)%pbcx(sindex) = pbcx 
	   nc_f(ncnum)%pbcy(sindex) = pbcy                     

	  !@ add to, if not in shadow previously 
	  Else
	   ver(vert)%shadownc(ncnum) = ncnum
	   nc_f(ncnum)%nshadow_vert = nc_f(ncnum)%nshadow_vert + 1 
	   If (nc_f(ncnum)%nshadow_vert .Gt. maxshadowvert) Call Throw_error_message(1)
	   nc_f(ncnum)%shadow_vertices(nc_f(ncnum)%nshadow_vert) = vert              
	   nc_f(ncnum)%pbcx(nc_f(ncnum)%nshadow_vert) = pbcx    
	   nc_f(ncnum)%pbcy(nc_f(ncnum)%nshadow_vert) = pbcy  
	  Endif
        Endif

	Call ncshadow_r(ncnum)                                                                         ! Compute the mean z and mean curvature in the shadow
	Call ncshadow_mcurv(ncnum)             
	End subroutine Compute_nc_shadow_vertices_local
!----------------------------------------------------------------------------------------------------------------------------
!                                Throw an error message     
!----------------------------------------------------------------------------------------------------------------------------
       Subroutine Throw_error_message(messageno) 
        Implicit None
        Integer :: messageno
         If (messageno .Eq. 1) Then
          Print*,"================================>==================>"
          Print*,str_red
          Print*,"Detected a setup problem"
          Print*,'Found an occurrence of number of shadow vertices exceeding ',maxshadowvert
          Print*,'Increase the maximum of maxshadowvertices in module_datastruct.f and restart'
          Print*,'Exiting'
          Print*
          Stop
         Endif
       End Subroutine Throw_error_message
!----------------------------------------------------------------------------------------------------------------------------
!                                Return the area of the membrane
!----------------------------------------------------------------------------------------------------------------------------
       Function get_membrane_area() Result(memarea)
       Implicit None 
       Real(Kind=8) :: memarea
       memarea = Sum(ver(1:nver)%totarea)
       Return 
       End function get_membrane_area

!----------------------------------------------------------------------------------------------------------------------------
!                                COMPUTE MEAN Z-POSITION FOR THE MEMBRANE       
!----------------------------------------------------------------------------------------------------------------------------
       Function Compute_mean_membrane_z() Result(meanz)
       Implicit None 
       Real(Kind=8) :: meanz
       meanz = Sum(ver(1:nver)%vcoord(3,1))/nver
       Return
       End Function Compute_mean_membrane_z
       
!----------------------------------------------------------------------------------------------------------------------------
!                                COMPUTE MEAN Z-POSITION FOR THE MEMBRANE       
!----------------------------------------------------------------------------------------------------------------------------
       Function Compute_closest_vertex(ncnum,xpos,ypos) Result(meanz)
       Use module_randomnumber
       Implicit None 
       Real(Kind=8) :: meanz,xpos,ypos,xll,xlu,yll,ylu
       Integer :: ncnum, listvertex(nver),i,ncounter,rannum

       listvertex=0
       ncounter=0

       xll = xpos - nc_f(ncnum)%soft_radius
       if (xll<0) xll=0
       xlu = xpos + nc_f(ncnum)%soft_radius
       if (xlu>periodic_box_length) xlu = periodic_box_length
       yll = ypos - nc_f(ncnum)%soft_radius
       if (yll<0) yll=0
       ylu = ypos + nc_f(ncnum)%soft_radius
       if (ylu>periodic_box_length) ylu = periodic_box_length

       Do i=1,nver,1
         If ((ver(i)%vcoord(1,1) .Gt. xll) .And. (ver(i)%vcoord(1,1) .Lt. xlu)) Then
          If ((ver(i)%vcoord(2,1) .Gt. yll) .And. (ver(i)%vcoord(2,1) .Lt. ylu)) Then
             ncounter=ncounter+1
             listvertex(ncounter)=i
          Endif
        Endif
       Enddo

       rannum= int(ncounter*ran2(seed)) + 1
       if (rannum .Gt. ncounter) rannum = ncounter
       meanz = ver(listvertex(rannum))%vcoord(3,1)
       Return
       End Function Compute_closest_vertex

!----------------------------------------------------------------------------------------------------------------------------
!                                COMPUTE MEAN Z-POSITION FOR THE MEMBRANE UNDER     THE GIVEN NANOCARRIER
!----------------------------------------------------------------------------------------------------------------------------      
       Subroutine ncshadow_r(ncnum)
       Implicit None 
       Integer :: ncnum


        nc_f(ncnum)%meanr(1,1) = nc_f(ncnum)%coord(1,1)
        nc_f(ncnum)%meanr(2,1) = nc_f(ncnum)%coord(2,1)
        nc_f(ncnum)%meanr(3,1) = membrane_init_z
        nc_f(ncnum)%meandist = Sqrt(Sum((nc_f(ncnum)%coord - nc_f(ncnum)%meanr)**2))        ! PBC is not needed since it is taken care in shadow{x,y,z}
     
       End Subroutine ncshadow_r

!----------------------------------------------------------------------------------------------------------------------------
!                                COMPUTE MEAN Z-POSITION FOR THE MEMBRANE       
!----------------------------------------------------------------------------------------------------------------------------       
       Subroutine ncshadow_mcurv(ncnum)
       Implicit None 
       Integer :: ncnum
       nc_f(ncnum)%meancurv=0.0
       End Subroutine ncshadow_mcurv
!----------------------------------------------------------------------------------------------------------------------------
!                               Reset the shadow vertices
!----------------------------------------------------------------------------------------------------------------------------
       Subroutine  Reset_shadow_vertices(ncnum)
       Implicit None 
       Integer :: ncnum,i
        Do i=1,pbnum,1
         If (ver(i)%shadownc(ncnum) .Eq. ncnum) ver(i)%shadownc(ncnum) = 0
        Enddo
       End Subroutine Reset_shadow_vertices

!----------------------------------------------------------------------------------------------------------------------------
!                                Store the shadow vertices
!----------------------------------------------------------------------------------------------------------------------------
       Subroutine  Store_Restore_shadowvertices(flag,ncnum)
       Implicit None 
       Integer :: flag,ncnum 
       If (flag .Eq. 0) mp%nc_f(ncnum) = nc_f(ncnum)                                                                                 ! Store if flag is 0
       If (flag .Eq. 1) nc_f(ncnum) = mp%nc_f(ncnum)                                                                                 ! Restore if flag is 1
       End Subroutine Store_Restore_shadowvertices      
!----------------------------------------------------------------------------------------------------------------------------
!                                UPDATE MEAN MEMBRANE CURVATURE IN THE SHADOW OF THE NC DURING VERTEX MOVE    
!----------------------------------------------------------------------------------------------------------------------------
       Function Update_meanH_vertex_move(oldmeanH,curmcur,oldmcur,nshadowvert) Result(meanH)
       Implicit None 
       Real(Kind=8) :: oldmeanH,dH,meanH,oldmcur,curmcur
       Integer :: nshadowvert
       dH = 0.0
       dH = curmcur - oldmcur
       meanH = oldmeanH + dH/nshadowvert
       Return
       End Function Update_meanH_vertex_move

!----------------------------------------------------------------------------------------------------------------------------
!                                UPDATE MEAN MEMBRANE CURVATURE IN THE SHADOW OF THE NC
!----------------------------------------------------------------------------------------------------------------------------
       Function Update_meanH(oldmeanH,curmcur,oldmcur,nshadowvert) Result(meanH)
       Implicit None 
       Real(Kind=8) :: oldmeanH,meanH,oldmcur,curmcur
       Integer :: nshadowvert
       meanH = oldmeanH + (curmcur-oldmcur)/nshadowvert
       Return
       End Function Update_meanH
       
!----------------------------------------------------------------------------------------------------------------------------
!                                COMPUTE THE CHANGE IN REACTION ENERGY
!----------------------------------------------------------------------------------------------------------------------------
       Function get_delE() Result(delE)
       Implicit None 
       Real(Kind=8) :: delE
       delE=reac_energy_change*1.38*1E-23*300                                                                            
       Return
       End Function get_delE       

!---------------------------------------------------------------------------------------------------------------------------       
!                 Subroutine to check if a triangle is within the periodic box
!---------------------------------------------------------------------------------------------------------------------------       
	Subroutine check_triangle_in_pbbox()
	Implicit None
	Integer :: i,j,v(3),flag
	tri(:)%boxtriangle=1                                                                                                        ! Initialize such that all triangle are inside
	Do i=1,ntr,1      
	 flag=0 ; v=tri(i)%vert ; j=1
  	 Do While ((j.Le.3) .And. (flag .Eq. 0))
	 If((ver(v(j))%vcoord(1,1) .Gt. periodic_box_length) .Or. (ver(v(j))%vcoord(2,1) .Gt. periodic_box_length))flag=1           ! Check for the coordinates of each vertex making up the triangle
	 j=j+1
	 Enddo
	 If (flag.Eq.1) tri(i)%boxtriangle=0                                                                                        ! If a vertex lies outside the box, the flag is set to zero 
	Enddo
	End Subroutine check_triangle_in_pbbox

!---------------------------------------------------------------------------------------------------------------------------       
!                 Subroutine to check if a vertex is within the periodic box
!---------------------------------------------------------------------------------------------------------------------------       
	Subroutine check_vert_in_pbbox()
	Implicit None
	Integer :: i
	ver(:)%boxvert=1                                                                                                            ! Initialize such that all triangle are inside
	Do i=1,pbnum,1      
	 If((ver(i)%vcoord(1,1) .Gt. periodic_box_length) .Or. (ver(i)%vcoord(2,1) .Gt. periodic_box_length)) ver(i)%boxvert=0      ! Check for the coordinates of each vertex making up the triangle
	Enddo
	End Subroutine check_vert_in_pbbox


!---------------------------------------------------------------------------------------------------------------------------       
!                 Subroutine to compute the volume of the box under the membrane
!---------------------------------------------------------------------------------------------------------------------------       
!       the volume of the individual triangle element is computed in the global coordinate system with origin at (0,0,0)

	Subroutine compute_tri_volume(trno)
	Implicit None
	Integer :: trno
        Real(Kind=8),Dimension(3,1) :: a,b,c    

	 If(tri(trno)%boxtriangle .Eq. 1) Then
          a=ver(tri(trno)%vert(1))%vcoord  
          b=ver(tri(trno)%vert(2))%vcoord  
          c=ver(tri(trno)%vert(3))%vcoord
          tri(trno)%vol=  (a(1,1)*(b(2,1)*c(3,1)-b(3,1)*c(2,1)) + &
                           a(2,1)*(b(3,1)*c(1,1)-b(1,1)*c(3,1)) + &
                           a(3,1)*(b(1,1)*c(2,1)-b(2,1)*c(1,1)))/6.0                   ! (a.(bxc))/6.0 is the volume 
         Endif
	End Subroutine compute_tri_volume

	
!---------------------------------------------------------------------------------------------------------------------------       
!                 Subroutine to compute the volume of the box under the membrane
!---------------------------------------------------------------------------------------------------------------------------       
!      (2/101/14, RN): The volume enclosed between the membrane and box boundary is computed by setting the coordinate system at 
!      the center of mass of the membrane. The cuboid volume is calculated seperately and the volume under the membrane
!      is calculated seperately and these contribution are added up later

	Subroutine compute_total_volume()
	Implicit None
        Integer :: i
        Real(Kind=8),Dimension(3,1) :: a,b,c,com(3,1)
	mp%vol=0.0

    	com(1,1)=Sum(ver(:)%vcoord(1,1))/pbnum	
    	com(2,1)=Sum(ver(:)%vcoord(2,1))/pbnum	
    	com(3,1)=Sum(ver(:)%vcoord(3,1))/pbnum	

        vol_calc:DO i=1,ntr
	 If(tri(i)%boxtriangle .Eq. 1) Then
         a=ver(tri(i)%vert(1))%vcoord-com
         b=ver(tri(i)%vert(2))%vcoord-com
         c=ver(tri(i)%vert(3))%vcoord-com
         mp%vol= mp%vol+ (a(1,1)*(b(2,1)*c(3,1)-b(3,1)*c(2,1)) + &
                          a(2,1)*(b(3,1)*c(1,1)-b(1,1)*c(3,1)) + &
                          a(3,1)*(b(1,1)*c(2,1)-b(2,1)*c(1,1)))/6.0             ! (a.(bxc))/6.0 is the volume 
         Endif
         ENDDO vol_calc
	 Call allocate_cubevertices(periodic_box_length)
	 mp%vol=mp%vol+compute_cuboid_vol(com,1,5,6)
	 mp%vol=mp%vol+compute_cuboid_vol(com,1,6,2)
	 mp%vol=mp%vol+compute_cuboid_vol(com,2,6,7)
	 mp%vol=mp%vol+compute_cuboid_vol(com,2,7,3)
	 mp%vol=mp%vol+compute_cuboid_vol(com,3,7,8)
	 mp%vol=mp%vol+compute_cuboid_vol(com,3,8,4)
	 mp%vol=mp%vol+compute_cuboid_vol(com,4,8,5)
	 mp%vol=mp%vol+compute_cuboid_vol(com,4,5,1)
	 mp%vol=mp%vol+compute_cuboid_vol(com,5,7,6)
	 mp%vol=mp%vol+compute_cuboid_vol(com,5,8,7)
	 mp%vol=mp%vol+compute_cuboid_vol(com,1,2,3)
	 mp%vol=mp%vol+compute_cuboid_vol(com,1,3,4)
	End Subroutine compute_total_volume


!--------------------------------------------------------------------------------------------------------------------------
!                 Subroutine to estimate the volume by dividing into a triangular tube + a tetrahedron
!--------------------------------------------------------------------------------------------------------------------------
	Subroutine compute_membrane_vol_enclosed()
	Implicit None
	Integer :: i,i1(1)
	Real(Kind=8),Dimension(3,1) :: a,b,c,a1,b1,c1,r21,r31
	Real(Kind=8) :: zcoor(3),minz,area

	mp%vol=0.0

	 DO i=1,ntr
	 If(tri(i)%boxtriangle .Eq. 1) Then	
         a=ver(tri(i)%vert(1))%vcoord
         b=ver(tri(i)%vert(2))%vcoord
         c=ver(tri(i)%vert(3))%vcoord
	 zcoor=[a(3,1),b(3,1),c(3,1)]
	 minz=Minval(zcoor) ; i1=Minloc(zcoor)
	 r21=b-a ; r31=c-a
	 area=Sqrt((r21(1,1)*r31(2,1)-r21(2,1)*r31(1,1))**2)*0.5
	 mp%vol=mp%vol+minz*area                                                                                               	       ! Volume of the triangular tube 

	 If (i1(1) .Eq. 1) Then                                                                                                        ! a corresponds to the vertex with minz and will be the reference point for the next calc.  
         a=ver(tri(i)%vert(1))%vcoord ; a1=a
         b=ver(tri(i)%vert(2))%vcoord ; b1=reshape([b(1,1),b(2,1),minz],(/3,1/))
         c=ver(tri(i)%vert(3))%vcoord ; c1=reshape([c(1,1),c(2,1),minz],(/3,1/))
	 Else If (i1(1) .Eq. 2) Then
         a=ver(tri(i)%vert(2))%vcoord ; a1=a
         b=ver(tri(i)%vert(3))%vcoord  ; b1=reshape([b(1,1),b(2,1),minz],(/3,1/))
         c=ver(tri(i)%vert(1))%vcoord  ; c1=reshape([c(1,1),c(2,1),minz],(/3,1/))
	 Else If (i1(1) .Eq. 3) Then
         a=ver(tri(i)%vert(3))%vcoord ; a1=a
         b=ver(tri(i)%vert(1))%vcoord ; b1=reshape([b(1,1),b(2,1),minz],(/3,1/))
         c=ver(tri(i)%vert(2))%vcoord ; c1=reshape([c(1,1),c(2,1),minz],(/3,1/))
	 Endif

	 b=b-a ; b1=b1-a ; c=c-a ; c1=c1-a

         mp%vol= mp%vol+(c(1,1)*(b1(2,1)*c1(3,1)-b1(3,1)*c1(2,1))+ &
		 c(2,1)*(b1(3,1)*c1(1,1)-b1(1,1)*c1(3,1))+c(3,1)*(b1(1,1)*c1(2,1)-b1(2,1)*c1(1,1)))/6.0                             ! (a.(bxc))/6.0 is the volume 

         mp%vol= mp%vol+ (c(1,1)*(b(2,1)*b1(3,1)-b(3,1)*b1(2,1))+ &
		 c(2,1)*(b(3,1)*b1(1,1)-b(1,1)*b1(3,1))+c(3,1)*(b(1,1)*b1(2,1)-b(2,1)*b1(1,1)))/6.0                                 ! (a.(bxc))/6.0 is the volume 

    	 Endif
         ENDDO

	End Subroutine compute_membrane_vol_enclosed

!---------------------------------------------------------------------------------------------------------------------------       
!                 Subroutine to compute the volume enclosed by the cuboid defined by the boundary vertices
!---------------------------------------------------------------------------------------------------------------------------       
	Function compute_cuboid_vol(com,v1,v2,v3) Result(cvol)
	Implicit None
	Real(Kind=8) :: com(3,1),cvol,a(3,1),b(3,1),c(3,1)
	Integer :: v1,v2,v3
        a=geomb(v1)%coord-com
        b=geomb(v2)%coord-com
        c=geomb(v3)%coord-com
        cvol= &
            (a(1,1)*(b(2,1)*c(3,1)-b(3,1)*c(2,1))+a(2,1)*(b(3,1)*c(1,1)-b(1,1)*c(3,1))+a(3,1)*(b(1,1)*c(2,1)-b(2,1)*c(1,1)))/6.0    ! (a.(bxc))/6.0 is the volume 
	Return
	End function compute_cuboid_vol

!-------------------------------------------------------------------------------------------------------------------------------
!                         Subroutine to allocate the position of the 8 vertices required for volume calculations
!-------------------------------------------------------------------------------------------------------------------------------
! eight vertices define the boundary of the box enclosed by the membrane (membrane + box under it)
! These eight vertices are taken and numbered as follows
!            1. (xmin,ymin,zmin)
!            2. (xmin,ymax,zmin)
!            3. (xmax,ymax,zmin)
!            4. (xmax,ymin,zmin)
!            5. (xmin,ymin,zmax)
!            6. (xmin,ymax,zmax)
!            7. (xmax,ymax,zmax)
!            8. (xmax,ymin,zmax)

!!! The triangles formed by these vertices with outward pointing normals are given by
!!!!!!!! (1,5,6) (1,6,2) (2,6,7) (2,7,3) (3,7,8) (3,8,4) (4,8,5) (4,5,1) (5,7,6) (5,8,7) (1,2,3) (1,3,4)

! This data will be used in  total volume calculation to determine the total
! volume enclosed between the membrane and box boundary

	Subroutine allocate_cubevertices(pbl)
	Implicit None
	Real(Kind=8) :: pbl
	geomb(1)%vert=0
	geomb(1)%coord=reshape([0.0,0.0,0.0],(/3,1/))
	geomb(2)%vert=0
	geomb(2)%coord=reshape([0.0,Real(pbl),0.0],(/3,1/))
	geomb(3)%vert=0
	geomb(3)%coord=reshape([Real(pbl),Real(pbl),0.0],(/3,1/))
	geomb(4)%vert=0
	geomb(4)%coord=reshape([Real(pbl),0.0,0.0],(/3,1/))
	geomb(5)%vert=gsize**2+2
	geomb(5)%coord=ver(geomb(5)%vert)%vcoord
	geomb(6)%vert=gsize*(gsize+1)+4
	geomb(6)%coord=ver(geomb(6)%vert)%vcoord
	geomb(7)%vert=gsize**2
	geomb(7)%coord=ver(geomb(7)%vert)%vcoord
	geomb(8)%vert=gsize*(gsize+3)+3
	geomb(8)%coord=ver(geomb(8)%vert)%vcoord
	End Subroutine allocate_cubevertices

      End MODULE module_datastruct


