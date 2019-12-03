
! For use in debugging if an antigen is present in the old and new linklist
! correctly. For use in module_mcsmoves.f/Antigen_diffusion_on_triangle()
    If (ant_no.Eq.32)Then
    Print*,str_red
    Print*,'Computed link cell for antigen ', ant_no
    Call print_cellnumber(ant_no,'c')
    old_curcell=return_linkcelldata('c',ant_no)
    Print*,str_black
    Print*,str_green
    if (check_ifin_linkcell(ant_no,old_curcell,'c') .Eq. 1)Then
	    Print*,'antigen ',ant_no,' is present in cell ',old_curcell
    Else
	    Print*,'antigen ',ant_no,' is not present in cell ',old_curcell
    Endif
    Print*,str_black
    Endif

    If (ant_no.Eq.32)Then
    Print*,str_blue
    Print*,'Computed link cell for antigen ',ant_no 
    Call print_cellnumber(ant_no,'c')
    curcell=return_linkcelldata('c',ant_no)
    Print*,str_black
    Print*,str_green
    if (check_ifin_linkcell(ant_no,curcell,'c') .Eq. 1)Then
	    Print*,'antigen ',ant_no,' is present in cell ',curcell
    Else
	    Print*,'antigen ',ant_no,' is not present in cell it is added to: ',curcell
	    Pause
    Endif
    if (check_ifin_linkcell(ant_no,old_curcell,'c') .Eq. 1)Then
	    Print*,'antigen ',ant_no,' is present in old cell ',old_curcell
    Else
	    Print*,'antigen ',ant_no,' is not present in old cell : ',old_curcell
	    Pause
    Endif

    Print*,str_black
    Endif


! use in link flip module
	  If(an_ver(i,2).Eq.34)Then
	 Print*,'updating link flip ',i,an_ver(i,1),an_ver(i,2)
	 Print*,antig(an_ver(i,2))%base_coord
	  Print*,'completing link flip ',i
	  Print*,'base_coord after'
	  Print*,antig(an_ver(i,2))%base_coord
	  Print*,'Antigen cell number after '
	  Call print_cellnumber(an_ver(i,2),'c')
	  Print*,'Linkcell data for antigen after',an_ver(i,2)
	  Call Print_linkcelldata('c',an_ver(i,2))
	 Print*,'============================================================================'
	  Endif


 	  If (an_ver(i,2) .Eq.34) Then
	  Print*
	  Print*,'============================================================================'
	  Print*,'Linkcell data for antigen before',an_ver(i,2)
	  Print*,'The code reports the current link cell for antigen as'
	  Call print_cellnumber(an_ver(i,2),'c')
	  Call Print_linkcelldata('c',an_ver(i,2))
	  Print*,'------------------------------------------------------------------------'
!	  Print*,'evaluating link flip ',i
!	  Print*,'Triangle vertices',new_tri_no
!	  Print*, tri(new_tri_no)%vert
!	  Print*,'vertex involved ',cver(i)
!	  Print*,'vertex normal'
!	  Print*,ver(cver(i))%vnor
!	  Print*,'Antigen disp ',antig(an_ver(i,2))%disp1,antig(an_ver(i,2))%disp2
	  Print*,'base_coord before'
	  Print*,antig(an_ver(i,2))%base_coord
!	  Print*,'Antigen cell number '
!	  Call print_cellnumber(an_ver(i,2),'c')
	  Endif


! debug the linkcell data before and after the call for antigen_diffusion_on_triangle. Originally used in module_mcsmoves.f

       If(arand  .Eq. 34) Then
       Print*,mcsstep,i, 'Antigen-cellnumber-before'
       Call print_cellnumber(34,'c')
       Call print_linkcelldata('c',34)
       Endif
       If(arand  .Eq. 34) Then
       Print*,mcsstep,i, 'Antigen-cellnumber'
       Call print_cellnumber(34,'c')
       Call print_linkcelldata('c',34)
       Pause
       Endif


!-->
	  If(an_ver(i,2).eq.34)Then
	  curcell=return_linkcelldata('c',an_ver(i,2))
	  Print*,str_red
	  Print*,str_hash
	  Print*,i,an_ver(i,1),an_ver(i,2)
	  Print*,'3 .New cell for antigen ',an_ver(i,2),' is ',old_curcell,int(antig(an_ver(i,2))%base_coord(1,1)*0.012)+&
	  int(antig(an_ver(i,2))%base_coord(2,1)*0.012)*3+int(antig(an_ver(i,2))%base_coord(3,1)*0.012)*9
	  if (check_ifin_linkcell(an_ver(i,2), curcell,'c') .Eq. 1)Then
	  Print*,'antigen ',an_ver(i,2),' is present in cell ',curcell
	  Else
	  Print*,'antigen ',an_ver(i,2),' is not present in cell ',curcell,str_black
	  Endif
  	  Print*,str_blue
	  if (check_ifin_linkcell(an_ver(i,2),old_curcell,'c') .Eq. 1)Then
	  Print*,'antigen ',an_ver(i,2),' is present in cell ',old_curcell
	  Else
	  Print*,'antigen ',an_ver(i,2),' is not present in cell ',old_curcell
	  Endif
	  Print*,str_hash,str_black
	  Endif
!-->

!-->
	  If(an_ver(i,2).Eq.34)Then
	 Print*,'updating link flip ',i,an_ver(i,1),an_ver(i,2)
	 Print*,antig(an_ver(i,2))%base_coord
	  Print*,'completing link flip ',i
	  Print*,'base_coord after'
	  Print*,antig(an_ver(i,2))%base_coord
	  Print*,'Antigen cell number after '
	  Call print_cellnumber(an_ver(i,2),'c')
	  Print*,'Linkcell data for antigen after',an_ver(i,2)
	  Call Print_linkcelldata('c',an_ver(i,2))
	 Print*,'============================================================================'
	  Endif
!-->	  
       

!-->
       If(ANY(tri(new_tri_no)%vert.Eq.0))Then
       Print*,str_green
       Print*,'Triangle new after rearrange ', new_tri_no,old_tri_no
       Print*, tri(new_tri_no)%vert
       Pause
      Endif
      If((old_tri_no.Ne.0) .and. (ANY(tri(old_tri_no)%vert.Eq.0)))Then
       Print*,str_green
       Print*,'Triangle old after rearrange ', old_tri_no,new_tri_no
       PrinT*,'At start ',old_vert_list
       Print*, tri(old_tri_no)%vert
       Pause
      Endif
!-->

!-->
!       make_trial_triangles: If (old_tri_no .Ne. 0) Then
!	 new_trial_triangle=old_tri_no	                                                                                      	    ! All the moves are made in the triangle the antigen is associated with
!	 old_trial_triangle=old_tri_no
!       Else
!	 Do i=1,num_bias_moves,1
!	 nei_no=Nint(ran2(seed)*ver(ver_no)%nonei)+1
!    	 If(nei_no .Gt. ver(ver_no)%nonei) nei_no=ver(ver_no)%nonei
!	 new_trial_triangle(i)=ver(ver_no)%vneitr(nei_no)                                                                            ! Generate num_bias_moves random triangle list
!	 nei_no=Nint(ran2(seed)*ver(ver_no)%nonei)+1
!    	 If(nei_no .Gt. ver(ver_no)%nonei) nei_no=ver(ver_no)%nonei
!	 old_trial_triangle(i)=ver(ver_no)%vneitr(nei_no)
!	 Enddo 
!       Endif make_trial_triangles
!         Do i=1,num_bias_moves,1                                              		                                            ! Choose one triangle randomly from the neighbour list
!      	 new_tri_no=new_trial_triangle(i)											    ! Store the exact triangle where the antigen is currently localized
!         fres=Generate_Antigen_Trial_Position(ver_no,new_tri_no,ant_no,disp)
!         If (fres .Eq. 0) Then
!         bias_disp1(i)=disp(1) ; bias_disp2(i)=disp(2)                                                                              ! Record the displacement made  if the move satisfies self-avoidance
!         bias_energy(i)= antigen_reaction_energy(ant_no-1,antig(ant_no)%tip_coord(1,1))
!         Else
!         bias_disp1(i)=mp%antig(ant_no)%disp1 ; bias_disp2(i)=mp%antig(ant_no)%disp2                                                ! The old position is maintained if the move violates self-avoidance
!         bias_energy(i)= energy_init_state ; 
!         Endif
!         trial_probability(i)=rosbluth_weight_new+exp(-beta*bias_energy(i))                                                         ! compute the rosenbluth factor for the new configuration 
!         rosbluth_weight_new=rosbluth_weight_new+trial_probability(i)                                                               ! compute the rosenbluth factor for the new configuration 
!         Enddo
!         rosbluth_weight_old=rosbluth_weight_old+exp(-beta*energy_init_state)							    ! Generate trial configurations starting from the old position
!         Do i=2,num_bias_moves,1
!      	 new_tri_no=old_trial_triangle(i)											    ! Store the exact triangle where the antigen is currently localized
!         fres=Generate_Antigen_Trial_Position(ver_no,new_tri_no,ant_no,disp)
!         If (fres .Eq. 0) Then
!         energy_trial_state = antigen_reaction_energy(ant_no-1,antig(ant_no)%tip_coord(1,1))
!         Else
!         energy_trial_state = energy_init_state
!         Endif
!         rosbluth_weight_new=rosbluth_weight_new+exp(-beta*energy_trial_state)                                                      ! compute the rosenbluth factor for the old configuration 
!         Enddo
!         selected_trial_move=Select_from_trial_moves(trial_probability,rosbluth_weight_new,num_bias_moves)
!         select_a_trial: If (selected_trial_move .Eq. 0) Then
!          move_rejected_flag=1
!         Else
!          avoid_zero_weight:If(rosbluth_weight_old .Eq. 0.0) Then
!          acc_prob=rosbluth_weight_new/(1.0E-20)
!          Else
!	  acc_prob=rosbluth_weight_new/rosbluth_weight_old
!	 Endif avoid_zero_weight
!	  chosen_triangle = new_trial_triangle(selected_trial_move)
!       	 rosbluth_Metropolis: If(ran2(seed) .Lt.  acc_prob)Then
!	  antig(ant_no)%disp1=bias_disp1(selected_trial_move)
!	  antig(ant_no)%disp2=bias_disp2(selected_trial_move)
!	  fres = Displace_Antigen_on_Triangle(ver_no,chosen_triangle,ant_no)
!          Call update_linkcells(ant_no,&
!	  mp%antig(ant_no)%base_coord(1,1),mp%antig(ant_no)%base_coord(2,1),mp%antig(ant_no)%base_coord(3,1),'c')   		    ! Update the link cell
!	 Else
!	  move_rejected_flag=1 
!         Endif rosbluth_metropolis
!	Endif select_a_trial
!      If(move_rejected_flag .Eq. 1) antig(ant_no)=mp%antig(ant_no) 

!-->
