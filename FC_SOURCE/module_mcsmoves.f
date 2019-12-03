!===========================================================================================================================
!===========================================================================================================================
!===========================================================================================================================
!      $$ modflip                          MODULE TO FLIP THE BONDS And DISPLACE THE CHOSEN VERTEX POINTS
!===========================================================================================================================
!===========================================================================================================================

       MODULE module_mcsmoves
       Integer :: norchk,amov,aflip,tmov,tflip,moveno
       Real(KIND=8):: prob,rnum,spexen                                        ! flen,mven,nemven are temp var for energy cal 
       Real(KIND=8):: dfang
       contains
!---------------------------------------------------------------------------------------------------------------------------
!                                           Subroutine to write the membrane area
!--------------------------------------------------------------------------------------------------------------------------- 
       Subroutine Membrane_MonteCarlo()
       Use Module_datastruct ; Use Module_randomnumber
       Integer :: i,mvinter
       mvinter=NINT(Real(tlink)/nver)
       
       Do i=1,tlink,1
       Call Link_flip()
       If(mod(i,mvinter).Eq.0)Then
       Call Vertex_Move()
       Endif 
       Enddo
       End Subroutine Membrane_MonteCarlo
!---------------------------------------------------------------------------------------------------------------------------
!                                           Subroutine to write the membrane area
!--------------------------------------------------------------------------------------------------------------------------- 
       Subroutine Membrane_MonteCarlo_biased()
       Use Module_datastruct ; Use Module_randomnumber
       Integer :: i,mvinter
       
       mvinter=NINT(Real(tlink)/nver)
       Do i=1,tlink,1
       Call Link_flip_biased()
       If(mod(i,mvinter).Eq.0)Then
       Call Vertex_Move_biased()
       Endif 
       Enddo
       End Subroutine Membrane_MonteCarlo_biased
!---------------------------------------------------------------------------------------------------------------------------
!                                           Subroutine TO FLIP A LINK
!--------------------------------------------------------------------------------------------------------------------------- 
      Subroutine Link_Flip()
      USE module_datastruct ; USE module_curvcalc 
      USE module_makesurface ; Use module_randomnumber 
      IMPLICIT NONE                                                                                                                 ! tdv and tdt are the local data structure
      Include '_mod_mcsmoves_interface.h'
      
      Real(KIND=8) :: delRE,delE,delPr
      Integer :: vt1,vt2,i,j,ep1,ep2,fv1,fv2,blch,i1                                                                                ! vt -->vert of triangle; norchk -->nor check  
      Integer :: fvp1,fvm1,fvp2,fvm2,tn1,tn2,rand,mchk,ant_chk_flag,tn11,tn21       
      Integer :: lp1,lp2,ch1,trn,llist(6)
      Integer,DIMENSION(4) :: cver                                                                                                  ! Temproary variables for case statement        
      Integer,DIMENSION(3) :: t1,t2,tmp,tmp1                                                                                        ! ep -->endpoint of link                        

      fv1=0 ; fv2=0; vt1=0; vt2=0
      
      rand=Nint((1.0-2.0*ran2(seed))*tlink)
      
      If ((abs(rand).GT.0) .And. (lin(rand)%boundary .Eq. 0))Then
       t1=tri(lin(rand)%tr)%vert                                                                                                    ! To read the vertices linked to the triangles 
       t2=tri(lin(-rand)%tr)%vert 
       tn1=lin(rand)%tr ; tn2=lin(-rand)%tr                                                                                         ! tn1 and tn2 are the names of the triangles    
       lp1=lin(rand)%sep(1) ; lp2=lin(rand)%sep(2)
       tn11=tn1 ; tn21=tn2
       ep1=0;ep2=0 ; blch=0 ;  ant_chk_flag=0

        llist=(/tri(tn1)%li,tri(tn2)%li/)

        Do i=1,3  
        If(t1(i).EQ.lp2)Then ; ep1=i ; Endif                                                                                        ! Position of vertex where the chosen bond ends   
        If(t2(i).EQ.lp1)Then ; ep2=i ; Endif                                                                                        ! The position of vert where conj chosen bond ends

        If((t1(i).NE.lp1).And.(t1(i).NE.lp2))then
        vt1=t1(i); fv1=i                                                                                                             ! Position of the free vertex in triang1
        Endif                                                                                                                       ! fv1 and fv2 are the position of the free vert

        If((t2(i).NE.lp2).And.(t2(i).NE.lp1))then
        vt2=t2(i); fv2=i                                                                                                             ! Position of the free vertex in triang2
        Endif 
        Enddo

        cver=(/vt1,vt2,lp1,lp2/)                                                                                                    ! All relevant vertice num in one array

        Do i=1,ver(vt1)%nonei,1 
        If(ver(vt1)%vneipt(i).EQ.vt2) RETURN                                                                                        ! Do not proceed if already connected
        Enddo
         
        blch=blcheck(vt2,vt1)                                                                                                       ! Check for the new bond length        
        If(blch .NE.0) RETURN                                                                                                       ! bond length chk ; blch=0 -> satsified

        If((ver(vt1)%nonei.GE.9).OR.(ver(vt2)%nonei.GE.9))  RETURN                                                                  ! Maximum limit on the number of neighbours
        If((ver(lp1)%nonei.LE.3).And.(ver(lp2)%nonei.LE.3)) RETURN                                                                  ! Minimum limit on the number of neighbours


	!@-----> Store all triangles, vertices and links that are affected
        mp%triangle(tn1)=tri(tn1) ; mp%triangle(tn2)=tri(tn2)
        ForAll(i=1:6) mp%link(llist(i))=lin(llist(i))
	Do i=1,4,1
           mp%vertex(cver(i))=ver(cver(i))                                                                                          ! Original state of involved vertices
	   Do j=1,ver(cver(i))%pbimno,1
	   mp%vertex(ver(cver(i))%pbmap(j)) = ver(ver(cver(i))%pbmap(j))                                                       	    ! All its images
	   Enddo
	   
           If(ver(cver(i))%antigen_flag .Eq.1)Then
	     Do i1=1,ver(cver(i))%nantigen
	        mp%antig(ver(cver(i))%antigen_list(i1))=antig(ver(cver(i))%antigen_list(i1))                                        ! store the old values related to antigens on vertex cver(1)
	     Enddo
	   Endif
	Enddo
	!@----->

        t1(ep1)=vt2 ; t2(ep2)=vt1                                                                                                   ! Free vertices are connected to new trian
        j=ver(cver(3))%nonei                                                                                                        ! lp2 is removed from list of lp1 
        Do i=1,j-1,1
        If(ver(cver(3))%vneipt(i).EQ.cver(4))Then
        ver(cver(3))%vneipt(i:j-1)=ver(cver(3))%vneipt(i+1:j)
        Endif
        Enddo  
        ver(cver(3))%vneipt(j:10)=0                                                                                                 ! All sites beyond the neighbour size=0 
        ver(cver(3))%nonei=j-1

        j=ver(cver(4))%nonei                                                                                                        ! lp1 is removed from the list of lp2 
        Do i=1,j-1,1
        If(ver(cver(4))%vneipt(i).EQ.cver(3))Then
        ver(cver(4))%vneipt(i:j-1)=ver(cver(4))%vneipt(i+1:j)
        Endif  
        Enddo  
        ver(cver(4))%vneipt(j:10)=0                                                                                                 ! All sites above neigh size=0 
        ver(cver(4))%nonei=j-1
  
        tmp1=0  
        Do i=1,3,1                                                                                                                  ! To find the the links whose triangles will change
        If(lin(tri(tn1)%li(i))%sep(2).EQ.vt1) tmp1(1)=tri(tn1)%li(i)
        If(lin(tri(tn2)%li(i))%sep(2).EQ.vt2) tmp1(2)=tri(tn2)%li(i)
        Enddo


        lin(rand)%sep=(/vt2,vt1/)                                                                                                   ! Update the start&end of the chosen random link 
        lin(-rand)%sep=(/vt1,vt2/)                                
        tri(tn1)%vert=t1                                                                                                            ! Update the Vertex of the triangle 
        tri(tn2)%vert=t2  

        fvp1=fv1+1 ; If(fv1.EQ.3)fvp1=1                                                                                             ! Circular boundary conditions
        fvm1=fv1-1 ; If(fv1.EQ.1)fvm1=3
        fvp2=fv2+1 ; If(fv2.EQ.3)fvp2=1
        fvm2=fv2-1 ; If(fv2.EQ.1)fvm2=3

        tmp=0 
        tmp(1)=tri(tn1)%li(fvm1)                                                                                                    ! tmp is a temp array that stores links
        tmp(2)=tri(tn2)%li(fvm2)

        tri(tn1)%li(fvm1)=tri(tn1)%li(fvp1)                                                                                         ! Updating the triangle --link
        tri(tn2)%li(fvm2)=tri(tn2)%li(fvp2)
        tri(tn1)%li(fvp1)=tmp(2)
        tri(tn2)%li(fvp2)=tmp(1)

        lin(tmp1(1))%tr=tn2                                                                                                         ! Updating the link - triangle
        lin(tmp1(2))%tr=tn1

        ver(cver(1))%nonei=ver(cver(1))%nonei+1                                                                                     ! Increment the neigh of ver1 by  1
        ver(cver(2))%nonei=ver(cver(2))%nonei+1                                                                                     ! Increment the neigh of ver2 by  1

        ch1=0 ; i=1
        f1_ver:Do WHILE((i.LE.mp%vertex(cver(1))%nonei).And.(ch1.EQ.0)) 
        If(mp%vertex(cver(1))%vneipt(i).EQ.cver(3))then
        Do j=i+1,mp%vertex(cver(1))%nonei
        ver(cver(1))%vneipt(j+1)=mp%vertex(cver(1))%vneipt(j)                                                                       ! The change in neighbouring  vertex is put in 
        Enddo 
        ver(cver(1))%vneipt(i+1)=cver(2) ; ch1=1
        Endif
        i=i+1
        Enddo f1_ver

        ch1=0 ; i=1
        f2_ver:Do WHILE((i.LE.mp%vertex(cver(2))%nonei).And.(ch1.EQ.0))                                                             ! The vertex order is changed
        If(mp%vertex(cver(2))%vneipt(i).EQ.cver(4))then
         Do j=i+1,mp%vertex(cver(2))%nonei                                                                                          ! Chosen link goes from cver(3) to cver(4)        
         ver(cver(2))%vneipt(j+1)=mp%vertex(cver(2))%vneipt(j)                                                                      !      cver(4)               4                    
         Enddo                                                                                                                      !         *                  *  After flip        
         ver(cver(2))%vneipt(i+1)=cver(1) ;ch1=1                                                                                    !       / | \              /   \ chosen link goes 
         Endif                                                                                                                      !      /  |  \            /     \ from 2 to 1     
         i=i+1                                                                                                                      ! vt1 /   |   \  vt2     /  tn2  \                
        Enddo f2_ver                                                                                                                !(or)*tn1 |tn2 *(or)  1 *---------* 2             
                                                                                                                                    !cver \   |   / cver(2)  \  tn1  /                
        ch1=0 ; i=1                                                                                                                 ! (1)  \  |  /            \     /                 
        f1_tr:Do WHILE((i.LE.mp%vertex(cver(1))%nonei).And.(ch1.EQ.0))                                                              !       \ | /              \   /                  
        If(mp%vertex(cver(1))%vneitr(i).EQ.tn1)Then                                                                                 !         *                  *                    
         Do j=i+1,mp%vertex(cver(1))%nonei                                                                                          !       cver(3)              3                    
         ver(cver(1))%vneitr(j+1)=mp%vertex(cver(1))%vneitr(j)
         Enddo 
         ver(cver(1))%vneitr(i+1)=tn2 ;ch1=1                                                                                        ! The first triangle free vertex is linked to tri 2
         Endif
         i=i+1
        Enddo f1_tr

        ch1=0 ; i=1
        f2_tr: Do WHILE((i.LE.mp%vertex(cver(2))%nonei).And.(ch1.EQ.0))
         If(mp%vertex(cver(2))%vneitr(i).EQ.tn2)Then
         Do j=i+1,mp%vertex(cver(2))%nonei
         ver(cver(2))%vneitr(j+1)=mp%vertex(cver(2))%vneitr(j)                                                                      ! The first triangle free vertex is linked to tri 2
         Enddo 
         ver(cver(2))%vneitr(i+1)=tn1 ;ch1=1
         Endif
         i=i+1
        Enddo f2_tr

        ch1=0 ; i=1
        f3_tr:Do WHILE(i.LE.mp%vertex(cver(3))%nonei .And.(ch1.EQ.0))                                                               ! Updating the neig trian list of the vertex where 
         If(ver(cver(3))%vneitr(i).EQ.tn2)Then                                                                                      ! the chosen bond started earlier(either t1 or t2) 
          Do j=i,ver(cver(3))%nonei                                                                                                 ! removed
          ver(cver(3))%vneitr(j)=ver(cver(3))%vneitr(j+1)
          Enddo                
          If(i.EQ.mp%vertex(cver(3))%nonei)Then                                                                                     ! If tn2 is last in the list it is rearranged for
          ver(cver(3))%vneitr(mp%vertex(cver(3))%nonei)=ver(cver(3))%vneitr(1)                                                      ! matching the link list
          ver(cver(3))%vneitr(1:ver(cver(3))%nonei)=ver(cver(3))%vneitr(2:mp%vertex(cver(3))%nonei)
          Endif
          ch1=1
         Endif 
         i=i+1
        Enddo f3_tr
        ver(cver(3))%vneitr(mp%vertex(cver(3))%nonei:10)=0 
     
        ch1=0 ; i=1
        f4_tr:Do WHILE((i.LE.mp%vertex(cver(4))%nonei).And.(ch1.EQ.0))                                                              ! Updating the neig trian list of the vertex 
         If(ver(cver(4))%vneitr(i).EQ.tn1) Then
         Do j=i,ver(cver(4))%nonei                                 
         ver(cver(4))%vneitr(j)=ver(cver(4))%vneitr(j+1)
         Enddo                
         If(i.EQ.mp%vertex(cver(4))%nonei)Then
         ver(cver(4))%vneitr(mp%vertex(cver(4))%nonei)=ver(cver(4))%vneitr(1)
         ver(cver(4))%vneitr(1:ver(cver(4))%nonei)=ver(cver(4))%vneitr(2:mp%vertex(cver(4))%nonei)
         Endif 
         ch1=1
         Endif
         i=i+1
        Enddo f4_tr
        ver(cver(4))%vneitr(mp%vertex(cver(4))%nonei:10)=0 

       Call onlyarea(tn1)                                                                                                            ! Calculate the area of the two new triangles 
       Call onlyarea(tn2)  

       Do i=1,4,1
       ver(cver(i))%totarea=0
       Do j=1,ver(cver(i))%nonei
       trn=ver(cver(i))%vneitr(j)
       ver(cver(i))%totarea=ver(cver(i))%totarea+tri(trn)%ar/3.0
       Enddo
       Enddo
        
       mchk=0

       Call faceangchk(tn1)                                                                                                         ! Check  angle between tri 1 and neigh faces       
       fang_tr1:If(norchk .EQ.0) Then                                                                                               ! Proceed further only if Yes (Max =150 \degrees)  
       Call faceangchk(tn2)                                                                                                         ! Check for angle between tri2 and its neigh faces 
       fang_tr2 : If(norchk .EQ.0) Then                             

        Do i=1,4,1
        Call normalcalc(cver(i))                                                                                                    ! Compute Curvature at the main vertex
	Call mapvertextoimages(cver(i))                                                                                             ! Propagate its value all over its images
          If(ver(cver(i))%antigen_flag .Eq.1)Then
	     If((i.Eq.1) .And. (mchk.Eq.0)) mchk = Reorient_Antigens_on_VerUTri(cver(1),tn1)                                        ! Since cver(1) has triangle tn1 before and after the move, reorient existing antigens 
             If((i.Eq.2) .And. (mchk.Eq.0)) mchk = Reorient_Antigens_on_VerUTri(cver(2),tn2)                                                 	    ! Same as above for cver(2) and tn2
	     If((i.Eq.3) .And. (mchk.Eq.0)) Then
	        mchk = Reorient_Antigens_on_VerUTri(cver(3),tn1)                                                                    ! same as above for cver(3) and tn1 
	        If(mchk .Eq. 0) mchk = Reorient_Antigens_on_VerUTri1_VerUTri2(cver(3),tn2,tn1)                                      ! tn2 is lost from the triangle list of cver(3). Move antigens (cver(3)Utn2) -> (cver(3)Utn1)
	     Endif
	     If((i.Eq.4) .And. (mchk.Eq.0)) Then
	        mchk = Reorient_Antigens_on_VerUTri(cver(4),tn2)                                                                    ! same as above for cver(4) and tn2
	        If(mchk .Eq. 0) mchk = Reorient_Antigens_on_VerUTri1_VerUTri2(cver(4),tn1,tn2)                                      ! tn1 is lost from the triangle list of cver(4). Move antigens (cver(4)Utn1) -> (cver(4)Utn2) 
	     Endif
	     If(mchk .Eq. 0) mchk=Selfavoidance_VertexAntigens_NeighVertAntigens(cver(i))                                           ! Do a self avoidance check between the antigens related to cver(i)  
          Endif
        Enddo

        If(mchk .Eq.0)Then
	  delPr = pressure*((tri(tn1)%vol+tri(tn2)%vol)-(mp%triangle(tn1)%vol+mp%triangle(tn2)%vol))
  	  delRE = TriangleAntigen_Reaction_Energy_Change(tn1)+TriangleAntigen_Reaction_Energy_Change(tn2)                           ! Compute the reaction energy due to antigens on each triangle 
	  delE = delRE + delPr + compute_linkflip_energies(cver)                                                                    ! the energy change associated with the move
          If((delE.GT.0) .And. (exp(-beta*delE).Lt.ran2(Seed))) mchk=1                                                              ! Metropolis Algorithm
	Endif

     Else
       	mchk=1 
     Endif fang_tr2

     Else
       	mchk=1 
     Endif fang_tr1

     flip_failed:If((mchk .Eq. 1))then
	!@----> Restore the old state	     
       tri(tn1)=mp%triangle(tn1) ; tri(tn2)=mp%triangle(tn2)
    	ForAll (i=1:6) lin(llist(i))=mp%link(llist(i))
        Do i=1,4,1
	   ver(cver(i))=mp%vertex(cver(i))
	   Do j=1,ver(cver(i))%pbimno,1
	   ver(ver(cver(i))%pbmap(j)) = mp%vertex(ver(cver(i))%pbmap(j))
	   Enddo

     	   If(mp%vertex(cver(i))%antigen_flag .Eq. 1) Then
     		Do i1=1,ver(cver(i))%nantigen,1
     	 	antig(ver(cver(i))%antigen_list(i1))=mp%antig(ver(cver(i))%antigen_list(i1))              		            ! Restore the antigen tip position if move is rejected
     		Enddo
     	   Endif
     	Enddo 
	!@---->
	     
    Endif flip_failed

     flip_accp: If(mchk==0)Then                                
	
	!@----> Update link cells for antigen     
        Do i=1,4,1
       	If(ver(cver(i))%antigen_flag .Eq. 1) Call VertexAntigen_Update_Linkcell(cver(i))
        Enddo
	!@---->

	!@----> Update the properties of the shadow
        Do i=1,num_nanocarrier,1
	Call ncshadow_mcurv(i)
        Enddo
	!@---->

      Endif flip_accp 
    Endif

    End Subroutine Link_flip

!---------------------------------------------------------------------------------------------------------------------------
!                                           Subroutine TO FLIP A LINK
!--------------------------------------------------------------------------------------------------------------------------- 
      Subroutine Link_Flip_biased()
      USE module_datastruct ; USE module_curvcalc 
      USE module_makesurface ; Use module_randomnumber 
      IMPLICIT NONE                                                                                                                 ! tdv and tdt are the local data structure
      Include '_mod_mcsmoves_interface.h'
       
      Real(KIND=8) :: delRE,delE,delPr,bias_ener
      Integer :: vt1,vt2,i,j,ep1,ep2,fv1,fv2,blch,i1                                                                                ! vt -->vert of triangle; norchk -->nor check  
      Integer :: fvp1,fvm1,fvp2,fvm2,tn1,tn2,rand,mchk,ant_chk_flag,tn11,tn21       
      Integer :: lp1,lp2,ch1,trn,llist(6)
      Integer,DIMENSION(4) :: cver                                                                                                  ! Temproary variables for case statement        
      Integer,DIMENSION(3) :: t1,t2,tmp,tmp1                                                                                        ! ep -->endpoint of link                        
      Character :: bias_mode


      fv1=0 ; fv2=0 ; bias_ener = 0.0
      vt1 = 0; vt2 =0
      
      rand=Nint((1.0-2.0*ran2(seed))*tlink)

      
      If ((abs(rand).GT.0) .And. (lin(rand)%boundary .Eq. 0))Then
       t1=tri(lin(rand)%tr)%vert                                                                                                    ! To read the vertices linked to the triangles 
       t2=tri(lin(-rand)%tr)%vert 
       tn1=lin(rand)%tr ; tn2=lin(-rand)%tr                                                                                         ! tn1 and tn2 are the names of the triangles    
       lp1=lin(rand)%sep(1) ; lp2=lin(rand)%sep(2)
       tn11=tn1 ; tn21=tn2
       ep1=0;ep2=0 ; blch=0 ;  ant_chk_flag=0

        llist=(/tri(tn1)%li,tri(tn2)%li/)

        Do i=1,3  
        If(t1(i).EQ.lp2)Then ; ep1=i ; Endif                                                                                        ! Position of vertex where the chosen bond ends   
        If(t2(i).EQ.lp1)Then ; ep2=i ; Endif                                                                                        ! The position of vert where conj chosen bond ends

        If((t1(i).NE.lp1).And.(t1(i).NE.lp2))then
        vt1=t1(i);fv1=i                                                                                                             ! Position of the free vertex in triang1
        Endif                                                                                                                       ! fv1 and fv2 are the position of the free vert

        If((t2(i).NE.lp2).And.(t2(i).NE.lp1))then
        vt2=t2(i);fv2=i                                                                                                             ! Position of the free vertex in triang2
        Endif 
        Enddo

        cver=(/vt1,vt2,lp1,lp2/)                                                                                                    ! All relevant vertice num in one array

        Do i=1,ver(vt1)%nonei,1 
        If(ver(vt1)%vneipt(i).EQ.vt2) RETURN                                                                                        ! Do not proceed if already connected
        Enddo
         
        blch=blcheck(vt2,vt1)                                                                                                       ! Check for the new bond length        
        If(blch .NE.0) RETURN                                                                                                       ! bond length chk ; blch=0 -> satsified

        If((ver(vt1)%nonei.GE.9).OR.(ver(vt2)%nonei.GE.9))  RETURN                                                                  ! Maximum limit on the number of neighbours
        If((ver(lp1)%nonei.LE.3).And.(ver(lp2)%nonei.LE.3)) RETURN                                                                  ! Minimum limit on the number of neighbours

        Do i=1,num_nanocarrier,1
        if (nc_f(i)%bias_mode == 'H')  mp%nc_f(i)%meancurv = nc_f(i)%meancurv                                                      ! Store only mean curvature under the shadow of each nanocarrier
        Enddo
        
        !@-----> Store all triangles, vertices and links that are affected
        mp%triangle(tn1)=tri(tn1) ; mp%triangle(tn2)=tri(tn2)
        ForAll(i=1:6) mp%link(llist(i))=lin(llist(i))
        Do i=1,4,1
         mp%vertex(cver(i))=ver(cver(i))                                                                                          ! Original state of involved vertices
         Do j=1,ver(cver(i))%pbimno,1
          mp%vertex(ver(cver(i))%pbmap(j)) = ver(ver(cver(i))%pbmap(j))                                                       	    ! All its images
         Enddo
   
         If(ver(cver(i))%antigen_flag .Eq.1)Then
          Do i1=1,ver(cver(i))%nantigen
           mp%antig(ver(cver(i))%antigen_list(i1))=antig(ver(cver(i))%antigen_list(i1))                                        ! store the old values related to antigens on vertex cver(1)
          Enddo
         Endif
        Enddo
       !@----->

        t1(ep1)=vt2 ; t2(ep2)=vt1                                                                                                   ! Free vertices are connected to new trian

        j=ver(cver(3))%nonei                                                                                                        ! lp2 is removed from list of lp1 
        Do i=1,j-1,1
        If(ver(cver(3))%vneipt(i).EQ.cver(4))Then
        ver(cver(3))%vneipt(i:j-1)=ver(cver(3))%vneipt(i+1:j)
        Endif
        Enddo  
        ver(cver(3))%vneipt(j:10)=0                                                                                                 ! All sites beyond the neighbour size=0 
        ver(cver(3))%nonei=j-1

        j=ver(cver(4))%nonei                                                                                                        ! lp1 is removed from the list of lp2 
        Do i=1,j-1,1
        If(ver(cver(4))%vneipt(i).EQ.cver(3))Then
        ver(cver(4))%vneipt(i:j-1)=ver(cver(4))%vneipt(i+1:j)
        Endif  
        Enddo  
        ver(cver(4))%vneipt(j:10)=0                                                                                                 ! All sites above neigh size=0 
        ver(cver(4))%nonei=j-1
  
        tmp1=0  
        Do i=1,3,1                                                                                                                  ! To find the the links whose triangles will change
        If(lin(tri(tn1)%li(i))%sep(2).EQ.vt1) tmp1(1)=tri(tn1)%li(i)
        If(lin(tri(tn2)%li(i))%sep(2).EQ.vt2) tmp1(2)=tri(tn2)%li(i)
        Enddo


        lin(rand)%sep=(/vt2,vt1/)                                                                                                   ! Update the start&end of the chosen random link 
        lin(-rand)%sep=(/vt1,vt2/)                                
        tri(tn1)%vert=t1                                                                                                            ! Update the Vertex of the triangle 
        tri(tn2)%vert=t2  

        fvp1=fv1+1 ; If(fv1.EQ.3)fvp1=1                                                                                             ! Circular boundary conditions
        fvm1=fv1-1 ; If(fv1.EQ.1)fvm1=3
        fvp2=fv2+1 ; If(fv2.EQ.3)fvp2=1
        fvm2=fv2-1 ; If(fv2.EQ.1)fvm2=3

        tmp=0 
        tmp(1)=tri(tn1)%li(fvm1)                                                                                                    ! tmp is a temp array that stores links
        tmp(2)=tri(tn2)%li(fvm2)

        tri(tn1)%li(fvm1)=tri(tn1)%li(fvp1)                                                                                         ! Updating the triangle --link
        tri(tn2)%li(fvm2)=tri(tn2)%li(fvp2)
        tri(tn1)%li(fvp1)=tmp(2)
        tri(tn2)%li(fvp2)=tmp(1)

        lin(tmp1(1))%tr=tn2                                                                                                         ! Updating the link - triangle
        lin(tmp1(2))%tr=tn1

        ver(cver(1))%nonei=ver(cver(1))%nonei+1                                                                                     ! Increment the neigh of ver1 by  1
        ver(cver(2))%nonei=ver(cver(2))%nonei+1                                                                                     ! Increment the neigh of ver2 by  1

        ch1=0 ; i=1
        f1_ver:Do WHILE((i.LE.mp%vertex(cver(1))%nonei).And.(ch1.EQ.0)) 
        If(mp%vertex(cver(1))%vneipt(i).EQ.cver(3))then
        Do j=i+1,mp%vertex(cver(1))%nonei
        ver(cver(1))%vneipt(j+1)=mp%vertex(cver(1))%vneipt(j)                                                                       ! The change in neighbouring  vertex is put in 
        Enddo 
        ver(cver(1))%vneipt(i+1)=cver(2) ; ch1=1
        Endif
        i=i+1
        Enddo f1_ver

        ch1=0 ; i=1
        f2_ver:Do WHILE((i.LE.mp%vertex(cver(2))%nonei).And.(ch1.EQ.0))                                                             ! The vertex order is changed
        If(mp%vertex(cver(2))%vneipt(i).EQ.cver(4))then
         Do j=i+1,mp%vertex(cver(2))%nonei                                                                                          ! Chosen link goes from cver(3) to cver(4)        
         ver(cver(2))%vneipt(j+1)=mp%vertex(cver(2))%vneipt(j)                                                                      !      cver(4)               4                    
         Enddo                                                                                                                      !         *                  *  After flip        
         ver(cver(2))%vneipt(i+1)=cver(1) ;ch1=1                                                                                    !       / | \              /   \ chosen link goes 
         Endif                                                                                                                      !      /  |  \            /     \ from 2 to 1     
         i=i+1                                                                                                                      ! vt1 /   |   \  vt2     /  tn2  \                
        Enddo f2_ver                                                                                                                !(or)*tn1 |tn2 *(or)  1 *---------* 2             
                                                                                                                                    !cver \   |   / cver(2)  \  tn1  /                
        ch1=0 ; i=1                                                                                                                 ! (1)  \  |  /            \     /                 
        f1_tr:Do WHILE((i.LE.mp%vertex(cver(1))%nonei).And.(ch1.EQ.0))                                                              !       \ | /              \   /                  
        If(mp%vertex(cver(1))%vneitr(i).EQ.tn1)Then                                                                                 !         *                  *                    
         Do j=i+1,mp%vertex(cver(1))%nonei                                                                                          !       cver(3)              3                    
         ver(cver(1))%vneitr(j+1)=mp%vertex(cver(1))%vneitr(j)
         Enddo 
         ver(cver(1))%vneitr(i+1)=tn2 ;ch1=1                                                                                        ! The first triangle free vertex is linked to tri 2
         Endif
         i=i+1
        Enddo f1_tr

        ch1=0 ; i=1
        f2_tr: Do WHILE((i.LE.mp%vertex(cver(2))%nonei).And.(ch1.EQ.0))
         If(mp%vertex(cver(2))%vneitr(i).EQ.tn2)Then
         Do j=i+1,mp%vertex(cver(2))%nonei
         ver(cver(2))%vneitr(j+1)=mp%vertex(cver(2))%vneitr(j)                                                                      ! The first triangle free vertex is linked to tri 2
         Enddo 
         ver(cver(2))%vneitr(i+1)=tn1 ;ch1=1
         Endif
         i=i+1
        Enddo f2_tr

        ch1=0 ; i=1
        f3_tr:Do WHILE(i.LE.mp%vertex(cver(3))%nonei .And.(ch1.EQ.0))                                                               ! Updating the neig trian list of the vertex where 
         If(ver(cver(3))%vneitr(i).EQ.tn2)Then                                                                                      ! the chosen bond started earlier(either t1 or t2) 
          Do j=i,ver(cver(3))%nonei                                                                                                 ! removed
          ver(cver(3))%vneitr(j)=ver(cver(3))%vneitr(j+1)
          Enddo                
          If(i.EQ.mp%vertex(cver(3))%nonei)Then                                                                                     ! If tn2 is last in the list it is rearranged for
          ver(cver(3))%vneitr(mp%vertex(cver(3))%nonei)=ver(cver(3))%vneitr(1)                                                      ! matching the link list
          ver(cver(3))%vneitr(1:ver(cver(3))%nonei)=ver(cver(3))%vneitr(2:mp%vertex(cver(3))%nonei)
          Endif
          ch1=1
         Endif 
         i=i+1
        Enddo f3_tr
        ver(cver(3))%vneitr(mp%vertex(cver(3))%nonei:10)=0 
     
        ch1=0 ; i=1
        f4_tr:Do WHILE((i.LE.mp%vertex(cver(4))%nonei).And.(ch1.EQ.0))                                                              ! Updating the neig trian list of the vertex 
         If(ver(cver(4))%vneitr(i).EQ.tn1) Then
         Do j=i,ver(cver(4))%nonei                                 
         ver(cver(4))%vneitr(j)=ver(cver(4))%vneitr(j+1)
         Enddo                
         If(i.EQ.mp%vertex(cver(4))%nonei)Then
         ver(cver(4))%vneitr(mp%vertex(cver(4))%nonei)=ver(cver(4))%vneitr(1)
         ver(cver(4))%vneitr(1:ver(cver(4))%nonei)=ver(cver(4))%vneitr(2:mp%vertex(cver(4))%nonei)
         Endif 
         ch1=1
         Endif
         i=i+1
        Enddo f4_tr
        ver(cver(4))%vneitr(mp%vertex(cver(4))%nonei:10)=0 

       Call onlyarea(tn1)                                                                                                            ! Calculate the area of the two new triangles 
       Call onlyarea(tn2)  

       Do i=1,4,1
       ver(cver(i))%totarea=0
       Do j=1,ver(cver(i))%nonei
       trn=ver(cver(i))%vneitr(j)
       ver(cver(i))%totarea=ver(cver(i))%totarea+tri(trn)%ar/3.0
       Enddo
       Enddo
        
       mchk=0

       Call faceangchk(tn1)                                                                                                         ! Check  angle between tri 1 and neigh faces       
       fang_tr1:If(norchk .EQ.0) Then                                                                                               ! Proceed further only if Yes (Max =150 \degrees)  
       Call faceangchk(tn2)                                                                                                         ! Check for angle between tri2 and its neigh faces 
       fang_tr2 : If(norchk .EQ.0) Then                             

        Do i=1,4,1
        Call normalcalc(cver(i))                                                                                                    ! Compute Curvature at the main vertex
	Call mapvertextoimages(cver(i))                                                                                             ! Propagate its value all over its images

          If(ver(cver(i))%antigen_flag .Eq.1)Then
	     If((i.Eq.1) .And. (mchk.Eq.0)) mchk = Reorient_Antigens_on_VerUTri(cver(1),tn1)                                        ! Since cver(1) has triangle tn1 before and after the move, reorient existing antigens 
             If((i.Eq.2) .And. (mchk.Eq.0)) mchk = Reorient_Antigens_on_VerUTri(cver(2),tn2)                                                 	    ! Same as above for cver(2) and tn2
	     If((i.Eq.3) .And. (mchk.Eq.0)) Then
	        mchk = Reorient_Antigens_on_VerUTri(cver(3),tn1)                                                                    ! same as above for cver(3) and tn1 
	        If(mchk .Eq. 0) mchk = Reorient_Antigens_on_VerUTri1_VerUTri2(cver(3),tn2,tn1)                                      ! tn2 is lost from the triangle list of cver(3). Move antigens (cver(3)Utn2) -> (cver(3)Utn1)
	     Endif
	     If((i.Eq.4) .And. (mchk.Eq.0)) Then
	        mchk = Reorient_Antigens_on_VerUTri(cver(4),tn2)                                                                    ! same as above for cver(4) and tn2
		If(mchk .Eq. 0) mchk = Reorient_Antigens_on_VerUTri1_VerUTri2(cver(4),tn1,tn2)                                      ! tn1 is lost from the triangle list of cver(4). Move antigens (cver(4)Utn1) -> (cver(4)Utn2) 
	     Endif
	     If(mchk .Eq. 0) mchk=Selfavoidance_VertexAntigens_NeighVertAntigens(cver(i))                                           ! Do a self avoidance check between the antigens related to cver(i)  
          Endif
        Enddo

        If(mchk .Eq.0) Then
          delPr = pressure*((tri(tn1)%vol+tri(tn2)%vol)-(mp%triangle(tn1)%vol+mp%triangle(tn2)%vol))
          delRE = TriangleAntigen_Reaction_Energy_Change(tn1)+TriangleAntigen_Reaction_Energy_Change(tn2)                           ! Compute the reaction energy due to antigens on each triangle 
  
         Do i=1,num_nanocarrier,1
          If(nc_f(i)%bias_mode == 'H') Then
           Call ncshadow_mcurv(i)
           bias_ener = bias_ener + meanH_biasing_potential(i, mp%nc_f(i)%meancurv, nc_f(i)%meancurv)
          Endif
         Enddo
 
         delE=compute_linkflip_energies(cver) + delRE + delPr + bias_ener                                                          ! the energy change associated with the move	  
         If((delE.GT.0) .And. (exp(-beta*delE).Lt.ran2(Seed))) mchk=1                                                              ! Metropolis Algorithm
        Endif

     Else
        mchk=1 
     Endif fang_tr2

     Else
        mchk=1 
     Endif fang_tr1

     flip_failed:If((mchk .Eq. 1))then
       !@----> Restore the old state	     
       tri(tn1)=mp%triangle(tn1) ; tri(tn2)=mp%triangle(tn2)
       ForAll (i=1:6) lin(llist(i))=mp%link(llist(i))
       Do i=1,4,1
        ver(cver(i))=mp%vertex(cver(i))
        Do j=1,ver(cver(i))%pbimno,1
         ver(ver(cver(i))%pbmap(j)) = mp%vertex(ver(cver(i))%pbmap(j))
        Enddo
        
        If(mp%vertex(cver(i))%antigen_flag .Eq. 1) Then
         Do i1=1,ver(cver(i))%nantigen,1
         antig(ver(cver(i))%antigen_list(i1))=mp%antig(ver(cver(i))%antigen_list(i1))              ! Restore the antigen tip position if move is rejected
         Enddo
        Endif
       Enddo 

       !@---->

     If(bias_mode == 'H') nc_f(1:num_nanocarrier)%meancurv = mp%nc_f(1:num_nanocarrier)%meancurv   ! Restore the curvature around each nanocarrier if the move is rejected

     Endif flip_failed

     flip_accp: If(mchk==0)Then                                
       Do i=1,4,1
        If(ver(cver(i))%antigen_flag .Eq. 1) Call VertexAntigen_Update_Linkcell(cver(i))
       Enddo
      Endif flip_accp 
    Endif

    End Subroutine Link_flip_biased

!---------------------------------------------------------------------------------------------------------------------------
!                  $vermov           Subroutine TO MOVE THE VERTEX
!---------------------------------------------------------------------------------------------------------------------------
! The move subroutine has been written in such a way that whenever a  boundary node is displaced its corrresponding periodic image 
! is also given the same displacement. While area of the periodic triangles are calculated by normal means the surface quantifiers 
! of the periodic vertex image are got from its original vertex

      Subroutine Vertex_Move()
      USE module_datastruct; USE module_curvcalc ; USE module_makesurface ; USE module_writedata ;   Use module_randomnumber
      Use, Intrinsic :: ISO_C_BINDING 
      IMPLICIT NONE 
      Include '_mod_mcsmoves_interface.h'
      Integer::vert,i,j,k,trno,blch,imgv,opbmap(3),mchk,i1,j1                                                                       ! i,j,k,trno -> temp var, *ch-> check vals
      Integer :: an_ver(11,4)
      Real(KIND=8)::inen,fien,delE,delRE,delPr                                                                                      ! init and fin ener, delE, del reaction ener, meanz before move, change in pressure contrib.        
      Real(KIND=8),DIMENSION(3,1)::dr                                                                                               
      Logical :: self_avoidance_flag

      self_avoidance_flag = .False.

      reac_energy_change=0.0
      vermov_attempt=vermov_attempt+1
      If(vermov_attempt .Eq. vermov_step_adjust) Call Reset_vertexmove_Stepsize()                                                   ! Recalculate step size for every 1000 attempts

      vert=NINT(nver*ran2(Seed))+1
      If(vert.Le.nver) Then
      inen=0; fien=0 ; an_ver=0 

      If(ver(vert)%boundary.Eq.0)Then 
      dr=1.0-2.0*reshape((/ran2(Seed),ran2(Seed),ran2(Seed)/),(/3,1/))                                                              ! A Small displacement vector
      Else
      If (fixed_frame_flag .Eqv. .True.) Then
	      Return
      Else
	      dr=reshape((/zero,zero,1.-2.*ran2(Seed)/),(/3,1/))                                                                    ! The boundary vertex moves only along z
      Endif ; Endif

      dr=dr*vermov_step_size

      !@----> Store the state of the chosen vertex and update the position
      mp%vertex(vert)=ver(vert)                                                                                                     ! Store the data for the chosen vertex
      Do i1=1,ver(vert)%nantigen,1
      mp%antig(ver(vert)%antigen_list(i1))=antig(ver(vert)%antigen_list(i1))
      Enddo
      ver(vert)%vcoord=ver(vert)%vcoord+dr                                                                                          ! New displaced position                     

       If(ver(vert)%vcoord(3,1) .Lt. 0) Then
       ver(vert)%vcoord=mp%vertex(vert)%vcoord
       RETURN 
       Endif

       nei_len_chk:Do i=1,ver(vert)%nonei
       blch=blcheck(vert,ver(vert)%vneipt(i))
       If(blch .NE.0 ) then                       
       ver(vert)%vcoord=mp%vertex(vert)%vcoord
       RETURN 
       Endif                                                                                                                        ! with all neighbouring vertices.
       Enddo nei_len_chk

       If(Check_minimum_distance(vert).NE.0) Then                                                                                   !Check for minimum self avoidance 
       ver(vert)%vcoord=mp%vertex(vert)%vcoord
       RETURN 
       Endif

       If(Selfavoidance_membrane_Nanocarrier(vert) .NE. 0) Then                                                                       ! Check Self avoidance of vertex with Nanocarrier 
       ver(vert)%vcoord=mp%vertex(vert)%vcoord
       RETURN 
       Endif

       vert_PBimage:If(ver(vert)%boundary .Eq.1) Then                                                                               ! If the chosen vertex is on the boundary      
       opbmap=ver(vert)%pbmap
       PB_ver:Do i=1,ver(vert)%pbimno                                                                                               ! All possible images are taken care of        
       imgv=ver(vert)%pbmap(i)                                                                                                      ! The image values are stored for further use  
       mp%vertex(imgv)=ver(imgv)
       ver(imgv)%vcoord=ver(imgv)%vcoord+dr                                                                                         ! Update the position of all periodic images   

       nei_len_imgchk:Do k=1,ver(imgv)%nonei,1                                                                                      ! Over all neighbouring vertices. Does not changes with time  
       blch=blcheck(imgv,ver(imgv)%vneipt(k))
       If(blch .NE.0 ) then                       
       Do j=1,i,1                                                                                                                   ! Restore the coordinates of all vertices that have been checked 
       ver(opbmap(j))%vcoord=mp%vertex(opbmap(j))%vcoord
       Enddo
       ver(vert)%vcoord=mp%vertex(vert)%vcoord                                                                                      ! Replace the original coordinates of the chosen vertex 
       RETURN 
       Endif                                                                                                                        ! with all neighbouring vertices.
       Enddo nei_len_imgchk

       If(Check_minimum_distance(imgv) .NE.0) Then                                                                                  ! Check for minimum length const with all verts
       Do j=1,i,1                                                                                                                   ! Restore the coordinates of all vertices that have been checked 
       ver(opbmap(j))%vcoord=mp%vertex(opbmap(j))%vcoord 
       Enddo  
       ver(vert)%vcoord=mp%vertex(vert)%vcoord                                                                                      ! Replace the original coordinates of the chosen vertex 
       RETURN  
       Endif                                                                                                                        !to proceed further else replace the original val

       Enddo PB_ver
       Endif vert_PBimage


	!@----> Store the starting state for the chosen vertex, neighbours and images

       vert_antigen:If((ver(vert)%antigen_flag .Eq.1)) Then             
       an_ver(1,1)=1 ;  an_ver(1,2)=vert                                                                                            ! Collect the antigen flag and antigen vertices 
       Do i1=1,ver(vert)%nantigen,1
       mp%antig(ver(vert)%antigen_list(i1))=antig(ver(vert)%antigen_list(i1))
       Enddo
       Endif  vert_antigen

       Do i=1,ver(vert)%nonei,1                                                                                                     ! Store all neighbouring vertices and triangles
       j=ver(vert)%vneitr(i) ; k=ver(vert)%vneipt(i)
       mp%vertex(k)=ver(k) ;  mp%triangle(j)=tri(j)
       nei_antigen: If((ver(k)%antigen_flag .Eq. 1))Then                 
           Do i1=1,ver(k)%nantigen,1
  	   mp%antig(ver(k)%antigen_list(i1))=antig(ver(k)%antigen_list(i1))
 	   Enddo
 	   an_ver(i+1,1)=1 ; an_ver(i+1,2)=k                                                                                ! Collect the antigen flag and antigen vertices 
       Else
 	   an_ver(i+1,:)=0
       Endif nei_antigen
       Enddo

       Do i=1,ver(vert)%pbimno,1
        imgv = ver(vert)%pbmap(i) 
        Do j=1,ver(imgv)%nonei,1                                                                                                       ! Store all neighbouring vertices and triangle for image vertex v1
        k=ver(imgv)%vneipt(j) ; mp%vertex(k)=ver(k) 
        If (j .Lt. ver(imgv)%nonei) Then                                                                               	            ! All image vertices have nonei-1 triangles
        j1=ver(imgv)%vneitr(j) ; mp%triangle(j1)=tri(j1)
        Endif
       Enddo
       Enddo
       !@---->       


       !@------> Area and volume calculations
       delPr=0.0
       upd_ntr_area:  Do i=1,ver(vert)%nonei                                                                                        ! Update the area of all neighbouring triangles  
       call areacalc(ver(vert)%vneitr(i))                                                                                           ! Calculate the new area of each face            
       delPr=delPr+tri(ver(vert)%vneitr(i))%vol-mp%triangle(ver(vert)%vneitr(i))%vol                                                ! Pressure contribution from each term
       Enddo upd_ntr_area


       Do i=1,ver(vert)%pbimno,1
       imgv = ver(vert)%pbmap(i) 
       upd_ntr_area_img:  Do j=1,ver(imgv)%nonei-1                                                                                  ! Update the area of all neighbouring triangles  
       call areacalc(ver(imgv)%vneitr(j))                                                                                           ! Calculate the new area of each face            
       delPr=delPr+tri(ver(imgv)%vneitr(j))%vol-mp%triangle(ver(imgv)%vneitr(j))%vol                                                ! Pressure contribution from each term
       Enddo upd_ntr_area_img
       Enddo
       !@------>

       
       !@-----> Face angle checks to prevent acute conformations
       norchk=0 ; mchk=0 ; trno=1

       face_angle_check: Do While((trno.Le.ver(vert)%nonei) .And. (norchk.Eq.0))
	       If(tri(ver(vert)%vneitr(trno))%pbflag.Eq.0) Then                                                                     ! For submarginal trinagles call faceangchk()  
		       Call faceangchk(ver(vert)%vneitr(trno)) 
	       Else
		       Call faceangchkPBtr(ver(vert)%vneitr(trno))                                                                  ! For the marginal triangle call faceangchkPBtr()
	       Endif
	       If(norchk.Ne.0) Then
		       mchk=1   ; Exit  
	       Endif                                                                                                                ! set mchk=0 if face angle check fails  
	       trno=trno+1
       Enddo face_angle_check
       !@----->

       !@-----> Curvature calculation at the chosen vertex and neighbours and the computed values are propogated to the images
       vert_update_antigen:If(mchk.Eq.0) Then
       Call normalcalc(vert)                                                                                                        ! mchk=0 is the condition for proceeding with further calculations  
       Call mapvertextoimages(vert)
       If(an_ver(1,1).Eq.1) mchk = Reorient_Vertex_linked_Antigens(an_ver(1,2))                                                     ! Try reorienting all antigens sitting on vertex v (prechecks self avoidance with NC)
       Endif vert_update_antigen

       i=1
       vertnei_update_antigen: Do While( (mchk.Eq.0) .And. (i.Le.ver(vert)%nonei) )
        j=ver(vert)%vneipt(i)       
        If(j .Gt. nver) j = ver(j)%imver                                                                                           ! Call only those vertices which are not periodic images 
         Call normalcalc(j)     
         Call mapvertextoimages(j)                                                                                                  ! map the original vertex props to the image        
         If(an_ver(i+1,1).Eq.1) mchk = Reorient_Vertex_linked_Antigens(an_ver(i+1,2))                                               ! Reorient all antigens on neighbouring vertex i 
         i=i+1
       Enddo vertnei_update_antigen
       

      i=1
      Do While((i.Le.ver(vert)%nonei+1) .And. (mchk.Eq.0))
      If(an_ver(i,1) .Eq. 1)  mchk= Selfavoidance_Antigens_In_VertexNeighbourhood(an_ver(i,2))                                      ! Check self avoidance if the vertex has an antigen
	 i=i+1
      Enddo
      !@----->        

      !@-----> Point where all extra constraints are found to be satisfied      
      self_avoidance_satisfied: If(mchk .Eq.0)Then
       fien = vermovener(vert)                                                                                                      ! Call for final energy calculation

       delRE=0.0
       Do i=1,ver(vert)%nonei+1,1
       If (an_ver(i,1).Eq.1) delRE = delRE+VertexAntigen_Reaction_Energy_Change(an_ver(i,2))                                        ! Compute the change in the reaction energy  (-1 takes care of indices in C++)
       Enddo

       delE=compute_vertexmove_energies(vert)+pressure*delPr+delRE                                                                  ! Total change in energy of the system due to vertex move                     

       move_metropolis: If(exp(-beta*delE) .Ge. ran2(Seed)) Then                                                                    ! Metropolis Scheme to reject a move                                          
	 vermov_accepted=vermov_accepted+1                                                                                   	    ! Append one to the accepted move counter
	 reac_energy_change=delRE
         Call Update_linkcells(vert,mp%vertex(vert)%vcoord(1,1),mp%vertex(vert)%vcoord(2,1),mp%vertex(vert)%vcoord(3,1),'m')        ! Update the link cells after the vertex move is accepted
         If(an_ver(1,1) .Eq. 1) Call VertexAntigen_Update_Linkcell(an_ver(1,2))
         Do i=1,ver(vert)%nonei,1
          If(an_ver(i+1,1) .Eq. 1) Call VertexAntigen_Update_Linkcell(an_ver(i+1,2))
         Enddo

         Do i=1,num_nanocarrier,1
         Call Compute_nc_shadow_vertices_local(i,vert)                                                                                         ! Recompute all shadow vertices again
         Enddo

       Else
       mchk=1                                                                                                                       ! Rejected from Metropolis 
       Endif move_metropolis
       Endif self_avoidance_satisfied

       If (mchk .Eq.1) Then
	 !@----> Restore the chosen vertex 
         ver(vert)=mp%vertex(vert)                                                                                                  !vertex

	 !@----> Restore the antigen at the chosen vertex
	 If(an_ver(1,1).Eq.1) Then
  	 Do i1=1,ver(an_ver(1,2))%nantigen,1                                                                   
  	 j1=ver(an_ver(1,2))%antigen_list(i1)
  	 antig(j1)=mp%antig(j1)                                                                                                     ! Restore the neighbour antigen base and tip position if move is rejected
  	 Enddo
         Endif
	 !@---->

	 !@----> Restore all neighbours and the associated antigens	 
	 rest_neig: Do i=1,ver(vert)%nonei,1                                                                                        ! Update all neighbours       
         j=ver(vert)%vneipt(i) ; k=ver(vert)%vneitr(i)                                                                              ! Neigh points  and triangles 
         ver(j)=mp%vertex(j)   ; tri(k)=mp%triangle(k)
         If(an_ver(i+1,1).Eq.1) Then                                       
  	 Do i1=1,ver(an_ver(i+1,2))%nantigen,1                                                                   
  	 j1=ver(an_ver(i+1,2))%antigen_list(i1)
  	 antig(j1)=mp%antig(j1)                                                                                                     ! Restore the neighbour antigen base and tip position if move is rejected
  	 Enddo
         Endif
         Enddo rest_neig
	 !@---->

	 !@----> Restore all image vertices and their neighbours
	 Do i=1,ver(vert)%pbimno,1
	    imgv = ver(vert)%pbmap(i)
	    ver(imgv) = mp%vertex(imgv)                                                                                             !images   
           Do k=1,ver(imgv)%nonei,1
	       ver(ver(imgv)%vneipt(k))=mp%vertex(ver(imgv)%vneipt(k)) 	                                                            !neighbours of images
       	       If(k .Lt. ver(imgv)%nonei) tri(ver(imgv)%vneitr(k))=mp%triangle(ver(imgv)%vneitr(k)) 
       	   Enddo
         Enddo
	 !@---->

       Endif
       Endif 
       End Subroutine Vertex_Move
    

!---------------------------------------------------------------------------------------------------------------------------
!                  $vermov           Subroutine TO MOVE THE VERTEX
!---------------------------------------------------------------------------------------------------------------------------
! The move subroutine has been written in such a way that whenever a  boundary node is displaced its corrresponding periodic image 
! is also given the same displacement. While area of the periodic triangles are calculated in the usual manner the surface quantifiers 
! of the perriodic vertex image are got from its original vertex

      Subroutine Vertex_Move_biased()
      USE module_datastruct; USE module_curvcalc ; USE module_makesurface ; USE module_writedata ;   Use module_randomnumber
      Use, Intrinsic :: ISO_C_BINDING 
      IMPLICIT NONE 
      Include '_mod_mcsmoves_interface.h'
      Integer :: vert,i,j,k,trno,blch,imgv,opbmap(3),mchk,i1,j1                                                                       ! i,j,k,trno -> temp var, *ch-> check vals
      Integer :: an_ver(11,4)
      Real(KIND=8):: delE,delRE,delPr,delbias_ener                                                                                   ! init and fin ener, delE, del reaction ener, meanz before move, change in pressure contrib.        
      Real(KIND=8),DIMENSION(3,1) :: dr
      Character :: bias_mode      
      Logical :: self_avoidance_flag

      self_avoidance_flag = .False.

      reac_energy_change=0.0 ; delbias_ener =0.0 
      vermov_attempt=vermov_attempt+1
      If(vermov_attempt .Eq. vermov_step_adjust) Call Reset_vertexmove_Stepsize()                                                   ! Recalculate step size for every 1000 attempts

      vert=NINT(nver*ran2(Seed))+1

      If(vert.Le.nver) Then
      an_ver=0 

      If(ver(vert)%boundary.Eq.0)Then 
      dr=1.0-2.0*reshape((/ran2(Seed),ran2(Seed),ran2(Seed)/),(/3,1/))                                                              ! A Small displacement vector
      Else
      If (fixed_frame_flag .Eqv. .True.) Then
      Return
      Else
      dr=reshape((/zero,zero,1.-2.*ran2(Seed)/),(/3,1/))                                                                           ! The boundary vertex moves only along z
      Endif ; Endif

      dr=dr*vermov_step_size

      !@----> Store the state of the chosen vertex and update the position
      mp%vertex(vert)=ver(vert)                                                                                                     ! Store the data for the chosen vertex
      Do i1=1,ver(vert)%nantigen,1
      mp%antig(ver(vert)%antigen_list(i1)) = antig(ver(vert)%antigen_list(i1))
      Enddo
      ver(vert)%vcoord = ver(vert)%vcoord+dr                                                                                          ! New displaced position                     

       If(ver(vert)%vcoord(3,1) .Lt. 0) Then
       ver(vert)%vcoord = mp%vertex(vert)%vcoord
       RETURN 
       Endif

       nei_len_chk:Do i=1,ver(vert)%nonei
       blch=blcheck(vert,ver(vert)%vneipt(i))
       If(blch .NE.0 ) then                       
       ver(vert)%vcoord = mp%vertex(vert)%vcoord
       RETURN 
       Endif                                                                                                                        ! with all neighbouring vertices.
       Enddo nei_len_chk

       If(Check_minimum_distance(vert).NE.0) Then                                                                                   !Check for minimum self avoidance 
       ver(vert)%vcoord = mp%vertex(vert)%vcoord
       RETURN 
       Endif

       If(Selfavoidance_membrane_Nanocarrier(vert) .NE. 0) Then                                                                       ! Check Self avoidance of vertex with Nanocarrier 
       ver(vert)%vcoord = mp%vertex(vert)%vcoord
       RETURN 
       Endif

       vert_PBimage:If(ver(vert)%boundary .Eq.1) Then                                                                               ! If the chosen vertex is on the boundary      
       opbmap=ver(vert)%pbmap
       PB_ver:Do i=1,ver(vert)%pbimno                                                                                               ! All possible images are taken care of        
       imgv=ver(vert)%pbmap(i)                                                                                                      ! The image values are stored for further use  
       mp%vertex(imgv)=ver(imgv)
       ver(imgv)%vcoord=ver(imgv)%vcoord+dr                                                                                         ! Update the position of all periodic images   

       nei_len_imgchk:Do k=1,ver(imgv)%nonei,1                                                                                      ! Over all neighbouring vertices. Does not changes with time  
       blch=blcheck(imgv,ver(imgv)%vneipt(k))
       If(blch .NE.0 ) then                       
       Do j=1,i,1                                                                                                                   ! Restore the coordinates of all vertices that have been checked 
       ver(opbmap(j))%vcoord=mp%vertex(opbmap(j))%vcoord
       Enddo
       ver(vert)%vcoord=mp%vertex(vert)%vcoord                                                                                      ! Replace the original coordinates of the chosen vertex 
       RETURN 
       Endif                                                                                                                        ! with all neighbouring vertices.
       Enddo nei_len_imgchk

       If(Check_minimum_distance(imgv) .NE.0) Then                                                                                  ! Check for minimum length const with all verts
       Do j=1,i,1                                                                                                                   ! Restore the coordinates of all vertices that have been checked 
       ver(opbmap(j))%vcoord=mp%vertex(opbmap(j))%vcoord 
       Enddo  
       ver(vert)%vcoord=mp%vertex(vert)%vcoord                                                                                      ! Replace the original coordinates of the chosen vertex 
       RETURN  
       Endif                                                                                                                        !to proceed further else replace the original val

       Enddo PB_ver
       Endif vert_PBimage


       !@----> Store the starting state for the chosen vertex, neighbours and images

       vert_antigen: If((ver(vert)%antigen_flag .Eq. 1)) Then             
       an_ver(1,1)=1 ;  an_ver(1,2)=vert                                                                                            ! Collect the antigen flag and antigen vertices 
       Do i1=1,ver(vert)%nantigen,1
       mp%antig(ver(vert)%antigen_list(i1))=antig(ver(vert)%antigen_list(i1))
       Enddo
       Endif  vert_antigen

       Do i=1,ver(vert)%nonei,1                                                                                                     ! Store all neighbouring vertices and triangles
        j=ver(vert)%vneitr(i) ; k=ver(vert)%vneipt(i)
        mp%vertex(k)=ver(k) ;  mp%triangle(j)=tri(j)
        nei_antigen: If((ver(k)%antigen_flag .Eq. 1))Then                 
         Do i1=1,ver(k)%nantigen,1
          mp%antig(ver(k)%antigen_list(i1)) = antig(ver(k)%antigen_list(i1))
         Enddo
        an_ver(i+1,1)=1 ; an_ver(i+1,2)=k                                                                                        ! Collect the antigen flag and antigen vertices 
        Else
         an_ver(i+1,:)=0
        Endif nei_antigen
       Enddo

       Do i=1,ver(vert)%pbimno,1
        imgv = ver(vert)%pbmap(i) 
        Do j=1,ver(imgv)%nonei,1                                                                                                    ! Store all neighbouring vertices and triangle for image vertex v1
        k=ver(imgv)%vneipt(j) ; mp%vertex(k)=ver(k) 
        If (j .Lt. ver(imgv)%nonei) Then                                                                               	            ! All image vertices have nonei-1 triangles
        j1=ver(imgv)%vneitr(j) ; mp%triangle(j1)=tri(j1)
        Endif
       Enddo
       Enddo
       !@---->       

       !@------> Area and volume calculations
       delPr=0.0
       upd_ntr_area:  Do i=1,ver(vert)%nonei                                                                                        ! Update the area of all neighbouring triangles  
       call areacalc(ver(vert)%vneitr(i))                                                                                           ! Calculate the new area of each face            
       delPr=delPr+tri(ver(vert)%vneitr(i))%vol-mp%triangle(ver(vert)%vneitr(i))%vol                                                ! Pressure contribution from each term
       Enddo upd_ntr_area

       Do i=1,ver(vert)%pbimno,1
       imgv = ver(vert)%pbmap(i) 
       upd_ntr_area_img:  Do j=1,ver(imgv)%nonei-1                                                                                  ! Update the area of all neighbouring triangles  
       call areacalc(ver(imgv)%vneitr(j))                                                                                           ! Calculate the new area of each face            
       delPr=delPr+tri(ver(imgv)%vneitr(j))%vol-mp%triangle(ver(imgv)%vneitr(j))%vol                                                ! Pressure contribution from each term
       Enddo upd_ntr_area_img
       Enddo
       !@------>
       

       !@-----> Face angle checks to prevent acute conformations
       norchk=0 ; mchk=0 ; trno=1
       face_angle_check: Do While((trno.Le.ver(vert)%nonei) .And. (norchk.Eq.0))
	If(tri(ver(vert)%vneitr(trno))%pbflag.Eq.0) Then                                                                     ! For submarginal trinagles call faceangchk()  
	 Call faceangchk(ver(vert)%vneitr(trno)) 
	Else
	 Call faceangchkPBtr(ver(vert)%vneitr(trno))                                                                  ! For the marginal triangle call faceangchkPBtr()
	Endif
	If(norchk.Ne.0) Then
	 mchk=1   ; Exit  
	Endif                                                                                                                ! set mchk=0 if face angle check fails  
	trno=trno+1
       Enddo face_angle_check
       !@----->

       !@-----> Curvature calculation at the chosen vertex and neighbours and the computed values are propogated to the images
       vert_update_antigen:If(mchk.Eq.0) Then
       Call normalcalc(vert)                                                                                                        ! mchk=0 is the condition for proceeding with further calculations  
       Call mapvertextoimages(vert)
       If(an_ver(1,1).Eq.1) mchk = Reorient_Vertex_linked_Antigens(an_ver(1,2))                                                     ! Try reorienting all antigens sitting on vertex v (prechecks self avoidance with NC)
       Endif vert_update_antigen

       i=1
       vertnei_update_antigen: Do While( (mchk .Eq. 0) .And. (i .Le. ver(vert)%nonei) )
         j=ver(vert)%vneipt(i)       
         If(j .Gt. nver) j = ver(j)%imver                                                                                           ! Call only those vertices which are not periodic images 
         Call normalcalc(j)     
         Call mapvertextoimages(j)                                                                                                  ! map the original vertex props to the image        
         If(an_ver(i+1,1) .Eq. 1) mchk = Reorient_Vertex_linked_Antigens(an_ver(i+1,2))                                               ! Reorient all antigens on neighbouring vertex i 
         i=i+1
       Enddo vertnei_update_antigen
     

      i=1
      Do While((i.Le.ver(vert)%nonei+1) .And. (mchk.Eq.0))
      If(an_ver(i,1) .Eq. 1)  mchk = Selfavoidance_Antigens_In_VertexNeighbourhood(an_ver(i,2))                                      ! Check self avoidance if the vertex has an antigen
       i=i+1
      Enddo
      !@----->  


      If ((mchk .Eq. 0) .And. (any(an_ver(:,1) .Eq. 1))) Then           
       If (Does_Move_breaks_bond(ver(vert)%nonei+1,an_ver(1:ver(vert)%nonei+1,:)) .Eqv. .True.) mchk = 1                            ! Check if a bond is broken during a vertex move (if an antigen is on the vertex)
      Endif


      !@-----> Point where all extra constraints are found to be satisfied
      self_avoidance_satisfied: If(mchk .Eq. 0)Then
       self_avoidance_flag = .True.      
       delRE=0.0

       Do i=1,ver(vert)%nonei+1,1
       If (an_ver(i,1) .Eq. 1) delRE = delRE + VertexAntigen_Reaction_Energy_Change(an_ver(i,2))                                    ! Compute the change in the reaction energy
       Enddo

       Do i=1, num_nanocarrier,1
         mp%nc_f(i) = nc_f(i)                                                                                                       ! Store all details of the nanocarrier shadow
         bias_mode = nc_f(i)%bias_mode;
         Call Compute_nc_shadow_vertices_local(i,vert)                                                                              ! Recompute all shadow vertices again
         If( (bias_mode .Eq. 'Z') .Or. (bias_mode .Eq. 'T') ) &
          delbias_ener = delbias_ener + meanr_biasing_potential(i, mp%nc_f(i)%meanr, nc_f(i)%meanr)                                   ! Compute the change in biasing energy 
         If (bias_mode .Eq. 'H') delbias_ener = delbias_ener + meanH_biasing_potential(i, mp%nc_f(i)%meancurv, nc_f(i)%meancurv)    ! (for H bias) 
       Enddo

       delE = compute_vertexmove_energies(vert)+delRE +pressure*delPr + delbias_ener                                                ! Total change in energy of the system due to vertex move                   
       
        move_metropolis: If(exp(-beta*delE) .Ge. ran2(Seed)) Then                                                                   ! Metropolis Scheme to reject a move                                          
          vermov_accepted=vermov_accepted+1                                                                                         ! Append one to the accepted move counter
          reac_energy_change=delRE
          Call Update_linkcells(vert,mp%vertex(vert)%vcoord(1,1),mp%vertex(vert)%vcoord(2,1),mp%vertex(vert)%vcoord(3,1),'m')       ! Update the link cells after the vertex move is accepted
          If(an_ver(1,1) .Eq. 1) Call VertexAntigen_Update_Linkcell(an_ver(1,2))
          Do i=1,ver(vert)%nonei,1
           If(an_ver(i+1,1) .Eq. 1) Call VertexAntigen_Update_Linkcell(an_ver(i+1,2))
          Enddo
        Else
          mchk=1                                                                                                                   ! Rejected from Metropolis 
        Endif move_metropolis

       Endif self_avoidance_satisfied

       If (mchk .Eq.1) Then

	 !@----> Restore the chosen vertex 
         ver(vert)=mp%vertex(vert)                                                                                                  !vertex

	 !@----> Restore the antigen at the chosen vertex
	 If(an_ver(1,1).Eq.1) Then
  	 Do i1=1,ver(an_ver(1,2))%nantigen,1                                                                   
  	 j1=ver(an_ver(1,2))%antigen_list(i1)
  	 antig(j1)=mp%antig(j1)                                                                                                     ! Restore the neighbour antigen base and tip position if move is rejected
  	 Enddo
         Endif
	 !@---->

	 !@----> Restore all neighbours and the associated antigens	 
	 rest_neig: Do i=1,ver(vert)%nonei,1                                                                                        ! Update all neighbours       
         j=ver(vert)%vneipt(i) ; k=ver(vert)%vneitr(i)                                                                              ! Neigh points  and triangles 
         ver(j)=mp%vertex(j)   ; tri(k)=mp%triangle(k)
         If(an_ver(i+1,1).Eq.1) Then                                       
  	 Do i1=1,ver(an_ver(i+1,2))%nantigen,1                                                                   
  	 j1=ver(an_ver(i+1,2))%antigen_list(i1)
  	 antig(j1)=mp%antig(j1)                                                                                                     ! Restore the neighbour antigen base and tip position if move is rejected
  	 Enddo
         Endif
         Enddo rest_neig
	 !@---->

	 !@----> Restore all image vertices and their neighbours
	 Do i=1,ver(vert)%pbimno,1
	    imgv = ver(vert)%pbmap(i)
	    ver(imgv) = mp%vertex(imgv)                                                                                             !images   
           Do k=1,ver(imgv)%nonei,1
	       ver(ver(imgv)%vneipt(k))=mp%vertex(ver(imgv)%vneipt(k)) 	                                                            !neighbours of images
       	       If(k .Lt. ver(imgv)%nonei) tri(ver(imgv)%vneitr(k))=mp%triangle(ver(imgv)%vneitr(k)) 
       	   Enddo
         Enddo
	 !@---->

	 !@----> Restore the state of the nanocarrier (shadow data)
	 If (self_avoidance_flag) Then                                                                                              ! The ncshadow data is stored only when the self avoidance of antigens is satisfied 
         Do i=1,num_nanocarrier,1
         nc_f(i) = mp%nc_f(i)
         Enddo
	 Endif
	 !@---->

       Endif
       Endif
       End Subroutine Vertex_Move_biased

!---------------------------------------------------------------------------------------------------------------------------
!                                   Subroutine to compute vertex move energies
!---------------------------------------------------------------------------------------------------------------------------
    Function compute_vertexmove_energies(vert) Result(energy)
	Use module_datastruct
	Implicit None
	Integer :: i,j,k,vert,k1
	Real (Kind =8 )::energy

	energy = 0.0
	energy = energy + kappa*((ver(vert)%mcur-ver(vert)%czero)**2*ver(vert)%totarea - &
	                (mp%vertex(vert)%mcur-mp%vertex(vert)%czero)**2*mp%vertex(vert)%totarea)

	Do i=1,ver(vert)%pbimno,1
	   j=ver(vert)%pbmap(i)
	   energy = energy + kappa*((ver(j)%mcur-ver(j)%czero)**2*ver(j)%totarea - &
	                            (mp%vertex(j)%mcur-mp%vertex(j)%czero)**2*mp%vertex(j)%totarea)
        Enddo

	Do k=1,ver(vert)%nonei,1
	k1 = ver(vert)%vneipt(k)
	if (k1 .Gt. nver) k1 = ver(k1)%imver
	energy = energy + kappa*((ver(k1)%mcur-ver(k1)%czero)**2*ver(k1)%totarea - &
	                            (mp%vertex(k1)%mcur-mp%vertex(k1)%czero)**2*mp%vertex(k1)%totarea)
	Do i=1,ver(k1)%pbimno,1
	   j=ver(k1)%pbmap(i)
	   energy = energy + kappa*((ver(j)%mcur-ver(j)%czero)**2*ver(j)%totarea - &
	                            (mp%vertex(j)%mcur-mp%vertex(j)%czero)**2*mp%vertex(j)%totarea)
        Enddo
	Enddo

	Return
	End Function compute_vertexmove_energies


!---------------------------------------------------------------------------------------------------------------------------
!                                   Subroutine to compute vertex move energies
!---------------------------------------------------------------------------------------------------------------------------
    Function compute_linkflip_energies(vert) Result(energy)
	Use module_datastruct
	Implicit None
	Integer :: i,j,k,vert(4),k1
	Real (Kind =8 )::energy
	energy = 0.0

	Do k=1,4,1
	k1 = vert(k)
	energy = energy + kappa*((ver(k1)%mcur-ver(k1)%czero)**2*ver(k1)%totarea -&
	                            (mp%vertex(k1)%mcur-mp%vertex(k1)%czero)**2*mp%vertex(k1)%totarea)
	Do i=1,ver(k1)%pbimno,1
	   j=ver(k1)%pbmap(i)
	   energy = energy + kappa*((ver(j)%mcur-ver(j)%czero)**2*ver(j)%totarea -&
	                            (mp%vertex(j)%mcur-mp%vertex(j)%czero)**2*mp%vertex(j)%totarea)
        Enddo
	Enddo

	Return

	End Function compute_linkflip_energies
       
!---------------------------------------------------------------------------------------------------------------------------
!           Compute change in reaction energy for all antigens linked to a   vertex                 
!---------------------------------------------------------------------------------------------------------------------------
    Subroutine VertexAntigen_Update_Linkcell(vertex_no)
    Use module_datastruct
    Implicit None
    Include '_mod_mcsmoves_interface.h'
    Integer :: vertex_no,antigen_no
    Real(Kind=8) :: delRE
    Integer :: i
    
    delRE=0.0
    Do i=1,ver(vertex_no)%nantigen,1
    antigen_no=ver(vertex_no)%antigen_list(i)
    Call Update_linkcells(antigen_no, &
    mp%antig(antigen_no)%base_coord(1,1),mp%antig(antigen_no)%base_coord(2,1),mp%antig(antigen_no)%base_coord(3,1),'c')             ! Update the link cell for antigens if the moved vertex has an antigen  
    Enddo
    Return
    
    End Subroutine VertexAntigen_Update_Linkcell
!---------------------------------------------------------------------------------------------------------------------------
!           Compute change in reaction energy for all antigens linked to a   vertex                 
!---------------------------------------------------------------------------------------------------------------------------
    Function VertexAntigen_Reaction_Energy_Change(vertex_no)Result(delRE)
    Use module_datastruct
    Implicit None
    Include '_mod_mcsmoves_interface.h'
    Integer :: vertex_no,antigen_no
    Real(Kind=8) :: delRE
    Integer :: i
    
    delRE=0.0
    Do i=1,ver(vertex_no)%nantigen,1
    antigen_no=ver(vertex_no)%antigen_list(i)
    If((bondedstate(antigen_no-1) .Eqv. .TRUE.)) delRE=delRE +&
     antigen_reaction_E_change(antigen_no-1,antig(antigen_no)%tip_coord(1,1),mp%antig(antigen_no)%tip_coord(1,1))
    Enddo
    Return
    End Function VertexAntigen_Reaction_Energy_Change
!---------------------------------------------------------------------------------------------------------------------------
!           Check if a vertex move breaks the antigen-antibody bond
!---------------------------------------------------------------------------------------------------------------------------
   Function Does_Move_breaks_bond(nvertices,an_ver) Result(breakflag)
    Use module_datastruct
    Implicit None
    Include '_mod_mcsmoves_interface.h'
    Integer :: nvertices,antigen_no,an_ver(nvertices,4)
    Logical :: breakflag
    Integer :: i,i1
    
    breakflag = .False.
    Do i=1,nvertices,1    
     If (an_ver(i,1) .Eq. 1) Then
      Do i1 = 1,ver(an_ver(i,2))%nantigen,1	
       antigen_no=ver(an_ver(i,2))%antigen_list(i1)-1
       If (bondedstate(antigen_no)) Then	 
        If((does_bond_breaks(antigen_no) .Eqv. .True.))Then
         breakflag = .True.
         Return
        Endif
       Endif
      Enddo
     Endif
    Enddo
    Return
   End Function Does_Move_breaks_bond
!---------------------------------------------------------------------------------------------------------------------------
!           Compute change in reaction energy for all antigens linked to a   vertex                 
!---------------------------------------------------------------------------------------------------------------------------
    Function TriangleAntigen_Reaction_Energy_Change(triang_no)Result(delRE)
    Use module_datastruct
    Implicit None
    Include '_mod_mcsmoves_interface.h'
    Integer :: triang_no,antigen_no
    Real(Kind=8) :: delRE
    Integer :: i
    
    delRE=0.0
    Do i=1,tri(triang_no)%nantigen,1
     antigen_no = tri(triang_no)%antigen_list(i)
     If((bondedstate(antigen_no-1) .Eqv. .TRUE.)) delRE = delRE + &
     antigen_reaction_E_change(antigen_no-1,antig(antigen_no)%tip_coord(1,1),mp%antig(antigen_no)%tip_coord(1,1))
    Enddo
    Return
    End Function TriangleAntigen_Reaction_Energy_Change

!---------------------------------------------------------------------------------------------------------------------------
!                  Reorient all antigens linked to a vertex depending on whether it is on a vertex or triangle
!---------------------------------------------------------------------------------------------------------------------------
    Function Reorient_Vertex_linked_Antigens(vertex_no)Result(fres)
    Use module_datastruct
    Implicit None
    Integer :: vertex_no,antigen_no,fres,diffus_tri
    Integer :: i
    i=1 ; fres=0
    
    Do While ((i .Le. ver(vertex_no)%nantigen) .And. (fres.Eq.0))
      antigen_no = ver(vertex_no)%antigen_list(i)
      diffus_tri = antig(antigen_no)%diffus_tri
      fres = Compute_Antigen_Base_Tip_Coord(vertex_no,diffus_tri,antigen_no)
      i=i+1
    Enddo
    
    If(fres.Eq.0) fres = Selfavoidance_VertexAntigens_Nanocarrier(vertex_no)                                                        ! Check antigen selfavoidance with nanocarriers
    Return
    End Function Reorient_Vertex_linked_Antigens

!---------------------------------------------------------------------------------------------------------------------------
!                  Reorient all antigens linked to a vertex but not on triangles
!---------------------------------------------------------------------------------------------------------------------------
    Function Reorient_Antigens_on_Vertex(vertex_no)Result(fres)
    Use module_datastruct
    Implicit None
    Integer :: vertex_no,antigen_no,fres,diffus_tri
    Integer :: i
    fres=0 ; i=1
    Do While((i .Le. ver(vertex_no)%nantigen) .And. (fres.Eq.0))
     antigen_no = ver(vertex_no)%antigen_list(i)
     diffus_tri = antig(antigen_no)%diffus_tri
     If(diffus_tri .Eq. 0) fres = Compute_Antigen_Base_Tip_Coord(vertex_no,diffus_tri,antigen_no)
     i=i+1
    Enddo
    If(fres .Eq. 0) fres = Selfavoidance_VertexAntigens_Nanocarrier(vertex_no)                                                      ! Check antigen selfavoidance with nanocarriers
    Return
    End Function Reorient_Antigens_on_Vertex

!---------------------------------------------------------------------------------------------------------------------------
!                                 Reorient all antigens linked to a vertex and a triangle 
!              (Read as all antigens that are on triangle(triang_no)and  also linked to vertex(vertex_no)
!---------------------------------------------------------------------------------------------------------------------------
    Function Reorient_Antigens_on_VerUTri(vertex_no,triang_no) Result(fres)
    Use module_datastruct
    Implicit None
    Integer :: i,vertex_no,triang_no,antigen_no,fres,diffus_tri
    fres=0 ; i=1

    Do While((i .Le.ver(vertex_no)%nantigen) .And. (fres .Eq. 0))
     antigen_no=ver(vertex_no)%antigen_list(i)
     diffus_tri=antig(antigen_no)%diffus_tri
     If(diffus_tri .Eq. triang_no) fres = Compute_Antigen_Base_Tip_Coord(vertex_no,triang_no,antigen_no)
     i=i+1
    Enddo
    If (fres .Eq. 0) fres=Selfavoidance_TriangleAntigens_Nanocarrier(triang_no)                                                    ! Check antigen selfavoidance with nanocarriers
    Return
    End Function Reorient_Antigens_on_VerUTri
!---------------------------------------------------------------------------------------------------------------------------
!              Reorient all antigens linked to a vertex on triangle 1 (old) to vertex on triangle 2 (new)
!---------------------------------------------------------------------------------------------------------------------------
    Function Reorient_Antigens_on_VerUTri1_VerUTri2(vertex_no,triang_no1,triang_no2) Result(fres)
    Use module_datastruct
    Implicit None
    Integer :: vertex_no,triang_no1,triang_no2,antigen_no,fres,diffus_tri
    Integer :: i,i1,ant_index
    
    fres=0 ; i=1 ; ant_index=0
    
    Do While((i .Le. ver(vertex_no)%nantigen) .And. (fres .Eq. 0))
     antigen_no=ver(vertex_no)%antigen_list(i)
     diffus_tri=antig(antigen_no)%diffus_tri
     If(diffus_tri .Eq. triang_no1)Then
      Do i1=1,tri(triang_no1)%nantigen,1
       If(tri(triang_no1)%antigen_list(i1).Eq. antigen_no) ant_index=i1                                                    ! Identify position of chosen antigen in the antigen list of triangle
      Enddo
      
      Do i1=ant_index,tri(triang_no1)%nantigen-1,1
       tri(triang_no1)%antigen_list(i1)=tri(triang_no1)%antigen_list(i1+1)                                                ! Remove the chosen antigen from the list of triang_no1  
       tri(triang_no1)%antigen_list(i1+1)=0
      Enddo
      
      tri(triang_no1)%nantigen=tri(triang_no1)%nantigen-1						   		    ! Decrease the number of antigens in triang_no1 by 1
      tri(triang_no2)%nantigen=tri(triang_no2)%nantigen+1	                                           	            ! Increment the number of antigen in triang_no2 by 1 and add antigen_no to its list
      If(tri(triang_no2)%nantigen .Gt. maxantigen_tri)Then
       fres=1 ; Return ; Endif
       tri(triang_no2)%antigen_list(tri(triang_no2)%nantigen)=antigen_no 
       antig(antigen_no)%diffus_tri=triang_no2									   	    ! Move antigen to new triangle
       fres=Compute_Antigen_Base_Tip_Coord(vertex_no,triang_no2,antigen_no)
      Endif
      i=i+1
    Enddo
    If(fres .Eq. 0) fres = Selfavoidance_TriangleAntigens_Nanocarrier(triang_no2)                                                   ! Check antigen selfavoidance with nanocarriers
    Return
    
    End Function Reorient_Antigens_on_VerUTri1_VerUTri2
!---------------------------------------------------------------------------------------------------------------------------
!                  Reorient all antigens linked to a vertex depending on whether it is on a vertex or triangle
!---------------------------------------------------------------------------------------------------------------------------
    Function Reorient_Antigens_on_triangle(triang_no) Result(fres)
    Use module_datastruct
    Implicit None
    Integer :: i,vertex_no,triang_no,antigen_no,fres
    fres=0 ;  i=1

    Do While((i.Le.tri(triang_no)%nantigen) .And. (fres.Eq.0))
     antigen_no=tri(triang_no)%antigen_list(i)
     vertex_no=antig(antigen_no)%vertex
     fres=Compute_Antigen_Base_Tip_Coord(vertex_no,triang_no,antigen_no)
     i=i+1
    Enddo
    If(fres .Eq. 0) fres = Selfavoidance_TriangleAntigens_Nanocarrier(triang_no)                                                    ! Check antigen selfavoidance with nanocarriers
    Return
    
    End Function Reorient_Antigens_on_triangle
!---------------------------------------------------------------------------------------------------------------------------
!   Subroutine to move an unbound antigen(antigen_no) freely on triangle(triang_no) with vertex(vertex_no) as reference point
!---------------------------------------------------------------------------------------------------------------------------
    Function Displace_Antigen_on_Triangle(vertex_no,triang_no,antigen_no) Result(fres)
    Use module_datastruct
    Implicit None
    Integer :: vertex_no,triang_no,antigen_no,fres
    Real(Kind=8) :: tr_vec1(3,1),tr_vec2(3,1),ant_length
    fres=0
    
    ant_length = antig(antigen_no)%length
    If (triang_no > 0) Then                                                                                                         ! If the antigen_no is associated with a triang_nole move on it with disp1 and disp2 
       Call Compute_triangle_basis_vectors(vertex_no,triang_no,tr_vec1,tr_vec2)                                                     ! compute the two unit vectors of the sides of a triangle with vertex_no as the reference    
       antig(antigen_no)%base_coord=ver(vertex_no)%vcoord+antig(antigen_no)%disp1*tr_vec1+antig(antigen_no)%disp2*tr_vec2                     
       antig(antigen_no)%tip_coord=antig(antigen_no)%base_coord+ant_length*tri(triang_no)%fnor                                      ! The antigen_no takes the face normal when on a triang_nole
   Else
       antig(antigen_no)%base_coord=ver(vertex_no)%vcoord                                                                           ! If the antigen_no is not associated with a triang_no reorient it on the vertex_noex it resides on
       antig(antigen_no)%tip_coord=antig(antigen_no)%base_coord+ant_length*ver(vertex_no)%vnor                                      ! Right now the flexure of the antigen_no is not restored (Needs some thought)
    Endif

       If(antig(antigen_no)%tip_coord(1,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen_no,1,-periodic_box_length)          ! Apply the periodic boundary condition for x comp
       If(antig(antigen_no)%tip_coord(1,1).Lt.0) Call PBC_Antigen_tip(antigen_no,1,periodic_box_length)                             ! Apply the periodic boundary condition for x comp
       If(antig(antigen_no)%tip_coord(2,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen_no,2,-periodic_box_length)          ! Apply the periodic boundary condition for y comp
       If(antig(antigen_no)%tip_coord(2,1).Lt.0) Call PBC_Antigen_tip(antigen_no,2,periodic_box_length)                             ! Apply the periodic boundary condition for y comp
       If((Antigen_selfavoidance(antigen_no).Eq.1) .Or. (Selfavoidance_Antigen_Nanocarrier(antigen_no).Eq.1)) fres=1                ! Check antigen selfavoidance with other antigens and other nanocarriers 
    Return
    End Function Displace_Antigen_on_Triangle 

!---------------------------------------------------------------------------------------------------------------------------
!                     Subroutine to move the antigen and compute the tip position based on theta and phi
!---------------------------------------------------------------------------------------------------------------------------
    Function Reposition_Antigen_on_Triangle(vertex_no,triang_no,antigen_no) Result(fres)
    Use module_datastruct ; Use module_curvcalc
    Implicit None
    Include '_mod_mcsmoves_interface.h'
    Integer :: vertex_no,triang_no,antigen_no,fres
    Real(Kind=8) :: tr_vec1(3,1),tr_vec2(3,1),normal(3,1),theta,phi,ant_vec_z(3,1),ant_length
    fres=0
    ant_length = antig(antigen_no)%length

    If (triang_no > 0) Then                                                                                                         ! If the antigen_no is associated with a triang_nole move on it with disp1 and disp2 
       normal=tri(triang_no)%fnor
       Call Compute_triangle_basis_vectors(vertex_no,triang_no,tr_vec1,tr_vec2)
       antig(antigen_no)%base_coord=ver(vertex_no)%vcoord+antig(antigen_no)%disp1*tr_vec1+antig(antigen_no)%disp2*tr_vec2                     
       theta=antig(antigen_no)%theta ; phi=antig(antigen_no)%phi
       ant_vec_z(1,1)=ant_length*sin(theta)*cos(phi)                                                                                ! Compute the antigen with respect to the z direction for given theta and phi 
       ant_vec_z(2,1)=ant_length*sin(theta)*sin(phi)
       ant_vec_z(3,1)=ant_length*cos(theta)
       Call Compute_Triangle_Householdermatrix(triang_no)                                                                           ! update the Householder matrix for the triangle 
       ant_vec_z=Matmul(tri(triang_no)%HHM,ant_vec_z)                                                                               ! transform the antigen vector such that the z direction points along the triangle normal 
       antig(antigen_no)%tip_coord=antig(antigen_no)%base_coord+ant_vec_z                                                           ! translate the computed antigen orientation to the antigen base
   Else
       normal=ver(vertex_no)%vnor
       antig(antigen_no)%base_coord=ver(vertex_no)%vcoord                                                                           ! If the antigen_no is not associated with a triang_no reorient it on the vertex_noex it resides on
       theta=antig(antigen_no)%theta ; phi=antig(antigen_no)%phi 
       ant_vec_z(1,1)=ant_length*sin(theta)*cos(phi)                                                                                ! Compute the antigen with respect to the z direction for given theta and phi 
       ant_vec_z(2,1)=ant_length*sin(theta)*sin(phi)
       ant_vec_z(3,1)=ant_length*cos(theta)
       Call Compute_Vertex_Householdermatrix(vertex_no)                                                                             ! update the Householder matrix for the triangle 
       ant_vec_z=Matmul(ver(vertex_no)%HHM,ant_vec_z)                                                                               ! transform the antigen vector such that the z direction points along the triangle normal 
       antig(antigen_no)%tip_coord=antig(antigen_no)%base_coord+ant_vec_z                                                           ! translate the computed antigen orientation to the antigen base
    Endif

       If(antig(antigen_no)%tip_coord(1,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen_no,1,-periodic_box_length)          ! Apply the periodic boundary condition for x comp
       If(antig(antigen_no)%tip_coord(1,1).Lt.0) Call PBC_Antigen_tip(antigen_no,1,periodic_box_length)                             ! Apply the periodic boundary condition for x comp
       If(antig(antigen_no)%tip_coord(2,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen_no,2,-periodic_box_length)          ! Apply the periodic boundary condition for y comp
       If(antig(antigen_no)%tip_coord(2,1).Lt.0) Call PBC_Antigen_tip(antigen_no,2,periodic_box_length)                             ! Apply the periodic boundary condition for y comp
       If((Antigen_selfavoidance(antigen_no).Eq.1) .Or. (Selfavoidance_Antigen_Nanocarrier(antigen_no).Eq.1)) fres=1                ! Check self avoidance of antigen with other antigens and other Nanocarriers 
       
       If ((fres .Eq. 0) .And. (bondedstate(antigen_no-1) .Eqv. .True.))Then                                                        ! Additional check if the antigen displaced is bound to an antibody on a vesicle  
        If(does_bond_breaks(antigen_no-1) .Eqv. .True. ) fres=1                                                              ! If the move breaks a bond the function returns .True.
       Endif
    Return
    End Function Reposition_Antigen_on_Triangle

!---------------------------------------------------------------------------------------------------------------------------
!                     Subroutine to move the antigen and compute the tip position based on theta and phi
!---------------------------------------------------------------------------------------------------------------------------
    Subroutine Compute_flexed_base_tip_coord(antigen_no)
    Use module_datastruct ; Use module_curvcalc
    Implicit None
    Include '_mod_mcsmoves_interface.h'
    Integer :: vertex_no,triang_no,antigen_no
    Real(Kind=8) :: tr_vec1(3,1),tr_vec2(3,1),normal(3,1),theta,phi,ant_vec_z(3,1),ant_length
    
    triang_no=antig(antigen_no)%diffus_tri
    vertex_no=antig(antigen_no)%vertex
    ant_length = antig(antigen_no)%length
    
    If (triang_no > 0) Then                                                                                                         ! If the antigen_no is associated with a triang_nole move on it with disp1 and disp2 
       normal=tri(triang_no)%fnor
       Call Compute_triangle_basis_vectors(vertex_no,triang_no,tr_vec1,tr_vec2)
       antig(antigen_no)%base_coord = ver(vertex_no)%vcoord + antig(antigen_no)%disp1*tr_vec1 + antig(antigen_no)%disp2*tr_vec2                     
       theta=antig(antigen_no)%theta ; phi=antig(antigen_no)%phi
       ant_vec_z(1,1)=ant_length*sin(theta)*cos(phi)                                                                                ! Compute the antigen with respect to the z direction for given theta and phi 
       ant_vec_z(2,1)=ant_length*sin(theta)*sin(phi)
       ant_vec_z(3,1)=ant_length*cos(theta)
       Call Compute_Triangle_Householdermatrix(triang_no)                                                                           ! update the Householder matrix for the triangle 
       ant_vec_z=Matmul(tri(triang_no)%HHM,ant_vec_z)                                                                               ! transform the antigen vector such that the z direction points along the triangle normal 
       antig(antigen_no)%tip_coord=antig(antigen_no)%base_coord+ant_vec_z                                                           ! translate the computed antigen orientation to the antigen base
   Else
       normal = ver(vertex_no)%vnor
       antig(antigen_no)%base_coord=ver(vertex_no)%vcoord                                                                           ! If the antigen_no is not associated with a triang_no reorient it on the vertex_noex it resides on
       theta = antig(antigen_no)%theta ; phi=antig(antigen_no)%phi
       ant_vec_z(1,1) = ant_length*sin(theta)*cos(phi)                                                                              ! Compute the antigen with respect to the z direction for given theta and phi 
       ant_vec_z(2,1) = ant_length*sin(theta)*sin(phi)
       ant_vec_z(3,1) = ant_length*cos(theta)
       Call Compute_Vertex_Householdermatrix(vertex_no)                                                                             ! update the Householder matrix for the triangle 
       ant_vec_z=Matmul(ver(vertex_no)%HHM,ant_vec_z)                                                                               ! transform the antigen vector such that the z direction points along the triangle normal 
       antig(antigen_no)%tip_coord=antig(antigen_no)%base_coord+ant_vec_z                                                           ! translate the computed antigen orientation to the antigen base
    Endif

       If(antig(antigen_no)%tip_coord(1,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen_no,1,-periodic_box_length)          ! Apply the periodic boundary condition for x comp
       If(antig(antigen_no)%tip_coord(1,1).Lt.0) Call PBC_Antigen_tip(antigen_no,1,periodic_box_length)                             ! Apply the periodic boundary condition for x comp
       If(antig(antigen_no)%tip_coord(2,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen_no,2,-periodic_box_length)          ! Apply the periodic boundary condition for y comp
       If(antig(antigen_no)%tip_coord(2,1).Lt.0) Call PBC_Antigen_tip(antigen_no,2,periodic_box_length)                             ! Apply the periodic boundary condition for y comp

    End Subroutine Compute_flexed_base_tip_coord
!---------------------------------------------------------------------------------------------------------------------------
!                     Subroutine to move the antigen and compute the tip position based on theta and phi
!---------------------------------------------------------------------------------------------------------------------------
    Function Compute_Antigen_Base_Tip_Coord(vertex_no,triang_no,antigen_no) Result(fres)
    Use module_datastruct ; Use module_curvcalc
    Implicit None
    Include '_mod_mcsmoves_interface.h'
    Integer :: vertex_no,triang_no,antigen_no,fres
    Real(Kind=8) :: tr_vec1(3,1),tr_vec2(3,1),normal(3,1),theta,phi,ant_vec_z(3,1),ant_length
    fres=0

    ant_length = antig(antigen_no)%length

    If (triang_no > 0) Then                                                                                                         ! Antigen is on a triangle
       normal=tri(triang_no)%fnor
       Call Compute_triangle_basis_vectors(vertex_no,triang_no,tr_vec1,tr_vec2)
       antig(antigen_no)%base_coord=ver(vertex_no)%vcoord+antig(antigen_no)%disp1*tr_vec1+antig(antigen_no)%disp2*tr_vec2           ! Base position is computed
       If (bondedstate(antigen_no-1) .Eqv. .True.) Then	
	       theta=antig(antigen_no)%theta ; phi=antig(antigen_no)%phi
	       ant_vec_z(1,1)=ant_length*sin(theta)*cos(phi)                                                                        ! Compute the antigen with respect to the z direction for given theta and phi 
	       ant_vec_z(2,1)=ant_length*sin(theta)*sin(phi)
	       ant_vec_z(3,1)=ant_length*cos(theta)
	       Call Compute_Triangle_Householdermatrix(triang_no)                                                                   ! update the Householder matrix for the triangle 
	       ant_vec_z=Matmul(tri(triang_no)%HHM,ant_vec_z)                                                                       ! transform the antigen vector such that the z direction points along the triangle normal 
	       antig(antigen_no)%tip_coord=antig(antigen_no)%base_coord+ant_vec_z                                                   ! translate the computed antigen orientation to the antigen base
	       Call Impose_PBC_Antigen_coords(antigen_no)
	       If(does_bond_breaks(antigen_no-1) .Eqv. .True.) fres=1
       Else
	       antig(antigen_no)%tip_coord=antig(antigen_no)%base_coord+ant_length*normal                                           ! The unbonded antigen follows the face normal
	       Call Impose_PBC_Antigen_coords(antigen_no)
       Endif   
   Else
       normal=ver(vertex_no)%vnor
       antig(antigen_no)%base_coord=ver(vertex_no)%vcoord                                                                           ! If the antigen_no is not associated with a triang_no reorient it on the vertex_nor it resides on
       If (bondedstate(antigen_no-1) .Eqv. .True.) Then	
	       theta=antig(antigen_no)%theta ; phi=antig(antigen_no)%phi
      	       ant_vec_z(1,1)=ant_length*sin(theta)*cos(phi)                                                                        ! Compute the antigen with respect to the z direction for given theta and phi 
	       ant_vec_z(2,1)=ant_length*sin(theta)*sin(phi)
	       ant_vec_z(3,1)=ant_length*cos(theta)
	       Call Compute_Vertex_Householdermatrix(vertex_no)                                                                     ! update the Householder matrix for the triangle 
	       ant_vec_z=Matmul(ver(vertex_no)%HHM,ant_vec_z)                                                                       ! transform the antigen vector such that the z direction points along the triangle normal 
       	       antig(antigen_no)%tip_coord=antig(antigen_no)%base_coord+ant_vec_z                                                   ! translate the computed antigen orientation to the antigen base
	       Call Impose_PBC_Antigen_coords(antigen_no)
	       If(does_bond_breaks(antigen_no-1) .Eqv. .True.) fres=1
       Else
       	       antig(antigen_no)%tip_coord=antig(antigen_no)%base_coord+ant_length*normal                                           ! translate the computed antigen orientation to the antigen base
	       Call Impose_PBC_Antigen_coords(antigen_no)
       Endif
    Endif
    Return
    End Function Compute_Antigen_Base_Tip_Coord

!---------------------------------------------------------------------------------------------------------------------------
!                              Subroutine to Impose PBC on antigen base and tip coords
!---------------------------------------------------------------------------------------------------------------------------
    Subroutine Impose_PBC_Antigen_coords(antigen_no)
    Use module_datastruct
    Implicit None
    Integer :: antigen_no
       If(antig(antigen_no)%tip_coord(1,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen_no,1,-periodic_box_length)          ! Apply the periodic boundary condition for x comp
       If(antig(antigen_no)%tip_coord(1,1).Lt.0) Call PBC_Antigen_tip(antigen_no,1,periodic_box_length)                             ! Apply the periodic boundary condition for x comp
       If(antig(antigen_no)%tip_coord(2,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen_no,2,-periodic_box_length)          ! Apply the periodic boundary condition for y comp
       If(antig(antigen_no)%tip_coord(2,1).Lt.0) Call PBC_Antigen_tip(antigen_no,2,periodic_box_length)                             ! Apply the periodic boundary condition for y comp
    End Subroutine Impose_PBC_Antigen_coords
!---------------------------------------------------------------------------------------------------------------------------
!                              Subroutine to rearrange antigens between an old and new triangle
!---------------------------------------------------------------------------------------------------------------------------
    Subroutine Rearrange_Antigen_on_triangles(old_tri_no,new_tri_no,ant_no)
    Use module_datastruct
    Implicit None
    Integer :: i,ant_no,new_tri_no,old_tri_no,array_index
    array_index=0

    If(old_tri_no .Ne. new_tri_no) Then
      nonzero_oldtri:If(old_tri_no .Ne. 0) Then
      Do i=1,tri(old_tri_no)%nantigen,1                                                                                             ! find the index of antigen in the triangle old_triang_no
	      If (tri(old_tri_no)%antigen_list(i) .Eq. ant_no) Then
	      array_index=i 
      	      Endif 
      Enddo

      tri(old_tri_no)%antigen_list(array_index)=tri(old_tri_no)%antigen_list(tri(old_tri_no)%nantigen)
      tri(old_tri_no)%antigen_list(tri(old_tri_no)%nantigen) = 0 
      tri(old_tri_no)%nantigen = tri(old_tri_no)%nantigen-1                                                                         ! decrease number of antigen in old triangle by 1
      Endif nonzero_oldtri

      tri(new_tri_no)%nantigen = tri(new_tri_no)%nantigen+1
      tri(new_tri_no)%antigen_list(tri(new_tri_no)%nantigen)=ant_no                                                                 ! Add antigen number to the new triangle it has moved to 
      antig(ant_no)%diffus_tri=new_tri_no
    Endif
    End Subroutine Rearrange_Antigen_on_triangles

!---------------------------------------------------------------------------------------------------------------------------
!                 Subroutine to transfer an antigen from V1UT1 to V2UT2
!---------------------------------------------------------------------------------------------------------------------------
    Subroutine Transfer_Antigen_from_VT1_to_VT2(ant_no,old_vert_no,old_tri_no,new_vert_no,new_tri_no)
    Use module_datastruct
    Implicit None
    Integer :: i,ant_no,old_vert_no,new_vert_no,new_tri_no,old_tri_no,array_index,nantigen
      array_index = 0
      nantigen=ver(old_vert_no)%nantigen
      Do i=1,nantigen,1
      If(ver(old_vert_no)%antigen_list(i) .Eq. ant_no) array_index=i
      Enddo	
  
      ver(old_vert_no)%antigen_list(array_index)=ver(old_vert_no)%antigen_list(nantigen)                                            ! transfer contents of nantigen'th' list element to ith list element
      ver(old_vert_no)%antigen_list(nantigen)=0                                                                                     ! set the nantigen'th' list element to 0   
      ver(old_vert_no)%nantigen=ver(old_vert_no)%nantigen-1
      If (ver(old_vert_no)%nantigen .Eq. 0)  ver(old_vert_no)%antigen_flag=0                                                        ! Set the antigen flag to zero if all antigens are removed

  
      If(old_tri_no .Ne. 0) Then
  	    nantigen=tri(old_tri_no)%nantigen
  	    Do i=1,nantigen,1                                                                                                       ! find the index of antigen in the triangle old_triang_no
              If (tri(old_tri_no)%antigen_list(i) .Eq. ant_no) array_index = i 
  	    Enddo
      	    tri(old_tri_no)%antigen_list(array_index)=tri(old_tri_no)%antigen_list(nantigen)
              tri(old_tri_no)%antigen_list(nantigen)=0
  	    tri(old_tri_no)%nantigen = tri(old_tri_no)%nantigen-1                                                                   ! decrease number of antigen in old triangle by 1
      Endif
  
      ver(new_vert_no)%nantigen = ver(new_vert_no)%nantigen+1
      ver(new_vert_no)%antigen_list(ver(new_vert_no)%nantigen)=ant_no                                                               ! Add antigen number to the new triangle it has moved to 
      ver(new_vert_no)%antigen_flag=1                                                                                               ! By default set all antigen flags to 1 when an antigen is added to the vertex  

      If(new_tri_no .Ne. 0) Then
      tri(new_tri_no)%nantigen = tri(new_tri_no)%nantigen+1
      tri(new_tri_no)%antigen_list(tri(new_tri_no)%nantigen)=ant_no                                                                 ! Add antigen number to the new triangle it has moved to 
      antig(ant_no)%diffus_tri=new_tri_no
      antig(ant_no)%vertex = new_vert_no
      Endif
    End Subroutine Transfer_Antigen_from_VT1_to_VT2
!---------------------------------------------------------------------------------------------------------------------------
!                              Subroutine to move an antigen freely on the triangle surface
!---------------------------------------------------------------------------------------------------------------------------
    Subroutine Antigen_diffusion_on_triangle()
    Use module_datastruct ; Use module_randomnumber
    Implicit None
    Include '_mod_mcsmoves_interface.h'
    Include '_mod_writer_interface.h'
    Integer :: ant_no,ver_no,new_tri_no,nei_no,old_tri_no,move_rejected_flag
    Real(Kind=8) :: ran_disp1,ran_disp2

    ant_no=Nint(ran2(seed)*num_antigens)+1
    If(ant_no .Le. num_antigens) Then
    move_rejected_flag=0 ; ver_no=antig(ant_no)%vertex ; reac_energy_change=0.0
    old_tri_no=antig(ant_no)%diffus_tri ; mp%antig(ant_no)=antig(ant_no)

    bonded_unbonded: If (bondedstate(ant_no-1) .Eqv. .False.) Then                                                                  ! If the antigen is not bonded, perform an unbiased Monte Carlo move. 
      If(old_tri_no .Ne. 0) Then
	       new_tri_no=old_tri_no
		   ran_disp1= (1.0-2.0*ran2(seed))*0.1
	       antig(ant_no)%disp1=antig(ant_no)%disp1+ran_disp1
	       If (antig(ant_no)%disp1 .Gt.1.0) antig(ant_no)%disp1=antig(ant_no)%disp1-ran_disp1                                   ! reverse the direction of movement if the displacement moves the antigen out of triangle
	       ran_disp2= (1.0-2.0*ran2(seed))*0.1
	       antig(ant_no)%disp2=antig(ant_no)%disp2+ran_disp2
	      If (antig(ant_no)%disp2 .Gt.1.0) antig(ant_no)%disp2=antig(ant_no)%disp2-ran_disp2                                    ! reverse the direction of movement if the displacement moves the antigen out of triangle
	       If ((antig(ant_no)%disp1 .Lt. 0) .Or. (antig(ant_no)%disp2 .Lt. 0))Then
		       nei_no=Nint(ran2(seed)*ver(ver_no)%nonei)+1
			   If(nei_no .Gt. ver(ver_no)%nonei) nei_no=ver(ver_no)%nonei                                                   ! Choose one triangle randomly from the neighbour list
		       new_tri_no=ver(ver_no)%vneitr(nei_no)
		       If((tri(new_tri_no)%nantigen .Ge. maxantigen_tri) .Or. (tri(new_tri_no)%pbflag .Eq. 1)) Then                 !  No diffusion on triangles containing periodic image vertices
			       antig(ant_no)%disp1=mp%antig(ant_no)%disp1
			       antig(ant_no)%disp2=mp%antig(ant_no)%disp2
			       Return
		       Endif
		       antig(ant_no)%diffus_tri=new_tri_no                                                                          ! Store the exact triangle where the antigen is currently localized
		       antig(ant_no)%disp1=ran2(seed)*0.1 ;  antig(ant_no)%disp2=ran2(seed)*0.1
	       Endif
      Else
	       nei_no=Nint(ran2(seed)*ver(ver_no)%nonei)+1
	       If(nei_no .Gt. ver(ver_no)%nonei) nei_no=ver(ver_no)%nonei                                                           ! Choose one triangle randomly from the neighbour list
	       new_tri_no=ver(ver_no)%vneitr(nei_no) 
	       If((tri(new_tri_no)%pbflag .Eq. 1) .Or. (tri(new_tri_no)%nantigen .Ge. maxantigen_tri)) Return                       ! Diffusion on triangles containing peridic vertices is not allowed 
	       antig(ant_no)%diffus_tri=new_tri_no                                                                                  ! Store the exact triangle where the antigen is currently localized
	       antig(ant_no)%disp1=ran2(seed)*0.1 ;  antig(ant_no)%disp2=ran2(seed)*0.1
      Endif

      If(Displace_Antigen_on_Triangle(ver_no,new_tri_no,ant_no) .Eq. 1) Then
	      antig(ant_no)=mp%antig(ant_no) ;   Return 
      Else
	      Call Rearrange_Antigen_on_triangles(old_tri_no,new_tri_no,ant_no)                                                     ! Reorder the antigen list on triangles and vice-versa
	      Call update_linkcells(ant_no,&
	      mp%antig(ant_no)%base_coord(1,1),mp%antig(ant_no)%base_coord(2,1),mp%antig(ant_no)%base_coord(3,1),'c')               ! Update the link cell
      Endif

    Else
    
      If(Rosenbluth_Sampling_Bound_Antigens(ant_no) .Eq. 1) antig(ant_no) = mp%antig(ant_no)                                          ! Call to Perform Rosenbluth sampling on the bound antigens  

    Endif bonded_unbonded
    Endif

    End Subroutine Antigen_diffusion_on_triangle

!---------------------------------------------------------------------------------------------------------------------------
!                              Subroutine to randomly displace an antigen from a vertex to another vertex
!---------------------------------------------------------------------------------------------------------------------------
    Subroutine Antigen_Hopping_on_vertex()
    Use module_datastruct ; Use module_randomnumber
    Implicit None
    Include '_mod_mcsmoves_interface.h'
    Include '_mod_writer_interface.h'
    Integer :: ant_no,nei_no,move_rejected_flag,move_accp_flag
    Integer :: old_vert_no,new_vert_no,old_tri_no,new_tri_no

    ant_no=min(Nint(ran2(seed)*num_antigens)+1,num_antigens)
    move_rejected_flag=0 ; move_accp_flag = 0
    old_vert_no=antig(ant_no)%vertex 
    old_tri_no=antig(ant_no)%diffus_tri 

    If (bondedstate(ant_no-1) .Eqv. .False.) Then                                                                                   ! If the antigen is not bonded, perform an unbiased Monte Carlo move. 
	       nei_no=min(Nint(ran2(seed)*ver(old_vert_no)%nonei)+1,ver(old_vert_no)%nonei)
	       new_vert_no=ver(old_vert_no)%vneipt(nei_no)                                                                          ! Choose a new vertex in the neighbourhood of old_vertex
	       If(new_vert_no .Gt. nver) new_vert_no = ver(new_vert_no)%imver							    ! If the chosen vertex is a phantom, we use its image

	       If(ver(new_vert_no)%nantigen .Lt. maxantigen_ver) Then
		  nei_no=min(Nint(ran2(seed)*ver(new_vert_no)%nonei)+1,ver(new_vert_no)%nonei)
	 	  new_tri_no=ver(new_vert_no)%vneitr(nei_no)                                                                        ! Choose a new triangle in the neighbourhood of new_vertex

		  If( (tri(new_tri_no)%pbflag .Ne. 1) .And.(tri(new_tri_no)%nantigen .Lt. maxantigen_tri)) Then
			  mp%antig(ant_no)=antig(ant_no)                                                                            ! Store the original state before peforming the operations  
			  antig(ant_no)%disp1 = ran2(Seed)
			  antig(ant_no)%disp2 = ran2(Seed)
			  move_accp_flag = Displace_Antigen_on_Triangle(new_vert_no,new_tri_no,ant_no)
			  If (move_accp_flag .Eq. 1) Then                                                                           ! If the move violates self avoidance then restore the original state  
				  antig(ant_no) = mp%antig(ant_no)
			  Else
		             Call Transfer_Antigen_from_VT1_to_VT2(ant_no,old_vert_no,old_tri_no,new_vert_no,new_tri_no)            ! Else perform some bookkeeping operations and
			     Call  update_linkcells(ant_no,&
	                     mp%antig(ant_no)%base_coord(1,1),mp%antig(ant_no)%base_coord(2,1),mp%antig(ant_no)%base_coord(3,1),'c')! Update the link cell
		          Endif
		  Endif
	       Endif
    Endif 
    End Subroutine Antigen_Hopping_on_vertex
!---------------------------------------------------------------------------------------------------------------------------
!                  Subroutine to generate a large number of antigen position with fixed flexure and bonded state
!---------------------------------------------------------------------------------------------------------------------------
    Function Rosenbluth_Sampling_Bound_Antigens(antigen_no) Result(fres)
    Use module_datastruct ; Use module_randomnumber
    Implicit None
    Include '_mod_mcsmoves_interface.h'    
    Integer :: vertex_no,antigen_no,tr_vertex_no(3),fres,old_triangle_no,new_triangle_no
    Real(Kind=8) :: rosbluth_weight_old,rosbluth_weight_new,acc_prob    
    Real(Kind=8) :: ran_disp1,ran_disp2,odisp1,odisp2
    Real(Kind=8) :: bias_disp1(num_bias_moves),bias_disp2(num_bias_moves),bias_energy(num_bias_moves)
    Real(Kind=8) :: trial_probability(num_bias_moves),disp_mag
    Real(Kind=8) :: energy_init_state,energy_trial_state
    Integer :: i,n,new_triangle_flag,neigh_no,selected_trial_move,trial_disp_triang(num_bias_moves)

    bias_disp1=0.0 ; bias_disp2=0.0 
    odisp1 = antig(antigen_no)%disp1; odisp2 = antig(antigen_no)%disp2
    bias_energy=0.0 ; trial_probability=0.0 ;  trial_disp_triang=0
    disp_mag=0.05 ; rosbluth_weight_new=0.0 ; rosbluth_weight_old=0.0 ; tr_vertex_no=0

    old_triangle_no = antig(antigen_no)%diffus_tri ;  vertex_no = antig(antigen_no)%vertex
    energy_init_state = antigen_reaction_energy(antigen_no-1,antig(antigen_no)%tip_coord(1,1))  ! Compute the reaction energy (-1 takes care of indices in C++)
    n=0
    Do i=1,num_bias_moves,1
     new_triangle_flag=0
     get_tri_flag: If(old_triangle_no .Eq. 0) Then                                              ! Find a new triangle if the  antigen is on a vertex
      new_triangle_flag=1
     Else
      new_triangle_no=old_triangle_no                                                           ! start with the existing triangle
      ran_disp1=disp_mag*(1-2.0*ran2(seed))
      antig(antigen_no)%disp1=odisp1+ran_disp1
      If(antig(antigen_no)%disp1 .Gt. 1.0) antig(antigen_no)%disp1 = odisp1 - ran_disp1  ! reverse translation if beyond the extent of the triangle 
       ran_disp2=disp_mag*(1-2.0*ran2(seed))  
       antig(antigen_no)%disp2=odisp2+ran_disp2
       If(antig(antigen_no)%disp2 .Gt. 1.0) antig(antigen_no)%disp2=odisp2-ran_disp2
       If((antig(antigen_no)%disp1 .Lt. 0.0) .Or. (antig(antigen_no)%disp2 .Lt. 0.0)) new_triangle_flag=1      ! If the translation takes it to another triangle on the same vertex sample for new_triangles 
     Endif get_tri_flag

     get_new_triangle: If (new_triangle_flag .Eq.1) Then
      neigh_no = Nint(ver(vertex_no)%nonei*ran2(Seed))+1
      If(neigh_no .Gt. ver(vertex_no)%nonei) neigh_no=ver(vertex_no)%nonei
      new_triangle_no = ver(vertex_no)%vneitr(neigh_no)
      antig(antigen_no)%disp1 = disp_mag*ran2(Seed)                                              		    ! only positive displacements since we are moving to a new triangle
      antig(antigen_no)%disp2 = disp_mag*ran2(Seed)
     Endif get_new_triangle


     If ((tri(new_triangle_no)%pbflag .Eq. 0) .And. (tri(new_triangle_no)%nantigen .Lt. maxantigen_tri)) Then
      If(Reposition_Antigen_on_Triangle(vertex_no,new_triangle_no,antigen_no) .Eq. 0) Then
       n=n+1
       bias_disp1(n) = antig(antigen_no)%disp1;                                                                    ! Store trial displacement 1 and 2                            
       bias_disp2(n) = antig(antigen_no)%disp2;
       trial_disp_triang(n) = new_triangle_no;
       bias_energy(n) = antigen_reaction_energy(antigen_no-1,antig(antigen_no)%tip_coord(1,1))
       trial_probability(n) = exp(-beta*bias_energy(n))
       rosbluth_weight_new = rosbluth_weight_new + trial_probability(n)
      Endif 
     Endif
    Enddo

    nonzero_trial: If (n .Gt. 0) Then
     selected_trial_move = Select_from_trial_moves(trial_probability(1:n),rosbluth_weight_new,n)                               ! Select a move from the 
     rosbluth_weight_old = rosbluth_weight_old+exp(-beta*energy_init_state)		  			            ! Generate trial configurations starting from the old position

     Do i=2,num_bias_moves,1
       new_triangle_flag=0
       get_tri_flag2: If(old_triangle_no .Eq. 0) Then                                                                      ! Find a new triangle if the  antigen is on a vertex
        new_triangle_flag=1
       Else
        new_triangle_no = old_triangle_no                                                                               ! start with the existing triangle
        ran_disp1 = disp_mag*(1-2.0*ran2(seed))
        antig(antigen_no)%disp1=odisp1+ran_disp1
        If(antig(antigen_no)%disp1 .Gt. 1.0) antig(antigen_no)%disp1 = odisp1-ran_disp1                ! As before reverse translation if it takes beyond the extent of the triangle 
        ran_disp2 = disp_mag*(1-2.0*ran2(seed))  
        antig(antigen_no)%disp2 = odisp2+ran_disp2
        If(antig(antigen_no)%disp2 .Gt. 1.0) antig(antigen_no)%disp2=odisp2-ran_disp2
        If((antig(antigen_no)%disp1 .Lt. 0.0) .Or. (antig(antigen_no)%disp2 .Lt. 0.0)) new_triangle_flag=1            ! If the translation takes it to another triangle on the same vertex sample for new_triangles	
       Endif get_tri_flag2
    
      get_new_triangle2: If (new_triangle_flag .Eq.1) Then
       neigh_no = Nint(ver(vertex_no)%nonei*ran2(Seed))+1
       If(neigh_no .Gt. ver(vertex_no)%nonei) neigh_no=ver(vertex_no)%nonei
        new_triangle_no = ver(vertex_no)%vneitr(neigh_no)
        antig(antigen_no)%disp1 = disp_mag*ran2(Seed)
        antig(antigen_no)%disp2 = disp_mag*ran2(Seed)
       Endif get_new_triangle2

       If ((tri(new_triangle_no)%pbflag .Eq. 0) .And. (tri(new_triangle_no)%nantigen .Lt. maxantigen_tri)) Then
        If (Reposition_Antigen_on_Triangle(vertex_no,new_triangle_no,antigen_no) .Eq. 0) Then                           ! Trial move is performed only on non boundary triangles
         energy_trial_state = antigen_reaction_energy(antigen_no-1,antig(antigen_no)%tip_coord(1,1))
         rosbluth_weight_old = rosbluth_weight_old + exp(-beta*energy_trial_state)
        Endif 
       Endif
     Enddo
    
    If(rosbluth_weight_old .Gt. 0.0) Then
     acc_prob=rosbluth_weight_new/rosbluth_weight_old
    Else
     acc_prob=1.0
    Endif
    
     rosbluth_Metropolis: If(acc_prob .Gt. ran2(Seed))Then
      new_triangle_no = trial_disp_triang(selected_trial_move)
      antig(antigen_no)%disp1=bias_disp1(selected_trial_move)
      antig(antigen_no)%disp2=bias_disp2(selected_trial_move)
      antig(antigen_no)%diffus_tri = new_triangle_no
      fres = Reposition_Antigen_on_Triangle(vertex_no,new_triangle_no,antigen_no)                                              ! Moving antigens to the prescribed position 
      If (old_triangle_no .Ne. new_triangle_no) Call Rearrange_Antigen_on_Triangles(old_triangle_no,new_triangle_no,antigen_no)! Book keeping step
      Call update_linkcells(antigen_no,&
      mp%antig(antigen_no)%base_coord(1,1),mp%antig(antigen_no)%base_coord(2,1),mp%antig(antigen_no)%base_coord(3,1),'c')      ! Update the link cell
      reac_energy_change=bias_energy(selected_trial_move)-energy_init_state                                                    ! Change in the reaction energy ( passed to module datastruct and accessed if needed)   
      fres=0
     Else
      fres=1
     Endif rosbluth_metropolis
   Else
    fres=1
   Endif nonzero_trial
   
   Return
   End Function Rosenbluth_Sampling_Bound_Antigens
!--------------------------------------------------------------------------------------------------------------------------
!   Subroutine to compute unit vector from a given vertex to other two vertices forming the given triangle
!--------------------------------------------------------------------------------------------------------------------------
      Subroutine Compute_triangle_basis_vectors(vertex_no,triang_no,tr_vec1,tr_vec2)
      Use module_datastruct
      Implicit None
      Integer :: i,n,tr_vertex_no(3),vertex_no,triang_no
      Real(Kind=8) :: tr_vec1(3,1),tr_vec2(3,1)
      n=1 ; tr_vec1=0.0 ; tr_vec2=0.0 ; tr_vertex_no=0
       tr_vertex_no(n)=vertex_no
       Do i=1,3,1
	       If(tri(triang_no)%vert(i).Ne. vertex_no) Then
	       n=n+1 ; tr_vertex_no(n)=tri(triang_no)%vert(i)
       	       Endif 
       Enddo
       tr_vec1=0.5*(ver(tr_vertex_no(2))%vcoord-ver(tr_vertex_no(1))%vcoord)
       tr_vec2=0.5*(ver(tr_vertex_no(3))%vcoord-ver(tr_vertex_no(1))%vcoord)
       End Subroutine Compute_triangle_basis_vectors

!---------------------------------------------------------------------------------------------------------------------------
!      Function that selects a trial function for  Rosenbluth sampling  (algorithm 41, Appendix J, Frenkel and Smit)
!---------------------------------------------------------------------------------------------------------------------------    
    Function Select_from_trial_moves(weight,rosbluth_weight_fac,number_moves) Result(selected_trial)
    Use module_randomnumber ; Use module_datastruct
    Implicit None
    Integer :: number_moves,selected_trial
    Real(Kind=8) ::weight(number_moves),rosbluth_weight_fac,curr_weight,rand_weight
    rand_weight = ran2(Seed)*rosbluth_weight_fac
    selected_trial=1 ; curr_weight=weight(selected_trial)
    Do While(curr_weight .Lt. rand_weight)
    selected_trial=selected_trial+1
    curr_weight = curr_weight+weight(selected_trial)
    Enddo
    Return
    End Function Select_from_trial_moves

!---------------------------------------------------------------------------------------------------------------------------
!	Function  to print the triangle list for a chosen triangle
!---------------------------------------------------------------------------------------------------------------------------    
     Subroutine Print_triangle_vert(tri_no)
     Use module_datastruct
    Implicit None
     Integer :: tri_no
     Print*,str_green,'Vertex list for triangle ',tri_no,' is ',tri(tri_no)%vert
     End Subroutine Print_triangle_vert

!---------------------------------------------------------------------------------------------------------------------------
!     Function to print the image vertex for a chosen vertex
!---------------------------------------------------------------------------------------------------------------------------    
     Subroutine Print_image_vert(ver_no)
     Use module_datastruct
     Implicit None
     Integer :: ver_no
     Print*,str_green,'Image vertex for vertex ',ver_no, ' is ',ver(ver_no)%imver
     End Subroutine Print_image_vert
     
!---------------------------------------------------------------------------------------------------------------------------
!                               Subroutine TO CALCULATE THE ENERGY FOR FLIPPING MOVE
!---------------------------------------------------------------------------------------------------------------------------
    Function flipener(cver) Result(flen)
    USE module_datastruct ; USE module_curvcalc
    IMPLICIT NONE
    Integer,DIMENSION(4) :: cver
    Integer :: i
    Real(Kind=8) :: flen
    flen=0
    four_points: Do i=1,4,1
    flen=flen+((ver(cver(i))%mcur-ver(cver(i))%czero)**2)*kappa*ver(cver(i))%totarea                                                                     !  Helfrich contrib by vertex vert       
    Enddo four_points
    End Function flipener     
!------------------------------------------------------------------------------------------------------------------------
!                               Subroutine TO CALCULATE THE ENERGY FOR VERTEX MOVES
!------------------------------------------------------------------------------------------------------------------------
    Function vermovener(vno) Result(mven)
    USE module_datastruct ; USE module_curvcalc
    IMPLICIT NONE
    Integer :: v,i,vno
    Real(Kind=8) :: mven
    mven=0
    mv_neigh: Do i=1,ver(vno)%nonei,1
    v=ver(vno)%vneipt(i)   
    mven=mven+((ver(v)%mcur-ver(v)%czero)**2)*kappa*ver(v)%totarea                                                                                 !  Helfrich  Energy    
    Enddo mv_neigh 
    mven=mven+((ver(vno)%mcur-ver(vno)%czero)**2)*kappa*ver(vno)%totarea                                                                      
    End Function vermovener
!------------------------------------------------------------------------------------------------------------------------
!                            Subroutine to Implement the Periodic boundary conditions for antigen tip
!------------------------------------------------------------------------------------------------------------------------
    Subroutine PBC_Antigen_tip(ano,indexx,shift)                                                                                     ! ano -antigen number, indexx-component to shift, shift= ±periodic_box_length  
    USE module_datastruct ; USE module_curvcalc
    IMPLICIT NONE
    Integer :: ano,indexx
    Real(Kind=8) :: shift
    antig(ano)%tip_coord(indexx,1) = antig(ano)%tip_coord(indexx,1) + shift
    End Subroutine PBC_Antigen_tip
!------------------------------------------------------------------------------------------------------------------------
!                                Subroutine TO CHECK for the ANGLE BETWEEN THE GIVEN FACES
!------------------------------------------------------------------------------------------------------------------------
    Subroutine faceangchk(trno)
    USE module_datastruct
    IMPLICIT NONE 
    Integer :: trno,j,k                                                                                                              ! trno --> Triangle to be chked for
    Real(KIND=8) :: npro(1,1)                                                                                                        ! Dot product of the Normals       
    Integer,DIMENSION(3) :: ll
    norchk=0 ; ll=0  
    ll= tri(trno)%li
    Do k=1,3
    j=lin(-ll(k))%tr 
    npro=MATMUL(TRANSPOSE(tri(trno)%fnor),tri(j)%fnor)  
    If(npro(1,1).LT.fangle) Then
    norchk=1                                                                                                                         ! norchk --> 0 for the conf to be accepted
    RETURN
    Endif
    Enddo
    RETURN
    End Subroutine faceangchk
!------------------------------------------------------------------------------------------------------------------------
!                 Subroutine TO CHECK FOR THE ANGLE BETWEEN THE GIVEN FACES for triangles with periodic images
!------------------------------------------------------------------------------------------------------------------------
! The idea of having a different subroutine instead of using the ! faceangle check just above is that : 1) the links to triangle
! mapping and vice-versa have not been  implemented for the PB images vertices. This subroutine uses the known vertices of the 
! triangle and compares it to the triangles in the neighboourhood of these vertices. Find the triangles next to each other and 
! perform the face-angle check.

     Subroutine faceangchkPBtr(trno)
     Use module_datastruct
     Implicit None
     Integer:: trno,v(3),i,j,vntr,vn(3),neiflag
     Real(Kind=8):: npro(1,1)
     norchk=0 
     v=tri(trno)%vert                                                                                                               ! Read all the vertices of the given triangle

     Do i=1,3,1                                                                                                                     ! Run over all the 3 vertices making up 'trno'    
     Do j=1,ver(v(i))%nonei-1                                                                                                       ! Explore the neigh triangs of 'trno'            
     neiflag=0
     vntr=ver(v(i))%vneitr(j)                                                                                                       ! vntr is the 'j' th neigh tri of vertex 'v(i)'  
     If(vntr.Ne.trno)Then                                          
       vn=tri(ver(v(i))%vneitr(j))%vert                                                                                             ! Vertices making up the triangle 'vntr'             
       If(v(1).Eq.vn(1) .Or.v(1).Eq.vn(2) .Or. v(1).Eq.vn(3)) Then                                                                  ! Check if the chosen triangle and the triangle in $ 
       neiflag=neiflag+1                                                                                                            ! $ the neig of its 'i'th vertex are neighbours      
       Endif
       If(v(2).Eq.vn(1) .Or.v(2).Eq.vn(2) .Or. v(2).Eq.vn(3)) Then  
       neiflag=neiflag+1
       Endif
       If(v(3).Eq.vn(1) .Or.v(3).Eq.vn(2) .Or. v(3).Eq.vn(3)) Then
       neiflag=neiflag+1
       Endif

       If(neiflag.Eq.2) Then                                                                                                        ! A value of 2 for neiflag is a sign of neig of trno & vntr 
       npro=MATMUL(TRANSPOSE(tri(trno)%fnor),tri(vntr)%fnor)  
       If(npro(1,1).LT.fangle) Then
       norchk=1                                                                                                                     ! norchk --> 0 for the conf to be accepted
       RETURN
       Endif ; Endif ; Endif  
      Enddo ;  Enddo

     End Subroutine faceangchkPBtr
!------------------------------------------------------------------------------------------------------------------
!                 Subroutine TO CALCULATE THE AREA OF THE GIVEN TRIANGLE
!------------------------------------------------------------------------------------------------------------------
     Subroutine onlyarea(tr)
     USE module_datastruct  
     IMPLICIT NONE
     Integer::i,j,k,tr
     Real(KIND=8) ::ax,ay,az,area
     Real(KIND=8),DIMENSION(3,1)::r1,r2,r3,r21,r31

     i=tri(tr)%vert(1) ; j=tri(tr)%vert(2) ; k=tri(tr)%vert(3)                                                                      ! The vertices that make up the triangle tr 
     r1=ver(i)%vcoord ; r2=ver(j)%vcoord ; r3=ver(k)%vcoord                                                                         ! Their corresponding coordinates
     r21=r2-r1 ; r31=r3-r1

     ax=(r21(2,1)*r31(3,1)-r21(3,1)*r31(2,1))*0.5
     ay=(r21(3,1)*r31(1,1)-r21(1,1)*r31(3,1))*0.5
     az=(r21(1,1)*r31(2,1)-r21(2,1)*r31(1,1))*0.5
     area=SQRT(ax**2+ay**2+az**2)

     tri(tr)%ar=area ; ax=ax/area ; ay=ay/area ; az=az/area
     tri(tr)%fnor=RESHAPE((/ax,ay,az/),SHAPE(tri(tr)%fnor))

     tri(tr)%vol=(r1(1,1)*(r2(2,1)*r3(3,1)-r2(3,1)*r3(2,1))+r1(2,1)*(r2(3,1)*r3(1,1)-r2(1,1)*r3(3,1))+ &
       &r1(3,1)*(r2(1,1)*r3(2,1)-r2(2,1)*r3(1,1)))/6.0

     RETURN
     End Subroutine onlyarea
!-------------------------------------------------------------------------------------------------------------------------
!                                   Function TO CHECK THE MAX And MIN BONDLENGTH
!-------------------------------------------------------------------------------------------------------------------------
     Function blcheck(ver1,ver2) result(fres) 
     USE module_datastruct ; USE module_makesurface
     IMPLICIT NONE
     Real(KIND=8)::bl
     Integer::ver1,ver2,fres                                                                                                        !fres=1 ===>false ,0===>true
     fres=0 
     bl=SUM((ver(ver1)%vcoord-ver(ver2)%vcoord)**2)
     If ((bl.GT.blUcut) .OR. (bl.LT.blLcut)) fres=1            
     RETURN
     End Function blcheck
!-------------------------------------------------------------------------------------------------------------------------
!                                   Function TO CHECK THE MINIMUM BONDLENGTH
!-------------------------------------------------------------------------------------------------------------------------
     Function blcheckmin(ver1,ver2) result(fres) 
     USE module_datastruct ; USE module_makesurface
     IMPLICIT NONE
     Real(KIND=8)::bl
     Integer::ver1,ver2,fres                                                                                                        ! fres=1 ===>false ,0===>true
     fres=0 
     bl=SUM((ver(ver1)%vcoord-ver(ver2)%vcoord)**2)
     If(bl.LT.blLcut) fres=1            
     RETURN 
     End Function blcheckmin
!-------------------------------------------------------------------------------------------------------------------------
!                                   Function TO CHECK THE MINIMUM BONDLENGTH
!-------------------------------------------------------------------------------------------------------------------------
     Function Check_minimum_distance(v) result(fres) 
     USE module_datastruct ; USE module_makesurface
     Use, Intrinsic :: ISO_C_BINDING 
     IMPLICIT NONE
     Include '_mod_linkcells_interface.h'

     Real(KIND=8)::bl
     Integer::i,v,fres,mem                                                                                                          ! fres=1 ===>false ,0===>true
     Integer,Pointer :: ring(:),lscl(:),head(:)
     fres=0

     If(v.Le.nver)Then
      Call C_F_Pointer(getcells(1,v,'m'),ring,[1])                                                                                  ! convert the C-pointer returned by neighcells to an one-D fortran pointer ring     
      Call C_F_Pointer(get_lscllist('m'),lscl,[1])                                                                                  ! convert the C-pointer returned by get_lscllist to an one-D fortran pointer lscl   
      Call C_F_Pointer(get_headlist('m'),head,[1])                                                                                  ! convert the C-pointer returned by get_headlist to an one-D fortran pointer head   

      Do i=1,27,1                                                                                                                   ! Perform the self avoidance check over the one ring cells around the cell containing v 
       mem=head(ring(i)+1)                                                                                                          ! get the first membrane particle in cell i 
       Do While(mem .Ne. -1)
        bl=SUM((ver(v)%vcoord-ver(mem+1)%vcoord)**2)                                                                                ! Compute distance if the membrane particle is present 
        If((bl .Lt. blLcut) .And. ((mem+1) .Ne. v))Then                         
        fres=1 ; Return                                                                                                             ! Return if minimum bondlength condition is violated even once 
        Endif
        mem=lscl(mem+1)                                                                                                             ! Move to the next particle 
       Enddo
       Enddo
       Return
     Else                                                                                                                           ! Link Cells have not been implemented for the Phantoms 
       i=1
       Do while(i.Lt.pbnum)                                                                                                         ! Hence we use the traditional method to check over all vertices
       bl=SUM((ver(v)%vcoord-ver(i)%vcoord)**2)
        If((bl .Lt. blLcut) .And. (i .Ne. v) ) Then
          fres=1  ; Return                                                                                                          ! Return fres=1 even if self avoidance is violated once 
        Endif
       i=i+1 
       Enddo
       RETURN 
     Endif
     Return
     End Function Check_minimum_distance

!-------------------------------------------------------------------------------------------------------------------------
!                     Check Self Avoidance of membrane with Nanocarrier
!-------------------------------------------------------------------------------------------------------------------------
      Function Selfavoidance_membrane_Nanocarrier(v) result(fres) 
      USE module_datastruct ; USE module_makesurface
      IMPLICIT NONE
      
      Real(KIND=8) :: bl
      Integer::i, v, fres                                                                                                             ! fres=1 ===>false ,0===>true      
      fres=0
      Do i = 1, num_Nanocarrier, 1                                                                                           
       bl = Sum((ver(v)%vcoord - nc_f(i)%coord)**2)
       If( bl .Lt. sad%vesmem(i))Then                                                   ! A membrane particle cannot be closer than soft radius 
        fres=1 ; Return
       Endif
      Enddo
      Return
      End Function Selfavoidance_membrane_Nanocarrier
!-------------------------------------------------------------------------------------------------------------------------
!                     Check Self Avoidance of Antigen with Nanocarrier
!-------------------------------------------------------------------------------------------------------------------------
      Function Selfavoidance_Antigen_Nanocarrier(antigen_no) result(fres) 
      USE module_datastruct ; USE module_makesurface
      IMPLICIT NONE
      
      Real(KIND=8) :: tip_dist, deltip(3,1)
      Integer::i, antigen_no, fres                                                                     ! fres=1 ===>false ,0===>true      
      fres=0
      Do i=1,num_Nanocarrier,1 
        deltip=antig(antigen_no)%tip_coord - nc_f(i)%coord
        deltip(1,1) = deltip(1,1) - Nint(deltip(1,1)/periodic_box_length)*periodic_box_length
        deltip(2,1) = deltip(2,1) - Nint(deltip(2,1)/periodic_box_length)*periodic_box_length
        tip_dist = Sum(deltip**2)                                    
        If(tip_dist .Lt. sad%antves(antigen_no,i)) Then         
        fres=1 ; Return 
        Endif
      Enddo
      Return
      End Function Selfavoidance_Antigen_Nanocarrier

!-------------------------------------------------------------------------------------------------------------------------
!                     Check Self Avoidance of all antigens on a given vertex (vertex_no) with the Nanocarriers
!-------------------------------------------------------------------------------------------------------------------------
      Function Selfavoidance_VertexAntigens_Nanocarrier(vertex_no) result(fres) 
      USE module_datastruct ; USE module_makesurface
      IMPLICIT NONE
      Real(KIND=8) :: tip_dist, deltip(3,1)
      Integer::i, j, antigen_no, fres, vertex_no
      fres=0
      Do i=1,num_Nanocarrier,1
       Do j=1,ver(vertex_no)%nantigen,1
        antigen_no=ver(vertex_no)%antigen_list(j)
        deltip=antig(antigen_no)%tip_coord - nc_f(i)%coord
        deltip(1,1)=deltip(1,1)-Nint(deltip(1,1)/periodic_box_length)*periodic_box_length
        deltip(2,1)=deltip(2,1)-Nint(deltip(2,1)/periodic_box_length)*periodic_box_length
        tip_dist=Sum(deltip**2)
        If(tip_dist .Lt. sad%antves(antigen_no,i))Then                                                                             ! Membrane-Antigen distance cutoff change from soft radius to radius to test freezing 
        fres=1 ; Return ; Endif
       Enddo
      Enddo
      Return
      End Function Selfavoidance_VertexAntigens_Nanocarrier
!-------------------------------------------------------------------------------------------------------------------------
!                     Check Self Avoidance of all antigens on a given vertex (vertex_no) with the Nanocarriers
!-------------------------------------------------------------------------------------------------------------------------
      Function Selfavoidance_TriangleAntigens_Nanocarrier(triang_no) result(fres) 
      USE module_datastruct ; USE module_makesurface
      IMPLICIT NONE
      
      Real(KIND=8):: tip_dist,deltip(3,1)
      Integer::i, j, antigen_no,fres,triang_no
      fres=0
      Do i=1,num_Nanocarrier,1 
        Do j=1, tri(triang_no)%nantigen, 1
        antigen_no = tri(triang_no)%antigen_list(j)
        deltip = antig(antigen_no)%tip_coord - nc_f(i)%coord
        deltip(1,1)=deltip(1,1) - Nint(deltip(1,1)/periodic_box_length) * periodic_box_length
        deltip(2,1)=deltip(2,1) - Nint(deltip(2,1)/periodic_box_length) * periodic_box_length
        tip_dist=Sum(deltip**2)                         
        If(tip_dist .Lt.  sad%antves(antigen_no,i))Then                                                                             ! Membrane-Antigen distance cutoff change from soft radius to radius to test freezing 
        fres=1 ; Return ; Endif
        Enddo
      Enddo
      Return
      End Function Selfavoidance_TriangleAntigens_Nanocarrier
!-------------------------------------------------------------------------------------------------------------------------
!                     Check Self Avoidance of Antigen with another antigen
!-------------------------------------------------------------------------------------------------------------------------
      Function Selfavoidance_Antigen_Antigen(a1,a2) result(fres) 
      USE module_datastruct ; USE module_makesurface
      IMPLICIT NONE
      Integer::a1,a2,fres                                                                             ! fres=1 ===>false ,0===>true      
      Real(KIND=8):: delbase(3,1),deltip(3,1),base_dist,tip_dist
      
      fres=0
      delbase = antig(a1)%base_coord - antig(a2)%base_coord
      delbase(1,1)=delbase(1,1)-Nint(delbase(1,1)/periodic_box_length)*periodic_box_length
      delbase(2,1)=delbase(2,1)-Nint(delbase(2,1)/periodic_box_length)*periodic_box_length
      base_dist=Sum(delbase**2)                                                                       ! Distance between the base of antigen a1 and a2
      
      If(base_dist .Lt. sad%antant(a1,a2))Then
       fres=1 ; Return                                                                              ! Antigen-Antigen distance cutoff 
      Endif
      
      deltip=antig(a1)%tip_coord-antig(a2)%tip_coord
      deltip(1,1)=deltip(1,1)-Nint(deltip(1,1)/periodic_box_length)*periodic_box_length
      deltip(2,1)=deltip(2,1)-Nint(deltip(2,1)/periodic_box_length)*periodic_box_length
      tip_dist=Sum(deltip**2)                                                                        ! Distance between the base of antigen a1 and a2
      If(tip_dist .Lt. sad%antant(a1,a2)) fres = 1                                                ! Antigen-Antigen distance cutoff 
      Return
      End Function Selfavoidance_Antigen_Antigen

!-------------------------------------------------------------------------------------------------------------------------
!                     Check Self Avoidance of Antigen with another antigen
!-------------------------------------------------------------------------------------------------------------------------
      Function  Selfavoidance_Antigens_In_VertexNeighbourhood(v1) result(fres) 
      USE module_datastruct ; USE module_makesurface
      IMPLICIT NONE
      Integer::v1,v2,fres,i1,i,j,a1,a2,j1,k,t2,a3,nneigh                                                                                            ! fres=1 ===>false ,0===>true      
      Real(KIND=8):: delbase(3,1),deltip(3,1),base_dist,tip_dist
      fres=0

      all_antigens_v1: Do i=1,ver(v1)%nantigen,1                        ! Self avoidance within the antigens on the same vertex
         a1=ver(v1)%antigen_list(i)                                     ! check every pair of antigens that are on vertex v1
 	  If(i .Lt. ver(v1)%nantigen) Then                              ! Avoid multiple calculations for the same pair 
              Do j=i+1,ver(v1)%nantigen
              a2=ver(v1)%antigen_list(j)
	      delbase=antig(a1)%base_coord-antig(a2)%base_coord
	      delbase(1,1)=delbase(1,1)-Nint(delbase(1,1)/periodic_box_length)*periodic_box_length
	      delbase(2,1)=delbase(2,1)-Nint(delbase(2,1)/periodic_box_length)*periodic_box_length
	      base_dist=Sum(delbase**2)                                                       ! Distance between the base of antigen a1 and a2
	      If(base_dist .Lt. sad%antant(a1,a2))Then
	      fres=1 ; Return                                                                 ! Antigen-Antigen distance cutoff 
	      Endif
	      deltip=antig(a1)%tip_coord-antig(a2)%tip_coord
	      deltip(1,1)=deltip(1,1)-Nint(deltip(1,1)/periodic_box_length)*periodic_box_length
	      deltip(2,1)=deltip(2,1)-Nint(deltip(2,1)/periodic_box_length)*periodic_box_length
	      tip_dist=Sum(deltip**2)                                                          ! Distance between the base of antigen a1 and a2
	      If(tip_dist .Lt. sad%antant(a1,a2)) Then
	      fres=1 ; Return
	      Endif
	      Enddo
          Endif

     allneigh_vertices: Do i1=1,ver(v1)%nonei,1                                  ! Self Avoidance over all vertices v2 which are neighbours of vertex v1
      v2=ver(v1)%vneipt(i1)
      all_antigens_v2: Do j=1,ver(v2)%nantigen,1                                 ! with all antigens of v2
	   a2=ver(v2)%antigen_list(j)
	      delbase = antig(a1)%base_coord-antig(a2)%base_coord
	      delbase(1,1) = delbase(1,1)-Nint(delbase(1,1)/periodic_box_length)*periodic_box_length
	      delbase(2,1) = delbase(2,1)-Nint(delbase(2,1)/periodic_box_length)*periodic_box_length
	      base_dist = Sum(delbase**2)                                                            ! Distance between the base of antigen a1 and a2
	      If(base_dist .Lt. sad%antant(a1,a2)) Then
 	      fres=1 ; Return                                                                      ! Antigen-Antigen distance cutoff 
	      Endif
	      deltip=antig(a1)%tip_coord-antig(a2)%tip_coord
	      deltip(1,1)=deltip(1,1)-Nint(deltip(1,1)/periodic_box_length)*periodic_box_length
	      deltip(2,1)=deltip(2,1)-Nint(deltip(2,1)/periodic_box_length)*periodic_box_length
	      tip_dist=Sum(deltip**2)                                                             ! Distance between the base of antigen a1 and a2
	      If(tip_dist .Lt. sad%antant(a1,a2)) Then
	      fres=1 ; Return
	      Endif
        Enddo all_antigens_v2

         If (v2 .Le. nver) Then                                                           ! Boundary vertices  have one triangle less than number of vertices
            nneigh=ver(v2)%nonei
         Else
            nneigh=ver(v2)%nonei-1
         Endif


        Do k=1,nneigh,1
         t2=ver(v2)%vneitr(k)
	 Do j1=1,tri(t2)%nantigen,1
	  a3 = tri(t2)%antigen_list(j1)
	  If(a1 .Ne. a3) Then
	      delbase=antig(a1)%base_coord-antig(a3)%base_coord
	      delbase(1,1)=delbase(1,1)-Nint(delbase(1,1)/periodic_box_length)*periodic_box_length
	      delbase(2,1)=delbase(2,1)-Nint(delbase(2,1)/periodic_box_length)*periodic_box_length
	      base_dist=Sum(delbase**2)                                                                ! Distance between the base of antigen a1 and a2 
	      If(base_dist .Lt. sad%antant(a1,a3))Then
 	      fres=1 ; Return                                                                          ! Antigen-Antigen distance cutoff 
	      Endif
	      deltip=antig(a1)%tip_coord-antig(a3)%tip_coord
	      deltip(1,1)=deltip(1,1)-Nint(deltip(1,1)/periodic_box_length)*periodic_box_length
	      deltip(2,1)=deltip(2,1)-Nint(deltip(2,1)/periodic_box_length)*periodic_box_length
	      tip_dist=Sum(deltip**2)                                                                   ! Distance between the base of antigen a1 and a2
	      If(tip_dist .Lt. sad%antant(a1,a3)) Then
	      fres=1 ; Return
	      Endif
	 Endif
	 Enddo ; Enddo
	
      Enddo allneigh_vertices
      Enddo all_antigens_v1
      End Function  Selfavoidance_Antigens_In_VertexNeighbourhood


!-------------------------------------------------------------------------------------------------------------------------
!                     Check Self Avoidance of Antigen with another antigen
!-------------------------------------------------------------------------------------------------------------------------
      Function Selfavoidance_VertexAntigens_NeighVertAntigens(v1) result(fres) 
      USE module_datastruct ; USE module_makesurface
      IMPLICIT NONE
      Integer::v1,v2,t2,fres,i1,i,j,j1,k,a1,a2,a3,nneigh                                             ! fres=1 ===>false ,0===>true      
      Real(KIND=8):: delbase(3,1),deltip(3,1),base_dist,tip_dist
      fres=0

      Do i=1,ver(v1)%nantigen-1,1                                                                    ! Self avoidance within the antigens on the same vertex
         a1=ver(v1)%antigen_list(i)
         Do j=i+1,ver(v1)%nantigen
           a2=ver(v1)%antigen_list(j)
	      delbase=antig(a1)%base_coord-antig(a2)%base_coord
	      delbase(1,1)=delbase(1,1)-Nint(delbase(1,1)/periodic_box_length)*periodic_box_length
	      delbase(2,1)=delbase(2,1)-Nint(delbase(2,1)/periodic_box_length)*periodic_box_length
	      base_dist=Sum(delbase**2)                                                                     ! Distance between the base of antigen a1 and a2
	      If(base_dist .Lt. sad%antant(a1,a2))Then
 	      fres=1 ; Return                                                                                ! Antigen-Antigen distance cutoff 
	      Endif
	      deltip=antig(a1)%tip_coord-antig(a2)%tip_coord
	      deltip(1,1)=deltip(1,1)-Nint(deltip(1,1)/periodic_box_length)*periodic_box_length
	      deltip(2,1)=deltip(2,1)-Nint(deltip(2,1)/periodic_box_length)*periodic_box_length
	      tip_dist=Sum(deltip**2)                                                                    ! Distance between the base of antigen a1 and a2
	      If(tip_dist .Lt. sad%antant(a1,a2)) Then
	      fres=1 ; Return
	      Endif
       Enddo
     Enddo

     allneigh_vertices: Do i1=1,ver(v1)%nonei,1                                       ! Self Avoidance over all vertices v2 which are neighbours of vertex v1
      v2 = ver(v1)%vneipt(i1)
         Do i=1,ver(v2)%nantigen-1,1                                                   ! Self avoidance within the antigens on the same vertex
          a1 = ver(v2)%antigen_list(i)
          Do j = i+1,ver(v2)%nantigen
           a2 = ver(v2)%antigen_list(j)
           delbase=antig(a1)%base_coord-antig(a2)%base_coord
           delbase(1,1)=delbase(1,1)-Nint(delbase(1,1)/periodic_box_length)*periodic_box_length
           delbase(2,1)=delbase(2,1)-Nint(delbase(2,1)/periodic_box_length)*periodic_box_length
           base_dist=Sum(delbase**2)                                                       ! Distance between the base of antigen a1 and a2
           If(base_dist .Lt. sad%antant(a1,a2))Then
            fres=1 ; Return                                                         ! Antigen-Antigen distance cutoff 
           Endif
           deltip=antig(a1)%tip_coord-antig(a2)%tip_coord
           deltip(1,1)=deltip(1,1)-Nint(deltip(1,1)/periodic_box_length)*periodic_box_length
           deltip(2,1)=deltip(2,1)-Nint(deltip(2,1)/periodic_box_length)*periodic_box_length
           tip_dist=Sum(deltip**2)                                                         ! Distance between the base of antigen a1 and a2
           If(tip_dist .Lt. sad%antant(a1,a2)) Then
             fres=1 ; Return
           Endif
         Enddo
        Enddo

      Do i=1,ver(v1)%nantigen,1                                                                                                     ! over all antigens of v1 
        a1 = ver(v1)%antigen_list(i)
        Do j = 1,ver(v2)%nantigen,1                                                                                                   ! with all antigens of v2
         a2 = ver(v2)%antigen_list(j)
         delbase=antig(a1)%base_coord-antig(a2)%base_coord
         delbase(1,1)=delbase(1,1)-Nint(delbase(1,1)/periodic_box_length)*periodic_box_length
         delbase(2,1)=delbase(2,1)-Nint(delbase(2,1)/periodic_box_length)*periodic_box_length
         base_dist=Sum(delbase**2)                                                                                             ! Distance between the base of antigen a1 and a2
         If(base_dist .Lt. sad%antant(a1,a2)) Then
           fres=1 ; Return                                                           ! Antigen-Antigen distance cutoff 
         Endif
         deltip=antig(a1)%tip_coord-antig(a2)%tip_coord
         deltip(1,1)=deltip(1,1)-Nint(deltip(1,1)/periodic_box_length)*periodic_box_length
         deltip(2,1)=deltip(2,1)-Nint(deltip(2,1)/periodic_box_length)*periodic_box_length
         tip_dist=Sum(deltip**2)                                                                                                ! Distance between the base of antigen a1 and a2
         If(tip_dist .Lt. sad%antant(a1,a2)) Then
          fres=1 ; Return
         Endif
        Enddo 

        If (v2 .Le. nver) Then                                                                                                    ! Boundary vertices  have one triangle less than number of vertices
           nneigh=ver(v2)%nonei
        Else
           nneigh=ver(v2)%nonei-1
        Endif

      Do k=1,nneigh,1
       t2=ver(v2)%vneitr(k)
        Do j1=1,tri(t2)%nantigen,1
          a3 = tri(t2)%antigen_list(j1)
          If(a1 .Ne. a3) Then
            delbase=antig(a1)%base_coord-antig(a3)%base_coord
            delbase(1,1)=delbase(1,1)-Nint(delbase(1,1)/periodic_box_length)*periodic_box_length
            delbase(2,1)=delbase(2,1)-Nint(delbase(2,1)/periodic_box_length)*periodic_box_length
            base_dist=Sum(delbase**2)                                                          ! Distance between the base of antigen a1 and a2
            If(base_dist .Lt. sad%antant(a1,a3)) Then
             fres=1 ; Return                                                                   ! Antigen-Antigen distance cutoff 
            Endif
            deltip=antig(a1)%tip_coord-antig(a3)%tip_coord
            deltip(1,1)=deltip(1,1)-Nint(deltip(1,1)/periodic_box_length)*periodic_box_length
            deltip(2,1)=deltip(2,1)-Nint(deltip(2,1)/periodic_box_length)*periodic_box_length
            tip_dist=Sum(deltip**2)                                                                                                ! Distance between the base of antigen a1 and a2
            If(tip_dist .Lt. sad%antant(a1,a3)) Then
             fres=1 ; Return
            Endif
          Endif
        Enddo 
      Enddo

      Enddo
      Enddo allneigh_vertices
      End Function Selfavoidance_VertexAntigens_NeighVertAntigens
           
!-------------------------------------------------------------------------------------------------------------------------
!                     Check Self Avoidance of Antigen with another antigen
!-------------------------------------------------------------------------------------------------------------------------
      Function Selfavoidance_Allantigens() result(fres) 
      USE module_datastruct ; USE module_makesurface
      IMPLICIT NONE
      Integer::a1,a2,fres                                                                            ! fres=1 ===>false ,0===>true      
      Real(KIND=8):: delbase(3,1),base_dist
      fres=0
      Do a1=1,num_antigens-1
       Do a2=a1+1,num_antigens,1
        delbase=antig(a1)%base_coord-antig(a2)%base_coord
        delbase(1,1)=delbase(1,1) - Nint(delbase(1,1)/periodic_box_length)*periodic_box_length
        delbase(2,1)=delbase(2,1) - Nint(delbase(2,1)/periodic_box_length)*periodic_box_length
        base_dist=Sum(delbase**2)                                                                     ! Distance between the base of antigen a1 and a2
        If(base_dist .Lt. sad%antant(a1,a2)) Then
          fres=1 ; Return                                                                            ! Antigen-Antigen distance cutoff 
        Endif
       Enddo 
      Enddo
      End Function Selfavoidance_Allantigens
!-------------------------------------------------------------------------------------------------------------------------
!                     Check Self Avoidance of Antigen with another antigen
!-------------------------------------------------------------------------------------------------------------------------
      Function Selfavoidance_Allantigens1() result(fres) 
      USE module_datastruct ; USE module_makesurface
      IMPLICIT NONE
      Integer::a1,a2,fres                                                                                                           ! fres=1 ===>false ,0===>true      
      Real(KIND=8):: delbase(3,1),base_dist
      fres=0
      Do a1=1,num_antigens-1
        Do a2=a1+1,num_antigens,1
         delbase=antig(a1)%base_coord-antig(a2)%base_coord
         delbase(1,1)=delbase(1,1)-Nint(delbase(1,1)/periodic_box_length)*periodic_box_length
         delbase(2,1)=delbase(2,1)-Nint(delbase(2,1)/periodic_box_length)*periodic_box_length
         base_dist=Sum(delbase**2)                                                                                                     ! Distance between the base of antigen a1 and a2
        If(base_dist .Lt. sad%antant(a1,a2))Then
          Print*,'Self avoidance violated for antigens ',a1,' and ',a2
          Print*,str_red
          Print*,'Antigen 1 details ',antig(a1)%diffus_tri,antig(a1)%vertex
          Print*,'Antigen 2 details ',antig(a2)%diffus_tri,antig(a2)%vertex
          Print*,'Antigen 1 coords ',antig(a1)%base_coord
          Print*,'Antigen 2 coords ',antig(a2)%base_coord
          Print*,'base_dist ',base_dist
          Print*,str_blue
          Print*,'Antigen 1 old coords ',mp%antig(a1)%base_coord
          Print*,'Antigen 2 old coords ',mp%antig(a2)%base_coord
          Print*,str_black
          fres=1 ; Return                                                             ! Antigen-Antigen distance cutoff 
        Endif
       Enddo 
      Enddo
      Return
      End Function Selfavoidance_Allantigens1 
!-------------------------------------------------------------------------------------------------------------------------
!                                   Function to check selfavoidance of one antigen with all
!-------------------------------------------------------------------------------------------------------------------------
      Function Antigen_selfavoidance(a1) result(fres) 
      USE module_datastruct ; USE module_makesurface
      Use, Intrinsic :: ISO_C_BINDING 
      IMPLICIT NONE
      Include '_mod_linkcells_interface.h'

      Real(KIND=8) :: base_dist, tip_dist, delbase(3,1), deltip(3,1)
      Integer :: i, fres, mem, a1                                       ! fres=1 ===>false ,0===>true
      Integer,Pointer :: ring(:), lscl(:), head(:)
      fres=0

       Call C_F_Pointer(getcells(1,a1,'c'),ring,[1])                ! convert the C-pointer returned by neighcells to an one-D fortran pointer ring     
       Call C_F_Pointer(get_lscllist('c'),lscl,[1])                 ! convert the C-pointer returned by get_lscllist to an one-D fortran pointer lscl   
       Call C_F_Pointer(get_headlist('c'),head,[1])                 ! convert the C-pointer returned by get_headlist to an one-D fortran pointer head   

       Do i=1,27,1                                                  ! Perform the self avoidance check over the one ring cells around the cell containing a1
        mem=head(ring(i)+1)                                         ! get the first antigen particles in cell mem ( ring contain cells from 0 to N_cells)
        Do While(mem .Ne. -1)
         If((mem+1) .Ne. a1) Then                                   ! head object + 1 to make it accessible to fortran
            delbase=antig(a1)%base_coord-antig(mem+1)%base_coord
	    delbase(1,1)=delbase(1,1)-Nint(delbase(1,1)/periodic_box_length)*periodic_box_length   ! Periodic boundary conditions applied to check the distance 
	    delbase(2,1)=delbase(2,1)-Nint(delbase(2,1)/periodic_box_length)*periodic_box_length
	    base_dist=Sum(delbase**2)                                                              ! Distance between the base of antigen a1 and a2
	    If(base_dist .Lt. sad%antant(a1,mem+1))Then
	    fres=1 ; Return 
	    Endif                                                                ! Antigen-Antigen distance cutoff 
            deltip=antig(a1)%tip_coord-antig(mem+1)%tip_coord
	    deltip(1,1)=deltip(1,1)-Nint(deltip(1,1)/periodic_box_length)*periodic_box_length
	    deltip(2,1)=deltip(2,1)-Nint(deltip(2,1)/periodic_box_length)*periodic_box_length
	    tip_dist=Sum(deltip**2)                                                               ! Distance between the base of antigen a1 and a2
	    If(tip_dist .Lt. sad%antant(a1,mem+1)) Then
	    fres=1 ; Return ; Endif                                                               ! Antigen-Antigen distance cutoff 
        Endif
         mem=lscl(mem+1)                                                                            ! Move to the next antigen 
        Enddo
        Enddo
        Return
      End Function Antigen_selfavoidance
!-------------------------------------------------------------------------------------------------------------------------
!                                   Function to check selfavoidance of one antigen with all
!-------------------------------------------------------------------------------------------------------------------------
      Subroutine Print_linkcelldata(ch,ant_no)
      Use, Intrinsic :: ISO_C_BINDING 
      IMPLICIT NONE
      Include '_mod_linkcells_interface.h'
      Integer:: ant_no,i,mem
      Character :: ch 
      Integer,Pointer :: ring(:),lscl(:),head(:)
      
      Call C_F_Pointer(getcells(1,ant_no,ch),ring,[1])                                                                              ! convert the C-pointer returned by neighcells to an one-D fortran pointer ring     
      Call C_F_Pointer(get_lscllist(ch),lscl,[1])                                                                                   ! convert the C-pointer returned by get_lscllist to an one-D fortran pointer lscl   
      Call C_F_Pointer(get_headlist(ch),head,[1])                                                                                   ! convert the C-pointer returned by get_headlist to an one-D fortran pointer head   

        Print*,'The antigen number is checked for is ',ant_no
        Do i=1,27,1                                                                                                                 ! Perform the self avoidance check over the one ring cells around the cell containing a1
        mem=head(ring(i)+1)                                                                                                         ! get the first antigen particles in cell mem
        Do While(mem .Ne. -1)
        if ((mem+1) .Eq. ant_no) Then
        Print*, 'mem+1 and ant_no ',i,mem+1,ant_no
        Print*, 'antigen number found in cell ', ring(i)-1
        Print*
        Endif
        mem=lscl(mem+1)                                                                                                             ! Move to the next antigen 
        Enddo ;  Enddo 
      End Subroutine Print_linkcelldata      
!-------------------------------------------------------------------------------------------------------------------------
!                                   Function to check selfavoidance of one antigen with all
!-------------------------------------------------------------------------------------------------------------------------
      Function return_linkcelldata(ch,ant_no) Result(fres)
      Use, Intrinsic :: ISO_C_BINDING 
      IMPLICIT NONE
      Include '_mod_linkcells_interface.h'
      Integer:: ant_no,i,mem,fres
      Character :: ch 
      Integer,Pointer :: ring(:),lscl(:),head(:)
      fres=-1
      
      Call C_F_Pointer(getcells(1,ant_no,ch),ring,[1])                                                                              ! convert the C-pointer returned by neighcells to an one-D fortran pointer ring     
      Call C_F_Pointer(get_lscllist(ch),lscl,[1])                                                                                   ! convert the C-pointer returned by get_lscllist to an one-D fortran pointer lscl   
      Call C_F_Pointer(get_headlist(ch),head,[1])                                                                                   ! convert the C-pointer returned by get_headlist to an one-D fortran pointer head   

        Do i=1,27,1                                                                                                                 ! Perform the self avoidance check over the one ring cells around the cell containing a1
        mem=head(ring(i)+1)                                                                                                         ! get the first antigen particles in cell mem
        Do While(mem .Ne. -1)
        if ((mem+1) .Eq. ant_no) Then
        fres=ring(i)-1
        Endif
        mem=lscl(mem+1)                                                                                                             ! Move to the next antigen 
        Enddo ;  Enddo 
        Return
      End function return_linkcelldata 
!-------------------------------------------------------------------------------------------------------------------------
!                                   Function to check if an object is present in a cell or not
!-------------------------------------------------------------------------------------------------------------------------
      Function check_ifin_linkcell(obj_no,cellno,ch) Result(fres)
      Use, Intrinsic :: ISO_C_BINDING 
      IMPLICIT NONE
      Include '_mod_linkcells_interface.h'
      Integer:: cellno,obj_no,mem,fres
      Character :: ch 
      Integer,Pointer :: lscl(:),head(:)
      fres=0
      
      Call C_F_Pointer(get_lscllist(ch),lscl,[1])                                                                                   ! convert the C-pointer returned by get_lscllist to an one-D fortran pointer lscl   
      Call C_F_Pointer(get_headlist(ch),head,[1])                                                                                   ! convert the C-pointer returned by get_headlist to an one-D fortran pointer head   

      mem=head(cellno+1)
      Do while(mem .Ne. -1) 
      If ((mem+1) .Eq. obj_no) Then                                                                                                 ! check if an object of type ch is present in cell cellno  
      fres=1                                                                                                                        ! If present return 1 else return 0
      Return
      Endif
      mem=lscl(mem+1)
      Enddo
      End function check_ifin_linkcell


!--------------------------------------------------------------------------------------------------------------------------
!                                 Function to compute the change in biasing potential when the membrane is moved
!--------------------------------------------------------------------------------------------------------------------------
      Function meanr_biasing_potential(ncnum,oldmeanr,currmeanr) Result(bias_ener)
      Use module_datastruct
      Implicit None
      Include '_mod_getcoordinates_interface.h'
      Real(Kind=8) :: oldmeanr(3,1),currmeanr(3,1),bias_ener
      Real(Kind=8) :: refdisold, refdisnew
      Integer :: ncnum
      bias_ener=0
      refdisold = Sqrt(Sum((nc_f(ncnum)%coord - oldmeanr)**2)) - nc_f(ncnum)%biasref
      refdisnew = Sqrt(Sum((nc_f(ncnum)%coord - currmeanr)**2)) - nc_f(ncnum)%biasref
      bias_ener = 0.5*nc_f(ncnum)%kbias*(refdisnew**2 - refdisold**2)                      ! change in the biaising energy
      Return
      End Function meanr_biasing_potential
!--------------------------------------------------------------------------------------------------------------------------
!                                 Function to compute the change in biasing potential when the membrane is moved
!--------------------------------------------------------------------------------------------------------------------------
      Function meanH_biasing_potential(ncnum, oldmeanH, currmeanH) Result(bias_ener)
      Use module_datastruct
      Implicit None
      Include '_mod_getcoordinates_interface.h'
      Real(Kind=8) :: oldmeanH, currmeanH, bias_ener, refold, refnew
      Integer :: ncnum
      refold = oldmeanH - nc_f(ncnum)%biasref
      refnew = currmeanH - nc_f(ncnum)%biasref
      bias_ener = 0.5*nc_f(ncnum)%kbias*(refnew**2-refold**2)                                                                       ! sum of all bias_ener   ( kh =0.5*k_curv)
      Return
      End Function meanH_biasing_potential

      End MODULE module_mcsmoves
