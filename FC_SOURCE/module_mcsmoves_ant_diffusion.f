!===========================================================================================================================
!===========================================================================================================================
!===========================================================================================================================
!      $$ modflip                          MODULE TO FLIP THE BONDS And DISPLACE THE CHOSEN VERTEX POINTS
!===========================================================================================================================
!===========================================================================================================================

       MODULE module_mcsmoves
       Integer ::norchk,amov,aflip,tmov,tflip
       Real(KIND=8)::flen,mven,dmven,spien,prob,rnum,spexen                                                                         ! flen,mven,nemven are temp var for energy cal 
       Real(KIND=8)::dfang
       contains

       Subroutine Membrane_MonteCarlo(mcsstep)
       Use Module_datastruct ; Use Module_randomnumber
       Use Module_makesurface, ONLY: antigenzero
       Implicit None
       Include '_mod_mcsmoves_interface.h'
       Integer :: i,lrand,vrand,arand,mvinter,mcsstep
       mvinter=NINT(Real(tlink)/nver)
       mcsmovecounter=mcsmovecounter+1

       Do i=1,num_antigens,1  
       arand=Nint(ran2(seed)*num_antigens)+1
       If(arand.Le.num_antigens) Call Antigen_diffusion_on_triangle(arand)
       Enddo

       If(mod(mcsmovecounter,datainterval).Eq.0)Then
       Open(01,File='./Membrane-area.dat',Position='Append')
       Write(01,*),mcsmovecounter/datainterval,Sum(ver(1:nver)%totarea)
       Endif
       End Subroutine Membrane_MonteCarlo

!---------------------------------------------------------------------------------------------------------------------------
!                                           Subroutine TO FLIP A LINK
!--------------------------------------------------------------------------------------------------------------------------- 

      Subroutine flipping(rand)
      USE module_datastruct ; USE module_curvcalc 
      USE module_makesurface ; Use module_randomnumber
      IMPLICIT NONE                                                                                                                 ! tdv and tdt are the local data structure
      Include '_mod_mcsmoves_interface.h'
       
      Real(KIND=8) :: inen,fien,delRE,delE                                                                                            ! blcheck,blcheckmin -->functions to be called  
      Integer :: vt1,vt2,i,j,ep1,ep2,fv1,fv2,blch,old_tri_no,new_tri_no                                                              ! vt -->vert of triangle; norchk -->nor check  
      Integer :: fvp1,fvm1,fvp2,fvm2,tn1,tn2,conch,rand,mchk,ant_chk_flag       
      Integer :: lp1,lp2,ch1,trn,llist(6),antigen_ind,an_ver(4,2),curcell,old_curcell
      Integer,DIMENSION(4) :: cver                                                                                                    ! Temproary variables for case statement        
      Integer,DIMENSION(3) :: t1,t2,tmp,tmp1                                                                                          ! ep -->endpoint of link                        

             
       t1=tri(lin(rand)%tr)%vert                                                                                                    ! To read the vertices linked to the triangles 
       t2=tri(lin(-rand)%tr)%vert 
       tn1=lin(rand)%tr ; tn2=lin(-rand)%tr                                                                                         ! tn1 and tn2 are the names of the triangles    
       lp1=lin(rand)%sep(1) ; lp2=lin(rand)%sep(2)
       ep1=0;ep2=0 ; blch=0 ; an_ver=0 ; ant_chk_flag=0

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

        mp%triangle(tn1)=tri(tn1) ; mp%triangle(tn2)=tri(tn2)
        ForAll(i=1:6) mp%link(llist(i))=lin(llist(i))
        ForAll(i=1:4) mp%vertex(cver(i))=ver(cver(i))                                                                               ! Original state of involved vertices
      
        flen=0 ; call flipener(cver)                                                                                                ! Call for initial energy calculation
        inen=flen 

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
        
        
       flen=0 ; mchk=0
                                                      
       Call faceangchk(tn1)                                                                                                         ! Check  angle between tri 1 and neigh faces       
       fang_tr1:If(norchk .EQ.0) Then                                                                                               ! Proceed further only if Yes (Max =150 \degrees)  
       Call faceangchk(tn2)                                                                                                         ! Check for angle between tri2 and its neigh faces 
       fang_tr2 : If(norchk .EQ.0) Then                             

        Do i=1,4,1
        Call normalcalc(cver(i))                                                                                                    ! Update the curvature at the chosen vertices 
          If(ver(cver(i))%antigen_flag .Eq.1 .And. mchk.Eq.0)Then
          an_ver(i,1)=1 ; an_ver(i,2)=ver(cver(i))%antigen_index                                                                    ! Store the antigen flag and antigen index for use in restoration
          mp%antig(an_ver(i,2))=antig(an_ver(i,2))			                                                            ! Store the antigen tip position
	  old_tri_no=antig(an_ver(i,2))%diffus_tri                                                                                  ! Choose the triangle the antigen is associated to before the flip move   
	  new_tri_no=old_tri_no                                                                                                     ! If the triangle has the same vertex after flip the same association is maintained 
	  If (i .Eq.3 .And. old_tri_no .Eq. tn2) new_tri_no=tn1			          			  	            ! cver(3) and cver(4) will lose one triangle each during flip moves
	  If (i .Eq.4 .And. old_tri_no .Eq. tn1) new_tri_no=tn2								  	    ! If the antigen is in the lost triangle move it to the existing one  
          If(Displace_Antigen_on_Triangle(cver(i),new_tri_no,an_ver(i,2)).Eq.1) mchk=1                                              ! displaces the antigen and checks for self-avoidance with other antigens
	  If((Selfavoidance_Antigen_Nanocarrier(an_ver(i,2)).Eq.1) .And. (mchk.Eq.0)) mchk=1                                        ! check for self-avoidance of the antigen with the Nanocarriers 
          If(mchk.Eq.0) Call Rearrange_Antigen_on_Triangles(old_tri_no,new_tri_no,an_ver(i,2))                                      ! Book keeping step for antigen list in triangles and triangle list in antigens 

          If(bondedstate(an_ver(i,2)-1) .Eqv. .TRUE. .And. mchk.Eq.0) delRE=delRE +&
          antigen_reaction_E_change(an_ver(i,2)-1,antig(an_ver(i,2))%tip_coord(1,1),mp%antig(an_ver(i,2))%tip_coord(1,1))             ! Compute the reaction energy if the antigen is already bonded to an antibody
          Endif
        Enddo

	If(mchk .Eq.0)Then
        call flipener(cver)                                                                                                         ! Call for final energy calculation             
        fien=flen                                                                                                           
        delE=delRE+fien-inen                                                                                                        ! the energy change associated with the move
        If((delE.GT.0) .And. (exp(-beta*delE).Lt.ran2(Seed))) mchk=1                                                                ! Metropolis Algorithm
	Endif

        Else
        mchk=1 
        Endif fang_tr2
        Else
        mchk=1 
        Endif fang_tr1

        flip_failed:If((mchk .Eq. 1))then
        tri(tn1)=mp%triangle(tn1) ; tri(tn2)=mp%triangle(tn2)
        ForAll (i=1:6) lin(llist(i))=mp%link(llist(i))
        Do i=1,4,1
        ver(cver(i))=mp%vertex(cver(i))
        If(an_ver(i,1).Eq.1) antig(an_ver(i,2))=mp%antig(an_ver(i,2))                                                               ! Restore the antigen tip position if move is rejected
        Enddo 
        Endif flip_failed

	flip_accp: If(mchk==0)Then                                
        Do i=1,4,1
 	 If((an_ver(i,1).Eq.1))Then
         Call update_linkcells(an_ver(i,2),&
         mp%antig(an_ver(i,2))%base_coord(1,1),mp%antig(an_ver(i,2))%base_coord(2,1),mp%antig(an_ver(i,2))%base_coord(3,1),'c')     ! Update the linkcells for the antigen at vertex cver(i)  
         If (bondedstate(an_ver(i,2)-1) .Eqv. .TRUE.) Call update_antigenbondedstate(an_ver(i,2)-1)                                 ! Check if the bonded antigen tip is within a cutoff
         Endif
 	 Do j=1,ver(cver(i))%pbimno,1                                                                                               ! Update the surface quantifiers of the corresponding periodic images
        call mapimagetovertex(ver(cver(i))%pbmap(j))
        Enddo
        Enddo
        Endif flip_accp 

      End Subroutine flipping 
!---------------------------------------------------------------------------------------------------------------------------
!                  $vermov           Subroutine TO MOVE THE VERTEX
!---------------------------------------------------------------------------------------------------------------------------
! The move subroutine has been written in such a waay that whenever a  boundary node is displaced its corrresponding periodic image 
! is also given the same displacement. While area of the periodic triangles are calculated by normal means the surface quantifiers 
! of the perriodic vertex image are got from its original vertex

      Subroutine movevertex(vert)
      USE module_datastruct; USE module_curvcalc 
      USE module_makesurface ; USE module_writedata
      Use module_randomnumber
      Use, Intrinsic :: ISO_C_BINDING 
      IMPLICIT NONE 
      Include '_mod_mcsmoves_interface.h'
      Integer::vert,i,ip1,j,k,trno,blch,imgv,opbmap(3),mchk,v1                                                                      ! i,j,k,trno -> temp var, *ch-> check vals
      Integer :: antigen_ind,totantigen,an_ver(11,4)
      Real(KIND=8)::inen,fien,delE,delRE,prob,len                                                                                   ! init and fin ener, probability          
      Real(KIND=8),DIMENSION(3,1)::dr                                                                                               

      inen=0; fien=0 ; an_ver=0 

      If(ver(vert)%boundary.Eq.0)Then 
      dr=1.0-2.0*reshape((/ran2(Seed),ran2(Seed),ran2(Seed)/),(/3,1/))                                                              ! A Small displacement vector
      Else
      dr=reshape((/zero,zero,1.-2.*ran2(Seed)/),(/3,1/))                                                                            ! The boundary vertex moves only along z
      Endif

      dr=dr*dispmag;
      mp%vertex(vert)=ver(vert)                                                                                                     ! Store the data for the chosen vertex

       If(SUM(dr**2).GT.maxdisp) RETURN                                                                                             ! Max allowed displacement is (0.0025*blen)  
       ver(vert)%vcoord=ver(vert)%vcoord+dr                                                                                         ! New displaced position                     

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

       If(Selfavoidance_membrane_Nanocarrier(vert).NE.0) Then                                                                       ! Check Self avoidance of vertex with Nanocarrier 
       ver(vert)%vcoord=mp%vertex(vert)%vcoord
       RETURN 
       Endif
                                                                                                                                    ! do the same min and max check for the images 
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
       ver(opbmap(j))=mp%vertex(opbmap(j))
       Enddo
       ver(vert)%vcoord=mp%vertex(vert)%vcoord                                                                                      ! Replace the original coordinates of the chosen vertex 
       RETURN 
       Endif                                                                                                                        ! with all neighbouring vertices.
       Enddo nei_len_imgchk

       If(Check_minimum_distance(imgv) .NE.0) Then                                                                                  ! Check for minimum blength const with all verts
       Do j=1,i,1                                                                                                                   ! Restore the coordinates of all vertices that have been checked 
       ver(opbmap(j))=mp%vertex(opbmap(j))
       Enddo
       ver(vert)%vcoord=mp%vertex(vert)%vcoord                                                                                      ! Replace the original coordinates of the chosen vertex 
       RETURN  
       Endif                                                                                                                        !to proceed further else replace the original val

       Do k=1,ver(imgv)%nonei
       mp%vertex(ver(imgv)%vneipt(k))=ver(ver(imgv)%vneipt(k))
       If(k.Lt.ver(imgv)%nonei)Then 
       mp%triangle(ver(imgv)%vneitr(k))=tri(ver(imgv)%vneitr(k)) 
       Endif
       Enddo

       Do k=1,ver(imgv)%nonei-1                                                                                                     ! 7 Mar 12  | Areacalc(j) has been moved to seperate loop   
       Call areacalc(ver(imgv)%vneitr(k))                                                                                           ! Area of all triangles surrouding an image is updated      
       Enddo
       Enddo PB_ver
       Endif vert_PBimage

       vert_antigen:If((ver(vert)%antigen_flag .Eq.1)) Then             
       an_ver(1,1)=1 ; an_ver(1,2)=ver(vert)%antigen_index  ; an_ver(1,3)=vert                                                      ! Collect the antigen flag and antigen vertices 
       mp%antig(an_ver(1,2))=antig(an_ver(1,2))
       Endif  vert_antigen

       mven=0 ; Call vermovener(vert)                                                                                               ! Call for the initial energy calculation 
       inen=mven 

       store_ntr_area:  Do i=1,ver(vert)%nonei                                                                                      ! Update the area of all neighbouring triangles  
       j=ver(vert)%vneitr(i) ; k=ver(vert)%vneipt(i)
       mp%vertex(k)=ver(k) ;  mp%triangle(j)=tri(j)
       If(k .Gt.nver) k=ver(k)%imver
       nei_antigen: If((ver(k)%antigen_flag .Eq.1))Then                 
       an_ver(i+1,1)=1 ; an_ver(i+1,2)=ver(k)%antigen_index  ;an_ver(i+1,3)=k                                                       ! Collect the antigen flag and antigen vertices 
       mp%antig(an_ver(i+1,2))=antig(an_ver(i+1,2))                                                                                 ! Store the antigen to  mp%antig (i+1 -> neighbour_number+1)
       Else
       an_ver(i+1,1)=0
       Endif nei_antigen
       Enddo store_ntr_area 

       upd_ntr_area:  Do i=1,ver(vert)%nonei                                                                                        ! Update the area of all neighbouring triangles  
       call areacalc(ver(vert)%vneitr(i))                                                                                           ! Calculate the new area of each face            
       Enddo upd_ntr_area
       norchk=0 ; mchk=0 ; trno=1

       Do While(trno.Le.ver(vert)%nonei .And. norchk.Eq.0)
       If(tri(ver(vert)%vneitr(trno))%pbflag.Eq.0) Then                                                                             ! For submarginal trinagles call faceangchk()  
       call faceangchk(ver(vert)%vneitr(trno)) 
       Else
       Call faceangchkPBtr(ver(vert)%vneitr(trno))                                                                                  ! For the marginal triangle call faceangchkPBtr()
       Endif
       If(norchk.Ne.0) Then
       mchk=1   ; Exit  ; Endif                                                                                                     ! set mchk=0 if face angle check fails  
       trno=trno+1
       Enddo

       vert_update_antigen:If(mchk.Eq.0) Then
       Call normalcalc(vert)                                                                                                        ! mchk=0 is the condition for proceeding with further calculations  
       If(an_ver(1,1).Eq.1)Then
       If(Displace_Antigen_on_Triangle(an_ver(1,3),antig(an_ver(1,2))%diffus_tri,an_ver(1,2)).Eq.1) mchk=1                          ! Triangle the antigen is currently sitting on 	
       If((mchk.Eq.0) .And. (Selfavoidance_Antigen_Nanocarrier(an_ver(1,2)).Eq.1)) mchk=1                                           ! Set flag to 1 with the antigen at vert violates self avoidance  
       Endif
       Endif vert_update_antigen

       i=1
       vertnei_update_antigen: Do While((mchk.Eq.0) .And. (i.Le.ver(vert)%nonei))
       j=ver(vert)%vneipt(i) 
       If(j.Gt. nver)Then                                                                                                           ! Call only those vertices which are not periodic images 
       Call normalcalc(ver(j)%imver)                                                                                                ! For periodic images, update their original vertex 
       Call mapimagetovertex(j)                                                                                                     ! map the original vertex props to the image        
       Else 
       Call normalcalc(j)     
       Endif

       If(an_ver(i+1,1).Eq.1)Then
       If(Displace_Antigen_on_Triangle(an_ver(i+1,3),antig(an_ver(i+1,2))%diffus_tri,an_ver(i+1,2)).Eq.1) mchk=1                    ! Triangle the antigen is currently sitting on 	
       If((mchk.Eq.0) .And. (Selfavoidance_Antigen_Nanocarrier(an_ver(1,2)).Eq.1)) mchk=1                                           ! Set flag to 1 with the antigen at vert violates self avoidance  
       If ((mchk.Eq.0) .And. (Selfavoidance_Antigen_Nanocarrier(an_ver(i+1,2)).Eq.1)) mchk=1
       Endif  
       i=i+1
       Enddo vertnei_update_antigen

      self_avoidance_satisfied: If(mchk .Eq.0)Then
       mven=0 ;  call vermovener(vert)                                                                                              ! Call for final energy calculation
       fien=mven    

       delRE=0.0
       Do i=1,ver(vert)%nonei+1,1
       If ((an_ver(i,1).Eq.1) .And. (bondedstate(an_ver(i,2)-1).Eqv. .TRUE.))Then
       delRE = delRE+antigen_reaction_E_change(an_ver(i,2)-1,antig(an_ver(i,2))%tip_coord(1,1),mp%antig(an_ver(i,2))%tip_coord(1,1))! Compute the change in the reaction energy  (-1 takes care of indices in C++)
       Endif
       Enddo

       delE=fien+delRE-inen                                                                                                         ! Total change in energy of the system due to vertex move                     

       move_metropolis: If(exp(-beta*delE) .Ge. ran2(Seed)) Then                                                                    ! Metropolis Scheme to reject a move                                          
        Call Update_linkcells(vert,mp%vertex(vert)%vcoord(1,1),mp%vertex(vert)%vcoord(2,1),mp%vertex(vert)%vcoord(3,1),'m')         ! Update the link cells after the vertex move is accepted
        If(an_ver(1,1) .Eq. 1) Then
        Call Update_linkcells(an_ver(1,2),&
        mp%antig(an_ver(1,2))%base_coord(1,1),mp%antig(an_ver(1,2))%base_coord(2,1),mp%antig(an_ver(1,2))%base_coord(3,1),'c')      ! Update the link cell for antigens if the moved vertex has an antigen  
	If (bondedstate(an_ver(1,2)-1).Eqv. .TRUE.) Call update_antigenbondedstate(an_ver(1,2)-1)                                   ! Check if the antigen has moved beyond the maximum stretch cutoff
	Endif

        Do i=1,ver(vert)%nonei,1
        If ((an_ver(i+1,1).Eq.1) )Then
        Call Update_linkcells(an_ver(i+1,2),&
        mp%antig(an_ver(i+1,2))%base_coord(1,1),mp%antig(an_ver(i+1,2))%base_coord(2,1),mp%antig(an_ver(i+1,2))%base_coord(3,1),'c')! Update the link cell for antigens if the moved vertex has an antigen  
    	If (bondedstate(an_ver(i+1,2)-1).Eqv. .TRUE.) Call update_antigenbondedstate(an_ver(i+1,2)-1)       		            ! Check if the antigen has moved beyond the maximum stretch cutoff
	Endif
        k=ver(vert)%vneipt(i)
        Do j=1,ver(k)%pbimno,1                                                                                                      ! Images of boundary vertices are updated
        call mapimagetovertex(ver(k)%pbmap(j))
        Enddo ; Enddo
        Do j=1,ver(vert)%pbimno,1                                                                                                    ! Images of vert is updated 
        call mapimagetovertex(ver(vert)%pbmap(j))
        Enddo
       Else
       mchk=1                                                                                                                        ! Rejected from Metropolis 
       Endif move_metropolis
       Endif self_avoidance_satisfied

       If (mchk .Eq.1) Then
       ver(vert)=mp%vertex(vert)
       If(an_ver(1,1).Eq.1) antig(an_ver(1,2))=mp%antig(an_ver(1,2))                                                                ! Restore the antigen base and tip position if move is rejected
       rest_neig: Do i=1,ver(vert)%nonei,1                                                                                          ! Update all neighbours       
       j=ver(vert)%vneipt(i) ; k=ver(vert)%vneitr(i)                                                                                ! Neigh points  and triangles 
       ver(j)=mp%vertex(j)   ; tri(k)=mp%triangle(k)
       If(an_ver(i+1,1).Eq.1) antig(an_ver(i+1,2))=mp%antig(an_ver(i+1,2))                                                          ! Restore the neighbour antigen base and tip position if move is rejected
       Enddo rest_neig

       PBimage_rest:If(ver(vert)%boundary .Eq.1) Then                                                                               ! If the chosen vertex is on the boundary      
       PBver_rest:Do i=1,ver(vert)%pbimno                                                                                           ! All possible images are taken care of        
       imgv=ver(vert)%pbmap(i)                    
       ver(imgv)=mp%vertex(imgv)                                                                                                    ! restore the old values
       Nver_rest:Do k=1,ver(imgv)%nonei,1
       j=ver(imgv)%vneipt(k) ; ver(j)=mp%vertex(j) 
       If(k .Lt. ver(imgv)%nonei)Then
       tri(ver(imgv)%vneitr(k))=mp%triangle(ver(imgv)%vneitr(k))
       Endif
       Enddo Nver_rest
       Enddo PBVer_rest
       Endif PBimage_rest
       Endif
       End Subroutine movevertex

!---------------------------------------------------------------------------------------------------------------------------
!                              Subroutine to move an antigen freely on the triangle surface
!---------------------------------------------------------------------------------------------------------------------------
    Function Generate_Antigen_Trial_Position(vert,triang,antigen,disp) Result(fres)
    Use module_datastruct ; Use module_randomnumber
    Implicit None
    Include '_mod_mcsmoves_interface.h'    
    Integer :: vert,triang,antigen,tr_vert(3),fres
    Real(Kind=8) :: tr_vec1(3,1),tr_vec2(3,1),disp(2)
    Integer :: i,j,n
       n=1 ; fres=0
       tr_vert(n)=vert
       Do i=1,3,1
       If(tri(triang)%vert(i).Ne. vert) Then
       n=n+1 ; tr_vert(n)=tri(triang)%vert(i)
       Endif ;  Enddo
       tr_vec1=0.5*(ver(tr_vert(2))%vcoord-ver(tr_vert(1))%vcoord)
       tr_vec2=0.5*(ver(tr_vert(3))%vcoord-ver(tr_vert(1))%vcoord)
       disp(1)=ran2(seed) ; disp(2)=ran2(seed)
       antig(antigen)%base_coord=ver(vert)%vcoord+disp(1)*tr_vec1+disp(2)*tr_vec2                     
       antig(antigen)%tip_coord=antig(antigen)%base_coord+eqm_bond_dist*ver(vert)%vnor                                                 ! Right now the flexure of the antigen is not restored (Needs some thought)
       If(antig(antigen)%tip_coord(1,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen,1,-periodic_box_length)                   ! Apply the periodic boundary condition for x comp
       If(antig(antigen)%tip_coord(1,1).Lt.0) Call PBC_Antigen_tip(antigen,1,periodic_box_length)                                      ! Apply the periodic boundary condition for x comp
       If(antig(antigen)%tip_coord(2,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen,2,-periodic_box_length)                   ! Apply the periodic boundary condition for y comp
       If(antig(antigen)%tip_coord(2,1).Lt.0) Call PBC_Antigen_tip(antigen,2,periodic_box_length)                                      ! Apply the periodic boundary condition for y comp
       If(Antigen_selfavoidance(antigen).Eq.1) fres=1	                                    				               ! check for the self avoidance of the antigen with others in its vicinity
    Return
    End Function Generate_Antigen_Trial_Position

!---------------------------------------------------------------------------------------------------------------------------
!                              Subroutine to move an antigen freely on the triangle surface
!---------------------------------------------------------------------------------------------------------------------------
    Function Displace_Antigen_on_Triangle(vert,triang,antigen) Result(fres)
    Use module_datastruct
    Implicit None
    Include '_mod_mcsmoves_interface.h'    
    Integer :: vert,triang,antigen,tr_vert(3),fres
    Real(Kind=8) :: tr_vec1(3,1),tr_vec2(3,1)
    Integer :: i,j,n
    If (triang > 0) Then                                                                                                              ! If the antigen is associated with a triangle move on it with disp1 and disp2 
       n=1 ; fres=0
       tr_vert(n)=vert
       Do i=1,3,1
       If(tri(triang)%vert(i).Ne. vert) Then
       n=n+1 ; tr_vert(n)=tri(triang)%vert(i)
       Endif ;  Enddo
       tr_vec1=0.5*(ver(tr_vert(2))%vcoord-ver(tr_vert(1))%vcoord)
       tr_vec2=0.5*(ver(tr_vert(3))%vcoord-ver(tr_vert(1))%vcoord)
       antig(antigen)%base_coord=ver(vert)%vcoord+antig(antigen)%disp1*tr_vec1+antig(antigen)%disp2*tr_vec2                     
       antig(antigen)%tip_coord=antig(antigen)%base_coord+eqm_bond_dist*ver(vert)%vnor                                                 ! Right now the flexure of the antigen is not restored (Needs some thought)
       If(antig(antigen)%tip_coord(1,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen,1,-periodic_box_length)                   ! Apply the periodic boundary condition for x comp
       If(antig(antigen)%tip_coord(1,1).Lt.0) Call PBC_Antigen_tip(antigen,1,periodic_box_length)                                      ! Apply the periodic boundary condition for x comp
       If(antig(antigen)%tip_coord(2,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen,2,-periodic_box_length)                   ! Apply the periodic boundary condition for y comp
       If(antig(antigen)%tip_coord(2,1).Lt.0) Call PBC_Antigen_tip(antigen,2,periodic_box_length)                                      ! Apply the periodic boundary condition for y comp
       If(Antigen_selfavoidance(antigen).Eq.1) fres=1	                                    				               ! check for the self avoidance of the antigen with others in its vicinity
   Else
       antig(antigen)%base_coord=ver(vert)%vcoord                                                                                      ! If the antigen is not associated with a triangle reorient it on the vertex it resides on
       antig(antigen)%tip_coord=antig(antigen)%base_coord+eqm_bond_dist*ver(vert)%vnor                                                 ! Right now the flexure of the antigen is not restored (Needs some thought)
       If(antig(antigen)%tip_coord(1,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen,1,-periodic_box_length)                   ! Apply the periodic boundary condition for x comp
       If(antig(antigen)%tip_coord(1,1).Lt.0) Call PBC_Antigen_tip(antigen,1,periodic_box_length)                                      ! Apply the periodic boundary condition for x comp
       If(antig(antigen)%tip_coord(2,1).Gt.periodic_box_length) Call PBC_Antigen_tip(antigen,2,-periodic_box_length)                   ! Apply the periodic boundary condition for y comp
       If(antig(antigen)%tip_coord(2,1).Lt.0) Call PBC_Antigen_tip(antigen,2,periodic_box_length)                                      ! Apply the periodic boundary condition for y comp
       If(Antigen_selfavoidance(antigen).Eq.1) fres=1
    Endif
    Return
    End Function Displace_Antigen_on_Triangle 
    

!---------------------------------------------------------------------------------------------------------------------------
!                              Subroutine to rearrange antigens between an old and new triangle
!---------------------------------------------------------------------------------------------------------------------------
    Subroutine Rearrange_Antigen_on_triangles(old_tri_no,new_tri_no,ant_no)
    Use module_datastruct ; Use module_randomnumber
    Implicit None
    Include '_mod_mcsmoves_interface.h'
    Integer :: i,ant_no,ver_no,new_tri_no,nei_no,old_tri_no,array_index
    If(old_tri_no .Ne. new_tri_no) Then
    Do i=1,3,1                                                                                                                      ! find the index of antigen in the triangle old_triang_no
    If (tri(old_tri_no)%antigen(i) .Eq. ant_no) array_index=i
    Enddo
    Do i= array_index,2,1                                                                                                           ! Remove the antigen from the antigen list for old triangles
    tri(old_tri_no)%antigen(i)=tri(old_tri_no)%antigen(i+1)
    tri(old_tri_no)%antigen(i+1)=0
    Enddo
    tri(old_tri_no)%nantigen = tri(old_tri_no)%nantigen-1					   	                            ! decrease number of antigen in old triangle by 1
    tri(new_tri_no)%nantigen = tri(new_tri_no)%nantigen+1
    tri(new_tri_no)%antigen(tri(new_tri_no)%nantigen)=ant_no                                                                        ! Add antigen number to the new triangle it has moved to 
    antig(ant_no)%diffus_tri=new_tri_no
    Endif
    End Subroutine Rearrange_Antigen_on_triangles
      
!---------------------------------------------------------------------------------------------------------------------------
!                              Subroutine to move an antigen freely on the triangle surface
!---------------------------------------------------------------------------------------------------------------------------
    Subroutine Antigen_diffusion_on_triangle(ant_no)
    Use module_datastruct ; Use module_randomnumber
    Implicit None
    Include '_mod_mcsmoves_interface.h'
    Integer :: i,ant_no,ver_no,new_tri_no,nei_no,old_tri_no,array_index,fres,selected_trial_move,move_rejected_flag
    Real(Kind=8) :: rosbluth_weight_old,rosbluth_weight_new,acc_prob
    Real(Kind=8) :: tri_vec1(3,1),tri_vec2(3,1),base_coord(3,1),tip_coord(3,1),energy_init_state,energy_trial_state
    Real(Kind=8) :: bias_disp1(num_bias_moves),bias_disp2(num_bias_moves),bias_energy(num_bias_moves),disp(2)
    Real(Kind=8) :: trial_probability(num_bias_moves)
    Integer ::  new_trial_triangle(num_bias_moves),old_trial_triangle(num_bias_moves),chosen_triangle


    move_rejected_flag=0
    ver_no=antig(ant_no)%vertex ; old_tri_no=antig(ant_no)%diffus_tri
    mp%antig(ant_no)=antig(ant_no)												    ! Store the old values	
    bonded_unbonded: If (bondedstate(ant_no-1).Eq.0) Then                                                                           ! If the antigen is not bonded perform a unbiased Monte Carlo move. 
      nei_no=Nint(ran2(seed)*ver(ver_no)%nonei)+1
      If(nei_no .Gt. ver(ver_no)%nonei) nei_no=ver(ver_no)%nonei                                		                    ! Choose one triangle randomly from the neighbour list
      new_tri_no=ver(ver_no)%vneitr(nei_no) ;  antig(ant_no)%diffus_tri=new_tri_no                                                  ! Store the exact triangle where the antigen is currently localized
      antig(ant_no)%disp1=ran2(seed) ; antig(ant_no)%disp2=ran2(seed)                                                               ! current random position with respect to the triangle coordinates

      If(Displace_Antigen_on_Triangle(ver_no,new_tri_no,ant_no).Eq.1)Then
      antig(ant_no)=mp%antig(ant_no) ;  Return 
      Else
      Call Rearrange_Antigen_on_triangles(old_tri_no,new_tri_no,ant_no)                                                             ! Reorder the antigen list on triangles and vice-versa
      Call update_linkcells(ant_no,&
      mp%antig(ant_no)%base_coord(1,1),mp%antig(ant_no)%base_coord(2,1),mp%antig(ant_no)%base_coord(3,1),'c')   	            ! Update the link cell
      Endif

    Else
       rosbluth_weight_old=0.0 ; rosbluth_weight_new=0.0                                                                            ! Initialize the rosenbluth weights 
       energy_init_state=antigen_reaction_energy(ant_no-1,antig(ant_no)%tip_coord(1,1))                                             ! Compute the reaction energy  (-1 takes care of indices in C++, passing the pointer)

       make_trial_triangles: If (old_tri_no .Ne. 0) Then
	 new_trial_triangle=old_tri_no	                                                                                      	    ! All the moves are made in the triangle the antigen is associated with
	 old_trial_triangle=old_tri_no
       Else
	 Do i=1,num_bias_moves,1
	 nei_no=Nint(ran2(seed)*ver(ver_no)%nonei)+1
    	 If(nei_no .Gt. ver(ver_no)%nonei) nei_no=ver(ver_no)%nonei
	 new_trial_triangle(i)=ver(ver_no)%vneitr(nei_no)                                                                            ! Generate num_bias_moves random triangle list
	 nei_no=Nint(ran2(seed)*ver(ver_no)%nonei)+1
    	 If(nei_no .Gt. ver(ver_no)%nonei) nei_no=ver(ver_no)%nonei
	 old_trial_triangle(i)=ver(ver_no)%vneitr(nei_no)
	 Enddo 
       Endif make_trial_triangles

         Do i=1,num_bias_moves,1                                              		                                            ! Choose one triangle randomly from the neighbour list
      	 new_tri_no=new_trial_triangle(i)											    ! Store the exact triangle where the antigen is currently localized
         fres=Generate_Antigen_Trial_Position(ver_no,new_tri_no,ant_no,disp)
         If (fres .Eq. 0) Then
         bias_disp1(i)=disp(1) ; bias_disp2(i)=disp(2)                                                                              ! Record the displacement made  if the move satisfies self-avoidance
         bias_energy(i)= antigen_reaction_energy(ant_no-1,antig(ant_no)%tip_coord(1,1))
         Else
         bias_disp1(i)=mp%antig(ant_no)%disp1 ; bias_disp2(i)=mp%antig(ant_no)%disp2                                                ! The old position is maintained if the move violates self-avoidance
         bias_energy(i)= energy_init_state ; 
         Endif
         trial_probability(i)=rosbluth_weight_new+exp(-beta*bias_energy(i))                                                         ! compute the rosenbluth factor for the new configuration 
         rosbluth_weight_new=rosbluth_weight_new+trial_probability(i)                                                               ! compute the rosenbluth factor for the new configuration 
         Enddo
  
         rosbluth_weight_old=rosbluth_weight_old+exp(-beta*energy_init_state)							    ! Generate trial configurations starting from the old position
         Do i=2,num_bias_moves,1
      	 new_tri_no=old_trial_triangle(i)											    ! Store the exact triangle where the antigen is currently localized
         fres=Generate_Antigen_Trial_Position(ver_no,new_tri_no,ant_no,disp)
         If (fres .Eq. 0) Then
         energy_trial_state = antigen_reaction_energy(ant_no-1,antig(ant_no)%tip_coord(1,1))
         Else
         energy_trial_state = energy_init_state
         Endif
         rosbluth_weight_new=rosbluth_weight_new+exp(-beta*energy_trial_state)                                                        ! compute the rosenbluth factor for the old configuration 
         Enddo

         selected_trial_move=Select_from_trial_moves(trial_probability,rosbluth_weight_new,num_bias_moves)

         select_a_trial: If (selected_trial_move .Eq. 0) Then
          move_rejected_flag=1
         Else
          avoid_zero_weight:If(rosbluth_weight_old .Eq. 0.0) Then
          acc_prob=rosbluth_weight_new/(1.0E-20)
          Else
	  acc_prob=rosbluth_weight_new/rosbluth_weight_old
	 Endif avoid_zero_weight
	  chosen_triangle = new_trial_triangle(selected_trial_move)

       	 rosbluth_Metropolis: If(ran2(seed) .Lt.  acc_prob)Then
	  antig(ant_no)%disp1=bias_disp1(selected_trial_move)
	  antig(ant_no)%disp2=bias_disp2(selected_trial_move)
	  fres = Displace_Antigen_on_Triangle(ver_no,chosen_triangle,ant_no)
          Call update_linkcells(ant_no,&
	  mp%antig(ant_no)%base_coord(1,1),mp%antig(ant_no)%base_coord(2,1),mp%antig(ant_no)%base_coord(3,1),'c')   		    ! Update the link cell
	 Else
	  move_rejected_flag=1 
         Endif rosbluth_metropolis
	Endif select_a_trial

      If(move_rejected_flag .Eq. 1) antig(ant_no)=mp%antig(ant_no)                                                          	    ! Restore the old attributes if the move fails 
    Endif bonded_unbonded
    End Subroutine Antigen_diffusion_on_triangle

!---------------------------------------------------------------------------------------------------------------------------
!	Function that selects a trial function for  Rosenbluth sampling  (algorithm 41, Appendix J, Frenkel and Smit)
!---------------------------------------------------------------------------------------------------------------------------    
	Function Select_from_trial_moves(weight,rosbluth_weight_fac,number_moves) Result(selected_trial)
	Use module_randomnumber ; Use module_datastruct
	Implicit None
	Integer :: number_moves,selected_trial,i
	Real(Kind=8) ::	weight(number_moves),rosbluth_weight_fac,curr_weight,rand_weight
        rand_weight = ran2(Seed)*rosbluth_weight_fac
	selected_trial=0 ; curr_weight=0
	Do i=1,number_moves,1
	curr_weight = curr_weight+weight(i)
	If(curr_weight.Lt.rand_weight) Then
		selected_trial=i ; Return
	Endif
	Enddo
        Return
	End Function Select_from_trial_moves
!---------------------------------------------------------------------------------------------------------------------------
!                               Subroutine TO CALCULATE THE ENERGY FOR FLIPPING MOVE
!---------------------------------------------------------------------------------------------------------------------------
    Subroutine flipener(cver)
    USE module_datastruct ; USE module_curvcalc
    IMPLICIT NONE
    Integer,DIMENSION(4) :: cver
    Integer :: i
    flen=0
    four_points: Do i=1,4,1
    flen=flen+(ver(cver(i))%mcur**2)*kappa*ver(cver(i))%totarea                                                                !  Helfrich contrib by vertex vert       
    Enddo four_points
    End Subroutine flipener     
!------------------------------------------------------------------------------------------------------------------------
!                               Subroutine TO CALCULATE THE ENERGY FOR VERTEX MOVES
!------------------------------------------------------------------------------------------------------------------------
    Subroutine vermovener(vno)
    USE module_datastruct ; USE module_curvcalc
    IMPLICIT NONE
    Integer :: v,i,vno,t 
    mven=0
    mv_neigh: Do i=1,ver(vno)%nonei,1
    v=ver(vno)%vneipt(i)   
    mven=mven+(ver(v)%mcur**2)*kappa*ver(v)%totarea                                                                             !  Helfrich  Energy    
    Enddo mv_neigh 
    mven=mven+(ver(vno)%mcur**2)*kappa*ver(vno)%totarea                                                                      
    End Subroutine vermovener
!------------------------------------------------------------------------------------------------------------------------
!                            Subroutine to Implement the Periodic boundary conditions for antigen tip
!------------------------------------------------------------------------------------------------------------------------
    Subroutine PBC_Antigen_tip(ano,indexx,shift)                                                                                ! ano -antigen number, indexx-component to shift, shift= ±periodic_box_length  
    USE module_datastruct ; USE module_curvcalc
    IMPLICIT NONE
    Integer :: ano,indexx
    Real(Kind=8) :: shift
    antig(ano)%tip_coord(indexx,1)=antig(ano)%tip_coord(indexx,1)+shift
    End Subroutine PBC_Antigen_tip
!------------------------------------------------------------------------------------------------------------------------
!                                Subroutine TO CHECK FOR THE ANGLE BETWEEN THE GIVEN FACES
!------------------------------------------------------------------------------------------------------------------------
   Subroutine faceangchk(trno)
   USE module_datastruct
   IMPLICIT NONE
   Integer :: trno,j,k                                                                                                         ! trno --> Triangle to be chked for
   Real(KIND=8) :: npro(1,1)                                                                                                   ! Dot product of the Normals       
   Integer,DIMENSION(3) :: ll
   norchk=0 ; ll=0  
   ll= tri(trno)%li
   Do k=1,3
   j=lin(-ll(k))%tr 
   npro=MATMUL(TRANSPOSE(tri(trno)%fnor),tri(j)%fnor)  
   If(npro(1,1).LT.fangle) Then
   norchk=1                                                                                                                    ! norchk --> 0 for the conf to be accepted
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
       v=tri(trno)%vert                                                                                                             ! Read all the vertices of the given triangle
  
       Do i=1,3,1                                                                                                                   ! Run over all the 3 vertices making up 'trno'    
       Do j=1,ver(v(i))%nonei-1                                                                                                     ! Explore the neigh triangs of 'trno'            
       neiflag=0
       vntr=ver(v(i))%vneitr(j)                                                                                                     ! vntr is the 'j' th neigh tri of vertex 'v(i)'  
       If(vntr.Ne.trno)Then                                          
         vn=tri(ver(v(i))%vneitr(j))%vert                                                                                           ! Vertices making up the triangle 'vntr'             
         If(v(1).Eq.vn(1) .Or.v(1).Eq.vn(2) .Or. v(1).Eq.vn(3)) Then                                                                ! Check if the chosen triangle and the triangle in $ 
         neiflag=neiflag+1                                                                                                          ! $ the neig of its 'i'th vertex are neighbours      
         Endif
         If(v(2).Eq.vn(1) .Or.v(2).Eq.vn(2) .Or. v(2).Eq.vn(3)) Then  
         neiflag=neiflag+1
         Endif
         If(v(3).Eq.vn(1) .Or.v(3).Eq.vn(2) .Or. v(3).Eq.vn(3)) Then
         neiflag=neiflag+1
         Endif

         If(neiflag.Eq.2) Then                                                                                                       ! A value of 2 for neiflag is a sign of neig of trno & vntr 
         npro=MATMUL(TRANSPOSE(tri(trno)%fnor),tri(vntr)%fnor)  
         If(npro(1,1).LT.fangle) Then
         norchk=1                                                                                                                    ! norchk --> 0 for the conf to be accepted
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

       i=tri(tr)%vert(1) ; j=tri(tr)%vert(2) ; k=tri(tr)%vert(3)                                                                       ! The vertices that make up the triangle tr 
       r1=ver(i)%vcoord ; r2=ver(j)%vcoord ; r3=ver(k)%vcoord                                                                          ! Their corresponding coordinates
       r21=r2-r1 ; r31=r3-r1

       ax=(r21(2,1)*r31(3,1)-r21(3,1)*r31(2,1))*0.5
       ay=(r21(3,1)*r31(1,1)-r21(1,1)*r31(3,1))*0.5
       az=(r21(1,1)*r31(2,1)-r21(2,1)*r31(1,1))*0.5
       area=SQRT(ax**2+ay**2+az**2)

       tri(tr)%ar=area ; ax=ax/area ; ay=ay/area ; az=az/area
       tri(tr)%fnor=RESHAPE((/ax,ay,az/),SHAPE(tri(tr)%fnor))

       RETURN
       End Subroutine onlyarea
!-------------------------------------------------------------------------------------------------------------------------
!                                   Function TO CHECK THE MAX And MIN BONDLENGTH
!-------------------------------------------------------------------------------------------------------------------------
      Function blcheck(ver1,ver2) result(fres) 
      USE module_datastruct ; USE module_makesurface
      IMPLICIT NONE
      Real(KIND=8)::bl
      Integer::ver1,ver2,fres                                                                                                           !fres=1 ===>false ,0===>true
      fres=0 
      bl=SUM((ver(ver1)%vcoord-ver(ver2)%vcoord)**2)
      If (bl.GT.blUcut.OR. bl.LE.blLcut) fres=1            
      RETURN
      End Function blcheck

!-------------------------------------------------------------------------------------------------------------------------
!                                   Function TO CHECK THE MINIMUM BONDLENGTH
!-------------------------------------------------------------------------------------------------------------------------
      Function blcheckmin(ver1,ver2) result(fres) 
      USE module_datastruct ; USE module_makesurface
      IMPLICIT NONE
      Real(KIND=8)::bl
      Integer::ver1,ver2,fres                                                                                                           ! fres=1 ===>false ,0===>true
      fres=0 
      bl=SUM((ver(ver1)%vcoord-ver(ver2)%vcoord)**2)
      If(bl.LE.blLcut) fres=1            
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
      Integer::i,v,fres,mem                                                                                                         ! fres=1 ===>false ,0===>true
      Integer,Pointer :: ring(:),lscl(:),head(:)
      fres=0

      If(v.Le.nver)Then
       Call C_F_Pointer(getcells(1,v,'m'),ring,[1])                                                                                 ! convert the C-pointer returned by neighcells to an one-D fortran pointer ring     
       Call C_F_Pointer(get_lscllist('m'),lscl,[1])                                                                                 ! convert the C-pointer returned by get_lscllist to an one-D fortran pointer lscl   
       Call C_F_Pointer(get_headlist('m'),head,[1])                                                                                 ! convert the C-pointer returned by get_headlist to an one-D fortran pointer head   

       Do i=1,27,1                                                                                                                  ! Perform the self avoidance check over the one ring cells around the cell containing v 
        mem=head(ring(i))                                                                                                           ! get the first membrane particle in cell i 
        Do While(mem .Ne. -1)
         blen=SUM((ver(v)%vcoord-ver(mem+1)%vcoord)**2)                                                                             ! Compute distance if the membrane particle is present 
         If((blen.Le.blLcut) .And. (mem+1.Ne.v))Then                         
         fres=1 ; Return                                                                                                            ! Return if minimum bondlength condition is violated even once 
         Endif
         mem=lscl(mem+1)                                                                                                            ! Move to the next particle 
        Enddo
        Enddo
        Return
      Else                                                                                                                          ! Link Cells have not been implemented for the Phantoms 
        i=1
        Do while(i.Lt.pbnum)                                                                                                        ! Hence we use the traditional method to check over all vertices
        bl=SUM((ver(v)%vcoord-ver(i)%vcoord)**2)
         If((bl.LE.blLcut) .And. (i.Ne.v) ) Then
           fres=1  ; Return                                                                                                         ! Return fres=1 even if self avoidance is violated once 
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
      Use, Intrinsic :: ISO_C_BINDING 
      IMPLICIT NONE
      Include '_mod_getcoordinates_interface.h'

      Real(KIND=8)::bl,obj_coord(3,1)
      Integer::i,v,fres                                                                                                             ! fres=1 ===>false ,0===>true      
      fres=0
      Do i=0,num_Nanocarrier-1,1                                                                                                    ! start from zero since we are accessing data from c++
      obj_coord=Reshape((/getx(i,'v'),gety(i,'v'),getz(i,'v')/),(/3,1/))                                                            ! Coordinate of the Nanocarrier i  
      blen=Sum((ver(v)%vcoord-obj_coord)**2)
      If(blen .Lt. ves_soft_rad_sq)Then                                                                                             ! A membrane particle cannot be closer than soft radius 
      fres=1 ; Return
      Endif
      Enddo
      Return
      End Function Selfavoidance_membrane_Nanocarrier

!-------------------------------------------------------------------------------------------------------------------------
!                     Check Self Avoidance of Antigen with Nanocarrier
!-------------------------------------------------------------------------------------------------------------------------
      Function Selfavoidance_Antigen_Nanocarrier(v) result(fres) 
      USE module_datastruct ; USE module_makesurface
      Use, Intrinsic :: ISO_C_BINDING 
      IMPLICIT NONE
      Include '_mod_getcoordinates_interface.h'

      Real(KIND=8)::bl,obj_coord(3,1)
      Integer::i,v,fres                                                                                                             ! fres=1 ===>false ,0===>true      
      fres=0
      Do i=0,num_Nanocarrier-1,1
      obj_coord=Reshape((/getx(i,'v'),gety(i,'v'),getz(i,'v')/),(/3,1/))                                                            ! Coordinate of the Nanocarrier i  
      blen=Sum((antig(v)%tip_coord-obj_coord)**2)
      If(blen .Lt. ves_rad_sq)Then                                                                                                  ! Membrae-Antigen distance cutoff change from soft radius to radius to test freezing 
      fres=1 ; Return
      Endif
      Enddo
      Return
      End Function Selfavoidance_Antigen_Nanocarrier

!-------------------------------------------------------------------------------------------------------------------------
!                     Check Self Avoidance of Antigen with another antigen
!-------------------------------------------------------------------------------------------------------------------------
      Function Selfavoidance_Antigen_Antigen(a1,a2) result(fres) 
      USE module_datastruct ; USE module_makesurface
      IMPLICIT NONE
      Integer::a1,a2,fres                                                                                                             ! fres=1 ===>false ,0===>true      
      Real(KIND=8):: base_dist,tip_dist
      fres=0
      base_dist=Sum((antig(a1)%base_coord-antig(a2)%base_coord)**2)                                                                 ! Distance between the base of antigen a1 and a2
      tip_dist=Sum((antig(a1)%tip_coord-antig(a2)%tip_coord)**2)                                                                    ! Distance between the tip of antigen a1 and a2
      If((base_dist .Lt. antigen_overlap_dist) .Or. (tip_dist .Lt. antigen_overlap_dist)) fres=1                                    ! Antigen-Antigen distance cutoff 
      Return
      End Function Selfavoidance_Antigen_Antigen

!-------------------------------------------------------------------------------------------------------------------------
!                                   Function to check selfavoidance of one antigen with all
!-------------------------------------------------------------------------------------------------------------------------
      Function Antigen_selfavoidance(a1) result(fres) 
      USE module_datastruct ; USE module_makesurface
      Use, Intrinsic :: ISO_C_BINDING 
      IMPLICIT NONE
      Include '_mod_linkcells_interface.h'

      Real(KIND=8)::base_dist,tip_dist
      Integer::i,v,fres,mem,a1                                                                                                      ! fres=1 ===>false ,0===>true
      Integer,Pointer :: ring(:),lscl(:),head(:)
      fres=0

       Call C_F_Pointer(getcells(1,a1,'c'),ring,[1])                                                                                 ! convert the C-pointer returned by neighcells to an one-D fortran pointer ring     
       Call C_F_Pointer(get_lscllist('c'),lscl,[1])                                                                                 ! convert the C-pointer returned by get_lscllist to an one-D fortran pointer lscl   
       Call C_F_Pointer(get_headlist('c'),head,[1])                                                                                 ! convert the C-pointer returned by get_headlist to an one-D fortran pointer head   

       Do i=1,27,1                                                                                                                  ! Perform the self avoidance check over the one ring cells around the cell containing a1
        mem=head(ring(i))                                                                                                           ! get the first antigen particles in cell mem
        Do While(mem .Ne. -1)
         If(mem+1 .Ne. a1)Then
         base_dist=SUM((antig(a1)%base_coord-antig(mem+1)%base_coord)**2)                                                           ! Compute distance between base of two antigens
         tip_dist=SUM((antig(a1)%tip_coord-antig(mem+1)%tip_coord)**2)                                                              ! Compute distance between tip of two antigens
         If((base_dist.Lt.antigen_overlap_dist).Or.(tip_dist .Lt.antigen_overlap_dist))Then                         
         fres=1 ; Return                                                                                                            ! Return if self avoidance condition is violated even once 
         Endif
	 Endif
         mem=lscl(mem+1)                                                                                                            ! Move to the next antigen 
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
      Integer:: cellno,ant_no,i,mem
      Character :: ch 
      Integer,Pointer :: ring(:),lscl(:),head(:)
      
      Call C_F_Pointer(getcells(1,ant_no,ch),ring,[1])                                                                              ! convert the C-pointer returned by neighcells to an one-D fortran pointer ring     
      Call C_F_Pointer(get_lscllist(ch),lscl,[1])                                                                                   ! convert the C-pointer returned by get_lscllist to an one-D fortran pointer lscl   
      Call C_F_Pointer(get_headlist(ch),head,[1])                                                                                   ! convert the C-pointer returned by get_headlist to an one-D fortran pointer head   

        Print*,'The antigen number is checked for is ',ant_no
        Do i=1,27,1    														    ! Perform the self avoidance check over the one ring cells around the cell containing a1
        mem=head(ring(i))                                                                                                           ! get the first antigen particles in cell mem
        Do While(mem .Ne. -1)
	if (mem+1 .Eq. ant_no) Then
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
      Integer:: cellno,ant_no,i,mem,fres
      Character :: ch 
      Integer,Pointer :: ring(:),lscl(:),head(:)
      fres=-1
      
      Call C_F_Pointer(getcells(1,ant_no,ch),ring,[1])                                                                              ! convert the C-pointer returned by neighcells to an one-D fortran pointer ring     
      Call C_F_Pointer(get_lscllist(ch),lscl,[1])                                                                                   ! convert the C-pointer returned by get_lscllist to an one-D fortran pointer lscl   
      Call C_F_Pointer(get_headlist(ch),head,[1])                                                                                   ! convert the C-pointer returned by get_headlist to an one-D fortran pointer head   

        Do i=1,27,1    														    ! Perform the self avoidance check over the one ring cells around the cell containing a1
        mem=head(ring(i))                                                                                                           ! get the first antigen particles in cell mem
        Do While(mem .Ne. -1)
	if (mem+1 .Eq. ant_no) Then
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
      Integer:: cellno,obj_no,i,mem,fres
      Character :: ch 
      Integer,Pointer :: ring(:),lscl(:),head(:)
      fres=0
      
      Call C_F_Pointer(get_lscllist(ch),lscl,[1])                                                                                   ! convert the C-pointer returned by get_lscllist to an one-D fortran pointer lscl   
      Call C_F_Pointer(get_headlist(ch),head,[1])                                                                                   ! convert the C-pointer returned by get_headlist to an one-D fortran pointer head   

      mem=head(cellno+1)
      Do while(mem .Ne. -1) 
	If (mem+1 .Eq. obj_no) Then                                                                                                 ! check if an object of type ch is present in cell cellno  
	fres=1															    ! If present return 1 else return 0
	Return
	Endif	
	mem=lscl(mem+1)
      Enddo
      End function check_ifin_linkcell
     
      End MODULE module_mcsmoves
      
