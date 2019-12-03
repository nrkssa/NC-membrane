!===========================================================================================================================
!===========================================================================================================================
!        $modcurv                    MODULE TO CALCULATE THE CURVATURE And OTHER RELATED QUANTITIES
!===========================================================================================================================
!===========================================================================================================================
      MODULE  module_curvcalc
      IMPLICIT NONE
      Real(KIND=8):: unitmat(3,3),zdir(3,1)
      Real(KIND=8),DIMENSION(10,3):: tanproj 
      
      contains

!-----------------------------------------------------------------------------------------------------------------------------
!                       SUBROUTINE TO INITIALIZE CONSTANTS FOR THIS MODULE
!-----------------------------------------------------------------------------------------------------------------------------
      Subroutine Initialize_curvcalc()
      Implicit None
      unitmat=reshape([1.,0.,0.,0.,1.,0.,0.,0.,1.],(/3,3/))
      zdir=unitmat(1:3,3:3)
      End Subroutine Initialize_curvcalc
!-----------------------------------------------------------------------------------------------------------------------------
!                       SUBROUTINE TO COMPUTE THE CURVATURE AT A VERTEX VERNU
!-----------------------------------------------------------------------------------------------------------------------------

       Subroutine normalcalc(vertex_no)
       USE module_datastruct 
       IMPLICIT NONE
       Integer :: i,ip1,im1,j,vertex_no,tr1,tr2                    								    ! Loop integers
       Real(KIND=8),DIMENSION(3,1)::nor,rtan,amat,smat    
       Real(KIND=8),DIMENSION(3,3)::pmat,evec,cmat,hmat                         					            ! g2l is the global to local frame conversion matrix
       Real(KIND=8),DIMENSION(1,1)::wei,xc,an  
       Real(KIND=8)::fwei,p,q,r,s,insq,e1,e2,elen,dih,emcu,si 
       Real(KIND=8),DIMENSION(3,1) :: en,euv,v1,v2,cpm,fn1,fn2       
       nor=0.0 ; cmat=0.0 ; hmat=0.0 ; evec=0.0 ; fwei=0.0

       weightcalc:Do i=1,ver(vertex_no)%nonei,1
       j=ver(vertex_no)%vneitr(i)
       fwei=fwei+tri(j)%ar                                             								    ! Total area of all neigh triangles
       Enddo weightcalc


       normal : Do i=1,ver(vertex_no)%nonei
       j=ver(vertex_no)%vneitr(i)
       nor=nor+(tri(j)%ar/fwei)*tri(j)%fnor
       Enddo normal


       nor=nor/SQRT(SUM(nor**2))                                       								    ! Normalize the calculated normal.( Global Frame)
       ver(vertex_no)%vnor=nor
       pmat=(unitmat-MATMUL(nor,TRANSPOSE(nor)))                      								    ! Tangent plane corresponding to normal at vertex_no 

       Do i=1,ver(vertex_no)%nonei
       im1=i-1 ; If(i.EQ.1) im1=ver(vertex_no)%nonei
       ip1=i+1 ; If(i.EQ.ver(vertex_no)%nonei) ip1=1

       v2=ver(ver(vertex_no)%vneipt(i))%vcoord-ver(vertex_no)%vcoord

       elen=SQRT(SUM(v2**2))                                       								    ! Length of the edge from vertex_no --> i    
       v2=v2/elen ; euv=v2                                        								    ! Unit vector along the edge and along i 
       tr1=ver(vertex_no)%vneitr(im1) ; tr2=ver(vertex_no)%vneitr(i)    							    ! Two triangles sharing the edge 
       fn1=tri(tr1)%fnor ; fn2=tri(tr2)%fnor   
       an=MATMUL(TRANSPOSE(fn1),fn2)
       If(an(1,1).GT.1.0000000) an(1,1)=1.0000000000000
       
       v1(1,1)=fn1(2,1)*fn2(3,1)-fn1(3,1)*fn2(2,1)              								    ! Calculate (N1 X N2) 
       v1(2,1)=fn1(3,1)*fn2(1,1)-fn1(1,1)*fn2(3,1)
       v1(3,1)=fn1(1,1)*fn2(2,1)-fn1(2,1)*fn2(1,1)

       si=1.000 
       xc=MATMUL(TRANSPOSE(euv),v1)                                								    ! E.(N1 X N2)           
       si=SIGN(si,xc(1,1))                                        								    ! Sign of E.(N1 X N2)   
       dih=pi+si*acos(an(1,1))                                     								    ! Signed Dihedral angle 
        
       en=tri(tr1)%fnor+tri(tr2)%fnor                            								    ! Normal along the above edge    
       en=en/SQRT(SUM(en**2))                                      								    ! Normalized                     
       wei=MATMUL(TRANSPOSE(nor),en)
       emcu=2*elen*cos(dih*0.5)                                   								    ! Mean curvature along a the edge 
       cpm(1,1)=(euv(2,1)*en(3,1)-euv(3,1)*en(2,1))
       cpm(2,1)=(euv(3,1)*en(1,1)-euv(1,1)*en(3,1))               								    ! Direction orthogonal to edge and edge normal 
       cpm(3,1)=(euv(1,1)*en(2,1)-euv(2,1)*en(1,1))
       rtan=MATMUL(pmat,cpm) ; rtan=rtan/SQRT(SUM(rtan**2))
       cmat=cmat+0.5*wei(1,1)*emcu*MATMUL(rtan,TRANSPOSE(rtan))
       Enddo

        amat=zdir+nor ; smat=zdir-nor  
        householder:If(SUM(amat**2) .GT. SUM(smat**2))Then      								    ! Use Householder trans to reduce to a 2X2 
        amat=amat/SQRT(SUM(amat**2))
        hmat=-(unitmat-2*MATMUL(amat,TRANSPOSE(amat)))
        Else
        smat=smat/SQRT(SUM(smat**2))
        hmat=(unitmat-2*MATMUL(smat,TRANSPOSE(smat)))
        Endif householder

        pmat=0
        pmat=MATMUL(TRANSPOSE(hmat),MATMUL(cmat,hmat))           								    ! Diagonalize the constructed matrix 

        p=pmat(1,1) ; q=pmat(1,2)                               								    ! Components of the 2X2 minor 
        r=pmat(2,1) ; s=pmat(2,2) 


        non_diagonal:If(q.NE.0.0 .And. r.NE.0.0)Then            								    ! Non diagonal matrices eigen values 
        insq=(p+s)**2-4*(p*s-q*r)

         If((insq .GT. 0.0))Then                               									    ! Complex values are avoided 
          If((p+s).LT.0)Then                                  									    ! Pick up the largest eigenvalue 
          e1=((p+s)-SQRT(insq))*0.5                          									    ! Eigenvalues corresponding to largest eigenvalue
          e2=((p+s)+SQRT(insq))*0.5
          Else
          e1=((p+s)+SQRT(insq))*0.5                          									    ! Eigenvalues 
          e2=((p+s)-SQRT(insq))*0.5                         								            ! Set up such that e1 is the largest eigenvalue
          Endif
         Else
         e1=(p+s)*0.5                                       									    ! Degenrate eigenvalues 
         e2=e1
         Endif
 
        Else
         If(p.GT.s)Then                                          								    ! Picking up the largest eigenvalue for diagonal matrix 
         e1=p ; e2=s
         evec=unitmat                                            								    ! The eigenvector is same as the unit matrix 
         Else 
         e1=s ; e2=p
         Endif
        Endif non_diagonal

         If(abs(e1).LT.10E-10) e1=0
         If(abs(e2).LT.10E-10) e2=0

         ver(vertex_no)%cur1=e1/ver(vertex_no)%totarea                   							    ! Principal curvature 1                       
         ver(vertex_no)%cur2=e2/ver(vertex_no)%totarea                      							    ! Principal curvature 2                       
         ver(vertex_no)%mcur=(ver(vertex_no)%cur1+ver(vertex_no)%cur2)								    ! Mean curvature (theorema egregium of Gauss) (2H)
        End Subroutine normalcalc

!-----------------------------------------------------------------------------------------------------------------------------
!                       SUBROUTINE TO COMPUTE THE CURVATURE AT A VERTEX VERNU
!-----------------------------------------------------------------------------------------------------------------------------
       Subroutine Compute_Antigen_orientation(antigen_num)
       USE module_datastruct 
       IMPLICIT NONE
       Integer :: vertex_no,antigen_num                          							            ! Loop integers
       Real(KIND=8) :: proj_mat(3,3),antigen_vec(3,1),proj_antigen_vec(3,1),costheta(1,1),cosphi1(1,1),cosphi2(1,1)
       Real(Kind=8) :: proj_vec_len,ant_length

       proj_mat=0.0
       vertex_no=antig(antigen_num)%vertex ; ant_length = antig(antigen_num)%length
       Call Compute_Principal_directions(vertex_no,proj_mat)                   	    						    ! update principal direction and projection matrix
       antigen_vec=(antig(antigen_num)%tip_coord-antig(antigen_num)%base_coord)/ant_length                                          ! Normalized Antigen vector 
       proj_antigen_vec=Matmul(proj_mat,antigen_vec)                                                                                ! Its projection on the tangent plane
       proj_vec_len=Sum(proj_antigen_vec**2)
       if (proj_vec_len .Gt.0.0) proj_antigen_vec=proj_antigen_vec/Sqrt(proj_vec_len)                                               ! Normalize the projection vector
       costheta = Matmul(Transpose(ver(vertex_no)%vnor),antigen_vec)                                                                ! Cosine of Polar angle wrt  the normal
       cosphi1 = Matmul(Transpose(ver(vertex_no)%t1),proj_antigen_vec)                                                              ! Cosine of Azimuthal angle wrt t1
       cosphi2 = Matmul(Transpose(ver(vertex_no)%t2),proj_antigen_vec)                                                              ! Cosine of Azimuthal angle wrt t2 
       if (costheta(1,1).Gt.1.0) costheta=1.0
       if (costheta(1,1).Lt.-1.0) costheta=-1.0
       antig(antigen_num)%theta=acos(costheta(1,1))
       if(cosphi2(1,1).Ge.0) Then                                                                                                   ! Projection is in the first and second quadrants 
	       if (cosphi1(1,1).Gt.1.0) cosphi1=1.0
	       if (cosphi1(1,1).Lt.-1.0) cosphi1=-1.0 
	       antig(antigen_num)%phi=acos(cosphi1(1,1))
       Else                                                                                                                         ! In the third and fourth quadrants 
	       if (cosphi1(1,1).Gt.1.0) cosphi1=1.0
	       if (cosphi1(1,1).Lt.-1.0) cosphi1=-1.0
	       antig(antigen_num)%phi=2*pi-acos(cosphi1(1,1))
       Endif

       End Subroutine Compute_Antigen_orientation

!-----------------------------------------------------------------------------------------------------------------------------
!                       SUBROUTINE TO COMPUTE THE CURVATURE AT A VERTEX VERNU
!-----------------------------------------------------------------------------------------------------------------------------

       Subroutine Compute_Principal_directions(vertex_no,proj_mat)
       USE module_datastruct 
       IMPLICIT NONE
       Integer :: i,ip1,im1,vertex_no,tr1,tr2                                                                                       ! Loop integers
       Real(KIND=8),DIMENSION(3,1)::nor,rtan,amat,smat    
       Real(KIND=8),DIMENSION(3,3)::pmat,evec,cmat,hmat,proj_mat                                                                    ! g2l is the global to local frame conversion matrix
       Real(KIND=8),DIMENSION(1,1)::wei,xc,an  
       Real(KIND=8)::fwei,p,q,r,s,insq,e1,e2,elen,dih,emcu,si 
       Real(KIND=8),DIMENSION(3,1) :: en,euv,v1,v2,cpm,fn1,fn2       
       nor=0 ; cmat=0 ; hmat=0 ; evec=0 ; fwei=0

       pmat=(unitmat-MATMUL(ver(vertex_no)%vnor,TRANSPOSE(ver(vertex_no)%vnor)))                                                    ! Tangent plane corresponding to normal at vertex_no 

       Do i=1,ver(vertex_no)%nonei
       im1=i-1 ; If(i.EQ.1) im1=ver(vertex_no)%nonei
       ip1=i+1 ; If(i.EQ.ver(vertex_no)%nonei) ip1=1

       v2=ver(ver(vertex_no)%vneipt(i))%vcoord-ver(vertex_no)%vcoord

       elen=SQRT(SUM(v2**2)) 		                                   ! Length of the edge from vertex_no --> i    
       v2=v2/elen ; euv=v2              		                   ! Unit vector along the edge and along i 
       tr1=ver(vertex_no)%vneitr(im1) ; tr2=ver(vertex_no)%vneitr(i)  		   ! Two triangles sharing the edge 
       fn1=tri(tr1)%fnor ; fn2=tri(tr2)%fnor   
       an=MATMUL(TRANSPOSE(fn1),fn2)
       If(an(1,1).GT.1.0000000) an(1,1)=1.0000000000000
       
       v1(1,1)=fn1(2,1)*fn2(3,1)-fn1(3,1)*fn2(2,1)              	   ! Calculate (N1 X N2) 
       v1(2,1)=fn1(3,1)*fn2(1,1)-fn1(1,1)*fn2(3,1)
       v1(3,1)=fn1(1,1)*fn2(2,1)-fn1(2,1)*fn2(1,1)

       si=1.000 
       xc=MATMUL(TRANSPOSE(euv),v1)                            		    ! E.(N1 X N2)           
       si=SIGN(si,xc(1,1))                                         ! Sign of E.(N1 X N2)   
       dih=pi+si*acos(an(1,1))                                     ! Signed Dihedral angle 
        
       en=tri(tr1)%fnor+tri(tr2)%fnor                              ! Normal along the above edge    
       en=en/SQRT(SUM(en**2))                                      ! Normalized                     
       wei=MATMUL(TRANSPOSE(nor),en)
       emcu=2*elen*cos(dih*0.5)                                    ! Mean curvature along a the edge 
       cpm(1,1)=(euv(2,1)*en(3,1)-euv(3,1)*en(2,1))
       cpm(2,1)=(euv(3,1)*en(1,1)-euv(1,1)*en(3,1))                ! Direction orthogonal to edge and edge normal 
       cpm(3,1)=(euv(1,1)*en(2,1)-euv(2,1)*en(1,1))
       rtan=MATMUL(pmat,cpm) ; rtan=rtan/SQRT(SUM(rtan**2))
       cmat=cmat+0.5*wei(1,1)*emcu*MATMUL(rtan,TRANSPOSE(rtan))
       Enddo

        amat=zdir+nor ; smat=zdir-nor  
        householder:If(SUM(amat**2) .GT. SUM(smat**2))Then       ! Use Householder trans to reduce to a 2X2 
        amat=amat/SQRT(SUM(amat**2))
        hmat=-(unitmat-2*MATMUL(amat,TRANSPOSE(amat)))
        Else
        smat=smat/SQRT(SUM(smat**2))
        hmat=(unitmat-2*MATMUL(smat,TRANSPOSE(smat)))
        Endif householder

        proj_mat=pmat
	pmat=0.0
        pmat=MATMUL(TRANSPOSE(hmat),MATMUL(cmat,hmat))           ! Diagonalize the constructed matrix 

        p=pmat(1,1) ; q=pmat(1,2)                                ! Components of the 2X2 minor 
        r=pmat(2,1) ; s=pmat(2,2) 


        non_diagonal:If(q.NE.0.0 .And. r.NE.0.0)Then             ! Non diagonal matrices eigen values 
        insq=(p+s)**2-4*(p*s-q*r)

         If((insq .GT. 0.0))Then                                  ! Complex values are avoided 
          If((p+s).LT.0)Then                                  ! Pick up the largest eigenvalue 
          e1=((p+s)-SQRT(insq))*0.5                          ! Eigenvalues corresponding to largest eigenvalue
          e2=((p+s)+SQRT(insq))*0.5
          Else
          e1=((p+s)+SQRT(insq))*0.5                          ! Eigenvalues 
          e2=((p+s)-SQRT(insq))*0.5                          ! Set up such that e1 is the largest eigenvalue
          Endif
         Else
         e1=(p+s)*0.5                                          ! Degenrate eigenvalues 
         e2=e1
         Endif
 
        evec(1,1)=q/SQRT(q**2+(e1-p)**2)
        evec(2,1)=(e1-p)/SQRT(q**2+(e1-p)**2)
        evec(1,2)=-evec(2,1)
        evec(2,2)=evec(1,1)
        evec(3,3)=1      

        Else
         If(p.GT.s)Then                                          ! Picking up the largest eigenvalue for diagonal matrix 
         e1=p ; e2=s
         evec=unitmat                                            ! The eigenvector is same as the unit matrix 
         Else
         e1=s ; e2=p
        evec(:,1:1)=unitmat(:,2:2) 
        evec(:,2:2)=unitmat(:,1:1) 
        evec(:,3:3)=unitmat(:,3:3) 
        Endif
        Endif non_diagonal

        If(abs(e1).LT.10E-10) e1=0
        If(abs(e2).LT.10E-10) e2=0

        ver(vertex_no)%t1(1:3,1:1)=Matmul(hmat,evec(1:3,1:1))
        ver(vertex_no)%t2(1:3,1:1)=Matmul(hmat,evec(1:3,2:2))           ! eigenvectors
        End Subroutine Compute_Principal_directions
	

!-------------------------------------------------------------------------------------------------------------------------
!               Compute the householder matrix that rotates an arbitrary vector to z direction
!--------------------------------------------------------------------------------------------------------------------------
        Subroutine Compute_Vertex_Householdermatrix(vertex_no) 
        Use module_datastruct
        Implicit None
        Real(Kind=8) :: amat(3,1),smat(3,1)
        Integer :: vertex_no 
        amat=zdir+ver(vertex_no)%vnor ; smat=zdir-ver(vertex_no)%vnor
        If(SUM(amat**2) .GT. SUM(smat**2))Then                          ! Use Householder trans to reduce to a 2X2 
        amat=amat/SQRT(SUM(amat**2))
        ver(vertex_no)%HHM=-(unitmat-2*MATMUL(amat,TRANSPOSE(amat)))
        Else
        smat=smat/SQRT(SUM(smat**2))
        ver(vertex_no)%HHM=(unitmat-2*MATMUL(smat,TRANSPOSE(smat)))
        Endif
        End Subroutine Compute_Vertex_HouseholderMatrix

!-------------------------------------------------------------------------------------------------------------------------
!               Compute the householder matrix that rotates an arbitrary vector to z direction
!--------------------------------------------------------------------------------------------------------------------------
        Subroutine compute_Triangle_Householdermatrix(triangle_no) 
        Use module_datastruct
        Implicit None
        Real(Kind=8) :: amat(3,1),smat(3,1)
        Integer :: triangle_no
        amat=zdir+tri(triangle_no)%fnor ; smat=zdir-tri(triangle_no)%fnor
        If(SUM(amat**2) .GT. SUM(smat**2))Then                          ! Use Householder trans to reduce to a 2X2 
        amat=amat/SQRT(SUM(amat**2))
        tri(triangle_no)%HHM=-(unitmat-2*MATMUL(amat,TRANSPOSE(amat)))
        Else
        smat=smat/SQRT(SUM(smat**2))
        tri(triangle_no)%HHM=(unitmat-2*MATMUL(smat,TRANSPOSE(smat)))
        Endif

        End Subroutine Compute_Triangle_HouseholderMatrix

!-------------------------------------------------------------------------------------------------------------------------
!               Compute the householder matrix that rotates an arbitrary vector to z direction
!--------------------------------------------------------------------------------------------------------------------------
        Subroutine compute_Householdermatrix(objno,ch) 
        Use module_datastruct
        Implicit None
        Integer :: objno
	Character :: ch
	If (ch .Eq. 'n') Then
		Call Compute_Vertex_Householdermatrix(objno)
	Else If (ch .Eq. 'f') Then
		Call Compute_Triangle_Householdermatrix(objno)
	Else
		Print*,str_red,'Passed a invalid character ',ch,' to compute_householder'
		Stop
	Endif
	Return
        End Subroutine Compute_HouseholderMatrix
	
        End MODULE module_curvcalc 
