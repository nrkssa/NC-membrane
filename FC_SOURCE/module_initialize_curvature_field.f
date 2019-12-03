      Module module_Initialize_curvature_field
      Integer,Allocatable,Dimension(:,:) :: Vneighlist
      Contains

!===========================================================================================================================
!                  Subroutine to initialize epsins and create phantom epsins 
!===========================================================================================================================
        Subroutine Initialize_spontaneous_curvature(ver_no,nver_curzero,curzero,curzero_mode,ann_begin,ann_end)                  
        Use module_datastruct                         
        Implicit None    
	Integer, Optional :: ver_no,ann_begin,ann_end
	Integer :: nver_curzero	
	Real(Kind=8) :: curzero
	Character(100) :: curzero_mode
        ver(:)%czero=0  ; ver(:)%czero_flag=0                                                                                       ! Initialize all vertex to epsin mapping
	If(Trim(Adjustl(curzero_mode)).Eq. 'CIRCULAR') Call Initialize_CircularPatch(ver_no,nver_curzero,curzero)
	If(Trim(Adjustl(curzero_mode)).Eq. 'RANDOM') Call Initialize_Randomly(nver_curzero,curzero)
	If(Trim(Adjustl(curzero_mode)).Eq. 'ANNULUS') Call Initialize_Annulus(ver_no,ann_begin,ann_end,curzero)
	If(Trim(Adjustl(curzero_mode)).Eq. 'ANNULUS_CONC') Call Initialize_Annulus_conc(nver_curzero,ver_no,ann_begin,curzero)
        End subroutine Initialize_Spontaneous_curvature	
!---------------------------------------------------------------------------------------------------------------------------
!                                             Subroutine to make a patch containing epsins
!---------------------------------------------------------------------------------------------------------------------------
        Subroutine Initialize_CircularPatch(ver_no,nver_curzero,curzero)                                                               ! Makes a patch with scalar field 1 and czero_flag field 1 
        Use module_datastruct
        Implicit None
        Integer :: i,j,ncount,ncount1,verlist(nver)
	Integer :: ver_no,nver_curzero
	Real(Kind=8) :: curzero

        Ver(:)%czero=0.0  ; ver(:)%czero_flag=-1 ;   verlist=0 
        ncount=1 ; ncount1=1 ; verlist(ncount1)=ver_no               
        ver(ver_no)%czero_flag=1 ; ver(ver_no)%czero=curzero   
         Do While((ncount.Lt.nver_curzero) .And. (ncount1 .Lt. nver))
         Do i=1,ver(verlist(ncount1))%nonei
         j=ver(verlist(ncount1))%vneipt(i)
         If(ver(j)%czero_flag .Eq. -1)Then
         ncount=ncount+1 ; ver(j)%czero_flag=1 ; verlist(ncount)=j
         ver(j)%czero=curzero
         Endif
         Enddo
         ncount1=ncount1+1
         Enddo
        End Subroutine Initialize_CircularPatch

!---------------------------------------------------------------------------------------------------------------------------
!                 Subroutine to distribute the epsins at random locations
!---------------------------------------------------------------------------------------------------------------------------
        Subroutine Initialize_Randomly(nver_curzero,curzero)                                                                          ! Makes a patch with scalar field 1 and czero_flag field 1 
        Use module_datastruct ; Use module_randomnumber
        Implicit None
        Integer :: ncount,rand
	Integer :: nver_curzero
	Real(Kind=8) :: curzero

        ver(:)%czero=0.0 ; ver(:)%czero_flag=-1 ; ncount=0
        Do While((ncount.Lt.nver_curzero))
          rand=Nint(ran2(Seed)*nver)+1
          If(rand.Gt.nver) rand=nver
          If(ver(rand)%czero_flag .Eq. -1)Then
          ver(rand)%czero_flag=1
          ncount=ncount+1
          ver(rand)%czero=curzero
          Endif
        Enddo
        End Subroutine Initialize_Randomly
!---------------------------------------------------------------------------------------------------------------------------
!                              Initialize  vector field in an annulus
!---------------------------------------------------------------------------------------------------------------------------
        Subroutine Initialize_Annulus(ver_no,ann_begin,ann_end,curzero)                                                    ! Makes a patch with scalar field 1 and czero_flag field 1 
        Use module_datastruct 
        Implicit None
        Integer :: i,ann_no
        Real(Kind=8) :: curzero
        Integer :: ver_no,ann_begin,ann_end

        Call Complete_neighbourhood(ver_no,'RETURN_LIST')
        ver(:)%czero_flag=-1 ; ver(:)%czero=0.0
	ann_no=ann_begin

	Do while(ann_no .Le. ann_end) 
	i=1
	Do while ((i.Le.nver))
        	If(Vneighlist(i,2).Eq. ann_no )Then                                         	                                    ! Put vector field in the annulus between N and N+vicinity 
	        ver(vneighlist(i,1))%czero_flag=1 
	        ver(Vneighlist(i,1))%czero=curzero
	        Endif
		i=i+1
        Enddo
	ann_no = ann_no+1
	Enddo
        DeAllocate(Vneighlist)
        End Subroutine Initialize_Annulus

!---------------------------------------------------------------------------------------------------------------------------
!                              Initialize  vector field in an annulus
!---------------------------------------------------------------------------------------------------------------------------
        Subroutine Initialize_Annulus_conc(nver_curzero,ver_no,ann_begin,curzero)                                                       ! Makes a patch with scalar field 1 and czero_flag field 1 
        Use module_datastruct 
        Implicit None
        Integer :: i,ann_no,ncount
        Real(Kind=8) :: curzero
        Integer :: ver_no,nver_curzero,ann_begin

        Call Complete_neighbourhood(ver_no,'RETURN_LIST')
        ver(:)%czero_flag=-1 ;  ncount=0 ; ver(:)%czero=0.0
	ann_no=ann_begin

	Do while(ncount .Le. nver_curzero) 
	i=1
	Do while ((i.Le.nver).And.(ncount .Le. nver_curzero))
        	If(Vneighlist(i,2).Eq. ann_no )Then                                                                         	    ! Put vector field in the annulus between N and N+vicinity 
	        ver(vneighlist(i,1))%czero_flag=1 
	        ver(Vneighlist(i,1))%czero=curzero
		ncount=ncount+1
	        Endif
		i=i+1
        Enddo
	ann_no = ann_no+1
	Enddo
        DeAllocate(Vneighlist)
        End Subroutine Initialize_Annulus_conc
	

!---------------------------------------------------------------------------------------------------------------------------
!                             Compute the curvature of an Epsin in its neighbourhood after being moved
!---------------------------------------------------------------------------------------------------------------------------
        Subroutine Complete_neighbourhood(v,retflag)                    
        Use module_datastruct
        Implicit None
        Integer :: v
        Character(Len=9),Allocatable :: accounted(:)
        Character(Len=9), Intent(In), Optional :: retflag
        Integer :: i,j,ncount,ncount1
        Integer,Allocatable :: neilist(:,:)
        ver(:)%neigh=0

        Allocate(neilist(nver,2)) ; Allocate(accounted(nver))
        accounted='NOTACC' ; neilist=0
        neilist(1,1)=v ; neilist(1,2)=0                                  ! The first entry is the vertex itself  
        ncount=ver(v)%nonei+1
        neilist(2:ncount,1)=ver(v)%vneipt(1:ver(v)%nonei)                ! Make a list with vertex v as a starting point             
        neilist(2:ncount,2)=1                                            ! Set the number of the neighbourhood we are accounting for 
        Forall (i=1:ncount) accounted(neilist(i,1))='ACCOUNTED'          ! this maintains a list that avoids overcounting            
        Forall (i=2:ncount) ver(neilist(i,1))%neigh=neilist(i,2)         ! this maintains a list that avoids overcounting            

        ncount1=2
        Do While(ncount.Lt.nver)
           Do i=1,ver(neilist(ncount1,1))%nonei,1
              j=ver(neilist(ncount1,1))%vneipt(i)
              If(accounted(j).Eq.'NOTACC')Then
              ncount=ncount+1 
              neilist(ncount,1)=j
              neilist(ncount,2)=neilist(ncount1,2)+1
              ver(j)%neigh=neilist(ncount,2)
              accounted(j)='ACCOUNTED'
              Endif
           Enddo
           ncount1=ncount1+1
        Enddo

        If(Present(retflag)) Then                                        ! Return the list only if asked for 
        Allocate(Vneighlist(nver,2))        
        Vneighlist=neilist
        Endif

        DeAllocate(neilist) ; DeAllocate(accounted)
        End Subroutine Complete_neighbourhood      

      End Module module_Initialize_curvature_field


