        Module Module_WriteData
        Contains

!------------------------------------------------------------------------------------------------------------------------
!                                    SUBROUTINE TO STORE THE STATE OF THE SYSTEM
!------------------------------------------------------------------------------------------------------------------------
        Subroutine membrane_dump1(ensno)
        Use module_datastruct
        Integer :: ensno
        Character(100):: filename, ensname
        
        Write(ensname,*)ensno
		filename='../SYS_STATE/DUMP/membrane_state-'//Trim(Adjustl(ensname))//'.dump'
		Call membrane_dump(filename)
		
        End Subroutine membrane_dump1
!----------------------------------- ------------------------------------------------------------------------------------
!                                    SUBROUTINE TO STORE THE STATE OF THE SYSTEM
!------------------------------------------------------------------------------------------------------------------------
        Subroutine membrane_dump2(ensno,time)
        Use module_datastruct
        Integer :: ensno,time
        Character(100):: filename,tname,ensname
        
        Write(ensname,*)ensno
        Write(tname,*)time
        filename='../SYS_STATE/DUMP/membrane_state-ENS-'//Trim(Adjustl(ensname))//'-'//Trim(Adjustl(tname))//'.dump'
        Call membrane_dump(filename)
        
        End Subroutine membrane_dump2
!------------------------------------------------------------------------------------------------------------------------
!                                    SUBROUTINE TO STORE THE STATE OF THE SYSTEM
!------------------------------------------------------------------------------------------------------------------------
        Subroutine membrane_dump(filename)
        Use module_datastruct
        Integer :: i
        Character(100):: filename
        
        Open(01,File=Trim(Adjustl(filename)),Form='Formatted')
        Write(01,*)blen,blLcut,blUcut                      	    ! Write the bond length, upper cutoff and lowercutoff used to generate the given configuration
        Write(01,'(4(I8,1x))')nver,pbnum,ntr,tlink
        Write(01,*)'---BeginVertexData'
        Do i=1,pbnum,1
        Write(01,'(3(F14.6,1x))')ver(i)%vnor
        Write(01,'(3(F14.6,1x))')ver(i)%vcoord
        Write(01,'(10(I8,1x))')ver(i)%vneipt
        Write(01,'(10(I8,1x))')ver(i)%vneitr
        Write(01,'(10(I8,1x))')ver(i)%PBCver
        Write(01,'(9(F14.6,1x))')ver(i)%L2G(:,1),ver(i)%L2G(:,2),ver(i)%L2G(:,3)
        Write(01,'(9(F14.6,1x))')ver(i)%HHM(:,1),ver(i)%HHM(:,2),ver(i)%HHM(:,3)
        Write(01,'(9(I9,1x))')ver(i)%nonei,ver(i)%boundary,ver(i)%pbimno,ver(i)%pbmap,ver(i)%imver,ver(i)%nover,ver(i)%boxvert
        Write(01,'(4(F14.6,1x))')ver(i)%mcur,ver(i)%cur1,ver(i)%cur2,ver(i)%totarea
        Write(01,'(I4,1x,F14.6,1x)')ver(i)%czero_flag,ver(i)%czero
        Write(01,'(I8,1x)')ver(i)%cellno
        Enddo
        Write(01,*)'---EndofVertexData'
        
        Write(01,*)'---BeginAntigenData'
        Do i=1,num_antigens,1
        Write(01,*)antig(i)%vertex
        Write(01,'(6(F14.6,1x))')antig(i)%base_coord,antig(i)%tip_coord
        Write(01,*)antig(i)%diffus_tri,antig(i)%disp1,antig(i)%disp2
        Write(01,*)antig(i)%theta,antig(i)%phi
        Enddo
        Write(01,*)'---EndofAntigenData'
        
        Write(01,*)'---BeginTriangleData'
        Do i=1,ntr,1
        Write(01,'(F14.6,1x,I2,1x,I2,1x,6(I8,1x),3(F14.6,1x))') &
        tri(i)%ar,tri(i)%pbflag,tri(i)%boxtriangle,tri(i)%li,tri(i)%vert,tri(i)%fnor
        Enddo 
        Write(01,*)'---EndofTriangleData'
        
        Write(01,*)'---BeginLinkData'
        Do i=1,tlink
        Write(01,'(6(I8,1x))')lin(i)%tr,lin(-i)%tr,lin(i)%boundary,lin(i)%pbflag,lin(i)%sep
        Enddo
        Write(01,*)'---End of Link Data'
        Close(01)
        
        End Subroutine membrane_dump
                
!------------------------------------------------------------------------------------------------------------------------
!                                    SUBROUTINE TO STORE XYZ Coords of  MEMBRANE  SYSTEM
!------------------------------------------------------------------------------------------------------------------------
      Subroutine Write_membrane_xyz(rank,window,time)
       Use module_datastruct
        Integer :: time,i,rank,window
         Character(100):: filename,tname,wname,rname
         Write(tname,*)time
	 	 Write(wname,*)window
		 Write(rname,*)rank
        filename='../SYS_STATE/DUMP/memb_xyz-'//Trim(Adjustl(rname))//'-'//Trim(Adjustl(wname))//'-'//Trim(Adjustl(tname))//'.xyz'
         Open(01,File=Trim(Adjustl(filename)),Form='Formatted')
	 	 Write(01,*)nver,pbnum,periodic_box_length
         Do i=1,pbnum,1
         Write(01,'(4(F14.6,1x))')ver(i)%vcoord,ver(i)%czero
         Enddo
         Close(01)
         End Subroutine Write_membrane_xyz
	 
 !------------------------------------------------------------------------------------------------------------------------
 !                                    SUBROUTINE TO STORE THE AREA OF THE MEMBRANE
 !------------------------------------------------------------------------------------------------------------------------
      Subroutine Write_membrane_area(rank,window,time)
       Use module_datastruct
        Integer :: time,rank,window
		Character(50):: filename,tname,rankname,winname
        Write(tname,*)time
	    Write(rankname,*)rank
	    Write(winname,*)window
	    Call compute_projectedarea()
        filename='./Membrane_area_ENS-'//Trim(Adjustl(rankname))//'_Frame-'//Trim(Adjustl(winname))//'.dat'
        Open(01,File=filename,Form='Formatted',Position='Append')
		Write(01,*)time,Sum(ver(:)%totarea),mp%projarea
        Close(01)
      End Subroutine Write_membrane_area
!------------------------------------------------------------------------------------------------------------------------
!                                    SUBROUTINE TO WRITE THE FILES FOR VTK VIEWERS
!------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE vtkformat(time)
        Use module_datastruct
        INTEGER:: time,num,i,j,nonpertr
        Integer,Dimension(:),Allocatable:: foff,ftype 
        CHARACTER(1000) ::temfi,timename
        nonpertr=0
2000    FORMAT(A23,I5,A1,1X,A15,I5,A2)

        Call check_triangle_in_pbbox()
        Call check_vert_in_pbbox()

		Write(timename,*)time
		Write(temfi,*)'conf-'//Trim(Adjustl(timename))//'.vtu'

        OPEN(11,FILE=Trim(Adjustl(temfi)),FORM='FORMATTED')
        WRITE(11,*)'<VTKFile type="UnstructuredGrid" version="0.1"','  byte_order="BigEndian">'
        WRITE(11,*)'<UnstructuredGrid>'
        WRITE(11,2000)'<Piece NumberOfPoints="',pbnum,'"','NumberOfCells="',ntr,'">'
        WRITE(11,*)'<PointData Scalars="scalars">'

        WRITE(11,*)'<DataArray type="Float32" Name="fl" Format="ascii">'                                                            ! Mean curvature
        Do i=1,nver,1
        Write(11,*) 0
        Enddo
        Do i=nver+1,pbnum,1
        Write(11,*) 1
        Enddo
        WRITE(11,*)'</DataArray>'

		WRITE(11,*)'<DataArray type="Int32" Name="vinbox" Format="ascii">'                                                            ! Mean curvature
        Do i=1,pbnum,1
        Write(11,*) ver(i)%boxvert
        Enddo
        WRITE(11,*)'</DataArray>'


        WRITE(11,*)'<DataArray type="Float32" Name="H" Format="ascii">'                                                             ! Principal curvature 1
        Do i=1,nver,1
        Write(11,*) ver(i)%mcur
        Enddo
        Do i=nver+1,pbnum,1
        j=ver(i)%imver
        Write(11,*) ver(j)%mcur
        Enddo
        WRITE(11,*)'</DataArray>'

        WRITE(11,*)'<DataArray type="Int32" Name="Antigenflag" Format="ascii">'                                                     ! Principal curvature 1
        Do i=1,nver,1
        Write(11,*)ver(i)%antigen_flag
        Enddo
        Do i=nver+1,pbnum,1
        j=ver(i)%imver
        Write(11,*) ver(j)%antigen_flag
        Enddo
        WRITE(11,*)'</DataArray>'

		WRITE(11,*)'<DataArray type="Int32" Name="ncshadow" Format="ascii">'                                                     ! Principal curvature 1
        Do i=1,pbnum,1
        Write(11,*)ver(i)%shadownc
        Enddo
        WRITE(11,*)'</DataArray>'


        WRITE(11,*)'<DataArray type="Float32" Name="c1" Format="ascii">'                                                            ! Principal curvature 1
        Do i=1,nver,1
        Write(11,*) ver(i)%cur1
        Enddo
        Do i=nver+1,pbnum,1
        j=ver(i)%imver
        Write(11,*) ver(j)%cur1
        Enddo
        WRITE(11,*)'</DataArray>'

		WRITE(11,*)'<DataArray type="Float32" Name="czero" Format="ascii">'                                                            ! Principal curvature 1
        Do i=1,pbnum,1
        Write(11,*) ver(i)%czero
        Enddo
        WRITE(11,*)'</DataArray>'


		WRITE(11,*)'<DataArray type="Int32" Name="czero-flag" Format="ascii">'                                                            ! Principal curvature 1
        Do i=1,pbnum,1
        Write(11,*) ver(i)%czero_flag
        Enddo
        WRITE(11,*)'</DataArray>'

        WRITE(11,*)'<DataArray type="Float32" Name="c2" Format="ascii">'                                                            ! Principal curvature 2
        Do i=1,nver,1
        Write(11,*) ver(i)%cur2
        Enddo
        Do i=nver+1,pbnum,1
        j=ver(i)%imver
        Write(11,*) ver(j)%cur2
        Enddo
        WRITE(11,*)'</DataArray>'
        WRITE(11,*)'</PointData>'

		Write(11,*)'<CellData>'
		WRITE(11,*)'<DataArray type="Int32" Name="tinbox" Format="ascii">'                                                            ! Principal curvature 2
        Do i=1,ntr,1
        Write(11,*) tri(i)%boxtriangle
        Enddo
        WRITE(11,*)'</DataArray>'
		Write(11,*)'</CellData>'


        WRITE(11,*)'<Points>'
        WRITE(11,*)'<DataArray type="Float32"','  NumberOfComponents="3" Format="ascii">'
        Do i=1,pbnum,1
        Write(11,*)ver(i)%vcoord
        Enddo
        WRITE(11,*)'</DataArray>'
        WRITE(11,*)'</Points>'

        WRITE(11,*)'<Cells>'
        WRITE(11,*)'<DataArray type="Int32"','  Name="connectivity" Format="ascii">'
        Do i=1,ntr,1
        Write(11,*)tri(i)%vert-1
        Enddo
        WRITE(11,*)'</DataArray>'

        Allocate(foff(ntr)); Allocate(ftype(ntr))
        num=0
        Do i=1,ntr,1
        num=num+1
        foff(num)=3*num
        Enddo
        ftype=5

        WRITE(11,*)'<DataArray type="Int32" Name="offsets" Format="ascii">' 
        Write(11,*)foff
        Write(11,*)'</DataArray>'

        WRITE(11,*)'<DataArray type="Int32" Name="types" Format="ascii">' 
        Write(11,*)ftype
        Write(11,*)'</DataArray>'

        WRITE(11,*)'</Cells>'

        WRITE(11,*)'</Piece>'
        WRITE(11,*)'</UnstructuredGrid>'
        WRITE(11,*)'</VTKFile>'

        Close(11)
        End Subroutine vtkformat

!------------------------------------------------------------------------------------------------------------------------
!                                    SUBROUTINE TO WRITE THE FILES FOR VTK VIEWERS
!------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE Create_periodicbox()
        Use module_datastruct
        INTEGER::num,i
        Integer,Dimension(:),Allocatable:: foff,ftype 
        CHARACTER(1000) :: filename

        filename='PeriodicBox.vtu'         
        OPEN(11,FILE=filename,FORM='FORMATTED')
        WRITE(11,*)'<VTKFile type="UnstructuredGrid" version="0.1"', '  byte_order="BigEndian">'
        WRITE(11,*)'<UnstructuredGrid>'
        WRITE(11,*)'<Piece NumberOfPoints="8"   NumberOfCells="6">'

        WRITE(11,*)'<Points>'
        WRITE(11,*)'<DataArray type="Float32"','  NumberOfComponents="3" Format="ascii">'
        Write(11,*)zero,zero,zero
        Write(11,*)periodic_box_length,zero,zero
        Write(11,*)periodic_box_length,periodic_box_length,zero
        Write(11,*)zero,periodic_box_length,zero
        Write(11,*)zero,zero,periodic_box_height
        Write(11,*)periodic_box_length,zero,periodic_box_height
        Write(11,*)periodic_box_length,periodic_box_length, periodic_box_height
        Write(11,*)zero,periodic_box_length,periodic_box_height
        WRITE(11,*)'</DataArray>'
        WRITE(11,*)'</Points>'

        WRITE(11,*)'<Cells>'
        WRITE(11,*)'<DataArray type="Int32"','  Name="connectivity" Format="ascii">'
        Write(11,*)'0  1  2  3'
        Write(11,*)'0  1  5  4'
        Write(11,*)'0  3  7  4'
        Write(11,*)'4  5  6  7'
        Write(11,*)'2  3  7  6'
        Write(11,*)'2  6  5  1'
        WRITE(11,*)'</DataArray>'

        Allocate(foff(6)); Allocate(ftype(6))
        num=0
        Do i=1,6,1
        num=num+1
        foff(num)=4*num
        Enddo
        ftype=5

        WRITE(11,*)'<DataArray type="Int32" Name="offsets"  Format="ascii">' 
        Write(11,*)foff
        Write(11,*)'</DataArray>'

        WRITE(11,*)'<DataArray type="Int32" Name="types"  Format="ascii">' 
        Write(11,*)ftype
        Write(11,*)'</DataArray>'

        WRITE(11,*)'</Cells>'

        WRITE(11,*)'</Piece>'
        WRITE(11,*)'</UnstructuredGrid>'
        WRITE(11,*)'</VTKFile>'

        Close(11)
        End Subroutine Create_periodicbox

!===========================================================================================================================
!                  Subroutine to calculate the area of the given face and update the total area linked to the vertices
!===========================================================================================================================
       Subroutine compute_projectedarea()
       USE module_datastruct
       IMPLICIT NONE
       Integer::i,j,k,tr
       Real(KIND=8),DIMENSION(3,1)::r1,r2,r3,r21,r31
       mp%projarea = 0.0
       Do tr=1,ntr,1
        i=tri(tr)%vert(1) ; j=tri(tr)%vert(2) ; k=tri(tr)%vert(3)
        r1=ver(i)%vcoord ; r2=ver(j)%vcoord ; r3=ver(k)%vcoord                                                                       ! The coordinates of the 3 vertices
        r21=r2-r1 ; r31=r3-r1                                                                                                        ! Relative position vectors for area calc
        mp%projarea = mp%projarea + (r21(1,1)*r31(2,1)-r21(2,1)*r31(1,1))*0.5
       Enddo
       End Subroutine compute_projectedarea

        

        End Module Module_WriteData
