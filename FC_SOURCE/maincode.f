!===========================================================================================================================
!                                  $mcode                        MAIN PROGRAM
!===========================================================================================================================

        PROGRAM triangulation
        USE module_makesurface ; USE module_curvcalc 
        USE module_rerun ; USE module_datastruct 
        USE module_mcsmoves ; USE module_dataformat  
        IMPLICIT NONE  
        Integer::rand,i,j,k,dataint
        Real(KIND=8)::ran2,tlen,tstart,tstop                             ! sarea --> Surface area of the surface
        Character (Len=10) :: sys_geom
        fangle=-0.5 ;  maxtl=SQRT(3.0)                                   ! gsize --> size of the grid  
        pi=acos(-1.0000000000000)                                    

        Open(10,FILE='parameters.in')
        Read(10,*),kappa
        Read(10,*),gsize
        Read(10,*),sys_geom
        Read(10,*),depth
        Read(10,*),period
        Read(10,*),dataint
        Close(10)


        If(sys_geom .Eq. 'PLANAR')Then
        Call makeplane()
        Else If (sys_geom .Eq. 'SINUSOIDAL')Then
        Call makesinusoidal()
        Endif

        Allocate(mp%vertex(pbnum)) ; Allocate(mp%triangle(ntr))
        Allocate(mp%link(-tlink:tlink))

        ver(:)%cur1=0.0 ; ver(:)%cur2=0.0 ; ver(:)%mcur=0.0

         vertice_call: Do i=1,nver                                       ! Normal calculation over each vertex 
         call normalcalc(i)
         Enddo vertice_call

         Do i=1,nver
         If(ver(i)%boundary.Eq.0)Then
         Do j=1,ver(i)%nonei
         k=ver(i)%vneipt(j)
         tlen=Sqrt(sum((ver(i)%vcoord-ver(k)%vcoord)**2)) 
         If(tlen.Le.1.0 .or.tlen.Ge.sqrt(3.0))Then
         Print*,i,k,tlen
         Endif
         Enddo
         Endif
         Enddo

         CALL vtkformat(0)

          Do i=nver+1,pbnum,1                                            ! Initialize the periodic images 
          call mapimagetovertex(i)
          Enddo

         mcs_loop: Do mcs=itime,ftime,1                                  ! Monte Carlo loop 
           inner_loop: Do iloop=1,nver,1                                 ! Inner loop       

           rand=nint((1-2*ran2(seed))*tlink)                             ! Flipping of the links
           If (abs(rand).GT.0.0 .And. lin(rand)%boundary.EQ.0)Then
           call flipping(rand)
           Endif 

           rand=nint(ran2(seed)*nver)+1                                  ! Vertex Move  
           If (rand.LE.nver) Then
           Call movevertex(rand)
           Endif
           Enddo inner_loop

            If(mod(mcs,dataint).Eq.0)Then
            Call vtkformat(mcs/dataint)
            Endif

             Enddo mcs_loop
        End PROGRAM triangulation

!---------------------------------------------------------------------------------------------------------------------------
!                           Subroutine TO CALCULATE ALL ANALYTIC  QUANTITIES
!---------------------------------------------------------------------------------------------------------------------------
         Subroutine analqtys()  
         USE module_curvcalc ; USE module_datastruct 
         USE module_mcsmoves ; USE module_dataformat
         IMPLICIT NONE
         Integer :: i,j,k,v,ilen,ichk
         Real(KIND=8) :: vol,cvol                                       
         Real(KIND=8),DIMENSION(3,1) :: a,b,c    
         Real(KIND=8) :: xcm,ycm,zcm,elen,spspen,bl                      ! COM, Elastic,nematic-nematic,nematic-curvature 
         Real(KIND=8) :: r1(3,1),r2(3,1)  
         
!          mp%rg=0 ; mp%area=0 ; ilen=0 ; cvol=0; vol=0
!          xcm=SUM(ver(1:nver)%vcoord(1,1))/nver                          ! X component of the center of mass
!          ycm=SUM(ver(1:nver)%vcoord(2,1))/nver                          ! Y component                      
!          zcm=SUM(ver(1:nver)%vcoord(3,1))/nver                          ! Z component                      
!
!         mp%rg=(SUM((ver(1:nver)%vcoord(1,1)-xcm)**2)+                   ! Rg square calculation 
!     $          SUM((ver(1:nver)%vcoord(2,1)-ycm)**2)+
!     $          SUM((ver(1:nver)%vcoord(3,1)-zcm)**2))/nver
!
!         mp%area=SUM(tri(1:ntr)%ar)                                      ! Surface area of the membrane
!
!            vol_calc:Do i=1,ntr
!            a=0; b=0; c=0
!            a=ver(tri(i)%vert(1))%vcoord  
!            b=ver(tri(i)%vert(2))%vcoord  
!            c=ver(tri(i)%vert(3))%vcoord
!
!            cvol=(a(1,1)*(b(2,1)*c(3,1)-b(3,1)*c(2,1))+                   ! (a.(bxc))/6.0 is the volume 
!     $            a(2,1)*(b(3,1)*c(1,1)-b(1,1)*c(3,1))+
!     $            a(3,1)*(b(1,1)*c(2,1)-b(2,1)*c(1,1)))/6.0
!            vol=vol+cvol 
!            Enddo vol_calc
!            mp%vol=vol

            elen=0                                                       ! Energy calculation
            Do i=1,nver
            elen=elen+kappa*(ver(i)%mcur)**2*ver(i)%totarea
            Enddo
            mp%enel=elen

       End Subroutine analqtys

        


!-------------------------------------------------------------------------------------------------------------------------
!                                   Function TO CHECK THE MAX And MIN BONDLENGTH
!-------------------------------------------------------------------------------------------------------------------------
      Function blcheck(ver1,ver2) result(fres) 
      USE module_datastruct ; USE module_makesurface
      USE module_mcsmoves
      IMPLICIT NONE
      Real(KIND=8)::bl
      Real(KIND=8),DIMENSION(3,1)::co1,co2,co3
      Integer::ver1,ver2,fres                                            !fres=1 ===>false ,0===>true
      fres=0 
      co1=ver(ver1)%vcoord
      co2=ver(ver2)%vcoord
      co3=co1-co2  
      bl=sqrt(co3(1,1)**2+co3(2,1)**2+co3(3,1)**2)
      If (bl.GT.SQRT(3.000).OR. bl.LE.1.00000000) fres=1            
      RETURN
      End Function blcheck

!-------------------------------------------------------------------------------------------------------------------------
!                                   Function TO CHECK THE MINIMUM BONDLENGTH
!-------------------------------------------------------------------------------------------------------------------------
        Function blcheckmin(ver1,ver2) result(fres) 
        USE module_datastruct ; USE module_makesurface
        USE module_mcsmoves
        IMPLICIT NONE
        Real(KIND=8)::bl
        Real(KIND=8),DIMENSION(3,1)::co1,co2,co3                 
        Integer::ver1,ver2,fres                                           ! fres=1 ===>false ,0===>true
        fres=0 
        co1=ver(ver1)%vcoord ; co2=ver(ver2)%vcoord
        co3=co1-co2  
        bl=SQRT(SUM(co3**2))
        If(bl.LE.1.0000000000) fres=1            
        RETURN 
        End Function blcheckmin

!------------------------------------------------------------------------------------------------------------------------
!                                  Function to generate a random number
!------------------------------------------------------------------------------------------------------------------------  
      FUNCTION ran2(idum)
      IMPLICIT NONE
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real*8  ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1 
     *   ,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791
     *   ,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/

      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue  
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END  FUNCTION ran2
!---------------------------------------------------------------------------------------------------------------------------



















