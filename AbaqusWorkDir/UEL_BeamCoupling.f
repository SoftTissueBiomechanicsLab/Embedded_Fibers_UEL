c Beam to Solid Coupling Element (Compatible with 8 node Brick Elements)
c By: Sotirios Kakaletsis
c October 12th, 2021
c
c References
c 1) Steinbrecher, Ivo, et al. "A mortar-type finite element approach 
c for embedding 1D beams into 3D solid volumes." Computational Mechanics 
c  66.6 (2020): 1377-1398.
c =====================================================================
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
c
      INCLUDE 'ABA_PARAM.INC'
c
      INTEGER :: nTA, nMA, nTB, nMB, nTC, nMC
      PARAMETER (nTA=3,nMA=3,nTB=3,nMB=3)
      PARAMETER (nTC=3,nMC=3)
      DIMENSION Ag(nTA,nMA), Bg(nTB,nMB), Cg(nTC,nMC)
      DOUBLE PRECISION Ag, Bg, Cg
      CHARACTER (len = 120) :: ADir, BDir, CDir
      COMMON /KMortar/ Ag, Bg, Cg 
c
      IF(LOP==0)THEN
c     Load mortar matrices	  
        ADir = '/home/datastore/Sotiris/Mglobal.txt'
        BDir = '/home/datastore/Sotiris/Mglobal.txt'
        CDir = '/home/datastore/Sotiris/Mglobal.txt'
        write(*,*) 'START ANALYSIS, Reading Mortar Matrices...'
        call ReadMatrix(1, ADir, nTA, nMA, Ag)
        call ReadMatrix(2, BDir, nTB, nMB, Bg)
        call ReadMatrix(3, CDir, nTC, nMC, Cg)
        write(*,*) 'Mortar Matrices import: SUCCESS'
      ELSEIF(LOP == 3)THEN
	write(*,*) 'FINISHED ANALYISIS'
      ENDIF
      RETURN
      END
c =====================================================================
      SUBROUTINE UEL (RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     1PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     2DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     3PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS,
     4NJPROP, PERIOD)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
c Inherit common block for mortar matrices=============================
      INTEGER :: nTA, nMA, nTB, nMB, nTC, nMC, nG
      PARAMETER (nTA=3,nMA=3,nTB=3,nMB=3)
      PARAMETER (nTC=3,nMC=3)
      DIMENSION Ag(nTA,nMA), Bg(nTB,nMB), Cg(nTC,nMC)
      DOUBLE PRECISION Ag, Bg, Cg
      COMMON /KMortar/ Ag, Bg, Cg
c Internal Important Variables ========================================
      INTEGER :: elID
      REAL:: e_pen
      DIMENSION Al(3*nG, 3*nMA), Bl(3*nG, 3*nMB), Cl(3*nG, 3*nMC)
c	  
      DOUBLE PRECISION Al, Bl, Cl
c
c~~~~~~~~~~~~~~A PART~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      if (Jtype .eq. 1) then 
c
c Define Inputs========================================================
c
      e_pen=PROPS(1)   
      elID=PROPS(2)	  
c
c Initialize Matrices and Vectors======================================
c
      call k_matrix_zero(rhs,ndofel,nrhs)
      call k_matrix_zero(amatrx,ndofel,ndofel)
c
c Recover Local Coupling Matrices =====================================
c
      call k_matrix_zero(Al, 3*nG, 3*nMA)
c
      call AssembleA(Ag, Al, elID, nTA, nG, nMA)
c	  write(*,*) 'A matrix', A
c
c AMATRIX
      do ia=1, 3*nG
        do ib=1, 3*nMA
         amatrx(ia,ib) = amatrx(ia,ib) + e_pen*Al(ia,ib)
        end do
      end do
c	   
c RHS
c
      do ia = 1, 3*nG
	 do ib = 1, 3*nMA
	  rhs(ia,1)=rhs(ia,1)-e_pen*Al(ia,ib)*U(ib)
	 end do
      end do
c
c~~~~~~~~~~~~~~B PART~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      elseif (Jtype .eq. 2) then 
c Define Inputs========================================================
c
      e_pen=PROPS(1)   
      elID=PROPS(2)	  
c
c Initialize Matrices and Vectors======================================
c
      call k_matrix_zero(rhs,ndofel,nrhs)
      call k_matrix_zero(amatrx,ndofel,ndofel)
c
c Recover Local Coupling Matrices =====================================
c
      call k_matrix_zero(Bl, 3*nG, 3*nMB)
c
      call AssembleA(Bg, Bl, elID, nTB, nG, nMB)
c	  write(*,*) 'A matrix', A
c
c AMATRIX
      do ia=1, 3*nG
        do ib=1, 3*nMB
         amatrx(ia,ib) = amatrx(ia,ib) + e_pen*Bl(ia,ib)
        end do
      end do
c	   
c RHS
c
      do ia = 1, 3*nG
	   do ib = 1, 3*nMB
	    rhs(ia,1)=rhs(ia,1)-e_pen*Bl(ia,ib)*U(ib)
	   end do
      end do
c
      elseif (Jtype .eq. 3) then 
c Define Inputs========================================================
c
      e_pen=PROPS(1)   
      elID=PROPS(2)	  
c
c Initialize Matrices and Vectors======================================
c
      call k_matrix_zero(rhs,ndofel,nrhs)
      call k_matrix_zero(amatrx,ndofel,ndofel)
c
c Recover Local Coupling Matrices =====================================
c
      call k_matrix_zero(Cl, 3*nG, 3*nMC)
c
      call AssembleA(Cg, Cl, elID, nTC, nG, nMC)
c	  write(*,*) 'A matrix', A
c
c AMATRIX
      do ia=1, 3*nG
        do ib=1, 3*nMC
         amatrx(ia,ib) = amatrx(ia,ib) + e_pen*Cl(ia,ib)
        end do
      end do
c	   
c RHS
c
      do ia = 1, 3*nG
	   do ib = 1, 3*nMC
	    rhs(ia,1)=rhs(ia,1)-e_pen*Cl(ia,ib)*U(ib)
	   end do
      end do
c
      endif
      return
      end
c======================================================================
      subroutine k_matrix_zero(A,n,m)
      INCLUDE 'ABA_PARAM.INC'
      dimension A(n,m)
c
      do i = 1, n
         do j = 1, m
            A(i,j)=0.d0
         end do
      end do
c
      return
      end
c======================================================================
      subroutine k_vector_zero(A,n)
      INCLUDE 'ABA_PARAM.INC'
      dimension A(n)
c
      do i = 1, n
         A(i)=0.d0
      end do
c
      return
      end
c======================================================================    
      subroutine ReadMatrix(IDx, MatDir, n, m, A)
      INCLUDE 'ABA_PARAM.INC'
c
      DIMENSION A(n,m)
c
      INTEGER IDx, n, m, f
      DOUBLE PRECISION A
      CHARACTER (len = 120) :: MatDir
c
      IF(m>3)THEN
	N4rows = m/4
	N4mod = MOD(m,4)
	f = 4
      ELSE
	N4rows = 1
	N4mod = 0
	f = m
      ENDIF
c	  
      open (unit=IDx, file=MatDir, status='old', action='read')	  
c
      I = 1
      do while (I<=n)
	do j=1, N4rows
	   col1 = (j-1)*f + 1
	   col2 = col1 +(f-1)
	   read(IDx,*), A(I,col1:col2)
         end do
         if (N4mod>0) then
	   col1 = N4rows*f + 1
	   col2 = col1 + N4mod - 1
	   read(IDx,*), A(I,col1:col2)
	 endif
	 I = I + 1
        end do
      close(IDx)
c
      return
      end
c
c======================================================================
      subroutine AssembleA(Ag, A, elID, nT, nG, nM)
      INCLUDE 'ABA_PARAM.INC'
c
      DIMENSION Ag(nT,nM), A(3*nG,3*nM)
c
      INTEGER elID, nT, nG, nM, counter, col_id
      DOUBLE PRECISION Ag, A, Ajk
c	 
      do i=1,nG
       counter = 0
       row_id = (i-1)*3+1
       g_row=(elID-1)*nG+i
       do j = 1,nM    
	 counter = counter+1
	 Ajk = Ag(g_row, counter)	
	 col_id = (j-1)*3+1
	 do ii=0,2
	   A(row_id+ii,col_id+ii)=Ajk
	 end do		
       end do	
      end do  
c
      return
      end	  
c=================================END==================================