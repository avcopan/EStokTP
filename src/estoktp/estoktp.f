c EStokTP
c Electronic Structure to k(T,P) 
c This is a program to generate rate constants relying on calls to external codes to perform
c electronic structure calculations and master equation simulations
c Copyright (C) 2018  by Carlo Cavallotti and Matteo Pelucchi, Politecnico di MIlano, Italy
c and Stephen J. Klippenstein, Argonne National Laboratories, USA.
c This program is free software: you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation, version 3 of the License
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.

c The procedure is as follows:
c (i) tau optimization of reactants: molecule and radical; and products
c (ii) constrained optimization on grid to locate approximate ts
c (iii) ts from maximum on grid
c (iv) improved ts from prelim ts
c (v) tau optimization of ts 
c (vi) higher level E's at ts, reactants, and products
c (vii) 1D sampling of tau for ts, reactants, and products from tau_opt 
c (viii) IRC from ts for tau_opt
c (ix) preparation of me input
c (x) run me to produce rate
c (xi) extract me output and fit to modified Arrhenius form
c For now all of these procedures are fixed to some very specific appraoch
c Ultimately should allow various branchings to different methodologies
c (i) employs method_reac_opt
c (ii) employs method_grid_opt
c (iii) employs method_ts_0
c (iv) employs method_ts_1
c (v) employs method_ts_tauo
c (vi) employs method_ts_hl
c (vii) employs method_1dtau
c (viii) employs method_irc
c     
c data for main program consists of 
c (i) reactant specifications in internal coordinates
c (ii) logical flags for each part of the code 
c (iii) temperatures and pressures for calculation
c then there is separate data for each subcode consisting of the methods to be used etc.
c complications to be dealt with later
c (a) must specify which vibrational modes to remove
c (b) for pressure dependence need data for complex - e.g. vdw for abstraction, adduct for addition
c
c Immediate changes to make - 
c (ii) continue on to optimized torsions in ts
c	for this will need to have some specificiation of additional torsional degrees of freedom

      program estoktp

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      LOGICAL leof,lsec,ltit
      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7
      character*100 command1
      character*30 gmemlow,gmemhigh
      character*30 gmem

      include 'filcomm.f'

cc initialize directory structure

      command1='mkdir -p output'
      call commrun(command1)
      command1='mkdir -p geoms'
      call commrun(command1)
      command1='mkdir -p irc_files'
      call commrun(command1)
      command1='mkdir -p me_files'
      call commrun(command1)
      command1='mkdir -p me_blocks'
      call commrun(command1)
      command1='mkdir -p hl_logs'
      call commrun(command1)
      command1='mkdir -p hr_geoms'
      call commrun(command1)
      command1='mkdir -p md_tunn'
      call commrun(command1)

c random number generator

      call mtinit()

c initializations

      iopt_reac1=0
      iopt_reac2=0
      iopt_prod1=0
      iopt_prod2=0
      iopt_wellr=0
      igeom_wellr=0
      iopt_wellp=0
      igeom_wellp=0
      i2dgrid_opt_ts=0
      igrid_opt_ts=0
      iopt_ts_0=0
      itauo_ts=0
      iopt_reac1_1=0
      iopt_reac2_1=0
      iopt_prod1_1=0
      iopt_prod2_1=0
      iopt_wellr_1=0
      iopt_wellp_1=0
      iopt_ts_1=0
      ibmat_reac1=0
      ibmat_reac2=0
      ibmat_prod1=0
      ibmat_prod2=0
      ibmat_wellr=0
      ibmat_wellp=0
      ibmat_ts=0
      ibmat_test=0
      i1dtau_reac1=0
      i1dtau_reac2=0
      i1dtau_prod1=0
      i1dtau_prod2=0
      i1dtau_wellr=0
      i1dtau_wellp=0
      i1dtau_ts=0
      imdtau_reac1=0
      imdtau_reac2=0
      imdtau_prod1=0
      imdtau_prod2=0
      imdtau_wellr=0
      imdtau_wellp=0
      imdtau_ts=0
      isymm_reac1=0
      isymm_reac2=0
      isymm_prod1=0
      isymm_prod2=0
      isymm_wellr=0
      isymm_wellp=0
      isymm_ts=0
      ihl_reac1=0
      ihl_reac2=0
      ihl_reacs=0
      ihl_prod1=0
      ihl_prod2=0
      ihl_prods=0
      ihl_wellr=0
      ihl_wellp=0
      ihl_ts=0
      iirc=0
      iktp=0
      imodarr=0
      idebug=0
      iabs=0
      iadd=0
      ivar=0
      iresirc=0
      ifixdis=0
      imdtunn=0
      iprojrcoo=1
      inotunn=0
      irw=0
      ipw=0
      ip1=0
      ip2=0
      nts=1
      ifrozrts = 0
      irecov = 0
      itotcalc=0
      intfreq=0

      open (unit=5,file='./data/estoktp.dat',status='UNKNOWN')
      open (unit=6,file='./output/estoktp.out',status='UNKNOWN')
      open (unit=7,file='./output/warnings.out',status='UNKNOWN')

  100 continue
      CALL LineRead (0)
      CALL LineRead (5)
      if (WORD.EQ.'REACTIONTYPE') then
         if (WORD2.EQ.'ABSTRACTION') iabs = 1
         if (WORD2.EQ.'ADDITION') iadd = 1
         if (WORD2.EQ.'ISOMERIZATION') iiso = 1
         if (WORD2.EQ.'BETASCISSION') ibeta = 1
         if (WORD3.EQ.'1TS') nts=1
         if (WORD3.EQ.'2TS') nts=2
         if (WORD3.EQ.'3TS') nts=3
      endif
      if (WORD.EQ.'VARIATIONAL') then
         ivar = 1
         if(WORD2.EQ.'INTERNAL')intfreq=1
      endif
      if (WORD.EQ.'RESIRC') then
         iresirc = 1
      endif
      if (WORD.EQ.'RESIRCALL') then
         iresirc = 2
      endif
      if (WORD.EQ.'FROZRTS') then
         ifrozrts = 1
         open (unit=99,status='unknown') 
         write (99,*) WORD2
         REWIND (99)
         read (99,*) frozcoo
         close (99)
c         write(*,*)'fr coo is ',frozcoo
c         stop
      endif
      if (WORD.EQ.'RECOVER') then
         irecov = 1
      endif
      if (WORD.EQ.'MDTUNNEL') then
         imdtunn = 1
         iprojrcoo = 0
         if(WORD2.EQ.'INTERNAL')intfreq=1
      endif
      if (WORD.EQ.'PROJRCOO') then
         iprojrcoo = 0
      endif

      if (WORD.EQ.'NOTUNNEL') then
         inotunn = 1
      endif
      if (WORD.EQ.'PRODS') then
         ip1 = 1
         ip2 = 1
      endif
      if (WORD.EQ.'PROD') then
         ipr1 = 1
      endif
      if (WORD.EQ.'WELLR') then
         irw = 1
         if (WORD2.EQ.'FINDGEOM'.AND.WORD3.EQ.'LEVEL0') then
            igeom_wellr=1
         endif
         if (WORD2.EQ.'FINDGEOM'.AND.WORD3.EQ.'LEVEL1') then
            igeom_wellr=2
         endif
      endif
      if (WORD.EQ.'WELLP') then
         ipw = 1
         if (WORD2.EQ.'FINDGEOM'.AND.WORD3.EQ.'LEVEL0') then
            igeom_wellp=1
         endif
         if (WORD2.EQ.'FINDGEOM'.AND.WORD3.EQ.'LEVEL1') then
            igeom_wellp=2
         endif
      endif
c      if (WORD.EQ.'STOICHIOMETRY') then
c         stoich=word2
c      endif
      if (WORD.EQ.'DEBUG') then
         open (unit=99,status='unknown') 
         REWIND (99)
         write (99,*) WORD2
         REWIND (99)
         read (99,*) idebug
         close (99)
         if (idebug.ge.1) write (6,*) 'idebug test',idebug
      endif
      if (WORD.EQ.'OPT_REAC1') then
         iopt_reac1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'OPT_REAC2') then
         iopt_reac2=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'OPT_PROD1') then
         iopt_prod1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'OPT_PROD2') then
         iopt_prod2=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'OPT_WELLR') then
         iopt_wellr=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'OPT_WELLP') then
         iopt_wellp=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'GRID_OPT_TS') then
         igrid_opt_ts=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'2DGRID_OPT_TS') then
         i2dgrid_opt_ts=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'OPT_TS_0') then
         iopt_ts_0=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'TAUO_TS') then
         itauo_ts=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'OPT_REAC1_1') then
         iopt_reac1_1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'OPT_REAC2_1') then
         iopt_reac2_1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'OPT_PROD1_1') then
         iopt_prod1_1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'OPT_PROD2_1') then
         iopt_prod2_1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'OPT_WELLR_1') then
         iopt_wellr_1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'OPT_WELLP_1') then
         iopt_wellp_1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'OPT_TS_1') then
         iopt_ts_1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'BMAT_REAC1') then
         ibmat_reac1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'BMAT_REAC2') then
         ibmat_reac2=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'BMAT_PROD1') then
         ibmat_prod1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'BMAT_PROD2') then
         ibmat_prod2=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'BMAT_WELLR') then
         ibmat_wellr=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'BMAT_WELLP') then
         ibmat_wellp=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'BMAT_TS') then
         ibmat_ts=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'BMAT_TEST') then
         ibmat_test=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'1DTAU_REAC1') then
         i1dtau_reac1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'1DTAU_REAC2') then
         i1dtau_reac2=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'1DTAU_PROD1') then
         i1dtau_prod1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'1DTAU_PROD2') then
         i1dtau_prod2=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'1DTAU_WELLR') then
         i1dtau_wellr=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'1DTAU_WELLP') then
         i1dtau_wellp=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'1DTAU_TS') then
         i1dtau_ts=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'MDTAU_REAC1') then
         imdtau_reac1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'MDTAU_REAC2') then
         imdtau_reac2=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'MDTAU_PROD1') then
         imdtau_prod1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'MDTAU_PROD2') then
         imdtau_prod2=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'MDTAU_WELLR') then
         imdtau_wellr=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'MDTAU_WELLP') then
         imdtau_wellp=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'MDTAU_TS') then
         imdtau_ts=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'SYMM_REAC1') then
         isymm_reac1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'SYMM_REAC2') then
         isymm_reac2=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'SYMM_PROD1') then
         isymm_prod1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'SYMM_PROD2') then
         isymm_prod2=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'SYMM_WELLR') then
         isymm_wellr=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'SYMM_WELLP') then
         isymm_wellp=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'SYMM_TS') then
         isymm_ts=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'IRC') then
         iirc=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'HL_REAC1') then
         ihl_reac1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'HL_REAC2') then
         ihl_reac2=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'HL_REACS') then
         ihl_reacs=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'HL_PROD1') then
         ihl_prod1=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'HL_PROD2') then
         ihl_prod2=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'HL_PRODS') then
         ihl_prods=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'HL_WELLR') then
         ihl_wellr=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'HL_WELLP') then
         ihl_wellp=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'HL_TS') then
         ihl_ts=1
         itotcalc=itotcalc+1
      endif
      if (WORD.EQ.'KTP') then
         iktp=1
      endif
      if (WORD.EQ.'MODARR') then
         imodarr=1
      endif
      if (WORD.EQ.'END') go to 200
      goto 100
  200 continue

      if(itotcalc.gt.1.and.irecov.eq.1)then
         write(*,*)'the recover options works only'
         write(*,*)'for just one calculation at a time'
         write(*,*)'please modify the input accordingly and restart'
         stop
      endif

c how many processors to use for low level and high level calculations
      read (5,*) numprocll,numprochl
      read (5,*)
c memory for low level and high level calculations
      call LineRead(5)

      gmemlow=word
      gmemhigh=word2

cc intro closure of unit 5
      close (unit=5,status='keep')

c ispecies = 0, 1, 2, 3, 4, 5, 6, 11, 12 => ts, reac1, reac2, prod1, prod2, wellr, wellp, reactants, products
c determine torsional minimum for reactants

c first part of the calculations, low level of nodes and memory      
      numproc=numprocll
      open (unit=99,status='unknown') 
      REWIND (99)
      write (99,*) gmemlow
      REWIND (99)
      read (99,*) gmem
      close (99)
c      gmem=gmemlow

      if (iopt_reac1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting tau_opt_reac1'
         ispecies=1
         call tau_opt(ispecies)
      endif
      if (iopt_reac2.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting tau_opt_reac2'
         ispecies=2
         call tau_opt(ispecies)
      endif
c determine torsional minimum for products
      if (iopt_prod1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting tau_opt_prod1'
         ispecies=3
         call tau_opt(ispecies)
      endif
      if (iopt_prod2.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting tau_opt_prod2'
         ispecies=4
         call tau_opt(ispecies)
      endif
c start ts search with constrained coordinate grid search 
      if (igrid_opt_ts.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting grid_opt'
         if(iabs.eq.1.or.iadd.eq.1) then
            call grid_opt
         else if (iiso.eq.1) then
            call grid_opt_iso   
         else if (ibeta.eq.1) then
            call grid_opt_beta   
         else
            write(*,*) 'supported reacion types for grid search'
            write(*,*) 'abstraction'
            write(*,*) 'addition'
            write(*,*) 'isomerization'
            write(*,*) 'beta-scission'
            write(*,*) 'make an adequate selection and restart the code'
            stop
         endif
      endif
c perform can of TS on 2D grid 
      if (i2dgrid_opt_ts.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting 2d grid scan'
         call grid_opt_2d
      endif
c proceed to optimization from maximum on grid search
      if (iopt_ts_0.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting ts_0'
         call ts_0
      endif
c determine torsional minimum for ts
      if (itauo_ts.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting tauo_ts'
         call tauo_ts
      endif
c determine torsional minimum for wells
      if (iopt_wellr.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting tau_opt_wellr'
         ispecies=5
         if(igeom_wellr.eq.1) ispecies=51
         call tau_opt(ispecies)
      endif
      if (iopt_wellp.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting tau_opt_wellp'
         ispecies=6
         if(igeom_wellp.eq.1) ispecies=61
         call tau_opt(ispecies)
      endif
c reoptimize at higher level
      if (iopt_reac1_1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting level1_reac1'
         ispecies=1
         call level1(ispecies)
      endif
      if (iopt_reac2_1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting level1_reac2'
         ispecies=2
         call level1(ispecies)
      endif
      if (iopt_prod1_1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting level1_prod1'
         ispecies=3
         call level1(ispecies)
      endif
      if (iopt_prod2_1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting level1_prod2'
         ispecies=4
         call level1(ispecies)
      endif
      if (iopt_ts_1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting level1_ts'
         ispecies=0
         call level1(ispecies)
      endif
      if (iopt_wellr_1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting level1_wellr'
         ispecies=5
         if(igeom_wellr.eq.2.or.igeom_wellr.eq.1) ispecies=51
         call level1(ispecies)
      endif
      if (iopt_wellp_1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting level1_wellp'
         ispecies=6
         if(igeom_wellp.eq.2.or.igeom_wellp.eq.1) ispecies=61
         call level1(ispecies)
      endif

c determine B matrix for slected species
      if (ibmat_reac1.eq.1) then
         if (idebug.ge.1) write (6,*) 'computing reac1 bmat'
         ispecies=1
         call bmatrix(ispecies,0)
      endif
      if (ibmat_reac2.eq.1) then
         if (idebug.ge.1) write (6,*) 'comuting reac2 bmat'
         ispecies=2
         call bmatrix(ispecies,0)
      endif
      if (ibmat_prod1.eq.1) then
         if (idebug.ge.1) write (6,*) 'computing prod1 bmat'
         ispecies=3
         call bmatrix(ispecies,0)
      endif
      if (ibmat_prod2.eq.1) then
         if (idebug.ge.1) write (6,*) 'computing prod2 bmat'
         ispecies=4
         call bmatrix(ispecies,0)
      endif
      if (ibmat_ts.eq.1) then
         if (idebug.ge.1) write (6,*) 'computing ts bmat'
         ispecies=0
         call bmatrix(ispecies,0)
      endif
      if (ibmat_wellr.eq.1) then
         if (idebug.ge.1) write (6,*) 'computing wellr bmat'
         ispecies=5
         if(igeom_wellr.eq.2.or.igeom_wellr.eq.1) ispecies=51
         call bmatrix(ispecies,0)
      endif
      if (ibmat_wellp.eq.1) then
         if (idebug.ge.1) write (6,*) 'computing wellp bmat'
         ispecies=6
         if(igeom_wellp.eq.2.or.igeom_wellp.eq.1) ispecies=61
         call bmatrix(ispecies,0)
      endif
      if (ibmat_test.eq.1) then
         if (idebug.ge.1) write (6,*) 'computing test bmat'
         call bmatrix(0,1)
      endif

c determine 1-dimensional torsional potentials
      if (i1dtau_reac1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting 1dtau_reac1'
         ispecies=1
         call onedtau(ispecies)
      endif
      if (i1dtau_reac2.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting 1dtau_reac2'
         ispecies=2
         call onedtau(ispecies)
      endif
      if (i1dtau_prod1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting 1dtau_prod1'
         ispecies=3
         call onedtau(ispecies)
      endif
      if (i1dtau_prod2.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting 1dtau_prod2'
         ispecies=4
         call onedtau(ispecies)
      endif
      if (i1dtau_wellr.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting 1dtau_wellr'
         ispecies=5
         if(igeom_wellr.eq.1.or.igeom_wellr.eq.2) ispecies=51
         call onedtau(ispecies)
      endif
      if (i1dtau_wellp.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting 1dtau_wellp'
         ispecies=6
         if(igeom_wellp.eq.1.or.igeom_wellp.eq.2) ispecies=61
         call onedtau(ispecies)
      endif
      if (i1dtau_ts.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting 1dtau_ts'
         ispecies=0
         call onedtau(ispecies)
      endif
c
c determine multi-dimensional torsional potentials
      if (imdtau_reac1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting mdtau_reac1'
         ispecies=1
         call mdtau(ispecies)
      endif
      if (imdtau_reac2.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting mdtau_reac2'
         ispecies=2
         call mdtau(ispecies)
      endif
      if (imdtau_prod1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting mdtau_prod1'
         ispecies=3
         call mdtau(ispecies)
      endif
      if (imdtau_prod2.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting mdtau_prod2'
         ispecies=4
         call mdtau(ispecies)
      endif
      if (imdtau_wellr.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting mdtau_wellr'
         ispecies=5
         if(igeom_wellr.eq.1.or.igeom_wellr.eq.2) ispecies=51
         call mdtau(ispecies)
      endif
      if (imdtau_wellp.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting mdtau_wellp'
         ispecies=6
         if(igeom_well.eq.1.or.igeom_wellp.eq.2) ispecies=61
         call mdtau(ispecies)
      endif
      if (imdtau_ts.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting mdtau_ts'
         ispecies=0
         call mdtau(ispecies)
      endif

c determine symmetry factors

c determine multi-dimensional torsional potentials
      if (isymm_reac1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting symm_reac1'
         ispecies=1
         call symmetry(ispecies)
      endif
      if (isymm_reac2.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting symm_reac2'
         ispecies=2
         call symmetry(ispecies)
      endif
      if (isymm_prod1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting symm_prod1'
         ispecies=3
         call symmetry(ispecies)
      endif
      if (isymm_prod2.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting symm_prod2'
         ispecies=4
         call symmetry(ispecies)
      endif
      if (isymm_wellr.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting symm_wellr'
         ispecies=5
         call symmetry(ispecies)
      endif
      if (isymm_wellp.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting symm_wellp'
         ispecies=6
         call symmetry(ispecies)
      endif
      if (isymm_ts.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting symm_ts'
         ispecies=0
         call symmetry(ispecies)
      endif

c determine higher level energies for ts, reactants, and products

c second part of the calculations, high level of nodes and memory      

      numproc=numprochl
      open (unit=99,status='unknown') 
      REWIND (99)
      write (99,*) gmemhigh
      REWIND (99)
      read (99,*) gmem
      close (99)

c      gmem=gmemhigh

      if (ihl_reac1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting hl_reac1'
         ispecies=1
         call hl(ispecies)
      endif
      if (ihl_reac2.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting hl_reac2'
         ispecies=2
         call hl(ispecies)
      endif
      if (ihl_reacs.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting hl_reacs'
         ispecies=11
         call hl(ispecies)
      endif
      if (ihl_prod1.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting hl_prod1'
         ispecies=3
         call hl(ispecies)
      endif
      if (ihl_prod2.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting hl_prod2'
         ispecies=4
         call hl(ispecies)
      endif
      if (ihl_prods.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting hl_prods'
         ispecies=12
         call hl(ispecies)
      endif
      if (ihl_wellr.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting hl_wellr'
         ispecies=5
         if(igeom_wellr.eq.1.or.igeom_wellr.eq.2) ispecies=51
         call hl(ispecies)
      endif
      if (ihl_wellp.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting hl_wellp'
         ispecies=6
         if(igeom_wellp.eq.1.or.igeom_wellp.eq.2) ispecies=61
         call hl(ispecies)
      endif
      if (ihl_ts.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting hl_ts'
         ispecies=0
         call hl(ispecies)
      endif

cc perform IRC of TS
      numproc=numprocll
      open (unit=99,status='unknown') 
      REWIND (99)
      write (99,*) gmemlow
      REWIND (99)
      read (99,*) gmem
      close (99)
c      gmem=gmemlow

      if (iirc.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting irc'
         call irc
      endif

c evaluate rate constant
      numproc=numprocll
      if (iktp.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting ktp'
         call ktp
      endif
c produced modified Arrhenius fit and output for CHEMKIN
      if (imodarr.eq.1) then
         if (idebug.ge.1) write (6,*) 'starting modarr'
         call modarr
      endif

      close(5)
      close(6)
      close(7)

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tau_opt(ispecies)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      real *8 mtrandom

      dimension tau(ntaumx),taumn(ntaumx),taumx(ntaumx),freq(nmdmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim) 
      dimension vref(noptmx),tauo(ntaumx,noptmx),
     $ xinto(3*natommx,noptmx),dist(noptmx)
      dimension tauopt(ntaumx)

      character*70 comline1,comline2
      character*60 atomlabel(natommx)
      character*30 intcoor(3*natommx)
      character*20 bislab(ntaumx)
      character*5 nameout
      character*100 command1
      character*100 commandcopy
      character*160 word_l0
      character*30 gmem
      character*30 distname

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

c initialize word, word2, word3, word4, word5
      call LineRead (0)
      commandcopy='cp -f ./data/level0_molpro.dat 
     $                    level0_molpro.dat'            

c input data
cc this call to ts is not really used here, we have tauo_ts
cc the file ts_opt is written by the TS opt routines
cs added option for species dependent electronic structure methods

      open (unit=21,file='./data/theory.dat',status='unknown')
      if (ispecies.eq.1) then
         open (unit=15,file='./data/reac1.dat',status='old')
         open (unit=16,file='./output/reac1.out',status='unknown')
         open (unit=17,file='./output/reac1_opt.out',status='unknown')
         nameout='reac1'
         rewind(21)
         WORD_L0='LEVEL0_1'
         do while (WORD.NE.'LEVEL0_1')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L0 = 'LEVEL0'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level0_reac1_molpro.dat 
     $                    level0_molpro.dat'            
         endif
       endif
      if (ispecies.eq.2) then
         open (unit=15,file='./data/reac2.dat',status='old')
         open (unit=16,file='./output/reac2.out',status='unknown')
         open (unit=17,file='./output/reac2_opt.out',status='unknown')
         nameout='reac2'
         rewind(21)
         WORD_L0='LEVEL0_2'
         do while (WORD.NE.'LEVEL0_2')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L0 = 'LEVEL0'
               go to 900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level0_reac2_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      endif
      if (ispecies.eq.3) then
         open (unit=15,file='./data/prod1.dat',status='old')
         open (unit=16,file='./output/prod1.out',status='unknown')
         open (unit=17,file='./output/prod1_opt.out',status='unknown')
         nameout='prod1'
         rewind(21)
         WORD_L0='LEVEL0_3'
         do while (WORD.NE.'LEVEL0_3')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L0 = 'LEVEL0'
               go to 900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level0_prod1_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      endif
      if (ispecies.eq.4) then
         open (unit=15,file='./data/prod2.dat',status='old')
         open (unit=16,file='./output/prod2.out',status='unknown')
         open (unit=17,file='./output/prod2_opt.out',status='unknown')
         nameout='prod2'
         rewind(21)
         WORD_L0='LEVEL0_4'
         do while (WORD.NE.'LEVEL0_4')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L0 = 'LEVEL0'
               go to 900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level0_prod2_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      endif
      if (ispecies.eq.5) then
         open (unit=15,file='./data/wellr.dat',status='old')
         open (unit=16,file='./output/wellr.out',status='unknown')
         open (unit=17,file='./output/wellr_opt.out',status='unknown')
         nameout='wellr'
         rewind(21)
         WORD_L0='LEVEL0_5'
         do while (WORD.NE.'LEVEL0_5')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L0 = 'LEVEL0'
               go to 900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level0_wellr_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      endif
      if (ispecies.eq.51) then
         open (unit=15,file='./data/wellr.dat',status='old')
         open (unit=16,file='./output/wellr.out',status='unknown')
         open (unit=17,file='./output/wellr_opt.out',status='unknown')
         open (unit=18,file='./output/ts_opt.out',status='unknown')
         nameout='wellr_l0'
c         inp_type=2
         rewind(21)
         WORD_L0='LEVEL0_51'
         do while (WORD.NE.'LEVEL0_51')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L0 = 'LEVEL0'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level0_wellr_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      endif
      if (ispecies.eq.6) then
         open (unit=15,file='./data/wellp.dat',status='old')
         open (unit=16,file='./output/wellp.out',status='unknown')
         open (unit=17,file='./output/wellp_opt.out',status='unknown')
         nameout='wellp'
         rewind(21)
         WORD_L0='LEVEL0_6'
         do while (WORD.NE.'LEVEL0_6')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L0 = 'LEVEL0'
               go to 900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level0_wellp_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      endif
      if (ispecies.eq.61) then
         open (unit=15,file='./data/wellp.dat',status='old')
         open (unit=16,file='./output/wellp.out',status='unknown')
         open (unit=17,file='./output/wellp_opt.out',status='unknown')
         open (unit=18,file='./output/ts_opt.out',status='unknown')
         nameout='wellp_l0'
         rewind(21)
         WORD_L0='LEVEL0_61'
         do while (WORD.NE.'LEVEL0_61')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L0 = 'LEVEL0'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level0_wellp_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      endif
      if (ispecies.eq.11) then
         open (unit=15,file='./data/reacs.dat',status='old')
         open (unit=16,file='./output/reacs.out',status='unknown')
         open (unit=17,file='./output/reacs_opt.out',status='unknown')
         nameout='reaso'
         rewind(21)
         WORD_L0='LEVEL0_11'
         do while (WORD.NE.'LEVEL0_11')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L0 = 'LEVEL0'
               go to 900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level0_reacs_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      endif
      if (ispecies.eq.12) then
         open (unit=15,file='./data/prods.dat',status='old')
         open (unit=16,file='./output/prods.out',status='unknown')
         open (unit=17,file='./output/prods_opt.out',status='unknown')
         nameout='proso'
         rewind(21)
         WORD_L0='LEVEL0_12'
         do while (WORD.NE.'LEVEL0_12')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L0 = 'LEVEL0'
               go to 900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level0_prods_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      endif
  900 continue
      rewind(21)

      ifreq=0
      ires=0
      ilev0code=1
      ilin=0
      ilin1=0
      ilin2=0

c print extra info during testing to see where things die
c nosmp is number of monte carlo sampling points
c dthresh is threshold for same geometry
c ethresh is threshold for same energy 

      do while (WORD.NE.'NOSMP')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (6,*) 'number of sampling points must be defined'
            stop
         endif
      enddo
      read (15,*) nosmp,dthresh,ethresh
      rewind(15)
c

      do while (WORD.NE.'NTAU')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (6,*) 'sampling coordinates must be defined'
            stop
         endif
      enddo
      read (15,*) ntau
      if (ntau.gt.ntaumx) then
         write (16,*) 'ntau too large',ntau,ntaumx
         stop
      endif
      if (ntau.ne.0) then
         read (15,*)
         do itau = 1 , ntau
            read (15,*) bislab(itau),taumn(itau),taumx(itau)
         enddo
      endif
      rewind(15)

c
c     natomt is to account for dummy atoms
      if(ispecies.ne.51.and.ispecies.ne.61)then
         do while (WORD.NE.'NATOM')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'natom must be defined'
               stop
            endif
         enddo
         read (15,*) natom,natomt,ilin
         if (natomt.gt.natommx) then
            write (16,*) 'natomt too large',natomt,natommx
            stop
         endif
         rewind(15)
      endif
c
      if(ispecies.eq.51.or.ispecies.eq.61)then
         open(unit=99,file='./data/reac1.dat',status='unknown')
         do while (WORD.NE.'NATOM')
            call LineRead (99)
            if (WORD.EQ.'END') then
               write (16,*) 'natom must be defined'
               write (16,*) 'in reac1 when using findgeom from level0'
               stop
            endif
         enddo
         read (99,*) natom1,natomt1,ilin1
         close(99)
         WORD=''
         if(iabs.eq.1.or.iadd.eq.1)then
            open(unit=99,file='./data/reac2.dat',status='unknown')
            do while (WORD.NE.'NATOM')
               call LineRead (99)
               if (WORD.EQ.'END') then
                  write (16,*) 'natom must be defined'
                  write (16,*) 'in reac2 with findgeom from level0'
                  stop
               endif
            enddo
            read (99,*) natom2,natomt2,ilin2
            close(99)
         endif
         natom=natom1+natom2
         natomt=natomt1+natomt2
         if(iabs.eq.1)natomt=natomt+1
cc if linear molecule reacting with atom in abs assume TS is linear
         if(iabs.eq.1.and.natom2.eq.1.
     $        and.ilin1.eq.1) then
            ilin=1
         endif

cc if linear molecule reacting with linear radical in abs assume TS is linear
         if(iabs.eq.1.and.ilin2.eq.1.and
     $        .ilin1.eq.1)then
            ilin=1
         endif
      endif

      do while (WORD.NE.'CHARGE')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (16,*) 'charge and spin must be defined'
            stop
         endif
      enddo
      read (15,*) icharge,ispin
      if(ispecies.ne.51.and.ispecies.ne.61)then
         do iatom = 1 , natomt
            read (15,'(A60)') atomlabel(iatom)
         enddo
      endif
      rewind(15)

      do while (WORD.NE.WORD_L0)
         call LineRead (21)
         if (WORD.EQ.'END') then
            write (16,*) 'level0 of theory must be defined'
            write (16,*) 'in theory.dat'
            stop
         endif
      enddo
      if(word2.eq.'G09') then
         ilev0code=1
      else if (word2.eq.'MOLPRO') then
         ilev0code=2
         call commrun(commandcopy)
      else
         write (16,*) 'level0 of theory must be either'
         write (16,*) 'g09 or molpro in theory.dat'
         stop
      endif
      if(word3.eq.'RESTART') then
         open (unit=99,status='unknown')
         write (99,*) word4
         rewind (99)
         read (99,*) ires
         close (unit=99,status='keep')
      endif
      if(ilev0code.eq.1)then
         read (21,'(A70)') comline1
         read (21,'(A70)') comline2
      endif
      close(21)


      if (idebug.ge.2) write (16,*) 'level of theory ',comline1,comline2

      if (natom.ne.1.and.ispecies.ne.51.and.ispecies.ne.61) then
         do while (WORD.NE.'INTCOOR')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'internal coordinates must be defined'
               stop
            endif
         enddo
         ncoord = 3*natom-6-ntau
         if (natom.eq.1) ncoord = 0
         if (natom.eq.2) ncoord = 1
         do icoord = 1 , ncoord
            call LineRead (15)
            intcoor(icoord) = word
            if(intcoor(icoord).eq.'RTS'.and.ispecies.eq.1.
     $ and.ibeta.ne.1)then
               write(*,*)'no coord name can be RTS in reac1.dat'
               write(*,*)'change the name and restart'
               stop
            endif
            OPEN (unit=99,status='unknown')
            REWIND (99)
            WRITE (99,1000) WORD2
 1000       FORMAT (A30)
            REWIND (99)
            READ (99,*) xinti(icoord)
            close (unit=99,status='keep')
         enddo
         rewind(15)
      endif
c gaussian com file data

c     if (idebug.ge.2) write (6,*) 'word_l0 test',word_l0
c     if (idebug.ge.2) write (6,*) 'past z-matrix'
      if (idebug.ge.2) write (16,*) 'past z-matrix'

c initializations
c     write (16,*) ' starting intializations'
      nint = 3*natom-ntau-6
      if (natom.eq.1) nint=0
      if (natom.eq.2) nint=1

c start sampling

      if (idebug.ge.2) write (16,*) 'number of sampling points',nosmp
      if (idebug.ge.2) write (16,*) 'number of sampled coordinates',ntau
      if (idebug.ge.2) write (16,*) 'number of other coordinates',nint
      noptg = 0
      do ismp = 1 , nosmp
         if (idebug.ge.2) write (16,*) '*******************************'
         if (idebug.ge.2) write (16,*) 'Random number point ',ismp,
     $                                   ' out  of ',nosmp

c generate random tau

         do itau = 1 , ntau
            tau(itau) = taumn(itau)+(taumx(itau)-taumn(itau))*
     $       mtrandom()
         enddo
         write(16,*)'Sampled Tau coordinates'
         do itau = 1 , ntau
            write(16,*)bislab(itau),tau(itau)
         enddo


cc determine distance to check for the findgeom option for wellp and wellr

         if(ispecies.eq.51.or.ispecies.eq.61) then
            open (unit=99,file='data/ts.dat',status='old')
            do while (WORD.NE.'ISITE')
               call LineRead (99)
               if (WORD.EQ.'END') then
                  write (16,*) 'when using findgeom option'
                  write (16,*) 'from level 0'
                  write (16,*) 'reaction site must be defined'
                  write (16,*) 'in file ts.dat'
                  stop
               endif
            enddo
            read (99,*) isite,jsite,ksite
c check if reverse index
            ireverse=0
            if(iabs.eq.1)then
               rewind(99)
               do while (WORD.NE.'RMIN')
                  call LineRead (99)
                  if (WORD.EQ.'END') then
                     write (16,*) 'when using findgeom option'
                     write (16,*) 'from level 0'
                     write (16,*) 'rmin must be defined'
                     write (16,*) 'in file ts.dat'
                     stop
                  endif
               enddo
               if(WORD2.eq.'REVERSE')ireverse=1
            endif
            write(*,*)'ireverse is',ireverse

c determine number of reacting atom for beta-scission reactions

            if(ibeta.eq.1)then
               rewind(99)
               do while (WORD.NE.'RMAX1')
                  call LineRead (99)
                  if (WORD.EQ.'END') then
                     write (16,*) 'when using findgeom option'
                     write (16,*) 'from level 0'
                     write (16,*) 'grid coords must be defined'
                     write (16,*) 'in file ts.dat'
                     stop
                  endif
               enddo
               read (99,*) rmax1,rstp1,rmax2,rstp2,ireact
            endif
            close(99)
            
            write(*,*)'isite is',isite

c for beta-decomp determine label of breaking bond      
            distname='RTS'
            write(16,*)'ireact is ', ireact
            if(ibeta.eq.1)then
               open(unit=99,status='unknown')
               do iatom = 1 , natomt
                  rewind (99)
                  write(99,*)atomlabel(iatom)
                  rewind (99)
                  call LineRead(99)
                  if(iatom.eq.ireact)then
                     distname=word3
                  endif
               enddo
               close(99)
               write(16,*)'distname is ', distname
            endif
         endif
c         stop
c        write (6,*) 'starting g09pot'
         if(ispecies.ne.51.and.ispecies.ne.61)then
            do iint = 1 , nint
               xint(iint) = xinti(iint)
            enddo
         endif

         if(ispecies.eq.51.or.ispecies.eq.61)then
            read (18,*)
            do iatom = 1 , natomt
               read (18,'(A60)') atomlabel(iatom)
            enddo
c for abstraction in reverse direction label of breaking bond
            if(iabs.eq.1.and.ireverse.eq.1)then
               open(unit=99,status='unknown')
               do iatom = 1 , natomt
                  rewind (99)
                  write(99,*)atomlabel(iatom)
                  rewind (99)
                  call LineRead(99)
                  if(iatom.eq.isite)then
                     distname=word3
                  endif
               enddo
            endif
            write(16,*)'distname is ', distname
c            stop
            do iint = 1 , nint
               read (18,*) intcoor(iint),xint(iint)
               write (*,*) intcoor(iint),xint(iint)
               if(ispecies.eq.51)then
                  if(intcoor(iint).eq.distname) then
                     if(distcheck.gt.1.4.and.iabs.eq.1)then
                        fgdisp51=0.5
                     else if(distcheck.gt.1.4.and.iadd.eq.1)then
                        fgdisp51=0.5
                     else
                        fgdisp51=0.2
                     endif
                     if(ibeta.eq.1)fgdisp51=-fgdisp51
                     if(iabs.eq.1.and.ireverse.eq.1)fgdisp51=-fgdisp51
                     xint(iint)=xint(iint)+fgdisp51
                  endif
               endif
               if(ispecies.eq.61)then
                  if(intcoor(iint).eq.distname) then
                     if(xint(iint).gt.1.4) then
                        fgdisp61=0.5
                     else
                        fgdisp61=0.2
                     endif
                     if(ibeta.eq.1)fgdisp61=-fgdisp61
                     if(iabs.eq.1.and.ireverse.eq.1)fgdisp61=-fgdisp61
                     xint(iint)=xint(iint)-fgdisp61
                  endif
               endif
            enddo
            close(18)
         endif
         ircons=0
         ixyz=0
         ired=0
         ilev=0
         if(ispecies.eq.51.or.ispecies.eq.61)then
            ilin=0
            if(ilin1.eq.0) then
               ired=1
            else if(ilin1.eq.1)then
               ired=2
            endif
         endif
         if(iabs.eq.1) ireact=natom1+1
         if(iadd.eq.1) ireact=natom1+1

         write (16,*) 
         if(ilev0code.eq.1) then
            write (16,*) 'starting g09pot'
            call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $           coord,vtot_0,vtot,freq,ifreq,ilin,ismp,
     $           comline1,comline2,icharge,ispin,ircons,
     $           atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires
     $           ,ixyz,ired)
            write (16,*) 'finished g09pot'

         else if (ilev0code.eq.2) then
            write (16,*) 'starting molprofopt'         
            call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $           coord,vtot,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $           atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,
     $           ilev,ispecies)
            write (16,*) 'finished molprofopt'         
         endif

         write (16,*) 
         write (16,*) 'optimized tau coordinates,energy,sampling point'
         write (16,*) (tauopt(itau),itau=1,ntau),vtot,ismp
         write (16,*) 'optimized coordinates'
         write (16,1111) (xint(it),it=1,nint)
 1111    format (100(1X,f7.2))

         if (vtot.gt.1.d10) go to 5000
cc if natom=1 the above condition is satisfied
c determine if new geometry
         if (idebug.ge.2) write (16,*) 'geom diff for geom iopt 
     $ and threshold'
         dmin = 1.d20
         do iopt = 1 , noptg
            dist(iopt) = 0.0d0
            do itau = 1 , ntau
               dref = abs(tauopt(itau)-tauo(itau,iopt))
               dplus = abs(tauopt(itau)+360.d0-tauo(itau,iopt))
               dminus = abs(tauopt(itau)-360.d0-tauo(itau,iopt))
               dd = min(dref,dplus,dminus)**2
               dist(iopt) = dist(iopt) + dd
            enddo
            dist(iopt)=sqrt(dist(iopt))
            if (idebug.ge.2) write (16,*) dist(iopt),iopt,dthresh
            if (dist(iopt).lt.dthresh) go to 5000
         enddo
c enter this point only if it is a new geometry (thus also at step 1)
c save tau, internal coordinates and vibrational frequencies for each new geometry
cc and save in output xyz file for each new geometry
         noptg = noptg+1
         vref(noptg) = vtot
         do itau = 1 , ntau
            tauo(itau,noptg) = tauopt(itau)
         enddo 
         if (idebug.ge.2) write (6,*) noptg
         do iint = 1 , nint
            xinto(iint,noptg) = xint(iint)
         enddo
         write (16,*) 'new geom'
         write (16,*) 'tau ' 
         write (16,*) (tauopt(itau),itau=1,ntau)
         write (16,*) 'internal coordinates ' 
         do iint = 1 , nint
            write (16,*) intcoor(iint),xint(iint)
         enddo
         write (16,*) 'coordinates '
         do iatom = 1 ,natomt
            write (16,*) (coord(iatom,idim),idim=1,ndim)
         enddo
         write (16,*) 'rotational constants'
         write (16,*) (abcrot(idim),idim=1,ndim)
         write (16,*) 'energy and its difference with minimum energy'
         write (16,*) vtot, (vref(iopt)-vref(1))
         write (16,*)

c         command1='newzmat       
c     +           -ichk -oxyz tmp.chk geom.xyz           '
c         call commrun(command1)

         open (unit=98,file='headgeom.tmp',status='unknown')
         write(98,*)natom
         write(98,*)vtot
         close(98)
         command1='cat headgeom.tmp geom.xyz > geom1.xyz '
         call commrun(command1)
         
         open (unit=99,status='unknown')
         rewind (99)
         write (99,1112) nameout,noptg
         rewind (99)
         read (99,1020) command1
         call commrun(command1)         
         close(99)

         open (unit=98,file='optgeom.tmp',status='unknown')

         write (98,*) 'opt geom',noptg
         do iint = 1 , nint
            write (98,*) xinto(iint,noptg)
         enddo
         do itau=1,ntau
            write (98,*) tauo(itau,noptg)
         enddo
         write (98,*) vref(noptg) 
         close(98)

         open (unit=99,status='unknown')
         rewind (99)
         write (99,1113) nameout,noptg
         rewind (99)
         read (99,1020) command1
         call commrun(command1)         
         close(99)

 1112    format ('cp -f geom1.xyz geoms/'A5,'_',I0.2,'.xyz',60X)
 1113    format ('cp -f optgeom.tmp output/'A5,'_opt_',I0.2,'.out',60X)
 1020    format (A70)

         continue


c check to see if this is also a new energy
         inew = 1
         if (noptg.gt.1) then
cc         if (noptg.gt.1) then
cc            do iopt = 1 , noptg-1
            do iopt = 1 , noptg
               dv = abs(vref(iopt)-vref(noptg))
               if (dv.lt.ethresh) then
                  inew = 0
               endif
            enddo 
         endif
         if (inew.eq.1) then 
c            do itau = 1 , ntau
c               tauo2(itau,noptg) = tauo(itau,noptg)
c            enddo 
c            do iint = 1 , nint
c               xinto2(iint,noptg) = xinto(iint,noptg)
c            enddo
            write (16,*) 'new energy also'
            write (16,*)
         endif
         write (16,*) '*******************************'
 5000    continue
      enddo

c for the optimum geometry print out the optimized coordinates
c
      ioptmin = 1
      voptmin= vref(ioptmin)
      do iopt = 1 , noptg
         if (vref(iopt).lt.voptmin) then
            ioptmin=iopt
            voptmin= vref(ioptmin)
         endif
      enddo
      if (idebug.ge.2) write (16,*) 'noptg voptmin ioptmin',
     $                noptg,voptmin,ioptmin
      write (17,*) 'opt geom',ioptmin
      if(ispecies.eq.51.or.ispecies.eq.61)then
c     $   or.ispecies.eq.5.or.ispecies.eq.6)then
         do iatom = 1 , natomt
            write (17,'(A60)') atomlabel(iatom)
         enddo
         do iint = 1 , nint
            write (17,*)intcoor(iint), xinto(iint,ioptmin)
         enddo
         do itau=1,ntau
            write (17,*)bislab(itau), tauo(itau,ioptmin)
         enddo
      else
         do iint = 1 , nint
            write (17,*) xinto(iint,ioptmin)
         enddo
         do itau=1,ntau
            write (17,*) tauo(itau,ioptmin)
         enddo
      endif
c
cc    modified as not consistent with the read statement in the ts_zmat routine
c      write (17,*) (tauo(itau,ioptmin),itau=1,ntau),
c     $ vref(iopt)
      write (17,*) vref(ioptmin) 
      close (unit=15,status='keep')
      close (unit=16,status='keep')
      close (unit=17,status='keep')

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine grid_opt

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension tau(ntaumx),taumn(ntaumx),taumx(ntaumx),freq(nmdmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim)
      dimension xints(3*natommx)
      dimension xintsave(3*natommx)
      dimension tauopt(ntaumx)

      character*70 comline1,comline2
      character*60 atomlabel(natommx)
      character*30 intcoor(3*natommx)
      character*30 intcoors(3*natommx)
      character*30 intcoorsave(3*natommx)
      character*20 bislab(ntaumx)
      character*30 gmem
      character*30 revcoo
      character*100 commandcopy

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

cc parameter initialization
      ires=0
cc we start with a non linear TS
      ilin=0
      ilin1=0
      ilin2=0
c input data

      open (unit=25,file='./data/ts.dat',status='old')
      open (unit=26,file='./output/grid_opt.out',status='unknown')
      open (unit=21,file='./data/theory.dat',status='unknown')


c     write (26,*) 'idebug test',idebug

cc get data from react1 file

      open (unit=15,file='./data/reac1.dat',status='old')

      do while (WORD.NE.'NTAU')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (6,*) 'sampling coordinates of reac1 must be defined'
            stop
         endif
      enddo
      read (15,*) ntau1
      rewind 15

c     natomt is to account for dummy atoms
      do while (WORD.NE.'NATOM')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (16,*) 'natom must be defined'
            stop
         endif
      enddo
      read (15,*) natom1,natomt1,ilin1
      close (15)

cc get data from react2 file

      open (unit=15,file='./data/reac2.dat',status='old')

      do while (WORD.NE.'NTAU')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (6,*) 'sampling coordinates of reac2 must be defined'
            stop
         endif
      enddo
      read (15,*) ntau2
      rewind 15

c     natomt is to account for dummy atoms
      do while (WORD.NE.'NATOM')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (26,*) 'natom must be defined'
            stop
         endif
      enddo
      read (15,*) natom2,natomt2,ilin2
      close (15)

c which atom to add or abstract and atoms to connect to it
c     do while (WORD.NE.'REACTIONTYPE')
c        call LineRead (25)
c        if (WORD.EQ.'END') then
c           write(26,*) 'reaction type must be defined'
c           stop
c        endif
c     enddo
c     iabs=0
c     iadd=0
c     if (WORD2.EQ.'ABSTRACTION') iabs = 1
c     if (WORD2.EQ.'ADDITION') iadd = 1
c     rewind 25

      do while (WORD.NE.'ISITE')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (26,*) 'reaction site must be defined'
            stop
         endif
      enddo
      read (25,*) isite,jsite,ksite
      rewind 25

      ireverse=0
      do while (WORD.NE.'RMIN')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (26,*) 'grid coords must be defined'
            stop
         endif
      enddo
      if(word2.eq.'REVERSE')then
         ireverse=1
      endif
      read (25,*) rmin,rmax,nr
      read (25,*)
      read (25,*) aabs1,babs1,aabs2,babs2,babs3
      rewind 25

      do while (WORD.NE.'CHARGE')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (26,*) 'charge and spin must be defined'
            stop
         endif
      enddo
      read (25,*) icharge,ispin
      rewind(25)

      do while (WORD.NE.'LEVEL0')
         call LineRead (21)
         if (WORD.EQ.'END') then
            write (26,*) 'gaussian level0_ts of theory must be defined'
            write (26,*) 'in file theory.dat'
            stop
         endif
      enddo
      if(word2.eq.'G09') then
         ilev0code=1
      else if (word2.eq.'MOLPRO') then
         ilev0code=2
         commandcopy='cp -f ./data/level0_molpro.dat 
     $                    level0_molpro.dat'            
         call commrun(commandcopy)
      else
         write (26,*) 'level0 of theory must be either'
         write (26,*) 'g09 or molpro in theory.dat'
         stop
      endif
      if(word3.eq.'RESTART') then
         open (unit=99,status='unknown')
         write (99,*) word4
         rewind (99)
         read (99,*) ires
         close (unit=99,status='keep')
      endif
      if(ilev0code.eq.1)then
         read (21,'(A70)') comline1
         read (21,'(A70)') comline2
      endif
      close(21)

      if (idebug.ge.2) write (26,*) ' comline test',comline1,comline2


      natom = natom1+natom2
c     natomt = natomt1+natomt2
      ntau = ntau1+ntau2
      rstp = 0.0d0
      if (nr.gt.1) then
         rstp = (rmax-rmin)/(float(nr-1))
      endif

cc if linear molecule reacting with atom in abs assume TS is linear
      if(iabs.eq.1.and.natom2.eq.1.
     $     and.ilin1.eq.1) then
         ilin=1
      endif

cc if linear molecule reacting with linear radical in abs assume TS is linear
      if(iabs.eq.1.and.ilin2.eq.1.and
     $     .ilin1.eq.1)then
         ilin=1
      endif

c build z-mat input

      rts = 0.d0
      ismp=0
      
      if (idebug.ge.2) write (26,*) 'entering ts_zmat'

      call ts_zmat(atomlabel,xinti,
     $ intcoor,natom1,natom2,natomt1,natomt2,ntau1,ntau2,bislab)

      if (idebug.ge.2) write (26,*) 'finished ts_zmat'

      vrimax = -1.0d20

      if (iadd.eq.1) natomt = natomt1+natomt2
      if (iabs.eq.1) natomt = natomt1+natomt2+1
      if(ireverse.eq.1)then
         open (unit=99,status='unknown')
         write (99,*) atomlabel(isite)
         rewind (99)
         call lineread (99)
         close (unit=99,status='keep')
         revcoo=word3
         write(*,*)'revcoo is', revcoo
c         stop
      endif
      ntau = 0
      nint = 3*natom-ntau-6
c      if(natom1.eq.2)nint=nint+1
      do iint = 1 , nint
         xint(iint) = xinti(iint)
      enddo
c if reverse mode find intcoor and move to the end
      revcoo_st=0.
      if(ireverse.eq.1)then
         do iint = 1 , nint
            xintsave(iint)=xint(iint)
            intcoorsave(iint)=intcoor(iint)
            if(intcoor(iint).eq.revcoo)then
               revcoo_st=xint(iint)
            endif
         enddo
         do iint = 1 , nint
            if(intcoorsave(iint).eq.revcoo)then
               intcoor(iint)='RTS'
            endif
         enddo
         intcoor(nint)=revcoo
      endif
      do ir = 1 , nr
         rts = rmin+float(ir-1)*rstp
c initializations
c     write (16,*) ' starting intializations'
c        natomt = natomt1+natomt2+1
         xinti(nint) = rts
         ijump=0
         if(ireverse.eq.1.and.xinti(nint).lt.revcoo_st)then
            ijump=1
         endif
c call ab initio potential

c        write (6,*) 'starting g09pot'
         if(ifrozrts.ne.1)then
            do iint = 1 , nint
               xint(iint) = xinti(iint)
            enddo
         else
            xint(nint) = xinti(nint)
         endif

c        if (natom2.eq.1) xinti(natomtp+1) = rts
c        if (natom2.eq.2) xinti(natomtp+1) = rts
c        if (natom2.ge.3) xinti(natomtp+1) = rts
         ircons=1
         ifreq=0
         ixyz=0
         ired=0
         ilev=0
         if(ijump.eq.1) then
            vtotr=1.0d20
            goto 4999
         endif
         if(ilev0code.eq.1) then
            if (idebug.ge.2) write (26,*) 'entering g09fopt'
            call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $           coord,vtot_0,vtotr,freq,ifreq,ilin,ismp,
     $           comline1,comline2,icharge,ispin,ircons,
     $           atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,
     $           ixyz,ired)
            write (26,*) 'finished g09pot'
         else if (ilev0code.eq.2) then
            write (26,*) 'starting molprofopt'         
            call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $           coord,vtotr,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $           atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,
     $           ilev,ispecies)
            write (26,*) 'finished molprofopt'         
         endif
 4999    continue
         write (26,*) ir,rts,vtotr
c find maximum on grid of rts values
         if (vtotr.gt.1.d10) go to 5000
         if (vtotr.gt.vrimax) then 
            irs = ir
            vrimax = vtotr
            do iint = 1 , nint
               intcoors(iint) = intcoor(iint)
               xints(iint) = xint(iint)
            enddo
         endif
 5000    continue
      enddo

cc if search form reverse direction put RTS back at end of coord vector

      if(ireverse.eq.1)then
         do iint = 1 , nint
            xintsave(iint)=xints(iint)
            intcoorsave(iint)=intcoors(iint)
         enddo
         iprog=0
         do iint = 1 , nint
            if(intcoorsave(iint).ne.'RTS')then
               iprog=iprog+1
               intcoors(iprog)=intcoorsave(iint)
               xints(iprog)=xintsave(iint)
            else
               intcoors(nint)='RTS'
               xints(nint)=xintsave(iint)
            endif
         enddo
      endif

c save key data

      write (26,*) 'grid optimized z-matrix'
      do iatom = 1 , natomt
         write (26,'(A60)') atomlabel(iatom)
      enddo
      do iint = 1 , nint
         write (26,*) intcoors(iint),xints(iint)
         intcoor(iint) = intcoors(iint)
         xinti(iint) = xints(iint)
      enddo
      natomtp = 3*(natom1+natom2)-12
      aabs1 = xints(natomtp+1)
      babs1 = xints(natomtp+2)
      aabs2 = xints(natomtp+3)
      babs2 = xints(natomtp+4)
      babs3 = xints(natomtp+5)
      rts = xints(natomtp+6)

      close (unit=25,status='keep')
      close (unit=26,status='keep')
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine grid_opt_iso

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension tau(ntaumx),taumn(ntaumx),taumx(ntaumx),freq(nmdmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim),drea(ntaumx)
      dimension xints(3*natommx)
      dimension tauopt(ntaumx)

      character*70 comline1,comline2
      character*60 atomlabel(natommx)
      character*60 atomlabel1(natommx)
      character*30 intcoor(3*natommx)
      character*30 intcoors(3*natommx)
      character*20 bislab(ntaumx)
      character*30 gmem
      character*70 command1
      character*40 filename
      character*4  cjunk
      character*100 commandcopy
      character*30 dih_rot

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

cc parameter initialization
      ires=1

c input data

      open (unit=25,file='./data/ts.dat',status='old')
      open (unit=26,file='./output/grid_opt.out',status='unknown')
      open (unit=21,file='./data/theory.dat',status='unknown')

cc get data from react1 file

      open (unit=15,file='./data/reac1.dat',status='old')

      do while (WORD.NE.'NTAU')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (6,*) 'sampling coordinates of reac1 must be defined'
            stop
         endif
      enddo
      read (15,*) ntau1
      rewind 15

c     natomt is to account for dummy atoms
      do while (WORD.NE.'NATOM')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (26,*) 'natom must be defined'
            stop
         endif
      enddo
      read (15,*) natom1,natomt1
      close (15)

      do while (WORD.NE.'ISITE')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (26,*) 'reaction site must be defined'
            stop
         endif
      enddo
      read (25,*) isite,jsite,ksite
      rewind 25

      do while (WORD.NE.'RMIN1')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (26,*) 'grid coords must be defined'
            stop
         endif
      enddo
      read (25,*) rmin1,rstp1,rmin2,rstp2,ireact
      rewind 25

      do while (WORD.NE.'CHARGE')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (26,*) 'charge and spin must be defined'
            stop
         endif
      enddo
      read (25,*) icharge,ispin
      rewind(25)

      do while (WORD.NE.'LEVEL0')
         call LineRead (21)
         if (WORD.EQ.'END') then
            write (26,*) 'gaussian level0_ts of theory must be defined'
            write (26,*) 'in file theory.dat'
            stop
         endif
      enddo
      if(word2.eq.'G09') then
         ilev0code=1
      else if (word2.eq.'MOLPRO') then
         ilev0code=2
         commandcopy='cp -f ./data/level0_molpro.dat 
     $                    level0_molpro.dat'            
         call commrun(commandcopy)
      else
         write (26,*) 'level0 of theory must be either'
         write (26,*) 'g09 or molpro in theory.dat'
         stop
      endif
      if(word3.eq.'RESTART') then
         open (unit=99,status='unknown')
         write (99,*) word4
         rewind (99)
         read (99,*) ires
         close (unit=99,status='keep')
      endif
      if(ilev0code.eq.1)then
         read (21,'(A70)') comline1
         read (21,'(A70)') comline2
      endif
      close(21)

cc now determine the structure to use in the calculations

cc first determine how many structures were determined with the rotational scan

      command1='ls output/reac1_opt_* > temp.log '
      call commrun(command1)
      command1='grep -c ^ temp.log > num.log'
      call commrun(command1)
      open (unit=99,file='num.log',status='unknown')
      read (99,*) numstruct
      close (unit=99,status='keep')

cc now open the xyz files and search the one with the minimum distance
cc between the reacting atoms
      write(*,*)'NUMSTR IS', numstruct
      do j=1,numstruct
         open (unit=99,status='unknown')
         write(99,1201)j
         rewind (99)
         read (99,'(A40)') filename
         close (99)
         write(*,*) 'filename is ',filename
         open(unit=23,file=filename,status='unknown')
         read (23,*)
         read (23,*)
         do k=1,natom1
            read (23,*)cjunk,coox,cooy,cooz
            if(k.eq.ireact)then
               coox1=coox
               cooy1=cooy
               cooz1=cooz
            endif
            if(k.eq.isite)then
               coox2=coox
               cooy2=cooy
               cooz2=cooz
            endif
         enddo
         close(23)
         drea(j)=sqrt((coox1-coox2)**2+(cooy1-cooy2)**2
     $           +(cooz1-cooz2)**2)
         write(26,*)'irea-isite dist for str ',j,' is ',drea(j)
      enddo

cc now determine the structure with the minimum distance
      iopt=1
      do j=1,numstruct
         if(drea(j).lt.drea(iopt))then
            iopt=j
         endif
      enddo
      write(26,*) 'the structure with the min dist is',iopt

c      stop

      natom = natom1
      ntau = ntau1

c build z-mat input

      rts = 0.d0
      ismp=0
c      iopt=0

      if (idebug.ge.2) write (26,*) 'entering ts_zmat'

      call ts_zmat_iso(atomlabel,atomlabel1,xinti,
     $ intcoor,natom1,natomt1,ntau1,iopt,dstart)

      if (idebug.ge.2) write (26,*) 'finished ts_zmat'

      do iatom = 1 , natomt1
         atomlabel(iatom) = atomlabel1(iatom)
      enddo

      vrimax = -1.0d20
c      nr2=(dstart-rmin2)/rstp2+0.5
c      nr1=(rmin2-rmin1)/rstp1+0.5

      nr2=abs(dstart-rmin1)/rstp1+0.5
      nr1=abs(rmin1-rmin2)/rstp2+0.5

      nr=nr1+nr2

      write(26,*)'nr1 is',nr1
      write(26,*)'nr2 is',nr2
      write(26,*)'nr is',nr

      nint = 3*natom-6
      do iint = 1 , nint
         xint(iint) = xinti(iint)
      enddo

cc determine dihedral of reacting atom
      
      if(ireact.gt.3)then
         open (unit=99,status='unknown')
         write(99,*)atomlabel(ireact)
         rewind(99)
         read(99,*)cjunk,cjunk,cjunk,cjunk,cjunk,cjunk,dih_rot
         close(99)
c      write(*,*)'dname is ',dih_rot

cc update dihed variable to avoid planar guess

         do iint = 1 , nint
            if(dih_rot.eq.intcoor(iint))then
               check1=abs(xint(iint)-180.)
               check2=abs(xint(iint)-0.)
c            write(*,*)'check1 is ',check1
c            write(*,*)'check2 is ',check2
               if(check1.lt.1.)then
                  xint(iint)=xint(iint)+10.
                  check1=10.
               endif
               if(check2.lt.1.)then
                  xint(iint)=xint(iint)+10.
                  check2=10.
               endif
            endif
         enddo
      else if (ireact.eq.3) then
         open (unit=99,status='unknown')
         write(99,*)atomlabel(ireact+1)
         rewind(99)
         read(99,*)cjunk,cjunk,cjunk,cjunk,cjunk,cjunk,dih_rot
         rewind(99)
         write(99,*)dih_rot
         rewind(99)
         call lineread(99)
         dih_rot=word
         
c         read(99)dih_rot
c         read(99,*)cjunk,cjunk,cjunk,cjunk,cjunk,cjunk,dih_rot
         close(99)
c      write(*,*)'dname is ',dih_rot
c      stop
cc update dihed variable to avoid planar guess

         do iint = 1 , nint
            if(dih_rot.eq.intcoor(iint))then
               check1=abs(xint(iint)-180.)
               check2=abs(xint(iint)-0.)
            write(*,*)'check1 is ',check1
            write(*,*)'check2 is ',check2
               if(check1.lt.1.)then
                  xint(iint)=xint(iint)+10.
                  check1=10.
               endif
               if(check2.lt.1.)then
                  xint(iint)=xint(iint)+10.
                  check2=10.
               endif
            endif
         enddo
      endif
c     stop

cc start cycling

      do ir = 1 , nr
         if(ir.le.nr2) then
            if(rmin2.lt.rmin1)then
               rts = dstart-float(ir)*rstp1
            else
               rts = dstart+float(ir)*rstp1
            endif
         else
            if(rmin2.lt.rmin1)then
               rts = rmin1-float(ir-nr2)*rstp2
            else
               rts = rmin1+float(ir-nr2)*rstp2
            endif
         endif

c         write(26,*)'rts is',rts
c      enddo
c      stop
c      do ir = 1 , nr

c initializations
c     write (16,*) ' starting intializations'
         ntau = 0
         natomt = natomt1
c        natomt = natomt1+natomt2+1
         xint(nint) = rts

cc update dihed variable to avoid planar guess

      do iint = 1 , nint
         if(dih_rot.eq.intcoor(iint))then
            check1=abs(xint(iint)-180.)
            check2=abs(xint(iint)-0.)
c            write(*,*)'check1 is ',check1
c            write(*,*)'check2 is ',check2
            if(check1.lt.1.)then
               xint(iint)=xint(iint)+10.
               check1=10.
            endif
            if(check2.lt.1.)then
               xint(iint)=xint(iint)+10.
               check2=10.
            endif
         endif
      enddo

c call ab initio potential

c        write (6,*) 'starting g09pot'
c        if (natom2.eq.1) xinti(natomtp+1) = rts
c        if (natom2.eq.2) xinti(natomtp+1) = rts
c        if (natom2.ge.3) xinti(natomtp+1) = rts
cc we assume the TS is not linear         
         ilin=0
         ircons=1
         ifreq=0
         ixyz=0
         ired=0
         ilev=0

         if(ilev0code.eq.1) then
            if (idebug.ge.2) write (26,*) 'entering g09fopt'
            call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $           coord,vtot_0,vtotr,freq,ifreq,ilin,ismp,
     $           comline1,comline2,icharge,ispin,ircons,
     $           atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires
     $           ,ixyz,ired)
            write (26,*) 'finished g09pot'
         else if (ilev0code.eq.2) then
            write (26,*) 'starting molprofopt'         
            call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $           coord,vtotr,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $           atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,
     $           ilev,ispecies)
            write (26,*) 'finished molprofopt'         
         endif

         write (26,*) ir,rts,vtotr
c         stop
c find maximum on grid of rts values
         if (vtotr.gt.1.d10) go to 5000
         if (vtotr.gt.vrimax) then 
            irs = ir
            vrimax = vtotr
            do iint = 1 , nint
               intcoors(iint) = intcoor(iint)
               xints(iint) = xint(iint)
            enddo
         endif
 5000    continue
      enddo
c save key data
      write (26,*) 'grid optimized z-matrix'
      do iatom = 1 , natomt
         write (26,*) atomlabel(iatom)
      enddo
      do iint = 1 , nint
         write (26,*) intcoors(iint),xints(iint)
         intcoor(iint) = intcoors(iint)
         xinti(iint) = xints(iint)
      enddo
c      natomtp = 3*(natom1+natom2)-12
c      aabs1 = xints(natomtp+1)
c      babs1 = xints(natomtp+2)
c      aabs2 = xints(natomtp+3)
c      babs2 = xints(natomtp+4)
c      babs3 = xints(natomtp+5)
c      rts = xints(natomtp+6)

      close (unit=25,status='keep')
      close (unit=26,status='keep')

 1201 format ("./geoms/reac1_"I0.2".xyz")

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine grid_opt_beta

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension tau(ntaumx),taumn(ntaumx),taumx(ntaumx),freq(nmdmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim),drea(ntaumx)
      dimension xints(3*natommx)
      dimension tauopt(ntaumx)

      character*70 comline1,comline2
      character*60 atomlabel(natommx)
      character*60 atomlabel1(natommx)
      character*30 intcoor(3*natommx)
      character*30 intcoors(3*natommx)
      character*20 bislab(ntaumx)
      character*30 gmem
      character*70 command1
      character*40 filename
      character*4  cjunk
      character*100 commandcopy

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

cc parameter initialization
      ires=1

c input data
      write(*,*) 'passing from here'

      open (unit=25,file='./data/ts.dat',status='old')
      open (unit=26,file='./output/grid_opt.out',status='unknown')
      open (unit=21,file='./data/theory.dat',status='unknown')

cc get data from react1 file

      open (unit=15,file='./data/reac1.dat',status='old')

      do while (WORD.NE.'NTAU')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (6,*) 'sampling coordinates of reac1 must be defined'
            stop
         endif
      enddo
      read (15,*) ntau1
      rewind 15

c     natomt is to account for dummy atoms
      do while (WORD.NE.'NATOM')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (16,*) 'natom must be defined'
            stop
         endif
      enddo
      read (15,*) natom1,natomt1
      close (15)

      do while (WORD.NE.'ISITE')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (26,*) 'reaction site must be defined'
            stop
         endif
      enddo
      read (25,*) isite,jsite,ksite
      rewind 25

      do while (WORD.NE.'RMAX1')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (26,*) 'grid coords must be defined'
            stop
         endif
      enddo
      read (25,*) rmax1,rstp1,rmax2,rstp2,ireact
      rewind 25

      do while (WORD.NE.'CHARGE')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (26,*) 'charge and spin must be defined'
            stop
         endif
      enddo
      read (25,*) icharge,ispin
      rewind(25)

      do while (WORD.NE.'LEVEL0')
         call LineRead (21)
         if (WORD.EQ.'END') then
            write (26,*) 'level0_ts of theory must be defined'
            write (26,*) 'in file theory.dat'
            stop
         endif
      enddo
      if(word2.eq.'G09') then
         ilev0code=1
      else if (word2.eq.'MOLPRO') then
         ilev0code=2
         commandcopy='cp -f ./data/level0_molpro.dat 
     $                    level0_molpro.dat'            
         call commrun(commandcopy)
      else
         write (26,*) 'level0 of theory must be either'
         write (26,*) 'g09 or molpro in theory.dat'
         stop
      endif
      if(word3.eq.'RESTART') then
         open (unit=99,status='unknown')
         write (99,*) word4
         rewind (99)
         read (99,*) ires
         close (unit=99,status='keep')
      endif
      if(ilev0code.eq.1)then
         read (21,'(A70)') comline1
         read (21,'(A70)') comline2
      endif
      close(21)

      if (idebug.ge.2) write (26,*) ' comline test',comline1,comline2

c build z-mat input

      natom = natom1
      ntau = ntau1
      ismp=0
      iopt=0

      if (idebug.ge.2) write (26,*) 'entering ts_zmat'

      call ts_zmat_beta(atomlabel,xinti,
     $ intcoor,natom1,natomt1,ntau1,iopt,dstart)

      if (idebug.ge.2) write (26,*) 'finished ts_zmat'

c      do iatom = 1 , natomt1
c         atomlabel(iatom) = atomlabel1(iatom)
c      enddo

      vrimax = -1.0d20
      nr1=(rmax1-dstart)/rstp1+0.5
      nr2=(rmax2-rmax1)/rstp2+0.5

      nr=nr1+nr2

      write(26,*)'nr1 is',nr1
      write(26,*)'nr2 is',nr2
      write(26,*)'nr is',nr

      nint = 3*natom-6
      do iint = 1 , nint
         xint(iint) = xinti(iint)
      enddo

      do ir = 1 , nr
         if(ir.le.nr1) then
            rts = dstart+float(ir)*rstp1
         else
            rts = rmax1+float(ir-nr1)*rstp2
         endif

c         write(26,*)'rts is',rts
c      enddo
c      stop
c      do ir = 1 , nr

c initializations
c     write (16,*) ' starting intializations'
         ntau = 0
         natomt = natomt1
c        natomt = natomt1+natomt2+1
         xint(nint) = rts

c call ab initio potential

c        write (6,*) 'starting g09pot'
c        if (natom2.eq.1) xinti(natomtp+1) = rts
c        if (natom2.eq.2) xinti(natomtp+1) = rts
c        if (natom2.ge.3) xinti(natomtp+1) = rts
cc we assume the TS is not linear         
         ilin=0
         ircons=1
         ifreq=0
         ixyz=0
         ired=0
         ilev=0

         if(ilev0code.eq.1) then
            if (idebug.ge.2) write (26,*) 'entering g09fopt'
            call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $           coord,vtot_0,vtotr,freq,ifreq,ilin,ismp,
     $           comline1,comline2,icharge,ispin,ircons,
     $           atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,
     $           ixyz,ired)
            write (26,*) 'finished g09pot'
         else if (ilev0code.eq.2) then
            write (26,*) 'starting molprofopt'         
            call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $           coord,vtotr,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $           atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,
     $           ilev,ispecies)
            write (26,*) 'finished molprofopt'         
         endif

         write (26,*) ir,rts,vtotr
c         stop
c find maximum on grid of rts values
         if (vtotr.gt.1.d10) go to 5000
         if (vtotr.gt.vrimax) then 
            irs = ir
            vrimax = vtotr
            do iint = 1 , nint
               intcoors(iint) = intcoor(iint)
               xints(iint) = xint(iint)
            enddo
         endif
 5000    continue
      enddo
c save key data
      write (26,*) 'grid optimized z-matrix'
      do iatom = 1 , natomt
         write (26,*) atomlabel(iatom)
      enddo
      do iint = 1 , nint
         write (26,*) intcoors(iint),xints(iint)
         intcoor(iint) = intcoors(iint)
         xinti(iint) = xints(iint)
      enddo
c      natomtp = 3*(natom1+natom2)-12
c      aabs1 = xints(natomtp+1)
c      babs1 = xints(natomtp+2)
c      aabs2 = xints(natomtp+3)
c      babs2 = xints(natomtp+4)
c      babs3 = xints(natomtp+5)
c      rts = xints(natomtp+6)

      close (unit=25,status='keep')
      close (unit=26,status='keep')

 1201 format ("./geoms/reac1_"I0.2".xyz")

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine grid_opt_2d

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'


      dimension tau(ntaumx),taumn(ntaumx),taumx(ntaumx),freq(nmdmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim),drea(ntaumx)
      dimension xints(3*natommx)
      dimension tauopt(ntaumx)
      dimension xint2ds(3*natommx)

      character*70 comline1,comline2
      character*60 atomlabel(natommx)
      character*60 atomlabel1(natommx)
      character*30 intcoor(3*natommx)
      character*30 intcoori(3*natommx)
      character*30 intcoors(3*natommx)
      character*20 bislab(ntaumx)
      character*20 intlabel
      character*30 gmem
      character*70 command1
      character*40 filename
      character*4  cjunk
      character*100 commandcopy
      character*30 dih_rot

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

cc parameter initialization

      ires=1
      if (iiso.ne.1)then
         write(26,*)'2d grid option supported only for isomerizations'
         write(26,*)'EStokTP stops here'
         stop
      endif
         
c input data

      open (unit=25,file='./data/ts.dat',status='old')
      open (unit=26,file='./output/grid_opt.out',status='unknown')
      open (unit=21,file='./data/theory.dat',status='unknown')

cc get data from react1 file

      open (unit=15,file='./data/reac1.dat',status='old')

      do while (WORD.NE.'NTAU')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (26,*) 'sampling coordinates of reac1 must be defined'
            stop
         endif
      enddo
      read (15,*) ntau1
      rewind 15

c     natomt is to account for dummy atoms
      do while (WORD.NE.'NATOM')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (26,*) 'natom must be defined'
            stop
         endif
      enddo
      read (15,*) natom1,natomt1
      close (15)

      do while (WORD.NE.'ISITE')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (26,*) 'reaction site must be defined'
            stop
         endif
      enddo
      read (25,*) isite,jsite,ksite
      rewind 25

      do while (WORD.NE.'RMIN1_1')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (26,*) '2d grid coords for site1 must be defined'
            stop
         endif
      enddo
      read (25,*) rmin1_1,rstp1_1,rmin2_1,rstp2_1,ireact_1
      rewind 25

cc initialize ireact_1 to the RTS pulled coordinate
      ireact=ireact_1
cc

      do while (WORD.NE.'RMIN1_2')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (26,*) '2d grid coords for site 2 must be defined'
            stop
         endif
      enddo
      read (25,*) rmin1_2,rstp1_2,rmin2_2,rstp2_2,intlabel
      rewind 25

      do while (WORD.NE.'CHARGE')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (26,*) 'charge and spin must be defined'
            stop
         endif
      enddo
      read (25,*) icharge,ispin
      rewind(25)

      do while (WORD.NE.'LEVEL0')
         call LineRead (21)
         if (WORD.EQ.'END') then
            write (26,*) 'gaussian level0_ts of theory must be defined'
            write (26,*) 'in file theory.dat'
            stop
         endif
      enddo
      if(word2.eq.'G09') then
         ilev0code=1
      else if (word2.eq.'MOLPRO') then
         ilev0code=2
         commandcopy='cp -f ./data/level0_molpro.dat 
     $                    level0_molpro.dat'            
         call commrun(commandcopy)
      else
         write (26,*) 'level0 of theory must be either'
         write (26,*) 'g09 or molpro in theory.dat'
         stop
      endif
      if(word3.eq.'RESTART') then
         open (unit=99,status='unknown')
         write (99,*) word4
         rewind (99)
         read (99,*) ires
         close (unit=99,status='keep')
      endif
      if(ilev0code.eq.1)then
         read (21,'(A70)') comline1
         read (21,'(A70)') comline2
      endif
      close(21)

cc now determine the structure to use in the calculations

cc first determine how many structures were determined with the rotational scan

      command1='ls output/reac1_opt_* > temp.log '
      call commrun(command1)
      command1='grep -c ^ temp.log > num.log'
      call commrun(command1)
      open (unit=99,file='num.log',status='unknown')
      read (99,*) numstruct
      close (unit=99,status='keep')

cc now open the xyz files and search the one with the minimum distance
cc between the reacting atoms
      write(*,*)'NUMSTR IS', numstruct
      do j=1,numstruct
         open (unit=99,status='unknown')
         write(99,1201)j
         rewind (99)
         read (99,'(A40)') filename
         close (99)
         write(*,*) 'filename is ',filename
         open(unit=23,file=filename,status='unknown')
         read (23,*)
         read (23,*)
         do k=1,natom1
            read (23,*)cjunk,coox,cooy,cooz
            if(k.eq.ireact)then
               coox1=coox
               cooy1=cooy
               cooz1=cooz
            endif
            if(k.eq.isite)then
               coox2=coox
               cooy2=cooy
               cooz2=cooz
            endif
         enddo
         close(23)
         drea(j)=sqrt((coox1-coox2)**2+(cooy1-cooy2)**2
     $           +(cooz1-cooz2)**2)
         write(26,*)'irea-isite dist for str ',j,' is ',drea(j)
      enddo

cc now determine the structure with the minimum distance
      iopt=1
      do j=1,numstruct
         if(drea(j).lt.drea(iopt))then
            iopt=j
         endif
      enddo
      write(26,*) 'the structure with the min dist is',iopt

c      stop

      natom = natom1
      ntau = ntau1

c build z-mat input

      rts = 0.d0
      ismp=0
c      iopt=0

      if (idebug.ge.2) write (26,*) 'entering ts_zmat'

      call ts_zmat_iso(atomlabel,atomlabel1,xinti,
     $ intcoor,natom1,natomt1,ntau1,iopt,dstart)

      if (idebug.ge.2) write (26,*) 'finished ts_zmat'

      do iatom = 1 , natomt1
         atomlabel(iatom) = atomlabel1(iatom)
c         write(*,*)'atom label is',atomlabel(iatom)
      enddo

      nint = 3*natom-6
      do iint = 1 , nint
         xint(iint) = xinti(iint)
      enddo
c
c     now move the second coordinate to the end of the z-matrix
c

      write(*,*)'int lab is ',intlabel
c      stop

      do iatom = 1 , ncoord
         xinti(index)=0.
         intcoori(index)=''
      enddo
      index=0
      ncoord = 3*natom-6
      do iatom = 1 , ncoord-ntau-1
         if (intcoor(iatom).ne.intlabel) then
            index=index+1
            xinti(index)=xint(iatom)
            intcoori(index)=intcoor(iatom)
         else
c                  xinti(ncoord)=xint(iatom)
            xinti(ncoord-1)=xint(iatom)
            intcoori(ncoord-1)=intlabel              
            write(*,*)'int nc-1 lab is ',intlabel
         endif
      enddo
      index=0
      do iatom = ncoord-1-ntau,ncoord-2
         index=index+1
         xinti(iatom)=tau(index)
         intcoori(iatom)=bislab(index)
      enddo
      intcoori(ncoord)='RTS'
      do iint = 1 , ncoord
         xint(iint) = xinti(iint)
         intcoor(iint)=intcoori(iint)
         write(*,*)'coord and value',xint(iint),intcoor(iint)
      enddo
      
      vrimax = -1.0d20
c      nr2=(dstart-rmin2)/rstp2+0.5
c      nr1=(rmin2-rmin1)/rstp1+0.5

      nr2_1=(dstart-rmin1_1)/rstp1_1+0.5
      nr1_1=(rmin1_1-rmin2_1)/rstp2_1+0.5

      nr_1=nr1_1+nr2_1

      write(26,*)'dstart is',dstart
      write(26,*)'nr1_1 is',nr1_1
      write(26,*)'nr2_1 is',nr2_1
      write(26,*)'nr_1 is',nr_1

      nr2_2=abs((xint(ncoord-1)-rmin1_2)/rstp1_2)+0.5
      nr1_2=abs((rmin1_2-rmin2_2)/rstp2_2)+0.5

      nr_2=nr1_1+nr2_1

      write(26,*)'nr1_2 is',nr1_2
      write(26,*)'nr2_2 is',nr2_2
      write(26,*)'nr_2 is',nr_2

c      stop

cc determine dihedral of reacting atom
      
      if(ireact.gt.3)then
         open (unit=99,status='unknown')
         write(99,*)atomlabel(ireact)
         rewind(99)
         read(99,*)cjunk,cjunk,cjunk,cjunk,cjunk,cjunk,dih_rot
         close(99)
c      write(*,*)'dname is ',dih_rot

cc update dihed variable to avoid planar guess

         do iint = 1 , nint
            if(dih_rot.eq.intcoor(iint))then
               check1=abs(xint(iint)-180.)
               check2=abs(xint(iint)-0.)
c            write(*,*)'check1 is ',check1
c            write(*,*)'check2 is ',check2
               if(check1.lt.1.)then
                  xint(iint)=xint(iint)+10.
                  check1=10.
               endif
               if(check2.lt.1.)then
                  xint(iint)=xint(iint)+10.
                  check2=10.
               endif
            endif
         enddo
      else if (ireact.eq.3) then
         open (unit=99,status='unknown')
         write(99,*)atomlabel(ireact+1)
         rewind(99)
         read(99,*)cjunk,cjunk,cjunk,cjunk,cjunk,cjunk,dih_rot
         rewind(99)
         write(99,*)dih_rot
         rewind(99)
         call lineread(99)
         dih_rot=word
         
c         read(99)dih_rot
c         read(99,*)cjunk,cjunk,cjunk,cjunk,cjunk,cjunk,dih_rot
         close(99)
c      write(*,*)'dname is ',dih_rot
c      stop
cc update dihed variable to avoid planar guess

         do iint = 1 , nint
            if(dih_rot.eq.intcoor(iint))then
               check1=abs(xint(iint)-180.)
               check2=abs(xint(iint)-0.)
c            write(*,*)'check1 is ',check1
c            write(*,*)'check2 is ',check2
               if(check1.lt.1.)then
                  xint(iint)=xint(iint)+10.
                  check1=10.
               endif
               if(check2.lt.1.)then
                  xint(iint)=xint(iint)+10.
                  check2=10.
               endif
            endif
         enddo
      endif
c     stop


cc start 2D cycling
c            write(*,*)'ok 1'
c            stop


      do iint=1,nint
         xint2ds(iint)=0.
      enddo
      xint2start=xint(nint-1)

      do ir_1 = 1 , nr_1

         ir2_=0
         if(ir_2.ne.0.and.ir_1.ne.1)then
            do iint = 1 , nint
               xint(iint) = xint2ds(iint)
            enddo
         endif

         if(ir_1.le.nr2_1) then
c            rts = dstart-float(ir)*rstp2
            rts = dstart-float(ir_1)*rstp1_1
         else
c            rts = rmin2-float(ir-nr2)*rstp1
            rts = rmin1_1-float(ir_1-nr2_1)*rstp2_1
         endif

         if(ir_1.le.nr2_1) then
            nr_2_cyc=1
         else
            nr_2_cyc=nr_2
         endif
         do ir_2 = 1 , nr_2_cyc
            if(ir_2.le.nr_2) then
               if(rmin2_2.lt.rmin1_2)then
                  xint(nint-1) = xint2start-float(ir_2)*rstp1_2
               else
                  xint(nint-1) = xint2start+float(ir_2)*rstp1_2
               endif
            else
               if(rmin2_2.lt.rmin1_2)then
                  xint(nint-1) = rmin1_2-float(ir_2-nr2_2)*rstp2_2
               else
                  xint(nint-1) = rmin1_2+float(ir_2-nr2_2)*rstp2_2
               endif
            endif

c            write(*,*)'ok 2'
c      enddo
c            write (26,*) ir_1,ir_2,rts,xint(nint-1),vtotr

c            goto 5000
c      do ir = 1 , nr

c initializations
c     write (16,*) ' starting intializations'
            ntau = 0
            natomt = natomt1
c        natomt = natomt1+natomt2+1
            xint(nint) = rts

cc update dihed variable to avoid planar guess

            do iint = 1 , nint
               if(dih_rot.eq.intcoor(iint))then
                  check1=abs(xint(iint)-180.)
                  check2=abs(xint(iint)-0.)
c            write(*,*)'check1 is ',check1
c            write(*,*)'check2 is ',check2
                  if(check1.lt.1.)then
                     xint(iint)=xint(iint)+10.
                     check1=10.
                  endif
                  if(check2.lt.1.)then
                     xint(iint)=xint(iint)+10.
                     check2=10.
                  endif
               endif
            enddo

c call ab initio potential

cc we assume the TS is not linear         
            ilin=0
            ircons=2
            ifreq=0
            ixyz=0
            ired=0
            ilev=0

            if(ilev0code.eq.1) then
c               if (idebug.ge.2) write (26,*) 'entering g09fopt'
               call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $           coord,vtot_0,vtotr,freq,ifreq,ilin,ismp,
     $           comline1,comline2,icharge,ispin,ircons,
     $           atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires
     $           ,ixyz,ired)
c               write (26,*) 'finished g09pot'
            else if (ilev0code.eq.2) then
               write (26,*) 'starting molprofopt'         
               call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $           coord,vtotr,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $           atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,
     $           ilev,ispecies)
               write (26,*) 'finished molprofopt'         
            endif

            write (26,*) ir_1,ir_2,rts,xint(nint-1),vtotr

c
cc save intermediate data
c
            if(ir_2.eq.1)then
               do iint = 1 , nint
                  xint2ds(iint) = xint(iint)
               enddo
            endif

c         stop
c find maximum on grid of rts values
            if (vtotr.gt.1.d10) go to 5000
            if (vtotr.gt.vrimax) then 
               irs = ir_1
               vrimax = vtotr
               do iint = 1 , nint
                  intcoors(iint) = intcoor(iint)
                  xints(iint) = xint(iint)
               enddo
            endif

 5000       continue

         enddo
      enddo


c save key data
      write (26,*) 'grid optimized z-matrix'
      do iatom = 1 , natomt
         write (26,*) atomlabel(iatom)
      enddo
      do iint = 1 , nint
         write (26,*) intcoors(iint),xints(iint)
         intcoor(iint) = intcoors(iint)
         xinti(iint) = xints(iint)
      enddo
c      natomtp = 3*(natom1+natom2)-12
c      aabs1 = xints(natomtp+1)
c      babs1 = xints(natomtp+2)
c      aabs2 = xints(natomtp+3)
c      babs2 = xints(natomtp+4)
c      babs3 = xints(natomtp+5)
c      rts = xints(natomtp+6)

      close (unit=25,status='keep')
      close (unit=26,status='keep')

 1201 format ("./geoms/reac1_"I0.2".xyz")

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ts_0

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension tau(ntaumx),taumn(ntaumx),taumx(ntaumx),freq(nmdmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim)
      dimension tauopt(ntaumx)

      character*70 comline1,comline2
      character*60 atomlabel(natommx)
      character*30 intcoor(3*natommx)
      character*20 bislab(ntaumx)
      character*30 gmem
      character*100 commandcopy

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

cc initialize parameters
      ires=0
      ilin1=0
      ilin2=0
      ilin=0
c input data

      open (unit=35,file='./data/ts.dat',status='old')
      open (unit=36,file='./output/ts_opt_0.out',status='unknown')

      open (unit=15,file='./data/reac1.dat',status='old')
      open (unit=21,file='./data/theory.dat',status='unknown')

      do while (WORD.NE.'NTAU')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (6,*) 'sampling coordinates of reac1 must be defined'
            stop
         endif
      enddo
      read (15,*) ntau1
      rewind 15

c     natomt is to account for dummy atoms
      do while (WORD.NE.'NATOM')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (16,*) 'natom must be defined'
            stop
         endif
      enddo
      read (15,*) natom1,natomt1,ilin1
      close (15)

cc get data from react2 file
      if(iadd.eq.1.or.iabs.eq.1) then

         open (unit=15,file='./data/reac2.dat',status='old')

         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (6,*) 'sampling coordinates of reac2 
     $                      must be defined'
               stop
            endif
         enddo
         read (15,*) ntau2
         rewind 15

c     natomt is to account for dummy atoms
         do while (WORD.NE.'NATOM')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'natom must be defined'
               stop
            endif
         enddo
         read (15,*) natom2,natomt2,ilin2
c         write(*,*)'natom2 is',natom2
c         write(*,*)'natomt2 is',natomt2
c         stop
         close (15)
      endif


c which atom to abstract and atoms to connect to it
      if (iadd.eq.1.or.iabs.eq.1) natom = natom1+natom2
      if (iiso.eq.1.or.ibeta.eq.1) natom = natom1
      if (iadd.eq.1) natomt = natomt1+natomt2
      if (iabs.eq.1) natomt = natomt1+natomt2+1
      if (iiso.eq.1.or.ibeta.eq.1) natomt = natomt1
c     natomt = natomt1+natomt2+1
      if (iadd.eq.1.or.iabs.eq.1) then
         ntau = ntau1+ntau2
      else if (iiso.eq.1.or.ibeta.eq.1) then
         ntau = ntau1
      endif
c      write(*,*)'natomt1 is',natomt1
c      write(*,*)'natomt2 is',natomt2
c      write(*,*)'natomt is',natomt
c      stop
cc if linear molecule reacting with atom in abs assume TS is linear
         if(iabs.eq.1.and.natom2.eq.1.
     $      and.ilin1.eq.1) then
            ilin=1
         endif

cc if linear molecule reacting with linear radical in abs assume TS is linear
         if(iabs.eq.1.and.ilin2.eq.1.and
     $      .ilin1.eq.1)then
            ilin=1
         endif


c input template data

      do while (WORD.NE.'NOSMP')
         call LineRead (35)
         if (WORD.EQ.'END') then
            write (6,*) 'number of sampling points must be defined'
            stop
         endif
      enddo
      read (35,*) nosmp,dthresh,ethresh
      rewind(35)
c
      do while (WORD.NE.'CHARGE')
         call LineRead (35)
         if (WORD.EQ.'END') then
            write (36,*) 'charge and spin must be defined'
            stop
         endif
      enddo
      read (35,*) icharge,ispin
      rewind(35)
c
      do while (WORD.NE.'LEVEL0_TS')
         call LineRead (21)
         if (WORD.EQ.'END') then
            write (36,*) ' level0 of theory must be defined for ts'
            write (36,*) 'in file theory.dat'
            stop
         endif
      enddo
      if(word2.eq.'G09') then
         ilev0code=1
      else if (word2.eq.'MOLPRO') then
         ilev0code=2
         commandcopy='cp -f ./data/level0_ts_molpro.dat 
     $                    level0_molpro.dat'            
         call commrun(commandcopy)
      else
         write (36,*) 'level0 of theory must be either'
         write (36,*) 'g09 or molpro in theory.dat'
         stop
      endif
      if(word3.eq.'RESTART') then
         open (unit=99,status='unknown')
         write (99,*) word4
         rewind (99)
         read (99,*) ires
         close (unit=99,status='keep')
      endif
      if(ilev0code.eq.1)then
         read (21,'(A70)') comline1
         read (21,'(A70)') comline2
      endif
      close(21)
      if (idebug.ge.2) write (36,*) ' comline test',comline1,comline2

      open (unit=26,file='./output/grid_opt.out',status='unknown')
 1000 continue
      CALL LineRead (26)
      if (WORD.eq.'GRID') go to 1500
      go to 1000
 1500 continue 
      do iatom = 1 , natomt
         read (26,'(A60)') atomlabel(iatom)
c        write (6,*) 'atomlabel',atomlabel(iatom)
      enddo 
      nint = 3*natom-6
      do iint = 1 , nint
         read (26,*) intcoor(iint),xint(iint)
c         write (6,*) 'intcoor, xint ',intcoor(iint),xint(iint)
      enddo
      close (unit=26,status='keep')
c      stop

c     do iint = 1 , nint
c        xint(iint) = xinti(iint)
c     enddo
      ntau = 0
      ircons=0
      ifreq=0
      ismp=0
      ixyz=0
      ired=0
      ilev=0
      if(ifrozrts.eq.1) then 
         ircons=1
         xint(nint)=frozcoo
      endif

c         write(*,*)'ilin is',ilin
c         stop

      if(ilev0code.eq.1) then
         call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $        coord,vtot_0,vtotr,freq,ifreq,ilin,ismp,
     $        comline1,comline2,icharge,ispin,ircons,
     $        atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires
     $        ,ixyz,ired)

      else if (ilev0code.eq.2) then
         call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $        coord,vtotr,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $        atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,
     $        ilev,ispecies)
      endif

cc added to save optimized TS z-matrix coordinates

      write (36,*) 'energy '
      write (36,*) vtotr
      write (36,*) 'TS_z-matrix'
      do iatom = 1 , natomt
         write (36,'(A60)') atomlabel(iatom)
      enddo
      do iint = 1 , nint
         write (36,*) intcoor(iint),xint(iint)
      enddo

      close (unit=35,status='keep')
      close (unit=36,status='keep')

cc added to have a file with highest level TS z-matrix parameters
cc this is updated by the higher level TS search and random torsional scan subs

      open (unit=37,file='./output/ts_opt.out',status='unknown')

      write (37,*) 'TS_z-matrix'
      do iatom = 1 , natomt
         write (37,'(A60)') atomlabel(iatom)
      enddo
      do iint = 1 , nint
         write (37,*) intcoor(iint),xint(iint)
      enddo
      close (unit=37,status='keep')

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tauo_ts

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'


      dimension tau(ntaumx),taumn(ntaumx),taumx(ntaumx),freq(nmdmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim),xint_opt0(ntaumx)
      dimension tauopt(ntaumx)
      dimension vref(noptmx),tauo(ntaumx,noptmx),
     $ xinto(3*natommx,noptmx),dist(noptmx)
      dimension rtsopt(noptmx),bbopt(noptmx)

      real *8 mtrandom

      character*70 comline1,comline2
      character*60 atomlabel(natommx)
      character*30 intcoor(3*natommx)
      character*30 intcoori(3*natommx)
      character*20 bislab(ntaumx)
      character*70 command1
      character*30 gmem,distcheck,bbcheck
      character*100 commandcopy
      character*40 filename
      character*40 cjunk

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

      noptg = 0
      ires = 0
      ilin = 0
      ilin1=0
      ilin2=0
cc      call mtinit()

c     open (unit=55,file='tauo_ts.dat',status='old')
      open (unit=55,file='./data/ts.dat',status='old')
      open (unit=56,file='./output/tauo_ts.out',status='unknown')
      open (unit=21,file='./data/theory.dat',status='unknown')

c input data
      do while (WORD.NE.'NOSMP')
         call LineRead (55)
         if (WORD.EQ.'END') then
            write (56,*) 'number of sampling points must be defined'
            stop
         endif
      enddo
      read (55,*) nosmp,dthresh,ethresh
      rewind(55)

c determine number of reacting number for beta-scission reactions
      if(ibeta.eq.1)then
         do while (WORD.NE.'RMAX1')
            call LineRead (55)
            if (WORD.EQ.'END') then
               write (56,*) 'grid coords must be defined'
               stop
            endif
         enddo
         read (55,*) rmax1,rstp1,rmax2,rstp2,ireact
         rewind(55)
      endif

c
      open (unit=15,file='./data/reac1.dat',status='old')

      do while (WORD.NE.'NTAU')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (56,*) 'sampling coordinates of reac1 must be defined'
            stop
         endif
      enddo
      read (15,*) ntau1
      rewind 15

c     natomt is to account for dummy atoms
      do while (WORD.NE.'NATOM')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (56,*) 'natom must be defined'
            stop
         endif
      enddo
      read (15,*) natom1,natomt1,ilin1
      close(15)

cc get data from react2 file

      if (iadd.eq.1.or.iabs.eq.1)then
         open (unit=15,file='./data/reac2.dat',status='old')

         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (56,*) 'sampling coord of reac2 must be defined'
               stop
            endif
         enddo
         read (15,*) ntau2
         rewind 15

c     natomt is to account for dummy atoms
         do while (WORD.NE.'NATOM')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (56,*) 'natom must be defined'
               stop
            endif
         enddo
         read (15,*) natom2,natomt2,ilin2
         close (15)
      endif

c total number of atoms initialization
      if (iadd.eq.1.or.iabs.eq.1)then
         natom = natom1+natom2
      else if (iiso.eq.1.or.ibeta.eq.1)then
         natom = natom1
      endif
      if (iadd.eq.1) natomt = natomt1+natomt2
      if (iabs.eq.1) natomt = natomt1+natomt2+1
      if (iiso.eq.1) natomt = natomt1
      if (ibeta.eq.1) natomt = natomt1

cc if linear molecule reacting with atom in abs assume TS is linear
         if(iabs.eq.1.and.natom2.eq.1.
     $      and.ilin1.eq.1) then
            ilin=1
         endif

cc if linear molecule reacting with linear radical in abs assume TS is linear
         if(iabs.eq.1.and.ilin2.eq.1.and
     $      .ilin1.eq.1)then
            ilin=1
         endif


c     natomt = natomt1+natomt2+1
cc get isite number from ts file
      do while (WORD.NE.'ISITE')
         call LineRead (55)
         if (WORD.EQ.'END') then
            write (56,*) 'reaction site must be defined'
            stop
         endif
      enddo
      read (55,*) isite,jsite,ksite
      rewind 55

c number and name of sampled coordinates

      do while (WORD.NE.'NTAU')
         call LineRead (55)
         if (WORD.EQ.'END') then
            write (6,*) 'sampling coordinates must be defined'
            stop
         endif
      enddo
      read (55,*) ntau
      if (ntau.gt.ntaumx) then
         write (56,*) 'ntau too large',ntau,ntaumx
         stop
      endif
      write(56,*)'sampled coordinates'
      if (ntau.ne.0) then
         read (55,*)
         do itau = 1 , ntau
            read (55,*) bislab(itau),taumn(itau),taumx(itau)
            write (56,*) bislab(itau),taumn(itau),taumx(itau)
            open (unit=99,status='unknown')
            rewind (99)
            write (99,*)bislab(itau)
            rewind (99)
            call LineRead (99)
            bislab(itau)=WORD
            close (unit=99,status='keep')
         enddo
      endif
      rewind(55)

c input template data

      do while (WORD.NE.'CHARGE')
         call LineRead (55)
         if (WORD.EQ.'END') then
            write (56,*) 'charge and spin must be defined'
            stop
         endif
      enddo
      read (55,*) icharge,ispin
      rewind(55)

      do while (WORD.NE.'LEVEL0_TS')
         call LineRead (21)
         if (WORD.EQ.'END') then
            write (56,*) 'gaussian level0 of theory must be defined'
            stop
         endif
      enddo
      if(word2.eq.'G09') then
         ilev0code=1
      else if (word2.eq.'MOLPRO') then
         ilev0code=2
         commandcopy='cp -f ./data/level0_ts_molpro.dat 
     $                    level0_molpro.dat'            
         call commrun(commandcopy)
      else
         write (56,*) 'level0 of theory must be either'
         write (56,*) 'g09 or molpro in theory.dat'
         stop
      endif
      if(word3.eq.'RESTART') then
         open (unit=99,status='unknown')
         write (99,*) word4
         rewind (99)
         read (99,*) ires
         close (unit=99,status='keep')
      endif
      if(ilev0code.eq.1)then
         read (21,'(A70)') comline1
         read (21,'(A70)') comline2
      endif
      close(21)
      write(56,*)
      if (idebug.ge.2) write (56,*) ' command lines',comline1,comline2

      open (unit=36,file='./output/ts_opt.out',status='unknown')
 1000 continue
      CALL LineRead (36)
      if (WORD.eq.'TS_Z-MATRIX') go to 1500
      write (6,*) 'WORD',WORD
      go to 1000
 1500 continue 
      do iatom = 1 , natomt
         read (36,'(A60)') atomlabel(iatom)
         write (6,*) 'atomlabel',atomlabel(iatom)
      enddo 

cc for abstraction determine label of breaking bond and initialize bbcheck
      if(iabs.eq.1)then
         open(unit=99,status='unknown')
         do iatom = 1 , natomt
            rewind (99)
            write(99,*)atomlabel(iatom)
            rewind (99)
            call LineRead(99)
            if(iatom.eq.isite)then
               bbcheck=word3
            endif
         enddo
         write(56,*)'bbcheck is ', bbcheck
c         write(*,*)'bbcheck is ', bbcheck
      endif
c      stop


c for beta-decomp determine label of breaking bond      
      if(ibeta.eq.1)then
         open(unit=99,status='unknown')
         do iatom = 1 , natomt
            rewind (99)
            write(99,*)atomlabel(iatom)
            rewind (99)
            call LineRead(99)
            if(iatom.eq.ireact)then
               distcheck=word3
            endif
         enddo
         close(99)
         write(56,*)'distcheck is ', distcheck
      endif
c      stop
      nint = 3*natom-6
cc
      rtscheck=0.
      bbdischeck=0.
      do iint = 1 , nint
         read (36,*) intcoori(iint),xinti(iint)
         write (6,*) 'intcoor, xint ',intcoori(iint),xinti(iint)
         if(intcoori(iint).eq.'RTS') then
            rtscheck=xinti(iint)
         endif
         if(ibeta.eq.1) then
            if(intcoori(iint).eq.distcheck) then
               rtscheck=xinti(iint)
            endif
         endif
         if(iabs.eq.1) then
            if(intcoori(iint).eq.bbcheck) then
               bbdischeck=xinti(iint)
            endif
         endif

      enddo
      close (unit=36,status='keep')
      if(irecov.ne.1)write(56,*)'rtscheck is ',rtscheck
      if(iabs.eq.1)write(56,*)'bbcheck is ',bbdischeck

c      write(*,*)'rtscheck is ',rtscheck
c      if(iabs.eq.1)write(*,*)'bbcheck is ',bbdischeck

c      stop
cc    re-order z-matrix, taking out torsions      

      write (6,*) 're-ordered z-matrix of TS'

      index=0
      ncoord = 3*natom-6
      do iatom = 1 , ncoord
         do itau=1,ntau
            if (intcoori(iatom).eq.bislab(itau)) then
               xint_opt0(itau)=xinti(iatom)
               goto 998
               endif
c            if (intcoori(iatom).eq.bislab(itau)) goto 998
         enddo
         index=index+1
         xint(index)=xinti(iatom)
         intcoor(index)=intcoori(iatom)
         write (6,*) 'intcoor, xint ',intcoor(index),xint(index)
 998     continue
      enddo

      nint = 3*natom-ntau-6

cc re-initialize initial coordinates

      do iint = 1 , nint
         xinti(iint) = xint(iint)
      enddo

c start sampling

      noptg=0
      ntau_save=ntau

      if(irecov.ne.1)then
         do ismp = 1 , nosmp


            write (56,*) '*******************************'
            write (56,*) 'Random number point ',ismp,' out  of ',nosmp

c generate random tau

            do itau = 1 , ntau
               tau(itau) = taumn(itau)+(taumx(itau)-taumn(itau))*
     $              mtrandom()
            enddo
            write(56,*)'Sampled Tau coordinates'
            if(ismp.ne.1) then
               do itau = 1 , ntau
                  write(56,*)bislab(itau),tau(itau)
               enddo
            else
               do itau = 1 , ntau
                  write(56,*)bislab(itau),xint_opt0(itau)
               enddo
            endif

cc for the first point start from level0

            if(ismp.eq.1) then
               do itau = 1 , ntau
                  tau(itau) =xint_opt0(itau)
               enddo
            endif

c initialize TS coordinates to structure optimized in sub ts_opt

            do iint = 1 , nint
               xint(iint) = xinti(iint)
            enddo

            ircons=0
            ifreq=0
            if(ifrozrts.eq.1) then
cc if using a frozen coordinate, then
cc move rts to last position and freeze it
               do iatom = 1 , ncoord
                  xinti(index)=0.
                  intcoori(index)=''
               enddo
               index=0
               ncoord = 3*natom-6
               do iatom = 1 , ncoord-ntau
                  if (intcoor(iatom).ne.'RTS') then
                     index=index+1
                     xinti(index)=xint(iatom)
                     intcoori(index)=intcoor(iatom)
                  else
c                  xinti(ncoord)=xint(iatom)
                     xinti(ncoord)=frozcoo
                     intcoori(ncoord)=intcoor(iatom)                  
                  endif
               enddo
               index=0
               do iatom = ncoord-ntau,ncoord-1
                  index=index+1
                  xinti(iatom)=tau(index)
                  intcoori(iatom)=bislab(index)
               enddo

               do iint = 1 , ncoord
                  xint(iint) = xinti(iint)
                  intcoor(iint)=intcoori(iint)
               enddo
               ircons=1
               ntau=0
c            write(6,*)'tau sampling for TS not working '
c            write(6,*)'with frozen coordinate option'
c            stop
            endif
            ixyz=0
            ired=0
            ilev=0
            ispecies=0

            if(ilev0code.eq.1) then
               write (56,*) 'starting g09pot'
               call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $           coord,vtot_0,vtot,freq,ifreq,ilin,ismp,
     $           comline1,comline2,icharge,ispin,ircons,
     $           atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires
     $           ,ixyz,ired)
               write (56,*) 'finished g09pot'

            else if (ilev0code.eq.2) then
               write (56,*) 'starting molprofopt'         
               call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $           coord,vtot,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $           atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,
     $           ilev,ispecies)
               write (56,*) 'finished molprofopt'         
            endif

            write (56,*) 
          write (56,*) 'optimized tau coordinates,energy,sampling point'
            write (56,*) (tauopt(itau),itau=1,ntau),vtot,ismp
            write (56,*) 'optimized coordinates'
            write (56,1111) (xint(it),it=1,nint)
 1111       format (100(1X,f7.2))

            if(ifrozrts.eq.1) ntau=ntau_save
            if (vtot.gt.1.d10) go to 5000

c determine if new geometry
            dmin = 1.d30
            do iopt = 1 , noptg
               dist(iopt) = 0.0d0
               if(ifrozrts.eq.0)then
                  do itau = 1 , ntau
                     dref = abs(tauopt(itau)-tauo(itau,iopt))
                     dplus = abs(tauopt(itau)+360.d0-tauo(itau,iopt))
                     dminus = abs(tauopt(itau)-360.d0-tauo(itau,iopt))
                     dd = min(dref,dplus,dminus)**2
                     dist(iopt) = dist(iopt) + dd
                  enddo
               else if(ifrozrts.eq.1)then
                  icoofr=ncoord-1-ntau
                  do itau = icoofr , ncoord-1
                     dref = abs(xint(itau)-xinto(itau,iopt))
                     dplus = abs(xint(itau)+360.d0-xinto(itau,iopt))
                     dminus = abs(xint(itau)-360.d0-xinto(itau,iopt))
                     dd = min(dref,dplus,dminus)**2
                     dist(iopt) = dist(iopt) + dd
                  enddo
               endif
               dist(iopt)=sqrt(dist(iopt))
c           write (6,*) 'iopt test'
               if (idebug.ge.2) write (56,*) dist(iopt),iopt,dthresh
               if (dist(iopt).lt.dthresh) go to 5000
            enddo
c enter this point only if it is a new geometry
c save tau, internal coordinates and vibrational frequencies for each new geometry
            noptg = noptg+1
            vref(noptg) = vtot
            if(ifrozrts.ne.1) then
               do itau = 1 , ntau
                  tauo(itau,noptg) = tauopt(itau)
               enddo 
               do iint = 1 , nint
                  xinto(iint,noptg) = xint(iint)
                  if(intcoor(iint).eq.'RTS') then
                     rtsopt(noptg)=xint(iint)
                     write (56,*) 'rtsopt is : ',rtsopt(noptg)
                  endif
                  if(ibeta.eq.1) then
                     if(intcoor(iint).eq.distcheck) then
                        rtsopt(noptg)=xint(iint)
                        write (56,*) 'rtsopt is : ',rtsopt(noptg)
                     endif
                  endif
                  if(iabs.eq.1) then
                     if(intcoor(iint).eq.bbcheck) then
                        bbopt(noptg)=xint(iint)
                        write (56,*) 'bbopt is : ',bbopt(noptg)
                     endif
                  endif
               enddo
               write (56,*) 'new geom at nopt : ',noptg
               write (56,*) 'tau ' 
               write (56,*) (tauopt(itau),itau=1,ntau)
               write (56,*) 'internal coordinates ' 
               do iint = 1 , nint
                  write (56,*) intcoor(iint),xint(iint)
               enddo
               write (56,*) 'coordinates '
               do iatom = 1 ,natomt
                  write (56,*) (coord(iatom,idim),idim=1,ndim)
               enddo
               write (56,*) 'rotational constants'
               write (56,*) (abcrot(idim),idim=1,ndim)
               write (56,*) 'energy, energy diff with structure 1 '
               write (56,*) vtot, (vref(iopt)-vref(1))
            else
               do iint = 1 , ncoord
                  xinto(iint,noptg) = xint(iint)
                  if(intcoor(iint).eq.'RTS') then
                     rtsopt(noptg)=xint(iint)
                     write (56,*) 'rtsopt is : ',rtsopt(noptg)
                  endif
                  if(ibeta.eq.1) then
                     if(intcoor(iint).eq.distcheck) then
                        rtsopt(noptg)=xint(iint)
                        write (56,*) 'rtsopt is : ',rtsopt(noptg)
                     endif
                  endif
                  if(iabs.eq.1) then
                     if(intcoor(iint).eq.bbcheck) then
                        bbopt(noptg)=xint(iint)
                        write (56,*) 'bbopt is : ',bbopt(noptg)
                     endif
                  endif
               enddo
               write (56,*) 'new geom at nopt : ',noptg
               write (56,*) 'tau ' 
               write (56,*) (tauopt(itau),itau=1,ntau)
               write (56,*) 'internal coordinates ' 
               do iint = 1 , ncoord
                  write (56,*) intcoor(iint),xint(iint)
               enddo
               write (56,*) 'coordinates '
               do iatom = 1 ,natomt
                  write (56,*) (coord(iatom,idim),idim=1,ndim)
               enddo
               write (56,*) 'rotational constants'
               write (56,*) (abcrot(idim),idim=1,ndim)
               write (56,*) 'energy, energy diff with structure 1 '
               write (56,*) vtot, (vref(iopt)-vref(1))*627.5
               write (56,*) 'rtsopt is', rtsopt(noptg)
               if(iabs.eq.1) write (56,*) 'bbopt is', bbopt(noptg)
            endif


c         command1='newzmat       
c     +           -ichk -oxyz tmp.chk geom.xyz           '
c         call commrun(command1)
            open (unit=98,file='headgeom.tmp',status='unknown')
            write(98,*)natom
            write(98,*)vtot
            close(98)
            command1='cat headgeom.tmp geom.xyz > geom1.xyz '
            call commrun(command1)

            open (unit=99,status='unknown')
            rewind (99)
            write (99,1112) noptg
            rewind (99)
            read (99,1020) command1
            close (99)
            call commrun(command1)         
 1112       format ('cp -f geom1.xyz geoms/ts_',I0.2,'.xyz')

            open (unit=98,file='tsgeom.tmp',status='unknown')
            write (98,*) 'TS_z-matrix'
            do iatom = 1 , natomt
               write (98,'(A60)') atomlabel(iatom)
            enddo
            if(ifrozrts.ne.1)then
               do iint = 1 , nint
                  write (98,*) intcoor(iint),xinto(iint,noptg)
               enddo
               do itau=1,ntau
                  write (98,*) bislab(itau),tauo(itau,noptg)
               enddo
            else
               do iint = 1 , ncoord
                  write (98,*) intcoor(iint),xinto(iint,noptg)
               enddo
            endif
            close(98)
            open (unit=99,status='unknown')
            rewind (99)
            write (99,1114) noptg
            rewind (99)
            read (99,1020) command1
            close (99)
            call commrun(command1)         
 1114       format ('cp -f tsgeom.tmp output/ts_opt_',I0.2,'.out')

            if(idebug.ge.2)then
               open (unit=99,status='unknown')
               rewind (99)
               write (99,1113) noptg
               rewind (99)
               read (99,1020) command1
               close (99)
               call commrun(command1)         
 1113          format ('cp -f geom.log geoms/ts_',I0.2,'.log')
            endif

 1020       format (A70)
            
            continue
c check to see if this is also a new energy
            inew = 1
            if (noptg.gt.1) then
               do iopt = 1 , noptg-1
                  dv = abs(vref(iopt)-vref(noptg))
                  if (dv.lt.ethresh) then
                     inew = 0
                  endif
               enddo 
            endif
            if (inew.eq.1) then 
c            do itau = 1 , ntau
c               tauo2(itau,noptg) = tauo(itau,noptg)
c            enddo 
c            do iint = 1 , nint
c               xinto2(iint,noptg) = xinto(iint,noptg)
c            enddo
               write (56,*) 'new energy also'
            endif
 5000       continue
         enddo
      else

cc recover results from previous calculation

cc first determine how many structures were determined with the rotational scan

         command1='ls output/ts_opt_*.out > temp.log '
         call commrun(command1)
         command1='grep -c ^ temp.log > num.log'
         call commrun(command1)
         open (unit=99,file='num.log',status='unknown')
         read (99,*) numstruct
         close (unit=99,status='keep')
c         numstruct=numstruct-1
         write(*,*)'num struct to choose from is: ',numstruct

cc now open the xyz files and get the energies
c
         write(*,*)'NUMSTR IS', numstruct
         do j=1,numstruct
            open (unit=99,status='unknown')
            write(99,1201)j
            rewind (99)
            read (99,'(A40)') filename
            close (99)
c            write(*,*) 'filename is ',filename
            open(unit=23,file=filename,status='unknown')
            read (23,*)
            read (23,*)vref(j)
            close(23)
            write(*,*)'the energy of stru ',j,'is ',vref(j)
         enddo

         do j=1,numstruct
            open (unit=99,status='unknown')
            write(99,1202)j
            rewind (99)
            read (99,'(A40)') filename
            close (99)
c            write(*,*) 'filename is ',filename
            open(unit=23,file=filename,status='unknown')
            read (23,*)
            do ik=1,natomt
               read (23,*)
            enddo
            if(ifrozrts.ne.1)then
               do iint = 1 , nint
                  read (23,*)cjunk, xinto(iint,j)
               enddo
               do itau=1,ntau
                  read (23,*)cjunk, tauo(itau,j)
               enddo
            else
               do iint = 1 , ncoord
                  read (23,*)cjunk, xinto(iint,j)
               enddo
            endif
            close(23)
            do iint = 1 , ncoord
               if(intcoor(iint).eq.'RTS') then
                  rtsopt(j)=xinto(iint,j)
               endif
               if(iabs.eq.1)then
                  if(intcoor(iint).eq.bbcheck) then
                     bbopt(j)=xinto(iint,j)
                  endif
               endif
            enddo
         write (*,*) 'rtsopt of stru j is : ',rtsopt(j)

         enddo
cc now determine rts coord for each structure
         rtscheck=rtsopt(1)
         if(iabs.eq.1)bbdischeck=bbopt(1)
         write(56,*)'rtscheck is ',rtscheck
      endif

 1201 format ("./geoms/ts_"I0.2".xyz")
 1202 format ("./output/ts_opt_"I0.2".out")

c      stop


c for the optimum geometry print out the optimized coordinates
      if(irecov.eq.1)noptg=numstruct
      ioptmin = 1
      voptmin= vref(ioptmin)
c     write (16,*) 'opt test',noptg,voptmin,ioptmin
      write (56,*) 'rtscheck ',rtscheck
      write (56,*) 'numstruct ',noptg
      write (56,*) rtscheck
      do iopt = 1 , noptg
         if(rtscheck.ne.0)then
            if (vref(iopt).lt.voptmin) then
               if(iabs.eq.1)then
                  if(abs(rtscheck-rtsopt(iopt)).lt.(0.4))then
                     if(abs(bbdischeck-bbopt(iopt)).lt.(0.4))then
                        ioptmin=iopt
                        voptmin= vref(ioptmin)
                     endif
                  endif
               else if(iadd.eq.1)then
                  if(abs(rtscheck-rtsopt(iopt)).lt.(0.15))then
                     ioptmin=iopt
                     voptmin= vref(ioptmin)
                  endif
               else if(abs(rtscheck-rtsopt(iopt)).lt.(0.4))then
                  ioptmin=iopt
                  voptmin= vref(ioptmin)
               endif
            endif
         else
           if (vref(iopt).lt.voptmin) then
              ioptmin=iopt
              voptmin= vref(ioptmin)
           endif
         endif
      enddo
      write (56,*) 'energy '
      write (56,*) voptmin
      write (56,*) 'ref energy '
      write (56,*) vref(ioptmin)
      write (56,*) 'opt geom',ioptmin
      if (ifrozrts.ne.1)then
         do iint = 1 , nint
            write (56,*) xinto(iint,ioptmin)
         enddo
         do itau=1,ntau
            write (56,*) tauo(itau,ioptmin)
         enddo
      else
         do iint = 1 , ncoord
            write (56,*) xinto(iint,ioptmin)
         enddo
      endif
      close (unit=55,status='keep')
      close (unit=56,status='keep')

cc added to have a file with highest level TS z-matrix parameters
cc this is updated by the higher level TS search and random dihedral scans

      open (unit=37,file='./output/ts_opt.out',status='unknown')

      write (37,*) 'TS_z-matrix'
      do iatom = 1 , natomt
         write (37,'(A60)') atomlabel(iatom)
      enddo
      if (ifrozrts.ne.1)then
         do iint = 1 , nint
            write (37,*) intcoor(iint),xinto(iint,ioptmin)
         enddo
         do itau=1,ntau
            write (37,*) bislab(itau),tauo(itau,ioptmin)
         enddo
      else
         do iint = 1 , ncoord
            write (37,*) intcoor(iint),xinto(iint,ioptmin)
         enddo
      endif
      close (unit=37,status='keep')

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine level1(ispecies)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      character*70 comline1,comline2
      character*70 comline3,comline4
      character*85 geomline
      character*30 intcoor(3*natommx)
      character*60 atomlabel(natommx)
      character*100 command1
      character*100 commandcopy
      character*2 atomname
      character*2 cjunk
      character*8 nameout
      character*20 bislab(ntaumx)
      character*160 word_l1
      character*30 gmem
      character*30 distname

      dimension freq(nmdmx),freqp(nmdmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim) 
      dimension gelec(nelecmx),eelec(nelecmx)
      dimension taumn(ntaumx),taumx(ntaumx)

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'


      open (unit=66,file='./output/level1.out',status='unknown')
      open (unit=21,file='./data/theory.dat',status='unknown')

      call LineRead(0)
      commandcopy='cp -f ./data/level1_molpro.dat 
     $                    level1_molpro.dat'            

cs added option for species dependent levels of theory

      if (ispecies.eq.1) then
         open (unit=15,file='./data/reac1.dat',status='old')
         open (unit=16,file='./me_files/reac1_ge.me',status='unknown')
         open (unit=19,file='./me_files/reac1_fr.me',status='unknown')
         open (unit=20,file='./me_files/reac1_zpe.me',status='unknown')
         open (unit=17,file='./output/reac1_opt.out',status='unknown')
c         open (unit=116,file='./me_files/rpst_ge1.me',status='unknown')
c         open (unit=119,file='./me_files/rpst_fr1.me',status='unknown')
c         open (unit=120,file='./me_files/rpst_zpe1.me',status='unknown')
         nameout='reac1_l1'
         inp_type=1
         rewind(21)
         WORD_L1='LEVEL1_1'
         do while (WORD.NE.'LEVEL1_1')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L1 = 'LEVEL1'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level1_reac1_molpro.dat 
     $                    level1_molpro.dat'            
         endif
      endif
      if (ispecies.eq.2) then
         open (unit=15,file='./data/reac2.dat',status='old')
         open (unit=16,file='./me_files/reac2_ge.me',status='unknown')
         open (unit=19,file='./me_files/reac2_fr.me',status='unknown')
         open (unit=20,file='./me_files/reac2_zpe.me',status='unknown')
         open (unit=17,file='./output/reac2_opt.out',status='unknown')
c         open (unit=116,file='./me_files/rpst_ge2.me',status='unknown')
c         open (unit=119,file='./me_files/rpst_fr2.me',status='unknown')
c         open (unit=120,file='./me_files/rpst_zpe2.me',status='unknown')
         nameout='reac2_l1'
         inp_type=1
         rewind(21)
         WORD_L1='LEVEL1_2'
         do while (WORD.NE.'LEVEL1_2')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L1 = 'LEVEL1'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level1_reac2_molpro.dat 
     $                    level1_molpro.dat'            
         endif
      endif
      if (ispecies.eq.3) then
         open (unit=15,file='./data/prod1.dat',status='old')
         open (unit=16,file='./me_files/prod1_ge.me',status='unknown')
         open (unit=19,file='./me_files/prod1_fr.me',status='unknown')
         open (unit=20,file='./me_files/prod1_zpe.me',status='unknown')
         open (unit=17,file='./output/prod1_opt.out',status='unknown')
c         open (unit=116,file='./me_files/ppst_ge1.me',status='unknown')
c         open (unit=119,file='./me_files/ppst_fr1.me',status='unknown')
c         open (unit=120,file='./me_files/ppst_zpe1.me',status='unknown')
         nameout='prod1_l1'
         inp_type=1
         rewind(21)
         WORD_L1='LEVEL1_3'
         do while (WORD.NE.'LEVEL1_3')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L1 = 'LEVEL1'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level1_prod1_molpro.dat 
     $                    level1_molpro.dat'            
         endif
      endif
      if (ispecies.eq.4) then
         open (unit=15,file='./data/prod2.dat',status='old')
         open (unit=16,file='./me_files/prod2_ge.me',status='unknown')
         open (unit=19,file='./me_files/prod2_fr.me',status='unknown')
         open (unit=20,file='./me_files/prod2_zpe.me',status='unknown')
         open (unit=17,file='./output/prod2_opt.out',status='unknown')
c         open (unit=116,file='./me_files/ppst_ge2.me',status='unknown')
c         open (unit=119,file='./me_files/ppst_fr2.me',status='unknown')
c         open (unit=120,file='./me_files/ppst_zpe2.me',status='unknown')
         nameout='prod2_l1'
         inp_type=1
         rewind(21)
         WORD_L1='LEVEL1_4'
         do while (WORD.NE.'LEVEL1_4')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L1 = 'LEVEL1'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level1_prod2_molpro.dat 
     $                    level1_molpro.dat'            
         endif
      endif
      if (ispecies.eq.5) then
         open (unit=15,file='./data/wellr.dat',status='old')
         open (unit=16,file='./me_files/wellr_ge.me',status='unknown')
         open (unit=19,file='./me_files/wellr_fr.me',status='unknown')
         open (unit=20,file='./me_files/wellr_zpe.me',status='unknown')
         open (unit=17,file='./output/wellr_opt.out',status='unknown')
         nameout='wellr_l1'
         inp_type=1
         rewind(21)
         WORD_L1='LEVEL1_5'
         do while (WORD.NE.'LEVEL1_5')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L1 = 'LEVEL1'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level1_wellr_molpro.dat 
     $                    level1_molpro.dat'            
         endif
      endif
      if (ispecies.eq.51) then
         open (unit=15,file='./data/wellr.dat',status='old')
         open (unit=16,file='./me_files/wellr_ge.me',status='unknown')
         open (unit=19,file='./me_files/wellr_fr.me',status='unknown')
         open (unit=20,file='./me_files/wellr_zpe.me',status='unknown')
         if(igeom_wellr.eq.2) then
            open (unit=17,file='./output/ts_opt.out',status='unknown')
         else
            open(unit=17,file='./output/wellr_opt.out',status='unknown')
         endif
         nameout='wellr_l1'
         inp_type=2
         rewind(21)
         WORD_L1='LEVEL1_51'
         do while (WORD.NE.'LEVEL1_51')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L1 = 'LEVEL1'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level1_wellr_molpro.dat 
     $                    level1_molpro.dat'            
         endif
      endif
      if (ispecies.eq.6) then
         open (unit=15,file='./data/wellp.dat',status='old')
         open (unit=16,file='./me_files/wellp_ge.me',status='unknown')
         open (unit=19,file='./me_files/wellp_fr.me',status='unknown')
         open (unit=20,file='./me_files/wellp_zpe.me',status='unknown')
         open (unit=17,file='./output/wellp_opt.out',status='unknown')
         nameout='wellp_l1'
         inp_type=1
         rewind(21)
         WORD_L1='LEVEL1_6'
         do while (WORD.NE.'LEVEL1_6')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L1 = 'LEVEL1'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level1_wellp_molpro.dat 
     $                    level1_molpro.dat'            
         endif
      endif
      if (ispecies.eq.61) then
         open (unit=15,file='./data/wellp.dat',status='old')
         open (unit=16,file='./me_files/wellp_ge.me',status='unknown')
         open (unit=19,file='./me_files/wellp_fr.me',status='unknown')
         open (unit=20,file='./me_files/wellp_zpe.me',status='unknown')
         if(igeom_wellp.eq.2)then
            open (unit=17,file='./output/ts_opt.out',status='unknown')
         else
            open(unit=17,file='./output/wellp_opt.out',status='unknown')
         endif
         nameout='wellp_l1'
         inp_type=2
         rewind(21)
         WORD_L1='LEVEL1_61'
         do while (WORD.NE.'LEVEL1_61')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L1 = 'LEVEL1'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level1_wellp_molpro.dat 
     $                    level1_molpro.dat'            
         endif
      endif
      if (ispecies.eq.11) then
         open (unit=15,file='./data/reacs.dat',status='old')
         open (unit=16,file='./me_files/reacs_ge.me',
     +         status='unknown')
         open (unit=19,file='./me_files/reacs_fr.me',
     +         status='unknown')
         open (unit=20,file='./me_files/reacs_zpe.me',status='unknown')
         open (unit=17,file='./output/reacs_opt.out',status='unknown')
         nameout='reacs_l1'
         inp_type=1
         rewind(21)
         WORD_L1='LEVEL1_11'
         do while (WORD.NE.'LEVEL1_11')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L1 = 'LEVEL1'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level1_reacs_molpro.dat 
     $                    level1_molpro.dat'            
         endif
      endif
      if (ispecies.eq.12) then
         open (unit=15,file='./data/prods.dat',status='old')
         open (unit=16,file='./me_files/prods_ge.me',
     +         status='unknown')
         open (unit=19,file='./me_files/prods_fr.me',
     +         status='unknown')
         open (unit=20,file='./me_files/prods_zpe.me',status='unknown')
         open (unit=17,file='./output/prods_opt.out',status='unknown')
         nameout='prods_l1'
         inp_type=1
         rewind(21)
         WORD_L1='LEVEL1_12'
         do while (WORD.NE.'LEVEL1_12')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L1 = 'LEVEL1'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level1_prods_molpro.dat 
     $                    level1_molpro.dat'            
         endif
      endif
      if (ispecies.eq.0) then
         do while (WORD.NE.'LEVEL1_TS')
            call LineRead(21)
            if (WORD.EQ.'END') then
               write(66,*)'The TS level1 of theory must be specified'
               write(66,*)'in file theory.dat'
               stop
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/level1_ts_molpro.dat 
     $                    level1_molpro.dat'            
         endif
         WORD=''
         rewind(21)

         open (unit=15,file='./data/ts.dat',status='old')
         do while (WORD.NE.'TS_GEOM')
            call LineRead(15)
            if (WORD.EQ.'END') then
               rewind(15)
               go to  901
            endif
         enddo
         if(word2.ne.'HR')then
            read(15,*)its
            open (unit=99,status='unknown')
            rewind (99)
            write (99,1114) its
            rewind (99)
            read (99,'(A100)') command1
            close (99)
            call commrun(command1)         
 1114       format ('cp -f  output/ts_opt_',I0.2,
     $           '.out output/ts_opt.out')
         else
            read(15,*)ihind,its
            write(*,*)'taking geometry from hind rot ',ihind
            write(*,*)'point ',its
            open (unit=99,status='unknown')
            write (99,1115) ihind,its
            rewind (99)
            read (99,'(A100)') command1
            close (99)
            call commrun(command1)         
 1115       format ('cp -f  hr_geoms/geom_isp00_hr',I0.2,'_hpt',I0.2,
     $           '.out output/ts_opt.out')
c            stop
         endif
         rewind(15)
 901     continue

         open (unit=16,file='./me_files/ts_ge.me',status='unknown')
         open (unit=19,file='./me_files/ts_fr.me',status='unknown')
         open (unit=20,file='./me_files/ts_zpe.me',status='unknown')
         open (unit=17,file='./output/ts_opt.out',status='unknown')
         nameout='tsgta_l1'
         inp_type=2
         rewind(21)
         WORD_L1='LEVEL1_TS'
         do while (WORD.NE.'LEVEL1_TS')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_L1 = 'LEVEL1'
               go to  900
            endif
         enddo
      endif
 900  continue
      rewind(21)

cc initialize paramters

      ires=1

cc first read blocks that are common for input type 1 and 2

      if (idebug.ge.2) write (66,*) ' starting zmat input'
      call LineRead (0)

      do while (WORD.NE.'SYMMETRYFACTOR')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (66,*) 'symmetry factor must be defined'
            stop
         endif
      enddo
      read (15,*) symf
      rewind(15)

      do while (WORD.NE.'NHIND')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (66,*) 'hind rotors must be defined'
            stop
         endif
      enddo
      read (15,*) nhind
      rewind(15)

      do while (WORD.NE.'NELEC')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (66,*) 'electronic states must be defined'
            stop
         endif
      enddo
      read (15,*) nelec
c      write (66,*) 'test4',nelec
      do ielec = 1 , nelec
         read (15,*) eelec(ielec),gelec(ielec)
      enddo
      rewind(15)

cc read scaling factor for frequencies.

      sclfr=1.0
      iscl=0
      rewind(21)
      call LineRead (0)
      do while (WORD.NE.'END')
         call LineRead (21)
         if (WORD.EQ.'SCALEFREQ'.and.iscl.eq.0) then
            read(21,*)sclfr
         endif
         if (WORD.EQ.'SCALEFREQ_TS'.and.iscl.eq.0.and.
     +        ispecies.eq.0) then
            read(21,*)sclfr
            iscl=1
         endif
         if (WORD.EQ.'SCALEFREQ_REAC1'.and.iscl.eq.0.and.
     +        ispecies.eq.1) then
            read(21,*)sclfr
            iscl=1
         endif
         if (WORD.EQ.'SCALEFREQ_REAC2'.and.iscl.eq.0.and.
     +        ispecies.eq.2) then
            read(21,*)sclfr
            iscl=1
         endif
         if (WORD.EQ.'SCALEFREQ_PROD1'.and.iscl.eq.0.and.
     +        ispecies.eq.3) then
            read(21,*)sclfr
            iscl=1
         endif
         if (WORD.EQ.'SCALEFREQ_PROD2'.and.iscl.eq.0.and.
     +        ispecies.eq.4) then
            read(21,*)sclfr
            iscl=1
         endif
         if (WORD.EQ.'SCALEFREQ_WELLR'.and.iscl.eq.0.and.
     +        ispecies.eq.5) then
            read(21,*)sclfr
            iscl=1
         endif
         if (WORD.EQ.'SCALEFREQ_WELLP'.and.iscl.eq.0.and.
     +        ispecies.eq.6) then
            read(21,*)sclfr
            iscl=1
         endif
         if (WORD.EQ.'SCALEFREQ_WELLR'.and.iscl.eq.0.and.
     +        ispecies.eq.51) then
            read(21,*)sclfr
            iscl=1
         endif
         if (WORD.EQ.'SCALEFREQ_WELLP'.and.iscl.eq.0.and.
     +        ispecies.eq.61) then
            read(21,*)sclfr
            iscl=1
         endif
      enddo
      rewind(21)
c      write(*,*)'ok cc1'
c      stop

cc read specific input for input type 1

      write (6,*) 'inp_type test',inp_type,ispecies

      if (inp_type.eq.1) then

c     natomt is to account for dummy atoms
         do while (WORD.NE.'NATOM')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'natom must be defined'
               stop
            endif
         enddo
         read (15,*) natom,natomt,ilin
         if (natomt.gt.natommx) then
            write (66,*) 'natomt too large',natomt,natommx
            stop
         endif
         rewind(15)

         write (66,*) 'test0',natomt,atomlabel(1)
         do while (WORD.NE.'CHARGE')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'charge and spin must be defined'
               stop
            endif
         enddo
         read (15,*) icharge,ispin
         do iatom = 1 , natomt
            read (15,'(A60)') atomlabel(iatom)
         enddo
         rewind(15)

         do while (WORD.NE.WORD_L1)
            call LineRead (21)
            if (WORD.EQ.'END') then
               write (66,*) 'level1 of theory must be defined'
               write (66,*) 'in file theory.dat'
               stop
            endif
         enddo
         if(word2.eq.'G09') then
            ilev1code=1
            read (21,'(A70)') comline1
            read (21,'(A70)') comline2
            call comline34_g09(ispecies,comline1,comline2,
     $                         comline3,comline4)
            if (idebug.ge.2) write (66,*) 'coml 12 ',comline1,comline2
            if (idebug.ge.2) write (66,*) 'coml 34 ',comline3,comline4
         else if(word2.eq.'MOLPRO') then
            ilev1code=2
            call commrun(commandcopy)
         else
            write (66,*) 'level1 of theory must be either'
            write (66,*) 'g09 or molpro in theory.dat'
            stop
         endif
         if(word3.eq.'RESTART') then
            open (unit=99,status='unknown')
            write (99,*) word4
            rewind (99)
            read (99,*) ires
            close (unit=99,status='keep')
         endif
         close(21)


         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'sampling coordinates must be defined'
               stop
            endif
         enddo
         read (15,*) ntau
         rewind(15)
         call  LineRead (0)

         if (natom.ne.1) then
            do while (WORD.NE.'INTCOOR')
               call LineRead (15)
               if (WORD.EQ.'END') then
                  write (66,*) 'internal coordinates must be defined'
                  stop
               endif
            enddo
            ncoord = 3*natom-6-ntau
            if (natom.eq.1) ncoord = 0
            if (natom.eq.2) ncoord = 1
            do icoord = 1 , ncoord
               call LineRead (15)
               intcoor(icoord) = word
            enddo
            rewind(15)
         endif

         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'sampling coordinates must be defined'
               stop
            endif
         enddo
         read (15,*) ntau
         if (ntau.gt.ntaumx) then
            write (66,*) 'ntau too large',ntau,ntaumx
            stop
         endif
         if (ntau.ne.0) then
            read (15,*)
            do itau = 1 , ntau
               read (15,*) bislab(itau),taumn(itau),taumx(itau)
               write (66,*) bislab(itau),taumn(itau),taumx(itau)
               open (unit=99,status='unknown')
               rewind (99)
               write (99,*)bislab(itau)
               rewind (99)
               call LineRead (99)
               icoord=ncoord+itau
               intcoor(icoord)=WORD
               close (unit=99,status='keep')
            enddo
         endif
         rewind(15)

cc now read optimized geometry parameters
         ncoord = 3*natom-6
         if(natom.eq.2) ncoord = 1
         if(natom.eq.1) ncoord = 0
         read (17,*)
         do icoord = 1 , ncoord
            call LineRead (17)
c     xinti(icoord) = word
c     intcoori(icoord) = word
            OPEN (unit=99,status='unknown')
            REWIND (99)
            WRITE (99,1000) WORD
 1000       FORMAT (A30)
            REWIND (99)
            READ (99,*) xint(icoord)
            close (unit=99,status='keep')
         enddo
         close (unit=17,status='keep')
         close (unit=15,status='keep')

cc get data for wellp, necessary for species 5 or 6 optimization

         if(ispecies.eq.5.or.ispecies.eq.6) then
            open (unit=99,file='data/ts.dat',status='old')
            do while (WORD.NE.'ISITE')
               call LineRead (99)
               if (WORD.EQ.'END') then
                  write (66,*) 'reaction site must be defined'
                  write (66,*) 'in file ts.dat'
                  stop
               endif
            enddo
            read (99,*) isite,jsite,ksite

            if(ibeta.eq.1)then
               rewind(99)
               do while (WORD.NE.'RMAX1')
                  call LineRead (99)
                  if (WORD.EQ.'END') then
                     write (66,*) 'grid coords must be defined'
                     stop
                  endif
               enddo
               read (99,*) rmax1,rstp1,rmax2,rstp2,ireact
            endif

            if(iiso.eq.1)then
               rewind(99)
               do while (WORD.NE.'RMIN1')
                  call LineRead (99)
                  if (WORD.EQ.'END') then
                     write (66,*) 'grid coords must be defined'
                     stop
                  endif
               enddo
               read (99,*) rmin1,rstp1,rmin2,rstp2,ireact
            endif
            close(99)
         endif

      else if (inp_type.eq.2) then

cc here we assume that the TS is not linear

         ilin=0
         ilin1=0
         ilin2=0

cc now read input of type 2

         do while (WORD.NE.'CHARGE')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'charge and spin must be defined'
               stop
            endif
         enddo
         read (15,*) icharge,ispin
         rewind(15)

c      write(*,*)'ok cc1'
c      stop
         if(ispecies.eq.0) then
            do while (WORD.NE.'LEVEL1_TS')
               call LineRead (21)
               if (WORD.EQ.'END') then
                  write (66,*) 'level1 of theory must be defined for TS'
                  write (66,*) 'in file theory.dat'
                  stop
               endif
            enddo
         else if(ispecies.eq.51.or.ispecies.eq.61) then
            do while (WORD.NE.WORD_L1)
               call LineRead (21)
               if (WORD.EQ.'END') then
                  write (66,*) 'level1 of theory must be defined'
                  write (66,*) 'in file theory.dat'
                  stop
               endif
            enddo
         endif
         if(word2.eq.'G09') then
            ilev1code=1
            read (21,'(A70)') comline1
            read (21,'(A70)') comline2
            call comline34_g09(ispecies,comline1,comline2,
     $                         comline3,comline4)
            if (idebug.ge.2) write (66,*)' coml 12',comline1,comline2
            if (idebug.ge.2) write (66,*)' coml 34',comline3,comline4
         else if(word2.eq.'MOLPRO') then
            ilev1code=2
            call commrun(commandcopy)
         else
            write (66,*) 'level1 of theory must be either'
            write (66,*) 'g09 or molpro in theory.dat'
            stop
         endif
         if(word3.eq.'RESTART') then
            open (unit=99,status='unknown')
            write (99,*) word4
            rewind (99)
            read (99,*) ires
            close (unit=99,status='keep')
         endif
         close(21)

         open (unit=25,file='./data/reac1.dat',status='old')

         do while (WORD.NE.'NTAU')
            call LineRead (25)
            if (WORD.EQ.'END') then
               write (66,*) 'sampling coords of reac1 must be defined'
               stop
            endif
         enddo
         read (25,*) ntau1
         rewind 25

c     natomt is to account for dummy atoms
         do while (WORD.NE.'NATOM')
            call LineRead (25)
            if (WORD.EQ.'END') then
               write (66,*) 'natom in reac1 must be defined'
               stop
            endif
         enddo
         read (25,*) natom1,natomt1,ilin1
         close (25)

cc get data from react2 file

         if(iabs.eq.1.or.iadd.eq.1)then
            open (unit=25,file='./data/reac2.dat',status='old')
            
            do while (WORD.NE.'NTAU')
               call LineRead (25)
               if (WORD.EQ.'END') then
                  write (66,*) 'samp coords of reac2 must be defined'
                  stop
               endif
            enddo
            read (25,*) ntau2
            rewind 25

c     natomt is to account for dummy atoms
            do while (WORD.NE.'NATOM')
               call LineRead (25)
               if (WORD.EQ.'END') then
                  write (66,*) 'natom must be defined'
                  stop
               endif
            enddo
            read (25,*) natom2,natomt2,ilin2
            close (25)
         endif
cc now we can determine the total number of atoms for the TS/wellr/wellp
         if(iadd.eq.1.or.iabs.eq.1)then
            natom = natom1+natom2
         else if (iiso.eq.1.or.ibeta.eq.1)then
            natom = natom1
         endif
cc modified         if (iadd.eq.1) natomt = natomt1+natomt2
         if (iadd.eq.1) natomt = natomt1+natomt2
         if (iabs.eq.1) natomt = natomt1+natomt2+1
         if (iiso.eq.1) natomt = natomt1
         if (ibeta.eq.1) natomt = natomt1

cc if linear molecule reacting with atom in abs assume TS is linear
         if(iabs.eq.1.and.natom2.eq.1.and.ispecies.eq.0.
     $      and.ilin1.eq.1) then
            ilin=1
         endif

cc if linear molecule reacting with linear radical in abs assume TS is linear
         if(iabs.eq.1.and.ilin2.eq.1.and.ispecies.eq.0.and
     $      .ilin1.eq.1)then
            ilin=1
         endif
c        natomt = natomt1+natomt2+1
         
         close (unit=15,status='keep')

c natomt is to account for dummy atoms
         if (natomt.gt.natommx) then
            write (66,*) 'natomt too large',natomt,natommx
            stop
         endif

c gaussian com file data
         read (17,*)
         if (idebug.ge.2) write (6,*) ' starting gaussian input'
         do iatom = 1 , natomt
            read (17,'(A60)') atomlabel(iatom)
         enddo

cc determine distance to check for the findgeom option for wellp and wellr

         if(ispecies.eq.51.or.ispecies.eq.61)then
            open (unit=99,file='data/ts.dat',status='old')
            do while (WORD.NE.'ISITE')
               call LineRead (99)
               if (WORD.EQ.'END') then
                  write (66,*) 'reaction site must be defined'
                  write (66,*) 'in file ts.dat'
                  stop
               endif
            enddo
            read (99,*) isite,jsite,ksite

c determine number of reacting atom for beta-scission 
c and isomerization reactions
            if(ibeta.eq.1)then
               rewind(99)
               do while (WORD.NE.'RMAX1')
                  call LineRead (99)
                  if (WORD.EQ.'END') then
                     write (66,*) 'grid coords must be defined'
                     stop
                  endif
               enddo
               read (99,*) rmax1,rstp1,rmax2,rstp2,ireact
            endif

            if(iiso.eq.1)then
               rewind(99)
               do while (WORD.NE.'RMIN1')
                  call LineRead (99)
                  if (WORD.EQ.'END') then
                     write (66,*) 'grid coords must be defined'
                     stop
                  endif
               enddo
               read (99,*) rmin1,rstp1,rmin2,rstp2,ireact
            endif
            close(99)
            

            write(*,*)'isite is',isite
            open (unit=99,file='geoms/tsgta_l1.xyz',status='old')
            read(99,*)
            read(99,*)
            do j=1,natom
               read(99,*)cjunk,coox,cooy,cooz
               if (j.eq.isite)then
                  xsite=coox
                  ysite=cooy
                  zsite=cooz
                  write(*,*)'xsite is',xsite
               endif
               if (j.eq.jsite)then
                  xji=coox
                  yji=cooy
                  zji=cooz
               endif
            enddo
            close(99)
            distcheck=sqrt((xsite-xji)**2+(ysite-yji)**2
     $                +(zsite-zji)**2)

c for beta-decomp determine label of breaking bond      
            distname='RTS'
            write(66,*)'ireact is ', ireact
            if(ibeta.eq.1)then
               open(unit=99,status='unknown')
               do iatom = 1 , natomt
                  rewind (99)
                  write(99,*)atomlabel(iatom)
                  rewind (99)
                  call LineRead(99)
                  if(iatom.eq.ireact)then
                     distname=word3
                  endif
               enddo
               write(66,*)'distname is ', distname
            endif
         endif

cc read coordinate names

         ncoord = 3*natom-6

         do iint = 1 , ncoord
            read (17,*) intcoor(iint),xint(iint)
            if(ispecies.eq.51.and.igeom_wellr.eq.2)then
               if(intcoor(iint).eq.distname) then
                  if(distcheck.gt.1.4.and.iabs.eq.1)then
                     fgdisp51=0.5
                  else if(distcheck.gt.1.4.and.iadd.eq.1)then
                     fgdisp51=0.5
                  else
                     fgdisp51=0.2
                  endif
                  if(ibeta.eq.1)fgdisp51=-fgdisp51
                  xint(iint)=xint(iint)+fgdisp51
               endif
            endif
            if(ispecies.eq.61.and.igeom_wellp.eq.2)then
               if(intcoor(iint).eq.distname) then
                  if(xint(iint).gt.1.4) then
                     fgdisp61=0.5
                  else
                     fgdisp61=0.2
                  endif
                  if(ibeta.eq.1)fgdisp61=-fgdisp61
                  xint(iint)=xint(iint)-fgdisp61
               endif
            endif
         enddo
         close (unit=17,status='keep')
      endif

      ntau = 0
      ircons=0
      ismp=0
      ifreq=1
      ilev=1
      if(ifrozrts.eq.1.and.ispecies.eq.0) then
         ircons=1
         xint(ncoord)=frozcoo
      endif
      ixyz=0
      ired=0
      if(iabs.eq.1) ireact=natom1+1
      if(iadd.eq.1) ireact=natom1+1
      if(ispecies.eq.51.or.ispecies.eq.61)then
         if(ilin1.eq.0)then
            ired=1
         else if(ilin1.eq.1)then
            ired=2
         endif
      endif
      if(ispecies.eq.5.or.ispecies.eq.6)then
         if(ilin1.eq.0)then
            ired=1
         else if(ilin1.eq.1)then
            ired=2
         endif
      endif
      if(ispecies.eq.6.and.iadd.eq.1)ired=0
      if(ispecies.eq.61.and.iadd.eq.1)ired=0
      
      write(*,*)'ilin is',ilin

      if(ilev1code.eq.1) then
         call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $ coord,vtot_0,vtotr,freq,ifreq,ilin,ismp,comline1,
     $ comline2,icharge,ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ired)

      else if (ilev1code.eq.2) then
         numproc=numprochl
         
         call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $ coord,vtotr,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ilev
     $ ,ispecies)
         numproc=numprocll
      endif

      if (inp_type.eq.1) then

cc write me input
         if (ispecies.eq.1) then
            if(iadd.eq.1.or.iabs.eq.1)then
               write (16,*) 'Bimolecular REACS'
               write (16,*) 'Fragment REACT1'
            else if (iiso.eq.1.or.ibeta.eq.1)then
               write (16,*) 'Well REACS'
               write (16,*) 'Species'
            else
               write (16,*) 'Well REACS'
               write (16,*) 'Species'
            endif
c            write (116,*) 'Barrier	B1  REACT WR'
c            write (116,*) '  RRHO '
c            write (116,*) '    Stoichiometry  ',stoich
c            write (116,*) '    Core PhaseSpaceTheory'
c            write (116,*) '      FragmentGeometry[angstrom]    ',natom
         endif
         if (ispecies.eq.2) then
            write (16,*) 'Fragment REACT2'
c            write (116,*) '      FragmentGeometry[angstrom]    ',natom
         endif
cc for calculation of backward rate         
         if (ispecies.eq.3.and.ipr1.eq.0) then
            write (16,*) 'Bimolecular PRODS'
            write (16,*) 'Fragment PROD1'
c            write (116,*) 'Barrier	B3  WP PRODS'
c            write (116,*) '  RRHO '
c            write (116,*) '    Stoichiometry  ',stoich
c            write (116,*) '    Core PhaseSpaceTheory'
c            write (116,*) '      FragmentGeometry[angstrom]    ',natom
         else if (ispecies.eq.3.and.ipr1.eq.1) then
            write (16,*) 'Well PRODS'
            write (16,*) 'Species'
         endif
         if (ispecies.eq.4) then
            write (16,*) 'Fragment PROD2'
c            write (116,*) '      FragmentGeometry[angstrom]    ',natom
         endif
         if (ispecies.eq.5) then
            write (16,*) 'Well WR'
            write (16,*) 'Species'
         endif
         if (ispecies.eq.6) then
            write (16,*) 'Well WP'
            write (16,*) 'Species'
         endif
         if (natom.ne.1) then
            write (16,*) 'RRHO     '
            write (16,*) 'Geometry[angstrom]    ',natom
c            command1='newzmat       
c     +           -ichk -oxyz tmp.chk geom.xyz           '
c            call commrun(command1)
            open (unit=98,file='headgeom.tmp',status='unknown')
            write(98,*)natom
            write(98,*)vtotr
            close(98)
            command1='cat headgeom.tmp geom.xyz > geom1.xyz '
            call commrun(command1)

            open (unit=99,status='unknown')
            rewind (99)
            write (99,1112) nameout
            rewind (99)
            read (99,1020) command1
            close(99)
            call commrun(command1)         
 1112       format ('cp -f geom1.xyz geoms/'A8,'.xyz')

            open (unit=99,status='unknown')
            rewind (99)
            write (99,1122) nameout
            rewind (99)
            read (99,1020) command1
            close(99)
            call commrun(command1)         
 1122       format ('cp -f geom.log geoms/'A8,'.log')
            if (ilev1code.eq.2) then
               open (unit=99,status='unknown')
               rewind (99)
               write (99,1123) nameout
               rewind (99)
               read (99,1020) command1
               close(99)
               call commrun(command1)         
 1123          format ('cp -f molpro.molden geoms/'A8,'.molden')
            endif

 1020       format (A70)

c         if (natom.ne.1) then
            open (unit=100,FILE='geom.xyz',status='unknown')
            do j=1,natom
               read (100,1001) geomline
               write (16,*) geomline
c               if (ispecies.le.4) write (116,*) geomline
            enddo
 1001       FORMAT (A85)
            close (unit=100,status='keep')
c     command1='cat reac_1.me geom.xyz > temp.me'
c     call commrun(command1)
c     command1='cp -f temp.me reac_1.me'
c     call commrun(command1)
            write (16,*) '	Core 	RigidRotor'
            write (16,*) '		SymmetryFactor ',symf
            write (16,*) '	End'
         endif
c        if (ispecies.eq.1) symfr1 = symf
c        if (ispecies.eq.2) symfr2 = symf
c        symfr = symfr1*symfr2
c        if (ispecies.eq.3) symfp1 = symf
c        if (ispecies.eq.4) symfp2 = symf
c        symfp = symfp1*symfp2
c these will create problems if reactant 1 and 2 are not run together or if
c products 1 and 2 are not run together
c an alternative would be to read the symmetry factors from the appropriate me files
         if (natom.eq.1) then

            open (unit=98,file='headgeom.tmp',status='unknown')
            write(98,*)natom
            write(98,*)vtotr
            close(98)
            command1='cat headgeom.tmp geom.xyz > geom1.xyz '
            call commrun(command1)

            open (unit=99,status='unknown')
            rewind (99)
            write (99,3112) nameout
            rewind (99)
            read (99,1020) command1
            close(99)
            call commrun(command1)         
 3112       format ('cp -f geom1.xyz geoms/'A8,'.xyz')

            open (unit=99,status='unknown')
            rewind (99)
            write (99,3122) nameout
            rewind (99)
            read (99,1020) command1
            close(99)
            call commrun(command1)         
 3122       format ('cp -f geom.log geoms/'A8,'.log')
            if (ilev1code.eq.2) then
               open (unit=99,status='unknown')
               rewind (99)
               write (99,3123) nameout
               rewind (99)
               read (99,1020) command1
               close(99)
               call commrun(command1)         
 3123          format ('cp -f molpro.molden geoms/'A8,'.molden')
            endif


            write (16,*) ' Atom'
c            command1='newzmat       
c     +           -ichk -oxyz tmp.chk geom.xyz           '
c            call commrun(command1)
            open (unit=100,FILE='geom.xyz',status='unknown')
            read (100,*) atomname
            write (16,*) ' Name ',atomname
c           write (116,*) atomname, '  0.  0.  0. '
            close (100)
         endif
         if (ispecies.eq.2) then 
            open (unit=126,file='./me_files/reac1_ge.me',status=
     $       'unknown')
            do while (WORD.NE.'SYMMETRYFACTOR')
               call LineRead (126)
            enddo
            open (unit=99)
            rewind (99)
            write (99,*) word2
            rewind (99)
            read (99,*) symfr1
            close (unit=99)
            symfr = symfr1*symf
c            write (116,*) '		SymmetryFactor ',symfr
            close (unit=126,status='keep')
         endif
         if (ispecies.eq.4) then 
            open (unit=126,file='./me_files/prod1_ge.me',status=
     $       'unknown')
            do while (WORD.NE.'SYMMETRYFACTOR')
               call LineRead (126)
            enddo
            open (unit=99)
            rewind (99)
            write (99,*) word2
            rewind (99)
            read (99,*) symfp1
            close (unit=99)
            symfp = symfp1*symf
c            write (116,*) '		SymmetryFactor ',symfp
            close (unit=126,status='keep')
         endif
c         if ((ispecies.eq.2).or.(ispecies.eq.4)) then
c            write (116,*) '          PotentialPrefactor[au] 	10.' 
c            write (116,*) '          PotentialPowerExponent 	6' 
c            write (116,*) '	    End'
c         endif
         nfreq=3*natom-6
         if (natom.eq.2) nfreq=1
         if (natom.eq.1) nfreq=0
         if (idebug.ge.2) write (6,*) 'freq test',natom,ispecies,
     $    nfreq,freq(1)
         if (natom.ne.1) then
            if(ilin.eq.1)nfreq=nfreq+1
            write (19,*) '    Frequencies[1/cm] ',nfreq
c            if ((ispecies.eq.1).or.(ispecies.eq.3)) 
c     $        write (119,*) '    Frequencies[1/cm] $nfreq'
c            if (ispecies.eq.1) nfreqr = nfreq
c            if (ispecies.eq.2) nfreqr = nfreqr+nfreq
c           write (6,*) 'nfreqr test',nfreqr,ispecies,nfreq
c            if (ispecies.eq.3) nfreqp = nfreq
c            if (ispecies.eq.4) nfreqp = nfreqp+nfreq
            write (19,8010) (freq(j),j=1,nfreq)
c            write (119,8010) (freq(j),j=1,nfreq)
8010        format (1x,10G12.5)
            zpe=0
            do j=1,nfreq
               zpe=zpe+freq(j)
            enddo
cc       divide by two and convert to kcal/mol
            zpe=zpe/(cautoicm*2.d0)
c           zpe=zpe*29979200000*6.626E-34/2.*6.022E23/4.184/1000/627.503
            write (20,*) zpe
c            write (120,*) zpe
         endif
         if (natom.eq.1)  then 
            write (20,*) ' 0.'
c            write (120,*) ' 0.'
         endif
         if (natom.ne.1) then
            if (ispecies.le.4) then
               if(ipr1.eq.0) then
                  write (19,*) ' ZeroEnergy[kcal/mol]             0.'
               else if(ipr1.eq.1.and.ispecies.eq.3) then
                  write (19,*) ' ZeroEnergy[kcal/mol]      $proden'
               else if(ipr1.eq.1) then
                  write (19,*) ' ZeroEnergy[kcal/mol]             0.'
               endif
            endif
            if (ispecies.eq.5) 
     $       write (19,*) ' ZeroEnergy[kcal/mol]            $wellren'
            if (ispecies.eq.6) 
     $       write (19,*) ' ZeroEnergy[kcal/mol]            $wellpen'
         endif
         write (19,*) ' ElectronicLevels[1/cm]           ',nelec
         do ielec = 1, nelec
            write (19,*) eelec(ielec),gelec(ielec)
         enddo
         write (19,*) 'End '

cc here all possible reaction types should be listed
         if (ispecies.eq.1.and.iiso.eq.1) then
            write (19,*) 'End'
         else if (ispecies.eq.1.and.ibeta.eq.1) then
            write (19,*) 'End'
         else if (ispecies.eq.3.and.ipr1.eq.0) then
            write (19,*) ''
         else if (ispecies.eq.3.and.ipr1.eq.1) then
            write (19,*) 'End'
         else if (ispecies.eq.4) then
            write (19,*) ''
         else if (ispecies.eq.1.and.iadd.eq.1) then
            write (19,*) ''
         else if (ispecies.eq.1.and.iabs.eq.1) then
            write (19,*) ''
         else if (ispecies.eq.2.and.iabs.eq.1) then
            write (19,*) ''
         else if (ispecies.eq.2.and.iadd.eq.1) then
            write (19,*) ''
         else
            write (19,*) 'End'
         endif

         if (ispecies.eq.2.or.ispecies.eq.4) then
            if (ispecies.eq.2) write (19,*) 
     $       'GroundEnergy[kcal/mol] 0.0'
            if (ispecies.eq.4) write (19,*) 
     $       'GroundEnergy[kcal/mol] $proden'
            write (19,*) 'End'
c            if (ispecies.eq.2) write (119,*) 
c     $       ' ZeroEnergy[kcal/mol]             0.'
c            if (ispecies.eq.4) write (119,*) 
c     $       ' ZeroEnergy[kcal/mol]             $proden'
c            write (119,*) ' ElectronicLevels[1/cm]           ',nelec
c            do ielec = 1, nelec
c               write (119,*) eelec(ielec),gelec(ielec)
c            enddo
c            write (119,*) 'End '
         endif
c         if ((ispecies.eq.5).or.(ispecies.eq.6)) 
c     $    write (19,*) 'End'
         write (19,*) '!************************************'
c        write (119,*) '!************************************'

cc now update saved geometries at level1 of theory

         if (ispecies.eq.1) then
          open (unit=17,file='./output/reac1_opt.out',status='unknown')
         endif
         if (ispecies.eq.2) then
          open (unit=17,file='./output/reac2_opt.out',status='unknown')
         endif
         if (ispecies.eq.3) then
          open (unit=17,file='./output/prod1_opt.out',status='unknown')
         endif
         if (ispecies.eq.4) then
          open (unit=17,file='./output/prod2_opt.out',status='unknown')
         endif
         if (ispecies.eq.5) then
          open (unit=17,file='./output/wellr_opt.out',status='unknown')
         endif
         if (ispecies.eq.6) then
          open (unit=17,file='./output/wellp_opt.out',status='unknown')
         endif

         ncoord = 3*natom-6
         if(natom.eq.2) ncoord = 1
         if(natom.eq.1) ncoord = 0
         write (17,*) 'opt level1 0'
         do iint = 1 , ncoord
            write (17,*) xint(iint)
         enddo
         close(17)

c     compute fc matrix for HR analysis
         ifreq=0
         ires=0
         ixyz=0
         ired=0

         if (nhind.ne.0) then
            comline1=comline3
            comline2=comline4
            atomlabel(1)=' '
            atomlabel(2)=' '


            if(ilev1code.eq.1) then
               call g09fopt(tau,ntau,natom,natom,numproc,gmem,
     $ coord,vtot_0,vtotr,freq,ifreq,ilin,ismp,comline1,
     $ comline2,icharge,ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ired)

               if (ispecies.eq.1) then
                  command1='cp -f geom.log ./output/reac1_fcmat.log '
               endif
               if (ispecies.eq.2) then
                  command1='cp -f geom.log ./output/reac2_fcmat.log '
               endif
               if (ispecies.eq.3) then
                  command1='cp -f geom.log ./output/prod1_fcmat.log '
               endif
               if (ispecies.eq.4) then
                  command1='cp -f geom.log ./output/prod2_fcmat.log '
               endif
               if (ispecies.eq.5) then
                  command1='cp -f geom.log ./output/wellr_fcmat.log '
               endif
               if (ispecies.eq.6) then
                  command1='cp -f geom.log ./output/wellp_fcmat.log '
               endif
               call commrun(command1)
            else if (ilev1code.eq.2) then
               if (ispecies.eq.1) then
                  command1='cp -f fcmat.log ./output/reac1_fcmat.log '
               endif
               if (ispecies.eq.2) then
                  command1='cp -f fcmat.log ./output/reac2_fcmat.log '
               endif
               if (ispecies.eq.3) then
                  command1='cp -f fcmat.log ./output/prod1_fcmat.log '
               endif
               if (ispecies.eq.4) then
                  command1='cp -f fcmat.log ./output/prod2_fcmat.log '
               endif
               if (ispecies.eq.5) then
                  command1='cp -f fcmat.log ./output/wellr_fcmat.log '
               endif
               if (ispecies.eq.6) then
                  command1='cp -f fcmat.log ./output/wellp_fcmat.log '
               endif
               call commrun(command1)
            endif            
         endif

      else if (inp_type.eq.2) then

c     save fc matrix for HR analysis

         if(ispecies.eq.0)then
            open (unit=18,file='./me_files/wellr_fake.me',
     +           status='unknown') 
            open (unit=28,file='./me_files/wellp_fake.me',
     +           status='unknown') 
         endif
         if(ispecies.eq.0)then
            if (iabs.eq.1) then
               if (nts.eq.1) write (16,*) 'Barrier TS REACS WR'
               if (nts.eq.2) write (16,*) 'Barrier B2 WR PRODS'
               if (nts.eq.3) write (16,*) 'Barrier B2 WR WP'
            endif
            if (iadd.eq.1) then
               if (nts.eq.1) write (16,*) 'Barrier TS REACS WR'
               if (nts.eq.2) write (16,*) 'Barrier B2 WR WP'
            endif
            if (iiso.eq.1) then
               if (nts.eq.1.and.ipr1.eq.0) then
                  write (16,*) 'Barrier TS REACS WP'
               else if(nts.eq.1.and.ipr1.eq.1)then
                  write (16,*) 'Barrier TS REACS PRODS'
               endif
            endif
            if (ibeta.eq.1) then
               if (nts.eq.1) write (16,*) 'Barrier TS REACS PRODS'
            endif
         endif
         if (ispecies.eq.51.or.ispecies.eq.5) then
            write (16,*) 'Well WR'
            write (16,*) 'Species'
         endif
         if (ispecies.eq.61.or.ispecies.eq.6) then
            write (16,*) 'Well WP'
            write (16,*) 'Species'
         endif

         if(ispecies.eq.0)then
            write (18,*) 'Well WR'
            write (28,*) 'Well WP'
            write (18,*) 'Species'
            write (28,*) 'Species'
            write (18,*) 'RRHO     ! fake well'
            write (28,*) 'RRHO     ! fake well'
            write (18,*) 'Geometry[angstrom]    ',natom
            write (28,*) 'Geometry[angstrom]    ',natom
         endif
         write (16,*) 'RRHO     ! transition state'
         write (16,*) 'Geometry[angstrom]    ',natom

c         command1='newzmat       
c     +    -ichk -oxyz tmp.chk geom.xyz           '
c         call commrun(command1)
         open (unit=98,file='headgeom.tmp',status='unknown')
         write(98,*)natom
         write(98,*)vtotr
         close(98)
         command1='cat headgeom.tmp geom.xyz > geom1.xyz '
         call commrun(command1)

         open (unit=99,status='unknown')
         rewind (99)
         write (99,2112) nameout
         rewind (99)
         read (99,2020) command1
         call commrun(command1)         
 2112    format ('cp -f geom1.xyz geoms/'A8,'.xyz')

         open (unit=99,status='unknown')
         rewind (99)
         write (99,2122) nameout
         rewind (99)
         read (99,1020) command1
         close(99)
         call commrun(command1)         
 2122    format ('cp -f geom.log geoms/'A8,'.log')
         if (ilev1code.eq.2) then
            open (unit=99,status='unknown')
            rewind (99)
            write (99,2123) nameout
            rewind (99)
            read (99,2020) command1
            close(99)
            call commrun(command1)         
 2123       format ('cp -f molpro.molden geoms/'A8,'.molden')
         endif

 2020    format (A70)

         open (unit=100,FILE='geom.xyz',status='unknown')
         do j=1,natom
            read (100,1002) geomline
            write (16,*) geomline
            if(ispecies.eq.0)write (18,*) geomline
            if(ispecies.eq.0)write (28,*) geomline
         enddo
 1002    FORMAT (A85)
         close (unit=100,status='keep')
c         command1='cat reac_1.me geom.xyz > temp.me'
c         call commrun(command1)
c         command1='cp -f temp.me reac_1.me'
c         call commrun(command1)
         write (16,*) '	Core 	RigidRotor'
         write (16,*) '		SymmetryFactor ',symf
         write (16,*) '	End'
         if(ispecies.eq.0) then
            write (18,*) '	Core 	RigidRotor'
            write (18,*) '		SymmetryFactor ',symf
            write (18,*) '	End'
            write (28,*) '	Core 	RigidRotor'
            write (28,*) '		SymmetryFactor ',symf
            write (28,*) '	End'
         endif

         nfreq=3*natom-6
         if(ilin.eq.1)nfreq=nfreq+1
         if (natom.eq.2) nfreq=1

         jt = 0
         do j = 1, nfreq
            if (freq(j).gt.0.) then 
               jt = jt + 1
               freqp(jt) = freq(j)*sclfr
            endif
         enddo

         write (19,*) '    Frequencies[1/cm] ',jt

         if(ispecies.eq.0)then
            write (18,*) '    Frequencies[1/cm] ',jt
            write (28,*) '    Frequencies[1/cm] ',jt
         endif

         write (19,8010) (freqp(j),j=1,jt)

         if(ispecies.eq.0)then
            write (18,8010) (freqp(j),j=1,jt)
            write (28,8010) (freqp(j),j=1,jt)
         endif

         zpe=0
         do j=1,jt
            zpe=zpe+freqp(j)
         enddo
cc       divide by two and convert to kcal/mol
         zpe=zpe/(cautoicm*2.d0)
c        zpe=zpe*29979200000*6.626E-34/2.*6.022E23/4.184/1000/627.503
         write (20,*) zpe

         if(ispecies.eq.0)then
            write (18,*) ' ZeroEnergy[kcal/mol]            0. '
            write (28,*) ' ZeroEnergy[kcal/mol]            0. '
            write (18,*)' ElectronicLevels[1/cm]           1'
            write (18,*)'     0.0           ',ispin
            write (28,*)' ElectronicLevels[1/cm]           1'
            write (28,*)'     0.0           ',ispin
            write (18,*)'End '
            write (18,*)'End '
            write (18,*)'!************************************'
            write (28,*)'End '
            write (28,*)'End '
            write (28,*)'!************************************'
         endif
         
         if(ispecies.eq.0)
     $      write (19,*) ' ZeroEnergy[kcal/mol]            $tsen '
         if (ispecies.eq.51.or.ispecies.eq.5) 
     $       write (19,*) ' ZeroEnergy[kcal/mol]            $wellren'
         if (ispecies.eq.61.or.ispecies.eq.6) 
     $       write (19,*) ' ZeroEnergy[kcal/mol]            $wellpen'
c         endif

         write (19,*) ' ElectronicLevels[1/cm]           ',nelec
         do ielec = 1, nelec
            write (19,*) eelec(ielec),gelec(ielec)
         enddo
         if(ispecies.eq.0.and.ifrozrts.ne.1)then
            if(imdtunn.ne.1)then
               write (19,*)'      Tunneling    Eckart'
               write (19,*)'   ImaginaryFrequency[1/cm]   ',abs(freq(1))
               write (19,*)'   WellDepth[kcal/mol]       $wdepfor'
               write (19,*)'   WellDepth[kcal/mol]       $wdepback'
            else
               write (19,*)'      Tunneling    Read'
               write (19,*)'   CutoffEnergy[1/cm]       2500'
               write (19,*)'   ImaginaryFrequency[1/cm]   ',abs(freq(1))
               write (19,*)'   File imactint.dat'
            endif
         endif
         write (19,*)'End '
         write (19,*)'End '
c         write (19,*)'End '
c         write (19,*)'End '
         write (19,*)'!************************************'

c save updated TS geometry

         if(ispecies.eq.0)then
            open (unit=37,file='./output/ts_opt.out',status='unknown')
            write (37,*) 'TS_z-matrix level1'
            do iatom = 1 , natomt
               write (37,*) atomlabel(iatom)
            enddo
            ncoord = 3*natom-6
            do iint = 1 , ncoord
               write (37,*) intcoor(iint),xint(iint)
            enddo
            close (unit=37,status='keep')
         endif
         if(ispecies.eq.51)then
            open (unit=37,file='./output/wellr_opt.out',
     $            status='unknown')
            write (37,*) 'wellr_z-matrix level1'
            do iatom = 1 , natomt
               write (37,*) atomlabel(iatom)
            enddo
            ncoord = 3*natom-6
            do iint = 1 , ncoord
               write (37,*) intcoor(iint),xint(iint)
            enddo
            close (unit=37,status='keep')
         endif
         if(ispecies.eq.61)then
            open (unit=37,file='./output/wellp_opt.out',
     $            status='unknown')
            write (37,*) 'wellp_z-matrix level1'
            do iatom = 1 , natomt
               write (37,*) atomlabel(iatom)
            enddo
            ncoord = 3*natom-6
            do iint = 1 , ncoord
               write (37,*) intcoor(iint),xint(iint)
            enddo
            close (unit=37,status='keep')
         endif

c     compute fc matrix for HR analysis

         ires=0
         ixyz=0
         ifreq=0
         ired=0

         if (nhind.ne.0) then
            comline1=comline3
            comline2=comline4
            atomlabel(1)=' '
            atomlabel(2)=' '

            if(ilev1code.eq.1) then
               call g09fopt(tau,ntau,natom,natom,numproc,gmem,
     $ coord,vtot_0,vtotr,freq,ifreq,ilin,ismp,comline1,
     $ comline2,icharge,ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ired)

               if (ispecies.eq.51) then
                  command1='cp -f geom.log ./output/wellr_fcmat.log '
               endif
               if (ispecies.eq.61) then
                  command1='cp -f geom.log ./output/wellp_fcmat.log '
               endif
               if (ispecies.eq.0) then
                  command1='cp -f geom.log ./output/ts_fcmat.log '
               endif
               call commrun(command1)
            else if (ilev1code.eq.2) then
               if (ispecies.eq.51) then
                  command1='cp -f fcmat.log ./output/wellr_fcmat.log '
               endif
               if (ispecies.eq.61) then
                  command1='cp -f fcmat.log ./output/wellp_fcmat.log '
               endif
               if (ispecies.eq.0) then
                  command1='cp -f fcmat.log ./output/ts_fcmat.log '
               endif
               call commrun(command1)
            endif
         endif
      endif

      close (unit=16,status='keep')
      close (unit=19,status='keep')
      close (unit=20,status='keep')
      close (unit=66,status='keep')

cc save me_files of frequencies so that they are not overwritten by 1Dhscan


      if (ispecies.eq.1) then
         command1='cp -f me_files/reac1_ge.me
     $            me_files/reac1_1dge.me '
         call commrun(command1)
         command1='cp -f ./me_files/reac1_fr.me 
     $                   ./me_files/reac1_unpfr.me'
      endif
      if (ispecies.eq.2) then
         command1='cp -f me_files/reac2_ge.me
     $            me_files/reac2_1dge.me '
         call commrun(command1)
         command1='cp -f ./me_files/reac2_fr.me 
     $                   ./me_files/reac2_unpfr.me'
      endif
      if (ispecies.eq.3) then
         command1='cp -f me_files/prod1_ge.me
     $            me_files/prod1_1dge.me '
         call commrun(command1)
         command1='cp -f ./me_files/prod1_fr.me 
     $                   ./me_files/prod1_unpfr.me'
      endif
      if (ispecies.eq.4) then
         command1='cp -f me_files/prod2_ge.me
     $            me_files/prod2_1dge.me '
         call commrun(command1)
         command1='cp -f ./me_files/prod2_fr.me 
     $                   ./me_files/prod2_unpfr.me'
      endif
      if (ispecies.eq.5) then
         command1='cp -f me_files/wellr_ge.me
     $            me_files/wellr_1dge.me '
         call commrun(command1)
         command1='cp -f ./me_files/wellr_fr.me 
     $                   ./me_files/wellr_unpfr.me'
      endif
      if (ispecies.eq.6) then
         command1='cp -f me_files/wellp_ge.me
     $            me_files/wellp_1dge.me '
         call commrun(command1)
         command1='cp -f ./me_files/wellp_fr.me 
     $                   ./me_files/wellp_unpfr.me'
      endif
      if (ispecies.eq.51) then
         command1='cp -f me_files/wellr_ge.me
     $            me_files/wellr_1dge.me '
         call commrun(command1)
         command1='cp -f ./me_files/wellr_fr.me 
     $                   ./me_files/wellr_unpfr.me'
      endif
      if (ispecies.eq.61) then
         command1='cp -f me_files/wellp_ge.me
     $            me_files/wellp_1dge.me '
         call commrun(command1)
         command1='cp -f ./me_files/wellp_fr.me 
     $                   ./me_files/wellp_unpfr.me'
      endif
      if (ispecies.eq.0) then
         command1='cp -f me_files/ts_ge.me
     $            me_files/ts_1dge.me '
         call commrun(command1)
         command1='cp -f ./me_files/ts_fr.me 
     $                   ./me_files/ts_unpfr.me'
      endif
      
      call commrun(command1)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine onedtau(ispecies)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

cc    this is the hindered rotor section

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'


      dimension dhindmn(nhindmx),dhindmx(nhindmx),
     $ freq(nmdmx),nhindsteps(nhindmx),nhindsymm(nhindmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim),xintt(3*natommx),xint_save(3*natommx)
      dimension vref(noptmx),tauo(ntaumx,noptmx),
     $ xinto(3*natommx,noptmx),dist(noptmx)
      dimension tauopt(ntaumx)
      dimension tau(ntaumx)
      dimension iatomtype(natommx)
      dimension itop_ind(natommx,nhindmx),
     $          itop_ind_nodum(natommx,nhindmx)
      dimension idummy(natommx),
     $          ngroup(nhindmx,natommx)
      dimension igrouptot(nhindmx),ipivotA(nhindmx),ipivotB(nhindmx)
      dimension taumn(ntaumx),taumx(ntaumx)
      dimension gelec(nelecmx),eelec(nelecmx)

      character*70 comline1,comline2,comline3,comline4,comsave1,comsave2
      character*60 atomlabel(natommx)
      character*4 atomname(natommx)
      character*4 atomconn(natommx)
      character*4 atomcoo(natommx)
      character*4 pivotA(nhindmx),pivotB(nhindmx)
      character*4 dicheck1(nhindmx),dicheck2(nhindmx)
      character*30 intcoor(3*natommx)
      character*30 intcoori(3*natommx)
      character*30 intcoort(3*natommx)
      character*20 bislab(ntaumx)
      character*20 hindlab(nhindmx)
      character*20 wordlab
      character*1 atfirstlet
      character*160 command1
      character*160 commandcopy
      character*160 word_hr
      character*30 gmem

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

cc initialize index for only reactant calculation
      ireactonly=0
      if(iabs.eq.0.and.iadd.eq.0.and.iiso.eq.0.and.ibeta.eq.0)then
         ireactonly=1
      endif
c input data
cs added option for species dependent hr level

      call Lineread(0)
      commandcopy='cp -f ./data/onedtau_molpro.dat 
     $       onedtau_molpro.dat'

      open (unit=21,file='./data/theory.dat',status='unknown')

      if (ispecies.eq.1) then
         open (unit=15,file='./data/reac1.dat',status='old')
         open (unit=16,file='./output/reac1_hr.out',status='unknown')
         open (unit=17,file='./output/reac1_opt.out',status='unknown')
         if(irecov.ne.1) then
            open (unit=19,file='./me_files/reac1_hr.me',
     $ status='unknown')
         endif
c         open (unit=119,file='./me_files/rpst_hr1.me',status='unknown')
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_1'
         do while (WORD.NE.'HIND_ROTOR_1')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/onedtau_reac1_molpro.dat 
     $       onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.2) then
         open (unit=15,file='./data/reac2.dat',status='unknown')
         open (unit=16,file='./output/reac2_hr.out',status='unknown')
         open (unit=17,file='./output/reac2_opt.out',status='unknown')
         if(irecov.ne.1) then
            open (unit=19,file='./me_files/reac2_hr.me',
     $status='unknown')
         endif
c         open (unit=119,file='./me_files/rpst_hr2.me',status='unknown')
         inp_type=1
c        return
cs I removed this return - it would be fine for H,O,OH,O2,CH3,HO2 but
cs will create problems for larger radicals
cs I don't know if this removal will create other problems
         rewind(21)
         WORD_HR='HIND_ROTOR_2'
         do while (WORD.NE.'HIND_ROTOR_2')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/onedtau_reac2_molpro.dat 
     $       onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.3) then
         open (unit=15,file='./data/prod1.dat',status='old')
         open (unit=16,file='./output/prod1_hr.out',status='unknown')
         open (unit=17,file='./output/prod1_opt.out',status='unknown')
         if(irecov.ne.1) then
            open (unit=19,file='./me_files/prod1_hr.me',
     $ status='unknown')
         endif
c         open (unit=119,file='./me_files/ppst_hr1.me',status='unknown')
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_3'
         do while (WORD.NE.'HIND_ROTOR_3')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/onedtau_prod1_molpro.dat 
     $       onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.4) then
         open (unit=15,file='./data/prod2.dat',status='old')
         open (unit=16,file='./output/prod2_hr.out',status='unknown')
         open (unit=17,file='./output/prod2_opt.out',status='unknown')
         if(irecov.ne.1) then
            open (unit=19,file='./me_files/prod2_hr.me',
     $ status='unknown')
         endif
c         open (unit=119,file='./me_files/ppst_hr2.me',status='unknown')
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_4'
         do while (WORD.NE.'HIND_ROTOR_4')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/onedtau_prod2_molpro.dat 
     $       onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.5) then
         open (unit=15,file='./data/wellr.dat',status='old')
         open (unit=16,file='./output/wellr_hr.out',status='unknown')
         open (unit=17,file='./output/wellr_opt.out',status='unknown')
         if(irecov.ne.1) then
            open (unit=19,file='./me_files/wellr_hr.me',
     $ status='unknown')
         endif
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_5'
         do while (WORD.NE.'HIND_ROTOR_5')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/onedtau_wellr_molpro.dat 
     $       onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.51) then
         open (unit=15,file='./data/wellr.dat',status='old')
         open (unit=16,file='./output/wellr_hr.out',status='unknown')
         open (unit=17,file='./output/wellr_opt.out',status='unknown')
         if(irecov.ne.1) then
            open (unit=19,file='./me_files/wellr_hr.me',
     $    status='unknown')
         endif
         inp_type=2
         rewind(21)
         WORD_HR='HIND_ROTOR_51'
         do while (WORD.NE.'HIND_ROTOR_51')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/onedtau_wellr_molpro.dat 
     $       onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.6) then
         open (unit=15,file='./data/wellp.dat',status='old')
         open (unit=16,file='./output/wellp_hr.out',status='unknown')
         open (unit=17,file='./output/wellp_opt.out',status='unknown')
         if(irecov.ne.1) then
            open (unit=19,file='./me_files/wellp_hr.me',
     $  status='unknown')
         endif
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_6'
         do while (WORD.NE.'HIND_ROTOR_6')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/onedtau_wellp_molpro.dat 
     $       onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.61) then
         open (unit=15,file='./data/wellp.dat',status='old')
         open (unit=16,file='./output/wellp_hr.out',status='unknown')
         open (unit=17,file='./output/wellp_opt.out',status='unknown')
         if(irecov.ne.1) then
            open (unit=19,file='./me_files/wellp_hr.me',
     $  status='unknown')
         endif
         inp_type=2
         rewind(21)
         WORD_HR='HIND_ROTOR_61'
         do while (WORD.NE.'HIND_ROTOR_61')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/onedtau_wellp_molpro.dat 
     $       onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.11) then
         open (unit=15,file='./data/reacs.dat',status='old')
         open (unit=16,file='./output/reacs_hr.out',status='unknown')
         open (unit=17,file='./output/reacs_opt.out',status='unknown')
         if(irecov.ne.1) then
            open (unit=19,file='./me_files/reacs_hr.me'
     $        ,status='unknown')
         endif
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_11'
         do while (WORD.NE.'HIND_ROTOR_11')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/onedtau_reacs_molpro.dat 
     $       onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.12) then
         open (unit=15,file='./data/prods.dat',status='old')
         open (unit=16,file='./output/prods_hr.out',status='unknown')
         open (unit=17,file='./output/prods_opt.out',status='unknown')
         if(irecov.ne.1) then
            open (unit=19,file='./me_files/prods_hr.me',
     $ status='unknown')
         endif
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_12'
         do while (WORD.NE.'HIND_ROTOR_12')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/onedtau_prods_molpro.dat 
     $       onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.0) then
         open (unit=15,file='./data/ts.dat',status='unknown')
         open (unit=16,file='./output/ts_hr.out',status='unknown')
         open (unit=17,file='./output/ts_opt.out',status='unknown')
         if(irecov.ne.1) then
            open (unit=19,file='./me_files/ts_hr.me',status='unknown')
         endif
         inp_type=2
         rewind(21)
         WORD_HR='HIND_ROTOR_TS'
         do while (WORD.NE.'HIND_ROTOR_TS')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/onedtau_ts_molpro.dat 
     $       onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.100.or.ispecies.eq.101) then
         open (unit=14,file='./irc_files/hrcc.dat',status='unknown')
         open (unit=15,file='./data/ts.dat',status='unknown')
         open (unit=16,file='./output/ts_hr.out',status='unknown')
         open (unit=17,file='./output/ts_opt.out',status='unknown')
         if(irecov.ne.1) then
            open (unit=19,file='./irc_files/hrcc.me',status='unknown')
         endif
         inp_type=2
         rewind(21)
         WORD_HR='HIND_ROTOR_TS'
         do while (WORD.NE.'HIND_ROTOR_TS')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/onedtau_ts_molpro.dat 
     $       onedtau_molpro.dat'
         endif
      endif
 900  continue
      rewind(21)

      ifreq=0
      noptg = 0
      ires=0
cc we assume that a linear molecule has no rotational dof         
      ilin=0

      do while (WORD.NE.'NHIND')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (16,*) 'hind rotors must be defined'
            stop
         endif
      enddo
c      call LineRead (15)
      read (15,*) nhind
      read (15,*)
      ntau_fr = nhind
      do ihind = 1 , nhind
         read (15,*) hindlab(ihind),dhindmn(ihind),dhindmx(ihind),
     +              nhindsteps(ihind),nhindsymm(ihind)
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) hindlab(ihind)
         rewind (99)
         call LineRead (99)
         hindlab(ihind)=WORD
         close (unit=99,status='keep')
      enddo
      if (nhind.eq.0) then
         close (unit=15,status='keep')
         close (unit=16,status='keep')
         close (unit=17,status='keep')
         if(irecov.ne.1) close (unit=19,status='keep')
         word=' '
         return
      endif

      do while (WORD.NE.'NELEC')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (66,*) 'electronic states must be defined'
            stop
         endif
      enddo
      read (15,*) nelec
      do ielec = 1 , nelec
         read (15,*) eelec(ielec),gelec(ielec)
      enddo
      rewind(15)


c read input of type 1

      if (inp_type.eq.1) then

c read level of theory for hindered rotor scan

         do while (WORD.NE.'NATOM')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'natom must be defined'
               stop
            endif
         enddo
         read (15,*) natom,natomt
         rewind(15)

         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'sampling coordinates must be defined'
               stop
            endif
         enddo
         read (15,*) ntau
         rewind(15)


         do while (WORD.NE.'CHARGE')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'charge and spin must be defined'
               stop
            endif
         enddo
         read (15,*) icharge,ispin
         do iatom = 1 , natomt
            read (15,'(A60)') atomlabel(iatom)
         enddo
         rewind(15)

         do while (WORD.NE.WORD_HR)
            call LineRead (21)
            if (WORD.EQ.'END') then
               write (16,*) 'hind rotors must be defined'
               write (16,*) 'in file theory.dat'
               stop
            endif
         enddo
         if(word2.eq.'G09') then
            ilev1code=1
            read (21,'(A70)') comline1
            read (21,'(A70)') comline2
         else if(word2.eq.'MOLPRO')then
            ilev1code=2
            call commrun(commandcopy)
         else
            write (16,*) 'hind rot of theory must be either'
            write (16,*) 'g09 or molpro in theory.dat'
            stop
         endif
         if(word3.eq.'RESTART') then
            open (unit=99,status='unknown')
            write (99,*) word4
            rewind (99)
            read (99,*) ires
            close (unit=99,status='keep')
         endif

         if (idebug.ge.2) write (16,*) ' comline test',comline1,comline2
         close(21)

         if (natom.ne.1) then
            do while (WORD.NE.'INTCOOR')
               call LineRead (15)
               if (WORD.EQ.'END') then
                  write (16,*) 'internal coordinates must be defined'
                  stop
               endif
            enddo
            ncoord = 3*natom-6-ntau
            if (natom.eq.1) ncoord = 0
            if (natom.eq.2) ncoord = 1
            do icoord = 1 , ncoord
               call LineRead (15)
               intcoori(icoord) = word
            enddo
            rewind(15)
         endif

         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'sampling coordinates must be defined'
               stop
            endif
         enddo
         read (15,*) ntau
         if (ntau.gt.ntaumx) then
            stop
         endif
         if (ntau.ne.0) then
            read (15,*)
            do itau = 1 , ntau
               read (15,*) bislab(itau),taumn(itau),taumx(itau)
               write (16,*) bislab(itau),taumn(itau),taumx(itau)
               open (unit=99,status='unknown')
               rewind (99)
               write (99,*)bislab(itau)
               rewind (99)
               call LineRead (99)
               icoord=ncoord+itau
               intcoori(icoord)=WORD
               close (unit=99,status='keep')
            enddo
         endif
         rewind(15)


         ncoord = 3*natom-6
         read (17,*)
         do icoord = 1 , ncoord
            call LineRead (17)
            OPEN (unit=99,status='unknown')
            REWIND (99)
            WRITE (99,1000) WORD
 1000       FORMAT (A30)
            REWIND (99)
            READ (99,*) xinti(icoord)
            close (unit=99,status='keep')
         enddo
         close(17)

         if (idebug.ge.2) write (16,*) 'past z-matrix'
      endif


c now read input of type 2

      if (inp_type.eq.2) then


         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'sampling coords of reac1 must be defined'
               stop
            endif
         enddo
         read (15,*) ntau
         rewind 15

         open (unit=25,file='./data/reac1.dat',status='old')
         word=' '
         do while (WORD.NE.'NTAU')
            call LineRead (25)
            if (WORD.EQ.'END') then
               write (16,*) 'sampling coords of reac1 must be defined'
               stop
            endif
         enddo
         read (25,*) ntau1
         rewind 25

c     natomt is to account for dummy atoms
         do while (WORD.NE.'NATOM')
            call LineRead (25)
            if (WORD.EQ.'END') then
               write (16,*) 'natom in reac1 must be defined'
               stop
            endif
         enddo
         read (25,*) natom1,natomt1
         close (25)

cc get data from react2 file

         if(iabs.eq.1.or.iadd.eq.1)then
            open (unit=25,file='./data/reac2.dat',status='old')

            do while (WORD.NE.'NTAU')
               call LineRead (25)
               if (WORD.EQ.'END') then
                  write (66,*) 'sampl coords of reac2 must be defined'
                  stop
               endif
            enddo
            read (25,*) ntau2
            rewind 25

c     natomt is to account for dummy atoms
            do while (WORD.NE.'NATOM')
               call LineRead (25)
               if (WORD.EQ.'END') then
                  write (66,*) 'natom must be defined'
                  stop
               endif
            enddo
            read (25,*) natom2,natomt2
            close (25)
         endif

         if(iadd.eq.1.or.iabs.eq.1)then
            natom = natom1+natom2
         else if (iiso.eq.1.or.ibeta.eq.1) then
            natom = natom1
         endif
         if (iadd.eq.1) natomt = natomt1+natomt2
         if (iabs.eq.1) natomt = natomt1+natomt2+1
         if (iiso.eq.1) natomt = natomt1
         if (ibeta.eq.1) natomt = natomt1
c        natomt = natomt1+natomt2+1

        do while (WORD.NE.'CHARGE')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'charge and spin must be defined'
               stop
            endif
         enddo
         read (15,*) icharge,ispin
c         rewind(15)
         close (unit=15,status='keep')
         
         if(ispecies.eq.0) then
            do while (WORD.NE.WORD_HR)
               call LineRead (21)
               if (WORD.EQ.'END') then
                  write (16,*) 'hind rotors must be defined for TS'
                  write (16,*) 'in file theory.dat'
                  stop
               endif
            enddo
         else if((ispecies.eq.100).or.(ispecies.eq.51).or.
     $    (ispecies.eq.61).or.(ispecies.eq.101)) then
            do while (WORD.NE.WORD_HR)
               call LineRead (21)
               if (WORD.EQ.'END') then
                  write (16,*) 'hind rotors must be defined'
                  write (16,*) 'in file theory.dat'
                  stop
               endif
            enddo
         endif
         if(word2.eq.'G09') then
            ilev1code=1
            read (21,'(A70)') comline1
            read (21,'(A70)') comline2
            if(ispecies.eq.0.or.ispecies.eq.100.or.ispecies.eq.101) then
               read (21,'(A70)') comline3
               read (21,'(A70)') comline4
               if(ispecies.eq.100.or.ispecies.eq.101) then
                  comline1=comline3
                  comline2=comline4
               endif
            endif
         else if(word2.eq.'MOLPRO')then
            ilev1code=2
            call commrun(commandcopy)
         else
            write (16,*) 'hind rot code must be either'
            write (16,*) 'g09 or molpro in theory.dat'
            stop
         endif
         if(word3.eq.'RESTART') then
            open (unit=99,status='unknown')
            write (99,*) word4
            rewind (99)
            read (99,*) ires
            close (unit=99,status='keep')
         endif

         if (idebug.ge.2) write (16,*) ' comline test',comline1,comline2
         if (idebug.ge.2) write (16,*) ' test2       ',comline3,comline4
         close(21)

c now read z-mat input

         read (17,*)
         if (idebug.ge.2) write (6,*) ' starting gaussian input'
         do iatom = 1 , natomt
            read (17,'(A60)') atomlabel(iatom)
         enddo

cc read coordinate names

         ncoord = 3*natom-6
         do iint = 1 , ncoord
            read (17,*) intcoori(iint),xinti(iint)
         enddo
         close (unit=17,status='keep')

         if (idebug.ge.2) write (16,*) 'past z-matrix'

cc now read abstraction site from grid input file

         open (unit=18,file='./data/ts.dat',status='old')

         do while (WORD.NE.'ISITE')
            call LineRead (18)
            if (WORD.EQ.'END') then
               write (16,*) 'abstraction site must be defined'
               stop
            endif
         enddo
         read (18,*) isite,jsite,ksite

         close (unit=18,status='keep')

cc if called from IRC read additional input
         if(ispecies.eq.100) then
            read(14,*)ihr_at1,disthr_at1
            read(14,*)ihr_at2,disthr_at2
            close(14)
         else if (ispecies.eq.101) then
            read(14,*)ihr_at1,disthr_at1
            close(14)
         endif
      endif


cc scan one hindered rotor at a time

c      nint = 3*natom-1-6

cc identify pivot atoms and connectivity in z-mat
      
      OPEN (unit=99,status='unknown')
      REWIND (99)
      write (99,*) atomlabel(1)
      REWIND (99)
      read (99,*) atomname(1)
      close (99)
      atomconn(1)='0000'
      atomcoo(1)='0000'
      
      OPEN (unit=99,status='unknown')
      REWIND (99)
      write (99,*) atomlabel(2)
      REWIND (99)
      read (99,*) atomname(2),atomconn(2),atomcoo(2)
      close (99)

      OPEN (unit=99,status='unknown')
      REWIND (99)
      write (99,*) atomlabel(3)
      REWIND (99)
      read (99,*) atomname(3),atomconn(3),atomcoo(3)
      close (99)

      do iatom = 4 , natomt
         OPEN (unit=99,status='unknown')
         REWIND (99)
         write (99,*) atomlabel(iatom)
         REWIND (99)
         read (99,*) atomname(iatom),atomconn(iatom),atomcoo(iatom)
     $        ,word2,word,word,wordlab
         REWIND (99)
         write (99,*) wordlab
         close (99)
         OPEN (unit=99,status='unknown')
         call LineRead (99)
         close (99)
         do ihind = 1 , nhind
            if (word.eq.hindlab(ihind)) then 
               OPEN (unit=99,status='unknown')
               REWIND (99)
               write (99,*) atomlabel(iatom)
               REWIND (99)
               read (99,*) dicheck1(ihind),pivotA(ihind),word,
     $        pivotB(ihind),word,dicheck2(ihind),word
            endif   
         enddo
      enddo

cc  assign connectivity of atom1 to atom2

      atomconn(1)=atomname(2)

cc      now assign connectivities and PivotB of atoms after dummy for TS

      if (inp_type.eq.2) then

         if(iabs.eq.1) then
cc for abstraction
cc initialize dummy atom and following atom connectivity and pivotB name
            do ihind=1 , nhind
               if(pivotB(ihind).eq.atomconn(natomt1+1))then
                  pivotB(ihind)=atomname(isite)
               endif
               if(pivotB(ihind).eq.atomconn(natomt1+2))then
                  pivotB(ihind)=atomname(isite)
               endif
            enddo
            atomconn(natomt1+1)=atomname(isite)
            atomconn(natomt1+2)=atomname(isite)
         endif
cc for addition
cc no dummy atom in this case         
         if(iadd.eq.1) then
            do ihind=1 , nhind
               if(pivotB(ihind).eq.atomconn(natomt1+1))then
                  pivotB(ihind)=atomname(isite)
               endif
               if (hindlab(ihind).eq.'BABS1') then 
                  pivotA(ihind)=atomname(isite)
                  pivotB(ihind)=atomname(jsite)
               endif
            enddo
            atomconn(natomt1+1)=atomname(isite)
         endif
      endif
cc now convert all labels to uppercase

      do iatom = 1 , natomt

         OPEN (unit=99,status='unknown')
         REWIND (99)
         write (99,*) atomname(iatom),' ',atomconn(iatom),' ',
     $               atomcoo(iatom)
         REWIND (99)
         call LineRead(99)
         atomname(iatom)=word
         atomconn(iatom)=word2
         atomcoo(iatom)=word3
         close (99)
      enddo

      do ihind=1, nhind
         OPEN (unit=99,status='unknown')
         REWIND (99)
         write (99,*)pivotA(ihind),' ',pivotB(ihind)
         REWIND (99)
         call LineRead(99)
         pivotA(ihind)=word
         pivotB(ihind)=word2
         close (99)

         OPEN (unit=99,status='unknown')
         REWIND (99)
         write (99,*)dicheck1(ihind),' ',dicheck2(ihind)
         REWIND (99)
         call LineRead(99)
         dicheck1(ihind)=word
         dicheck2(ihind)=word2
         close (99)
      enddo

cc
cc now determine the atoms belonging to the top
cc

      do ihind = 1 , nhind
         do iatom = 1 , natomt
            itop_ind(iatom,ihind)=0
         enddo
      enddo

      do ihind = 1 , nhind

         do iatom = 1 , natomt
            if (atomname(iatom).eq.pivotA(ihind)) 
     +        itop_ind(iatom,ihind)=1
            if (atomname(iatom).eq.pivotB(ihind)) 
     +        itop_ind(iatom,ihind)=2
         enddo

         istep=0
 1010    continue
         do iatom = 1 , natomt
            counter=0
            do jatom =1, natomt
               if (itop_ind(jatom,ihind).ne.0) then
c               if (iatom.eq.1) counter=0.
                  counter=counter+1
               endif
            enddo
            if (counter.eq.natomt) goto 1020
c            if (itop_ind(iatom,ihind).ne.0) goto 1011

            do iconn = 1, natomt
               if (atomconn(iatom).eq.atomname(iconn)) then
                  if (itop_ind(iconn,ihind).ne.0) then
                     if (itop_ind(iatom,ihind).eq.0) then
                        itop_ind(iatom,ihind)=itop_ind(iconn,ihind)
                     endif
                  endif

cc do it also the other way around
                  if (itop_ind(iatom,ihind).ne.0) then
                     if (itop_ind(iconn,ihind).eq.0) then
                        itop_ind(iconn,ihind)=itop_ind(iatom,ihind)
                     endif
                  endif

               endif
            enddo
 1011       continue
         enddo
         istep=istep+1
         if (istep.eq.1000) goto 1021
         goto 1010 
 1021    continue
         write (7,*) 'something did not work here'
         write (7,*) 'in sub onedtau'
         write (7,*) 'in the assignment of atoms to tops'
         write (7,*) 'check the me input file'
 1020    continue
         if (idebug.ge.2) then
            write (16,*)'pivA pivB at '
     $           ,pivotA(ihind),pivotB(ihind)
            do iatom = 1 , natomt
               write (16,*)'iatom atomconn',iatom,atomconn(iatom)
            enddo
            do iatom = 1 , natomt
               write (16,*)'iatom top ind',iatom,itop_ind(iatom,ihind)
            enddo
         endif
      enddo

c         stop

      do ihind = 1 , nhind
         itopcheck1=0
         itopcheck2=0
         do iatom = 1 , natomt
            if (dicheck1(ihind).eq.atomname(iatom)) then
               itopcheck1=itop_ind(iatom,ihind)
            endif
            if (dicheck2(ihind).eq.atomname(iatom)) then
               itopcheck2=itop_ind(iatom,ihind)
            endif
         enddo

         if (itopcheck1.eq.itopcheck2) then
           write (7,*) 'there is something wrong with this rotor'
           write (7,*) hindlab(ihind)
           write (7,*) 'ispecies is ', ispecies
           write (7,*) 'check the input file'
           write (7,*)'dicheck1 is ',dicheck1(ihind),'top is',
     $      itopcheck1
           write (7,*)'dicheck2 is ',dicheck2(ihind),'top is',
     $      itopcheck2
         endif
      enddo

cc
cc now inizialize the index for dummy atoms
cc
      do iatom = 1 , natomt
         idummy(iatom)=0
      enddo
      do iatom = 1 , natomt
         OPEN (unit=99,status='unknown')
         REWIND (99)
         write (99,1100) atomname(iatom)
         REWIND (99)
         read (99,1100) atfirstlet
         if (atfirstlet.eq.'X') idummy(iatom)=1
         close (99)
      enddo

 1100 FORMAT (A1)

cc
cc now update list of atoms in top taking out dummies
cc

      do ihind=1,nhind
         iatomup=0
         do iatom = 1 , natomt
            if (idummy(iatom).ne.1) then
               iatomup=iatomup+1
               itop_ind_nodum(iatomup,ihind)=itop_ind(iatom,ihind)
            endif
            if (atomname(iatom).eq.pivotA(ihind)) then
               ipivotA(ihind)=iatomup
            endif
            if (atomname(iatom).eq.pivotB(ihind)) then
               ipivotB(ihind)=iatomup
            endif

         enddo
      enddo
      if (iatomup.ne.natom) then
         write (7,*) 'there is a problem in hind rot section'
         write (7,*) 'with renumb to take out dummy atoms'
         write (7,*) 'ispecies is ', ispecies
      endif


cc
c now determine atoms in first top
cc

      do ihind=1,nhind
         do igroup=1,natom
            ngroup(ihind,igroup)=0
         enddo
         igrouptot(ihind)=0
      enddo

      do ihind=1,nhind
         igroup=0
         do iatom = 1 , natom
            if (iatom.ne.ipivotA(ihind).and.iatom.ne.ipivotB(ihind)) 
     $       then
               if (itop_ind_nodum(iatom,ihind).eq.1) then
                  igroup=igroup+1
                  ngroup(ihind,igroup)=iatom
               endif
            endif
         enddo
         igrouptot(ihind)=igroup
      enddo

c      do ihind=1,nhind
c      do igroup=1,igrouptot(ihind)
c      write (*,*)'hind ind atom',ihind,igroup,ngroup(ihind,igroup)
c      enddo
c      enddo

cc if ispecies is 100 move ihr_at1 and ihr_at2 bond coords to last two 
cc positions
      if(ispecies.eq.100) then

         index=0
         ncoord = 3*natom-6

         xintt(ncoord-1)=disthr_at1
         xintt(ncoord)=disthr_at2
         intcoort(ncoord-1)=atomcoo(ihr_at1)
         intcoort(ncoord)=atomcoo(ihr_at2)

c         write(*,*)'at index 1',ihr_at1
c         write(*,*)'at index 2',ihr_at2
c         write(*,*)'at coo 1 ',atomcoo(ihr_at1)
c         write(*,*)'at coo 2 ',atomcoo(ihr_at2)

         do icoo = 1 , ncoord
            if (intcoori(icoo).eq.atomcoo(ihr_at1).or.
     $          intcoori(icoo).eq.atomcoo(ihr_at2)) goto 98
            index=index+1
            xintt(index)=xinti(icoo)
            intcoort(index)=intcoori(icoo)
 98        continue
         enddo

cc rewrite vectors

         do iatom = 1 , ncoord
            intcoori(iatom)=intcoort(iatom)
            xinti(iatom)=xintt(iatom)
         enddo

      else if (ispecies.eq.101)then
         index=0
         ncoord = 3*natom-6

         xintt(ncoord)=disthr_at1
         intcoort(ncoord)=atomcoo(ihr_at1)

c         write(*,*)'at index 1',ihr_at1
c         write(*,*)'at index 2',ihr_at2
c         write(*,*)'at coo 1 ',atomcoo(ihr_at1)
c         write(*,*)'at coo 2 ',atomcoo(ihr_at2)

         do icoo = 1 , ncoord
            if (intcoori(icoo).eq.atomcoo(ihr_at1)) goto 88
            index=index+1
            xintt(index)=xinti(icoo)
            intcoort(index)=intcoori(icoo)
 88        continue
         enddo

cc rewrite vectors

         do iatom = 1 , ncoord
            intcoori(iatom)=intcoort(iatom)
            xinti(iatom)=xintt(iatom)
         enddo

c         do iatom = 1 , ncoord
c            write(*,*)'int coo xint ',intcoori(iatom),xinti(iatom)
c         enddo

c         stop

      endif


cc start hindered rotor scan: 
      comsave1=comline1
      comsave2=comline2

cc if recovering data skip PES re-evaluation
      if(irecov.eq.1)   goto 999

      do ihind = 1 , nhind

cc for each rotor start from first level theory

      comline1=comsave1
      comline2=comsave2

cc start rotational scan for hind rotor ihind

cc    re-order z-matrix, taking out torsion ihind      

         write (16,*) 're-ordered z-matrix of hind rotor'

         index=0
         ncoord = 3*natom-6

cc     first assign optimized hind parameter to last coordinate
         do iatom = 1 , ncoord
            if (intcoori(iatom).eq.hindlab(ihind)) then
               xint(ncoord)=xinti(iatom)
            endif
         enddo

cc then update all coordinates

         do iatom = 1 , ncoord
            if (intcoori(iatom).eq.hindlab(ihind)) goto 998
            index=index+1
            xint(index)=xinti(iatom)
            intcoor(index)=intcoori(iatom)
            write (16,*) 'intcoor, xint ',intcoor(index),xint(index)
 998        continue
         enddo
         intcoor(ncoord)=hindlab(ihind)

c        write (6,*) 'starting g09pot'
c
cc save coordinates

         do iint = 1 , ncoord
            xint_save(iint) = xint(iint)
         enddo
cc start scan on nhind scan points

         write (16,*) 'hindered_rotor ',ihind

         ihind_step = (dhindmx(ihind)-dhindmn(ihind))/nhindsteps(ihind)


         ircons=1
         irepeat=0
         xstart=xint(ncoord)
         ihalf=nhindsteps(ihind)/2
         iprog=0
         if(ifrozrts.eq.1.and.ispecies.eq.0) ircons=2

 111     continue

         do iscanhind=1,nhindsteps(ihind)

            if(iscanhind.le.ihalf+1)iprog=iscanhind
            if(iscanhind.gt.ihalf+1)iprog=nhindsteps(ihind)
     $        +ihalf-iscanhind+2

            if(ispecies.eq.100) ircons=3
            if(ispecies.eq.101) ircons=2
            ntau=0
            ismp=iscanhind
            ixyz=0
            ired=0
            ifreq=0
            ilev=10
            if(iabs.eq.1) ireact=natom1+1
            if(iadd.eq.1) ireact=natom1

            if(ilev1code.eq.1) then
               call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $ coord,vtot_0,vtot,freq,ifreq,ilin,ismp,comline1,
     $ comline2,icharge,ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ired)

            else if (ilev1code.eq.2) then
               numproc=numprochl
         
               call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $ coord,vtot,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ilev
     $ ,ispecies)
               numproc=numprocll
            endif

c            write(*,*)'before g09'

c      do ih=1,nhind
c      do igroup=1,igrouptot(ih)
c      write (*,*)'hind ind atom',ih,igroup,ngroup(ih,igroup)
c      enddo
c      enddo
c            call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
c     $   coord,vtot_0,vtot,freq,ifreq,ilin,ismp,gkeyword,igkey,comline1,
c     $    comline2,icharge,ispin,ircons,
c     $    atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ired)

c            write(*,*)'after g09'
c      do ih=1,nhind
c      do igr=1,igrouptot(ih)
c      write (*,*)'hind ind atom',ih,igr,ngroup(ih,igr)
c      enddo
c      enddo


            write (16,*) xint(ncoord),vtot

            if(vtot.gt.1.0e10) then

               if(ispecies.eq.0.and.irepeat.eq.0)then

c                  if(irepeat.eq.0)then

cc switch method to determine HR potential
                  comline1=comline3
                  comline2=comline4
                  
                  ncoord = 3*natom-6
c                 write (16,*) 'int-2s ',atomcoo(isite)
c                 write (16,*) 'int-1s ',atomcoo(natomt1+1)
                  do icoo = 1 , ncoord
                     if (intcoor(icoo).eq.atomcoo(isite))then
                        xintt(ncoord-2)=xint_save(icoo)
                        intcoort(ncoord-2)=intcoor(icoo)
                        if(ibeta.eq.1)then
                           xintt(ncoord-1)=xint_save(icoo)
                           intcoort(ncoord-1)=intcoor(icoo)
                        endif
c                        write (16,*) 'int-2 ',intcoort(ncoord-2)
                     endif
                     if (intcoor(icoo).eq.'RTS')then
                        xintt(ncoord-1)=xint_save(icoo)
                        intcoort(ncoord-1)=intcoor(icoo)
c                        write (16,*) 'int-1 ',intcoort(ncoord-1)
                     endif
                  enddo

                  index=0

                  do icoo = 1 , ncoord-1
                     if (intcoor(icoo).eq.intcoort(ncoord-2).or.
     $                    intcoor(icoo).eq.intcoort(ncoord-1))goto 198
                     index=index+1
c                     write(16,*)'index intcoort(index)',index,
c     $                           intcoort(index)
                     xintt(index)=xint_save(icoo)
                     intcoort(index)=intcoor(icoo)
 198                 continue
                  enddo

cc rewrite vectors
                  do icoo = 1 , ncoord-1
                     intcoor(icoo)=intcoort(icoo)
                     xint(icoo)=xintt(icoo)
                     xint_save(icoo)=xintt(icoo)
                 write (16,*) 'intcoor, xint ',intcoor(icoo),xint(icoo)
                  enddo
                  xint(ncoord)=xstart
                  irepeat=1
                  if(iabs.eq.1) then
                     ircons=3
                  else if (iadd.eq.1) then
                     ircons=2
                  else if (iiso.eq.1) then
                     ircons=2
                  else if (ibeta.eq.1) then
                     ircons=2
                  endif
                  goto 111
               else
                  write(7,*)'this rotor did not converge'
                  write(7,*)'rotor number ',ihind
                  write(7,*)'rotor step ',iscanhind
                  write(7,*)'using the last calculated energy'
                  write(7,*)'using the coordinates of first point'
                  write(7,*)'for next point'
                  if(ispecies.ne.0)then
                     do iint = 1 , ncoord-1
                        xint(iint) = xint_save(iint)
                     enddo
                  endif
                  open (unit=65,status='unknown')
                  read (65,*) vtot
                  close(65)
                  write (16,*) 'Not converged, using last energy'
                  write (16,*) xint(ncoord),vtot
               endif

            endif

            vref(iprog)=vtot


            if (iscanhind.ne.1) then
               vref(iprog)=(vref(iprog)-vref(1))*CAUTOKCAL
            endif

            if(vref(iprog).lt.-0.1.and.iscanhind.ne.1) then
               write(7,*)'there is something wrong with this rotor'
               write(7,*)'found a negative energy'
               write(7,*)'rotor number ',ihind
               write(7,*)'rotor step ',iscanhind
               write(7,*)'check the me input'
            endif

            write (16,*) xint(ncoord),vtot
            
cc save hindered rotor optimized geometry 

c         command1='newzmat       
c     +    -ichk -oxyz tmp.chk geom.xyz           '
c         call commrun(command1)
            open (unit=98,file='headgeom.tmp',status='unknown')
            write(98,*)natom
            write(98,*)vtot
            close(98)
            command1='cat headgeom.tmp geom.xyz > geom1.xyz '
            call commrun(command1)

            open (unit=99,status='unknown')
            rewind (99)
            write(99,1200)ispecies,ihind,iprog
            rewind (99)
            read (99,2001) command1
            close (99)
            call commrun(command1)

 1200       format (" cp -f geom1.xyz 
     $          hr_geoms/geom_isp"I0.2"_hr"I0.2"_hpt"I0.2".xyz  ")

cc if TS save also z-mat and coords

            if(ispecies.eq.0) then
               open (unit=98,file='geom_temp.out',status='unknown')
               write(98,*)'TS_z-matrix from HR ',ihind,' point ',iprog
            else 
               open (unit=98,file='geom_temp.out',status='unknown')
               write(98,*)'z-matrix from HR ',ihind,' point ',iprog
            endif
            do iatom = 1 , natomt
               write (98,*) atomlabel(iatom)
            enddo
            do icoord = 1 , ncoord
               write (98,*) intcoor(icoord),xint(icoord)
            enddo
            close(98)

            open (unit=99,status='unknown')
            write(99,1201)ispecies,ihind,iprog
            rewind (99)
            read (99,2001) command1
            close (99)
            call commrun(command1)

 1201       format (" cp -f geom_temp.out 
     $ hr_geoms/geom_isp"I0.2"_hr"I0.2"_hpt"I0.2".out  ")

            if(iscanhind.le.ihalf) then
               xint(ncoord)= xstart+ihind_step*iscanhind
            else if (iscanhind.gt.ihalf)then
               xint(ncoord)= xstart-ihind_step*(iscanhind-ihalf)
            endif

cc get back to starting point for reverse sweep

            if(iscanhind.eq.ihalf+1) then
               do iint = 1 , ncoord-1
                  xint(iint) = xint_save(iint)
               enddo
            endif
            if (xint(ncoord).gt.360.)then
               xint(ncoord)= xint(ncoord)-360.
            endif

         enddo
cc end of hindered rotors PES scan 

         vref(1)=0.
cc write results in ME format         
         write (19,*) 'Rotor                             Hindered'
         write (19,1112) (ngroup(ihind,igr),igr=1,
     $         igrouptot(ihind))
         write (19,*) 'Axis  ',ipivotA(ihind),ipivotB(ihind)
         write (19,*) 'Symmetry ',nhindsymm(ihind)
         write (19,*) 'Potential[kcal/mol] ',nhindsteps(ihind)
         write (19,1111) (vref(iscanhind),iscanhind=1,
     $        nhindsteps(ihind))
         write (19,*)'End'

cc build xyz file with all scanned structures, one for each rotor

         write(command1,2201)ispecies,ihind
         call commrun(command1)
         do j=2,nhindsteps(ihind) 
            write(command1,2202)ispecies,ihind,j
            call commrun(command1)
            command1='cp -f temp1.xyz temp.xyz'
            call commrun(command1)
         enddo
         write(command1,2203)
         call commrun(command1)
         write(command1,2204)ispecies,ihind
         call commrun(command1)

 2201 format (" cp -f hr_geoms/geom_isp"I0.2"_hr"I0.2"_hpt01.xyz
     $          temp.xyz  ")

 2202 format (" cat hr_geoms/geom_isp"I0.2"_hr"I0.2"_hpt"I0.2".xyz
     $          temp.xyz > temp1.xyz  ")

 2203 format (" sed -i '/^$/d' temp.xyz  ")

 2204 format (" cp -f temp.xyz hr_geoms/geom_isp"I0.2"_hr"I0.2"_all.xyz
     $        ")


      enddo

cc end of hindered rotors scan 

      close (unit=19,status='keep')

 999  continue

      close (unit=16,status='keep')



      if(ispecies.eq.100.or.ispecies.eq.101) goto 9999

cc now project fc matrix and re-eveluate frequencies
      open (unit=15,file='hind_rot_head.dat',status='unknown')
      write(15,*)'Act_energy(kcal/mol):       0. '
      write(15,*)'Initial_Temperature:        200'
      write(15,*)'Temperature_steps:          40'
      write(15,*)'Temperature_increment:      50'
      write(15,*)'Delta_Energy_rea:           0.'
      write(15,*)'Delta_Energy_pro:           0.'
      write(15,*)'Maxstep:                    1'
      write(15,*)'Npointsint:                 5 '
      write(15,*)'Maxtdev:                    0.5'
      write(15,*)'Rearrange(1=yes,0=no)       1'
      write(15,*)'SaddlePoint                 1'
      write(15,*)'ds(1=noexp,0=standard)      0'
      write(15,*)'isct_vtst(1=vtst_sct,0=sct) 1'
      write(15,*)'zerocurvature(1)            0'
      write(15,*)'reduced_mass                1.0'
      write(15,*)'minimum_frequency            50'
      write(15,*)'anim_freq(if_Maxstep=1)       2'
      write(15,*)'onlyrotors(0=yes,1=no)        0'
      close (unit=15,status='keep')

      open (unit=15,file='hrdata4proj.dat',status='unknown')
      write (15,*)'proj_rea_coo(0=yes(def),1=no) 1'
      write (15,*)'numrotors        ',nhind
      do ihind=1,nhind
         write (15,*)'pivotA           ',ipivotA(ihind)
         write (15,*)'pivotB           ',ipivotB(ihind)
         write (15,*)'atomsintopA      ',igrouptot(ihind)
         write (15,*)'topAatoms        ',(ngroup(ihind,igr),igr=1,
     $         igrouptot(ihind))
      enddo
      close (unit=15,status='keep')

      open (unit=99,status='unknown')
      rewind (99)
      if (ispecies.eq.1) then
         command1='cp -f hrdata4proj.dat
     $   output/hrdata4proj_reac1.dat '
         call commrun(command1)
         write (99,2000)'extract_data_projrot.sh 
     $   ./output/reac1_fcmat.log ',natom
      endif
      if (ispecies.eq.2) then
         command1='cp -f hrdata4proj.dat
     $   output/hrdata4proj_reac2.dat '
         call commrun(command1)
         write (99,2000)'extract_data_projrot.sh 
     $   ./output/reac2_fcmat.log ',natom
      endif
      if (ispecies.eq.3) then
         command1='cp -f hrdata4proj.dat
     $   output/hrdata4proj_prod1.dat '
         call commrun(command1)
         write (99,2000)'extract_data_projrot.sh 
     $   ./output/prod1_fcmat.log ',natom
      endif
      if (ispecies.eq.4) then
         command1='cp -f hrdata4proj.dat
     $   output/hrdata4proj_prod2.dat '
         call commrun(command1)
         write (99,2000)'extract_data_projrot.sh 
     $   ./output/prod2_fcmat.log ',natom
      endif
      if (ispecies.eq.5.or.ispecies.eq.51) then
         command1='cp -f hrdata4proj.dat
     $   output/hrdata4proj_wellr.dat '
         call commrun(command1)
         write (99,2000)'extract_data_projrot.sh 
     $   ./output/wellr_fcmat.log ',natom
      endif
      if (ispecies.eq.6.or.ispecies.eq.61) then
         command1='cp -f hrdata4proj.dat
     $   output/hrdata4proj_wellp.dat '
         call commrun(command1)
         write (99,2000)'extract_data_projrot.sh 
     $   ./output/wellp_fcmat.log ',natom
      endif
      if (ispecies.eq.11) then
         command1='cp -f hrdata4proj.dat
     $   output/hrdata4proj_reacs.dat '
         call commrun(command1)
         write (99,2000)'extract_data_projrot.sh 
     $   ./output/reacs_fcmat.log ',natom
c    $   ./output/reac1_fcmat.log ',natom
      endif
      if (ispecies.eq.12) then
         command1='cp -f hrdata4proj.dat
     $   output/hrdata4proj_prods.dat '
         call commrun(command1)
         write (99,2000)'extract_data_projrot.sh 
     $   ./output/prods_fcmat.log ',natom
c    $   ./output/reac1_fcmat.log ',natom  ! changed by sjk because it looks wrong
      endif
      if (ispecies.eq.0) then
         command1='cp -f hrdata4proj.dat
     $   output/hrdata4proj_ts.dat '
         call commrun(command1)
         write (99,2000)'extract_data_projrot.sh 
     $   ./output/ts_fcmat.log ',natom
      endif
      rewind (99)
      read (99,2001) command1
      close (99)
      call commrun(command1)
      command1='RPHt.exe'
      call commrun(command1)

      if (ispecies.eq.1) then
         command1='cp -f ./me_files/reac1_hr.me 
     $                   ./me_files/reac1_1dhr.me'
         call commrun(command1)         
         command1='cp -f hrproj_freq.dat ./output/hrproj_freq_reac1.out'
         open (unit=19,file='./me_files/reac1_fr.me',status='unknown')
         open (unit=29,file='./me_files/reac1_nofr.me',status='unknown')
      endif

      if (ispecies.eq.2) then
         command1='cp -f ./me_files/reac2_hr.me 
     $                   ./me_files/reac2_1dhr.me'
         call commrun(command1)         
         command1='cp -f hrproj_freq.dat ./output/hrproj_freq_reac2.out'
         open (unit=19,file='./me_files/reac2_fr.me',status='unknown')
         open (unit=29,file='./me_files/reac2_nofr.me',status='unknown')
      endif

      if (ispecies.eq.3) then
         command1='cp -f ./me_files/prod1_hr.me 
     $                   ./me_files/prod1_1dhr.me'
         call commrun(command1)         
         command1='cp -f hrproj_freq.dat ./output/hrproj_freq_prod1.out'
         open (unit=19,file='./me_files/prod1_fr.me',status='unknown')
         open (unit=29,file='./me_files/prod1_nofr.me',status='unknown')
      endif

      if (ispecies.eq.4) then
         command1='cp -f ./me_files/prod2_hr.me 
     $                   ./me_files/prod2_1dhr.me'
         call commrun(command1)         
         command1='cp -f hrproj_freq.dat ./output/hrproj_freq_prod2.out'
         open (unit=19,file='./me_files/prod2_fr.me',status='unknown')
         open (unit=29,file='./me_files/prod2_nofr.me',status='unknown')
      endif

      if (ispecies.eq.5.or.ispecies.eq.51) then
         command1='cp -f ./me_files/wellr_hr.me 
     $                   ./me_files/wellr_1dhr.me'
         call commrun(command1)         
         command1='cp -f hrproj_freq.dat ./output/hrproj_freq_wellr.out'
         open (unit=19,file='./me_files/wellr_fr.me',status='unknown')
         open (unit=29,file='./me_files/wellr_nofr.me',status='unknown')
      endif

      if (ispecies.eq.6.or.ispecies.eq.61) then
         command1='cp -f ./me_files/wellp_hr.me 
     $                   ./me_files/wellp_1dhr.me'
         call commrun(command1)         
         command1='cp -f hrproj_freq.dat ./output/hrproj_freq_wellp.out'
         open (unit=19,file='./me_files/wellp_fr.me',status='unknown')
         open (unit=29,file='./me_files/wellp_nofr.me',status='unknown')
      endif

      if (ispecies.eq.0) then
         command1='cp -f ./me_files/ts_hr.me 
     $                   ./me_files/ts_1dhr.me'
         call commrun(command1)         
         command1='cp -f hrproj_freq.dat ./output/hrproj_freq_ts.out'
         open (unit=19,file='./me_files/ts_fr.me',status='unknown')
         open (unit=29,file='./me_files/ts_nofr.me',status='unknown')
      endif
      call commrun(command1)         

      nfreq=3*natom-6-ntau_fr
      if (natom.eq.2) nfreq=1
      if (ispecies.eq.0) nfreq=3*natom-6-ntau_fr-1

      open (unit=20,file='hrproj_freq.dat',status='old')
      if (ispecies.ne.0) then
         do j=1,nfreq
            read (20,*) freq(j)
         enddo
      else
         do j=1,nfreq+6+ntau_fr+1
            read (20,*) freq(j)
         enddo
      endif
      close (20)

      write (19,*) '    Frequencies[1/cm] ',nfreq
      write (19,8010) (freq(j),j=1,nfreq)
 8010 format (1x,10G12.5)

      if (ispecies.ne.0) then
         if ((ispecies.eq.1).or.(ispecies.eq.2).or.(ispecies.eq.3).or.
     $      (ispecies.eq.4)) then 
            write (19,*) ' ZeroEnergy[kcal/mol]             0.'
            write (29,*) ' ZeroEnergy[kcal/mol]             0.'
         endif
         if ((ispecies.eq.5).or.(ispecies.eq.51)) then
            write (19,*) ' ZeroEnergy[kcal/mol]             $wellren'
            write (29,*) ' ZeroEnergy[kcal/mol]             $wellren'
         endif
         if ((ispecies.eq.6).or.(ispecies.eq.61)) then
            write (19,*) ' ZeroEnergy[kcal/mol]             $wellpen'
            write (29,*) ' ZeroEnergy[kcal/mol]             $wellpen'
         endif
         write (19,*) ' ElectronicLevels[1/cm]           ',nelec
         write (29,*) ' ElectronicLevels[1/cm]           ',nelec
         do ielec = 1, nelec
            write (19,*) eelec(ielec),gelec(ielec)
            write (29,*) eelec(ielec),gelec(ielec)
         enddo
         write (19,*) 'End '
         write (29,*) 'End '
         if (ispecies.eq.1.and.iiso.eq.1) then
            write (19,*) 'End '
            write (29,*) 'End '
         endif
         if (ispecies.eq.1.and.ibeta.eq.1) then
            write (19,*) 'End '
            write (29,*) 'End '
         endif
         if (ispecies.eq.1.and.ireactonly.eq.1) then
            write (19,*) 'End '
            write (29,*) 'End '
         endif
         if (ispecies.eq.2) then
            write (19,*) 'GroundEnergy[kcal/mol] 0.0'
            write (19,*) 'End'
            write (29,*) 'GroundEnergy[kcal/mol] 0.0'
            write (29,*) 'End'
         endif
         if (ispecies.eq.3.and.ipr1.eq.1) then
            write (19,*) 'End '
            write (29,*) 'End '
         endif
         if (ispecies.eq.4) then
            write (19,*) 'GroundEnergy[kcal/mol] $proden'
            write (19,*) 'End'
            write (29,*) 'GroundEnergy[kcal/mol] $proden'
            write (29,*) 'End'
         endif
         if ((ispecies.eq.5).or.(ispecies.eq.51)) then
            write (19,*) 'End'
            write (29,*) 'End'
         endif
         if ((ispecies.eq.6).or.(ispecies.eq.61)) then
            write (19,*) 'End'
            write (29,*) 'End'
         endif
         write (19,*) '!************************************'
         write (29,*) '!************************************'
      else
         write (19,*) ' ZeroEnergy[kcal/mol]            $tsen '
         write (29,*) ' ZeroEnergy[kcal/mol]            $tsen '
         write (19,*) ' ElectronicLevels[1/cm]           ',nelec
         write (29,*) ' ElectronicLevels[1/cm]           ',nelec
         do ielec = 1, nelec
            write (19,*) eelec(ielec),gelec(ielec)
            write (29,*) eelec(ielec),gelec(ielec)
         enddo
         if(ifrozrts.ne.1)then
            if(imdtunn.ne.1)then
               write (19,*)'      Tunneling    Eckart'
               write (29,*)'      Tunneling    Eckart'
               write (19,*)'   ImaginaryFrequency[1/cm]   ',
     $              abs(freq(nfreq+6+ntau_fr+1))
               write (29,*)'   ImaginaryFrequency[1/cm]   ',
     $              abs(freq(nfreq+6+ntau_fr+1))
               write (19,*)'   WellDepth[kcal/mol]       $wdepfor'
               write (29,*)'   WellDepth[kcal/mol]       $wdepfor'
               write (19,*)'   WellDepth[kcal/mol]       $wdepback'
               write (29,*)'   WellDepth[kcal/mol]       $wdepback'
            else
               write (19,*)'      Tunneling    Read'
               write (29,*)'      Tunneling    Read'
               write (19,*)'   CutoffEnergy[1/cm]       2500'
               write (29,*)'   CutoffEnergy[1/cm]       2500'
               write (19,*)'   ImaginaryFrequency[1/cm]   ',
     $              abs(freq(nfreq+6+ntau_fr+1))
               write (29,*)'   ImaginaryFrequency[1/cm]   ',
     $              abs(freq(nfreq+6+ntau_fr+1))
               write (19,*)'   File imactint.dat'
               write (29,*)'   File imactint.dat'
            endif
         endif
         write (19,*)'End '
         write (19,*)'End '
         write (29,*)'End '
         write (29,*)'End '
c         write (19,*)'End '
c         write (19,*)'End '
         write (19,*)'!************************************'
         write (29,*)'!************************************'
      endif

      close (unit=19,status='keep')
      close (unit=29,status='keep')

 9999 continue

 1111    format (100(1X,f7.2))
 1112    format (' Group ',100(1X,I3))
 2000    format (A80,1X,I10)
 2001    format (A160)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mdtau(ispecies)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

cc    this is the multi dimensional hindered rotor section

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'


      dimension dhindmn(nhindmx),dhindmx(nhindmx),
     $ freq1d(noptmdmx,nmdmx),
     $ freq2d(noptmdmx,noptmdmx,nmdmx),nhindsteps(nhindmx),
     $ nhindsymm(nhindmx),freq3d(noptmdmx,noptmdmx,noptmdmx,nmdmx)

      dimension dhindmn1Da(nhindmx),dhindmx1Da(nhindmx),
     $ nhindsteps1Da(nhindmx),nhindsymm1Da(nhindmx)

      dimension dhindmn2Da(nhindmx),dhindmx2Da(nhindmx),
     $ nhindsteps2Da(nhindmx),nhindsymm2Da(nhindmx)
      dimension dhindmn2Db(nhindmx),dhindmx2Db(nhindmx),
     $ nhindsteps2Db(nhindmx),nhindsymm2Db(nhindmx)

      dimension dhindmn3Da(nhindmx),dhindmx3Da(nhindmx),
     $ nhindsteps3Da(nhindmx),nhindsymm3Da(nhindmx)
      dimension dhindmn3Db(nhindmx),dhindmx3Db(nhindmx),
     $ nhindsteps3Db(nhindmx),nhindsymm3Db(nhindmx)
      dimension dhindmn3Dc(nhindmx),dhindmx3Dc(nhindmx),
     $ nhindsteps3Dc(nhindmx),nhindsymm3Dc(nhindmx)

      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim),xintt(3*natommx),xint_save(3*natommx),
     $ xint_int(3*natommx),xint_int2(3*natommx)
      dimension vref1D(noptmx)
      dimension vref2D(noptmx,noptmx),tauo(ntaumx,noptmx),
     $ xinto(3*natommx,noptmx),dist(noptmx)
      dimension vref3D(noptmx,noptmx,noptmx)
      dimension vpot(noptmx)
      dimension tauopt(ntaumx)
      dimension tau(ntaumx)
      dimension iatomtype(natommx)
      dimension itop_ind(natommx,nhindmx),
     $          itop_ind_nodum(natommx,nhindmx)
      dimension idummy(natommx),
     $          ngroup(nhindmx,natommx)
      dimension igrouptot(nhindmx),ipivotA(nhindmx),ipivotB(nhindmx)
      dimension i1DtopA(nhindmx)
      dimension i2DtopA(nhindmx),i2DtopB(nhindmx)
      dimension i3DtopA(nhindmx),i3DtopB(nhindmx),i3DtopC(nhindmx)
      dimension taumn(ntaumx),taumx(ntaumx)
      dimension gelec(nelecmx),eelec(nelecmx)

      character*70 comline1,comline2,comline3,comline4,comsave1
      character*70 comline5,comline6,comlineref1
      character*70 comlineref2,comsave2
      character*60 atomlabel(natommx),atomlabel_save(natommx)
      character*4 atomname(natommx)
      character*4 atomconn(natommx)
      character*4 atomcoo(natommx)
      character*4 pivotA(nhindmx),pivotB(nhindmx)
      character*4 dicheck1(nhindmx),dicheck2(nhindmx)
      character*30 intcoor(3*natommx)
      character*30 intcoori(3*natommx)
      character*30 intcoort(3*natommx)
      character*30 intcoor_int(3*natommx)
      character*30 intcoor_int2(3*natommx)
      character*20 bislab(ntaumx)
      character*20 hindlab(nhindmx)
      character*20 hindlab1Da(nhindmx)
      character*20 hindlab2Da(nhindmx)
      character*20 hindlab2Db(nhindmx)
      character*20 hindlab3Da(nhindmx)
      character*20 hindlab3Db(nhindmx)
      character*20 hindlab3Dc(nhindmx)
      character*20 wordlab
      character*35 namefile,namefile2
      character*1 atfirstlet
      character*140 command1
      character*140 commandcopy
      character*10 word_1dhrout
      character*10 word_2dhrout
      character*10 word_3dhrout
      character*160 word_hr
      character*30 gmem
      character*100 buffer

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

      call Lineread(0)

c if not specified otherwise below, we use 1d active space also for multiD scan
      if(irecov.eq.1)then
         write(*,*) 'recovery option not implemented for mdtau scans'
         write(*,*) 'please modify the input'
      endif

      write(*,*)'entered mdtau'

      commandcopy='cp -f ./data/onedtau_molpro.dat 
     $       mdtau_molpro.dat'

      open (unit=21,file='./data/theory.dat',status='unknown')

      if (ispecies.eq.1) then
         open (unit=15,file='./data/reac1.dat',status='old')
         open (unit=16,file='./output/reac1_mDhr.out',status='unknown')
         open (unit=17,file='./output/reac1_opt.out',status='unknown')
         open (unit=19,file='./me_files/reac1_1dhr.me',
     $            status='unknown')
         open (unit=20,file='./me_files/reac1_mdhr.me',
     $            status='unknown')
         open (unit=26,file='./me_files/reac1_mdhr_nofr.me',
     $        status='unknown')
         word_1dhrout='reac1_1dhr'
         word_2dhrout='reac1_2dhr'
         word_3dhrout='reac1_3dhr'
c         open (unit=119,file='./me_files/rpst_hr1.me',status='unknown')
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_1'
         do while (WORD.NE.'HIND_ROTOR_1')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/mdtau_reac1_molpro.dat 
     $ onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.2) then
         open (unit=15,file='./data/reac2.dat',status='unknown')
         open (unit=16,file='./output/reac2_mDhr.out',status='unknown')
         open (unit=17,file='./output/reac2_opt.out',status='unknown')
         open (unit=19,file='./me_files/reac2_1dhr.me',
     $ status='unknown')
         open (unit=20,file='./me_files/reac2_mdhr.me',
     $ status='unknown')
         open (unit=26,file='./me_files/reac2_mdhr_nofr.me',
     $        status='unknown')
         word_1dhrout='reac2_1dhr'
         word_2dhrout='reac2_2dhr'
         word_3dhrout='reac2_3dhr'
c         open (unit=119,file='./me_files/rpst_hr2.me',status='unknown')
         inp_type=1
c        return
cs I removed this return - it would be fine for H,O,OH,O2,CH3,HO2 but
cs will create problems for larger radicals
cs I don't know if this removal will create other problems
         rewind(21)
         WORD_HR='HIND_ROTOR_2'
         do while (WORD.NE.'HIND_ROTOR_2')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/mdtau_reac2_molpro.dat 
     $ onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.3) then
         open (unit=15,file='./data/prod1.dat',status='old')
         open (unit=16,file='./output/prod1_mDhr.out',status='unknown')
         open (unit=17,file='./output/prod1_opt.out',status='unknown')
         open (unit=19,file='./me_files/prod1_1dhr.me',
     $ status='unknown')
         open (unit=20,file='./me_files/prod1_mdhr.me',
     $ status='unknown')
         open (unit=26,file='./me_files/prod1_mdhr_nofr.me',
     $        status='unknown')
         word_1dhrout='prod1_1dhr'
         word_2dhrout='prod1_2dhr'
         word_3dhrout='prod1_3dhr'
c         open (unit=119,file='./me_files/ppst_hr1.me',status='unknown')
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_3'
         do while (WORD.NE.'HIND_ROTOR_3')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/mdtau_prod1_molpro.dat 
     $ onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.4) then
         open (unit=15,file='./data/prod2.dat',status='old')
         open (unit=16,file='./output/prod2_mDhr.out',status='unknown')
         open (unit=17,file='./output/prod2_opt.out',status='unknown')
         open (unit=19,file='./me_files/prod2_1dhr.me',
     $ status='unknown')
         open (unit=20,file='./me_files/prod2_mdhr.me',
     $ status='unknown')
         open (unit=26,file='./me_files/prod2_mdhr_nofr.me',
     $        status='unknown')
         word_1dhrout='prod2_1dhr'
         word_2dhrout='prod2_2dhr'
         word_3dhrout='prod2_3dhr'
c         open (unit=119,file='./me_files/ppst_hr2.me',status='unknown')
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_4'
         do while (WORD.NE.'HIND_ROTOR_4')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/mdtau_prod2_molpro.dat 
     $ onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.5) then
         open (unit=15,file='./data/wellr.dat',status='old')
         open (unit=16,file='./output/wellr_mDhr.out',status='unknown')
         open (unit=17,file='./output/wellr_opt.out',status='unknown')
         open (unit=19,file='./me_files/wellr_1dhr.me',
     $ status='unknown')
         open (unit=20,file='./me_files/wellr_mdhr.me',
     $ status='unknown')
         open (unit=26,file='./me_files/wellr_mdhr_nofr.me',
     $        status='unknown')
         word_1dhrout='wellr_1dhr'
         word_2dhrout='wellr_2dhr'
         word_3dhrout='wellr_3dhr'
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_5'
         do while (WORD.NE.'HIND_ROTOR_5')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/mdtau_wellr_molpro.dat 
     $ onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.51) then
         open (unit=15,file='./data/wellr.dat',status='old')
         open (unit=16,file='./output/wellr_mDhr.out',status='unknown')
         open (unit=17,file='./output/wellr_opt.out',status='unknown')
         open (unit=19,file='./me_files/wellr_1dhr.me',
     $ status='unknown')
         open (unit=20,file='./me_files/wellr_mdhr.me',
     $ status='unknown')
         open (unit=26,file='./me_files/wellr_mdhr_nofr.me',
     $        status='unknown')
         word_1dhrout='wellr_1dhr'
         word_2dhrout='wellr_2dhr'
         word_3dhrout='wellr_3dhr'
         inp_type=2
         rewind(21)
         WORD_HR='HIND_ROTOR_51'
         do while (WORD.NE.'HIND_ROTOR_51')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/mdtau_wellr_molpro.dat 
     $ onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.6) then
         open (unit=15,file='./data/wellp.dat',status='old')
         open (unit=16,file='./output/wellp_mDhr.out',status='unknown')
         open (unit=17,file='./output/wellp_opt.out',status='unknown')
         open (unit=19,file='./me_files/wellp_1dhr.me',
     $ status='unknown')
         open (unit=20,file='./me_files/wellp_mdhr.me',
     $ status='unknown')
         open (unit=26,file='./me_files/wellp_mdhr_nofr.me',
     $        status='unknown')
         word_1dhrout='wellp_1dhr'
         word_2dhrout='wellp_2dhr'
         word_3dhrout='wellp_3dhr'
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_6'
         do while (WORD.NE.'HIND_ROTOR_6')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/mdtau_wellp_molpro.dat 
     $ onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.61) then
         open (unit=15,file='./data/wellp.dat',status='old')
         open (unit=16,file='./output/wellp_mDhr.out',status='unknown')
         open (unit=17,file='./output/wellp_opt.out',status='unknown')
         open (unit=19,file='./me_files/wellp_1dhr.me',
     $ status='unknown')
         open (unit=20,file='./me_files/wellp_mdhr.me',
     $ status='unknown')
         open (unit=26,file='./me_files/wellp_mdhr_nofr.me',
     $        status='unknown')
         word_1dhrout='wellp_1dhr'
         word_2dhrout='wellp_2dhr'
         word_3dhrout='wellp_3dhr'
         inp_type=2
         rewind(21)
         WORD_HR='HIND_ROTOR_61'
         do while (WORD.NE.'HIND_ROTOR_61')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/mdtau_wellp_molpro.dat 
     $ onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.11) then
         open (unit=15,file='./data/reacs.dat',status='old')
         open (unit=16,file='./output/reacs_mDhr.out',status='unknown')
         open (unit=17,file='./output/reacs_opt.out',status='unknown')
         open (unit=19,file='./me_files/reacs_1dhr.me'
     $        ,status='unknown')
         open (unit=20,file='./me_files/reacs_mdhr.me',
     $ status='unknown')
         open (unit=26,file='./me_files/reacs_mdhr_nofr.me',
     $        status='unknown')
         word_1dhrout='reacs_1dhr'
         word_2dhrout='reacs_2dhr'
         word_3dhrout='reacs_3dhr'
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_11'
         do while (WORD.NE.'HIND_ROTOR_11')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/mdtau_reacs_molpro.dat 
     $ onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.12) then
         open (unit=15,file='./data/prods.dat',status='old')
         open (unit=16,file='./output/prods_mDhr.out',status='unknown')
         open (unit=17,file='./output/prods_opt.out',status='unknown')
         open (unit=19,file='./me_files/prods_1dhr.me',
     $ status='unknown')
         open (unit=20,file='./me_files/prods_mdhr.me',
     $ status='unknown')
         open (unit=26,file='./me_files/prods_mdhr_nofr.me',
     $        status='unknown')
         word_1dhrout='prods_1dhr'
         word_2dhrout='prods_2dhr'
         word_3dhrout='prods_3dhr'
         inp_type=1
         rewind(21)
         WORD_HR='HIND_ROTOR_12'
         do while (WORD.NE.'HIND_ROTOR_12')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/mdtau_prods_molpro.dat 
     $ onedtau_molpro.dat'
         endif
      endif
      if (ispecies.eq.0) then
         open (unit=15,file='./data/ts.dat',status='unknown')
         open (unit=16,file='./output/ts_mDhr.out',status='unknown')
         open (unit=17,file='./output/ts_opt.out',status='unknown')
         open (unit=19,file='./me_files/ts_1dhr.me',status='unknown')
         open (unit=20,file='./me_files/ts_mdhr.me',status='unknown')
         open (unit=26,file='./me_files/ts_mdhr_nofr.me',
     $        status='unknown')
         word_1dhrout='tsgta_1dhr'
         word_2dhrout='tsgta_2dhr'
         word_3dhrout='tsgta_3dhr'
         inp_type=2
         rewind(21)
         WORD_HR='HIND_ROTOR_TS'
         do while (WORD.NE.'HIND_ROTOR_TS')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_HR = 'HIND_ROTOR'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/mdtau_ts_molpro.dat 
     $ onedtau_molpro.dat'
         endif
      endif

 900  continue
      rewind(21)

      ifreq=0
      noptg = 0
      imhrfr=0
      ires=0
      nhind1d=0
      nhind2d=0
      nhind3d=0
cc we assume that a linear molecule has no rotational dof         
      ilin=0

      write(*,*)'ok1 mdtau'


cc initialize the potential vectors to 0.
      do i=1,noptmx
         vref1D(i)=0.
         do j=1,noptmx
            vref2D(i,j)=0.
            do k=1,noptmx
               vref3D(i,j,k)=0.
            enddo
         enddo
      enddo

cc read 1d hind rotor list
      do while (WORD.NE.'NHIND')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (16,*) 'hind rotors must be defined'
            stop
         endif
      enddo
c      call LineRead (15)
      read (15,*) nhind
      read (15,*)
      ntau_fr = nhind
      do ihind = 1 , nhind
         read (15,*) hindlab(ihind),dhindmn(ihind),dhindmx(ihind),
     +              nhindsteps(ihind),nhindsymm(ihind)
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) hindlab(ihind)
         rewind (99)
         call LineRead (99)
         hindlab(ihind)=WORD
         close (unit=99,status='keep')
      enddo
      if (nhind.eq.0) then
         close (unit=15,status='keep')
         close (unit=16,status='keep')
         close (unit=17,status='keep')
         close (unit=19,status='keep')
         word=' '
         return
      endif
      rewind(15)

c      write(*,*)'ok12 mdtau'

cc read 1d multihind rotor list
      mhind1=1
      do while (WORD.NE.'NHIND1D')
         call LineRead (15)
         if (WORD.EQ.'END') then
            mhind1=0
            goto 8999
         endif
      enddo
c      call LineRead (15)
      read (15,*) nhind1d
      read (15,*)
      do ihind = 1 , nhind1d
         read (15,*) hindlab1Da(ihind),dhindmn1Da(ihind),
     +        dhindmx1Da(ihind),nhindsteps1Da(ihind),nhindsymm1Da(ihind)
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) hindlab1Da(ihind)
         rewind (99)
         call LineRead (99)
         hindlab1Da(ihind)=WORD
         close (unit=99,status='keep')
      enddo
 8999  continue
      word=' '
      rewind(15)


c      write(*,*)'ok22 mdtau'
cc read 2d hind rotor list
      mhind2=1
      do while (WORD.NE.'NHIND2D')
         call LineRead (15)
         if (WORD.EQ.'END') then
            mhind2=0
            goto 9001
c            write (16,*) '2Dhind rotors must be defined'
c            stop
         endif
      enddo
c      call LineRead (15)
      read (15,*) nhind2d
      read (15,*)
      do ihind = 1 , nhind2d
         read (15,*) hindlab2Da(ihind),dhindmn2Da(ihind),
     +        dhindmx2Da(ihind),nhindsteps2Da(ihind),nhindsymm2Da(ihind)
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) hindlab2Da(ihind)
         rewind (99)
         call LineRead (99)
         hindlab2Da(ihind)=WORD
         close (unit=99,status='keep')

         read (15,*) hindlab2Db(ihind),dhindmn2Db(ihind),
     +        dhindmx2Db(ihind),nhindsteps2Db(ihind),nhindsymm2Db(ihind)
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) hindlab2Db(ihind)
         rewind (99)
         call LineRead (99)
         hindlab2Db(ihind)=WORD
         close (unit=99,status='keep')
      enddo
 9001 continue
      word=' '
      rewind(15)


cc read 3d hind rotor list
c      write(*,*)'ok112 mdtau'
      mhind3=1
      do while (WORD.NE.'NHIND3D')
         call LineRead (15)
c         write(16,*)'word is ',word
c         write(16,*)'mhind is ',mhind
         if (WORD.EQ.'END'.and.mhind1.eq.0.and.mhind2.eq.0) then
c            write(*,*)'from here'
            write (16,*) '1D, 2D or 3Dhind rotors must be defined'
            close(16)
            stop
         else if (WORD.EQ.'END') then
            mhind3=0
            goto 9002
         endif
      enddo
c      call LineRead (15)
      read (15,*) nhind3d
      read (15,*)
c      write(*,*)'ok32 mdtau'
c      stop
      do ihind = 1 , nhind3d
         read (15,*) hindlab3Da(ihind),dhindmn3Da(ihind),
     +        dhindmx3Da(ihind),nhindsteps3Da(ihind),nhindsymm3Da(ihind)
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) hindlab3Da(ihind)
         rewind (99)
         call LineRead (99)
         hindlab3Da(ihind)=WORD
c      write(*,*)'ok221 mdtau'
         close (unit=99,status='keep')
c
         read (15,*) hindlab3Db(ihind),dhindmn3Db(ihind),
     +        dhindmx3Db(ihind),nhindsteps3Db(ihind),nhindsymm3Db(ihind)
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) hindlab3Db(ihind)
         rewind (99)
         call LineRead (99)
         hindlab3Db(ihind)=WORD
         close (unit=99,status='keep')
c
         read (15,*) hindlab3Dc(ihind),dhindmn3Dc(ihind),
     +        dhindmx3Dc(ihind),nhindsteps3Dc(ihind),nhindsymm3Dc(ihind)
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) hindlab3Dc(ihind)
         rewind (99)
         call LineRead (99)
         hindlab3Dc(ihind)=WORD
         close (unit=99,status='keep')
      enddo
 9002 continue
      word=' '
      rewind(15)

      if (nhind1d.eq.0.and.nhind2d.eq.0.and.nhind3d.eq.0) then
         close (unit=15,status='keep')
         close (unit=16,status='keep')
         close (unit=17,status='keep')
         close (unit=19,status='keep')
         word=' '
         return
      endif
      nhind_check=nhind1d+nhind2d+nhind3d
      if(nhind_check.gt.1)then
         write (16,*) 'only one md rotor supported at present'
         write (16,*) 'the program will be stopped'
         close(16)
         stop
      endif
c      write(*,*)'ok2 mdtau'
      write(*,*)'nhind2d is ',nhind2d
      write(*,*)'nhind3d is ',nhind3d
c      stop

cc initialize vectors stating correspondence between 1D and 1Dmulti HRs
      do ihind = 1 , nhind1d
         i1DtopA(ihind)=0
      enddo
      do ihind = 1 , nhind1d
         do i = 1 , nhind
            if(hindlab1Da(ihind).eq.hindlab(i)) then
               i1DtopA(ihind)=i
            endif
         enddo
      enddo
      
      itestA=0
      iad=0
      do ihind = 1 , nhind1d
         if(i1DtopA(ihind).ne.0)iad=1
         itestA=itestA+iad
         iad=0
      enddo
      if(nhind1d.ne.0)then
         if(itestA.ne.1) then
            write (16,*) '1D multi rotors must be def. as 1D rotors'
            write (16,*) 'and 1D rotor calcs must preceed 1D mHRs'
            stop
         endif
      endif


cc initialize vectors stating correspondence between 1D and 2D HRs
      do ihind = 1 , nhind2d
         i2DtopA(ihind)=0
         i2DtopB(ihind)=0
      enddo
      do ihind = 1 , nhind2d
         do i = 1 , nhind
            if(hindlab2Da(ihind).eq.hindlab(i)) then
               i2DtopA(ihind)=i
            endif
            if(hindlab2Db(ihind).eq.hindlab(i)) then
               i2DtopB(ihind)=i
            endif
         enddo
      enddo
      
      itestA=0
      itestB=0
      iad=0
      do ihind = 1 , nhind2d
         if(i2DtopA(ihind).ne.0)iad=1
         itestA=itestA+iad
         iad=0
         if(i2DtopB(ihind).ne.0)iad=1
         itestB=itestB+iad
         iad=0
      enddo
      if(nhind2d.ne.0)then
         if(itestA.ne.1.or.itestB.ne.1) then
            write (16,*) '2D rotors must also be defined as 1D rotors'
            write (16,*) 'and 1D rotor calculations must preceed 2D HRs'
            stop
         endif
      endif

cc initialize vectors stating correspondence between 1D and 3D HRs
      do ihind = 1 , nhind3d
         i3DtopA(ihind)=0
         i3DtopB(ihind)=0
         i3DtopC(ihind)=0
      enddo
      do ihind = 1 , nhind3d
         do i = 1 , nhind
            if(hindlab3Da(ihind).eq.hindlab(i)) then
               i3DtopA(ihind)=i
            endif
            if(hindlab3Db(ihind).eq.hindlab(i)) then
               i3DtopB(ihind)=i
            endif
            if(hindlab3Dc(ihind).eq.hindlab(i)) then
               i3DtopC(ihind)=i
            endif
         enddo
      enddo
      
      itestA=0
      itestB=0
      itestC=0
      iad=0
      do ihind = 1 , nhind3d
         if(i3DtopA(ihind).ne.0)iad=1
         itestA=itestA+iad
         iad=0
         if(i3DtopB(ihind).ne.0)iad=1
         itestB=itestB+iad
         iad=0
         if(i3DtopC(ihind).ne.0)iad=1
         itestC=itestC+iad
         iad=0
      enddo
      if(nhind3d.ne.0)then
         if(itestA.ne.1.or.itestB.ne.1.or.itestC.ne.1) then
            write (16,*) '3D rotors must also be defined as 1D rotors'
            write (16,*) 'and 1D rotor calculations must preceed 3D HRs'
            stop
         endif
      endif

cc read input info for es calculations
      do while (WORD.NE.'NELEC')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (16,*) 'electronic states must be defined'
            stop
         endif
      enddo
      read (15,*) nelec
      do ielec = 1 , nelec
         read (15,*) eelec(ielec),gelec(ielec)
      enddo
      rewind(15)

c read input of type 1

      if (inp_type.eq.1) then

c read level of theory for hindered rotor scan

         do while (WORD.NE.'NATOM')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'natom must be defined'
               stop
            endif
         enddo
         read (15,*) natom,natomt
         rewind(15)

         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'sampling coordinates must be defined'
               stop
            endif
         enddo
         read (15,*) ntau
         rewind(15)


         do while (WORD.NE.'CHARGE')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'charge and spin must be defined'
               stop
            endif
         enddo
         read (15,*) icharge,ispin
         do iatom = 1 , natomt
            read (15,'(A60)') atomlabel(iatom)
            atomlabel_save(iatom)= atomlabel(iatom)
         enddo
         rewind(15)

         do while (WORD.NE.WORD_HR)
            call LineRead (21)
            if (WORD.EQ.'END') then
               write (16,*) 'hind rotors must be defined'
               write (16,*) 'in file theory.dat'
               stop
            endif
         enddo
         if(word2.eq.'G09') then
            ilev1code=1
            read (21,'(A70)') comline1
            read (21,'(A70)') comline2
         else if(word2.eq.'MOLPRO')then
            ilev1code=2
            call commrun(commandcopy)
         else
            write (16,*) 'multiD hind rot theory must be either'
            write (16,*) 'g09 or molpro in theory.dat'
            stop
         endif
         if(word3.eq.'RESTART') then
            open (unit=99,status='unknown')
            write (99,*) word4
            rewind (99)
            read (99,*) ires
            close (unit=99,status='keep')
         endif

         if (idebug.ge.2) write (16,*) ' comline test',comline1,comline2

cc now read theory lines for frequency calculations for multirotor analysis
         call LineRead (21)
         if(word.eq.'MHR_FREQS') then
            imhrfr=1
            if(ilev1code.eq.1)then
               call comline56_g09(ispecies,comline1,comline2,comline5,
     +                            comline6)
            endif
c            read (21,'(A70)') comline5
c            read (21,'(A70)') comline6
         endif
         close(21)

         if (natom.ne.1) then
            do while (WORD.NE.'INTCOOR')
               call LineRead (15)
               if (WORD.EQ.'END') then
                  write (16,*) 'internal coordinates must be defined'
                  stop
               endif
            enddo
            ncoord = 3*natom-6-ntau
            if (natom.eq.1) ncoord = 0
            if (natom.eq.2) ncoord = 1
            do icoord = 1 , ncoord
               call LineRead (15)
               intcoori(icoord) = word
            enddo
            rewind(15)
         endif

         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'sampling coordinates must be defined'
               stop
            endif
         enddo
         read (15,*) ntau
         if (ntau.gt.ntaumx) then
            stop
         endif
         if (ntau.ne.0) then
            read (15,*)
            do itau = 1 , ntau
               read (15,*) bislab(itau),taumn(itau),taumx(itau)
               write (16,*) bislab(itau),taumn(itau),taumx(itau)
               open (unit=99,status='unknown')
               rewind (99)
               write (99,*)bislab(itau)
               rewind (99)
               call LineRead (99)
               icoord=ncoord+itau
               intcoori(icoord)=WORD
               close (unit=99,status='keep')
            enddo
         endif
         rewind(15)


         ncoord = 3*natom-6
         read (17,*)
         do icoord = 1 , ncoord
            call LineRead (17)
            OPEN (unit=99,status='unknown')
            REWIND (99)
            WRITE (99,1000) WORD
 1000       FORMAT (A30)
            REWIND (99)
            READ (99,*) xinti(icoord)
            close (unit=99,status='keep')
         enddo
         
         if (idebug.ge.2) write (16,*) 'past z-matrix'
      endif


c now read input of type 2

      if (inp_type.eq.2) then

         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'sampling coords of reac1 must be defined'
               stop
            endif
         enddo
         read (15,*) ntau
         rewind 15

         open (unit=25,file='./data/reac1.dat',status='old')
         word=' '
         do while (WORD.NE.'NTAU')
            call LineRead (25)
            if (WORD.EQ.'END') then
               write (16,*) 'sampling coords of reac1 must be defined'
               stop
            endif
         enddo
         read (25,*) ntau1
         rewind 25

c     natomt is to account for dummy atoms
         do while (WORD.NE.'NATOM')
            call LineRead (25)
            if (WORD.EQ.'END') then
               write (16,*) 'natom in reac1 must be defined'
               stop
            endif
         enddo
         read (25,*) natom1,natomt1
         close (25)

cc get data from react2 file

         if(iabs.eq.1.or.iadd.eq.1)then
            open (unit=25,file='./data/reac2.dat',status='old')

            do while (WORD.NE.'NTAU')
               call LineRead (25)
               if (WORD.EQ.'END') then
                  write (66,*) 'define sampling coords of reac2'
                  stop
               endif
            enddo
            read (25,*) ntau2
            rewind 25

c     natomt is to account for dummy atoms
            do while (WORD.NE.'NATOM')
               call LineRead (25)
               if (WORD.EQ.'END') then
                  write (66,*) 'natom must be defined'
                  stop
               endif
            enddo
            read (25,*) natom2,natomt2
            close (25)
         else
            natom2=0
            natomt2=0
            ntau2=0
         endif
         natom = natom1+natom2
         if (iadd.eq.1) then
            natomt = natomt1+natomt2
         else if (iabs.eq.1) then
            natomt = natomt1+natomt2+1
         else
            natomt = natomt1+natomt2
         endif
c        natomt = natomt1+natomt2+1

        do while (WORD.NE.'CHARGE')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (16,*) 'charge and spin must be defined'
               stop
            endif
         enddo
         read (15,*) icharge,ispin
         close(15)
c         rewind(15)
         
         if(ispecies.eq.0) then
            do while (WORD.NE.WORD_HR)
               call LineRead (21)
               if (WORD.EQ.'END') then
                  write (16,*) 'hind rotors must be defined for TS'
                  write (16,*) 'in file theory.dat'
                  stop
               endif
            enddo
         else if((ispecies.eq.100).or.(ispecies.eq.51).or.
     $    (ispecies.eq.61)) then
            do while (WORD.NE.WORD_HR)
               call LineRead (21)
               if (WORD.EQ.'END') then
                  write (16,*) 'hind rotors must be defined'
                  write (16,*) 'in file theory.dat'
                  stop
               endif
            enddo
         endif
         if(word2.eq.'G09') then
            ilev1code=1
            read (21,'(A70)') comline1
            read (21,'(A70)') comline2
            if(ispecies.eq.0) then
               read (21,'(A70)') comline3
               read (21,'(A70)') comline4
            endif
            call LineRead (21)
            if(word.eq.'MHR_FREQS') then
c               read (21,'(A70)') comline5
c               read (21,'(A70)') comline6
               imhrfr=1
c               write(*,*)'imhrf is ',imhrfr
c               stop
               if(ilev1code.eq.1)then
                 call comline56_g09(ispecies,comline1,comline2,comline5,
     +                comline6)
               endif
            endif
         else if(word2.eq.'MOLPRO')then
            ilev1code=2
            call commrun(commandcopy)
         else
            write (16,*) 'multiD hind rot code must be either'
            write (16,*) 'g09 or molpro in theory.dat'
            stop
         endif
         if(word3.eq.'RESTART') then
            open (unit=99,status='unknown')
            write (99,*) word4
            rewind (99)
            read (99,*) ires
            close (unit=99,status='keep')
         endif
         if (idebug.ge.2) write (16,*) ' comline test',comline1,comline2
         if (idebug.ge.2) write (16,*) ' test2       ',comline3,comline4
         if (idebug.ge.2) write (16,*) ' test3       ',comline5,comline6
         close(21)

c now read z-mat input

         read (17,*)

         if (idebug.ge.2) write (6,*) ' starting gaussian input'
         do iatom = 1 , natomt
            read (17,'(A60)') atomlabel(iatom)
            atomlabel_save(iatom)= atomlabel(iatom)
c            write (6,*) atomlabel(iatom)
         enddo

cc read coordinate names

         ncoord = 3*natom-6

         do iint = 1 , ncoord
            read (17,*) intcoori(iint),xinti(iint)
            write (6,*) intcoori(iint),xinti(iint)
         enddo

         close (unit=17,status='keep')

         if (idebug.ge.2) write (16,*) 'past z-matrix'

cc now read abstraction site from grid input file

         open (unit=18,file='./data/ts.dat',status='old')

         do while (WORD.NE.'ISITE')
            call LineRead (18)
            if (WORD.EQ.'END') then
               write (16,*) 'abstraction site must be defined'
               stop
            endif
         enddo
         read (18,*) isite,jsite,ksite

         close (unit=18,status='keep')

      endif

cc initialize atomcoo vector

      OPEN (unit=99,status='unknown')
      REWIND (99)
      write (99,*) atomlabel(1)
      REWIND (99)
      read (99,*) atomname(1)
      close (99)
      atomconn(1)='0000'
      atomcoo(1)='0000'
      
      OPEN (unit=99,status='unknown')
      REWIND (99)
      write (99,*) atomlabel(2)
      REWIND (99)
      read (99,*) atomname(2),atomconn(2),atomcoo(2)
      close (99)

      OPEN (unit=99,status='unknown')
      REWIND (99)
      write (99,*) atomlabel(3)
      REWIND (99)
      read (99,*) atomname(3),atomconn(3),atomcoo(3)
      close (99)

      do iatom = 4 , natomt
         OPEN (unit=99,status='unknown')
         REWIND (99)
         write (99,*) atomlabel(iatom)
         REWIND (99)
         read (99,*) atomname(iatom),atomconn(iatom),atomcoo(iatom)
     $        ,word2,word,word,wordlab
         close (99)
      enddo
cc convert to upper cases
      do iatom = 1 , natomt
         OPEN (unit=99,status='unknown')
         REWIND (99)
         write (99,*) atomname(iatom),' ',atomconn(iatom),' ',
     $               atomcoo(iatom)
         REWIND (99)
         call LineRead(99)
         atomname(iatom)=word
         atomconn(iatom)=word2
         atomcoo(iatom)=word3
         close (99)
      enddo

cc 1D hindered rotor section
cc scan each 1D hindered rotor at a time

cc start hindered rotor scan: 


      comsave1=comline1
      comsave2=comline2
      comlineref1=comline1
      comlineref2=comline2

      write (16,*) '1D hindered_rotor ',ihind

      do ihind = 1 , nhind1d

cc for each 1Drotor start from first level theory

         comline1=comsave1
         comline2=comsave2

cc start rotational scan for hind rotor ihind

cc    re-order z-matrix, taking out torsion ihind      

         write (16,*) 're-ordered z-matrix of hind rotor'

         ncoord = 3*natom-6

         do iatom = 1 , ncoord
            xint_int(iatom)=0.
            intcoor_int(iatom)=''
         enddo

cc     first assign optimized hind parameter of the HR to last coordinate
         do iatom = 1 , ncoord
            if (intcoori(iatom).eq.hindlab1Da(ihind)) then
               xint(ncoord)=xinti(iatom)
            endif
         enddo

cc then update all coordinates for the rotor
         index=0
         do iatom = 1 , ncoord
            if (intcoori(iatom).eq.hindlab1Da(ihind)) goto 997
            index=index+1
            xint_int(index)=xinti(iatom)
            intcoor_int(index)=intcoori(iatom)
 997        continue
         enddo
         intcoor_int(ncoord)=hindlab1Da(ihind)

         do iatom = 1 , ncoord
            write (16,*) 'intc,xint ',intcoor_int(iatom),xint_int(iatom)
         enddo


cc         start first cycle of 1D rotors

         ihind_stepa = (dhindmx1Da(ihind)-dhindmn1Da(ihind))
     +      /nhindsteps1Da(ihind)

         xstarta=xint(ncoord)

         irepeat=0

         index=0
         do iatom = 1 , ncoord
            if (intcoor_int(iatom).eq.hindlab1Da(ihind)) goto 9198
            index=index+1
            xint(index)=xint_int(iatom)
            intcoor(index)=intcoor_int(iatom)
 9198       continue
         enddo
         intcoor(ncoord)=hindlab1Da(ihind)

cc save coordinates
         do iint = 1 , ncoord
            xint_save(iint) = xint(iint)
         enddo

 101     continue

         ircons=1
         refen=0.
         ihalf=nhindsteps(ihind)/2
         iprog=0

         do iscanhinda=1,nhindsteps1Da(ihind)

cc start scan on nhind scan points

            ntau=0
            ismp=iscanhinda
            ixyz=0
            ired=0
            ilev=10

            if(iscanhinda.le.ihalf+1)iprog=iscanhinda
            if(iscanhinda.gt.ihalf+1)iprog=nhindsteps1Da(ihind)
     $           +ihalf-iscanhinda+2



c               stop
c            write(*,*)'before g09'
c            write(*,*)'comline 1 ',comline1
c            write(*,*)'comline 2 ',comline2
c            stop

            if(ilev1code.eq.1) then

               call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $   coord,vtot_0,vtot,freq,ifreq,ilin,ismp,comline1,
     $    comline2,icharge,ispin,ircons,
     $    atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ired)

c            write(*,*)'after g09'
            else if (ilev1code.eq.2) then
               numproc=numprochl
               call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $ coord,vtot,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ilev
     $ ,ispecies)
               numproc=numprocll
            endif

            write (16,*) xint(ncoord),vtot
            
            if(vtot.gt.1.0e10) then
                  
               if(ispecies.eq.0.and.irepeat.eq.0)then

cc switch method to determine HR potential
                  comline1=comline3
                  comline2=comline4
                  
                  ncoord = 3*natom-6
 
                  do icoo = 1 , ncoord
                     if (intcoor(icoo).eq.atomcoo(isite))then
                        xintt(ncoord-2)=xint_save(icoo)
                        intcoort(ncoord-2)=intcoor(icoo)
                     endif
                     if (intcoor(icoo).eq.'RTS')then
                        xintt(ncoord-1)=xint_save(icoo)
                        intcoort(ncoord-1)=intcoor(icoo)
c     write (16,*) 'int-1 ',intcoort(ncoord-1)
                     endif
                  enddo

                  index=0

                  do icoo = 1 , ncoord-1
                     if (intcoor(icoo).eq.intcoort(ncoord-2).or.
     $                    intcoor(icoo).eq.intcoort(ncoord-1))goto 98
                     index=index+1
c     write(16,*)'index intcoort(index)',index,
c     $                           intcoort(index)
                     xintt(index)=xint_save(icoo)
                     intcoort(index)=intcoor(icoo)
 98                  continue
                  enddo

cc rewrite vectors
                  do icoo = 1 , ncoord-1
                     intcoor(icoo)=intcoort(icoo)
                     xint(icoo)=xintt(icoo)
                     xint_save(icoo)=xintt(icoo)
                     write (16,*) 'intcoor, xint ',intcoor(icoo)
     +                 ,xint(icoo)
                  enddo
                  xint(ncoord)=xstarta
                  irepeat=1
                  ircons=3
                  goto 101
               else
                  write(7,*)'this rotor did not converge'
                  write(7,*)'rotor number ',ihind
                  write(7,*)'rotor step ',iscanhinda
                  write(7,*)'using the last calculated energy'
                  write(7,*)'using the coordinates of first point'
                  write(7,*)'for next point'
                  if(ispecies.ne.0)then
                     do iint = 1 , ncoord-1
                        xint(iint) = xint_save(iint)
                     enddo
                  endif
                  open (unit=65,status='unknown')
                  read (65,*) vtot
                  close(65)
                  write (16,*) 'Not converged, using last energy'
                  write (16,*) xint(ncoord),vtot
               endif

            endif
         
            vref1D(iprog)=vtot

            if (iscanhinda.eq.1) then
               refen=vref1D(iscanhinda)
            else
               vref1D(iprog)=
     +              (vref1D(iprog)
     +              -refen)*CAUTOKCAL
            endif

            if(vref1D(iprog).lt.-0.1.
     +           and.iscanhinda.ne.1)then
               write(7,*)'there is something wrong with this rotor'
               write(7,*)'found a negative energy'
               write(7,*)'1Drotor number ',ihind
               write(7,*)'rotor step a',iscanhinda
               write(7,*)'rotor step b',iscanhindb
               write(7,*)'check the me input'
            endif

            if(imhrfr.eq.1) then
               ixyz=0
               ired=0

               if(ilev1code.eq.1) then
                  comline1=comline5
                  comline2=comline6
                  atomlabel(1)=' '
                  atomlabel(2)=' '

                  call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $   coord,vtot_0,vtot,freq,ifreq,ilin,ismp,comline1,
     $   comline2,icharge,ispin,ircons,
     $    atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ired)

                  atomlabel(1)= atomlabel_save(1)
                  atomlabel(2)= atomlabel_save(2)

               else if (ilev1code.eq.2) then
                  command1='cp -f fcmat.log ./geom.log '
                  call commrun(command1)
               endif

cc now project fc matrix and re-eveluate frequencies
cc use RPHt file from 1D rotor analysis
               open (unit=15,file='hind_rot_head.dat',status='unknown')
               write(15,*)'Act_energy(kcal/mol):       0. '
               write(15,*)'Initial_Temperature:        200'
               write(15,*)'Temperature_steps:          40'
               write(15,*)'Temperature_increment:      50'
               write(15,*)'Delta_Energy_rea:           0.'
               write(15,*)'Delta_Energy_pro:           0.'
               write(15,*)'Maxstep:                    1'
               write(15,*)'Npointsint:                 5 '
               write(15,*)'Maxtdev:                    0.5'
               write(15,*)'Rearrange(1=yes,0=no)       1'
               write(15,*)'SaddlePoint                 1'
               write(15,*)'ds(1=noexp,0=standard)      0'
               write(15,*)'isct_vtst(1=vtst_sct,0=sct) 1'
               write(15,*)'zerocurvature(1)            0'
               write(15,*)'reduced_mass                1.0'
               write(15,*)'minimum_frequency            50'
               write(15,*)'anim_freq(if_Maxstep=1)       2'
               write(15,*)'onlyrotors(0=yes,1=no)        0'
               close (unit=15,status='keep')

               open (unit=99,status='unknown')
               rewind (99)

               if (ispecies.eq.1) then
                  command1='cp -f output/hrdata4proj_reac1.dat
     $                  hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                  geom.log ',natom
               endif
               if (ispecies.eq.2) then
                  command1='cp -f output/hrdata4proj_reac2.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.3) then
                  command1='cp -f output/hrdata4proj_prod1.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.4) then
                  command1='cp -f output/hrdata4proj_prod2.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.5) then
                  command1='cp -f output/hrdata4proj_wellr.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                           geom.log ',natom
               endif
               if (ispecies.eq.6) then
                  command1='cp -f output/hrdata4proj_wellr.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.11) then
                  command1='cp -f output/hrdata4proj_reacs.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.12) then
                  command1='cp -f output/hrdata4proj_prods.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.0) then
                  command1='cp -f output/hrdata4proj_ts.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                           geom.log ',natom
               endif

               rewind (99)
               read (99,2001) command1
               close (99)
               call commrun(command1)
               command1='RPHt.exe'
               call commrun(command1)

               if(imhrfr.eq.1)then
                  nfreq=3*natom-6-ntau_fr
                  if (natom.eq.2) nfreq=1
                  if (ispecies.eq.0) nfreq=3*natom-6-ntau_fr-1
                  
                  open (unit=15,file='hrproj_freq.dat',status='old')
                  if (ispecies.ne.0) then
                     do k=1,nfreq
                        read (15,*) freq1d(iscanhinda,k)
                     enddo
                  else
                     do k=1,nfreq+6+ntau_fr+1
                        read (15,*) freq1d(iscanhinda,k)
                     enddo
                  endif
                  close (15)
               endif

c      write (19,*) '    Frequencies[1/cm] ',nfreq

               if(irepeat.eq.0)then
                  comline1=comlineref1
                  comline2=comlineref2
               else
                  comline1=comline3
                  comline2=comline4
               endif
            endif


c            write (16,*) xint(ncoord),vtot         

            if(iscanhinda.le.ihalf) then
               xint(ncoord)= xstarta+ihind_stepa*iscanhinda
            else if (iscanhinda.gt.ihalf)then
               xint(ncoord)= xstarta-ihind_stepa*(iscanhinda-ihalf)
            endif

cc get back to starting point for reverse sweep

            if(iscanhinda.eq.ihalf+1) then
               do iint = 1 , ncoord-1
                  xint(iint) = xint_save(iint)
               enddo
            endif
            if (xint(ncoord).gt.360.)then
               xint(ncoord)= xint(ncoord)-360.
            endif

         enddo

cc end of 1D hindered rotors PES scan

         write (16,*) 'End of Scan for 1D hind rotor ',ihind

cc end of scan over 1D hindered rotors

cc now save results of 1D PES rotational scan

         open (unit=99,status='unknown')
         rewind (99)
         write(99,1201)word_1dhrout,ihind
         rewind (99)
         read (99,'(A35)') namefile
         close (99)
      
         open (unit=99,status='unknown')
         rewind (99)
         write(99,1202)word_1dhrout,ihind
         rewind (99)
         read (99,'(A35)') namefile2
         close (99)
         
         ref=0.
         if(imhrfr.eq.1) open (unit=17,file=namefile,status='unknown')
         open (unit=14,file=namefile2,status='unknown')
         if(imhrfr.eq.1) write(17,*)nhindsteps1Da(ihind)
         write(14,*)nhindsteps1Da(ihind)
         if(imhrfr.eq.1)write(17,8011)(k,k=1,nfreq)
         write(14,*)'nofreq'
         if(imhrfr.eq.1)write(17,*)
         write(14,*)
         do i=1,nhindsteps1Da(ihind)
            if(i.eq.1)then
              if(imhrfr.eq.1)write(17,8009)i,ref,(freq1d(i,k),k=1,nfreq)
               write(14,*)i,ref
c     write (19,8010) (freq(j),j=1,nfreq)
c
            else
               if(imhrfr.eq.1) then
                  write(17,8009)i,vref1D(i),(freq1d(i,k),k=1,nfreq)
               endif
               write(14,*)i,vref1D(i)
            endif
         enddo
         close(14)
         close(17)
      enddo

      write (16,*) 'End of Scan for all 1D hind rotors '
      write (16,*) '********************************** '

cc 2D hindered rotor section
cc scan each 2D hindered rotor at a time

cc start hindered rotor scan: 

      comsave1=comline1
      comsave2=comline2
      comlineref1=comline1
      comlineref2=comline2

      do ihind = 1 , nhind2d

cc for each 2Drotor start from first level theory

         comline1=comsave1
         comline2=comsave2

cc start rotational scan for hind rotor ihind

cc    re-order z-matrix, taking out torsion ihind      

         write (16,*) 're-ordered z-matrix of hind rotor'

         ncoord = 3*natom-6

         do iatom = 1 , ncoord
            xint_int(iatom)=0.
            intcoor_int(iatom)=''
         enddo

cc     first assign optimized hind parameter of first HR to last coordinate
         do iatom = 1 , ncoord
            if (intcoori(iatom).eq.hindlab2Da(ihind)) then
               xint_int(ncoord)=xinti(iatom)
            endif
         enddo

cc then update all coordinates for first rotor
         index=0
         do iatom = 1 , ncoord
            if (intcoori(iatom).eq.hindlab2Da(ihind)) goto 998
            index=index+1
            xint_int(index)=xinti(iatom)
            intcoor_int(index)=intcoori(iatom)
 998        continue
         enddo
         intcoor_int(ncoord)=hindlab2Da(ihind)

         do iatom = 1 , ncoord
            write (16,*) 'intc,xint ',intcoor_int(iatom),xint_int(iatom)
         enddo

cc     then do the same for second rotor
         do iatom = 1 , ncoord
            if (intcoor_int(iatom).eq.hindlab2Db(ihind)) then
               xint(ncoord)=xint_int(iatom)
            endif
         enddo

cc then update all coordinates for second rotor

         index=0
         do iatom = 1 , ncoord
            if (intcoor_int(iatom).eq.hindlab2Db(ihind)) goto 1998
            index=index+1
            xint(index)=xint_int(iatom)
            intcoor(index)=intcoor_int(iatom)
 1998        continue
         enddo
         intcoor(ncoord)=hindlab2Db(ihind)

         write(16,*)'coords for optimization'
         do iatom = 1 , ncoord
            write (16,*) 'intcoo,xint ',intcoor(iatom),xint(iatom)
         enddo


cc         start first cycle of 2D rotors

         ihind_stepa = (dhindmx2Da(ihind)-dhindmn2Da(ihind))
     +      /nhindsteps2Da(ihind)

         xstarta=xint(ncoord-1)
         xstartb=xint(ncoord)
         irepeat=0


cc save coordinates
            do iint = 1 , ncoord
               xint_save(iint) = xint(iint)
            enddo


         refen=0.

         do iscanhinda=1,nhindsteps2Da(ihind)

cc start scan on nhind scan points

            write (16,*) '2D hindered_rotor ',ihind
            write (16,*) '2D hindered_rotor cycle A',iscanhinda
            write (16,*) '2D hindered_rotor cycle B'

            ihind_stepb = (dhindmx2Db(ihind)-dhindmn2Db(ihind))
     +      /nhindsteps2Db(ihind)

c            ihalf=nhindsteps(ihind)/2
            xint(ncoord)=xstartb
c            iprog=0

            comline1=comsave1
            comline2=comsave2
            ntau=0
            ismp=iscanhindb
            ixyz=0
            ired=0
            ircons=2
            ilev=10
            irepeat=0

            do iint = 1 , ncoord
               xint_save(iint) = xint(iint)
            enddo
 111     continue

            do iscanhindb=1,nhindsteps2Db(ihind)

c               if(iscanhind.le.ihalf+1)iprog=iscanhind
c               if(iscanhind.gt.ihalf+1)iprog=nhindsteps(ihind)
c     $              +ihalf-iscanhind+2


               if(ilev1code.eq.1) then

                  call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $                 coord,vtot_0,vtot,freq,ifreq,ilin,ismp,comline1,
     $                 comline2,icharge,ispin,ircons,
     $                 atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires
     $                 ,ixyz,ired)

               else if (ilev1code.eq.2) then
                  numproc=numprochl
         
                  call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $                 coord,vtot,freq,ifreq,ilin,ismp,icharge,ispin,
     $                 ircons,
     $                 atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,
     $                 ixyz,ilev,ispecies)
                  numproc=numprocll
               endif


c            write(*,*)'after g09'

               write (16,*) xint(ncoord-1),xint(ncoord),vtot

               if(vtot.gt.1.0e10) then
                  
                  if(ispecies.eq.0.and.irepeat.eq.0)then


c                  if(irepeat.eq.0)then

cc switch method to determine HR potential
                     comline1=comline3
                     comline2=comline4
                  
                     ncoord = 3*natom-6
c                 write (16,*) 'int-2s ',atomcoo(isite)
c                 write (16,*) 'int-1s ',atomcoo(natomt1+1)
 
                     do icoo = 1 , ncoord
                        if (intcoor(icoo).eq.atomcoo(isite))then
                           xintt(ncoord-3)=xint_save(icoo)
                           intcoort(ncoord-3)=intcoor(icoo)
c                        write (16,*) 'int-2 ',intcoort(ncoord-2)
                        endif
                        if (intcoor(icoo).eq.'RTS')then
                           xintt(ncoord-2)=xint_save(icoo)
                           intcoort(ncoord-2)=intcoor(icoo)
c                        write (16,*) 'int-1 ',intcoort(ncoord-1)
                        endif
                     enddo

                     index=0

                     do icoo = 1 , ncoord-2
                        if (intcoor(icoo).eq.intcoort(ncoord-3).or.
     $                  intcoor(icoo).eq.intcoort(ncoord-2))goto 198
                        index=index+1
c     write(16,*)'index intcoort(index)',index,
c     $                           intcoort(index)
                        xintt(index)=xint_save(icoo)
                        intcoort(index)=intcoor(icoo)
 198                    continue
                     enddo

cc rewrite vectors
                     do icoo = 1 , ncoord-2
                        intcoor(icoo)=intcoort(icoo)
                        xint(icoo)=xintt(icoo)
                        xint_save(icoo)=xintt(icoo)
                        write (16,*) 'intcoor, xint ',intcoor(icoo)
     +                               ,xint(icoo)
                     enddo
c                     xint(ncoord-1)=xstarta
                     xint(ncoord)=xstartb
                     irepeat=1
                     ircons=4
                     goto 111
                  else
                     write(7,*)'this rotor did not converge'
                     write(7,*)'rotor number ',ihind
                     write(7,*)'rotor step ',iscanhindb
                     write(7,*)'using the last calculated energy'
                     write(7,*)'using the coordinates of first point'
                     write(7,*)'for next point'
                     if(ispecies.ne.0)then
                        do iint = 1 , ncoord-1
                           xint(iint) = xint_save(iint)
                        enddo
                     endif
                     open (unit=65,status='unknown')
                     read (65,*) vtot
                     close(65)
                     write (16,*) 'Not converged, using last energy'
                     write (16,*) xint(ncoord),vtot
                  endif

               endif

               vref2D(iscanhinda,iscanhindb)=vtot

               if (iscanhindb.eq.1.and.iscanhinda.eq.1) then
                  refen=vref2D(iscanhinda,iscanhindb)
               else
                  vref2D(iscanhinda,iscanhindb)=
     +           (vref2D(iscanhinda,iscanhindb)
     +                 -refen)*CAUTOKCAL
               endif

               if(vref2D(iscanhinda,iscanhindb).lt.-0.1.
     +                   and.iscanhindb.ne.1)then
                  write(7,*)'there is something wrong with this rotor'
                  write(7,*)'found a negative energy'
                  write(7,*)'2Drotor number ',ihind
                  write(7,*)'rotor step a',iscanhinda
                  write(7,*)'rotor step b',iscanhindb
                  write(7,*)'check the me input'
               endif

               if(imhrfr.eq.1) then
                  ixyz=0
                  ired=0
 
                  if(ilev1code.eq.1) then
                     comline1=comline5
                     comline2=comline6
                     atomlabel(1)=' '
                     atomlabel(2)=' '

                     call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $    coord,vtot_0,vtot,freq,ifreq,ilin,ismp,comline1,
     $    comline2,icharge,ispin,ircons,
     $    atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ired)

                     atomlabel(1)= atomlabel_save(1)
                     atomlabel(2)= atomlabel_save(2)

                  else if (ilev1code.eq.2) then

                     command1='cp -f fcmat.log ./geom.log '
                     call commrun(command1)

                  endif

cc now project fc matrix and re-eveluate frequencies
cc use RPHt file from 1D rotor analysis
               open (unit=15,file='hind_rot_head.dat',status='unknown')
               write(15,*)'Act_energy(kcal/mol):       0. '
               write(15,*)'Initial_Temperature:        200'
               write(15,*)'Temperature_steps:          40'
               write(15,*)'Temperature_increment:      50'
               write(15,*)'Delta_Energy_rea:           0.'
               write(15,*)'Delta_Energy_pro:           0.'
               write(15,*)'Maxstep:                    1'
               write(15,*)'Npointsint:                 5 '
               write(15,*)'Maxtdev:                    0.5'
               write(15,*)'Rearrange(1=yes,0=no)       1'
               write(15,*)'SaddlePoint                 1'
               write(15,*)'ds(1=noexp,0=standard)      0'
               write(15,*)'isct_vtst(1=vtst_sct,0=sct) 1'
               write(15,*)'zerocurvature(1)            0'
               write(15,*)'reduced_mass                1.0'
               write(15,*)'minimum_frequency            50'
               write(15,*)'anim_freq(if_Maxstep=1)       2'
               write(15,*)'onlyrotors(0=yes,1=no)        0'
               close (unit=15,status='keep')

               open (unit=99,status='unknown')

               if (ispecies.eq.1) then
                  command1='cp -f output/hrdata4proj_reac1.dat
     $ hrdata4proj.dat '
                  call commrun(command1)
                  write (99,*)'extract_data_projrot.sh 
     $ geom.log ',natom
               endif

               if (ispecies.eq.2) then
                  command1='cp -f output/hrdata4proj_reac2.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.3) then
                  command1='cp -f output/hrdata4proj_prod1.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.4) then
                  command1='cp -f output/hrdata4proj_prod2.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.5) then
                  command1='cp -f output/hrdata4proj_wellr.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                           geom.log ',natom
               endif
               if (ispecies.eq.6) then
                  command1='cp -f output/hrdata4proj_wellr.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.11) then
                  command1='cp -f output/hrdata4proj_reacs.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.12) then
                  command1='cp -f output/hrdata4proj_prods.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.0) then
                  command1='cp -f output/hrdata4proj_ts.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                           geom.log ',natom
               endif

               rewind (99)
               read (99,2001) command1
               close (99)

               call commrun(command1)
               command1='RPHt.exe'
               call commrun(command1)

               if(imhrfr.eq.1)then
                  nfreq=3*natom-6-ntau_fr
                  if (natom.eq.2) nfreq=1
                  if (ispecies.eq.0) nfreq=3*natom-6-ntau_fr-1
               
                  open (unit=15,file='hrproj_freq.dat',status='old')
                  if (ispecies.ne.0) then
                     do k=1,nfreq
                        read (15,*) freq2d(iscanhinda,iscanhindb,k)
                     enddo
                  else
                     do k=1,nfreq+6+ntau_fr+1
                        read (15,*) freq2d(iscanhinda,iscanhindb,k)
                     enddo
                  endif
                  close (15)
               endif

c      write (19,*) '    Frequencies[1/cm] ',nfreq

               if(irepeat.eq.0)then
                  comline1=comlineref1
                  comline2=comlineref2
               else
                  if(ispecies.eq.0)then
                     comline1=comline3
                     comline2=comline4
                  else
                     comline1=comlineref1
                     comline2=comlineref2
                  endif
               endif
               endif

c               write (16,*) xint(ncoord-1),xint(ncoord),vtot
            
               xint(ncoord)= xstartb+ihind_stepb*iscanhindb
            
               if (xint(ncoord).gt.360.)then
                  xint(ncoord)= xint(ncoord)-360.
               endif


            enddo
cc end of 2D hindered rotors PES scan for rotor b

c            write (16,*) xint(ncoord-1),vtot
            xint(ncoord-1)= xstarta+ihind_stepa*iscanhinda
            if (xint(ncoord-1).gt.360.)then
               xint(ncoord-1)= xint(ncoord-1)-360.
            endif
cc now update saved coords coordinates
         do iint = 1 , ncoord
            xint_save(iint) = xint(iint)
         enddo

         enddo
cc end of 2D hindered rotors PES scan for rotor a

         write (16,*) 'End of Scan for 2D hind rotor ',ihind

cc end of scan over 2D hindered rotors

cc now save results of 2D PES rotational scan

         open (unit=99,status='unknown')
         rewind (99)
         write(99,1201)word_2dhrout,ihind
         rewind (99)
         read (99,'(A35)') namefile
         close (99)

         open (unit=99,status='unknown')
         rewind (99)
         write(99,1202)word_2dhrout,ihind
         rewind (99)
         read (99,'(A35)') namefile2
         close (99)
         
         ref=0.
         if(imhrfr.eq.1) open (unit=17,file=namefile,status='unknown')
         open (unit=14,file=namefile2,status='unknown')
         if(imhrfr.eq.1) write(17,*)nhindsteps2Da(ihind),
     $                   nhindsteps2Db(ihind)
         write(14,*)nhindsteps2Da(ihind),nhindsteps2Db(ihind)
         if(imhrfr.eq.1)write(17,8011)(k,k=1,nfreq)
         write(14,*)'nofreq'
         if(imhrfr.eq.1)write(17,*)
         write(14,*)
         do i=1,nhindsteps2Da(ihind)
            do j=1,nhindsteps2Db(ihind)
               if(i.eq.1.and.j.eq.1)then
                  if(imhrfr.eq.1)then
                     write(17,8010)i,j,ref,(freq2d(i,j,k),k=1,nfreq)
                  endif
                  write(14,*)i,j,ref
c      write (19,8010) (freq(j),j=1,nfreq)
c
               else
                  if(imhrfr.eq.1) then
                   write(17,8010)i,j,vref2D(i,j),
     $                    (freq2d(i,j,k),k=1,nfreq)
                  endif
                  write(14,*)i,j,vref2D(i,j)
               endif
            enddo
         enddo
         close(14)
         close(17)
      enddo

      write (16,*) 'End of Scan for all 2D hind rotors '
      write (16,*) '********************************** '


cc start scan for 3D hindered rotors

cc 3D hindered rotor section
cc scan each 3D hindered rotor at a time

cc start hindered rotor scan: 

      comsave1=comline1
      comsave2=comline2
      comlineref1=comline1
      comlineref2=comline2

      do ihind = 1 , nhind3d

cc for each 3Drotor start from first level theory

         comline1=comsave1
         comline2=comsave2

cc start rotational scan for hind rotor ihind

cc    re-order z-matrix, taking out torsion ihind      

         write (16,*) 're-ordered z-matrix of hind rotor'

         ncoord = 3*natom-6

         do iatom = 1 , ncoord
            xint_int(iatom)=0.
            intcoor_int(iatom)=''
            xint_int2(iatom)=0.
            intcoor_int2(iatom)=''
         enddo

cc     first assign optimized hind parameter of first HR to last coordinate
         do iatom = 1 , ncoord
            if (intcoori(iatom).eq.hindlab3Da(ihind)) then
               xint_int(ncoord)=xinti(iatom)
            endif
         enddo

cc then update all coordinates for first rotor
         index=0
         do iatom = 1 , ncoord
            if (intcoori(iatom).eq.hindlab3Da(ihind)) goto 5998
            index=index+1
            xint_int(index)=xinti(iatom)
            intcoor_int(index)=intcoori(iatom)
 5998       continue
         enddo
         intcoor_int(ncoord)=hindlab3Da(ihind)

         do iatom = 1 , ncoord
            write(16,*)'inta,xinta ',intcoor_int(iatom),xint_int(iatom)
         enddo

cc     then do the same for second rotor
         do iatom = 1 , ncoord
            if (intcoor_int(iatom).eq.hindlab3Db(ihind)) then
               xint_int2(ncoord)=xint_int(iatom)
            endif
         enddo

cc then update all coordinates for second rotor

         index=0
         do iatom = 1 , ncoord
            if (intcoor_int(iatom).eq.hindlab3Db(ihind)) goto 6998
            index=index+1
            xint_int2(index)=xint_int(iatom)
            intcoor_int2(index)=intcoor_int(iatom)
 6998        continue
         enddo
         intcoor_int2(ncoord)=hindlab3Db(ihind)

         do iatom = 1 , ncoord
            write(16,*)'intb,xintb ',intcoor_int(iatom),xint_int(iatom)
         enddo

cc     then do the same for the third rotor
         do iatom = 1 , ncoord
            if (intcoor_int2(iatom).eq.hindlab3Dc(ihind)) then
               xint(ncoord)=xint_int2(iatom)
            endif
         enddo

cc then update all coordinates for third rotor

         index=0
         do iatom = 1 , ncoord
            if (intcoor_int2(iatom).eq.hindlab3Dc(ihind)) goto 2998
            index=index+1
            xint(index)=xint_int2(iatom)
            intcoor(index)=intcoor_int2(iatom)
 2998        continue
         enddo
         intcoor(ncoord)=hindlab3Dc(ihind)

cc now we have all the coordinates in the desired order 

         write(16,*)'coords for optimization'
         do iatom = 1 , ncoord
            write (16,*) 'intcoo,xint ',intcoor(iatom),xint(iatom)
         enddo

cc         start first cycle of 3D rotors

         ihind_stepa = (dhindmx3Da(ihind)-dhindmn3Da(ihind))
     +      /nhindsteps3Da(ihind)

         xstarta=xint(ncoord-2)
         xstartb=xint(ncoord-1)
         xstartc=xint(ncoord)
         irepeat=0

cc save coordinates
         do iint = 1 , ncoord
            xint_save(iint) = xint(iint)
         enddo


         refen=0.

         do iscan_a=1,nhindsteps3Da(ihind)

cc start scan on nhind scan points

            ihind_stepb = (dhindmx3Db(ihind)-dhindmn3Db(ihind))
     +      /nhindsteps3Db(ihind)

c            ihalf=nhindsteps(ihind)/2
            xint(ncoord-1)=xstartb
c            iprog=0



            do iscan_b=1,nhindsteps3Db(ihind)

            write (16,*) '3D hindered_rotor ',ihind
            write (16,*) '3D hindered_rotor cycle A',iscan_a
            write (16,*) '3D hindered_rotor cycle B',iscan_b

c               if(iscanhind.le.ihalf+1)iprog=iscanhind
c               if(iscanhind.gt.ihalf+1)iprog=nhindsteps(ihind)
c     $              +ihalf-iscanhind+2


               ihind_stepc = (dhindmx3Dc(ihind)-dhindmn3Dc(ihind))
     +              /nhindsteps3Dc(ihind)

               xint(ncoord)=xstartc

               comline1=comsave1
               comline2=comsave2
               ntau=0
               ismp=iscan_c
               ixyz=0
               ired=0
               ircons=3
               ilev=10
               irepeat=0

               do iint = 1 , ncoord
                  xint_save(iint) = xint(iint)
               enddo

 1111     continue


               do iscan_c=1,nhindsteps3Dc(ihind)


c               stop
c            write(*,*)'before g09'

                  if(ilev1code.eq.1) then
 
                     call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $                    coord,vtot_0,vtot,freq,ifreq,ilin,ismp,
     $                    comline1,comline2,icharge,ispin,ircons,
     $                    atomlabel,intcoor,bislab,tauopt,xint,abcrot
     $                    ,ires,ixyz,ired)

                  else if (ilev1code.eq.2) then

                     numproc=numprochl
                     call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $                    coord,vtot,freq,ifreq,ilin,ismp,
     $                    icharge,ispin,ircons,
     $                    atomlabel,intcoor,bislab,tauopt,xint,abcrot,
     $                    ires,ixyz,ilev,ispecies)
                     numproc=numprocll

                  endif

c            write(*,*)'after g09'

                  write (16,1717) xint(ncoord-2),xint(ncoord-1)
     $                 ,xint(ncoord),vtot

                  if(vtot.gt.1.0e10) then
                  
                   if(ispecies.eq.0.and.irepeat.eq.0)then


c                  if(irepeat.eq.0)then

cc switch method to determine HR potential
                     comline1=comline3
                     comline2=comline4
                  
                     ncoord = 3*natom-6
c                 write (16,*) 'int-2s ',atomcoo(isite)
c                 write (16,*) 'int-1s ',atomcoo(natomt1+1)
 
                     do icoo = 1 , ncoord
                        if (intcoor(icoo).eq.atomcoo(isite))then
                           xintt(ncoord-4)=xint_save(icoo)
                           intcoort(ncoord-4)=intcoor(icoo)
c                        write (16,*) 'int-2 ',intcoort(ncoord-2)
                        endif
                        if (intcoor(icoo).eq.'RTS')then
                           xintt(ncoord-3)=xint_save(icoo)
                           intcoort(ncoord-3)=intcoor(icoo)
c                        write (16,*) 'int-1 ',intcoort(ncoord-1)
                        endif
                     enddo

                     index=0

                     do icoo = 1 , ncoord-3
                        if (intcoor(icoo).eq.intcoort(ncoord-4).or.
     $                  intcoor(icoo).eq.intcoort(ncoord-3))goto 1198
                        index=index+1
c     write(16,*)'index intcoort(index)',index,
c     $                           intcoort(index)
                        xintt(index)=xint_save(icoo)
                        intcoort(index)=intcoor(icoo)
 1198                   continue
                     enddo

cc rewrite vectors
                     do icoo = 1 , ncoord-3
                        intcoor(icoo)=intcoort(icoo)
                        xint(icoo)=xintt(icoo)
                        xint_save(icoo)=xintt(icoo)
                        write (16,*) 'intcoor, xint ',intcoor(icoo)
     +                               ,xint(icoo)
                     enddo
c                     xint(ncoord-2)=xstarta
c                     xint(ncoord-1)=xstartb
                     xint(ncoord)=xstartc
                     irepeat=1
                     ircons=5
                     goto 1111
                  else
                     write(7,*)'this rotor did not converge'
                     write(7,*)'rotor number ',ihind
                     write(7,*)'rotor step ',iscan_c
                     write(7,*)'using the last calculated energy'
                     write(7,*)'using the coordinates of first point'
                     write(7,*)'for next point'
                     if(ispecies.ne.0)then
                        do iint = 1 , ncoord-3
                           xint(iint) = xint_save(iint)
                        enddo
                     endif

                     open (unit=65,status='unknown')
                     read (65,*) vtot
                     close(65)
                     write (16,*) 'Not converged, using last energy'
                     write (16,*) xint(ncoord),vtot
                  endif

               endif

               vref3D(iscan_a,iscan_b,iscan_c)=vtot

               if (iscan_b.eq.1.and.iscan_a.eq.1.and.iscan_c.eq.1) then
                  refen=vref3D(iscan_a,iscan_b,iscan_c)
               else
                  vref3D(iscan_a,iscan_b,iscan_c)=
     +           (vref3D(iscan_a,iscan_b,iscan_c)
     +                 -refen)*CAUTOKCAL
               endif

               if(vref3D(iscan_a,iscan_b,iscan_c).lt.-0.1.
     +                   and.iscan_c.ne.1)then
                  write(7,*)'there is something wrong with this rotor'
                  write(7,*)'found a negative energy'
                  write(7,*)'3Drotor number ',ihind
                  write(7,*)'rotor step a',iscan_a
                  write(7,*)'rotor step b',iscan_b
                  write(7,*)'rotor step c',iscan_c
                  write(7,*)'check the me input'
               endif



c            if(iscanhind.le.ihalf) then
c               xint(ncoord)= xstart+ihind_step*iscanhind
c            else if (iscanhind.gt.ihalf)then
c               xint(ncoord)= xstart-ihind_step*(iscanhind-ihalf)
c            endif
cc reverse sweep not used in 2D scans
cc get back to starting point for reverse sweep
c            if(iscanhind.eq.ihalf+1) then
c               do iint = 1 , ncoord-1
c                  xint(iint) = xint_save(iint)
c               enddo
c            endif


cc save feature not activated for 2D scans
cc save hindered rotor optimized geometry zzzz

c         command1='newzmat       
c     +    -ichk -oxyz tmp.chk geom.xyz           '
c         call commrun(command1)
c         open (unit=98,file='headgeom.tmp',status='unknown')
c         write(98,*)natom
c         write(98,*)vtot
c         close(98)
c         command1='cat headgeom.tmp geom.xyz > geom1.xyz '
c         call commrun(command1)
c
c         open (unit=99,status='unknown')
c         rewind (99)
c         write(99,1200)ispecies,ihind,iscanhind
c         rewind (99)
c         read (99,2001) command1
c         close (99)
c         call commrun(command1)
c 1200 format (" cp -f geom1.xyz 
c     $          hr_geoms/geom_isp"I0.2"_hr"I0.2"_hpt"I0.2".xyz  ")

c     calculate hessian for optimized point

               if(imhrfr.eq.1) then
                  ixyz=0
                  ired=0
                  
                  if(ilev1code.eq.1) then
                     comline1=comline5
                     comline2=comline6
                     atomlabel(1)=' '
                     atomlabel(2)=' '

                     call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $                    coord,vtot_0,vtot,freq,ifreq,ilin,ismp,
     $                    comline1,comline2,icharge,ispin,ircons,
     $                    atomlabel,intcoor,bislab,tauopt,xint,abcrot,
     $                    ires,ixyz,ired)

                     atomlabel(1)= atomlabel_save(1)
                     atomlabel(2)= atomlabel_save(2)

                  else if (ilev1code.eq.2) then

                     command1='cp -f fcmat.log ./geom.log '
                     call commrun(command1)

                  endif
             

cc now project fc matrix and re-eveluate frequencies
cc use RPHt file from 1D rotor analysis
               open (unit=15,file='hind_rot_head.dat',status='unknown')
               write(15,*)'Act_energy(kcal/mol):       0. '
               write(15,*)'Initial_Temperature:        200'
               write(15,*)'Temperature_steps:          40'
               write(15,*)'Temperature_increment:      50'
               write(15,*)'Delta_Energy_rea:           0.'
               write(15,*)'Delta_Energy_pro:           0.'
               write(15,*)'Maxstep:                    1'
               write(15,*)'Npointsint:                 5 '
               write(15,*)'Maxtdev:                    0.5'
               write(15,*)'Rearrange(1=yes,0=no)       1'
               write(15,*)'SaddlePoint                 1'
               write(15,*)'ds(1=noexp,0=standard)      0'
               write(15,*)'isct_vtst(1=vtst_sct,0=sct) 1'
               write(15,*)'zerocurvature(1)            0'
               write(15,*)'reduced_mass                1.0'
               write(15,*)'minimum_frequency            50'
               write(15,*)'anim_freq(if_Maxstep=1)       2'
               write(15,*)'onlyrotors(0=yes,1=no)        0'
               close (unit=15,status='keep')

               open (unit=99,status='unknown')
               rewind (99)

               if (ispecies.eq.1) then
                  command1='cp -f output/hrdata4proj_reac1.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.2) then
                  command1='cp -f output/hrdata4proj_reac2.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.3) then
                  command1='cp -f output/hrdata4proj_prod1.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.4) then
                  command1='cp -f output/hrdata4proj_prod2.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.5) then
                  command1='cp -f output/hrdata4proj_wellr.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                           geom.log ',natom
               endif
               if (ispecies.eq.6) then
                  command1='cp -f output/hrdata4proj_wellr.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.11) then
                  command1='cp -f output/hrdata4proj_reacs.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.12) then
                  command1='cp -f output/hrdata4proj_prods.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                            geom.log ',natom
               endif
               if (ispecies.eq.0) then
                  command1='cp -f output/hrdata4proj_ts.dat
     $                      hrdata4proj.dat '
                  call commrun(command1)
                  write (99,2000)'extract_data_projrot.sh 
     $                           geom.log ',natom
               endif

               rewind (99)
               read (99,2001) command1
               close (99)
               call commrun(command1)
               command1='RPHt.exe'
               call commrun(command1)

               if(imhrfr.eq.1)then
                  nfreq=3*natom-6-ntau_fr
                  if (natom.eq.2) nfreq=1
                  if (ispecies.eq.0) nfreq=3*natom-6-ntau_fr-1
               
                  open (unit=15,file='hrproj_freq.dat',status='old')
                  if (ispecies.ne.0) then
                     do k=1,nfreq
                        read (15,*) freq3d(iscan_a,iscan_b,iscan_c,k)
                     enddo
                  else
                     do k=1,nfreq+6+ntau_fr+1
                        read (15,*) freq3d(iscan_a,iscan_b,iscan_c,k)
                     enddo
                  endif
                  close (15)
               endif

c      write (19,*) '    Frequencies[1/cm] ',nfreq

               if(irepeat.eq.0)then
                  comline1=comlineref1
                  comline2=comlineref2
               else
                  if(ispecies.eq.0)then
                     comline1=comline3
                     comline2=comline4
                  else
                     comline1=comlineref1
                     comline2=comlineref2
                  endif
               endif
               endif
               

c               write (16,1717) xint(ncoord-2),xint(ncoord-1),
c     $                      xint(ncoord),vtot
            
               xint(ncoord)= xstartc+ihind_stepc*iscan_c
            
               if (xint(ncoord).gt.360.)then
                  xint(ncoord)= xint(ncoord)-360.
               endif


            enddo
cc end of 3D hindered rotors PES scan for rotor c


            xint(ncoord-1)= xstartb+ihind_stepb*iscan_b
            if (xint(ncoord-1).gt.360.)then
               xint(ncoord-1)= xint(ncoord-1)-360.
            endif
         enddo

cc end of 3D hindered rotors PES scan for rotor b

c            write (16,*) xint(ncoord-1),vtot
            xint(ncoord-2)= xstarta+ihind_stepa*iscan_a
            if (xint(ncoord-2).gt.360.)then
               xint(ncoord-2)= xint(ncoord-2)-360.
            endif
         enddo
cc end of 3D hindered rotors PES scan for rotor a

         write (16,*) 'End of Scan for 3D hind rotor ',ihind

cc end of scan over 3D hindered multirotor ihind

cc now save results of 3D PES rotational scan for multirotor ihind

         open (unit=99,status='unknown')
         rewind (99)
         write(99,1201)word_3dhrout,ihind
         rewind (99)
         read (99,'(A35)') namefile
         close (99)

         open (unit=99,status='unknown')
         rewind (99)
         write(99,1202)word_3dhrout,ihind
         rewind (99)
         read (99,'(A35)') namefile2
         close (99)
         
         ref=0.
         if(imhrfr.eq.1) open (unit=17,file=namefile,status='unknown')
         open (unit=14,file=namefile2,status='unknown')
         if(imhrfr.eq.1) write(17,*)nhindsteps3Da(ihind),
     $                   nhindsteps3Db(ihind),nhindsteps3Dc(ihind)
         write(14,*)nhindsteps3Da(ihind),nhindsteps3Db(ihind),
     $              nhindsteps3Dc(ihind)
         if(imhrfr.eq.1)write(17,8011)(k,k=1,nfreq)
         write(14,*)'nofreq'
         if(imhrfr.eq.1)write(17,*)
         write(14,*)
         do i=1,nhindsteps3Da(ihind)
            do j=1,nhindsteps3Db(ihind)
               do ij=1,nhindsteps3Dc(ihind)
                  if(i.eq.1.and.j.eq.1.and.ij.eq.1)then
                     if(imhrfr.eq.1)write(17,9010)i,j,ij,ref,
     $                    (freq3d(i,j,ij,k),k=1,nfreq)
                     write(14,*)i,j,ij,ref
c      write (19,8010) (freq(j),j=1,nfreq)
c
                  else
                     if(imhrfr.eq.1) then
                        write(17,9010)i,j,ij,vref3D(i,j,ij),
     $                       (freq3d(i,j,ij,k),k=1,nfreq)
                     endif
                     write(14,*)i,j,ij,vref3D(i,j,ij)
                  endif
               enddo
            enddo
         enddo
         close(14)
         close(17)
      enddo

      write (16,*) 'End of Scan for all 3D hind rotors '
      write (16,*) '********************************** '


cc write mr input file for me

      open (unit=99,status='unknown')

cc first write 1D hindered rotor section
      do ihind=1,nhind1d
         write(20,*)"   Core MultiRotor"
         write(26,*)"   Core MultiRotor"
         write(20,*)"   SymmetryFactor  1.0"
         write(26,*)"   SymmetryFactor  1.0"
         write(20,*)"   InterpolationEnergyMax[kcal/mol]  100"
         write(26,*)"   InterpolationEnergyMax[kcal/mol]  100"
         rewind (99)
         write(99,1201)word_1dhrout,ihind
         rewind (99)
         read (99,'(A35)') namefile
         rewind (99)
         write(99,1202)word_1dhrout,ihind
         rewind (99)
         read (99,'(A35)') namefile2
         close (99)
         write(20,*)"   PotentialEnergySurface[kcal/mol] ",namefile
         write(26,*)"   PotentialEnergySurface[kcal/mol] ",namefile2
         write(20,*)"   QuantumLevelEnergyMax[kcal/mol]         5"
         write(26,*)"   QuantumLevelEnergyMax[kcal/mol]         5"
         do i=1,nhind
            if(i1DtopA(ihind).eq.i)then
               write(20,*)"   InternalRotation"
               write(26,*)"   InternalRotation"
            endif
            call LineRead(0)
            do while (WORD.NE.'END')
               read(19,'(A100)')buffer
               rewind (99)
               write(99,*)buffer
               rewind (99)
               call LineRead(99)
c               write(*,*)'test file 19 word is', word
               if(i1DtopA(ihind).eq.i)then
                  if(WORD.eq.'GROUP')write(20,'    (A60)')buffer
                  if(WORD.eq.'GROUP')write(26,'    (A60)')buffer
                  if(WORD.eq.'AXIS')write(20,'    (A60)')buffer
                  if(WORD.eq.'AXIS')write(26,'    (A60)')buffer
                  if(WORD.eq.'SYMMETRY')write(20,'    (A60)')buffer
                  if(WORD.eq.'SYMMETRY')write(26,'    (A60)')buffer
               endif
            enddo
            if(i1DtopA(ihind).eq.i)then
               write(20,*)"     MassExpansionSize 5"
               write(26,*)"     MassExpansionSize 5"
               write(20,*)"     HamiltonSizeMin  13"
               write(26,*)"     HamiltonSizeMin  13"
               write(20,*)"     HamiltonSizeMax  101"
               write(26,*)"     HamiltonSizeMax  101"
               write(20,*)"     GridSize         100"
               write(26,*)"     GridSize         100"
               write(20,*)"   End"
               write(26,*)"   End"
            endif
         enddo
         write(20,*)"   End"
         write(26,*)"   End"
      enddo

cc then write 2D hindered rotor section
      do ihind=1,nhind2d
         write(20,*)"   Core MultiRotor"
         write(26,*)"   Core MultiRotor"
         write(20,*)"   SymmetryFactor  1.0"
         write(26,*)"   SymmetryFactor  1.0"
         write(20,*)"   InterpolationEnergyMax[kcal/mol]  100"
         write(26,*)"   InterpolationEnergyMax[kcal/mol]  100"
         rewind (99)
         write(99,1201)word_2dhrout,ihind
         rewind (99)
         read (99,'(A35)') namefile
         rewind (99)
         write(99,1202)word_2dhrout,ihind
         rewind (99)
         read (99,'(A35)') namefile2
         close (99)
         write(20,*)"   PotentialEnergySurface[kcal/mol] ",namefile
         write(26,*)"   PotentialEnergySurface[kcal/mol] ",namefile2
         write(20,*)"   QuantumLevelEnergyMax[kcal/mol]         5"
         write(26,*)"   QuantumLevelEnergyMax[kcal/mol]         5"
         do i=1,nhind
            if(i2DtopA(ihind).eq.i.or.i2DtopB(ihind).eq.i)then
               write(20,*)"   InternalRotation"
               write(26,*)"   InternalRotation"
            endif
            call LineRead(0)
            do while (WORD.NE.'END')
               read(19,'(A100)')buffer
               rewind (99)
               write(99,*)buffer
               rewind (99)
               call LineRead(99)
c               write(*,*)'test file 19 word is', word
               if(i2DtopA(ihind).eq.i.or.i2DtopB(ihind).eq.i)then
                  if(WORD.eq.'GROUP')write(20,'    (A60)')buffer
                  if(WORD.eq.'GROUP')write(26,'    (A60)')buffer
                  if(WORD.eq.'AXIS')write(20,'    (A60)')buffer
                  if(WORD.eq.'AXIS')write(26,'    (A60)')buffer
                  if(WORD.eq.'SYMMETRY')write(20,'    (A60)')buffer
                  if(WORD.eq.'SYMMETRY')write(26,'    (A60)')buffer
               endif
            enddo
            if(i2DtopA(ihind).eq.i.or.i2DtopB(ihind).eq.i)then
               write(20,*)"     MassExpansionSize 5"
               write(26,*)"     MassExpansionSize 5"
               write(20,*)"     HamiltonSizeMin  13"
               write(26,*)"     HamiltonSizeMin  13"
               write(20,*)"     HamiltonSizeMax  101"
               write(26,*)"     HamiltonSizeMax  101"
               write(20,*)"     GridSize         100"
               write(26,*)"     GridSize         100"
               write(20,*)"   End"
               write(26,*)"   End"
            endif
         enddo
         write(20,*)"   End"
         write(26,*)"   End"
      enddo

cc then write 3D hindered rotor section
      do ihind=1,nhind3d
         write(20,*)"   Core MultiRotor"
         write(26,*)"   Core MultiRotor"
         write(20,*)"   SymmetryFactor 1.0"
         write(26,*)"   SymmetryFactor 1.0"
         write(20,*)"   InterpolationEnergyMax[kcal/mol]  100"
         write(26,*)"   InterpolationEnergyMax[kcal/mol]  100"
         rewind (99)
         write(99,1201)word_3dhrout,ihind
         rewind (99)
         read (99,'(A35)') namefile
         rewind (99)
         write(99,1202)word_3dhrout,ihind
         rewind (99)
         read (99,'(A35)') namefile2
         write(20,*)"   PotentialEnergySurface[kcal/mol] ",namefile
         write(26,*)"   PotentialEnergySurface[kcal/mol] ",namefile2
         write(20,*)"   QuantumLevelEnergyMax[kcal/mol]         5"
         write(26,*)"   QuantumLevelEnergyMax[kcal/mol]         5"
         do i=1,nhind
            if(i3DtopA(ihind).eq.i.or.i3DtopB(ihind).eq.i.or.
     $          i3DtopC(ihind).eq.i)then
               write(20,*)"   InternalRotation"
               write(26,*)"   InternalRotation"
            endif
            call LineRead(0)
            do while (WORD.NE.'END')
               read(19,'(A100)')buffer
               rewind (99)
               write(99,*)buffer
               rewind (99)
               call LineRead(99)
c               write(*,*)'test file 19 word is', word
               if(i3DtopA(ihind).eq.i.or.i3DtopB(ihind).eq.i.or.
     $              i3DtopC(ihind).eq.i)then
                  if(WORD.eq.'GROUP')write(20,'    (A60)')buffer
                  if(WORD.eq.'GROUP')write(26,'    (A60)')buffer
                  if(WORD.eq.'AXIS')write(20,'    (A60)')buffer
                  if(WORD.eq.'AXIS')write(26,'    (A60)')buffer
                  if(WORD.eq.'SYMMETRY')write(20,'    (A60)')buffer
                  if(WORD.eq.'SYMMETRY')write(26,'    (A60)')buffer
               endif
            enddo
            if(i3DtopA(ihind).eq.i.or.i3DtopB(ihind).eq.i.or.
     $           i3DtopC(ihind).eq.i)then
               write(20,*)"     MassExpansionSize 5"
               write(26,*)"     MassExpansionSize 5"
               write(20,*)"     HamiltonSizeMin  13"
               write(26,*)"     HamiltonSizeMin  13"
               write(20,*)"     HamiltonSizeMax  101"
               write(26,*)"     HamiltonSizeMax  101"
               write(20,*)"     GridSize         100"
               write(26,*)"     GridSize         100"
               write(20,*)"   End"
               write(26,*)"   End"
            endif
         enddo
         write(20,*)"   End"
         write(26,*)"   End"
      enddo

cc now write 1D hr potentials
      rewind(19)
      ipot=0
      do i=1,nhind
         iwrite=1
         do ij=1,nhind1d
            if(i1DtopA(ij).eq.i) iwrite=0
         enddo
         do ij=1,nhind2d
            if(i2DtopA(ij).eq.i) iwrite=0
            if(i2DtopB(ij).eq.i) iwrite=0
         enddo
         do ij=1,nhind3d
            if(i3DtopA(ij).eq.i) iwrite=0
            if(i3DtopB(ij).eq.i) iwrite=0
            if(i3DtopC(ij).eq.i) iwrite=0
         enddo
         call LineRead(0)
         do while (WORD.NE.'END')
            if(ipot.eq.0) then
               read(19,'(A100)')buffer
               rewind (99)
               write(99,*)buffer
               rewind (99)
               call LineRead(99)
               write(*,*)'test file 19bis word is', word
            else
               write(*,*)'ipot is', ipot
               read (19,*) (vpot(j),j=1,ipot)
            endif
            if(iwrite.eq.1)then
               if(ipot.eq.0) then
                  write(20,'(A60)')buffer
                  write(26,'(A60)')buffer
               else
                  write (20,7111) (vpot(j),j=1,ipot)
                  write (26,7111) (vpot(j),j=1,ipot)
                  ipot=0
                  WORD=''
               endif
               if(WORD.EQ.'POTENTIAL[KCAL/MOL]')then
                  rewind(99)
                  write(99,*)word2
                  rewind(99)
                  read(99,*)ipot
               endif
            endif
         enddo
      enddo
      close(99)
      close(19)
      close(20)
      close(26)
cc now copy mr file to hr file and update other files accordingly

      if (ispecies.eq.1) then
         if(imhrfr.eq.1)then
            command1='cp -f me_files/reac1_fr.me
     $            me_files/reac1_1dfr.me '
            call commrun(command1)
            command1='cp -f me_files/reac1_nofr.me
     $            me_files/reac1_fr.me '
            call commrun(command1)
            command1='cp -f me_files/reac1_mdhr.me
     $            me_files/reac1_hr.me '
         else
            command1='cp -f me_files/reac1_mdhr_nofr.me
     $            me_files/reac1_hr.me '
         endif
         call commrun(command1)

         command1='head -n -3  me_files/reac1_1dge.me
     $            > me_files/reac1_ge.me '
         call commrun(command1)

      endif
      if (ispecies.eq.2) then
         if(imhrfr.eq.1)then
            command1='cp -f me_files/reac2_fr.me
     $            me_files/reac2_1dfr.me '
            call commrun(command1)
            command1='cp -f me_files/reac2_nofr.me
     $            me_files/reac2_fr.me '
            call commrun(command1)

            command1='cp -f me_files/reac2_mdhr.me
     $            me_files/reac2_hr.me '
         else
            command1='cp -f me_files/reac2_mdhr_nofr.me
     $            me_files/reac2_hr.me '
         endif
         call commrun(command1)

         command1='head -n -3  me_files/reac2_1dge.me
     $            > me_files/reac2_ge.me '
         call commrun(command1)

      endif

      if (ispecies.eq.3) then
         if(imhrfr.eq.1)then
            command1='cp -f me_files/prod1_fr.me
     $            me_files/prod1_1dfr.me '
            call commrun(command1)
            command1='cp -f me_files/prod1_nofr.me
     $            me_files/prod1_fr.me '
            call commrun(command1)

            command1='cp -f me_files/prod1_mdhr.me
     $            me_files/prod1_hr.me '
         else
            command1='cp -f me_files/prod1_mdhr_nofr.me
     $            me_files/prod1_hr.me '
         endif
         call commrun(command1)

         command1='head -n -3  me_files/prod1_1dge.me
     $            > me_files/prod1_ge.me '
         call commrun(command1)
      endif

      if (ispecies.eq.4) then
         if(imhrfr.eq.1)then
            command1='cp -f me_files/prod2_fr.me
     $            me_files/prod2_1dfr.me '
            call commrun(command1)
            command1='cp -f me_files/prod2_nofr.me
     $            me_files/prod2_fr.me '
            call commrun(command1)

            command1='cp -f me_files/prod2_mdhr.me
     $            me_files/prod2_hr.me '
         else
            command1='cp -f me_files/prod2_mdhr_nofr.me
     $            me_files/prod2_hr.me '
         endif
         call commrun(command1)

         command1='head -n -3  me_files/prod2_1dge.me
     $            > me_files/prod2_ge.me '
         call commrun(command1)
      endif

      if (ispecies.eq.5) then
         if(imhrfr.eq.1)then
            command1='cp -f me_files/wellr_fr.me
     $            me_files/wellr_1dfr.me '
            call commrun(command1)
            command1='cp -f me_files/wellr_nofr.me
     $            me_files/wellr_fr.me '
            call commrun(command1)

            command1='cp -f me_files/wellr_mdhr.me
     $            me_files/wellr_hr.me '
         else
            command1='cp -f me_files/wellr_mdhr_nofr.me
     $            me_files/wellr_hr.me '
         endif
         call commrun(command1)
         command1='head -n -3  me_files/wellr_1dge.me
     $            > me_files/wellr_ge.me '
         call commrun(command1)
      endif

      if (ispecies.eq.6) then
         if(imhrfr.eq.1)then
            command1='cp -f me_files/wellp_fr.me
     $            me_files/wellp_1dfr.me '
            call commrun(command1)
            command1='cp -f me_files/wellp_nofr.me
     $            me_files/wellp_fr.me '
            call commrun(command1)

            command1='cp -f me_files/wellp_mdhr.me
     $            me_files/wellp_hr.me '
         else
            command1='cp -f me_files/wellp_mdhr_nofr.me
     $            me_files/wellp_hr.me '
         endif
         call commrun(command1)
         command1='head -n -3  me_files/wellp_1dge.me
     $            > me_files/wellp_ge.me '
         call commrun(command1)
      endif

      if (ispecies.eq.11) then
         if(imhrfr.eq.1)then
            command1='cp -f me_files/reacs_fr.me
     $            me_files/reacs_1dfr.me '
            call commrun(command1)
            command1='cp -f me_files/reacs_nofr.me
     $            me_files/reacs_fr.me '
            call commrun(command1)

            command1='cp -f me_files/reacs_mdhr.me
     $            me_files/reacs_hr.me '
         else
            command1='cp -f me_files/reacs_mdhr_nofr.me
     $me_files/reacs_hr.me '
         endif
         call commrun(command1)
         command1='cp -f me_files/reacs_ge.me
     $            me_files/reacs_1dge.me '
         call commrun(command1)
         command1='head -n -3  me_files/reacs_1dge.me
     $            > me_files/reacs_ge.me '
         call commrun(command1)
      endif

      if (ispecies.eq.12) then
         if(imhrfr.eq.1)then
            command1='cp -f me_files/prods_fr.me
     $            me_files/prods_1dfr.me '
            call commrun(command1)
            command1='cp -f me_files/prods_nofr.me
     $            me_files/prods_fr.me '
            call commrun(command1)

            command1='cp -f me_files/prods_mdhr.me
     $            me_files/prods_hr.me '
         else
            command1='cp -f me_files/prods_mdhr_nofr.me
     $            me_files/prods_hr.me '
         endif
         call commrun(command1)
         command1='cp -f me_files/prods_ge.me
     $            me_files/prods_1dge.me '
         call commrun(command1)
         command1='head -n -3  me_files/prods_1dge.me
     $            > me_files/prods_ge.me '
         call commrun(command1)
      endif

      if (ispecies.eq.0) then
         if(imhrfr.eq.1)then
            command1='cp -f me_files/ts_fr.me
     $            me_files/ts_1dfr.me '
            call commrun(command1)
            command1='cp -f me_files/ts_nofr.me
     $            me_files/ts_fr.me '
            call commrun(command1)
            command1='cp -f me_files/ts_mdhr.me
     $            me_files/ts_hr.me '
         else
            command1='cp -f me_files/ts_mdhr_nofr.me
     $            me_files/ts_hr.me '
         endif
         call commrun(command1)
         command1='head -n -3  me_files/ts_1dge.me
     $            > me_files/ts_ge.me '
         call commrun(command1)
      endif

 1201 format ( "./me_files/"A10"_"I0.2".dat")
 1202 format ( "./me_files/"A10"_nofr"I0.2".dat")
 7111 format (100(1X,f7.2))
 1717 format (3(1X,f7.2),1X,F10.5)
 2000 format (A80,1X,I10)
 2001 format (A160)
 8009 format (I3,1x,F7.3,1x,100G12.5)
 8010 format (I3,1x,I3,1x,F7.3,1x,100G12.5)
 8011 format (1x,100G10.3)
 9010 format (I3,1x,I3,1x,I3,1x,F7.3,1x,100G12.5)


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine symmetry(ispecies)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

cc    this is the hindered rotor section

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'


      dimension dhindmn(nhindmx),dhindmx(nhindmx),
     $ freq(nmdmx),nhindsteps(nhindmx),nhindsymm(nhindmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim),xintt(3*natommx),xint_save(3*natommx)
      dimension vref(noptmx),tauo(ntaumx,noptmx),
     $ xinto(3*natommx,noptmx),dist(noptmx)
      dimension tauopt(ntaumx)
      dimension energy(noptmx,nhindmx)
      dimension energymin(noptmx,nhindmx)
      dimension energyopt(noptmx,nhindmx)
      dimension istepmin(noptmx,nhindmx)
      dimension nmin(nhindmx)
      dimension nrotiso(nhindmx)
      dimension tau(ntaumx)
c      dimension iatomtype(natommx)
c      dimension itop_ind(natommx,nhindmx),
c     $          itop_ind_nodum(natommx,nhindmx)
c      dimension idummy(natommx),
c     $          ngroup(nhindmx,natommx)
c      dimension igrouptot(nhindmx),ipivotA(nhindmx),ipivotB(nhindmx)
c      dimension taumn(ntaumx),taumx(ntaumx)
c      dimension gelec(nelecmx),eelec(nelecmx)


      character*70 comline1,comline2,comline3,comline4,comsave1,comsave2
      character*60 atomlabel(natommx)
      character*4 atomname(natommx)
      character*4 atomconn(natommx)
      character*4 atomcoo(natommx)
      character*4 pivotA(nhindmx),pivotB(nhindmx)
      character*4 dicheck1(nhindmx),dicheck2(nhindmx)
      character*30 intcoor(3*natommx)
      character*30 intcoori(3*natommx)
      character*30 intcoort(3*natommx)
      character*20 bislab(ntaumx)
      character*20 hindlab(nhindmx)
      character*20 wordlab
      character*1 atfirstlet
      character*140 command1
      character*160 commandcopy
      character*30 gmem
      character*40 filename
      character*2 cjunk
      character*100 word_sym

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

cc energy threshold to check for optical isomers found 
cc through internal rotations

      enthreshold=0.00002

c input data
cs added option for species dependent hr level
cc assume level0 is used for molpro unless specified

      call Lineread(0)
      commandcopy='cp -f ./data/level0_molpro.dat 
     $       level0_molpro.dat'

      open (unit=21,file='./data/theory.dat',status='unknown')
      rewind(21)


c      WORD_SYM='SYMMETRY'
c      if(ispecies.eq.0) then
c         WORD_SYM='SYMMETRY_TS'
c      endif


      if (ispecies.eq.1) then
         open (unit=15,file='./data/reac1.dat',status='old')
         open (unit=16,file='./output/reac1_sym.out',status='unknown')
         open (unit=17,file='./geoms/reac1_l1.xyz',status='unknown')
         open (unit=18,file='./me_files/reac1_symm.me',status='unknown')
         WORD_SYM='SYMMETRY_1'
         do while (WORD.NE.'SYMMETRY_1')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_SYM = 'SYMMETRY'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/symm_reac1_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      else if (ispecies.eq.2) then
         open (unit=15,file='./data/reac2.dat',status='unknown')
         open (unit=16,file='./output/reac2_sym.out',status='unknown')
         open (unit=17,file='./geoms/reac2_l1.xyz',status='unknown')
         open (unit=18,file='./me_files/reac2_symm.me',status='unknown')
         WORD_SYM='SYMMETRY_2'
         do while (WORD.NE.'SYMMETRY_2')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_SYM = 'SYMMETRY'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/symm_reac2_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      else if (ispecies.eq.3) then
         open (unit=15,file='./data/prod1.dat',status='old')
         open (unit=16,file='./output/prod1_sym.out',status='unknown')
         open (unit=17,file='./geoms/prod1_l1.xyz',status='unknown')
         open (unit=18,file='./me_files/prod1_symm.me',status='unknown')
         WORD_SYM='SYMMETRY_3'
         do while (WORD.NE.'SYMMETRY_3')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_SYM = 'SYMMETRY'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/symm_prod1_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      else if (ispecies.eq.4) then
         open (unit=15,file='./data/prod2.dat',status='old')
         open (unit=16,file='./output/prod2_sym.out',status='unknown')
         open (unit=17,file='./geoms/prod2_l1.xyz',status='unknown')
         open (unit=18,file='./me_files/prod2_symm.me',status='unknown')
         WORD_SYM='SYMMETRY_4'
         do while (WORD.NE.'SYMMETRY_4')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_SYM = 'SYMMETRY'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/symm_prod2_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      else if (ispecies.eq.5.or.ispecies.eq.51) then
         open (unit=15,file='./data/wellr.dat',status='old')
         open (unit=16,file='./output/wellr_sym.out',status='unknown')
         open (unit=17,file='./geoms/wellr_l1.xyz',status='unknown')
         open (unit=18,file='./me_files/wellr_symm.me',status='unknown')
         WORD_SYM='SYMMETRY_5'
         do while (WORD.NE.'SYMMETRY_5')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_SYM = 'SYMMETRY'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/symm_wellr_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      else if (ispecies.eq.6.or.ispecies.eq.61) then
         open (unit=15,file='./data/wellp.dat',status='old')
         open (unit=16,file='./output/wellp_sym.out',status='unknown')
         open (unit=17,file='./geoms/wellp_l1.xyz',status='unknown')
         open (unit=18,file='./me_files/wellp_symm.me',status='unknown')
         WORD_SYM='SYMMETRY_6'
         do while (WORD.NE.'SYMMETRY_6')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_SYM = 'SYMMETRY'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/symm_wellp_molpro.dat 
     $                    level0_molpro.dat'            
         endif
      else if (ispecies.eq.0) then
         open (unit=15,file='./data/ts.dat',status='unknown')
         open (unit=16,file='./output/ts_sym.out',status='unknown')
         open (unit=17,file='./geoms/tsgta_l1.xyz',status='unknown')
         open (unit=18,file='./me_files/ts_symm.me',status='unknown')
         WORD_SYM='SYMMETRY_TS'
         do while (WORD.NE.'SYMMETRY_TS')
            call LineRead(21)
            if (WORD.EQ.'END') then
               WORD_SYM = 'SYMMETRY'
               go to  900
            endif
         enddo
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/symm_ts_molpro.dat 
     $ level0_molpro.dat'            
         endif
      else
         write(16,*)'ispecies not defined'
         write(16,*)'in sub symmetry'
         write(16,*)'change option and restart'
         stop
      endif

  900 continue
      rewind(21)

      ifreq=0
      noptg = 0
      ires=0
      isite=0
      imhind=0

cc first check if a multirotor scan has been performed
      call LineRead(0)
      do while (WORD.NE.'END')
         call LineRead (15)
         if (WORD.EQ.'NHIND1D') then
            read(15,*)nmhind
            if(nmhind.gt.0)imhind=1
         endif
      enddo
      rewind(15)
      call LineRead(0)
      do while (WORD.NE.'END')
         call LineRead (15)
         if (WORD.EQ.'NHIND2D') then
            read(15,*)nmhind
            if(nmhind.gt.0)imhind=1
         endif
      enddo
      rewind(15)
      call LineRead(0)
      do while (WORD.NE.'END')
         call LineRead (15)
         if (WORD.EQ.'NHIND3D') then
            read(15,*)nmhind
            if(nmhind.gt.0)imhind=1
         endif
      enddo
      rewind(15)

cc now check if the species is monoatomic  
cc in which case quit this subroutine

      if(ispecies.eq.1.or.ispecies.eq.2.or.
     $   ispecies.eq.3.or.ispecies.eq.4) then
         do while (WORD.NE.'NATOM')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (96,*) 'in symmetry routine'
               write (96,*) 'natom must be defined'
               stop
            endif
         enddo
         read (15,*) natom1,natomt1
         rewind(15)
         if(natom1.eq.1) goto 999
      endif

      do while (WORD.NE.'NHIND')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (16,*) 'hind rotors must be defined'
            stop
         endif
      enddo
c      call LineRead (15)
      read (15,*) nhind
      read (15,*)
      ntau_fr = nhind
      do ihind = 1 , nhind
         read (15,*) hindlab(ihind),dhindmn(ihind),dhindmx(ihind),
     +              nhindsteps(ihind),nhindsymm(ihind)
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) hindlab(ihind)
         rewind (99)
         call LineRead (99)
         hindlab(ihind)=WORD
         close (unit=99,status='keep')
      enddo
  
      do while (WORD.NE.'CHARGE')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (16,*) 'charge and spin must be defined'
            stop
         endif
      enddo
      read (15,*) icharge,ispin
      rewind 15

      if(ispecies.eq.0)then
         open (unit=25,file='./data/reac1.dat',status='old')
         do while (WORD.NE.'NATOM')
            call LineRead (25)
            if (WORD.EQ.'END') then
               write (96,*) 'natom in reac1 must be defined'
               stop
            endif
         enddo
         read (25,*) natom1,natomt1
         close (25)

         do while (WORD.NE.'ISITE')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (96,*) 'reaction site must be defined'
               stop
            endif
         enddo
         read (15,*) isite,jsite,ksite
         rewind 15
      endif

      close(15)
      
c read level of theory for symmetry determination

      do while (WORD.NE.WORD_SYM)
         call LineRead (21)
         if (WORD.EQ.'END') then
            if(nhind.ne.0) then
               write (16,*) 'symmetry level of theory must be defined'
               write (16,*) 'in file theory.dat'
               stop
            else 
               goto 9001
            endif
         endif
      enddo
      if(word2.eq.'G09')then
         ilev1code=1
         read (21,'(A70)') comline1
         read (21,'(A70)') comline2
      else if(word2.eq.'MOLPRO')then
         ilev1code=2
         call commrun(commandcopy)
      else
         write (16,*) 'symmetry level calculations either g09 or molpro'
         write (16,*) 'in file theory.dat'
      endif

      if (idebug.ge.2) write (16,*) ' comline test',comline1,comline2

 9001 continue
      close(21)


c now read xyz input

      read (17,*)natom
      read (17,*)enl1_ref
      do iatom = 1 , natom
         read(17,*)cjunk,coox,cooy,cooz
         open (unit=99,status='unknown')
         write(99,1201)cjunk,coox,cooy,cooz
         rewind(99)
         read (99,'(A60)') atomlabel(iatom)
         close(99)
      enddo
      if (idebug.ge.2) write (6,*) ' finished input section'
      close(17)

      open (unit=15,file='symm.inp',status='unknown')
      write(15,*)'Angstrom'
      write(15,*)'Tolerance 0.07'
      write(15,*)'Geometry ',natom
      do iatom = 1 , natom
         write (15,'(A60)') atomlabel(iatom)
      enddo
      close(15)

      command1='symmetry_number symm.inp > symm.out'
      call commrun(command1)

      open (unit=99,file='symm.out',status='unknown')
      read(99,*)nrotsymm,noptiso
      write(16,*)' Rotational Symmetry Number ',nrotsymm
      write(16,*)' Number of optical isomers ',noptiso
      close(99)

cc read potential for the ihind hind rotor
cc and identify minima 

      write(16,*)
      write(16,*)'Analysis for minima along hindered rotor path'
      write(16,*)

cc if we have a multiD hind rotor we skip this section
cc as it is assumed that optical isomers multiplicity is capture
cc by them md analysis
      if(imhind.eq.1)nhind=0

      if(nhind.ne.0) then
         do ihind=1,nhind
            write(16,*)'Analysis of hind rotor ',ihind
            do iscan=1,nhindsteps(ihind)
               open (unit=99,status='unknown')
               write(99,1200)ispecies,ihind,iscan
               rewind (99)
               read (99,'(A40)') filename
               close (99)
               open(unit=99,file=filename,status='unknown')
               read(99,*)
               read(99,*)energy(iscan,ihind)
               write(*,*)'energy is ',energy(iscan,ihind)
               close(99)
            enddo
            nm=1
            energymin(1,ihind)=energy(1,ihind)
            istepmin(nm,ihind)=1
            do iscan=1,nhindsteps(ihind)
               if(iscan.gt.2.and.iscan.le.nhindsteps(ihind)-1)then
                  if(energy(iscan,ihind).lt.energy(iscan-1,ihind).and.
     $                energy(iscan,ihind).lt.energy(iscan+1,ihind))then
                     nm=nm+1
                     energymin(nm,ihind)=energy(iscan,ihind)
                     istepmin(nm,ihind)=iscan
                  endif
               endif
            enddo
            nmin(ihind)=nm
            write(16,*)'nmin is ',nmin(ihind)
            do i=1,nmin(ihind)
               write(16,*)'energy is ',energymin(i,ihind)
               write(16,*)'step is ',istepmin(i,ihind)
            enddo
         enddo
      endif

cc 
cc determine energy and structure for each minima if nminima > 1

      do ihind=1,nhind
         if (nmin(ihind).gt.1)then
            do imin=1,nmin(ihind)
               open (unit=99,status='unknown')
               write(99,1200)ispecies,ihind,istepmin(imin,ihind)
               rewind (99)
               read (99,'(A40)') filename
               close (99)
               open(unit=99,file=filename,status='unknown')
               read(99,*)
               read(99,*)energy(iscan,ihind)
               write(*,*)'energy is ',energy(iscan,ihind)
               do iatom=1,natom
                  read(99,*)cjunk,coox,cooy,cooz
                  open (unit=199,status='unknown')
                  if(ifrozrts.eq.1.and.ispecies.eq.0)then
                     if(iatom.eq.isite.or.iatom.eq.(natom1+1))then
                        write(199,1204)cjunk,coox,cooy,cooz
                        rewind(199)
                        read (199,'(A60)') atomlabel(iatom)
                        close(199)
                     else
                        write(199,1205)cjunk,coox,cooy,cooz
                        rewind(199)
                        read (199,'(A60)') atomlabel(iatom)
                        close(199)
                     endif
                  else
                     write(199,1201)cjunk,coox,cooy,cooz
                     rewind(199)
                     read (199,'(A60)') atomlabel(iatom)
                     close(199)
                  endif
               enddo
               close(99)
               natomt=natom
               ircons=0
cc no rotational minima to investigate for linear molecules
               ilin=0
               ntau=0
               do i=1,natom
                  xint(i)=0
               enddo
               ncoord=3*natom-6
               do i=1,ncoord
                  intcoor(i)=''
               enddo
               ixyz=1
               ismp=imin
               ired=0
               ilev=0

               if(ilev1code.eq.1) then
                  call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $                 coord,vtot_0,vtot,freq,ifreq,ilin,ismp,comline1,
     $                 comline2,icharge,ispin,ircons,
     $                 atomlabel,intcoor,bislab,tauopt,xint,abcrot,
     $                 ires,ixyz,ired)

               else if (ilev1code.eq.2) then

                  numproc=numprochl         
                  call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $                 coord,vtot,freq,ifreq,ilin,ismp,
     $                 icharge,ispin,ircons,
     $                 atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,
     $                 ixyz,ilev,ispecies)
                  numproc=numprocll
               endif
            
               energyopt(imin,ihind)=vtot
               write(16,*)'En for hind rotor',ihind,' is ',imin,
     $              energyopt(imin,ihind)
               write(16,*)'Geometry '
               do iatom=1,natom
                  write(16,1202)iatom,(coord(iatom,idim),idim=1,3)
               enddo
            enddo
         endif
      enddo

      do ihind=1,nhind
         nrotiso(ihind)=1
         do imin=2,nmin(ihind)
            endiff=energyopt(1,ihind)-energyopt(imin,ihind)
            if(abs(endiff).lt.enthreshold)then
               nrotiso(ihind)=nrotiso(ihind)+1
            endif
         enddo
         write(16,*)'Opt Rot isomers of rotor ',ihind,nrotiso(ihind)
      enddo
      
      optisofin=noptiso
      do ihind=1,nhind
         optisofin=optisofin/nrotiso(ihind)
      enddo

cc when using a md rotor model the optical symm number is set to 1
cc as it is assumed that optical symmetry is captured by the HR analysis

      if(imhind.eq.1) then
         optisofin=1
         write(16,*)' setting optical symm to 1 as MD rotor is active  '
      endif

      globalsymm=float(nrotsymm)/optisofin

      write(16,*)' Final Results of symmetry analysis '
      write(16,*)' Rotational Symmetry Number ',nrotsymm
      write(16,*)' Number of optical isomers ',optisofin
      write(16,*)' Rot Symm/ Opt iso number ',globalsymm

      if(imhind.eq.0) then
         write(18,*)'        Core    RigidRotor'
         write(18,*)'               SymmetryFactor ',globalsymm
         write(18,*)'        End'
      else
         write(18,*)'   Core MultiRotor'
         write(18,*)'   SymmetryFactor ',globalsymm
      endif

      close(16)
      close(18)


      if (ispecies.eq.1) then
         if(imhind.eq.0)then
            command1='head -n -3  me_files/reac1_1dge.me
     $> temp.me '
            call commrun(command1)
            command1='cat temp.me  me_files/reac1_symm.me
     $> me_files/reac1_ge.me '
            call commrun(command1)
         else
            command1='tail -n +3  me_files/reac1_hr.me
     $> temp.me '
            call commrun(command1)
            command1='cat me_files/reac1_symm.me temp.me
     $> me_files/reac1_hr.me '
            call commrun(command1)            
         endif
      endif

      if (ispecies.eq.2) then
         if(imhind.eq.0)then
            command1='head -n -3  me_files/reac2_1dge.me
     $> temp.me '
            call commrun(command1)
            command1='cat temp.me  me_files/reac2_symm.me
     $> me_files/reac2_ge.me '
            call commrun(command1)
         else
            command1='tail -n +3  me_files/reac2_hr.me
     $> temp.me '
            call commrun(command1)
            command1='cat me_files/reac2_symm.me temp.me
     $> me_files/reac2_hr.me '
            call commrun(command1)            
         endif
      endif
      if (ispecies.eq.3) then
         if(imhind.eq.0)then
            command1='head -n -3  me_files/prod1_1dge.me
     $> temp.me '
            call commrun(command1)
            command1='cat temp.me  me_files/prod1_symm.me
     $> me_files/prod1_ge.me '
            call commrun(command1)
         else
            command1='tail -n +3  me_files/prod1_hr.me
     $> temp.me '
            call commrun(command1)
            command1='cat me_files/prod1_symm.me temp.me
     $> me_files/prod1_hr.me '
            call commrun(command1)            
         endif
      endif
      if (ispecies.eq.4) then
         if(imhind.eq.0)then
            command1='head -n -3  me_files/prod2_1dge.me
     $> temp.me '
            call commrun(command1)
            command1='cat temp.me  me_files/prod2_symm.me
     $> me_files/prod2_ge.me '
            call commrun(command1)
         else
            command1='tail -n +3  me_files/prod2_hr.me
     $> temp.me '
            call commrun(command1)
            command1='cat me_files/prod2_symm.me temp.me
     $> me_files/prod2_hr.me '
            call commrun(command1)            
         endif
      endif
      if (ispecies.eq.5) then
         if(imhind.eq.0)then
            command1='head -n -3  me_files/wellr_1dge.me
     $> temp.me '
            call commrun(command1)
            command1='cat temp.me  me_files/wellr_symm.me
     $> me_files/wellr_ge.me '
            call commrun(command1)
         else
            command1='tail -n +3  me_files/wellr_hr.me
     $> temp.me '
            call commrun(command1)
            command1='cat me_files/wellr_symm.me temp.me
     $> me_files/wellr_hr.me '
            call commrun(command1)            
         endif
      endif
      if (ispecies.eq.6) then
         if(imhind.eq.0)then
            command1='head -n -3  me_files/wellp_1dge.me
     $> temp.me '
            call commrun(command1)
            command1='cat temp.me  me_files/wellp_symm.me
     $> me_files/wellp_ge.me '
            call commrun(command1)
         else
            command1='tail -n +3  me_files/wellp_hr.me
     $> temp.me '
            call commrun(command1)
            command1='cat me_files/wellp_symm.me temp.me
     $> me_files/wellp_hr.me '
            call commrun(command1)            
         endif
      endif

      if (ispecies.eq.0) then
         if(imhind.eq.0)then
            command1='head -n -3  me_files/ts_1dge.me
     $> temp.me '
            call commrun(command1)
            command1='cat temp.me  me_files/ts_symm.me
     $> me_files/ts_ge.me '
            call commrun(command1)
         else
            command1='tail -n +3  me_files/ts_hr.me
     $> temp.me '
            call commrun(command1)
            command1='cat me_files/ts_symm.me temp.me
     $> me_files/ts_hr.me '
            call commrun(command1)            
         endif
      endif


 1200 format ("./hr_geoms/geom_isp"I0.2"_hr"I0.2"_hpt"I0.2".xyz       ")
 1201 format(1X,A2,1X,F9.5,1X,F9.5,1X,F9.5,40x' ')
 1202 format(1X,I2,1X,F9.5,1X,F9.5,1X,F9.5,1X)
 1203 format(1X,I2,1X,I2,1X,F9.5,1X,'F',40x' ')
 1204 format(1X,A2,1X,'-1',1X,F9.5,1X,F9.5,1X,F9.5,40x' ')
 1205 format(1X,A2,1X,'0',1X,F9.5,1X,F9.5,1X,F9.5,40x' ')

 999  continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine irc

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension dhindmn(nhindmx),dhindmx(nhindmx),
     $ freq(nmdmx),nhindsteps(nhindmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim) 
      dimension freqtot(3*natommx),freqtot1pro(3*natommx)
     $ ,imatch(3*natommx)
      dimension vref(noptmx),tauo(ntaumx,noptmx),
     $ xinto(3*natommx,noptmx),dist(noptmx)
      dimension tauopt(ntaumx)
      dimension tau(ntaumx)
      dimension igrouptot(nhindmx),ipivotA(nhindmx),ipivotB(nhindmx)
      dimension freqsubs(nhindmx,nircmx)
      dimension ngroup(nhindmx,natommx)
      dimension rc_ene(nircmx)
      dimension zpe_irc(nircmx)
      dimension rc_coord(nircmx)
      dimension rc_ene_kcal(nircmx)
      dimension rc_ene_hl(nircmx)
      dimension freqproj(3*natommx,nircmx)
      dimension freqintproj(3*natommx,nircmx)
      dimension freqRTproj(3*natommx,nircmx)
      dimension rotpot(noptmx)
      dimension gelec(nelecmx),eelec(nelecmx)
      dimension rotpot_ir(noptmx,nircmx,nhindmx)
      dimension numpot_ir(nhindmx)
      dimension nhrcc_points(nircmx)
      dimension ibconn(natommx),iaconn(natommx),idconn(natommx)

      character*20 hr_rescale,hr_rescale2f,hr_rescale2b,hr_rescale3
      character*70 comline1,comline2,comline3,comline4
c      character*70 comline3,comline4
      character*60 atomlabel(natommx)
      character*30 intcoor(3*natommx)
      character*30 intcoori(3*natommx)
      character*20 bislab(ntaumx)
      character*20 hindlab(nhindmx)
      character*30 cjunk
      character*60 rotline
      character*70 rot2dline
      character*180 command1
      character*80 atgeom_me(natommx,nircmx)
      character*100 atgeom_pr(natommx,nircmx)
cc dimension of this matrix should be double checked
      character*80 force_con(3*natommx*3*natommx/10,nircmx)
      character*80 grad(natommx,nircmx)
      character*100 gradts(natommx)
cc      character*20 rc_coord(nircmx)
      character*20 step(nircmx)
      character*20 filename
      character*2 aname
      character*30 gmem
      character*20 bname(natommx),anname(natommx),dname(natommx)
     $ ,atname(natommx),bconnt(natommx),aconnt(natommx),dconnt(natommx)


      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

      ires=0
      ionlyfor=0 
cc we assume the TS is not linear         
      ilin=0

c input data

      open (unit=96,file='./output/irc.out',status='unknown')
      open (unit=17,file='./output/ts_opt.out',status='unknown')
      open (unit=21,file='./data/theory.dat',status='unknown')
      open (unit=27,file='./me_files/ts_ge.me',status='unknown')

      call LineRead (0)

c read theory
      do while (WORD.NE.'IRC')
         call LineRead (21)
         if (WORD.EQ.'END') then
            write (96,*) 'IRC must be defined'
            write (96,*) 'in file theory.dat'
            stop
         endif
      enddo
      if(word2.eq.'G09') then
         ilev0code=1
         read (21,'(A70)') comline1
         read (21,'(A70)') comline2
         read (21,'(A70)') comline3
         read (21,'(A70)') comline4
      else
         write (96,*) 'level0 of theory must be defined'
         write (96,*) 'in theory.dat'
         stop
      endif
      iskiptheo=0
      iskiphl=0
      if(word3.eq.'SKIP') iskiptheo=1 
      if(irecov.eq.1) iskiptheo=1
      if(word3.eq.'SKIPALL')then 
         iskiptheo=1 
         iskiphl=1
      endif
c      write(*,*)'iskiptheo is ',iskiptheo
c      write(*,*)'irecover is ',irecov
c      stop

      if(word3.eq.'ONLYFOR'.or.word4.eq.'ONLYFOR') ionlyfor=1 

      call LineRead(21)
      hr_rescale=word
      hr_rescale2f=word2
      hr_rescale2b=word3
      hr_rescale3=word4
      close (21)

cc read structural parameters for TS


cc get data from react1 file

      open (unit=25,file='./data/reac1.dat',status='old')
      do while (WORD.NE.'NATOM')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (96,*) 'natom in reac1 must be defined'
            stop
         endif
      enddo
      read (25,*) natom1,natomt1
      close (25)

cc get data from react2 file
      if(iadd.eq.1.or.iabs.eq.1)then
         call LineRead (0)

         open (unit=25,file='./data/reac2.dat',status='old')

c     natomt is to account for dummy atoms
         do while (WORD.NE.'NATOM')
            call LineRead (25)
            if (WORD.EQ.'END') then
               write (96,*) 'natom must be defined'
               stop
            endif
         enddo
         read (25,*) natom2,natomt2
         close (25)
      endif

      if(iabs.eq.1.or.iadd.eq.1)then
         natom = natom1+natom2
      else if (iiso.eq.1.or.ibeta.eq.1) then
         natom = natom1
      endif

      if (iadd.eq.1) natomt = natomt1+natomt2
      if (iabs.eq.1) natomt = natomt1+natomt2+1
      if (iiso.eq.1) natomt = natomt1
      if (ibeta.eq.1) natomt = natomt1

      open (unit=15,file='./data/ts.dat',status='unknown')
      do while (WORD.NE.'NHIND')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (96,*) 'hind rotors must be defined'
            stop
         endif
      enddo
      read (15,*) nhind
      rewind 15

c determine number of reacting atom for beta-scission reactions
      if(ibeta.eq.1)then
         do while (WORD.NE.'RMAX1')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (96,*) 'grid coords must be defined'
               write (96,*) 'for isomerization of betascission reacts'
               write (96,*) 'for IRC scans'
               stop
            endif
         enddo
         read (15,*) rmax1,rstp1,rmax2,rstp2,ireact
         rewind(15)
      endif

      if(iiso.eq.1)then
         do while (WORD.NE.'RMIN1')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (96,*) 'grid coords must be defined'
               write (96,*) 'for isomerization of betascission reacts'
               write (96,*) 'for IRC scans'
               stop
            endif
         enddo
         read (15,*) rmax1,rstp1,rmax2,rstp2,ireact
         rewind(15)
      endif

cc read 2d hind rotor list
      do while (WORD.NE.'NHIND2D')
         call LineRead (15)
         if (WORD.EQ.'END') then
            nhind2d=0
            goto 1112
cc the nhind2d section must not necessarily be there
         endif
      enddo
c      call LineRead (15)
      read (15,*) nhind2d
 1112 continue  
      rewind 15

      do while (WORD.NE.'NHIND3D')
         call LineRead (15)
         if (WORD.EQ.'END') then
            nhind3d=0
            goto 1113
cc the nhind3d section must not necessarily be there
         endif
      enddo
c      call LineRead (15)
      read (15,*) nhind3d
 1113 continue  
      rewind 15

      nhindmd=nhind3d+nhind2d

      do while (WORD.NE.'ISITE')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (96,*) 'reaction site must be defined'
            stop
         endif
      enddo
      read (15,*) isite,jsite,ksite
      rewind 15

      do while (WORD.NE.'NTAU')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (96,*) 'sampling coords of reac2 must be defined'
            stop
         endif
      enddo
      read (15,*) ntau
      rewind(15)

      do while (WORD.NE.'CHARGE')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (96,*) 'charge and spin must be defined'
            stop
         endif
      enddo
      read (15,*) icharge,ispin
      rewind(15)


      do while (WORD.NE.'NELEC')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (96,*) 'electronic states must be defined'
            stop
         endif
      enddo
      read (15,*) nelec
      do ielec = 1 , nelec
         read (15,*) eelec(ielec),gelec(ielec)
      enddo
      rewind(15)

      close(15)


cc read symmetry from me file
cc necessary so that it includes the initialization from the symm sub

      if(nhindmd.eq.0) then
         do while (WORD.NE.'SYMMETRYFACTOR')
            call LineRead (27)
            if (WORD.EQ.'END') then
               write (96,*) 'symmetry factor must be defined'
               stop
            endif
         enddo
         open (unit=99,status='unknown')
         write (99,*)word2
         rewind (99)
         read (99,*) symf
         close (99)
      else
         symf=1
      endif
      close(27)

c  reading z-mat
      read (17,*)
      if (idebug.ge.2) write (6,*) ' reading z-mat input'
      do iatom = 1 , natomt
         read (17,'(A60)') atomlabel(iatom)
      enddo

cc read coordinate names

      ncoord = 3*natom-6

      do iint = 1 , ncoord
         read (17,*) intcoor(iint),xint(iint)
c         write (6,*) 'intcoor, xint ',intcoori(iint),xinti(iint)
      enddo

      close (unit=17,status='keep')


      if (idebug.ge.2) write (16,*) 'past z-matrix'



cc now perform IRC scan

      ircons=0
      ntau=0
      ifreq=0 
      ixyz=0
      ired=0
      ires=0

      if(ilev0code.eq.1) then
         write (16,*) 'starting g09pot'

         if(iskiptheo.ne.1) then
            call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $   coord,vtot_0,vtot,freq,ifreq,ilin,ismp,comline1,
     $    comline2,icharge,ispin,ircons,
     $    atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ired)

cc save the output

            command1='cp -f geom.log ./irc_files/irc_g09f.log'
            call commrun(command1)         

cc now do the backward step
            if(ionlyfor.ne.1) then
               comline1=comline3
               comline2=comline4

               call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $   coord,vtot_0,vtot,freq,ifreq,ilin,ismp,comline1,
     $    comline2,icharge,ispin,ircons,
     $    atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ired)

               command1='cp -f geom.log ./irc_files/irc_g09b.log'
               call commrun(command1)    
            endif
         endif


         call read_g09_ircout(force_con,natom,iread,numpointsf,
     $  numpointsb,
     $  atgeom_me,ifcread,rc_ene,rc_coord,grad)

c         write(*,*)'ok up to here cc2'
c         stop

      endif
      
c 
c  from the irc output we need to fill in the atgeom_me,atgeom_pr
c ,grad,force_con,rc_ene, and rc_coor vectors. 
c  the numpointsf, nummpoitsb and ifcread parameters are inizialited as well 
c

      if(ionlyfor.ne.1)then
         numpointstot=numpointsf+numpointsb+1
      else
         numpointstot=numpointsf+1
         numpointsb=0
      endif

cc first create the at_geom vector as necessary for freq projection

      do j=1,numpointstot
         if(j.ne.numpointsf+1)then
            do  iatom = 1, natom
               open (unit=99, status='unknown')
               write(99,*)grad(iatom,j)
               rewind(99)
               read(99,*)ind1,ind2
               ind3=0
               rewind(99)
               write(99,*)atgeom_me(iatom,j)
               rewind(99)
               read(99,*)cjunk,coox,cooy,cooz
               rewind(99)
               write(99,1199)ind1,ind2,ind3,coox,cooy,cooz
               rewind(99)
               read (99,'(A100)') atgeom_pr(iatom,j)
               close(99)
            enddo
         else
            open(unit=108,file='geoms/tsgta_l1.xyz',status='unknown')
            read(108,*)
            read(108,*)
            do iatom = 1, natom
               read(108,*)cjunk,coox,cooy,cooz
               open (unit=99, status='unknown')
               write(99,*)grad(iatom,numpointsf)
               rewind(99)
               read(99,*)ind1,ind2
               ind3=0
               close(99)
               open (unit=99, status='unknown')
               rewind(99)
               write(99,1199)ind1,ind2,ind3,coox,cooy,cooz
               rewind(99)
               read (99,'(A100)') atgeom_pr(iatom,numpointsf+1)
               close(99)
            enddo
            close(108)
         endif
      enddo
 1199 format (I3,1x,I3,1x,I3,1x,F9.5,1x,F9.5,1x,F9.5)


cc now write coordinates and energies file
c      
      open (unit=333,file='RPHt_coord_en.dat',status='unknown')

      write(333,*)'Point Coordinate Energy Bond1 Bond2'

      call read_zmat(atomlabel,natom,natomt,intcoor,bislab,ibconn,
     $ iaconn,idconn,bname,anname,dname,atname,idummy,isited,jsited,
     $ ksited,aconnt,bconnt,dconnt)


      do i=1,numpointsf+numpointsb+1
            open (unit=99, status='unknown')
            write(99,*)atgeom_me(isited,i)
            write(99,*)atgeom_me(jsited,i)
            if(iadd.eq.1.or.iabs.eq.1)then
               write(99,*)atgeom_me(natom1+1,i)
            else if (iiso.eq.1.or.ibeta.eq.1) then
               write(99,*)'X', 0.0, 0.0, 0.0
            endif
            rewind(99)
            read(99,*)cjunk,atcentx,atcenty,atcentz
            read(99,*)cjunk,atreax,atreay,atreaz
            read(99,*)cjunk,atprodx,atprody,atprodz
            close(99)
            dist_atc_rea=sqrt((atcentx-atreax)**2.
     $       +(atcenty-atreay)**2.
     $       +(atcentz-atreaz)**2.0)
            dist_atc_pro=sqrt((atcentx-atprodx)**2.0
     $       +(atcenty-atprody)**2.0
     $       +(atcentz-atprodz)**2.0)
 

            write(*,*)'energy = ',rc_ene(i)
         if(iadd.eq.1.or.iabs.eq.1)then
            write(333,1333)i,rc_coord(i),rc_ene(i),dist_atc_rea,
     $                  dist_atc_pro
            else if (iiso.eq.1.or.ibeta.eq.1) then
               write(333,1333)i,rc_coord(i),rc_ene(i),dist_atc_rea,
     $                  0.0
            endif
      enddo
c      do i=1,numpointsf+numpointsb+1
c         write(*,*)'FC = ',force_con(2,i)
c      enddo
      close(333)
 
 1333 format (I3,1X,F8.5,1X,F10.7,1X,F8.5,1X,F8.5)
 1032 format (A120)

cc
cc set to 0 the gradient at the saddle point

      do j=1,natom
         gradts(j)=grad(j,numpointsf+1)
         grad(j,numpointsf+1)=' 1 1 0. 0. 0.'
      enddo


cc at this point, all vectors are filled in 
cc now we can compute the projected frequencies

      open (unit=16,file='freqout.dat'
     $         ,status='unknown')
      open (unit=17,file='fresub.dat'
     $         ,status='unknown')
      open (unit=18,file='freqRTout.dat'
     $         ,status='unknown')

cc write input for multidimensional tunneling calculations

      open (unit=333,file='RPHt_all_data.dat',status='unknown')
      write(333,*)'Number_of_Atoms:            ',natom
      write(333,*)'Act_energy(kcal/mol):       0. '
      write(333,*)'Initial_Temperature:        50'
      write(333,*)'Temperature_steps:          40'
      write(333,*)'Temperature_increment:      40'
      write(333,*)'Delta_Energy_rea:           0.'
      write(333,*)'Delta_Energy_pro:           0.'
      write(333,*)'Maxstep:                   ',numpointstot
      write(333,*)'Npointsint:                 5 '
      write(333,*)'Maxtdev:                    0.5'
      write(333,*)'Rearrange(1=yes,0=no)       1'
      write(333,*)'SaddlePoint                ',numpointsf
         if(intfreq.eq.0)then
            write(333,*)'internalcoord(1=yes)     0'
         else if (intfreq.eq.1)then
            write(333,*)'internalcoord(1=yes)     1'            
         endif
      write(333,*)'isct_vtst(1=vtst_sct,0=sct) 1'
      write(333,*)'zerocurvature(1)            0'
      write(333,*)'reduced_mass                1.0'
      write(333,*)'minimum_frequency            50'
      write(333,*)'anim_freq(if_Maxstep=1)          2'
      write(333,*)'onlyrotors(0=yes,1=no)        1'
      write(333,*)'proj_rea_coo(0=yes(def),1=no) ',iprojrcoo
      write(333,*)'numrotors                     ',0

cc here starts the projection cycle over the total number of IRC points

      
      do inumpoints = 1, numpointstot
         open (unit=133,file='RPHt_input_data.dat',status='unknown')
c         open (unit=134,file='./data/hind_rot_head.dat',
c     +         status='unknown')

         write(133,*)'Number_of_Atoms: ',natom
         write(133,*)'Act_energy(kcal/mol):       0. '
         write(133,*)'Initial_Temperature:        200'
         write(133,*)'Temperature_steps:          40'
         write(133,*)'Temperature_increment:      40'
         write(133,*)'Delta_Energy_rea:           0.'
         write(133,*)'Delta_Energy_pro:           0.'
         write(133,*)'Maxstep:                    1'
         write(133,*)'Npointsint:                 5 '
         write(133,*)'Maxtdev:                    0.5'
         write(133,*)'Rearrange(1=yes,0=no)       1'
         write(133,*)'SaddlePoint                 1'
         if(intfreq.eq.0)then
            write(133,*)'internalcoord(1=yes)     0'
         else if (intfreq.eq.1)then
            write(133,*)'internalcoord(1=yes)     1'            
         endif
         write(133,*)'isct_vtst(1=vtst_sct,0=sct) 1'
         write(133,*)'zerocurvature(1)            0'
         write(133,*)'reduced_mass                1.0'
         write(133,*)'minimum_frequency            50'
         write(133,*)'anim_freq(if_Maxstep=1)       2'
         write(133,*)'onlyrotors(0=yes,1=no)        0'


c         do iatom = 1, 17
c            read (134,'(A60)') atomlabel(iatom) 
c            write (133,'(A60)') atomlabel(iatom) 
c         enddo
c         close(134)
         
         if(nhind.ne.0) then
            open (unit=15,file='./output/hrdata4proj_ts.dat'
     $           ,status='unknown')
            read (15,*)cjunk
            if(inumpoints.eq.numpointsf+1)then
               write (133,*)'proj_rea_coo(0=yes(def),1=no) ',1
            else
               write (133,*)'proj_rea_coo(0=yes(def),1=no) ',iprojrcoo
            endif
            read (15,*)cjunk,nhind
            write (133,*)cjunk,nhind

            do ir=1,nhind
               read (15,*)cjunk,ipivotA(ir)
               write (133,*)cjunk,ipivotA(ir)
               read (15,*)cjunk,ipivotB(ir)
               write (133,*)cjunk,ipivotB(ir)
               read (15,*)cjunk,igrouptot(ir)
               write (133,*)cjunk,igrouptot(ir)
               read (15,*)cjunk,(ngroup(ir,igr),igr=1,
     $              igrouptot(ir))
               write (133,*)cjunk,(ngroup(ir,igr),igr=1,
     $              igrouptot(ir))

            enddo
            close(15)
         else
            if(inumpoints.eq.numpointsf+1)then
               write (133,*)'proj_rea_coo(0=yes(def),1=no) ',1
            else
               write (133,*)'proj_rea_coo(0=yes(def),1=no) ',iprojrcoo
            endif
            write (133,*)'numrotors ',0
         endif

         write (133,*) 'Step', inumpoints
         write (333,*) 'Step', inumpoints
         write (133,*) 'geometry'
         write (333,*) 'geometry'
    
         do iatom = 1, natom
            write (133,'(A70)') atgeom_pr(iatom,inumpoints)
            write (333,'(A70)') atgeom_pr(iatom,inumpoints)
         enddo
         
         write (133,*) 'gradient'       
         write (333,*) 'gradient'       
         do iatom = 1, natom
            write (133,'(A70)') grad(iatom,inumpoints)
            if(inumpoints.eq.numpointsf+1)then
               write (333,'(A70)') gradts(iatom)
            else
               write (333,'(A70)') grad(iatom,inumpoints)
            endif
         enddo
        
         write (133,*) 'Hessian'
         write (333,*) 'Hessian'
         do inumlines = 1, ifcread
            write (133,'(A77)')  force_con(inumlines,inumpoints)
            write (333,'(A77)')  force_con(inumlines,inumpoints)
         enddo
          
         write (133,*) 'End'
       
         close(133)

         open (unit=108,status='unknown')
         write (108,1305)inumpoints
         rewind (108)
         read (108,1032) command1
         close (108)
 1305    format (" cp -f  RPHt_input_data.dat 
     +    ./irc_files/input"I0.3 )
         call commrun(command1)

         if(intfreq.eq.1)then
            open(unit=15,file='geom_bmat.xyz',status='unknown')
            write(15,*)natom
            write(15,*)'IRC geom num ',inumpoints
            open(unit=99,status='unknown')
            do j=1,natom
               write(99,*)atgeom_pr(j,inumpoints)
               rewind(99)
               read(99,*)ijunk,ijunk,ijunk,coox,cooy,cooz
               write (15,*)'at',coox,cooy,cooz
               rewind(99)
            enddo
            close(99)
            close(15)
            call bmatrix(0,1)
c            write(*,*)'ok up to here'
c            stop
         endif
         command1='RPHt.exe'
         call commrun(command1)

cc         now read projected frequencies

         open (unit=15,file='hrproj_freq.dat'
     $         ,status='unknown')
         
         nfreq=3*natom-nhind-1-6
         do j=1,nfreq
            read(15,*)freqproj(j,inumpoints)
            write(16,*)freqproj(j,inumpoints)
         enddo
         write(16,*)
         close(15)

         if(intfreq.eq.1)then
            open (unit=15,file='hrprojint_freq.dat'
     $         ,status='unknown')
            nfreq=3*natom-nhind-1-6
            do j=1,nfreq
               read(15,*)freqintproj(j,inumpoints),eigen
               write(16,*)freqintproj(j,inumpoints)
               if(eigen.lt.0)freqintproj(j,inumpoints)=0.
            enddo
            write(16,*)
            close(15)
         endif

cc         now read frequencies projected for rot trasl and RC but not HR

         open (unit=15,file='RTproj_freq.dat'
     $         ,status='unknown')
         
         nfreq=3*natom-1-6
         do j=1,nfreq
            read(15,*)freqRTproj(j,inumpoints)
            write(18,*)freqRTproj(j,inumpoints)
            if(j.ne.1)then
               if(freqRTproj(j,inumpoints).gt.
     $              freqRTproj(j-1,inumpoints))then
                  freqRTproj(j,inumpoints)=0.
               endif
            endif
         enddo
         write(18,*)
         close(15)


cc      now determine for each rotor which frequency has been projected

c    determine freq for one rotor at a time 

         do ir=1,nhind

            open (unit=133,file='RPHt_input_data.dat',status='unknown')
c            open (unit=134,file='./data/hind_rot_head.dat',
c     +         status='unknown')

            write(133,*)'Number_of_Atoms: ',natom
            write(133,*)'Act_energy(kcal/mol):       0. '
            write(133,*)'Initial_Temperature:        200'
            write(133,*)'Temperature_steps:          40'
            write(133,*)'Temperature_increment:      50'
            write(133,*)'Delta_Energy_rea:           0.'
            write(133,*)'Delta_Energy_pro:           0.'
            write(133,*)'Maxstep:                    1'
            write(133,*)'Npointsint:                 5 '
            write(133,*)'Maxtdev:                    0.5'
            write(133,*)'Rearrange(1=yes,0=no)       1'
            write(133,*)'SaddlePoint                 1'
            write(133,*)'ds(1=noexp,0=standard)      0'
            write(133,*)'isct_vtst(1=vtst_sct,0=sct) 1'
            write(133,*)'zerocurvature(1)            0'
            write(133,*)'reduced_mass                1.0'
            write(133,*)'minimum_frequency            50'
            write(133,*)'anim_freq(if_Maxstep=1)       2'
            write(133,*)'onlyrotors(0=yes,1=no)        0'


c            do iatom = 1, 17
c               read (134,'(A60)') atomlabel(iatom) 
c               write (133,'(A60)') atomlabel(iatom) 
c            enddo
c            close(134)


            write (133,*)'proj_rea_coo(0=yes(def),1=no) ',iprojrcoo
            write (133,*)'numrotors        ',1
            write (133,*)'pivotA           ',ipivotA(ir)
            write (133,*)'pivotB           ',ipivotB(ir)
            write (133,*)'atomsintopA      ',igrouptot(ir)
            write (133,*)'topAatoms        ',(ngroup(ir,igr),igr=1,
     $           igrouptot(ir))

            write (133,*) 'Step', inumpoints
            write (133,*) 'geometry'
    
            do iatom = 1, natom
               write (133,'(A70)') atgeom_pr(iatom,inumpoints)
            enddo
         
            write (133,*) 'gradient'       
            do iatom = 1, natom
               write (133,'(A70)') grad(iatom,inumpoints)
            enddo
        
            write (133,*) 'Hessian'
            do inumlines = 1, ifcread
               write (133,'(A77)')  force_con(inumlines,inumpoints)
            enddo
          
            write (133,*) 'End'
            close(133)
  
            command1='RPHt.exe'
            call commrun(command1)

c   now read 6D+1(RP) and 6D+1(RP)+1(HR) projected freqs     
            nfreq=3*natom-7
            nfreq1pro=3*natom-8

            open(unit=100,file='RTproj_freq.dat'
     $           ,status='unknown')
            do i=1,nfreq
               read(100,*) freqtot(i)
            enddo
            close(100)

            open(unit=100,file='hrproj_freq.dat'
     $           ,status='unknown')

            do i=1,nfreq1pro
               read(100,*) freqtot1pro(i)
            enddo
            close(100)

cc  now identify freq among freqtot that does not match fretotpro


            do i=1,nfreq
               imatch(i)=0
            enddo

            do j=1,nfreq1pro
               diff=0
               diffmin=10000.
               imin=0
               do i=1,nfreq
                  diff=abs(freqtot(i)-freqtot1pro(j))
                  if(diff.lt.diffmin.and.imatch(i).eq.0)then
                     diffmin=diff
                     imin=i
                  endif
               enddo
               imatch(imin)=1
c     write(*,*)'imatch is',imin
            enddo
         
            do i=1,nfreq
c     write(*,*)'imatch is',imatch(i)
               if(imatch(i).eq.0)then
                  freqsubs(ir,inumpoints)=freqtot(i)
c               write(*,*)'the not match freq is',freqtot(i)
               endif
            enddo
         enddo
         write(17,*)(freqsubs(ir,inumpoints),ir=1,nhind)

      enddo

cc here ends the cycle over the total number of IRC points

      close(16)
      close(17)
      close(18)
      close(333)


cc determine the scaling factors for the HR potentials

      if(hr_rescale.eq.'HRCC'.and.nhind.ne.0) then

         open (unit=99, status='unknown')
         write(99,*)hr_rescale2f
         rewind(99)
         read(99,*)num_hrcc_pointsf
         close(99)
         open (unit=99, status='unknown')
         write(99,*)hr_rescale2b
         rewind(99)
         read(99,*)num_hrcc_pointsb
         close(99)
  
cc initialize the vector with the sampled points

         istepf=(numpointsf)/(num_hrcc_pointsf+1)
         istepb=(numpointsb)/(num_hrcc_pointsb+1)
         if(istepf.eq.0.or.istepf.eq.1)then 
            istepf=1
            num_hrcc_pointsf=numpointsf-1
         endif
         if(istepb.eq.0)then
            istepb=1
            num_hrcc_pointsb=numpointsb-1
         endif
         if(istepb.eq.1)then
c            istepb=1
            num_hrcc_pointsb=numpointsb-1
         endif

c         num_hrcc_points_max=numpointsf+numpointsb+1
         num_hrcc_points_sampled=3+num_hrcc_pointsf+num_hrcc_pointsb
c         if(num_hrcc_points_sampled.gt.num_hrcc_points_max)then
c            num_hrcc_points_sampled=num_hrcc_points_max
c         endif

         num_hrcc_points=num_hrcc_points_sampled
         do j=1,num_hrcc_points
            if(j.eq.1) then 
               nhrcc_points(j)=1
            else if(j.eq.1+num_hrcc_pointsf+1) then
               nhrcc_points(j)=numpointsf+1
            else if(j.eq.3+num_hrcc_pointsf+num_hrcc_pointsb) then
               nhrcc_points(j)=numpointsf+numpointsb+1
            else if (j.lt.1+num_hrcc_pointsf+1) then
               nhrcc_points(j)=numpointsf-(num_hrcc_pointsf+1)*istepf
     $            +(j-1)*istepf+1
            else if (j.gt.1+num_hrcc_pointsf+1.and.
     $               j.lt.3+num_hrcc_pointsf+num_hrcc_pointsb) then
               nhrcc_points(j)=1+numpointsf+
     $                        (j-2-num_hrcc_pointsf)*istepb
            endif
         enddo

         do j=1,num_hrcc_points
            write(96,*)'HR sampled points is',nhrcc_points(j)
         enddo
c         stop

         if(iabs.eq.1.and.hr_rescale3.ne.'READ') then
cc calculate dist1 and dist2 of the forming and breaking bond
cc reacting atom: isite
cc isite is bound to ji
cc natom1+1 first atom on second reactant

            istep=(numpointsf+numpointsb)/(num_hrcc_points-1)
            index=0
            do j=1,num_hrcc_points
c               index=1+(j-1)*istep
               index=nhrcc_points(j)
               write(*,*)'index is',index
c            stop

               open (unit=99, status='unknown')
               write(99,*)atgeom_me(isite,index)
               write(99,*)atgeom_me(jsite,index)
               write(99,*)atgeom_me(natom1+1,index)
               rewind(99)
               read(99,*)cjunk,atcentx,atcenty,atcentz
               read(99,*)cjunk,atreax,atreay,atreaz
               read(99,*)cjunk,atprodx,atprody,atprodz
               close(99)
               dist_atc_rea=sqrt((atcentx-atreax)**2.
     $              +(atcenty-atreay)**2.
     $              +(atcentz-atreaz)**2.0)
               dist_atc_pro=sqrt((atcentx-atprodx)**2.0
     $              +(atcenty-atprody)**2.0
     $              +(atcentz-atprodz)**2.0)
               open (unit=107,file='./irc_files/hrcc.dat',
     $              status='unknown')
               write(107,*)isite,dist_atc_rea
               write(107,*)natomt1+2,dist_atc_pro
               close(107)
               call onedtau(100)
cc save HR potentials
               open (unit=108,status='unknown')
               write (108,1306)index
               rewind (108)
               read (108,1032) command1
               close (108)
 1306          format (" cp -f  ./irc_files/hrcc.me
     +           ./irc_files/hrcc.me_"I0.2 )
               call commrun(command1)
cc now read the output of HR scans
               open (unit=107, file='./irc_files/hrcc.me', status='old')
               do ir=1,nhind
                  read(107,*)
                  read(107,*)
                  read(107,*)
                  read(107,*)
                  call LineRead(107)
                  open (unit=99,status='unknown')
                  rewind (99)
                  write(99,*)word2
                  rewind(99)
                  read(99,*)numpot_ir(ir)
                  close(99)
                  read(107,*)(rotpot_ir(ik,index,ir),ik=1,numpot_ir(ir))
                  read(107,*)
               enddo
               close(107)
            enddo
         else if ((iadd.eq.1.or.ibeta.eq.1).and.
     $           (hr_rescale3.ne.'READ')) then
            istep=(numpointsf+numpointsb)/(num_hrcc_points-1)
            index=0
            do j=1,num_hrcc_points
               index=nhrcc_points(j)
               write(*,*)'test index is',index
               write(*,*)'istep is',istep
               open (unit=99, status='unknown')
               write(99,*)atgeom_me(isite,index)
               if(iadd.eq.1)then
                  write(99,*)atgeom_me(natom1+1,index)
               else if(ibeta.eq.1)then
                  write(99,*)atgeom_me(ireact,index)
               endif
               rewind(99)
               read(99,*)cjunk,atcentx,atcenty,atcentz
               read(99,*)cjunk,atreax,atreay,atreaz
               close(99)
               dist_atc_rea=sqrt((atcentx-atreax)**2.
     $              +(atcenty-atreay)**2.
     $              +(atcentz-atreaz)**2.0)
               open (unit=107,file='./irc_files/hrcc.dat',
     $              status='unknown')
c               rewind (107)
c               write(107,*)isite,dist_atc_rea
               if(iadd.eq.1)then
                  write(107,*)natomt1+1,dist_atc_rea
               else if(ibeta.eq.1)then
                  write(107,*)ireact,dist_atc_rea
               endif
               close(107)
c               goto 1010
               call onedtau(101)
cc save HR potentials
               open (unit=108,status='unknown')
               write (108,9306)index
               rewind (108)
               read (108,1032) command1
               close (108)
 9306          format (" cp -f  ./irc_files/hrcc.me
     +./irc_files/hrcc.me_"I0.2 )
               call commrun(command1)
cc now read the output of HR scans
               open(unit=107,file='./irc_files/hrcc.me',status='old')
               do ir=1,nhind
                  read(107,*)
                  read(107,*)
                  read(107,*)
                  read(107,*)
                  call LineRead(107)
                  open (unit=99,status='unknown')
                  rewind (99)
                  write(99,*)word2
                  rewind(99)
                  read(99,*)numpot_ir(ir)
                  close(99)
                  read(107,*)(rotpot_ir(ik,index,ir),ik=1,numpot_ir(ir))
                  read(107,*)
               enddo
               close(107)
            enddo
         else if (hr_rescale3.ne.'READ') then
            write(*,*)'HRCC option possible only '
            write(*,*)'for addition, abstraction and betascission'
            stop
         endif
cc now fit (or read) all the points along the MEP
         open (unit=107, file='./irc_files/hrcc_all.me', 
     $        status='unknown')
         if (hr_rescale3.ne.'READ') then
            do ir=1,nhind
               index=nhrcc_points(1)
c               index=1
               indexp1=0
               itest=1
               do j=1,numpointstot
c                  indexp1=index+istep
                  if(j.eq.nhrcc_points(itest).and.j.ne.numpointstot)then
                     jprog=0
                     ij=itest
                     index=nhrcc_points(ij)
                     indexp1=nhrcc_points(ij+1)
 1                   itest=itest+1
                  endif
                  if(j.ne.numpointstot)then
                     do ik=1,numpot_ir(ir)
                        rotpot_ir(ik,j,ir)= rotpot_ir(ik,index,ir)+
     $                       (rotpot_ir(ik,indexp1,ir)
     $                       -rotpot_ir(ik,index,ir))*
     $                       jprog/(nhrcc_points(ij+1)-nhrcc_points(ij))
                        write(107,1307)'pot irc HR ',
     $                       ik,j,ir,rotpot_ir(ik,j,ir)
                     enddo
                     jprog=jprog+1
                  else
                     do ik=1,numpot_ir(ir)
                        write(107,1307)'pot irc HR ',
     $                       ik,j,ir,rotpot_ir(ik,j,ir)
                     enddo
                     jprog=0
c                     index=index+istep
c                     index=nhrcc_points(j+1)
                  endif
               enddo
            enddo
         else
            open (unit=108, file='./irc_files/hrcc.me', status='old')
            do ir=1,nhind
               read(108,*)
               read(108,*)
               read(108,*)
               read(108,*)
               call LineRead(108)
               open (unit=99,status='unknown')
               rewind (99)
               write(99,*)word2
               rewind(99)
               read(99,*)numpot_ir(ir)
               close(99)
               read(108,*)
               read(108,*)
            enddo
            close(108)
c            write(*,*)'numpointstot is ',numpointstot
c            write(*,*)'nhind is ',nhind
c            write(*,*)'numpot_ir 1 is ',numpot_ir(1)
c            write(*,*)'numpot_ir 2 is ',numpot_ir(2)
c            stop
            do ir=1,nhind
               do j=1,numpointstot
                  do ik=1,numpot_ir(ir)
                     read(107,*)cjunk,cjunk,cjunk,
     $                    ijunk,jjunk,irjunk,rotpot_ir(ik,j,ir)
                  enddo
               enddo
            enddo
         endif
         close(107)
      endif

 1307    format (A12,1X,I5,1X,I5,1X,I5,1X,F9.2)

c         write(96,*)'numpotir1 is ',numpot_ir(1)
c         write(*,*)'numpotir1 is ',numpot_ir(1)
c         do j=1,numpot_ir(1)
c            write(96,*)'at ',j,' the pot I read is: ',rotpot_ir(j,101,1)
c         enddo
c         close(96)
c         stop
cc now we can write the ME input

cc first we determine the activation energy using HL data

      open (unit=107, file='./me_files/reac1_en.me', status='old')
      read (107,*) reac1_en
      close (107)
      open (unit=107, file='./me_files/reac1_zpe.me', status='old')
      read (107,*) reac1_zpe
      close (107)
      reac1_en=reac1_en+reac1_zpe

      if(iabs.eq.1.or.iadd.eq.1) then
         open (unit=107, file='./me_files/reac2_en.me', status='old')
         read (107,*) reac2_en
         close (107)
         open (unit=107, file='./me_files/reac2_zpe.me', status='old')
         read (107,*) reac2_zpe
         close (107)
         reac2_en=reac2_en+reac2_zpe
      else if (iiso.eq.1.or.ibeta.eq.1)then
         reac2_en=0.
      endif

      reac_en=reac1_en+reac2_en

      open (unit=107, file='./me_files/ts_en.me', status='old')
      read (107,*) ts_en
      close (107)
      open (unit=107, file='./me_files/ts_zpe.me', status='old')
      read (107,*) ts_zpe
      close (107)
      ts_en=ts_en+ts_zpe
      ts_en=ts_en-reac_en
      ts_en = ts_en*cautokcal

cc then we determine the energies along the reaction path

cc this is to be consistent with ME that does not include the 
cc lowest eigenvalue of the HR in calculating the HR partition function
      
      nfreq=3*natom-1-6

      do inumpoints = 1, numpointstot
         zpe=0.
         do j=1,nfreq
            zpe=zpe+freqRTproj(j,inumpoints)
         enddo
         zpe=zpe/(cautoicm*2.d0)
         zpe_irc(inumpoints)=zpe
         rc_ene_kcal(inumpoints)=(rc_ene(inumpoints)+zpe)*cautokcal
         if(inumpoints.eq.numpointsf+1)Ets=rc_ene_kcal(inumpoints)
      enddo
c now rescale with respect to TS energy
      do inumpoints = 1, numpointstot
         rc_ene_kcal(inumpoints)=rc_ene_kcal(inumpoints)-Ets
      enddo

c now determine energies with respect to reactants
c      write(*,*)'ts en',ts_en
      do inumpoints = 1, numpointstot
         rc_ene_kcal(inumpoints)=rc_ene_kcal(inumpoints)+ts_en
c         write(*,*)inumpoints,rc_ene_kcal(inumpoints)
      enddo

cc now rescale potential if requested

      if(iresirc.eq.1) then

cc first perform High level calculations, unless they have already been done
cc in which the skipall option is used to bypass this step

         if(iskiphl.ne.1) then
            numproc=numprochl
            open (unit=107, file='input.xyz', status='unknown')
            write (107,*)
            do iatom = 1, natom
               write (107,'(A80)') atgeom_me(iatom,1)
            enddo 
            close(107)
            call hl(100)
            if(numpointstot.ne.numpointsf+1) then
               open (unit=107, file='input.xyz', status='unknown')
               write (107,*)
               do iatom = 1, natom
                  write (107,'(A80)') atgeom_me(iatom,numpointstot)
               enddo
               close(107)
               call hl(101)
            endif
            numproc=numprocll
         endif
cc now read energies of starting point,TS, and arrival point
         open (unit=107, file='./me_files/ircst_en.me', status='old')
         read (107,*) en_ircst
         close(107)
         if(numpointstot.ne.numpointsf+1) then
            open (unit=107, file='./me_files/ircarr_en.me', 
     $            status='old')
            read (107,*) en_ircarr
            close(107)
         endif
         open (unit=107, file='./me_files/ts_en.me', status='old')
         read (107,*) en_ts
         close(107)
cc determine scaling factors for initial and final points
         destart_hl=en_ts-en_ircst
         destart_l1=rc_ene(numpointsf+1)-rc_ene(1)
         scalestart=destart_hl/destart_l1
         if(numpointstot.ne.numpointsf+1) then
            dearr_hl=en_ts-en_ircarr
            dearr_l1=rc_ene(numpointsf+1)-rc_ene(numpointstot)
            scalearr=dearr_hl/dearr_l1
         else
            scalearr=1
         endif

         write(*,*)'scalestart is',scalestart
         write(*,*)'scalearr is',scalearr
cc       now determine rescaled energies from start
         
         do inumpoints = 1, numpointstot
            rc_ene_kcal(inumpoints)=0
         enddo

         do inumpoints = 1, numpointsf+1
            rc_ene_kcal(inumpoints)=(rc_ene(inumpoints)
     $           -rc_ene(numpointsf+1))*scalestart
     $          +zpe_irc(inumpoints)-zpe_irc(numpointsf+1)
            rc_ene_kcal(inumpoints)=rc_ene_kcal(inumpoints)
     $           *cautokcal
         enddo
         if(numpointstot.ne.numpointsf+1) then
            do inumpoints = numpointsf+1,numpointstot
               rc_ene_kcal(inumpoints)=(rc_ene(inumpoints)
     $              -rc_ene(numpointsf+1))*scalearr
     $              +zpe_irc(inumpoints)-zpe_irc(numpointsf+1)
               rc_ene_kcal(inumpoints)=rc_ene_kcal(inumpoints)
     $              *cautokcal
            enddo
         endif
c     now determine energies with respect to reactants
c      write(*,*)'ts en',ts_en
         do inumpoints = 1, numpointstot
            rc_ene_kcal(inumpoints)=rc_ene_kcal(inumpoints)+ts_en
c            write(*,*)inumpoints,rc_ene_kcal(inumpoints)
         enddo

cc now rescale rc_ene vector and rewrite potential file for SCT
         do inumpoints = 1, numpointsf+1
            rc_ene(inumpoints)=rc_ene(inumpoints)*scalestart
         enddo
         if(numpointstot.ne.numpointsf+1) then
            do inumpoints = numpointsf+1,numpointstot
               rc_ene(inumpoints)=rc_ene(inumpoints)*scalearr
            enddo
         endif

         open (unit=333,file='RPHt_coord_en.dat',status='unknown')

         write(333,*)'Point Coordinate Energy Bond1 Bond2'

         do i=1,numpointsf+numpointsb+1
            open (unit=99, status='unknown')
               write(99,*)atgeom_me(isite,i)
               write(99,*)atgeom_me(jsite,i)
            if(iabs.eq.1.or.iadd.eq.1)then
               write(99,*)atgeom_me(natom1+1,i)
            else if(iiso.eq.1.or.ibeta.eq.1)then
               write(99,*)atgeom_me(ireact,i)
            endif
            rewind(99)
            read(99,*)cjunk,atcentx,atcenty,atcentz
            read(99,*)cjunk,atreax,atreay,atreaz
            read(99,*)cjunk,atprodx,atprody,atprodz
            close(99)
            dist_atc_rea=sqrt((atcentx-atreax)**2.
     $       +(atcenty-atreay)**2.
     $       +(atcentz-atreaz)**2.0)
            dist_atc_pro=sqrt((atcentx-atprodx)**2.0
     $       +(atcenty-atprody)**2.0
     $       +(atcentz-atprodz)**2.0)
 
c            write(*,*)'energy = ',rc_ene(i)
            write(333,1333)i,rc_coord(i),rc_ene(i),dist_atc_rea,
     $           dist_atc_pro
         enddo
c      do i=1,numpointsf+numpointsb+1
c         write(*,*)'FC = ',force_con(2,i)
c      enddo
         close(333)
c         stop   
      endif

cc now rescale potential if requested

      if(iresirc.eq.2) then

cc perform High level calculations for all points, unless they have already been done
cc in which case the skipall option is used to bypass this step

         if(iskiphl.ne.1) then
            numproc=numprochl
            do inump=1,numpointstot
               open (unit=107, file='input.xyz', status='unknown')
               write (107,*)
               do iatom = 1, natom
                  write (107,'(A80)') atgeom_me(iatom,inump)
               enddo 
               close(107)
               call hl(100)
               open (unit=107,file='./me_files/ircst_en.me',
     $               status='old')
               read (107,*) en_ircst
               write (96,*)'HL en of point ',inump,' is ', en_ircst
               close(107)
               rc_ene_hl(inump)=en_ircst
            enddo
            numproc=numprocll
            open (unit=107, file='./irc_files/hl_en.dat', 
     $            status='unknown')
            do inump=1,numpointstot
               write (107,*)inump,rc_ene_hl(inump)
            enddo
            close(107)
         else
            open (unit=107, file='./irc_files/hl_en.dat', 
     $            status='unknown')
            do inump=1,numpointstot
               read (107,*)numb,rc_ene_hl(inump)
            enddo
            close(107)
         endif

cc       now determine  HL energy profile         
         do inumpoints = 1, numpointstot
            rc_ene_kcal(inumpoints)=0
         enddo

         do inumpoints = 1, numpointstot
            rc_ene_kcal(inumpoints)=(rc_ene_hl(inumpoints)
     $           -rc_ene_hl(numpointsf+1))
     $          +zpe_irc(inumpoints)-zpe_irc(numpointsf+1)
            rc_ene_kcal(inumpoints)=rc_ene_kcal(inumpoints)
     $           *cautokcal
         enddo
c     now determine energies with respect to reactants
         do inumpoints = 1, numpointstot
            rc_ene_kcal(inumpoints)=rc_ene_kcal(inumpoints)+ts_en
         enddo

cc now rewrite rc_ene vector 

         open (unit=107, file='./me_files/ts_en.me', status='old')
         read (107,*) en_ts
         close(107)
c
         do inumpoints = 1, numpointstot
            rc_ene(inumpoints)=rc_ene_hl(inumpoints)-en_ts
         enddo

         open (unit=333,file='RPHt_coord_en.dat',status='unknown')

         write(333,*)'Point Coordinate Energy Bond1 Bond2'

         do i=1,numpointsf+numpointsb+1
            open (unit=99, status='unknown')
            write(99,*)atgeom_me(isite,i)
            write(99,*)atgeom_me(jsite,i)
            write(99,*)atgeom_me(natom1+1,i)
            rewind(99)
            read(99,*)cjunk,atcentx,atcenty,atcentz
            read(99,*)cjunk,atreax,atreay,atreaz
            read(99,*)cjunk,atprodx,atprody,atprodz
            close(99)
            dist_atc_rea=sqrt((atcentx-atreax)**2.
     $       +(atcenty-atreay)**2.
     $       +(atcentz-atreaz)**2.0)
            dist_atc_pro=sqrt((atcentx-atprodx)**2.0
     $       +(atcenty-atprody)**2.0
     $       +(atcentz-atprodz)**2.0)
 
c            write(*,*)'energy = ',rc_ene(i)
            write(333,1333)i,rc_coord(i),rc_ene(i),dist_atc_rea,
     $           dist_atc_pro
         enddo
c      do i=1,numpointsf+numpointsb+1
c         write(*,*)'FC = ',force_con(2,i)
c      enddo
         close(333)
c         stop   
      endif
      
cc now perform multidimensional tunneling calculations
cc and prepare input file for me tunneling

      open (unit=99,status='unknown')
      rewind (99)
      write (99,1021)
      rewind (99)
      read (99,1020) command1
      close (99)
      call commrun(command1)

 1020 format (A100)
 1021 format ("sed -i '/Cartesian/d' RPHt_all_data.dat")       

      command1=' cp -f RPHt_all_data.dat RPHt_input_data.dat'
      call commrun(command1)

      command1='RPHt.exe > ./md_tunn/mdtunn.out'
      call commrun(command1)

cc save relavent files for md tunn
      command1='cp -f trajec.xyz ./md_tunn'
      call commrun(command1)
      command1='cp -f VaG.txt ./md_tunn'
      call commrun(command1)
      command1='cp -f mueff.txt ./md_tunn'
      call commrun(command1)

cc finally we have all we need to write the variational input!
       
c      open (unit=112,file='./irc_files/test.out',status='unknown')
c      do iatom=1,natom
c         write (112,*) atgeom_me(iatom,1)
c      enddo
c      close(112)
c      stop

      open(unit=107,file='./me_files/de1_TSvar.me',
     +     status='unknown')
      write(107,*)rc_ene_kcal(1)-ts_en
      close(107)

      open(unit=107,file='./me_files/variational.me',
     +     status='unknown')
      if(intfreq.eq.1)then
         open(unit=109,file='./me_files/variational_xyz.me',
     +        status='unknown')
      endif


      if (iabs.eq.1) then
         if (nts.eq.1) write (107,*) 'Barrier TS REACS WR'
         if (nts.eq.2) write (107,*) 'Barrier B2 WR PRODS'
         if (nts.eq.3) write (107,*) 'Barrier B2 WR WP'
         if(intfreq.eq.1)then
            if (nts.eq.1) write (109,*) 'Barrier TS REACS WR'
            if (nts.eq.2) write (109,*) 'Barrier B2 WR PRODS'
            if (nts.eq.3) write (109,*) 'Barrier B2 WR WP'
         endif
      endif
      if (iadd.eq.1) then
         if (nts.eq.1) write (107,*) 'Barrier TS REACS WP'
         if (nts.eq.2) write (107,*) 'Barrier B2 WR WP'
         if(intfreq.eq.1)then
         if (nts.eq.1) write (109,*) 'Barrier TS REACS WP'
         if (nts.eq.2) write (109,*) 'Barrier B2 WR WP'
         endif
      endif
      if (iiso.eq.1) then
         if(ipr1.eq.0)then
            write (107,*) 'Barrier TS REACS WP'
            if(intfreq.eq.1)then
               write (109,*) 'Barrier TS REACS WP'
            endif
         else
            write (107,*) 'Barrier TS REACS PRODS'
            if(intfreq.eq.1)then
               write (109,*) 'Barrier TS REACS PRODS'
            endif
         endif
      endif
      if (ibeta.eq.1) then
         write (107,*) 'Barrier TS REACS PRODS'
            if(intfreq.eq.1)then
               write (109,*) 'Barrier TS REACS PRODS'
            endif
      endif

      write(107,*)'Variational'
      if(intfreq.eq.1)then
         write(109,*)'Variational'
      endif
      do inumpoints=1,numpointstot
         nfreqtest=3*natom-6-1-nhind
         nfreqw=0
         nfreqwxyz=0
         intfreqw=0
         if(intfreq.eq.1)then
            intfreqw=1
         endif
         if(intfreq.eq.0)then
            do j=1,nfreq
               if(freqproj(j,inumpoints).ne.0)then
                  nfreqw=nfreqw+1
               endif
            enddo
         else
            do j=1,nfreq
               if(freqintproj(j,inumpoints).ne.0)then
                  nfreqw=nfreqw+1
               endif
            enddo
            do j=1,nfreq
               if(freqproj(j,inumpoints).ne.0)then
                  nfreqwxyz=nfreqwxyz+1
               endif
            enddo
         endif
         if(nfreqw.ne.nfreqtest) goto 9999
         if(nfreqwxyz.ne.nfreqtest) intfreqw=0

         write(107,*)'RRHO            !  ',inumpoints
         write(107,*)'Geometry[angstrom] ',natom
         do j=1,natom
            write(107,'(A80)')atgeom_me(j,inumpoints)
         enddo
         if(nhindmd.eq.0)then
            write(107,*)'	Core 	RigidRotor'
            write(107,*)'          SymmetryFactor ',symf
            write(107,*)'End'
         endif
         if(intfreqw.eq.1)then
            write(109,*)'RRHO            !  ',inumpoints
            write(109,*)'Geometry[angstrom] ',natom
            do j=1,natom
               write(109,'(A80)')atgeom_me(j,inumpoints)
            enddo
            if(nhindmd.eq.0)then
               write(109,*)'	Core 	RigidRotor'
               write(109,*)'          SymmetryFactor ',symf
               write(109,*)'End'
            endif
         endif

         if(nhind.ne.0.or.nhindmd.ne.0)then
            if(nhindmd.eq.0)then
               open(unit=108,file='./me_files/ts_hr.me',
     +              status='unknown')
            endif
            if(nhindmd.ne.0)then
               open(unit=108,file='./me_files/ts_mdhr_nofr.me',
     +              status='unknown')
            endif
            do ir=1,nhindmd
               do ij=1,24
                  read(108,'(A70)')rot2dline
                  write(107,*)rot2dline
                  if(intfreqw.eq.1) write(109,*)rot2dline
               enddo
            enddo
            nhind_res=nhind-2*nhindmd
            do ir=1,nhind_res
               read(108,'(A60)')rotline
               write(107,*)rotline
               if(intfreqw.eq.1) write(109,*)rotline
               read(108,'(A60)')rotline
               write(107,*)rotline
               if(intfreqw.eq.1) write(109,*)rotline
               read(108,'(A60)')rotline
               write(107,*)rotline
               if(intfreqw.eq.1) write(109,*)rotline
               read(108,'(A60)')rotline
               write(107,*)rotline
               if(intfreqw.eq.1) write(109,*)rotline
               call LineRead(108)
               open (unit=99,status='unknown')
               rewind (99)
               write(99,*)word2
               rewind(99)
               read(99,*)numpot
               close(99)
             write(107,*)' Potential[kcal/mol] ',numpot 
             if(intfreqw.eq.1)write(109,*)' Potential[kcal/mol] ',numpot 
               read(108,*)(rotpot(j),j=1,numpot)
               if(hr_rescale.eq.'HRCC') then
                  do j=1,numpot
                    rotpot(j)= rotpot_ir(j,inumpoints,ir)
                 enddo
               endif
               write(107,1111)(rotpot(j),j=1,numpot)
               if(intfreqw.eq.1)write(109,1111)(rotpot(j),j=1,numpot)
               read(108,'(A60)')rotline
               write(107,*)rotline
               if(intfreqw.eq.1) write(109,*)rotmore hrline
            enddo
            close(108)
         endif
         nfreq=3*natom-1-6

cc before writing frequencies we determine how many are not zero
cc since along the reaction path the negative freqs that may emerge
cc are set to zero
      
         nfreqw=0
         if(intfreq.eq.0)then
            do j=1,nfreq
               if(freqproj(j,inumpoints).ne.0)then
                  nfreqw=nfreqw+1
               endif
            enddo
         else
            do j=1,nfreq
               if(freqintproj(j,inumpoints).ne.0)then
                  nfreqw=nfreqw+1
               endif
            enddo
         endif
         write(107,*)'    Frequencies[1/cm] ',nfreqw
         if(intfreqw.eq.1) write(109,*)'    Frequencies[1/cm] ',nfreqw
         if(intfreq.eq.0)then
            write(107,8010) (freqproj(j,inumpoints),j=1,nfreqw)            
         else
            write(107,8010) (freqintproj(j,inumpoints),j=1,nfreqw)
            if(intfreqw.eq.1)then
               write(109,8010) (freqproj(j,inumpoints),j=1,nfreqw)
            endif
         endif
         write(107,*)'ZeroEnergy[kcal/mol] ',rc_ene_kcal(inumpoints)
         write (107,*) ' ElectronicLevels[1/cm]           ',nelec
         do ielec = 1, nelec
            write (107,*) eelec(ielec),gelec(ielec)
         enddo
         write(107,*)'End'
         write(107,*)'!*********************************************'
         if(intfreqw.eq.1)then
            write(109,*)'ZeroEnergy[kcal/mol] ',rc_ene_kcal(inumpoints)
            write (109,*) ' ElectronicLevels[1/cm]           ',nelec
            do ielec = 1, nelec
               write (109,*) eelec(ielec),gelec(ielec)
            enddo
            write(109,*)'End'
            write(109,*)'!*********************************************'
         endif
 9999    continue
      enddo

      open(unit=108,file='./me_files/ts_fr.me',
     +           status='unknown')
      
      if(ifrozrts.ne.1)then
         do while (WORD.NE.'IMAGINARYFREQUENCY[1/CM]')
            call LineRead (108)
            if (WORD.EQ.'END') then
               write (96,*) 'Imag freq must be defined in 
     +./me_files/ts_fr.me'
               stop
            endif
         enddo
         read(word2,'(f9.3)') freq_imag
         close(108)
         if(imdtunn.ne.1) then
            write(107,*)'      Tunneling    Eckart'
            write(107,*)'      ImaginaryFrequency[1/cm] ',freq_imag
            write(107,*)'      WellDepth[kcal/mol]       $wdepfor'
            write(107,*)'      WellDepth[kcal/mol]       $wdepback'
         else
            write(107,*)'      Tunneling    Read'
            write(107,*)'   CutoffEnergy[1/cm]       2500'
            write(107,*)'   ImaginaryFrequency[1/cm]   ',freq_imag
            write(107,*)'   File imactint.dat'
         endif
         if(intfreqw.eq.1) then
            if(imdtunn.ne.1) then
               write(109,*)'      Tunneling    Eckart'
               write(109,*)'      ImaginaryFrequency[1/cm] ',freq_imag
               write(109,*)'      WellDepth[kcal/mol]       $wdepfor'
               write(109,*)'      WellDepth[kcal/mol]       $wdepback'
            else
               write(109,*)'      Tunneling    Read'
               write(109,*)'   CutoffEnergy[1/cm]       2500'
               write(109,*)'   ImaginaryFrequency[1/cm]   ',freq_imag
               write(109,*)'   File imactint.dat'
            endif
         endif
      endif
      write(107,*)'End '
      if(intfreqw.eq.1) write(109,*)'End '

c      write(107,*)'End '
c      write(107,*)'End '

      close (107)
      close (109)
      close (unit=96,status='keep')

      if(intfreq.eq.1)then
         command1='cp -f ./me_files/variational.me 
     + ./me_files/variational_int.me'
         call commrun(command1)
      endif

cc write frequencies along MEP
ccccc
      open(unit=17,file='freq_proj_xyz.out',status='unknown')
      do inumpoints=1,numpointstot
         write(17,8001)rc_coord(inumpoints),
     +                 (freqproj(j,inumpoints),j=1,nfreqw)
      enddo
      close(17)

      if(intfreq.eq.1)then
         open(unit=16,file='freq_proj_int.out',status='unknown')
         do inumpoints=1,numpointstot
            write(16,8001)rc_coord(inumpoints),
     +                    (freqintproj(j,inumpoints),j=1,nfreqw)
         enddo
         close(16)
      endif

      command1='cp -f ./freq_proj_xyz.out ./md_tunn'
      call commrun(command1)

      if(intfreq.eq.1)then
         command1='cp -f ./freq_proj_int.out ./md_tunn'
         call commrun(command1)
      endif


 2000 format (A80,1X,I10)
 2001 format (A160)
 1111 format (100(1X,f7.2))
 8001 format (1x,f7.2,1x,100G12.5)
 8010 format (1x,10G12.5)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine hl(ispecies)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      character*70 comline1,comline2
      character*30 intcoor(3*natommx)
      character*60 atomlabel(natommx)
      character*100 command1
      character*100 commandcopy
      character*20 bislab(ntaumx)
      character*30 gmem
      character*160 word_hl

      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim) 
      dimension taumn(ntaumx),taumx(ntaumx)

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      character *30 code_name

      include 'filcomm.f'

      
      call lineread(0)
      
c open files

      open (unit=66,file='./output/hl_en.out',status='unknown')
      open (unit=21,file='./data/theory.dat',status='unknown')

      do while (WORD.NE.'HLEVEL')
         call LineRead (21)
         if (WORD.EQ.'END') then
            write (66,*) 'code for HL must be defined'
            stop
         endif
      enddo
c      call LineRead (21)
      code_name=word2
      if(code_name.eq.'MOLPRO')then
         ilev1code=2
         commandcopy='cp -f ./data/hl_molpro.dat ./hl_molpro.dat'       
      else if (code_name.eq.'G09')then
         ilev1code=1
      else
         write(66,*) 'code not supported in high level calculations'
         write(66,*) code_name
         stop
      endif


      if (word3.eq.'RESCALE')then
         if(iadd.eq.1)then
            if(ispecies.eq.0)then
               open (unit=99,file='./me_files/reac1_en.me',status='old')
               read(99,*)hl_reac1_en
               close(99)
               open (unit=99,file='./me_files/reac2_en.me',status='old')
               read(99,*)hl_reac2_en
               close(99)
               open (unit=99,file='./me_files/wellp_en.me',status='old')
               read(99,*)hl_wp_en
               close(99)
               open (unit=99,file='./geoms/tsgta_l1.xyz',status='old')
               read(99,*)
               read(99,*)ts_en_l1
               close(99)
               open (unit=99,file='./geoms/reac1_l1.xyz',status='old')
               read(99,*)
               read(99,*)reac1_en_l1
               close(99)
               open (unit=99,file='./geoms/reac2_l1.xyz',status='old')
               read(99,*)
               read(99,*)reac2_en_l1
               close(99)
               open (unit=99,file='./geoms/wellp_l1.xyz',status='old')
               read(99,*)
               read(99,*)wp_en_l1
               close(99)

               scale_fact=(hl_wp_en-hl_reac1_en-hl_reac2_en)/
     $                 (wp_en_l1-reac1_en_l1-reac2_en_l1)

               hl_scaled_en=(ts_en_l1-reac1_en_l1-reac2_en_l1)
     $                      *scale_fact+hl_reac1_en+hl_reac2_en
               write(66,*)'scale fact is ',scale_fact
               write(66,*)' energy changes without ZPE corrections'
               write(66,*)'DE HL is         ',(hl_wp_en-hl_reac1_en
     $              -hl_reac2_en)*627.5
               write(66,*)'DE Level1 is     ',(wp_en_l1-reac1_en_l1
     $              -reac2_en_l1)*627.5
               write(66,*)'Unscaled Eact is ',(ts_en_l1-reac1_en_l1
     $              -reac2_en_l1)*627.5
               write(66,*)'Scaled Eact is   ',(hl_scaled_en-hl_reac1_en
     $               -hl_reac2_en)*627.5

               open (unit=99,file='./me_files/ts_en.me',
     $               status='unknown')
               write(99,*)hl_scaled_en
               close(99)
               close(21)
               goto 999
            endif
         endif
         if(ispecies.eq.0)then
            open (unit=99,file='./me_files/reac1_en.me',status='old')
            read(99,*)hl_react_en
            close(99)
            open (unit=99,file='./me_files/prod1_en.me',status='old')
            read(99,*)hl_prod1_en
            close(99)
            open (unit=99,file='./me_files/prod2_en.me',status='old')
            read(99,*)hl_prod2_en
            close(99)
            open (unit=99,file='./geoms/tsgta_l1.xyz',status='old')
            read(99,*)
            read(99,*)ts_en_l1
            close(99)
            open (unit=99,file='./geoms/reac1_l1.xyz',status='old')
            read(99,*)
            read(99,*)reac1_en_l1
            close(99)
            open (unit=99,file='../100/geoms/tsgta_l1.xyz',status='old')
            read(99,*)
            read(99,*)prod_en_l1
            close(99)
            scale_fact=(hl_prod1_en+hl_prod2_en-hl_react_en)/
     $                 (prod_en_l1-reac1_en_l1)
            hl_scaled_en=(ts_en_l1-reac1_en_l1)*scale_fact+hl_react_en
            write(66,*)'scale fact is ',scale_fact
            write(66,*)' energy changes without ZPE corrections'
            write(66,*)'DE HL is         ',(hl_prod1_en+hl_prod2_en
     $      -hl_react_en)*627.5
            write(66,*)'DE Level1 is     ',(prod_en_l1-reac1_en_l1)
     $                             *627.5
            write(66,*)'Unscaled Eact is ',(ts_en_l1-reac1_en_l1)
     $                             *627.5
            write(66,*)'Scaled Eact is   ',(hl_scaled_en-hl_react_en)
     $                             *627.5

            open (unit=99,file='./me_files/ts_en.me',status='unknown')
            write(99,*)hl_scaled_en
            close(99)
            close(21)
            goto 999
         endif
      endif

c         read(99,*)
c         read(99,*)l1_prod1_en
c         close(99)




cc now check for code specific/species specific input commands

      rewind(21)
      WORD_HL='HLEVEL'
      if(ispecies.eq.1)then
         do while (WORD.NE.'HLEVEL_1')
            call LineRead(21)
            if (WORD.EQ.'END') then
               go to  900
            endif
         enddo
         WORD_HL='HLEVEL_1'
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/hl_reac1_molpro.dat 
     $ hl_molpro.dat'            
         endif
      endif
      if(ispecies.eq.2)then
         do while (WORD.NE.'HLEVEL_2')
            call LineRead(21)
            if (WORD.EQ.'END') then
               go to  900
            endif
         enddo
         WORD_HL='HLEVEL_2'
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/hl_reac2_molpro.dat 
     $ hl_molpro.dat'            
         endif
      endif
      if(ispecies.eq.3)then
         do while (WORD.NE.'HLEVEL_3')
            call LineRead(21)
            if (WORD.EQ.'END') then
               go to  900
            endif
         enddo
         WORD_HL='HLEVEL_3'
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/hl_prod1_molpro.dat 
     $ hl_molpro.dat'            
         endif
      endif
      if(ispecies.eq.4)then
         do while (WORD.NE.'HLEVEL_4')
            call LineRead(21)
            if (WORD.EQ.'END') then
               go to  900
            endif
         enddo
         WORD_HL='HLEVEL_4'
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/hl_prod2_molpro.dat 
     $ hl_molpro.dat'            
         endif
      endif
      if(ispecies.eq.5)then
         do while (WORD.NE.'HLEVEL_5')
            call LineRead(21)
            if (WORD.EQ.'END') then
               go to  900
            endif
         enddo
         WORD_HL='HLEVEL_5'
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/hl_wellr_molpro.dat 
     $ hl_molpro.dat'            
         endif
      endif
      if(ispecies.eq.51)then
         do while (WORD.NE.'HLEVEL_51')
            call LineRead(21)
            if (WORD.EQ.'END') then
               go to  900
            endif
         enddo
         WORD_HL='HLEVEL_51'
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/hl_wellr_molpro.dat 
     $ hl_molpro.dat'            
         endif
      endif
      if(ispecies.eq.6)then
         do while (WORD.NE.'HLEVEL_6')
            call LineRead(21)
            if (WORD.EQ.'END') then
               go to  900
            endif
         enddo
         WORD_HL='HLEVEL_6'
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/hl_wellp_molpro.dat 
     $ hl_molpro.dat'            
         endif
      endif
      if(ispecies.eq.61)then
         do while (WORD.NE.'HLEVEL_61')
            call LineRead(21)
            if (WORD.EQ.'END') then
               go to  900
            endif
         enddo
         WORD_HL='HLEVEL_61'
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/hl_wellp_molpro.dat 
     $ hl_molpro.dat'            
         endif
      endif
      if(ispecies.eq.11)then
         do while (WORD.NE.'HLEVEL_11')
            call LineRead(21)
            if (WORD.EQ.'END') then
               go to  900
            endif
         enddo
         WORD_HL='HLEVEL_11'
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/hl_reacs_molpro.dat 
     $ hl_molpro.dat'            
         endif
      endif
      if(ispecies.eq.12)then
         do while (WORD.NE.'HLEVEL_12')
            call LineRead(21)
            if (WORD.EQ.'END') then
               go to  900
            endif
         enddo
         WORD_HL='HLEVEL_12'
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/hl_prods_molpro.dat 
     $ hl_molpro.dat'            
         endif
      endif
      if(ispecies.eq.0)then
         do while (WORD.NE.'HLEVEL_TS')
            call LineRead(21)
            if (WORD.EQ.'END') then
               go to  900
            endif
         enddo
         WORD_HL='HLEVEL_TS'
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/hl_ts_molpro.dat 
     $ hl_molpro.dat'            
         endif
      endif
      if(ispecies.eq.100.or.ispecies.eq.101)then
         do while (WORD.NE.'HLEVEL_TS')
            call LineRead(21)
            if (WORD.EQ.'END') then
               go to  900
            endif
         enddo
         WORD_HL='HLEVEL_TS'
         if(WORD2.eq.'MOLPRO')then
            commandcopy='cp -f ./data/hl_ts_molpro.dat 
     $ hl_molpro.dat'            
         endif
      endif

c      call commrun(commandcopy)

 900  continue
      close(21)

      if (ispecies.eq.1) then
         open (unit=15,file='./data/reac1.dat',status='old')
         open (unit=17,file='./output/reac1_opt.out',status='unknown')
         inp_type=1
      endif
      if (ispecies.eq.2) then
         open (unit=15,file='./data/reac2.dat',status='old')
         open (unit=17,file='./output/reac2_opt.out',status='unknown')
         inp_type=1
      endif
      if (ispecies.eq.3) then
         open (unit=15,file='./data/prod1.dat',status='old')
         open (unit=17,file='./output/prod1_opt.out',status='unknown')
         inp_type=1
      endif
      if (ispecies.eq.4) then
         open (unit=15,file='./data/prod2.dat',status='old')
         open (unit=17,file='./output/prod2_opt.out',status='unknown')
         inp_type=1
      endif
      if (ispecies.eq.5) then
         open (unit=15,file='./data/wellr.dat',status='old')
         open (unit=17,file='./output/wellr_opt.out',status='unknown')
         inp_type=1
      endif
      if (ispecies.eq.51) then
         open (unit=15,file='./data/wellr.dat',status='old')
         open (unit=17,file='./output/wellr_opt.out',status='unknown')
         inp_type=2
      endif
      if (ispecies.eq.6) then
         open (unit=15,file='./data/wellp.dat',status='old')
         open (unit=17,file='./output/wellp_opt.out',status='unknown')
         inp_type=1
      endif
      if (ispecies.eq.61) then
         open (unit=15,file='./data/wellp.dat',status='old')
         open (unit=17,file='./output/wellp_opt.out',status='unknown')
         inp_type=2
      endif
      if (ispecies.eq.11) then
         open (unit=15,file='./data/reacs.dat',status='old')
         open (unit=17,file='./output/reacs_opt.out',status='unknown')
         inp_type=1
      endif
      if (ispecies.eq.12) then
         open (unit=15,file='./data/prods.dat',status='old')
         open (unit=17,file='./output/prods_opt.out',status='unknown')
         inp_type=1
      endif
      if (ispecies.eq.0) then
         open (unit=15,file='./data/ts.dat',status='old')
         open (unit=17,file='./output/ts_opt.out',status='unknown')
         inp_type=2
      endif
      if (ispecies.eq.100.or.ispecies.eq.101) then
         open (unit=15,file='./data/ts.dat',status='old')
         open (unit=17,file='input.xyz',status='unknown')
         inp_type=2
      endif

      if (inp_type.eq.1) then

         do while (WORD.NE.'NATOM')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'natom must be defined'
               stop
            endif
         enddo
         read (15,*) natom,natomt,ilin
         rewind(15)

         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'sampling coordinates must be defined'
               stop
            endif
         enddo
         read (15,*) ntau
         rewind(15)

c read and write z mat input 

         if (idebug.ge.2) write (6,*) ' starting zmat input'

         do while (WORD.NE.'CHARGE')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'charge and spin must be defined'
               stop
            endif
         enddo
         read (15,*) icharge,ispin

         if(code_name.eq.'MOLPRO')then

            call commrun(commandcopy)
c            write (67,*)'geometry={angstrom'
            do iatom = 1 , natomt
               read (15,'(A60)') atomlabel(iatom)
c               write (67,*)atomlabel(iatom)
            enddo
c            write (67,*)'}'
            rewind(15)

         else if (code_name.eq.'G09')then   
            open (unit=21,file='./data/theory.dat',status='unknown')
            do while (WORD.NE.WORD_HL)
               call LineRead (21)
               if (WORD.EQ.'END') then
                  write (66,*) ' HL g09 theory must be defined'
                  write (66,*) ' in file theory.dat '
                  stop
               endif
            enddo
            read (21,'(A70)') comline1
            read (21,'(A70)') comline2
            close(21)

c            write(67,*) comline1
c            write(67,*) comline2
c            write(67,*)
c            write(67,*)' gaussian geom'
c            write(67,*)
c            write (67,*) icharge,ispin
            do iatom = 1 , natomt
               read (15,'(A60)') atomlabel(iatom)
c               write (67,*)atomlabel(iatom)
            enddo
c            write(67,*)
            rewind(15)
         else
            write(66,*) 'code not supported'
            write(66,*) code_name
            stop
         endif
           

cc read coordinate names

         if (natom.ne.1) then
            do while (WORD.NE.'INTCOOR')
               call LineRead (15)
               if (WORD.EQ.'END') then
                  write (66,*) 'internal coordinates must be defined'
                  stop
               endif
            enddo
            ncoord = 3*natom-6-ntau
            if (natom.eq.1) ncoord = 0
            if (natom.eq.2) ncoord = 1
            do icoord = 1 , ncoord
               call LineRead (15)
               intcoor(icoord) = word
            enddo
            rewind(15)
         endif

         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'sampling coordinates must be defined'
               stop
            endif
         enddo
         read (15,*) ntau
         ntau_fr=ntau
         if (ntau.gt.ntaumx) then
            write (66,*) 'ntau too large',ntau,ntaumx
            stop
         endif
         if (ntau.ne.0) then
            read (15,*)
            do itau = 1 , ntau
               read (15,*) bislab(itau),taumn(itau),taumx(itau)
               write (66,*) bislab(itau),taumn(itau),taumx(itau)
               open (unit=99,status='unknown')
               rewind (99)
               write (99,*) bislab(itau)
               rewind (99)
               call LineRead (99)
               icoord=ncoord+itau
               intcoor(icoord)=WORD
               close (unit=99,status='keep')
            enddo
         endif
         rewind(15)

         close (unit=15,status='keep')

cc now read optimized geometry parameters
         ncoord = 3*natom-6
         if (natom.eq.1) ncoord=0
         if (natom.eq.2) ncoord=1
         read (17,*)
         do icoord = 1 , ncoord
            call LineRead (17)
c         xinti(icoord) = word
c         intcoori(icoord) = word
            OPEN (unit=99,status='unknown')
            REWIND (99)
            WRITE (99,2000) WORD
 2000       FORMAT (A30)
            REWIND (99)
            READ (99,*) xint(icoord)
            close (unit=99,status='keep')
         enddo
         close (unit=17,status='keep')

cc now write z-matrix parameters
         
c         if(code_name.eq.'MOLPRO')then
c            do icoord = 1 , ncoord
c               write (67,*) intcoor(icoord),' = ',xint(icoord)
c            enddo
c            write (67,*)
c         else if (code_name.eq.'G09')then   
c            do icoord = 1 , ncoord
c               write (67,*) intcoor(icoord),xint(icoord)
c            enddo
c            write (67,*)
c         endif

cc now read input of type 2

      else if (inp_type.eq.2) then

cc we assume the TS is not linear         
         ilin=0

         open (unit=25,file='./data/reac1.dat',status='old')
         word=' '
         do while (WORD.NE.'NTAU')
            call LineRead (25)
            if (WORD.EQ.'END') then
               write (66,*) 'sampling coords of reac1 must be defined'
               stop
            endif
         enddo
         read (25,*) ntau1
         rewind 25

c     natomt is to account for dummy atoms
         do while (WORD.NE.'NATOM')
            call LineRead (25)
            if (WORD.EQ.'END') then
               write (66,*) 'natom in reac1 must be defined'
               stop
            endif
         enddo
         read (25,*) natom1,natomt1
         close (25)

cc get data from react2 file
         if(iadd.eq.1.or.iabs.eq.1) then
            open (unit=25,file='./data/reac2.dat',status='old')

            do while (WORD.NE.'NTAU')
               call LineRead (25)
               if (WORD.EQ.'END') then
                  write (66,*) 'sampl coords of reac2 must be defined'
                  stop
               endif
            enddo
            read (25,*) ntau2
            rewind 25

c     natomt is to account for dummy atoms
            do while (WORD.NE.'NATOM')
               call LineRead (25)
               if (WORD.EQ.'END') then
                  write (66,*) 'natom must be defined'
                  stop
               endif
            enddo
            read (25,*) natom2,natomt2
            close (25)
         endif
         if(iadd.eq.1.or.iabs.eq.1)then
            natom = natom1+natom2
         else if (iiso.eq.1.or.ibeta.eq.1) then
            natom = natom1
         endif

         if (iadd.eq.1) natomt = natomt1+natomt2
         if (iabs.eq.1) natomt = natomt1+natomt2+1
         if (iiso.eq.1) natomt = natomt1
         if (ibeta.eq.1) natomt = natomt1

         do while (WORD.NE.'CHARGE')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'charge and spin must be defined'
               stop
            endif
         enddo
         read (15,*) icharge,ispin
c         rewind(15)
         close(15)
c molpro z-mat file input from gaussian z-mat output

         if (code_name.eq.'MOLPRO')then   
            command1='cp -f ./data/hl_molpro.dat ./hl_molpro.dat'       
            call commrun(command1)
            call commrun(commandcopy)
            read (17,*)
            if(ispecies.eq.100.or.ispecies.eq.101)then
               if(iabs.eq.1) then
                  ntot=natomt-1
c               ntot=natom
               else
                  ntot=natomt
               endif
            else
               ntot=natomt
            endif
            do iatom = 1 , ntot
               read (17,'(A60)') atomlabel(iatom)
c               write (67,*)atomlabel(iatom)
            enddo
c            write (67,*)'}'
         else if (code_name.eq.'G09')then   
            open (unit=21,file='./data/theory.dat',status='unknown')
            do while (WORD.NE.WORD_HL)
               call LineRead (21)
               if (WORD.EQ.'END') then
                  write (66,*) ' HL g09 theory must be defined'
                  write (66,*) ' in file theory.dat '
                  stop
               endif
            enddo 
            read (21,'(A70)') comline1
            read (21,'(A70)') comline2
            close(21)
            read (17,*)
            if (idebug.ge.2) write (6,*) ' starting gaussian input'
            do iatom = 1 , natomt
               read (17,'(A60)') atomlabel(iatom)
            enddo
         endif

cc read coordinate names
         ncoord = 3*natom-6
         if (natom.eq.1) ncoord=0
         if (natom.eq.2) ncoord=1
         if((ispecies.ne.100).and.(ispecies.ne.101))then
            do iint = 1 , ncoord
               read (17,*) intcoor(iint),xint(iint)
c     write (6,*) 'intcoor, xint ',intcoori(iint),xinti(iint)
            enddo
         endif
         close (unit=17,status='keep')
      endif


      ircons=0
      ismp=0
      ifreq=0
      ilev=2
      ircons=0
      ixyz=0
      ired=0
      ntau=0
      if(iabs.eq.1) ireact=natom1+1
      if(iadd.eq.1) ireact=natom1
      if((ispecies.eq.100).or.(ispecies.eq.101))ixyz=1

cc now run the high level calculation

      if(ilev1code.eq.1) then
         call g09fopt(tau,ntau,natom,natomt,numproc,gmem,
     $ coord,vtot_0,vtotr,freq,ifreq,ilin,ismp,comline1,
     $ comline2,icharge,ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ired)

         vtotr=vtot_0

      else if (ilev1code.eq.2) then
         numproc=numprochl
         
         call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $ coord,vtotr,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ilev
     $ ,ispecies)

      endif

      open (unit=65,file='temp.log',status='unknown')
      rewind (65)
      write (65,*) vtotr
      close(65)

c        open (unit=99,status='unknown')
c        rewind (99)
c        write (99,1021)
c        rewind (99)
c        read (99,1020) command1
c        close (99)
c        call commrun(command1)
c1020    FORMAT (A100)
c1021    format ("sed -ie 's/D/E/g' temp.log")       
      

cc now save energy in me_files

      if (ispecies.eq.1) then
         command1='cp -f temp.log ./me_files/reac1_en.me'
      else if (ispecies.eq.2) then
         command1='cp -f temp.log ./me_files/reac2_en.me'
      else if (ispecies.eq.3) then
         command1='cp -f temp.log ./me_files/prod1_en.me'
      else if (ispecies.eq.4) then
         command1='cp -f temp.log ./me_files/prod2_en.me'
      else if (ispecies.eq.5.or.ispecies.eq.51) then
         command1='cp -f temp.log ./me_files/wellr_en.me'
      else if (ispecies.eq.6.or.ispecies.eq.61) then
         command1='cp -f temp.log ./me_files/wellp_en.me'
      else if (ispecies.eq.0) then
         command1='cp -f temp.log ./me_files/ts_en.me'
      else if (ispecies.eq.100) then
         command1='cp -f temp.log ./me_files/ircst_en.me'
      else if (ispecies.eq.101) then
         command1='cp -f temp.log ./me_files/ircarr_en.me'
      endif
      call commrun(command1)

cc now save the output log files in hl_logs directory
cc this is code dependent

      if (code_name.eq.'MOLPRO')then   
         open (unit=99,status='unknown')
         rewind (99)
         if (ispecies.eq.1) then
            command1='cp -f molpro.out ./hl_logs/reac1_molpro.out'
         endif
         if (ispecies.eq.2) then
            command1='cp -f molpro.out ./hl_logs/reac2_molpro.out'
         endif
         if (ispecies.eq.3) then
            command1='cp -f molpro.out ./hl_logs/prod1_molpro.out'
         endif
         if (ispecies.eq.4) then
            command1='cp -f molpro.out ./hl_logs/prod2_molpro.out'
         endif
         if (ispecies.eq.5.or.ispecies.eq.51) then
            command1='cp -f molpro.out ./hl_logs/wellr_molpro.out'
         endif
         if (ispecies.eq.6.or.ispecies.eq.61) then
            command1='cp -f molpro.out ./hl_logs/wellp_molpro.out'
         endif
         if (ispecies.eq.0) then
            command1='cp -f molpro.out ./hl_logs/ts_molpro.out'
         endif
         if (ispecies.eq.100) then
            command1='cp -f molpro.out ./hl_logs/ircst_molpro.out'
         endif
         if (ispecies.eq.101) then
            command1='cp -f molpro.out ./hl_logs/ircarr_molpro.out'
         endif
         call commrun(command1)
      endif

      if (code_name.eq.'G09')then   

         if (ispecies.eq.1) then
            command1='cp -f geom.log ./hl_logs/reac1_g09.out'
         else if (ispecies.eq.2) then
            command1='cp -f geom.log ./hl_logs/reac2_g09.out'
         else if (ispecies.eq.3) then
            command1='cp -f geom.log ./hl_logs/prod1_g09.out'
         else if (ispecies.eq.4) then
            command1='cp -f geom.log ./hl_logs/prod2_g09.out'
         else if (ispecies.eq.5.or.ispecies.eq.51) then
            command1='cp -f geom.log ./hl_logs/wellr_g09.out'
         else if (ispecies.eq.6.or.ispecies.eq.61) then
            command1='cp -f geom.log ./hl_logs/wellp_g09.out'
         else if (ispecies.eq.0) then
            command1='cp -f geom.log ./hl_logs/ts_g09.out'
         else if (ispecies.eq.100) then
            command1='cp -f geom.log ./hl_logs/ircst_g09.out'
         else if (ispecies.eq.101) then
            command1='cp -f geom.log ./hl_logs/ircarr_g09.out'
         endif
         call commrun(command1)
      endif
c     close (unit=65,status='keep')

 999  continue

      close (unit=66,status='keep')


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ktp

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      logical leof,lsec,ltit

      dimension freq(3*natommx)
      dimension gelec1(nelecmx),eelec1(nelecmx)
      dimension gelec2(nelecmx),eelec2(nelecmx)
      dimension gelec12(nelecmx),eelec12(nelecmx)
      dimension rotpot(noptmx)

      character*300 command1
      character*80 cread
      character*20 cread1,cread2
      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      character*30 gmem
      character*20 stoich_well
      character*20 filewell

      include 'filcomm.f'


cc before doing anything check and modify me files to have structures consistent with calculation type

cc check on products
      if(ip1.eq.1.and.ip2.eq.1)then
         numend=0
         command1='egrep End ./me_files/prod1_fr.me > tmp.log'
         call commrun(command1)
         command1='wc tmp.log > tmp1.log'
         call commrun(command1)
         open(unit=99,file='tmp1.log',status='unknown')
         read(99,*)numend
         close(99)
         if(numend.eq.2)then
            write(command1,5001)
            call commrun(command1)
         endif
         numend=0
         command1='egrep End ./me_files/prod2_fr.me > tmp.log'
         call commrun(command1)
         command1='wc tmp.log > tmp1.log'
         call commrun(command1)
         open(unit=99,file='tmp1.log',status='unknown')
         read(99,*)numend
         close(99)
         if(numend.eq.3)then
            write(command1,5002)
            call commrun(command1)
         endif
         write(command1,5004) 
         call commrun(command1)
         write(command1,5007) 
         call commrun(command1)
         write(command1,5008) 
         call commrun(command1)

 5001    format(" sed -ie '0,/End/{s/End/ /}' me_files/prod1_fr.me")
 5002    format(" sed -ie '0,/End/{s/End/ /}' me_files/prod2_fr.me")
c         stop
 5004    format(" sed -ie 's/Species/Fragment Prod2/g'
     +    me_files/prod2_ge.me")

 5007    format(" sed -ie 's/Well/Bimolecular/g'
     +    me_files/prod1_ge.me")
 5008    format(" sed -ie 's/Species/Fragment PROD1/g'
     +    me_files/prod1_ge.me")

      endif

cc check on single product

      if(ipr1.eq.1)then
         numend=0
         command1='egrep End ./me_files/prod1_fr.me > tmp.log'
         call commrun(command1)
         command1='wc tmp.log > tmp1.log'
         call commrun(command1)
         open(unit=99,file='tmp1.log',status='unknown')
         read(99,*)numend
         close(99)
         if(numend.eq.1)then
cc add a line at the end fo the file
            open(unit=99,file='tmp.log',status='unknown')
            write(99,*)'End'
            close(99)
            command1='cat ./me_files/prod1_fr.me tmp.log > tmp1.log'
            call commrun(command1)
            command1='cp -f tmp1.log ./me_files/prod1_fr.me'
            call commrun(command1)
         endif
         write(command1,5005)
         call commrun(command1)
         write(command1,5006)
         call commrun(command1)
         write(command1,5009)
         call commrun(command1)
         write(command1,5010)
         call commrun(command1)

      endif
 5003 format(" sed -ie '0,/End/{s/End/ /}' me_files/prod2_fr.me")
 5005 format(" sed -ie 's/Bimolecular/Well/g'
     +    me_files/prod1_ge.me")
 5006 format(" sed -ie 's/Fragment PROD1/Species/g'
     +    me_files/prod1_ge.me")
 5009 format(" sed -ie 's/TS REACS WP/TS REACS PRODS/g'
     +    me_files/ts_ge.me")
 5010 format(" sed -ie 's/TS REACS WP/TS REACS PRODS/g'
     +    me_files/variational.me")


c      endif

cc this is to be continued to solve for reac1 and reac2 ambiguities when imported
cc from previous calculations

c      stop


cc now open files with energy and estimate forward and backward energy barriers

      open (unit=107, file='./me_files/reac1_en.me', status='old')
      read (107,*) reac1_en
      close (107)
      open (unit=107, file='./me_files/reac1_zpe.me', status='old')
      read (107,*) reac1_zpe
      close (107)
      reac1_en=reac1_en+reac1_zpe

cc here we check if we are producing only the reactant block
      ireactonly=0
      if(iabs.eq.0.and.iadd.eq.0.and.iiso.eq.0.and.ibeta.eq.0)then
         ireactonly=1
      endif

      if(iabs.eq.1.or.iadd.eq.1)then
         open (unit=107, file='./me_files/reac2_en.me', status='old')
         read (107,*) reac2_en
         close (107)
         open (unit=107, file='./me_files/reac2_zpe.me', status='old')
         read (107,*) reac2_zpe
         close (107)
         reac2_en=reac2_en+reac2_zpe
      else if (iiso.eq.1.or.ibeta.eq.1) then
         reac2_en=0
      else
         reac2_en=0
      endif

      reac_en=reac1_en+reac2_en


cc here we copy to the block directory all files with *dhr*.dat 
cc in the me_files directory

      command1='cp -f ./me_files/*dhr*.dat ./me_blocks'
      call commrun(command1)

cc save absolute energy of reactants in me blocks

      if(iadd.eq.1.or.iabs.eq.1)then
         open (unit=107,file='./me_blocks/bimol1_en.me',
     $         status='unknown')  
         write(107,*)reac_en
         close(107)

         if (idebug.ge.2) write (6,*) 'ip test',ip1,ip2
         if ((ip1.eq.1).and.(ip2.eq.1)) then
            open (unit=107, file='./me_files/prod1_en.me', status='old')
            read (107,*) prod1_en
            close (107)
            open (unit=107, file='./me_files/prod1_zpe.me',status='old')
            read (107,*) prod1_zpe
            close (107)
            prod1_en=prod1_en+prod1_zpe
            
            open (unit=107, file='./me_files/prod2_en.me', status='old')
            read (107,*) prod2_en
            close (107)
            open (unit=107, file='./me_files/prod2_zpe.me',status='old')
            read (107,*) prod2_zpe
            close (107)
            prod2_en=prod2_en+prod2_zpe
            prod_en=prod1_en+prod2_en-reac_en
            prod_en_tot=prod1_en+prod2_en
            open (unit=107,file='./me_blocks/bimol2_en.me',
     $             status='unknown')  
            write(107,*) prod_en_tot
            close(107)
         else if (ipr1.eq.1)then
            open (unit=107, file='./me_files/prod1_en.me', status='old')
            read (107,*) prod1_en
            close (107)
            open (unit=107, file='./me_files/prod1_zpe.me',status='old')
            read (107,*) prod1_zpe
            close (107)
            prod1_en=prod1_en+prod1_zpe
            prod_en=prod1_en-reac_en
            prod_en_tot=prod1_en
         else
            prod_en=0.0
         endif
c      else if (iiso.eq.1)then
c         open (unit=107,file='./me_blocks/reac_iso_en.me',
c     $         status='unknown')  
c         write(107,*)reac_en
c         close(107)
c         prod_en=0.0
c         if ((ip1.eq.1).and.(ip2.eq.1)) then
c            write(*,*) ' the prods keyword is not'
c            write(*,*) ' compatible with the iso keyword '
c            stop
c         endif
      else if (ibeta.eq.1.or.iiso.eq.1)then
         open (unit=107,file='./me_blocks/reac_ibeta_en.me',
     $         status='unknown')  
         write(107,*)reac_en
         close(107)
         if (idebug.ge.2) write (6,*) 'ip test',ip1,ip2
         if ((ip1.eq.1).and.(ip2.eq.1)) then
            open (unit=107, file='./me_files/prod1_en.me', status='old')
            read (107,*) prod1_en
            close (107)
            open (unit=107, file='./me_files/prod1_zpe.me',status='old')
            read (107,*) prod1_zpe
            close (107)
            prod1_en=prod1_en+prod1_zpe
            
            open (unit=107, file='./me_files/prod2_en.me', status='old')
            read (107,*) prod2_en
            close (107)
            open (unit=107, file='./me_files/prod2_zpe.me',status='old')
            read (107,*) prod2_zpe
            close (107)
            prod2_en=prod2_en+prod2_zpe
            prod_en=prod1_en+prod2_en-reac_en
            prod_en_tot=prod1_en+prod2_en
            open (unit=107,file='./me_blocks/bimol2_en.me',
     $             status='unknown')  
            write(107,*) prod_en_tot
            close(107)
         else if (ipr1.eq.1)then
            open (unit=107, file='./me_files/prod1_en.me', status='old')
            read (107,*) prod1_en
            close (107)
            open (unit=107, file='./me_files/prod1_zpe.me',status='old')
            read (107,*) prod1_zpe
            close (107)
            prod1_en=prod1_en+prod1_zpe
            prod_en=prod1_en-reac_en
            prod_en_tot=prod1_en
         else
            prod_en=0.0
         endif
      else
cc if no reaction is specified, only the reactant block is written
         open (unit=107,file='./me_blocks/reac_en.me',
     $         status='unknown')  
         write(107,*)reac_en
         close(107)
      endif

      if(ireactonly.eq.1) goto 1010
c     write (6,*) 'prod en test',reac_en,prod1_en,prod2_en,prod_en

      if (irw.eq.1) then
         if(iiso.eq.1.or.ibeta.eq.1)then
            write(*,*)  ' the wellr keyword is not'
            write(*,*) ' compatible with the iso or beta keywords '
            stop
         endif
         open (unit=107, file='./me_files/wellr_en.me', status='old')
         read (107,*) wellr_en
         close (107)
         open (unit=107, file='./me_files/wellr_zpe.me', status='old')
         read (107,*) wellr_zpe
         close (107)
         wellr_el_en=wellr_en
         wellr_en=wellr_en+wellr_zpe
         wellr_en=wellr_en-reac_en
      else
         wellr_en=0.0
      endif


c      open (unit=107,file='./me_blocks/w1_en.me',status='unknown')  
      if (irw.eq.1) then
c         write(107,*)wellr_en+reac_en
         w1_en=wellr_en+reac_en
      else
c         write(107,*)'0.'
         w1_en=0.
      endif
c      close(107)

      if (ipw.eq.1) then
         open (unit=107, file='./me_files/wellp_en.me', status='old')
         read (107,*) wellp_en
         close (107)
         open (unit=107, file='./me_files/wellp_zpe.me', status='old')
         read (107,*) wellp_zpe
         close (107)
         wellp_el_en=wellp_en
         wellp_en=wellp_en+wellp_zpe
         wellp_en=wellp_en-reac_en
      else
         wellp_en=0.0
      endif
      write (6,*) 'ipw 2 test',ipw,wellp_en


c      open (unit=107,file='./me_blocks/w2_en.me',status='unknown')  
      if(ipw.eq.1)then
c         write(107,*)wellp_en+reac_en
         w2_en=wellp_en+reac_en
      else
c         write(107,*)'0'
         w2_en=0.
      endif
c      close(107)
      write (6,*) 'ipw 3 test',ipw,w2_en


      open (unit=107, file='./me_files/ts_en.me', status='old')
      read (107,*) ts_en
      close (107)
      open (unit=107, file='./me_files/ts_zpe.me', status='old')
      read (107,*) ts_zpe
      close (107)
      ts_el_en=ts_en
      ts_en=ts_en+ts_zpe
      ts_en_au=ts_en

c      open (unit=107,file='./me_blocks/ts0_en.me',status='unknown')  
c      write(107,*)ts_en
c      close(107)

      ts_en=ts_en-reac_en


      if (idebug.ge.2) write (6,*) 'en test',reac_en,prod_en,wellr_en,
     $  wellp_en,ts_en


      if (nts.eq.1) then
         if (iabs.eq.1) then 
            if(irw.eq.1) then
               wdepfor = (ts_en-wellr_en)*cautokcal
            else
               wdepfor = ts_en*cautokcal
            endif
            if(ipw.eq.1) then
               wdepback = (ts_en-wellp_en)*cautokcal
            else
               wdepback = (ts_en-prod_en)*cautokcal
            endif
            if (idebug.ge.2) write (6,*) 'wdep abs test',wdepfor,
     $       wdepback,ts_en,
     $       prod_en,reac_en
         endif
         if (iadd.eq.1) then 
            if(irw.eq.1) then
               wdepfor = (ts_en-wellr_en)*cautokcal
            else
               wdepfor = ts_en*cautokcal
            endif
            if(ipw.eq.1) then
               wdepback = (ts_en-wellp_en)*cautokcal
            else
               wdepback = (ts_en-prod_en)*cautokcal
            endif
            if (idebug.ge.2) write (6,*) 'wdep add test',wdepfor,
     $       wdepback,ts_en,
     $       wellp_en,reac_en
         endif
         if (iiso.eq.1) then
            wdepfor = ts_en*cautokcal
            wdepback = (ts_en-wellp_en)*cautokcal
            if ((ip1.eq.1).and.(ip2.eq.1)) then
               wdepback = (ts_en-prod_en)*cautokcal
               write(*,*)'wdep back is ',wdepback
               write(*,*)'ts en is ',ts_en
               write(*,*)'prod en is ',prod_en_tot
            else if (ipr1.eq.1)then
               wdepback = (ts_en-prod_en)*cautokcal
               write(*,*)'wdep back is ',wdepback
               write(*,*)'ts en is ',ts_en
               write(*,*)'prod en is ',prod_en_tot
            endif

            if(ipw.eq.1) then
               wdepback = (ts_en-wellp_en)*cautokcal
            endif
            if (idebug.ge.2) write (6,*) 'wdep add test',wdepfor,
     $           wdepback,ts_en,
     $           wellp_en,reac_en
         endif

         if (ibeta.eq.1) then
            wdepfor = ts_en*cautokcal
            wdepback = (ts_en-prod_en)*cautokcal
            if (idebug.ge.2) write (6,*) 'wdep add test',wdepfor,
     $           wdepback,ts_en,
     $           prod_en,reac_en
         endif

      endif
      if (nts.eq.2) then
         if (iabs.eq.1) then 
            wdepfor = (ts_en-wellr_en)*cautokcal
            wdepback = (ts_en-prod_en)*cautokcal
         endif
         if (iadd.eq.1) then 
            wdepfor = (ts_en-wellr_en)*cautokcal
            wdepback = (ts_en-wellp_en)*cautokcal
         endif
         if (iiso.eq.1.or.ibeta.eq.1) then 
            write(*,*)'only 1 TS supported for iso or beta reactions'
            stop
         endif
      endif
      if (nts.eq.3) then
         if (iabs.eq.1) then 
            wdepfor = (ts_en-wellr_en)*cautokcal
            wdepback = (ts_en-wellp_en)*cautokcal
         endif
         if (iadd.eq.1) then 
            wdepfor = (ts_en-wellr_en)*cautokcal
            wdepback = (ts_en-wellp_en)*cautokcal
         endif
         if (iiso.eq.1.or.ibeta.eq.1) then 
            write(*,*)'only 1 TS supported for iso or beta reactions'
            stop
         endif
      endif
      prod_en = prod_en*cautokcal
      wellr_en = wellr_en*cautokcal
      wellp_en = wellp_en*cautokcal
      ts_en = ts_en*cautokcal
      write (6,*) 'wellp_en test',wellp_en

cc it is possible that the electronic energy of the wells is smaller than
cc that of the TS, but that the the ZPE correction makes is higher. In this
cc case the eckart correction is set to 0

      check1=0.1
      check2=0.1
      if(irw.eq.1) then
         check1=ts_el_en-wellr_el_en
         check2=ts_el_en+ts_zpe-wellr_el_en - wellr_zpe
      endif
      if(ipw.eq.1) then
         check1=ts_el_en-wellp_el_en
         check2=ts_el_en+ts_zpe-wellp_el_en - wellp_zpe
      endif
      if(check1.gt.0.and.check2.lt.0) then
         write (6,*) 'setting Eckart contribution to negligible'
         write (6,*) 'since ZPE makes barrier height negative'
         write (*,*) 'setting Eckart contribution to negligible'
         write (*,*) 'since ZPE makes barrier height negative'
         wdepfor=0.01
         wdepback=0.01
      endif
      if(inotunn.eq.1) then
         write (6,*) 'setting Eckart contribution to negligible'
         write (6,*) 'as requested in input'
         write (*,*) 'setting Eckart contribution to negligible'
         write (*,*) 'as requested in input'
         wdepfor=0.01
         wdepback=0.01
      endif


cc now we seed the TS file of me produced by the HL subroutine with these data

      open (unit=99,status='unknown')
      if(ivar.eq.0) then
         write (99,999)
         rewind (99)
         read (99,1120) command1

 999     format ("cat ./me_files/ts_ge.me ./me_files/ts_hr.me 
     +  ./me_files/ts_fr.me > temp0.me ")

         call commrun(command1)

         rewind (99)
         write (99,1000) ts_en
         rewind (99)
         read (99,1120) command1
 1000    format(" sed -e 's/\$tsen/",1X,F7.2,1X,
     +    "/' temp0.me > temp1.me")

         call commrun(command1)
      else if (ivar.eq.1) then
         command1='cp -f ./me_files/variational.me ./temp1.me' 
         call commrun(command1)
      endif

      rewind (99)
      write (99,1001) wdepfor
      rewind (99)
      read (99,1120) command1
 1001 format (" sed -e 's/\$wdepfor/",1X,F7.2,1X,
     +    "/' ./temp1.me > temp2.me")

      call commrun(command1)

      rewind (99)
      write (99,1002) wdepback
      rewind (99)
      read (99,1120) command1
 1002 format (" sed -e 's/\$wdepback/",1X,F7.2,1X,
     +    "/' ./temp2.me > temp21.me")

      call commrun(command1)

cc now save the TS or variational input as blocks

      command1='echo End > temp211.me'
      call commrun(command1)

      if(ivar.eq.1)then
         command1='cat temp21.me temp211.me > temp21a.me'
         call commrun(command1)
      else if(ivar.eq.0)then
         command1='cp -f temp21.me temp21a.me'
         call commrun(command1)
      endif

      if(ivar.eq.0) then
         rewind (99)
         write (99,1003) 
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
 1003    format(" sed -ie 's/Barrier TS REACS WR/Barrier TS0/g'
     +    temp21a.me")
         rewind (99)
c         write(99,2003)ts_en_au
         write(99,2003)ts_en
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
 2003    format ("sed -i 's/ZeroEnergy.*/ZeroEnergy[kcal/mol]    ",
     $    f11.5," /' temp21a.me")
         command1='cp -f temp21a.me ./me_blocks/ts1.me'
         call commrun(command1)
         open (unit=107,file='./me_blocks/ts1_en.me',status='unknown')  
         write(107,*)ts_en_au
         close(107)
      else if (ivar.eq.1) then
cc save au value of maximum of variational energy
         open (unit=107,file='./me_files/de1_TSvar.me',status='unknown')  
         read(107,*)ts_en_var
         close(107)
         open (unit=107,file='./me_blocks/ts1_en.me',status='unknown')  
         write(107,*)ts_en_au+ts_en_var/CAUTOKCAL
         close(107)
         rewind (99)
         write (99,1005) 
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
 1005    format(" sed -ie 's/Barrier TS REACS WR/Barrier VARIATIONAL/g'
     +    temp21a.me")
         command1='cp -f temp21a.me ./me_blocks/variational.me'
         call commrun(command1)
      endif
      close(99)
cc now add  two lines to the TS file

      command1='echo End > temp211.me'
      call commrun(command1)
      command1='echo End > temp212.me'
      call commrun(command1)
      command1='cat temp211.me temp212.me > temp213.me'
      call commrun(command1)
      command1='cat temp21.me temp213.me > temp3.me'
      call commrun(command1)

c now formulate separate reactants, products, well, and PST TS files

 1010 continue

      open (unit=99,status='unknown')
      if(iabs.eq.1.or.iadd.eq.1)then
         write (99,150)
      else if (iiso.eq.1.or.ibeta.eq.1) then
         write (99,149)
      else
         write (99,149)
      endif

 150  format ("cat ./me_files/reac1_ge.me 
     +  ./me_files/reac1_hr.me ./me_files/reac1_fr.me 
     +  ./me_files/reac2_ge.me ./me_files/reac2_hr.me 
     +  ./me_files/reac2_fr.me > reactants.me")

 149  format ("cat ./me_files/reac1_ge.me 
     +  ./me_files/reac1_hr.me ./me_files/reac1_fr.me
     +  > reactants.me")

      rewind (99)
      read (99,1120) command1
      call commrun(command1)

cc save the reactants block as bimol1 for abstraction
cc and addition and as reac1_iso for isomerization

      if(iabs.eq.1.or.iadd.eq.1)then
         command1='cp -f reactants.me bimol1.me'
         call commrun(command1)
         rewind (99)
         write(99,151)
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
 151     format ("sed -ie 's/REACS/BIMOL1/g' bimol1.me")
         command1='cp -f bimol1.me ./me_blocks'
         call commrun(command1)
      else if (iiso.eq.1)then
         command1='cp -f reactants.me reac1_iso.me'
         call commrun(command1)
         rewind (99)
         write(99,161)
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
 161     format ("sed -ie 's/REACS/REAC1/g' reac1_iso.me")
         command1='cp -f reac1_iso.me ./me_blocks'
         call commrun(command1)
      else if (ibeta.eq.1)then
         command1='cp -f reactants.me reac1_beta.me'
         call commrun(command1)
         rewind (99)
         write(99,162)
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
 162     format ("sed -ie 's/REACS/REAC1/g' reac1_beta.me")
         command1='cp -f reac1_beta.me ./me_blocks'
         call commrun(command1)
      else
         command1='cp -f reactants.me reac1.me'
         call commrun(command1)
         rewind (99)
         write(99,163)
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
 163     format ("sed -ie 's/REACS/REAC1/g' reac1.me")
         command1='cp -f reac1.me ./me_blocks'
         call commrun(command1)
      endif

cc if no reactions is specified, we go the the end of the sub and exit

      if(ireactonly.eq.1) then
         close(99)
         goto 1020
      endif

cc create the products block
      
      if ((ip1.eq.1).and.(ip2.eq.1)) then
         rewind (99)
         write (99,250)
 250     format ("cat ./me_files/prod1_ge.me 
     +  ./me_files/prod1_hr.me ./me_files/prod1_fr.me 
     +  ./me_files/prod2_ge.me ./me_files/prod2_hr.me 
     +  ./me_files/prod2_fr.me > products0.me")
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
         rewind (99)
         write (99,260) prod_en
 260     format(" sed -e 's/\$proden/",1X,F7.2,1X,
     +        "/' products0.me > products.me")
         rewind (99)
         read (99,1120) command1
         call commrun(command1)

cc save the products block as bimol2

         command1='cp -f products.me bimol2.me'
         call commrun(command1)
         rewind (99)
         write(99,152)
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
 152     format ("sed -ie 's/PRODS/BIMOL2/g' bimol2.me")
         command1='cp -f bimol2.me ./me_blocks'
         call commrun(command1)
      else if (ipr1.eq.1) then
         rewind (99)
         write (99,1250)
 1250     format ("cat ./me_files/prod1_ge.me 
     +  ./me_files/prod1_hr.me ./me_files/prod1_fr.me 
     +   > products0.me")
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
         rewind (99)
         write (99,1261) 
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
         rewind (99)
         write (99,1260) prod_en 
         rewind (99)
         read (99,1120) command1
         call commrun(command1)

 1260     format(" sed -e 's/\$proden/",1X,F7.2,1X,
     +        "/' products0.me > products.me")
 1261    format(" sed -ie 's/ZeroEnergy\[kcal\/mol\]             0./
     + ZeroEnergy[kcal\/mol] \$proden/' products0.me ")

cc save the products block as prod1_iso 

         command1='cp -f products.me ./me_blocks/prod1_iso.me'
         call commrun(command1)
c         stop
      endif

      rewind (99)
      if (irw.eq.0) then
         write (99,311)
      else
         write (99,312)
      endif
 311  format ("cat ./me_files/wellr_fake.me 
     +   > wellr.me")
 312  format ("cat ./me_files/wellr_ge.me 
     +  ./me_files/wellr_hr.me ./me_files/wellr_fr.me 
     +  > wellr0.me")

      rewind (99)
      read (99,1120) command1
      call commrun(command1)

      if(irw.ne.0)then
         rewind (99)
         write (99,320) wellr_en
 320     format(" sed -e 's/\$wellren/",1X,F7.2,1X,
     +    "/' wellr0.me > wellr.me")

         rewind (99)
         read (99,1120) command1
         call commrun(command1)
      endif

cc save wellr if not fake

      call LineRead(0)
      if(irw.ne.0)then
         command1='cp -f wellr.me wellr1.me'
         call commrun(command1)
         rewind (99)
         write(99,153)
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
 153     format ("sed -ie 's/WR/W1/g' wellr1.me")
         rewind (99)
         write(99,1531)w1_en
         write(*,*)'w1_en is',w1_en
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
 1531    format ("sed -i 's/ZeroEnergy.*/ZeroEnergy[au]    ",
     $    f11.5," /' wellr1.me")
         command1='cp -f wellr1.me ./me_blocks/w1.me'
         call commrun(command1)
      endif


      rewind (99)
      if (ipw.eq.0) then
         write (99,331)
      else
         write (99,332)
      endif
 331  format ("cat ./me_files/wellp_fake.me 
     +   > wellp.me")
 332  format ("cat ./me_files/wellp_ge.me 
     +  ./me_files/wellp_hr.me ./me_files/wellp_fr.me 
     +  > wellp0.me")

      rewind (99)
      read (99,1120) command1
      call commrun(command1)

      if (ipw.ne.0) then
         write (6,*) 'wellp_en test',wellp_en
         rewind (99)
         write (99,340) wellp_en
         rewind (99)
         read (99,1120) command1
 340     format(" sed -e 's/\$wellpen/",1X,F7.2,1X,
     +    "/' wellp0.me > wellp.me")

         rewind (99)
         read (99,1120) command1
         call commrun(command1)
      endif

cc save wellp if not fake

      if (ipw.ne.0)then
         command1='cp -f wellp.me wellp0.me'
         call commrun(command1)
         rewind (99)
         write(99,154)
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
 154     format ("sed -ie 's/WP/W2/g' wellp0.me")
         rewind (99)
         write(99,1541)w2_en
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
 1541    format ("sed -i 's/ZeroEnergy.*/ZeroEnergy[au]    ",
     $    f11.5," /' wellp0.me")
         command1='cp -f wellp0.me ./me_blocks/w2.me'
         call commrun(command1)
      endif

      rewind (99)
      if(ipw.ne.0)then
         write (99,350)
      else
         write (99,351)
      endif
 350  format ("cat wellr.me wellp.me > wells.me")
 351  format ("cat wellr.me > wells.me")
      rewind (99)
      read (99,1120) command1
      call commrun(command1)

cc now create phase space theory blocks

cc start with reactant block

      if(iadd.eq.1.or.iabs.eq.1)then
         command1='echo 0 > temp1.log'
         call commrun(command1)
         rewind (99)
         write (99,371)
 371     format (" egrep Freq me_files/reac1_fr.me|awk '{print $2}'
     $   > temp.log")
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
         command1='cat temp.log temp1.log > temp2.log'
         call commrun(command1)
         open(unit=123,file='temp2.log',status='unknown')
         read(123,*)nfreqr1
         close(123)

         command1='echo 0 > temp1.log'
         call commrun(command1)
         rewind (99)
         write (99,372)
 372     format (" egrep Freq me_files/reac2_fr.me|awk '{print $2}'
     $> temp.log")
         rewind (99)
         read (99,1120) command1
         call commrun(command1)
         command1='cat temp.log temp1.log > temp2.log'
         call commrun(command1)
         open(unit=123,file='temp2.log',status='unknown')
         read(123,*)nfreqr2
         close(123)
         nfreqr=nfreqr1+nfreqr2

c     determine stoichiometry of wellr
         filewell='wellr.me'
         call stoichiometry(filewell,stoich_well)
cc
c         write(*,*)'stoich_wellr is ',stoich_wellr
         open(unit=123,file='pstr1.me',status='unknown')
         write (123,*) 'Barrier	B1  REACS WR'
         write (123,*) '  RRHO '
         write (123,*) '    Stoichiometry  ',stoich_well
         write (123,*) '    Core PhaseSpaceTheory'
         open(unit=124,file='me_files/reac1_ge.me',status='unknown')
         read(124,*)
         read(124,*)
         read(124,*)
         read(124,*)cread,natom
         write (123,*) '      FragmentGeometry[angstrom]    ',natom
         do j=1,natom
            read(124,'(A80)')cread
            write(123,'(A80)')cread
         enddo
         read(124,*)
         read(124,*)cread,symr1
         close(124)
cc now write geometry of second fragment
         open(unit=124,file='me_files/reac2_ge.me',status='unknown')
         read(124,*)
         read(124,*)cread
         if(cread.eq.'Atom')then
            write (123,*) '      FragmentGeometry[angstrom]    ',1
            read(124,*)cread1,cread2
            write(123,*)cread2,' 0. 0. 0.'
            symr2=1
         else
c         read(124,*)
            read(124,*)cread,natom
            write (123,*) '      FragmentGeometry[angstrom]    ',natom
            do j=1,natom
               read(124,'(A80)')cread
               write(123,'(A80)')cread
            enddo
            read(124,*)
            read(124,*)cread,symr2
         endif
         close(124)

         symr=symr1*symr2

         write (123,*) '		SymmetryFactor  ',symr
         write (123,*) '          PotentialPrefactor[au] 	10.' 
         write (123,*) '          PotentialPowerExponent 	6' 
         write (123,*) '	    End'

cc now write hindered rotor section
         call LineRead (0)
         open(unit=19,file='data/reac1.dat',status='unknown')
         do while (WORD.NE.'NHIND')
            call LineRead (19)
            if (WORD.EQ.'END') then
               write (6,*) ' error in reading input in sub ktp'
               stop
            endif
         enddo
         read (19,*) nhind
         read(19,*)
         if(nhind.ne.0)then
            do j=1,nhind
               open(unit=20,file='me_files/reac1_hr.me',
     $              status='unknown')
               read(20,'(A80)')cread
               write(123,'(A80)')cread
               open(unit=21,file='me_files/reac1_ge.me',
     $              status='unknown')
               read(21,*)
               read(21,*)
               read(21,*)
               read(21,*)cread,natom
               write(123,*)'Geometry[angstrom]    ',natom
               do i=1,natom
                  read(21,'(A80)')cread
                  write(123,'(A80)')cread
               enddo
               close(21)
               read(20,'(A80)')cread
               write(123,'(A80)')cread
               read(20,'(A80)')cread
               write(123,'(A80)')cread
               read(20,'(A80)')cread
               write(123,'(A80)')cread
               read(20,'(A80)')cread
               write(123,'(A80)')cread

               read(19,*)cread,start,arr,nptot
               read(20,*)(rotpot(ij),ij=1,nptot)
               write(123,8010)(rotpot(ij),ij=1,nptot)
               read(20,'(A80)')cread
               write(123,'(A80)')cread
            enddo
         endif
         close(19)
         close(20)

         call LineRead (0)
         open(unit=19,file='data/reac2.dat',status='unknown')
         do while (WORD.NE.'NHIND')
            call LineRead (19)
            if (WORD.EQ.'END') then
               write (6,*) ' error in reading input in sub ktp'
               stop
            endif
         enddo
         read (19,*) nhind
         read(19,*)
         if(nhind.ne.0)then
            do j=1,nhind
               open(unit=20,file='me_files/reac2_hr.me',
     $             status='unknown')
               read(20,'(A80)')cread
               write(123,'(A80)')cread
               open(unit=21,file='me_files/reac2_ge.me',
     $             status='unknown')
c            read(21,*)
               read(21,*)
               read(21,*)
               read(21,*)cread,natom
               write(123,*)'Geometry[angstrom]    ',natom
               do i=1,natom
                  read(21,'(A80)')cread
                  write(123,'(A80)')cread
               enddo
               close(21)
               read(20,'(A80)')cread
               write(123,'(A80)')cread
               read(20,'(A80)')cread
               write(123,'(A80)')cread
               read(20,'(A80)')cread
               write(123,'(A80)')cread
               read(20,'(A80)')cread
               write(123,'(A80)')cread
               read(19,*)cread,start,arr,nptot
               read(20,*)(rotpot(ij),ij=1,nptot)
               write(123,8010)(rotpot(ij),ij=1,nptot)
               read(20,'(A80)')cread
               write(123,'(A80)')cread
            enddo
         endif
         close(19)
         close(20)

         write (123,*) '    Frequencies[1/cm] ',nfreqr

         open(unit=124,file='me_files/reac1_fr.me',status='unknown')
         read(124,*)
         read(124,*)(freq(j),j=1,nfreqr1)
         write(123,8010)(freq(j),j=1,nfreqr1)
         close(124)
         if(nfreqr2.ne.0) then
            open(unit=124,file='me_files/reac2_fr.me',status='unknown')
            read(124,*)
            read(124,*)(freq(j),j=1,nfreqr2)
            write(123,8010)(freq(j),j=1,nfreqr2)
            close(124)
         endif
 8010    format (1x,10G12.5)

         open(unit=124,file='data/reac1.dat',status='unknown')
         do while (WORD.NE.'NELEC')
            call LineRead (124)
         enddo
         read (124,*) nelec1
         do ielec = 1 , nelec1
            read (124,*) eelec1(ielec),gelec1(ielec)
         enddo
         close(124)
         call Lineread(0)

         open(unit=124,file='data/reac2.dat',status='unknown')
         do while (WORD.NE.'NELEC')
            call LineRead (124)
         enddo
         read (124,*) nelec2
         do ielec = 1 , nelec2
            read (124,*) eelec2(ielec),gelec2(ielec)
         enddo
         close(124)
         
         nelec12=nelec1*nelec2
         index=0
         do j=1,nelec1
            do i=1,nelec2
               index=index+1
               gelec12(index)=gelec1(j)*gelec2(i)
               eelec12(index)=eelec1(j)+eelec2(i)
            enddo
         enddo

c     write(123,*)'!************************************'
         write(123,*)' ZeroEnergy[kcal/mol]            0. '
         write(123,*)' ElectronicLevels[1/cm]        ', nelec12
         do j=1,nelec12
            write(123,*)eelec12(j),gelec12(j)
         enddo
         write(123,*)'End'
         write(123,*)'!************************************'
         close(123)

cc now save the PST block

         command1='cp -f pstr1.me ./me_blocks/pstr1.me'
         call commrun(command1)

         command1='cp -f pstr1.me rpst.me'
         call commrun(command1)
      endif
c      rewind (99)
c      write (99,450)
c 450  format ("cat ./me_files/rpst_ge1.me ./me_files/rpst_ge2.me
c     +  ./me_files/rpst_hr1.me ./me_files/rpst_hr2.me 
c     +  ./me_files/rpst_fr1.me ./me_files/rpst_fr2.me 
c     +  > rpst0.me")
c      rewind (99)
c      read (99,1120) command1
c      call commrun(command1)
c      rewind (99)
c      write (99,470) nfreqr
c 470  format(" sed -e 's/\$nfreq/",1X,I5,1X,
c     +    "/' rpst0.me > rpst.me")
c      rewind (99)
c      read (99,1120) command1
c      call commrun(command1)

cc now prepare the PST block of products      

      if(iabs.eq.1) then

cc first get the product freqs
         if ((ip1.eq.1).and.(ip2.eq.1)) then

            command1='echo 0 > temp1.log'
            call commrun(command1)
            rewind (99)
            write (99,373)
 373        format (" egrep Freq me_files/prod1_fr.me|awk '{print $2}'
     $> temp.log")
            rewind (99)
            read (99,1120) command1
            call commrun(command1)
            command1='cat temp.log temp1.log > temp2.log'
            call commrun(command1)
            open(unit=123,file='temp2.log',status='unknown')
            read(123,*)nfreqp1
            close(123)

            command1='echo 0 > temp1.log'
            call commrun(command1)
            rewind (99)
            write (99,374)
 374        format (" egrep Freq me_files/prod2_fr.me|awk '{print $2}'
     $> temp.log")
            rewind (99)
            read (99,1120) command1
            call commrun(command1)
            command1='cat temp.log temp1.log > temp2.log'
            call commrun(command1)
            open(unit=123,file='temp2.log',status='unknown')
            read(123,*)nfreqp2
            close(123)

            nfreqp=nfreqp1+nfreqp2

            filewell='wellp.me'
            call stoichiometry(filewell,stoich_well)
cc

            open(unit=123,file='pstr2.me',status='unknown')
            write (123,*) 'Barrier	B3  WP PRODS'
            write (123,*) '  RRHO '
            write (123,*) '    Stoichiometry  ',stoich_well
            write (123,*) '    Core PhaseSpaceTheory'
            open(unit=124,file='me_files/prod1_ge.me',status='unknown')
            read(124,*)
            read(124,*)
            read(124,*)
            read(124,*)cread,natom
            write (123,*) '      FragmentGeometry[angstrom]    ',natom
            do j=1,natom
               read(124,'(A80)')cread
               write(123,'(A80)')cread
            enddo
            read(124,*)
            read(124,*)cread,symr1
            close(124)
cc now write geometry of second fragment
            open(unit=124,file='me_files/prod2_ge.me',status='unknown')
            read(124,*)
            read(124,*)cread
            if(cread.eq.'Atom')then
               write (123,*) '      FragmentGeometry[angstrom]    ',1
               read(124,*)cread1,cread2
               write(123,*)cread2,' 0. 0. 0.'
               symr2=1
            else
c     read(124,*)
               read(124,*)cread,natom
               write (123,*) '      FragmentGeometry[angstrom]    ',
     $                        natom
               do j=1,natom
                  read(124,'(A80)')cread
                  write(123,'(A80)')cread
               enddo
               read(124,*)
               read(124,*)cread,symr2
            endif
            close(124)
            
            symr=symr1*symr2

            write (123,*) '		SymmetryFactor  ',symr
            write (123,*) '          PotentialPrefactor[au] 	10.' 
            write (123,*) '          PotentialPowerExponent 	6' 
            write (123,*) '	    End'

cc now write hindered rotor section
            call LineRead (0)
            open(unit=19,file='data/prod1.dat',status='unknown')
            do while (WORD.NE.'NHIND')
               call LineRead (19)
               if (WORD.EQ.'END') then
                  write (6,*) ' error in reading input in sub ktp'
                  stop
               endif
            enddo
            read (19,*) nhind
            read(19,*)
            if(nhind.ne.0)then
               do j=1,nhind
                  open(unit=20,file='me_files/prod1_hr.me',
     $                 status='unknown')
                  read(20,'(A80)')cread
                  write(123,'(A80)')cread
                  open(unit=21,file='me_files/prod1_ge.me',
     $                 status='unknown')
                  read(21,*)
                  read(21,*)
                  read(21,*)
                  read(21,*)cread,natom
                  write(123,*)'Geometry[angstrom]    ',natom
                  do i=1,natom
                     read(21,'(A80)')cread
                     write(123,'(A80)')cread
                  enddo
                  close(21)
                  read(20,'(A80)')cread
                  write(123,'(A80)')cread
                  read(20,'(A80)')cread
                  write(123,'(A80)')cread
                  read(20,'(A80)')cread
                  write(123,'(A80)')cread
                  read(20,'(A80)')cread
                  write(123,'(A80)')cread

                  read(19,*)cread,start,arr,nptot
                  read(20,*)(rotpot(ij),ij=1,nptot)
                  write(123,8010)(rotpot(ij),ij=1,nptot)
                  read(20,'(A80)')cread
                  write(123,'(A80)')cread
               enddo
            endif
            close(19)
            close(20)

            call LineRead (0)
            open(unit=19,file='data/prod2.dat',status='unknown')
            do while (WORD.NE.'NHIND')
               call LineRead (19)
               if (WORD.EQ.'END') then
                  write (6,*) ' error in reading input in sub ktp'
                  stop
               endif
            enddo
            read (19,*) nhind
            read(19,*)
            if(nhind.ne.0)then
               do j=1,nhind
                  open(unit=20,file='me_files/prod2_hr.me',
     $                  status='unknown')
                  read(20,'(A80)')cread
                  write(123,'(A80)')cread
                  open(unit=21,file='me_files/prod2_ge.me',
     $                 status='unknown')
c                  read(21,*)
                  read(21,*)
                  read(21,*)
                  read(21,*)cread,natom1
                  write(123,*)'Geometry[angstrom]    ',natom1
                  do i=1,natom1
                     read(21,'(A80)')cread
                     write(123,'(A80)')cread
                  enddo
                  close(21)
                  read(20,'(A80)')cread
                  write(123,'(A80)')cread
                  read(20,'(A80)')cread
                  write(123,'(A80)')cread
                  read(20,'(A80)')cread
                  write(123,'(A80)')cread
                  read(20,'(A80)')cread
                  write(123,'(A80)')cread
                  read(19,*)cread,start,arr,nptot
                  read(20,*)(rotpot(ij),ij=1,nptot)
                  write(123,8010)(rotpot(ij),ij=1,nptot)
                  read(20,'(A80)')cread
                  write(123,'(A80)')cread
               enddo
            endif
            close(19)
            close(20)

            write (123,*) '    Frequencies[1/cm] ',nfreqp

            open(unit=124,file='me_files/prod1_fr.me',status='unknown')
            read(124,*)
            read(124,*)(freq(j),j=1,nfreqp1)
            write(123,8011)(freq(j),j=1,nfreqp1)
            close(124)
            if(nfreqp2.ne.0) then
               open(unit=124,file='me_files/prod2_fr.me',
     $              status='unknown')
               read(124,*)
               read(124,*)(freq(j),j=1,nfreqp2)
               write(123,8011)(freq(j),j=1,nfreqp2)
               close(124)
            endif
 8011       format (1x,10G12.5)

            call LineRead (0)
            open(unit=124,file='data/prod1.dat',status='unknown')
            do while (WORD.NE.'NELEC')
               call LineRead (124)
            enddo
            read (124,*) nelec1
            do ielec = 1 , nelec1
               read (124,*) eelec1(ielec),gelec1(ielec)
            enddo
            close(124)
            call Lineread(0)

            open(unit=124,file='data/prod2.dat',status='unknown')
            do while (WORD.NE.'NELEC')
               call LineRead (124)
            enddo
            read (124,*) nelec2
            do ielec = 1 , nelec2
               read (124,*) eelec2(ielec),gelec2(ielec)
            enddo
            close(124)
            
            nelec12=nelec1*nelec2
            index=0
            do j=1,nelec1
               do i=1,nelec2
                  index=index+1
                  gelec12(index)=gelec1(j)*gelec2(i)
                  eelec12(index)=eelec1(j)+eelec2(i)
               enddo
            enddo

c      write(123,*)'!************************************'
            write(123,*)' ZeroEnergy[kcal/mol]          ', prod_en
            write(123,*)' ElectronicLevels[1/cm]        ', nelec12
            do j=1,nelec12
               write(123,*)eelec12(j),gelec12(j)
            enddo
            write(123,*)'End'
            write(123,*)'!************************************'
            close(123)

cc now save the PST block

            command1='cp -f pstr2.me ./me_blocks/pstr2.me'
            call commrun(command1)

            command1='cp -f pstr2.me ppst.me'
            call commrun(command1)
         endif
      endif
c      stop


c      rewind (99)
c      write (99,550)
c 550  format ("cat ./me_files/ppst_ge1.me ./me_files/ppst_ge2.me
c     +  ./me_files/ppst_hr1.me ./me_files/ppst_hr2.me 
c     +  ./me_files/ppst_fr1.me ./me_files/ppst_fr2.me 
c     +  > ppst0.me")
c      rewind (99)
c      read (99,1120) command1
c      call commrun(command1)
c
c      rewind (99)
c      write (99,560) prod_en
c 560  format(" sed -e 's/\$proden/",1X,F7.2,1X,
c     +    "/' ppst0.me > ppst1.me")
c      rewind (99)
c      read (99,1120) command1
c      call commrun(command1)


ccccccccccccccccccccccc
cc old determination of freqs
c      rewind (99)
c      write (99,570) nfreqp
c 570  format(" sed -e 's/\$nfreq/",1X,I5,1X,
c     +    "/' ppst1.me > ppst.me")
c      rewind (99)
c      read (99,1120) command1
c      call commrun(command1)


cccccccccccccccccccccc

      if (iadd.eq.1.or.iiso.eq.1) then
         if(nts.eq.1) then
            rewind (99)
            write (99,1025) 
            rewind (99)
            read (99,1120) command1
            call commrun(command1)
1025  format(" sed -ie 's/Barrier TS REACS WR/ Barrier TS REACS WP/g'
     +   temp3.me")
            rewind (99)
            write (99,1925) 
            rewind (99)
            read (99,1120) command1
            call commrun(command1)
1925    format(" sed -ie 's/Barrier B2 WR WP/Barrier TS REACS WP/g'
     +    temp3.me")
            if(ip1.eq.1.and.ip2.eq.1)then
               rewind (99)
               write (99,1125) 
               rewind (99)
               read (99,1120) command1
               call commrun(command1) 
1125  format(" sed -ie 's/Barrier TS REACS WP/Barrier TS REACS PRODS/g'
     +    temp3.me")
            endif
         else if (nts.eq.2) then
            rewind (99)
            if(ip1.eq.1.and.ip2.eq.1)then
               write (99,1225) 
               rewind (99)
               read (99,1120) command1
               call commrun(command1)
               rewind(99)
               write (99,1325) 
               rewind (99)
               read (99,1120) command1
               call commrun(command1)
1225   format(" sed -ie 's/Barrier TS REACS WR/Barrier B2 WR PRODS/g'
     +   temp3.me")
1325   format(" sed -ie 's/Barrier B2 WR WP/Barrier B2 WR PRODS/g'
     +   temp3.me")
            else
               rewind (99)
               write (99,1026) 
               rewind (99)
               read (99,1120) command1
               call commrun(command1)
            endif
         endif
      endif

      if (iabs.eq.1) then
         if(nts.eq.1) then
            rewind (99)
            write (99,1025) 
            rewind (99)
            read (99,1120) command1
            call commrun(command1)
            rewind (99)
            write (99,3125) 
            rewind (99)
            read (99,1120) command1
            call commrun(command1)
            rewind (99)
            write (99,3126) 
            rewind (99)
            read (99,1120) command1
            call commrun(command1)
3125   format(" sed -ie 's/Barrier TS REACS WP/Barrier TS REACS WR/g'
     +   temp3.me")
3126   format(" sed -ie 's/Barrier B2 WR WP/Barrier TS REACS WR/g'
     +   temp3.me")
         else if (nts.eq.2) then
            rewind (99)
            if(ip1.eq.1.and.ip2.eq.1)then
               write (99,1225) 
               rewind (99)
               read (99,1120) command1
               call commrun(command1)
               rewind(99)
               write (99,1325) 
               rewind (99)
               read (99,1120) command1
               call commrun(command1)
            else if (ipr1.eq.1)then
               rewind (99)
               write (99,1027) 
               rewind (99)
               read (99,1120) command1
               call commrun(command1)
            else
               rewind (99)
               write (99,1026) 
               rewind (99)
               read (99,1120) command1
               call commrun(command1)
            endif
         else if (nts.eq.3) then
            rewind (99)
            write (99,1026) 
            rewind (99)
            read (99,1120) command1
            call commrun(command1)
1026  format(" sed -ie 's/Barrier TS REACS WP/Barrier B2 WR WP/g'
     +    temp3.me")
1027  format(" sed -ie 's/Barrier TS REACS WR/Barrier B2 WR PRODS/g'
     +   temp3.me")
         endif
      endif

c finally ready to write full me input and run mess

      rewind (99)
      if ((nts.eq.1).and.(iabs.eq.1))then
         write (99,1011)
      endif
      if ((nts.eq.1).and.(iadd.eq.1)) write (99,1015)
      if ((nts.eq.1).and.(iiso.eq.1)) then
         if(ipw.eq.1) then
            write (99,1015)
         else if (ip1.eq.1.and.ip2.eq.1)then
            write (99,1017)
         else if (ipr1.eq.1)then
            write (99,1017)
         else
            write (99,1015)
         endif
      endif
      if ((nts.eq.1).and.(ibeta.eq.1)) write (99,1016)
      if ((iabs.eq.1).and.(nts.eq.2)) write (99,1012)
      if ((iadd.eq.1).and.(nts.eq.2)) then
         if (ipr1.eq.1)then
            write (99,1018)
         else
            write (99,1013)
         endif
      endif
      if (nts.eq.3) write (99,1014)
      rewind (99)
      read (99,1120) command1

 1011 format ("cat ./data/me_head.dat ./reactants.me 
     +  ./wellr.me 
     +  ./temp3.me > me_ktp.inp")

 1015 format ("cat ./data/me_head.dat ./reactants.me 
     +  ./wellp.me 
     +  ./temp3.me > me_ktp.inp")
 1017 format ("cat ./data/me_head.dat ./reactants.me 
     +  ./products.me
     +  ./temp3.me > me_ktp.inp")

 1016 format ("cat ./data/me_head.dat ./reactants.me 
     +  ./products.me 
     +  ./temp3.me > me_ktp.inp")

 1012 format ("cat ./data/me_head.dat ./reactants.me 
     +  ./products.me ./wellr.me ./rpst.me 
     +  ./temp3.me > me_ktp.inp")

 1013 format ("cat ./data/me_head.dat ./reactants.me 
     +  ./wells.me ./rpst.me 
     +  ./temp3.me > me_ktp.inp")
 1018 format ("cat ./data/me_head.dat ./reactants.me 
     +  ./wells.me ./products.me ./rpst.me 
     +  ./temp3.me > me_ktp.inp")

 1014 format ("cat ./data/me_head.dat ./reactants.me 
     +  ./products.me ./wells.me ./rpst.me ./ppst.me
     +  ./temp3.me > me_ktp.inp")

      call commrun(command1)

c      stop
      
      command1='mess me_ktp.inp'

      call commrun(command1)

c      rewind (99)
c      write (99,1004)
c      rewind (99)
c      read (99,1120) command1

c 1004 format (" egrep -A25 'K\)    REACS->WR'
c     +     ./rate.out > ./output/kextract.out")
c      call commrun(command1)

      command1='cp -f me_ktp.inp output'
      call commrun(command1)
      command1='cp -f rate.out output'
      call commrun(command1)

      close (unit=99,status='keep')
c      close (unit=105,status='keep')
c      close (unit=106,status='keep')

 1120 FORMAT (A300)
 1020 continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine modarr

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      open (unit=125,file='modarr.dat',status='old')
      open (unit=126,file='modarr.out',status='unknown')
      close (unit=125,status='keep')
      close (unit=126,status='keep')
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ts_zmat(atomlabel,
     $ xinti,intcoor,natom1,natom2,natomt1,natomt2,ntau1,ntau2,bislab)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension tau(ntaumx),taumn(ntaumx),taumx(ntaumx),freq(nmdmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim) 

      dimension ibconn(natommx),iaconn(natommx),idconn(natommx)
      dimension idummy(natommx)
      character*20 bislab(ntaumx)
      character*20 bname(natommx),anname(natommx),dname(natommx)
     $ ,atname(natommx),bconnt(natommx),aconnt(natommx),dconnt(natommx)


      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      character*60 atomlabel(natommx)
      character*30 intcoor(3*natommx)
      character*30 gmem

      include 'filcomm.f'

      open (unit=116,file='./output/ts_zmat.out',status='unknown')
      write (116,*) 'entering ts_zmat'
      open (unit=21,file='./data/reac1.dat',status='old')

c read input from reac_1 for first part of ts
c first set of lines are just garbage from before
 1000 continue
      CALL LineRead (21)
      if (WORD.eq.'CHARGE') go to 1500
      write (116,*) 'WORD',WORD
      go to 1000
 1500 continue
      write (116,*) 'entering atom reads',natomt1
      read (21,*)
      do iatom = 1 , natomt1
         read (21,'(A60)') atomlabel(iatom)
      enddo
      rewind(21)

      natom1p = 3*natom1-6
      if(natom1.eq.2) natom1p=1

cc now read atom labels for react1
      if (natom1.ne.1) then
         do while (WORD.NE.'INTCOOR')
            call LineRead (21)
         enddo
         ncoord = 3*natom1-6-ntau1
         if (natom1.eq.1) ncoord = 0
         if (natom1.eq.2) ncoord = 1
         do icoord = 1 , ncoord
            call LineRead (21)
            intcoor(icoord) = word
         enddo
         rewind(21)
      endif

      do while (WORD.NE.'NTAU')
         call LineRead (21)
      enddo
      read (21,*) ntautemp
c      write (*,*) ntautemp
      
      if (ntautemp.ne.0) then
         read (21,*)
         do itau = 1 , ntautemp
            call LineRead (21)
cc modified intcoor vector index
            intcoor(natom1p+itau-ntau1) = word
c            read (21,*) intcoor(natom1p+itau-ntau1),taumn,taumx
c            write (*,*) intcoor(natom1p+itau-ntau1),taumn,taumx
         enddo
      endif
      close (unit=21,status='keep')

cc now read optimized coordinates

      open (unit=23,file='./output/reac1_opt.out',status='old')

      read (23,*)
      
      do icoord = 1 , natom1p
         read (23,*) xinti(icoord)
      enddo

      close (unit=23,status='keep')
  
c      write(*,*)'natomp1 ',natom1p
c      write(*,*)'intcoord ',intcoor(1)
c      stop

c    do icoord = 1 , natom1p-ntau1
c         call LineRead (21)
c         intcoor(icoord) = word
c        write (116,*) 'xint test',xinti(icoord),icoord
c      enddo


c      read (21,*)
c      do itau = 1 , ntau1
c         call LineRead (21)
cc modified intcoor vector index
cc         intcoor(natom1p+itau) = word
c         intcoor(natom1p+itau-ntau1) = word
c         read (23,*) xinti(natom1p-ntau1+itau)
c      enddo

c add lines connecting reac_1 and reac_2 via abs_site plus dummy atom
c to do so must have atom labels and some connections from file

c start with dummy atom for abstraction case
      if (iabs.eq.1) then
         if(natom1.ne.2) then
            OPEN (unit=99,status='unknown')
            REWIND (99)
            write (99,1901) 'X111',isite,'2.',jsite,'90.',ksite,
     $    '180.'
 1901       format (a4,1x,i4,1x,a2,1x,i4,1x,a3,1x,i4,1x,a4)
            REWIND (99)
            read (99,'(A60)') atomlabel(iatom)
            close (unit=99,status='keep')
         else if(natom1.eq.2.and.natomt1.eq.2)then
            OPEN (unit=99,status='unknown')
            REWIND (99)
            write (99,1902) 'X111',isite,'2.',jsite,'90.'
 1902       format (a4,1x,i4,1x,a2,1x,i4,1x,a3,1x)
            REWIND (99)
            read (99,'(A60)') atomlabel(iatom)
            close (unit=99,status='keep')
         else if(natom1.eq.2.and.natomt1.eq.3)then
            OPEN (unit=99,status='unknown')
            REWIND (99)
            write (99,1903) 'X111',isite,'2.',jsite,'90.',ksite,
     $    '180.'
 1903       format (a4,1x,i4,1x,a2,1x,i4,1x,a3,1x,i4,1x,a4)
            REWIND (99)
            read (99,'(A60)') atomlabel(iatom)
            close (unit=99,status='keep')
         endif
      endif

      open (unit=22,file='./data/reac2.dat',status='unknown')
 2000 continue
      CALL LineRead (22)
      if (WORD.eq.'CHARGE') go to 2500
      write (116,*) 'WORD reac_2a',WORD
      go to 2000
 2500 continue

c read remaining input from reac_2 for second part of ts
      read (22,*)
      if (iadd.eq.1) then
         iatomi = natomt1+1
         iatomf = natomt1+natomt2
      endif
      if (iabs.eq.1) then
         iatomi = natomt1+2
         iatomf = natomt1+natomt2+1
c        iatomi = iatomi+1
c        iatomf = iatomf+1
      endif
cc
cc
      do iatom = iatomi , iatomf
c     do iatom = natomt1+2 , natomt1+natomt2+1
         read (22,'(A60)') atomlabel(iatom)
         write (116,*) 'atomlabel reac_2b',atomlabel(iatom),iatom
         if (iatom.eq.iatomi) then
c        if (iatom.eq.natomt1+2) then
            OPEN (unit=99,status='unknown')
            REWIND (99)
            write (99,*) atomlabel(iatom)
            REWIND (99)
            call LineRead (99) 
            REWIND (99)
            if (iabs.eq.1) then
               if(natom1.eq.2.and.natomt1.eq.2) then
                  write (99,2601) word,isite,'rts',natomt1+1,
     $               'aabs1',jsite,'180.'
               else if(natom1.eq.2.and.natomt1.eq.3) then
                  write (99,2601) word,isite,'rts',natomt1+1,
     $               'aabs1',jsite,'180.'
               else
                  write (99,2601) word,isite,'rts',natomt1+1,
     $               'aabs1',jsite,'babs1'
               endif
            endif
            if (iadd.eq.1) then 
               if(natom1.eq.2.and.natomt1.eq.2) then
                  write (99,2601) word,isite,'rts',jsite,
     $       'aabs1',ksite,'180.'
               else if (natom1.eq.2.and.natomt1.eq.3) then
                  write (99,2601) word,isite,'rts',jsite,
     $       'aabs1',ksite,'180.'
               else
                  write (99,2601) word,isite,'rts',jsite,
     $       'aabs1',ksite,'babs1'
               endif
            endif
 2601       format (1x,a8,1x,i4,1x,a3,1x,i4,1x,a5,1x,i4,1x,a5)
            REWIND (99)
            read (99,'(A60)') atomlabel(iatom)
            write (116,*) 'atomlabel reac_3a',atomlabel(iatom)
            close (unit=99,status='keep')
         endif
         if (iatom.eq.iatomi+1) then
c        if (iatom.eq.natomt1+3) then
            OPEN (unit=99,status='unknown')
            REWIND (99)
            write (99,*) atomlabel(iatom)
            REWIND (99)
            call LineRead (99) 
            REWIND (99)
            if (iabs.eq.1) then
               if(natom2.eq.2.and.natomt2.eq.3)then
                  write (99,2602) word,word2,word3,isite,
     $       ' 90. ',natomt1+1,'babs2'
               else
                  write (99,2602) word,word2,word3,isite,
     $       'aabs2',natomt1+1,'babs2'
               endif
            endif
            if (iadd.eq.1) then
               if(natom2.eq.2.and.natomt2.eq.3)then
c                  write(*,*)'passing here'
c                  write(*,*)'stop'
                  write (99,2602) word,word2,word3,isite,
     $       ' 90. ',jsite,'babs2'
               else
                  write (99,2602) word,word2,word3,isite,
     $       'aabs2',jsite,'babs2'
               endif
            endif
c           write (99,2602) word,word2,word3,isite,'aabs2',
c    $       natomt1+1,'babs2'
 2602       format (1x,3a8,1x,i4,1x,a5,1x,i4,1x,a5)
            REWIND (99)
            read (99,'(A60)') atomlabel(iatom)
            write (116,*) 'atomlabel reac_3b',atomlabel(iatom)
            close (unit=99,status='keep')
         endif
         if (iatom.eq.iatomi+2) then
c        if (iatom.eq.natomt1+4) then
            OPEN (unit=99,status='unknown')
            REWIND (99)
            write (99,*) atomlabel(iatom)
            REWIND (99)
            call LineRead (99)
            REWIND (99)
c           if (iabs.eq.1) write (99,2603) word,word2,word3,word4,word5,
c    $       isite,'babs3'
c           if (iadd.eq.1) write (99,2603) word,word2,word3,word4,word5,
c    $       isite,'babs3'
            write (99,2603) word,word2,word3,word4,word5,isite,
     $       'babs3'
 2603       format (1x,5a8,1x,i4,1x,a5)
            REWIND (99)
            read (99,'(A60)') atomlabel(iatom)
            write (116,*) 'atomlabel reac_3c',atomlabel(iatom)
            close (unit=99,status='keep')
         endif
      enddo

      natomtp = 3*(natom1+natom2)-12
      if (natom2.eq.1) natomtp = natomtp+3
      if (natom2.eq.2) natomtp = natomtp+1
      if (natom1.eq.2) natomtp = natomtp+1
c      write(*,*)'natomtp is',natomtp
c      stop
cc now read remaining atom labels for react2
      if (natom2.ne.1) then
         do while (WORD.NE.'INTCOOR')
            call LineRead (22)
         enddo
         ncoord = 3*natom2-6-ntau1
         if (natom2.eq.1) ncoord = 0
         if (natom2.eq.2) ncoord = 1
         nstart=3*natom1-6+1
         if(natom1.eq.2)nstart=nstart+1
         do icoord=nstart , natomtp-ntau2
            call LineRead (22)
            intcoor(icoord) = word
         enddo
         rewind(22)
      endif

      rewind(22)
      do while (WORD.NE.'NTAU')
         call LineRead (22)
      enddo
      read (22,*) ntautemp

      if (ntautemp.ne.0) then
         read (22,*)
         do itau = 1 , ntautemp
            call LineRead (22)
cc modified intcoor vector index
            intcoor(natomtp+itau-ntau2) = word
c            read (21,*) intcoor(natom1p+itau-ntau1),taumn,taumx
c            write (*,*) intcoor(natom1p+itau-ntau1),taumn,taumx
         enddo
      endif

      open (unit=24,file='./output/reac2_opt.out',status='unknown')
      read (24,*)
      nstart=3*natom1-6+1
      if(natom1.eq.2)nstart=nstart+1
      do icoord = nstart , natomtp
         read (24,*) xinti(icoord)
      enddo

      if (natom2.eq.1) then
         if(natom1.eq.2.and.natomt1.eq.2)then
            xinti(natomtp+1) = aabs1
            xinti(natomtp+2) = rts
            intcoor(natomtp+1) = 'AABS1'
            intcoor(natomtp+2) = 'RTS'
         else if(natom1.eq.2.and.natomt1.eq.3)then
            xinti(natomtp+1) = aabs1
            xinti(natomtp+2) = rts
            intcoor(natomtp+1) = 'AABS1'
            intcoor(natomtp+2) = 'RTS'
         else
            xinti(natomtp+1) = aabs1
            xinti(natomtp+2) = babs1
            xinti(natomtp+3) = rts
            intcoor(natomtp+1) = 'AABS1'
            intcoor(natomtp+2) = 'BABS1'
            intcoor(natomtp+3) = 'RTS'
         endif
      endif
c      write(*,*)'test'
c      write(*,*) 'intcoor ', intcoor(natomtp+2)
c      STOP
      if (natom2.eq.2) then
         if(natom1.eq.2.and.natomt1.eq.2)then
            xinti(natomtp+1) = aabs1
            xinti(natomtp+2) = aabs2
            xinti(natomtp+3) = babs2
            xinti(natomtp+4) = rts
            intcoor(natomtp+1) = 'AABS1'
            intcoor(natomtp+2) = 'AABS2'
            intcoor(natomtp+3) = 'BABS2'
            intcoor(natomtp+4) = 'RTS'
         else if(natom1.eq.2.and.natomt1.eq.3)then
            xinti(natomtp+1) = aabs1
            xinti(natomtp+2) = aabs2
            xinti(natomtp+3) = babs2
            xinti(natomtp+4) = rts
            intcoor(natomtp+1) = 'AABS1'
            intcoor(natomtp+2) = 'AABS2'
            intcoor(natomtp+3) = 'BABS2'
            intcoor(natomtp+4) = 'RTS'
         else if(natomt2.eq.3)then
            xinti(natomtp+1) = aabs1
            xinti(natomtp+2) = babs1
            xinti(natomtp+3) = babs2
            xinti(natomtp+4) = babs3
            xinti(natomtp+5) = rts
            intcoor(natomtp+1) = 'AABS1'
            intcoor(natomtp+2) = 'BABS1'
            intcoor(natomtp+3) = 'BABS2'
            intcoor(natomtp+4) = 'BABS3'
            intcoor(natomtp+5) = 'RTS'
         else
            xinti(natomtp+1) = aabs1
            xinti(natomtp+2) = babs1
            xinti(natomtp+3) = aabs2
            xinti(natomtp+4) = babs2
            xinti(natomtp+5) = rts
            intcoor(natomtp+1) = 'AABS1'
            intcoor(natomtp+2) = 'BABS1'
            intcoor(natomtp+3) = 'AABS2'
            intcoor(natomtp+4) = 'BABS2'
            intcoor(natomtp+5) = 'RTS'
         endif
      endif
      if (natom2.ge.3) then
         if(natom1.eq.2.and.natomt1.eq.2)then
            xinti(natomtp+1) = aabs1
            xinti(natomtp+2) = aabs2
            xinti(natomtp+3) = babs2
            xinti(natomtp+4) = babs3
            xinti(natomtp+5) = rts
            intcoor(natomtp+1) = 'AABS1'
            intcoor(natomtp+2) = 'AABS2'
            intcoor(natomtp+3) = 'BABS2'
            intcoor(natomtp+4) = 'BABS3'
            intcoor(natomtp+5) = 'RTS'
         else
            xinti(natomtp+1) = aabs1
            xinti(natomtp+2) = babs1
            xinti(natomtp+3) = aabs2
            xinti(natomtp+4) = babs2
            xinti(natomtp+5) = babs3
            xinti(natomtp+6) = rts
            intcoor(natomtp+1) = 'AABS1'
            intcoor(natomtp+2) = 'BABS1'
            intcoor(natomtp+3) = 'AABS2'
            intcoor(natomtp+4) = 'BABS2'
            intcoor(natomtp+5) = 'BABS3'
            intcoor(natomtp+6) = 'RTS'
         endif
      endif
         
      close (unit=22,status='keep')
      close (unit=24,status='keep')
      close (unit=116,status='keep')

cc now  check special cases
cc check if dummy atom used as second atom in reac2 in species with natom>2
cc in which case the atomlabel should be changed for 

      natomtot=natom1+natom2
      call read_zmat(atomlabel,natomtot,iatomf,intcoor,bislab,ibconn,
     $ iaconn,idconn,bname,anname,dname,atname,idummy,isited,jsited,
     $ ksited,aconnt,bconnt,dconnt)

c      do j=1,iatomf
c         write(*,*)'idummy is ', j, idummy(j)
c      enddo
c      write(*,*)'idummy iatomi+1is ', idummy(iatomi+1)
 
cc search for coordinate assigned in remaining dof
      ntotcoo=3*natomtot-6

c      write(*,*)'iatomi is ',iatomi
c      stop
      if(idummy(iatomi+1).eq.1.and.natom2.gt.2)then
         do j=iatomi+2,iatomf
            if(idummy(j).ne.1)then
               icheck=0
               do ik=1,ntotcoo
                  if(intcoor(ik).eq.dname(j))icheck=1
               enddo
               if(icheck.eq.0)then
                  isub=j
                  write(*,*)'found unassigned coord',j,dname(j)
               endif
            endif
         enddo
cc now substitute variable dname(j) with aabs2 and change atom labels
         open (unit=99,status='unknown')
         write(99,*)dname(isub)
         rewind(99)
         read(99,*)xinti(natomtp+3)
         dname(isub)='AABS2'
c         write(*,*)'new param is',xinti(natomtp+3)
cc
         rewind(99)
         write (99,*) atomlabel(iatomi+1)
         rewind (99)
         call LineRead (99)
         rewind(99)
         if(iabs.eq.1)then
            write (99,2602) word,word2,word3,isite,
     $       ' 90. ',natomt1+1,'babs2'
         else if(iadd.eq.1)then
            write (99,2602) word,word2,word3,isite,
     $       ' 90. ',natomt1,'babs2'
         endif
         rewind(99)
         read (99,'(A60)') atomlabel(iatomi+1)
c         write(*,*)'new atom label 1', atomlabel(iatomi+1)
c         close(99)
c         open (unit=99,status='unknown')
         rewind(99)
         write (99,*) atomlabel(isub)
         REWIND(99)
         call LineRead3 (99)
         rewind(99)
         write (99,2604) word,word2,word3,word4,
     $       word5,word6,'aabs2'
         rewind(99)
         read (99,'(A60)') atomlabel(isub)         
      endif
 2604 format (1x,3a7,1x,a7,1x,a7,1x,a7,1x,a7)

c      do j=1,iatomf
c         write(*,*)atomlabel(j)
c      enddo
c      do j=1,ntotcoo
c         write(*,*)intcoor(j),xinti(j)
c      enddo

c      stop
      write (116,*) 'past z-matrix'
      return
      end

c (i) employs method_reac_opt
c (ii) employs method_grid_opt
c (iii) employs method_ts_0
c (iv) employs method_ts_1
c (v) employs method_ts_tauo
c (vi) employs method_ts_hl
c (vii) employs method_1dtau
c (viii) employs method_irc
      


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ts_zmat_iso(atomlabel,atomlabel1,
     $ xinti,intcoor,natom1,natomt1,ntau1,iopt,dstart)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension tau(ntaumx),taumn(ntaumx),taumx(ntaumx),freq(nmdmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim),xintt(3*natommx)

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      character*60 atomlabel(natommx)
      character*60 atomlabel1(natommx)
      character*30 intcoor(3*natommx)
      character*30 intcoort(3*natommx)
      character*30 gmem
      character*40 filename
      character*30 cjunk,adislab,aanglab,adihlab
      character*30 cname1,cname2,cname3,cname4
      character*1 xread

      include 'filcomm.f'

      open (unit=116,file='./output/ts_zmat.out',status='unknown')
      write (116,*) 'entering ts_zmat_iso'
      open (unit=21,file='./data/reac1.dat',status='old')

c read input from reac_1 for first part of ts
c first set of lines are just garbage from before
 1000 continue
      CALL LineRead (21)
      if (WORD.eq.'CHARGE') go to 1500
c      write (116,*) 'WORD',WORD
      go to 1000
 1500 continue
      write (116,*) 'entering atom reads',natomt1
      read (21,*)
      do iatom = 1 , natomt1
         read (21,'(A60)') atomlabel(iatom)
      enddo
      rewind(21)

      natom1p = 3*natom1-6

cc now check for dummy atoms and update ireact ixyz ij ik indexes

 1001 continue
      CALL LineRead (21)
      if (WORD.eq.'CHARGE') go to 1501
c      write (116,*) 'WORD',WORD
      go to 1001
 1501 continue
      write (116,*) 'entering atom reads',natomt1
      read (21,*)
      ndummy=0
      do iatom = 1 , natomt1
         read (21,*) xread
         if(xread.eq.'X')then
            ndummy=ndummy+1
         endif
         if(iatom.eq.ireact)ireact_xyz=ireact-ndummy
         if(iatom.eq.isite)isite_xyz=isite-ndummy
         if(iatom.eq.jsite)ji_xyz=jsite-ndummy
         if(iatom.eq.ksite)ki_xyz=ksite-ndummy
      enddo
      rewind(21)
c      write(*,*)'ireact is ', ireact
c      write(*,*)'ireact_xyz is ', ireact_xyz
c      write(*,*)'isite_xyz is ', isite_xyz
c      write(*,*)'ji_xyz is ', ji_xyz
c      write(*,*)'ki_xyz is ', ki_xyz
c      stop

cc  read coordinate labels for react1
      if (natom1.ne.1) then
         do while (WORD.NE.'INTCOOR')
            call LineRead (21)
         enddo
         ncoord = 3*natom1-6-ntau1
         if (natom1.eq.1) ncoord = 0
         if (natom1.eq.2) ncoord = 1
         do icoord = 1 , ncoord
            call LineRead (21)
            intcoor(icoord) = word
         enddo
         rewind(21)
      endif

      do while (WORD.NE.'NTAU')
         call LineRead (21)
      enddo
      read (21,*) ntautemp
c      write (*,*) ntautemp
      
      if (ntautemp.ne.0) then
         read (21,*)
         do itau = 1 , ntautemp
            call LineRead (21)
cc modified intcoor vector index
            intcoor(natom1p+itau-ntau1) = word
c            read (21,*) intcoor(natom1p+itau-ntau1),taumn,taumx
c            write (*,*) intcoor(natom1p+itau-ntau1),taumn,taumx
         enddo
      endif
      rewind(21)

cc read coordinate labels for reacting atom


      do while (WORD.NE.'CHARGE')
         call LineRead (21)
      enddo
      read (21,*)
      do iatom = 1 , natomt1
         if(iatom.eq.ireact)then
            read(21,*)cname1,cjunk,adislab,cjunk,aanglab,cjunk,adihlab
         else
            call LineRead (21)
         endif
      enddo
      close (unit=21,status='keep')

cc convert iatom parameters to upper case

      open (unit=99,status='unknown')
      write(99,*)adislab
      rewind (99)
      call LineRead (99)
      adislab=word
      close (99)

      open (unit=99,status='unknown')
      write(99,*)aanglab
      rewind (99)
      call LineRead (99)
      aanglab=word
      close (99)

      open (unit=99,status='unknown')
      write(99,*)adihlab
      rewind (99)
      call LineRead (99)
      adihlab=word
      close (99)

cc now read optimized coordinates

      if( iopt.ne.0) then
         open (unit=99,status='unknown')
         write(99,1200)iopt
         rewind (99)
         read (99,'(A40)') filename
         close (99)
      else
         filename='./output/reac1_opt.out'
      endif

 1200 format ("./output/reac1_opt_"I0.2".out")

      open(unit=23,file=filename,status='unknown')

      read (23,*)cjunk,cjunk,iminst
      
      do icoord = 1 , natom1p
         read (23,*) xinti(icoord)
      enddo

      close (unit=23,status='keep')
  
      write(*,*)'At params ', adislab,aanglab,adihlab

cc now determine new parameters for reacting atoms

      if( iopt.ne.0) then
         open (unit=99,status='unknown')
         write(99,1201)iopt
         rewind (99)
         read (99,'(A40)') filename
         close (99)
      else
         if(iminst.ne.0)then
            open (unit=99,status='unknown')
            write(99,1201)iminst
            rewind (99)
            read (99,'(A40)') filename
            close (99)
         else
            filename='reac1_l1.xyz'
         endif
      endif

 1201 format ("./geoms/reac1_"I0.2".xyz")

      open(unit=23,file=filename,status='unknown')
      read (23,*)
      read (23,*)
      do j=1,natom1
         read (23,*)cjunk,coox,cooy,cooz
         if(j.eq.ireact_xyz)then
            coox1=coox
            cooy1=cooy
            cooz1=cooz
         endif
         if(j.eq.isite_xyz)then
            coox2=coox
            cooy2=cooy
            cooz2=cooz
         endif
         if(j.eq.ji_xyz)then
            coox3=coox
            cooy3=cooy
            cooz3=cooz
         endif
         if(j.eq.ki_xyz)then
            coox4=coox
            cooy4=cooy
            cooz4=cooz
         endif
      enddo
      close(23)

      write(*,*)'coox1 ',coox1
      write(*,*)'coox2 ',coox2
      write(*,*)'coox3 ',coox3
      write(*,*)'coox4 ',coox4

      readis=sqrt((coox1-coox2)**2+(cooy1-cooy2)**2+(cooz1-cooz2)**2)
      dstart=readis

c      open(unit=23,file='dihed.dat',status='unknown')
c      write(23,*)coox1,cooy1,cooz1
c      write(23,*)coox2,cooy2,cooz2
c      write(23,*)coox3,cooy3,cooz3
c      write(23,*)coox4,cooy4,cooz4
c      close(23)
c      call dihedral
c      open(unit=23,file='dihed.res',status='unknown')
c      read(23,*)reaang
c      read(23,*)readih
c      close(23)

cc to be checked - changed on 19/10

      call xyz_to_zmat(coox1,cooy1,cooz1,coox2,cooy2,cooz2,
     +  coox3,cooy3,cooz3,coox4,cooy4,cooz4,reaang,readih)

      write(*,*)'readis is ',readis

cc now create new input

cc first update z-matrix, by connecting the reacting atom and isite

      do iatom=1,natomt1
         open(unit=99,status='unknown')
         write(99,*)atomlabel(iatom)
         rewind(99)
         call LineRead(99)
         if(iatom.eq.ireact)then
            cname1=word
         endif
         if(iatom.eq.isite)then
            cname2=word
         endif
         if(iatom.eq.jsite)then
            cname3=word
         endif
         if(iatom.eq.ksite)then
            cname4=word
         endif
         close(99)
      enddo

      write(*,*)'cname 1 is ',cname1
      write(*,*)'cname 2 is ',cname2
      write(*,*)'cname 3 is ',cname3
      write(*,*)'cname 4 is ',cname4
c      write(*,*)'ireact is',ireact
c      stop
cc now update atomlabels

      do iatom=1,natomt1
         if(iatom.eq.ireact.and.iatom.eq.3)then
            open(unit=99,status='unknown')
            write(99,1202)cname1,cname2,cname3,aanglab
            rewind(99)
            read(99,'(A60)')atomlabel1(iatom)
            close(99)
         else if(iatom.eq.ireact)then
            open(unit=99,status='unknown')
            write(99,1203)cname1,cname2,cname3,aanglab,cname4,
     $                    adihlab
            rewind(99)
            read(99,'(A60)')atomlabel1(iatom)
            close(99)
         else
            atomlabel1(iatom)=atomlabel(iatom)
         endif
      enddo

 1202 format (A3,1X,A3,1X,'RTS',1X,A3,1X,A5,1X)
 1203 format (A3,1X,A3,1X,'RTS',1X,A3,1X,A5,1X,A3,1X,A5,1X)

      do iatom=1,natomt1
         write(*,*)atomlabel1(iatom)
      enddo

cc now update xint vector

      do icoord = 1 , natom1p
         if(intcoor(icoord).eq.adislab)then
            xinti(icoord)=readis
         endif
         if(intcoor(icoord).eq.aanglab)then
            xinti(icoord)=reaang
         endif
         if(intcoor(icoord).eq.adihlab)then
            xinti(icoord)=readih
         endif
      enddo

c      do icoord = 1 , natom1p
c         write(*,*)intcoor(icoord),xinti(icoord)
c      enddo

cc      now move the bond coordinate to the last place and call it RTS

      index=0
      do icoord = 1 , natom1p
         if(intcoor(icoord).ne.adislab) then
            index=index+1
            intcoort(index)=intcoor(icoord)
            xintt(index)=xinti(icoord)
         endif
      enddo
      xintt(natom1p)=readis
      intcoort(natom1p)='RTS'

cc update coordinate names and print

      do iatom=1,natomt1
         write(*,*)atomlabel1(iatom)
      enddo
      do icoord = 1 , natom1p
            intcoor(icoord)=intcoort(icoord)
            xinti(icoord)=xintt(icoord)
            write(*,*)intcoor(icoord),xinti(icoord)
      enddo
         
      write (116,*) 'past z-matrix'

      close (unit=22,status='keep')
      close (unit=24,status='keep')
      close (unit=116,status='keep')

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ts_zmat_beta(atomlabel,
     $ xinti,intcoor,natom1,natomt1,ntau1,iopt,dstart)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension tau(ntaumx),taumn(ntaumx),taumx(ntaumx),freq(nmdmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim),xintt(3*natommx)

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      character*60 atomlabel(natommx)
      character*60 atomlabel1(natommx)
      character*30 intcoor(3*natommx)
      character*30 intcoort(3*natommx)
      character*30 gmem
      character*40 filename
      character*30 cjunk,adislab,aanglab,adihlab
      character*30 cname1,cname2,cname3,cname4

      include 'filcomm.f'

      open (unit=116,file='./output/ts_zmat.out',status='unknown')
      write (116,*) 'entering ts_zmat_beta'
      open (unit=21,file='./data/reac1.dat',status='old')

c read input from reac_1 for first part of ts
c first set of lines are just garbage from before

      write(6,*)'entering ts zmat beta'

 1000 continue
      CALL LineRead (21)
      if (WORD.eq.'CHARGE') go to 1500
      write (116,*) 'WORD',WORD
      go to 1000
 1500 continue
      write (116,*) 'entering atom reads',natomt1
      read (21,*)
      do iatom = 1 , natomt1
         read (21,'(A60)') atomlabel(iatom)
      enddo
      rewind(21)

      natom1p = 3*natom1-6

cc  read coordinate labels for react1
      if (natom1.ne.1) then
         do while (WORD.NE.'INTCOOR')
            call LineRead (21)
         enddo
         ncoord = 3*natom1-6-ntau1
         if (natom1.eq.1) ncoord = 0
         if (natom1.eq.2) ncoord = 1
         do icoord = 1 , ncoord
            call LineRead (21)
            intcoor(icoord) = word
         enddo
         rewind(21)
      endif

      do while (WORD.NE.'NTAU')
         call LineRead (21)
      enddo
      read (21,*) ntautemp

      if (ntautemp.ne.0) then
         read (21,*)
         do itau = 1 , ntautemp
            call LineRead (21)
cc modified intcoor vector index
            intcoor(natom1p+itau-ntau1) = word
c            read (21,*) intcoor(natom1p+itau-ntau1),taumn,taumx
c            write (*,*) intcoor(natom1p+itau-ntau1),taumn,taumx
         enddo
      endif
      rewind(21)

cc read coordinate labels for reacting atom

      do while (WORD.NE.'CHARGE')
         call LineRead (21)
      enddo
      read (21,*)
      do iatom = 1 , natomt1
         if(iatom.eq.ireact)then
            read(21,*)cname1,cjunk,adislab
         else
            call LineRead (21)
         endif
      enddo
      close (unit=21,status='keep')
c      write(*,*)'react at is ',adislab
c      write(*,*)'ireact at is ',ireact
c      stop

cc convert iatom parameters to upper case

      open (unit=99,status='unknown')
      write(99,*)adislab
      rewind (99)
      call LineRead (99)
      adislab=word
      close (99)

cc now read optimized coordinates

      if( iopt.ne.0) then
         open (unit=99,status='unknown')
         write(99,1200)iopt
         rewind (99)
         read (99,'(A40)') filename
         close (99)
      else
         filename='./output/reac1_opt.out'
      endif

 1200 format ("./output/reac1_opt_"I0.2".out")

      open(unit=23,file=filename,status='unknown')

      read (23,*)cjunk,cjunk,iminst
      
      do icoord = 1 , natom1p
         read (23,*) xinti(icoord)
      enddo

      close (unit=23,status='keep')

  
cc now move distance coordinate for reacting atom to last position

      do icoord = 1 , natom1p
         if(intcoor(icoord).eq.adislab)then
            icoordmove=icoord
         endif
      enddo

      index=0
      do icoord = 1 , natom1p
         if(icoord.ne.icoordmove)then
            index=index+1
            xintt(index) = xinti(icoord)
            intcoort(index) = intcoor(icoord)
         endif
         if(icoord.eq.icoordmove)then
            dstart=xinti(icoord)
         endif
      enddo
      xintt(natom1p)=dstart
      intcoort(natom1p)=adislab

c      stop

      do icoord = 1 , natom1p
         xinti(icoord)=xintt(icoord)
         intcoor(icoord) = intcoort(icoord)
      enddo

      do icoord = 1 , natom1p
         write(6,*)icoord, xinti(icoord)
         write(6,*)icoord, intcoor(icoord)
      enddo

      write (116,*) 'past z-matrix'

      close (unit=22,status='keep')
      close (unit=24,status='keep')
      close (unit=116,status='keep')

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine stoichiometry(filename,stoichname)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension natomnumb(natommx)

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      character*20 filename
      character*20 stoichname
      character*60 command1

      character*2 atomlabel(natommx)
      character*2 aname
      character*30 cjunk
      character*30 gmem

      include 'filcomm.f'

      write(*,*)'filename is ',filename
      open (unit=116,file=filename,status='old')
      read(116,*)cjunk
      read(116,*)cjunk
      read(116,*)cjunk
      read(116,*)cjunk,natom
c      write(*,*)'natom is ', natom
      numH=0
      numB=0
      numC=0
      numN=0
      numO=0
      numF=0
      numNe=0
      numAl=0
      numSi=0
      numP=0
      numS=0
      numCl=0
      numAr=0
      numGa=0
      numGe=0
      numAs=0
      numSe=0
      numBr=0
      numKr=0
      numIn=0
      numSn=0
      numSb=0
      numTe=0
      numI=0
      numXe=0
      do j=1,natom
         READ (116,*)aname,cjunk,cjunk,cjunk
         if (aname.eq.'H') numH=numH+1
         if (aname.eq.'B') numB=numB+1
         if (aname.eq.'C') numC=numC+1
         if (aname.eq.'N') numN=numN+1
         if (aname.eq.'O') numO=numO+1
         if (aname.eq.'F') numF=numF+1
         if (aname.eq.'Ne') numNe=numNe+1
         if (aname.eq.'Al') numAl=numAl+1
         if (aname.eq.'Si') numSi=numSi+1
         if (aname.eq.'P') numP=numP+1
         if (aname.eq.'S') numS=numS+1
         if (aname.eq.'Cl') numCl=numCl+1
         if (aname.eq.'Ar') numAr=numAr+1
         if (aname.eq.'Ga') numGa=numGa+1
         if (aname.eq.'Ge') numGe=numGe+1
         if (aname.eq.'As') numAs=numAs+1
         if (aname.eq.'Se') numSe=numSe+1
         if (aname.eq.'Br') numBr=numBr+1
         if (aname.eq.'Kr') numKr=numKr+1
         if (aname.eq.'In') numIn=numIn+1
         if (aname.eq.'Sn') numSn=numSn+1
         if (aname.eq.'Sb') numSb=numSb+1
         if (aname.eq.'Te') numTe=numTe+1
         if (aname.eq.'I') numI=numI+1
         if (aname.eq.'Xe') numXe=numXe+1
      enddo
      open(unit=117,file='stoich.res',status='unknown')

      indnum=1
      if(numH.ne.0) then
         atomlabel(indnum)='H'
         natomnumb(indnum)=numH
         indnum=indnum+1
      endif
      if(numB.ne.0) then
         atomlabel(indnum)='B'
         natomnumb(indnum)=numB
         indnum=indnum+1
      endif
      if(numC.ne.0) then
         atomlabel(indnum)='C'
         natomnumb(indnum)=numC
         indnum=indnum+1
      endif
      if(numN.ne.0) then
         atomlabel(indnum)='N'
         natomnumb(indnum)=numN
         indnum=indnum+1
      endif
      if(numO.ne.0) then
         atomlabel(indnum)='O'
         natomnumb(indnum)=numO
         indnum=indnum+1
      endif
      if(numF.ne.0) then
         atomlabel(indnum)='F'
         natomnumb(indnum)=numH
         indnum=indnum+1
      endif
      if(numNe.ne.0) then
         atomlabel(indnum)='Ne'
         natomnumb(indnum)=numNe
         indnum=indnum+1
      endif
      if(numAl.ne.0) then
         atomlabel(indnum)='Al'
         natomnumb(indnum)=numAl
         indnum=indnum+1
      endif
      if(numSi.ne.0) then
         atomlabel(indnum)='Si'
         natomnumb(indnum)=numSi
         indnum=indnum+1
      endif
      if(numP.ne.0) then
         atomlabel(indnum)='P'
         natomnumb(indnum)=numP
         indnum=indnum+1
      endif
      if(numS.ne.0) then
         atomlabel(indnum)='S'
         natomnumb(indnum)=numS
         indnum=indnum+1
      endif
      if(numCl.ne.0) then
         atomlabel(indnum)='Cl'
         natomnumb(indnum)=numCl
         indnum=indnum+1
      endif
      if(numAr.ne.0) then
         atomlabel(indnum)='Ar'
         natomnumb(indnum)=numAr
         indnum=indnum+1
      endif
      if(numGa.ne.0) then
         atomlabel(indnum)='Ga'
         natomnumb(indnum)=numGa
         indnum=indnum+1
      endif
      if(numGe.ne.0) then
         atomlabel(indnum)='Ge'
         natomnumb(indnum)=numGe
         indnum=indnum+1
      endif
      if(numAs.ne.0) then
         atomlabel(indnum)='As'
         natomnumb(indnum)=numAs
         indnum=indnum+1
      endif
      if(numSe.ne.0) then
         atomlabel(indnum)='Se'
         natomnumb(indnum)=numSe
         indnum=indnum+1
      endif
      if(numBr.ne.0) then
         atomlabel(indnum)='Br'
         natomnumb(indnum)=numBr
         indnum=indnum+1
      endif
      if(numKr.ne.0) then
         atomlabel(indnum)='Kr'
         natomnumb(indnum)=numKr
         indnum=indnum+1
      endif
      if(numIn.ne.0) then
         atomlabel(indnum)='In'
         natomnumb(indnum)=numIn
         indnum=indnum+1
      endif
      if(numSn.ne.0) then
         atomlabel(indnum)='Sn'
         natomnumb(indnum)=numSn
         indnum=indnum+1
      endif
      if(numTe.ne.0) then
         atomlabel(indnum)='Te'
         natomnumb(indnum)=numTe
         indnum=indnum+1
      endif
      if(numI.ne.0) then
         atomlabel(indnum)='I'
         natomnumb(indnum)=numI
         indnum=indnum+1
      endif
      if(numTe.ne.0) then
         atomlabel(indnum)='Te'
         natomnumb(indnum)=numTe
         indnum=indnum+1
      endif
      indnum=indnum-1
      write(*,*)'anumb(1)= ', natomnumb(1)
      write(*,*)'indnum= ', indnum
      open(unit=101,file='stoich.dat',status='unknown')  
      write(101,201)(atomlabel(j),natomnumb(j),j=1,indnum)
      close(101)
      write(command1,202)
      call commrun(command1)
      open(unit=101,file='stoich.dat',status='unknown')  
      read(101,*)stoichname
      close(101)

 201  format(10(A2,I2))
 202  format(" sed -ie 's/ //g' stoich.dat")

      close(116)
      close(117)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bmatrix(ispecies,ifile)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

c      dimension natomnumb(natommx)
      dimension tau(ntaumx),taumn(ntaumx),taumx(ntaumx),freq(nmdmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim),xintt(3*natommx)

      dimension xintp(3*natommx),xintn(3*natommx)
      dimension xintpp(3*natommx),xintpn(3*natommx)
      dimension xintnp(3*natommx),xintnn(3*natommx)

      dimension bmat(3*natommx,3*natommx)
      dimension cmat(3*natommx,3*natommx,3*natommx)

      dimension coox(natommx),cooy(natommx),cooz(natommx)
      dimension cooxn(natommx),cooyn(natommx),coozn(natommx)
      dimension cooxp(natommx),cooyp(natommx),coozp(natommx)

      dimension cooxpp(natommx),cooypp(natommx),coozpp(natommx)
      dimension cooxnp(natommx),cooynp(natommx),cooznp(natommx)
      dimension cooxpn(natommx),cooypn(natommx),coozpn(natommx)
      dimension cooxnn(natommx),cooynn(natommx),cooznn(natommx)

      dimension cooxyz(3*natommx)
      dimension ibconn(natommx),iaconn(natommx),idconn(natommx)
      dimension idummy(natommx),iangind(3*natommx)
      dimension ianginda(3*natommx),iangindd(3*natommx)
      dimension idihed(natommx)
      dimension bval(natommx),aval(natommx),dval(natommx)


      LOGICAL leof,lsec,ltit
      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      character*20 bname(natommx),anname(natommx),dname(natommx)
     $ ,atname(natommx),bconnt(natommx),aconnt(natommx),dconnt(natommx)

      character*20 filename
      character*20 stoichname
      character*60 command1

      character*60 atomlabel(natommx)
      character*2 aname
      character*1 aname1(natommx)
      character*30 cjunk
      character*30 angsub1,angsub2
      character*30 nameout,nameoutc
      character*30 gmem
      character*30 intcoor(3*natommx)
      character*30 intcoort(3*natommx)
      character*20 bislab(ntaumx)

      include 'filcomm.f'

      if (ispecies.eq.1) then
         open (unit=15,file='./data/reac1.dat',status='old')
         open (unit=17,file='./output/reac1_opt.out',status='unknown')
         if(ifile.eq.0)then
            open (unit=18,file='./geoms/reac1_l1.xyz',status='unknown')
         endif
         nameout='reac1_bmat'
         nameoutc='reac1_cmat'
         inp_type=1
      endif
      if (ispecies.eq.2) then
         open (unit=15,file='./data/reac2.dat',status='old')
         open (unit=17,file='./output/reac2_opt.out',status='unknown')
         if(ifile.eq.0)then
            open (unit=18,file='./geoms/reac2_l1.xyz',status='unknown')
         endif
         nameout='reac2_bmat'
         nameoutc='reac2_cmat'
         inp_type=1
      endif
      if (ispecies.eq.3) then
         open (unit=15,file='./data/prod1.dat',status='old')
         open (unit=17,file='./output/prod1_opt.out',status='unknown')
         if(ifile.eq.0)then
            open (unit=18,file='./geoms/prod1_l1.xyz',status='unknown')
         endif
         nameout='prod1_bmat'
         nameoutc='prod1_cmat'
         inp_type=1
      endif
      if (ispecies.eq.4) then
         open (unit=15,file='./data/prod2.dat',status='old')
         open (unit=17,file='./output/prod2_opt.out',status='unknown')
         if(ifile.eq.0)then
            open (unit=18,file='./geoms/prod2_l1.xyz',status='unknown')
         endif
         nameout='prod2_bmat'
         nameoutc='prod2_cmat'
         inp_type=1
      endif
      if (ispecies.eq.5) then
         open (unit=15,file='./data/wellr.dat',status='old')
         open (unit=17,file='./output/wellr_opt.out',status='unknown')
         if(ifile.eq.0)then
            open (unit=18,file='./geoms/wellr_l1.xyz',status='unknown')
         endif
         nameout='wellr_bmat'
         nameoutc='wellr_cmat'
         inp_type=2
      endif
      if (ispecies.eq.6) then
         open (unit=15,file='./data/wellp.dat',status='old')
         open (unit=17,file='./output/wellp_opt.out',status='unknown')
         if(ifile.eq.0)then
            open (unit=18,file='./geoms/wellp_l1.xyz',status='unknown')
         endif
         nameout='wellp_bmat'
         nameoutc='wellp_cmat'
         inp_type=1
      endif
      if (ispecies.eq.51) then
         open (unit=15,file='./data/wellp.dat',status='old')
         if(igeom_wellr.eq.2)then
            open (unit=17,file='./output/ts_opt.out',status='unknown')
         else
            open(unit=17,file='./output/wellp_opt.out',status='unknown')
         endif
         if(ifile.eq.0)then
            open (unit=18,file='./geoms/wellr_l1.xyz',status='unknown')
         endif
         nameout='wellr_bmat'
         nameoutc='wellr_cmat'
         inp_type=2
      endif
      if (ispecies.eq.61) then
         open (unit=15,file='./data/wellp.dat',status='old')
         if(igeom_wellp.eq.2)then
            open (unit=17,file='./output/ts_opt.out',status='unknown')
         else
            open(unit=17,file='./output/wellp_opt.out',status='unknown')
         endif
         if(ifile.eq.0)then
            open (unit=18,file='./geoms/wellp_l1.xyz',status='unknown')
         endif
         nameout='wellp_bmat'
         nameoutc='wellp_cmat'
         inp_type=2
      endif
      if (ispecies.eq.0) then
         open (unit=15,file='./data/ts.dat',status='old')
         open (unit=17,file='./output/ts_opt.out',status='unknown')
         if(ifile.eq.0)then
            open (unit=18,file='./geoms/tsgta_l1.xyz',status='unknown')
         endif
         nameout='tsgta_bmat'
         nameoutc='tsgta_cmat'
         inp_type=2
      endif

cc initialize parameters and files

      if(ifile.eq.1)then         
         open (unit=18,file='geom_bmat.xyz',status='unknown')
      endif

cc first read blocks that are common for input type 1 and 2

      if (idebug.ge.2) write (66,*) ' starting zmat input'
      call LineRead (0)

cc read specific input for input type 1

      if (inp_type.eq.1) then

c     natomt is to account for dummy atoms
         do while (WORD.NE.'NATOM')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'natom must be defined'
               stop
            endif
         enddo
         read (15,*) natom,natomt,ilin
         if (natomt.gt.natommx) then
            write (66,*) 'natomt too large',natomt,natommx
            stop
         endif
         rewind(15)

         write (66,*) 'test0',natomt,atomlabel(1)

         do while (WORD.NE.'CHARGE')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'charge and spin must be defined'
               stop
            endif
         enddo
         read (15,*) icharge,ispin
         do iatom = 1 , natomt
            read (15,'(A60)') atomlabel(iatom)
         enddo
         rewind(15)

         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'sampling coordinates must be defined'
               stop
            endif
         enddo
         read (15,*) ntau
         rewind(15)
         call  LineRead (0)

         if (natom.ne.1) then
            do while (WORD.NE.'INTCOOR')
               call LineRead (15)
               if (WORD.EQ.'END') then
                  write (66,*) 'internal coordinates must be defined'
                  stop
               endif
            enddo
            ncoord = 3*natom-6-ntau
            if (natom.eq.1) ncoord = 0
            if (natom.eq.2) ncoord = 1
            do icoord = 1 , ncoord
               call LineRead (15)
               intcoor(icoord) = word
            enddo
            rewind(15)
         endif

         do while (WORD.NE.'NTAU')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'sampling coordinates must be defined'
               stop
            endif
         enddo
         read (15,*) ntau
         if (ntau.gt.ntaumx) then
            write (66,*) 'ntau too large',ntau,ntaumx
            stop
         endif
         if (ntau.ne.0) then
            read (15,*)
            do itau = 1 , ntau
               read (15,*) bislab(itau),taumn(itau),taumx(itau)
               write (66,*) bislab(itau),taumn(itau),taumx(itau)
               open (unit=99,status='unknown')
               rewind (99)
               write (99,*)bislab(itau)
               rewind (99)
               call LineRead (99)
               icoord=ncoord+itau
               intcoor(icoord)=WORD
               close (unit=99,status='keep')
            enddo
         endif
         rewind(15)
cc now read optimized geometry parameters
         ncoord = 3*natom-6
         if(natom.eq.2) ncoord = 1
         if(natom.eq.1) ncoord = 0
         read (17,*)
         do icoord = 1 , ncoord
            call LineRead (17)
c     xinti(icoord) = word
c     intcoori(icoord) = word
            OPEN (unit=99,status='unknown')
            REWIND (99)
            WRITE (99,1000) WORD
 1000       FORMAT (A30)
            REWIND (99)
            READ (99,*) xint(icoord)
            close (unit=99,status='keep')
         enddo
         close (unit=17,status='keep')
         close (unit=15,status='keep')

      else if (inp_type.eq.2) then

cc here we assume that the TS is not linear

         ilin=0
         ilin1=0
         ilin2=0

cc now read input of type 2
         do while (WORD.NE.'CHARGE')
            call LineRead (15)
            if (WORD.EQ.'END') then
               write (66,*) 'charge and spin must be defined'
               stop
            endif
         enddo
         read (15,*) icharge,ispin
         rewind(15)
         call  LineRead (0)

         do while (WORD.NE.'END')
            call LineRead (15)
            if (WORD.EQ.'INTCOORD') then
               read (15,*) nintcoord
               nlbend=0
               do j=1,nintcoord
                  call LineRead (15)
                  if (WORD.EQ.'LBEND')then
                     nlbend=nlbend+1
                     angsub1=WORD2
                     angsub2=WORD3
                  endif
               enddo
            endif
         enddo
         rewind(15)
         call  LineRead (0)

c         write(*,*)'angsub1 ',angsub1
c         write(*,*)'angsub2 ',angsub2
c         stop

         open (unit=25,file='./data/reac1.dat',status='old')

         do while (WORD.NE.'NTAU')
            call LineRead (25)
            if (WORD.EQ.'END') then
               write (66,*) 'sampling coords of reac1 must be defined'
               stop
            endif
         enddo
         read (25,*) ntau1
         rewind 25

c     natomt is to account for dummy atoms
         do while (WORD.NE.'NATOM')
            call LineRead (25)
            if (WORD.EQ.'END') then
               write (66,*) 'natom in reac1 must be defined'
               stop
            endif
         enddo
         read (25,*) natom1,natomt1,ilin1
         close (25)

cc get data from react2 file

         if(iabs.eq.1.or.iadd.eq.1)then
            open (unit=25,file='./data/reac2.dat',status='old')
            
            do while (WORD.NE.'NTAU')
               call LineRead (25)
               if (WORD.EQ.'END') then
                  write (66,*) 'samp coords of reac2 must be defined'
                  stop
               endif
            enddo
            read (25,*) ntau2
            rewind 25

c     natomt is to account for dummy atoms
            do while (WORD.NE.'NATOM')
               call LineRead (25)
               if (WORD.EQ.'END') then
                  write (66,*) 'natom must be defined'
                  stop
               endif
            enddo
            read (25,*) natom2,natomt2,ilin2
            close (25)
         endif
cc now we can determine the total number of atoms for the TS/wellr/wellp
         if(iadd.eq.1.or.iabs.eq.1)then
            natom = natom1+natom2
         else if (iiso.eq.1.or.ibeta.eq.1)then
            natom = natom1
         endif
cc modified         if (iadd.eq.1) natomt = natomt1+natomt2
         if (iadd.eq.1) natomt = natomt1+natomt2
         if (iabs.eq.1) natomt = natomt1+natomt2+1
         if (iiso.eq.1) natomt = natomt1
         if (ibeta.eq.1) natomt = natomt1

cc if linear molecule reacting with atom in abs assume TS is linear
         if(iabs.eq.1.and.natom2.eq.1.and.ispecies.eq.0.
     $      and.ilin1.eq.1) then
            ilin=1
         endif

cc if linear molecule reacting with linear radical in abs assume TS is linear
         if(iabs.eq.1.and.ilin2.eq.1.and.ispecies.eq.0.and
     $      .ilin1.eq.1)then
            ilin=1
         endif
c        natomt = natomt1+natomt2+1
         
         close (unit=15,status='keep')

c natomt is to account for dummy atoms
         if (natomt.gt.natommx) then
            write (66,*) 'natomt too large',natomt,natommx
            stop
         endif

c gaussian com file data
         read (17,*)
         if (idebug.ge.2) write (6,*) ' starting gaussian input'
         do iatom = 1 , natomt
            read (17,'(A60)') atomlabel(iatom)
         enddo

cc read coordinate names

         ncoord = 3*natom-6

         do iint = 1 , ncoord
            read (17,*) intcoor(iint),xint(iint)
         enddo
         close (unit=17,status='keep')
      endif

      do j=1,natomt
         write(*,*)'atomlabel is ',atomlabel(j)
      enddo
      do j=1,ncoord
         write(*,*)'coord is ',intcoor(j),xint(j)
      enddo

cc initialize vectors
      do j=1,ncoord 
         bname(j)=''
         anname(j)=''
         dname(j)=''
         bval(j)=0.
         aval(j)=0.
         dval(j)=0.
      enddo
      do j=1,natomt
         ibconn(j)=0
         iaconn(j)=0
         idconn(j)=0
      enddo
cc
      call read_zmat(atomlabel,natom,natomt,intcoor,bislab,ibconn,
     $ iaconn,idconn,bname,anname,dname,atname,idummy,isited,jsited,
     $ ksited,bconnt,aconnt,dconnt)
 
c      write(*,*)'bname is ',bname(2)
c      write(*,*)'anname is ',anname(4)
c      do j=1,natomt
c         write(*,*) 'iaconn ',j,iaconn(j)
c         write(*,*) 'idconn ',j,idconn(j)
c      enddo
c      stop
cc this does not account for dummy atoms defined with  given quantities

      do k=1,natomt
         do j=1,ncoord
            if(intcoor(j).eq.bname(k))bval(k)=xint(j)
            if(intcoor(j).eq.anname(k))aval(k)=xint(j)
            if(intcoor(j).eq.dname(k))dval(k)=xint(j)
         enddo
      enddo
c      write(*,*)'aval4 is ',aval(4)


cc all that follows is not necessary. we read the intial geometry
cc I leave (coomented) it as it can be useful sooner or later

cc   define starting xyz geometry. 
cc   convention: atom 1 is is 0 0 0
cc   atom 2 bd 0 0
cc   atom 3 on xy plane 

c      coox(1)=0.
c      cooy(1)=0.
c      cooz(1)=0.
c      coox(2)=bval(2)
c      cooy(2)=0.
c      cooz(2)=0.
c      if(ibconn(3).eq.2) then
c         coox(3)=bval(2)-bval(3)*cos(aval(3)/180.*pigr)
c         cooy(3)=bval(3)*sin(aval(3)/180.*pigr)
c         cooz(3)=0.
c      else if(ibconn(3).eq.1)then
c         coox(3)=bval(3)*cos(aval(3)/180.*pigr)
c         cooy(3)=bval(3)*sin(aval(3)/180.*pigr)
c         cooz(3)=0.
c      else
c         write(*,*)'I am lost with this Z-Mat'
c         write(*,*)'in subroutine Bmatrix'
c         write(*,*)'for species ',ispecies
c         stop
c      endif
c      xa=0.
c      ya=0.
c      za=0.
c      do j=4,natomt
c         call zmat_to_xyz(xa,ya,za,coox(ibconn(j)),cooy(ibconn(j)),
c     $ cooz(ibconn(j)),coox(iaconn(j)),cooy(iaconn(j)),
c     $ cooz(iaconn(j)),
c     $ coox(idconn(j)),cooy(idconn(j)),cooz(idconn(j)),
c     $ bval(j),aval(j),dval(j))
c         coox(j)=xa
c         cooy(j)=ya
c         cooz(j)=za
c         write(*,*) j,coox(j),cooy(j),cooz(j)
c         write(*,*)
c      enddo


cc read the zmat from the level 1 output
      read(18,*)natomcar
      read(18,*)cjunk
      do j=1,natomcar
         read(18,*)cjunk,coox(j),cooy(j),cooz(j)
c         write(*,*)cjunk,coox(j),cooy(j),cooz(j)
      enddo
      close(18)
cc we do not need a dummy atom check if there are no dummy atoms in the cart matrix
c      if(natomcar.eq.natom)then
c         do j=1,natomcar
c            idummy(j)=0
c         enddo
c      endif
c      natomt=natom

cc this is a check that we can compute the xyz matrix and convert back to internal
c      write(*,*)'ok1 '
c      stop
      ntau=0
      ifilu=1
      call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,coox,cooy,cooz,xinti,tauopt,
     $     ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)
      if(xint(1).eq.0)then
         write(*,*)'failed in updating geometry'
         write(*,*)'probably error of g09, check output'
         write(*,*)'stopping now'
         stop
      endif

c      write(*,*)'ok2 '

      open(unit=10,file='temp1.xyz',status='unknown')
      do j=1,ncoord
         write(10,*)intcoor(j),xint(j),xinti(j),abs(xint(j)-xinti(j))
      enddo
      close(10)

      do j=1,ncoord
         xint(j)=xinti(j)
      enddo

c      stop

      open(unit=10,file='temp2.xyz',status='unknown')
      write(10,*)natom
      write(10,*)'test'
      inda=0
      do j=1,natom
         aname1(j)=atname(j)
         if(idummy(j).ne.1)then
            inda=inda+1
      endif
      write(10,*)aname1(inda),coox(j),cooy(j),cooz(j)
      enddo
      close(10)

c      stop

cc now we update the coordinates using those etransformed in our procedure for consistentcy
cc with dimension of dihedral angles

      do j=1,ncoord
         xint(j)=xinti(j)
      enddo

cc here I inizialize indexes to check if int coor is an angle or a dihedral
cc

      do k=1,ncoord
         iangind(k)=0
         ianginda(k)=0
         iangindd(k)=0
      enddo


      do k=1,ncoord
         do j=1,ncoord
            if(intcoor(k).eq.anname(j))iangind(k)=1
            if(intcoor(k).eq.anname(j))ianginda(k)=1
            if(intcoor(k).eq.dname(j))iangind(k)=1
            if(intcoor(k).eq.dname(j))iangindd(k)=1
         enddo
      enddo


cc this is to check
c      do ik=1,natom
c         cooxp(ik)=coox(ik)
c         cooyp(ik)=cooy(ik)
c         coozp(ik)=cooz(ik)
c      enddo
c      coozp(5)=coozp(5)+0.01

c      call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
c     $    ,idconn,bname,anname,dname,atname,cooxp,cooyp,coozp,xintp,
c     $      tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,
c     $      ifilu)

c      do j=1,ncoord
c         write(*,*)intcoor(j),xint(j),xintp(j)
c      enddo
c      stop





cc now we can compute the B matrix using central difference
cc Bik=dqi/dxk

      deltax=0.01
      deltay=0.01
      deltaz=0.01
      kind=0
      kaind=0
      ifilu=0

      do ij=1,ncoord
         do ik=1,3*natom
            bmat(ij,ik)=0.
         enddo
      enddo

c      do k=1,natomt
c         write(*,*)'idummy k atom is ', idummy(k)
c      enddo
c      stop

      do k=1,natomt
         if(idummy(k).ne.1)then
            kind=kind+1
            kaind=kaind+1
            write(*,*)'k atom is ', k
            write(*,*)'kaind is ', kaind
c           write(*,*)'idummy k atom is ', idummy(k)
            write(*,*)'kind is ', kind
            do ij=1,ncoord
               xintp(ij)=0.
            enddo

cc first perturb x + deltax
c            write(*,*)'ok1 '
c            do ik=1,natomt
            do ik=1,natom
               cooxp(ik)=coox(ik)
               cooyp(ik)=cooy(ik)
               coozp(ik)=cooz(ik)
            enddo
            cooxp(kaind)=cooxp(kaind)+deltax

            call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,cooxp,cooyp,coozp,xintp,
     $  tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)

cc then perturb x - deltax
c            write(*,*)'ok2 '
c            do ik=1,natomt
            do ik=1,natom
               cooxn(ik)=coox(ik)
               cooyn(ik)=cooy(ik)
               coozn(ik)=cooz(ik)
            enddo
            cooxn(kaind)=cooxn(kaind)-deltax

            call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,cooxn,cooyn,coozn,xintn,
     $  tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)

cc we can now compute the first component of the Bmatrix

            do i=1,ncoord
               if(abs(xintp(i)-xintn(i)).gt.300)then
                  if(xintn(i).lt.0)then
                     xintn(i)=xintn(i)+360.
                  else if(xintn(i).gt.0)then
                     xintn(i)=xintn(i)-360.
                  endif
                  if(abs(xintp(i)-xintn(i)).gt.300)then
                     write(*,*)'something did not work here'
                     write(*,*)'kind= ',kind
                     write(*,*)'intocoor= ',i
                     write(*,*)'xintp= ',xintp(i)
                     write(*,*)'xintn= ',xintn(i)
                     stop
                  endif
               endif
               bmat(i,kind)=(xintp(i)-xintn(i))/2./deltax
c               bmat(i,kind)=(xintp(i)-xint(i))/deltax
c                     write(*,*)'kind= ',kind
c                     write(*,*)'intocoor= ',i
c                     write(*,*)'xintp= ',xintp(i)
c                     write(*,*)'xintn= ',xintn(i)
            enddo

c            write(*,*)'ok3 '
cc we now replicate for y

cc first perturb y + deltay
c            do ik=1,natomt
            do ik=1,natom
               cooxp(ik)=coox(ik)
               cooyp(ik)=cooy(ik)
               coozp(ik)=cooz(ik)
            enddo
            cooyp(kaind)=cooyp(kaind)+deltay
c            write(*,*)'cooyp= ',cooyp(k)

            call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,cooxp,cooyp,coozp,xintp,
     $  tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)

cc then perturb y - deltay
c            do ik=1,natomt
            do ik=1,natom
               cooxn(ik)=coox(ik)
               cooyn(ik)=cooy(ik)
               coozn(ik)=cooz(ik)
            enddo
            cooyn(kaind)=cooyn(kaind)-deltay
c            write(*,*)'cooyn= ',cooyn(k)

            call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,cooxn,cooyn,coozn,xintn,
     $  tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)

cc we can now compute the second(y) component of the Bmatrix for atom k

            kind=kind+1
c           write(*,*)'kind= ',kind
           do i=1,ncoord
               if(abs(xintp(i)-xintn(i)).gt.300)then
                  if(xintn(i).lt.0)then
                     xintn(i)=xintn(i)+360.
                  else if(xintn(i).gt.0)then
                     xintn(i)=xintn(i)-360.
                  endif
                  if(abs(xintp(i)-xintn(i)).gt.300)then
                     write(*,*)'something did not work here'
                     write(*,*)'kind= ',kind
                     write(*,*)'intocoor= ',i
                     write(*,*)'xintp= ',xintp(i)
                     write(*,*)'xintn= ',xintn(i)
                     stop
                  endif
               endif
               bmat(i,kind)=(xintp(i)-xintn(i))/2./deltay
c               bmat(i,kind)=(xintp(i)-xint(i))/deltay
c               write(*,*)'intocoor= ',i
c               write(*,*)'xintp= ',xintp(i)
c               write(*,*)'xintn= ',xintn(i)
            enddo

cc we now replicate for z

cc first perturb z + deltaz
c            do ik=1,natomt
            do ik=1,natom
               cooxp(ik)=coox(ik)
               cooyp(ik)=cooy(ik)
               coozp(ik)=cooz(ik)
            enddo
            coozp(kaind)=coozp(kaind)+deltaz

            call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,cooxp,cooyp,coozp,xintp,
     $  tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)

cc then perturb z - deltaz
c            do ik=1,natomt
            do ik=1,natom
               cooxn(ik)=coox(ik)
               cooyn(ik)=cooy(ik)
               coozn(ik)=cooz(ik)
            enddo
            coozn(kaind)=coozn(kaind)-deltaz

            call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,cooxn,cooyn,coozn,xintn,
     $  tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)

cc we can now compute the third(z) component of the Bmatrix for atom k

            kind=kind+1
            do i=1,ncoord
               if(abs(xintp(i)-xintn(i)).gt.300)then
                  if(xintn(i).lt.0)then
                     xintn(i)=xintn(i)+360.
                  else if(xintn(i).gt.0)then
                     xintn(i)=xintn(i)-360.
                  endif
                  if(abs(xintp(i)-xintn(i)).gt.300)then
                     write(*,*)'something did not work here'
                     write(*,*)'kind= ',kind
                     write(*,*)'intocoor= ',i
                     write(*,*)'xintp= ',xintp(i)
                     write(*,*)'xintn= ',xintn(i)
                     stop
                  endif
               endif
c               write(*,*)'xintp= ',xintp(i)
c               write(*,*)'xintn= ',xintn(i)
               bmat(i,kind)=(xintp(i)-xintn(i))/2./deltaz
c               bmat(i,kind)=(xintp(i)-xint(i))/deltaz
            enddo
         endif
      enddo
c      write(*,*)'kind is ',kind
c      write(*,*)'kaind is ',kaind


c      do k=1,ncoord
c         write(*,*)'iangind(k)',k,iangind(k)
c      enddo

      do k=1,ncoord
         do j=1,3*natom
            if(iangind(k).eq.1)bmat(k,j)=bmat(k,j)/180.*pigr
         enddo
      enddo

      if(nintcoord.ne.0)then

         do j=1,ncoord
c         write(*,*)'intcoor(j) ',intcoor(j)
c         write(*,*)'angsub1 ',angsub1
            if(intcoor(j).eq.angsub1)ibmatsub1=j
            if(intcoor(j).eq.angsub2)ibmatsub2=j
         enddo

         iatsub1=0
         iatsub2=0
         do j=1,natomt
            if(anname(j).eq.angsub1)iatsub1=j
            if(dname(j).eq.angsub2)iatsub2=j
         enddo


c         write(*,*)'i angsub1 ',ibmatsub1
c         write(*,*)'i angsub2 ',ibmatsub2
c         write(*,*)'iatom angsub1 ',iatsub1
c         write(*,*)'iatom angsub2 ',iatsub2
c         write(*,*)'ibconn ',ibconn(iatsub1)
c         write(*,*)'iaconn ',iaconn(iatsub1)
c         
c     now update iangsub1 bmat component
         iatom=1
         iinda=0
         bcentx=0.
         do j=1,3*natom
            iinda=iinda+1
            if(iinda.eq.4)then
               iinda=1
               iatom=iatom+1
            endif
            if((iatom.ne.iatsub1).and.iatom.ne.ibconn(iatsub1).and.
     +          iatom.ne.iaconn(iatsub1))then
               bmat(ibmatsub1,j)=0.
            endif
         enddo
         bmat(ibmatsub1,(ibconn(iatsub1)-1)*3+1)=
     +       -(bmat(ibmatsub1,(iatsub1-1)*3+1)+
     +         bmat(ibmatsub1,(iaconn(iatsub1)-1)*3+1))/2.
         bmat(ibmatsub1,(ibconn(iatsub1)-1)*3+2)=
     +       -(bmat(ibmatsub1,(iatsub1-1)*3+2)+
     +         bmat(ibmatsub1,(iaconn(iatsub1)-1)*3+2))/2.
         bmat(ibmatsub1,(ibconn(iatsub1)-1)*3+3)=
     +       -(bmat(ibmatsub1,(iatsub1-1)*3+3)+
     +         bmat(ibmatsub1,(iaconn(iatsub1)-1)*3+3))/2.

cc
cc now determine displacement in orthogonal plane
cc
cc  first get components of vector on plane
         vpx=coox(iatsub1)-coox(ibconn(iatsub1))
         vpy=cooy(iatsub1)-cooy(ibconn(iatsub1))
         vpz=cooz(iatsub1)-cooz(ibconn(iatsub1))
         vpnorm=sqrt(vpx*vpx+vpy*vpy+vpz*vpz)
         vpx=vpx/vpnorm
         vpy=vpy/vpnorm
         vpz=vpz/vpnorm
cc now update components of of 2nd  bmat linear angle
cc starting values
         bm1x=bmat(ibmatsub1,(iatsub1-1)*3+1)
         bm1y=bmat(ibmatsub1,(iatsub1-1)*3+2)
         bm1z=bmat(ibmatsub1,(iatsub1-1)*3+3)
         bm2x=bmat(ibmatsub1,(ibconn(iatsub1)-1)*3+1)
         bm2y=bmat(ibmatsub1,(ibconn(iatsub1)-1)*3+2)
         bm2z=bmat(ibmatsub1,(ibconn(iatsub1)-1)*3+3)
         bm3x= bmat(ibmatsub1,(iaconn(iatsub1)-1)*3+1)
         bm3y= bmat(ibmatsub1,(iaconn(iatsub1)-1)*3+2)
         bm3z= bmat(ibmatsub1,(iaconn(iatsub1)-1)*3+3)
cc project to perpendicular plane
         bmn1x=bm1y*vpz-bm1z*vpy
         bmn1y=-bm1x*vpz+bm1z*vpx
         bmn1z=bm1x*vpy-bm1y*vpx
         bmn2x=bm2y*vpz-bm2z*vpy
         bmn2y=-bm2x*vpz+bm2z*vpx
         bmn2z=bm2x*vpy-bm2y*vpx
         bmn3x=bm3y*vpz-bm3z*vpy
         bmn3y=-bm3x*vpz+bm3z*vpx
         bmn3z=bm3x*vpy-bm3y*vpx

cc update vector.

         do j=1,3*natom
            bmat(ibmatsub2,j)=0.
         enddo
         bmat(ibmatsub2,(iatsub1-1)*3+1)=bmn1x
         bmat(ibmatsub2,(iatsub1-1)*3+2)=bmn1y
         bmat(ibmatsub2,(iatsub1-1)*3+3)=bmn1z
c         write(*,*)'test ',bmat(ibmatsub2,(iatsub1-1)*3+3)
c         write(*,*)'test2',(iatsub1-1)*3+3

         bmat(ibmatsub2,(ibconn(iatsub1)-1)*3+1)=bmn2x
         bmat(ibmatsub2,(ibconn(iatsub1)-1)*3+2)=bmn2y
         bmat(ibmatsub2,(ibconn(iatsub1)-1)*3+3)=bmn2z

         bmat(ibmatsub2,(iaconn(iatsub1)-1)*3+1)=bmn3x
         bmat(ibmatsub2,(iaconn(iatsub1)-1)*3+2)=bmn3y
         bmat(ibmatsub2,(iaconn(iatsub1)-1)*3+3)=bmn3z

      endif
c      stop
      


      open(unit=20,file='bmat.dat',status='unknown')
      write(20,*)ncoord,3*natom
      do k=1,ncoord
         write(20,102)(bmat(k,j),j=1,3*natom)
      enddo
      close(20)

      if(ifile.eq.0)then
         open (unit=99,status='unknown')
         rewind (99)
         write (99,100) nameout
         rewind (99)
         read (99,101) command1
         close(99)
         call commrun(command1)         
      endif

 100  format ('cp -f bmat.dat output/'A8,'.dat')
 101  format (A100)
 102  format (10000(1X,1PE11.4))


c      stop

cc now we can compute the C matrix using central difference
cc Cijk=d2qi/dxj/dxk

c      deltax=0.01
c      deltay=0.01
c      deltaz=0.01
      jind=0
      kind=0
      jaind=0
      kaind=0

      do j=1,natomt
         if(idummy(j).ne.1)then
 
c*******************************************************************
cc perturb xj
            jind=jind+1
            jaind=jaind+1
            kind=0
            kaind=0
            do k=1,natomt
            if(idummy(k).ne.1)then
               kaind=kaind+1

cc  perturb xj + deltaxj
cc  perturb xk + deltaxk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               cooxpp(jaind)=cooxpp(jaind)+deltax
               cooxpp(kaind)=cooxpp(kaind)+deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj - deltaxj
cc  perturb xk + deltaxk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooxnp(jaind)=cooxnp(jaind)-deltax
               cooxnp(kaind)=cooxnp(kaind)+deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj + deltaxj
cc  perturb xk - deltaxk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               cooxpn(jaind)=cooxpn(jaind)+deltax
               cooxpn(kaind)=cooxpn(kaind)-deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj - deltaxj
cc  perturb xk - deltaxk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooxnn(jaind)=cooxnn(jaind)-deltax
               cooxnn(kaind)=cooxnn(kaind)-deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)

cc we can now compute the first of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
c                  write(*,*)'xpp ',xintpp(i)
c                  write(*,*)'xnp ',xintnp(i)
c                  write(*,*)'xpn ',xintpn(i)
c                  write(*,*)'xnn ',xintnn(i)
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
c                  if(jind.eq.kind)then
c                     cmat(i,jind,kind)=(xintpp(i)-2*xint(i)+
c     $                    xintnn(i))/4./deltax/deltax
c                  else
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltax/deltax
c                  endif

c                  write(*,*)'cmat ', i,jind,kind,cmat(i,jind,kind)
c                  write(*,*)'xpp ',xintpp(i)
c                  write(*,*)'xnp ',xintnp(i)
c                  write(*,*)'xpn ',xintpn(i)
c                  write(*,*)'xnn ',xintnn(i)
               enddo
c               stop
cc we now replicate for y displacement of k
cc  perturb xj + deltaxj
cc  perturb yk + deltayk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               cooxpp(jaind)=cooxpp(jaind)+deltax
               cooypp(kaind)=cooypp(kaind)+deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj - deltaxj
cc  perturb yk + deltayk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooxnp(jaind)=cooxnp(jaind)-deltax
               cooynp(kaind)=cooynp(kaind)+deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj + deltaxj
cc  perturb yk - deltayk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               cooxpn(jaind)=cooxpn(jaind)+deltax
               cooypn(kaind)=cooypn(kaind)-deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj - deltaxj
cc  perturb yk - deltayk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooxnn(jaind)=cooxnn(jaind)-deltax
               cooynn(kaind)=cooynn(kaind)-deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc we can now compute the second of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
c                  if(jind.eq.kind)then
c                     cmat(i,jind,kind)=(xintpp(i)-2*xint(i)+
c     $                    xintnn(i))/4./deltax/deltay
c                  else
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltax/deltay
c                  endif
               enddo

cc we now replicate for z displacement of k
cc  perturb xj + deltaxj
cc  perturb zk + deltazk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               cooxpp(jaind)=cooxpp(jaind)+deltax
               coozpp(kaind)=coozpp(kaind)+deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj - deltaxj
cc  perturb zk + deltazk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooxnp(jaind)=cooxnp(jaind)-deltax
               cooznp(kaind)=cooznp(kaind)+deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj + deltaxj
cc  perturb zk - deltazk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               cooxpn(jaind)=cooxpn(jaind)+deltax
               coozpn(kaind)=coozpn(kaind)-deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj - deltaxj
cc  perturb zk - deltazk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooxnn(jaind)=cooxnn(jaind)-deltax
               cooznn(kaind)=cooznn(kaind)-deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc we can now compute the third of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
c                  if(jind.eq.kind)then
c                     cmat(i,jind,kind)=(xintpp(i)-2*xint(i)+
c     $                    xintnn(i))/4./deltax/deltaz
c                  else
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltax/deltaz
c                  endif
               enddo
            endif
            enddo
c*******************************************************************
cc perturb yj
            jind=jind+1
            kind=0
            kaind=0
            do k=1,natomt
            if(idummy(k).ne.1)then
               kaind=kaind+1
cc  perturb yj + deltayj
cc  perturb xk + deltaxk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               cooypp(jaind)=cooypp(jaind)+deltay
               cooxpp(kaind)=cooxpp(kaind)+deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj - deltayj
cc  perturb xk + deltaxk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooynp(jaind)=cooynp(jaind)-deltay
               cooxnp(kaind)=cooxnp(kaind)+deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj + deltayj
cc  perturb xk - deltaxk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               cooypn(jaind)=cooypn(jaind)+deltay
               cooxpn(kaind)=cooxpn(kaind)-deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj - deltayj
cc  perturb xk - deltaxk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooynn(jaind)=cooynn(jaind)-deltay
               cooxnn(kaind)=cooxnn(kaind)-deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)

cc we can now compute the fourth of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
c                  if(jind.eq.kind)then
c                     cmat(i,jind,kind)=(xintpp(i)-2*xint(i)+
c     $                    xintnn(i))/4./deltay/deltax
c                  else
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltay/deltax
c                  endif
               enddo

cc we now replicate for y displacement of k
cc  perturb yj + deltayj
cc  perturb yk + deltayk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               cooypp(jaind)=cooypp(jaind)+deltay
               cooypp(kaind)=cooypp(kaind)+deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj - deltayj
cc  perturb yk + deltayk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooynp(jaind)=cooynp(jaind)-deltay
               cooynp(kaind)=cooynp(kaind)+deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj + deltayj
cc  perturb yk - deltayk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               cooypn(jaind)=cooypn(jaind)+deltay
               cooypn(kaind)=cooypn(kaind)-deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj - deltayj
cc  perturb yk - deltayk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooynn(jaind)=cooynn(jaind)-deltay
               cooynn(kaind)=cooynn(kaind)-deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc we can now compute the fifth of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
c                  if(jind.eq.kind)then
c                     cmat(i,jind,kind)=(xintpp(i)-2*xint(i)+
c     $                    xintnn(i))/4./deltay/deltay
c                  else
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltay/deltay
c                  endif
               enddo

cc we now replicate for z displacement of k
cc  perturb yj + deltayj
cc  perturb zk + deltazk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               cooypp(jaind)=cooypp(jaind)+deltay
               coozpp(kaind)=coozpp(kaind)+deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj - deltayj
cc  perturb zk + deltazk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooynp(jaind)=cooynp(jaind)-deltay
               cooznp(kaind)=cooznp(kaind)+deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj + deltayj
cc  perturb zk - deltazk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               cooypn(jaind)=cooypn(jaind)+deltay
               coozpn(kaind)=coozpn(kaind)-deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj - deltayj
cc  perturb zk - deltazk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooynn(jaind)=cooynn(jaind)-deltay
               cooznn(kaind)=cooznn(kaind)-deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc we can now compute the sixth of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
c                  if(jind.eq.kind)then
c                     cmat(i,jind,kind)=(xintpp(i)-2*xint(i)+
c     $                    xintnn(i))/4./deltay/deltaz
c                  else
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltay/deltaz
c                  endif
               enddo
            endif
            enddo

c*******************************************************************
cc perturb zj
            jind=jind+1
            kind=0
            kaind=0
            do k=1,natomt
            if(idummy(k).ne.1)then
               kaind=kaind+1
cc  perturb zj + deltazj
cc  perturb xk + deltaxk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               coozpp(jaind)=coozpp(jaind)+deltaz
               cooxpp(kaind)=cooxpp(kaind)+deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj - deltazj
cc  perturb xk + deltaxk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooznp(jaind)=cooznp(jaind)-deltaz
               cooxnp(kaind)=cooxnp(kaind)+deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj + deltazj
cc  perturb xk - deltaxk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               coozpn(jaind)=coozpn(jaind)+deltaz
               cooxpn(kaind)=cooxpn(kaind)-deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj - deltazj
cc  perturb xk - deltaxk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooznn(jaind)=cooznn(jaind)-deltaz
               cooxnn(kaind)=cooxnn(kaind)-deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)

cc we can now compute the seventh of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
c                  if(jind.eq.kind)then
c                     cmat(i,jind,kind)=(xintpp(i)-2*xint(i)+
c     $                    xintnn(i))/4./deltaz/deltax
c                  else
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltaz/deltax
c                  endif
               enddo

cc we now replicate for y displacement of k
cc  perturb zj + deltazj
cc  perturb yk + deltayk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               coozpp(jaind)=coozpp(jaind)+deltaz
               cooypp(kaind)=cooypp(kaind)+deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj - deltazj
cc  perturb yk + deltayk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooznp(jaind)=cooznp(jaind)-deltaz
               cooynp(kaind)=cooynp(kaind)+deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj + deltazj
cc  perturb yk - deltayk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               coozpn(jaind)=coozpn(jaind)+deltaz
               cooypn(kaind)=cooypn(kaind)-deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj - deltazj
cc  perturb yk - deltayk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooznn(jaind)=cooznn(jaind)-deltaz
               cooynn(kaind)=cooynn(kaind)-deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc we can now compute the eigth of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
c                  if(jind.eq.kind)then
c                     cmat(i,jind,kind)=(xintpp(i)-2*xint(i)+
c     $                    xintnn(i))/4./deltaz/deltay
c                  else
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltaz/deltay
c                  endif
               enddo

cc we now replicate for z displacement of k
cc  perturb zj + deltazj
cc  perturb zk + deltazk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               coozpp(jaind)=coozpp(jaind)+deltaz
               coozpp(kaind)=coozpp(kaind)+deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj - deltazj
cc  perturb zk + deltazk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooznp(jaind)=cooznp(jaind)-deltaz
               cooznp(kaind)=cooznp(kaind)+deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj + deltazj
cc  perturb zk - deltazk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               coozpn(jaind)=coozpn(jaind)+deltaz
               coozpn(kaind)=coozpn(kaind)-deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj - deltazj
cc  perturb zk - deltazk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooznn(jaind)=cooznn(jaind)-deltaz
               cooznn(kaind)=cooznn(kaind)-deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc we can now compute the ninth of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
c                  if(jind.eq.kind)then
c                     cmat(i,jind,kind)=(xintpp(i)-2*xint(i)+
c     $                    xintnn(i))/4./deltaz/deltaz
c                  else
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltaz/deltaz
c                  endif
               enddo
            endif
            enddo
         endif
         write(*,*)'jind is ',jind
      enddo
      write(*,*)'kind is ',kind
      write(*,*)'jind is ',jind

      do k=1,ncoord
         do j=1,3*natom
            do i=1,3*natom
               if(iangind(k).eq.1)cmat(k,j,i)=cmat(k,j,i)/180.*pigr
            enddo
         enddo
      enddo


      open(unit=20,file='cmat.dat',status='unknown')
      write(20,*)ncoord,3*natom
      do i=1,ncoord
      write(20,*)i
         do j=1,3*natom
            write(20,105)(cmat(i,j,k),k=1,3*natom)
         enddo
      write(20,*)
      enddo
      close(20)

      if(ifile.eq.0)then
         open (unit=99,status='unknown')
         rewind (99)
         write (99,103) nameoutc
         rewind (99)
         read (99,104) command1
         close(99)
         call commrun(command1)         
      endif
 103  format ('cp -f cmat.dat output/'A8,'.dat')
 104  format (A100)
 105  format (10000(1X,1PE11.4))

      return
      end
