c molprofopt a program to start and read output from molpro 2010-2015
c Copyright (C) 2018  by Carlo Cavallotti and Matteo Pelucchi, Politecnico di MIlano, Italy
c and Stephen J. Klippenstein, Argonne National Laboratories, USA.
c This program is free software: you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation, version 3 of the License
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.

      SUBROUTINE molprofopt(tau,ntau,natom,natomt,nshared,gmemo,
     $ coord,vtot,freq,ifreq,ilin,ismp,icharge, ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ilev
     $ ,ispecies)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

c     parameter (ntaumx=10, nmdmx=300, natommx=100, ndim=3)

      integer natom,iatom,idim
      dimension freq(nmdmx),tau(ntaumx),tauopt(ntaumx)
      dimension coord(natommx,ndim),xint(3*natommx),abcrot(ndim)
      dimension coox(natommx),cooy(natommx),cooz(natommx)
      character*160 stringr
      character*70 comline1
      character*100 command1
      character*30 gmemo
      character*30 cjunk
      character*60 atomlabel(natommx)
      character*80 irconslabel
      character*60 label
      character*30 intcoor(3*natommx)
      character*30 cootoread
      character*20 bislab(ntaumx)
      character*2 aname
      character*2 atomtype(natommx)
      
c     character*160 cname

      logical leof,lsec,ltit

      character*1000 line,string
      character*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      character*20 ct1,ct2,ct3

c      COMMON /keyword/ leof,lsec,ltit
c      COMMON /key/ line,sename,string,word,word2,word3
c     $ ,word4,word5,title,title1

c      common /tsdat/ aabs1,babs1,aabs2,babs2,babs3,rts
c      common /itsdat/ iabs,iadd,ivar,isite,ji,ki
c      common /strucword/ stoich
      include 'filcomm.f'

111   continue

cccccccccccccc print out molpro input file  ccccccccccccc
c initialize word, word2, word3, word4, word5
      call LineRead (0)
      ichecken=0
      irepeat=0
c      ncoord = natom*3-6-ntau-ircons
      ncoord = natom*3-6-ntau
      if (natom.eq.2) ncoord = 1
c      write(*,*)'natomt is',natomt
c      stop
      write(*,*)'ircons is ',ircons

 999  continue

      open(unit=10,status='unknown',file='molpro.inp')
      if(ilev.eq.0) then
         open(unit=11,status='unknown',file='level0_molpro.dat')
      else if(ilev.eq.1) then
         open(unit=11,status='unknown',file='level1_molpro.dat')
      else if (ilev.eq.2) then
         open(unit=11,status='unknown',file='hl_molpro.dat')
      else if (ilev.eq.10) then
         open(unit=11,status='unknown',file='onedtau_molpro.dat')
      else
         write(*,*)'molpro call not implemented yet '
         write(*,*)'for this type of calculation'
         stop
      endif

cc assign memory for calculation
cc it is assumed the memory in teh input is given as MW

      igmem=100
      open(unit=99,file='temp.tmp',status='unknown')
      write(99,*)gmemo
      close(99)
      command1="sed -ie 's/MW/ /' temp.tmp"
      call commrun(command1)
      open(unit=99,file='temp.tmp',status='unknown')
      read(99,*)igmem
      close(99)
      write(10,*)'memory,',igmem,',m'

 100  continue
      read (11,'(A70)') comline1
      if (comline1.EQ.'End1') go to 200
      write (10,*) comline1
      goto 100
 200  continue

      write (10,*)'geometry={angstrom'
      if(ixyz.ne.1)then
         do iatom = 1 , natomt
            write (10,*)atomlabel(iatom)
         enddo
      else
         do iatom = 1 , natom
            open(unit=99,status='unknown')
            write(99,*)atomlabel(iatom)
            rewind(99)
            read(99,*)aname,tcoox,tcooy,tcooz
            write(10,*)aname,', ',tcoox,', ',tcooy,', ',tcooz
            close(99)
         enddo
      endif
      write (10,*)'}'
      
      if(ixyz.ne.1)then
         do icoord = 1 , ncoord
            write (10,*) intcoor(icoord),' = ',xint(icoord)
         enddo
         do itau = 1 , ntau
            write (10,*) bislab(itau),' = ',tau(itau)
         enddo
      endif
      write (10,*)
      write (10,*) 'SET,SPIN=',ispin-1

 101  continue

      read (11,'(A70)') comline1

      if (comline1.EQ.'End2') go to 201
      if (ispin.eq.1) write (10,*) comline1
      
      goto 101
 201  continue
      if (ispin.eq.1) go to 401
cc else read input for open shell
 301  continue

      read (11,'(A70)') comline1

      if (comline1.EQ.'End3') go to 401
      write (10,*) comline1
      goto 301
 401  continue
 
      close(11)
      write(10,*)'put,molden,molpro.molden'
      write(10,*)'---'
      close (unit=10,status='keep')

      ntotcoord = natom*3-6
      if(ircons.ne.0)then
         open(unit=99,file='temp.tmp',status='unknown')
         if (ircons.eq.1) then
            write (99,*) 'inactive,',intcoor(ntotcoord)
         else if (ircons.eq.2) then
            write (99,*)'inactive,',intcoor(ntotcoord),
     $           ',',intcoor(ntotcoord-1)
         else if (ircons.eq.3) then
            write (99,*)'inactive,',intcoor(ntotcoord),
     $           ',',intcoor(ntotcoord-1),
     $           ',',intcoor(ntotcoord-2)
         else if (ircons.eq.4) then
            write (99,*)'inactive,',intcoor(ntotcoord),
     $           ',',intcoor(ntotcoord-1),
     $           ',',intcoor(ntotcoord-2),
     $           ',',intcoor(ntotcoord-3)
         else if (ircons.eq.5) then
            write (99,*)'inactive,',intcoor(ntotcoord),
     $           ',',intcoor(ntotcoord-1),
     $           ',',intcoor(ntotcoord-2),
     $           ',',intcoor(ntotcoord-3),
     $           ',',intcoor(ntotcoord-4)
         endif
         close(99)
         command1="sed -ie 's/ //g' temp.tmp"
         call commrun(command1)

         open(unit=99,file='temp.tmp',status='unknown')
         read(99,'(A80)')irconslabel
         rewind(99)
         write(99,1200)irconslabel
         rewind(99)
         read(99,'(A80)')irconslabel
         close(99)
         command1=irconslabel
         call commrun(command1)

 1200    format (" sed -ie 's/inactive/"A35"/g' molpro.inp")
         write(*,*)'irconslab is ',irconslabel
      endif
c
      if(ispecies.eq.0.and.ircons.eq.3.and.ilev.eq.10)then
         open(unit=99,status='unknown')
         write(99,1202)
         rewind(99)
         read(99,'(A60)')label
         close(99)
         command1=label
         call commrun(command1)
 1202    format (" sed -ie 's/root=2/root=1/g' molpro.inp")
         open(unit=99,status='unknown')
         write(99,1203)
         rewind(99)
         read(99,'(A60)')label
         close(99)
         command1=label
         call commrun(command1)
 1203    format (" sed -ie 's/Root=2/root=1/g' molpro.inp")
      endif

c            write (*,*)'intcoord is ',intcoor(1)
c            stop

ccccccccccccc run molpro  cccccccccccccccc
      command1='rm -f molpro.xml'       
      call commrun(command1)
      command1='rm -f molpro.log_*'       
      call commrun(command1)

      if(irecov.ne.1)then
         command1='rm -f molpro.out'       
         call commrun(command1)
         call molprorun(numproc)
      endif

cc first copy output file to geom.log, for saving

      command1='cp -f molpro.out geom.log'       
      call commrun(command1)

cc now read results

      call LineRead (0)

cc first get energy

      command1='egrep CBSEN molpro.out > temp1.log'
      call commrun(command1)
      command1='tail -n 1 temp1.log > temp.log'
      call commrun(command1)

      open (unit=99,status='unknown',file='temp.log')
      read(99,*)cjunk,cjunk,cjunk,vtot
      close(99)
      write(*,*)'vtot is',vtot

cc then get coordinates if not an xyz calculations
      
c      open (unit=100,status='unknown',file='molpro.molden')
      if(ixyz.ne.1)then
         do j=1,ncoord+ntau
c         rewind(100)
            open (unit=99,status='unknown',file='temp.log')
            if(j.le.ncoord) write(99,*) intcoor(j)
            if(j.gt.ncoord) write(99,*) bislab(j-ncoord)
            close(99)
            call trimtemp
            open (unit=99,status='unknown',file='temp1.log')
            call lineread(99)
            cootoread=word
            write(*,*)'cootoread is ',cootoread
            rewind(99)
            write(99,1201)cootoread
 1201       format (" egrep "A20" molpro.molden > temp.log")
            rewind(99)
            read(99,'(A60)')command1
            call commrun(command1)
            close(99)         
            command1="sed -ie 's/=/ /g' temp.log"
            call commrun(command1)
            open (unit=99,status='unknown',file='temp.log')
            read(99,*)cjunk,xint(j)         
            close(99)
c         stop
c         do while (WORD.ne.cootoread)
c            CALL LineRead (100)
c            write(*,*)'word is ',word
c         enddo
c         open (unit=99,status='unknown')
c         write(99,*) word2
c         rewind(99)
c         read(99,*)xint(j)
c         close(99)
         enddo
c      close(100)
         iread=1
         itau=0
         do j=1,ncoord+ntau
            write(*,*)'opt coord is',xint(j)
            itau = iread-(natom*3-6-ntau)
            if (iread.gt.ncoord) then
               tauopt(itau)=xint(j)
               if (tauopt(itau).gt.360.0d0) tauopt(itau) = 
     $              tauopt(itau) - 360.0d0
               if (tauopt(itau).lt.0.0d0) tauopt(itau) = 
     $              tauopt(itau) + 360.0d0
            endif
            iread=iread+1
         enddo
      endif
c      stop
cc then get and write xyz coordinates

      open (unit=100,status='unknown',file='molpro.molden')
      open (unit=101,status='unknown',file='geom.xyz')

      do while (WORD.ne.'[ATOMS]')
         CALL LineRead (100)
      enddo
      do j=1,natom
         read(100,*)atomtype(j),djunk,djunk,coox(j),cooy(j),cooz(j)
         write(101,*)atomtype(j),coox(j),cooy(j),cooz(j)
         coord(j,1)=coox(j)
         coord(j,2)=cooy(j)
         coord(j,3)=cooz(j)
      enddo
      write(101,*)
      close(100)
      close(101)

cc now read frequencies, if computed

cc first do a check in the molden file to see if the frequencies key was there

      command1='grep -q "\[FREQ" molpro.molden && echo 1 > temp.log || 
     $  echo 0 > temp.log'
      call commrun(command1)

      open (unit=99,status='unknown',file='temp.log')
      read(99,*)ifreq
      close(99)

cc then read fcmat.log file, even if no freqs must be written

      open (unit=101,status='unknown',file='fcmat.log')
      write(101,*)'frequencies not found in molden file'
      close(101)

cc then read frequencies

      if(ifreq.eq.1)then
         open (unit=100,status='unknown',file='molpro.molden')

         do while (WORD.ne.'[FREQ]')
            CALL LineRead (100)
         enddo
         index=0.
         nread=3*natom
         do j=1,nread
            read(100,*)freqread
c            write(*,*)'freq ',j,' is ',freqread
            if(freqread.gt.0.5)then
               index=index+1
               freq(index)=freqread
cc now assign negavtive value to imaginary frequency
               if(freq(index).lt.freq(index-1))then
                  freq(index-1)=-freq(index-1)
c                  index=index-1
               endif
            endif
         enddo
         nfreq=3*natom-6
         if(ilin.eq.1)  nfreq=3*natom-5
         if(natom.eq.2)  nfreq=1
c         if(ispecies.eq.0)nfreq=nfreq-1
         if(nfreq.ne.index)then
            write(*,*)' there is disagreement between '
            write(*,*)' expected and read frequencies'
            write(*,*)' the program will be stopped'
            write(*,*)' read frequencies ',index
            write(*,*)' expected frequencies ',nfreq
            stop
         endif
         do j=1,nfreq
            write(*,*)'freq ',j,' is ',freq(j)
         enddo
         close(100)
cc now proceed to saving the hessian
         open (unit=100,status='unknown',file='molpro.out')
         open (unit=101,status='unknown',file='fcmat.log')
         write(101,*)'Force constants in Cartesian coordinates'
         numlines=0
         nummax=natom*3/5
         do j=0,nummax
            numlines=numlines+natom*3-5*j+1
         enddo
c         numlines=numlines-1
         write(*,*)'numlines =',numlines

         do while (stringr.ne.
     $  ' Force Constants (Second Derivatives of the Energy) in [a.u.]')
            read(100,'(A160)')stringr
            if(stringr.eq.' Variable memory released')then
               write(*,*)'did not found Hessian infos in the output'
               goto 500
            endif
         enddo
         do j=1,numlines
            read(100,'(A160)')stringr
            write(101,'(A100)')stringr
         enddo
         write(101,*)
         write(101,*)
 500     continue
         write(101,*)'Input orientation'
         write(101,*)'---'
         write(101,*)'Center'
         write(101,*)'Number'
         write(101,*)'---'

         do j=1,natom
            if(atomtype(j).eq.'H')ianumb=1
            if(atomtype(j).eq.'B')ianumb=5
            if(atomtype(j).eq.'C')ianumb=6
            if(atomtype(j).eq.'N')ianumb=7
            if(atomtype(j).eq.'O')ianumb=8
            if(atomtype(j).eq.'F')ianumb=9
            if(atomtype(j).eq.'Ne')ianumb=10
            if(atomtype(j).eq.'Al')ianumb=13
            if(atomtype(j).eq.'Si')ianumb=14
            if(atomtype(j).eq.'P')ianumb=15
            if(atomtype(j).eq.'S')ianumb=16
            if(atomtype(j).eq.'Cl')ianumb=17
            if(atomtype(j).eq.'Ar')ianumb=18
            if(atomtype(j).eq.'Ga')ianumb=31
            if(atomtype(j).eq.'Ge')ianumb=32
            if(atomtype(j).eq.'As')ianumb=33
            if(atomtype(j).eq.'Se')ianumb=34
            if(atomtype(j).eq.'Br')ianumb=35
            if(atomtype(j).eq.'Kr')ianumb=36
            if(atomtype(j).eq.'In')ianumb=49
            if(atomtype(j).eq.'Sn')ianumb=50
            if(atomtype(j).eq.'Sb')ianumb=51
            if(atomtype(j).eq.'Te')ianumb=52
            if(atomtype(j).eq.'I')ianumb=53
            if(atomtype(j).eq.'Xe')ianumb=54
            write(101,1002)j,ianumb,0,coox(j),cooy(j),cooz(j)
         enddo

         write(101,*)
         write(101,*)
         write(101,*)'Forces '
         write(101,*)'---'
         write(101,*)'Center'

         do j=1,natom
            if(atomtype(j).eq.'H')ianumb=1
            if(atomtype(j).eq.'B')ianumb=5
            if(atomtype(j).eq.'C')ianumb=6
            if(atomtype(j).eq.'N')ianumb=7
            if(atomtype(j).eq.'O')ianumb=8
            if(atomtype(j).eq.'F')ianumb=9
            if(atomtype(j).eq.'Ne')ianumb=10
            if(atomtype(j).eq.'Al')ianumb=13
            if(atomtype(j).eq.'Si')ianumb=14
            if(atomtype(j).eq.'P')ianumb=15
            if(atomtype(j).eq.'S')ianumb=16
            if(atomtype(j).eq.'Cl')ianumb=17
            if(atomtype(j).eq.'Ar')ianumb=18
            if(atomtype(j).eq.'Ga')ianumb=31
            if(atomtype(j).eq.'Ge')ianumb=32
            if(atomtype(j).eq.'As')ianumb=33
            if(atomtype(j).eq.'Se')ianumb=34
            if(atomtype(j).eq.'Br')ianumb=35
            if(atomtype(j).eq.'Kr')ianumb=36
            if(atomtype(j).eq.'In')ianumb=49
            if(atomtype(j).eq.'Sn')ianumb=50
            if(atomtype(j).eq.'Sb')ianumb=51
            if(atomtype(j).eq.'Te')ianumb=52
            if(atomtype(j).eq.'I')ianumb=53
            if(atomtype(j).eq.'Xe')ianumb=54
            write(101,1001)j,ianumb,0.000,0.000,0.000
         enddo
         write(101,*)
         write(101,*)
         close(100)
         close(101)
      endif

1001  format(1X,I2,1X,I2,1X,3(1x,F9.4))
1002  format(1X,I2,1X,I2,1X,I2,1X,3(1x,F9.4))

      write(6,*)'out of molprofopt' 

      RETURN
      END 




