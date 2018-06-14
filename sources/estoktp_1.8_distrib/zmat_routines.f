c a set of programa to read and convert z-matrixes
c Copyright (C) 2018  by Carlo Cavallotti and Matteo Pelucchi, Politecnico di MIlano, Italy
c and Stephen J. Klippenstein, Argonne National Laboratories, USA.
c This program is free software: you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation, version 3 of the License
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

      subroutine read_zmat(atomlabel,natom,natomt,intcoor,bislab,ibconn,
     $ iaconn,idconn,bname,anname,dname,atname,idummy,isited,jsited,
     $ ksited,bconnt,aconnt,dconnt)

      implicit double precision (a-h,o-z)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
 
c      integer ibconn
      dimension ibconn(natommx),iaconn(natommx),idconn(natommx)
      dimension idummy(natommx)
      character*60 atomlabel(natommx)
      character*30 intcoor(3*natommx)
      character*20 bislab(ntaumx)
      character*20 bname(natommx),anname(natommx),dname(natommx)
     $ ,atname(natommx),bconnt(natommx),aconnt(natommx),dconnt(natommx)
      character*1 xread
      character*5 ccheck
      LOGICAL leof,lsec,ltit

      CHARACTER*160 line,sename,string,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7
      character*20 ct1,ct2,ct3
      include 'filcomm.f'

      write(*,*)'natomt is ', natomt

      do j=1,natomt
         atname(j)=''
         bconnt(j)=''
         bname(j)=''
         aconnt(j)=''
         anname(j)=''
         dconnt(j)=''
         dname(j)=''
         ibconn(j)=0
         iaconn(j)=0
         idconn(j)=0
      enddo
c      write(*,*)'ok 1'
c      stop
c
      do j=1,natomt
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) atomlabel(j)
         rewind (99)
         call LineRead3 (99)
         close(99)
         atname(j)=word
         bconnt(j)=word2
         bname(j)=word3
         aconnt(j)=word4
         anname(j)=word5
         dconnt(j)=word6
         dname(j)=word7
      enddo

      do k=1,natomt
c         do j=1,natommx
         do j=1,19
            if(j.lt.10)write(ccheck,"(I1)")j
            if(j.lt.100.and.j.gt.9)write(ccheck,"(I2)")j
            if(j.lt.1000.and.j.gt.99)write(ccheck,"(I3)")j
c            write(*,*)'ccheck is ',ccheck
c            write(*,*)'bconnt is ',bconnt(k)
            if(bconnt(k).eq.ccheck)then
               bconnt(k)=atname(j)
            endif
            if(aconnt(k).eq.ccheck)then
               aconnt(k)=atname(j)
            endif
            if(dconnt(k).eq.ccheck)then
               dconnt(k)=atname(j)
            endif
         enddo
c         stop
      enddo
cc now determine connectivity in terms of progressive atom numbering

      do j=1,natomt
         do i=1,natomt
            if(bconnt(j).eq.atname(i))then
               ibconn(j)=i
            endif
            if(aconnt(j).eq.atname(i))then
               iaconn(j)=i
            endif
            if(dconnt(j).eq.atname(i))then
               idconn(j)=i
            endif
         enddo
      enddo

      do j=1,natomt
         write(*,*)'at ',j,' conn is: ',ibconn(j),iaconn(j),idconn(j)
      enddo

cc now check for dummy atoms
      ndummy=0
      do j=1,natomt
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) atname(j)
         rewind (99)
         read(99,*)xread
         close(99)
         if(xread.eq.'X') then
            idummy(j)=1
         else
            idummy(j)=0
         endif
c         write(*,*)'idummy is j: ',j, idummy(j)
      enddo
      isited=0
      jsited=0
      ksited=0
      do j=1,isite
         if(idummy(j).eq.0)isited=isited+1
      enddo
      do j=1,jsite
         if(idummy(j).eq.0)jsited=jsited+1
      enddo
      do j=1,ksite
         if(idummy(j).eq.0)ksited=ksited+1
      enddo

      write(*,*)'isite isited ',isite,isited
      write(*,*)'jsite jsited ',jsite,jsited
      write(*,*)'ksite ksited ',ksite,ksited

      write(*,*)'out of read_zmat'

      return
      end

c ******************************************************

      subroutine update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $ ,idconn,bname,anname,dname,atname,coox,cooy,cooz,xint,tauopt
     $ ,ntau,idummy,ilin,aconnt,bconnt,dconnt,atomlabel)

      implicit double precision (a-h,o-z)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
 
c     
      dimension ibconn(natommx),iaconn(natommx),idconn(natommx)
      dimension xint(3*natommx)
      dimension xintt(3*natommx)
      dimension coox(natommx),cooy(natommx),cooz(natommx)
      dimension cooxt(natommx),cooyt(natommx),coozt(natommx)
      dimension tauopt(ntaumx)
      dimension idummy(natommx)
      character*30 intcoor(3*natommx)
      character*20 bislab(ntaumx)
      character*20 cname
      character*20 bname(natommx),anname(natommx),dname(natommx)
     $ ,atname(natommx),bconnt(natommx),aconnt(natommx),dconnt(natommx)
      LOGICAL leof,lsec,ltit
      character*60 atomlabel(natommx)

      CHARACTER*160 line,sename,string,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'


      nint = natom*3-6-ntau
      if (natom.eq.2) nint=1
      
      do j=1,nint
         xintt(j)=xint(j)
         xint(j)=0.
      enddo
 
      do j=1,ntau
         tauopt(j)=0.
      enddo

      open (unit=98,status='unknown')
      do j=1,natom
         write(98,*)atname(j),coox(j),cooy(j),cooz(j)
      enddo
      close(98)
c      stop


c rototraslate the xyz matrix
      call rototrasl(natom,coox,cooy,cooz,ilin)
cc
cc write xyz matrix   

      open (unit=97,status='unknown')
      do j=1,natom
         write(97,*)atname(j),coox(j),cooy(j),cooz(j)
      enddo
      write(97,*)'ilin is',ilin
      close(97)
c      stop

cc first introduce fake coordinates for dummy atoms
cc to be consistent with numering in input z-mat

      ind=0
      do j=1,natomt
         iangname=0
         write(*,*)'idummy is ', j,idummy(j)
         if(idummy(j).eq.0) then
            ind=ind+1
            cooxt(j)=coox(ind)
            cooyt(j)=cooy(ind)
            coozt(j)=cooz(ind)
         else if (j.eq.3.and.idummy(j).eq.2) then
cc it is assumed a dummy atom is on the z axis with y=0 perp to atom2 or 1 
cc this has been deactivated
            open (unit=99,status='unknown')
            rewind (99)
            write (99,*) bname(j)
            write (99,*) anname(j)
            rewind (99)
            read(99,*)bd
            read(99,*)ang
            close(99)
            if(ibconn(j).eq.1)then
               cooxt(j)=0.
               cooyt(j)=0.
               coozt(j)=bd
            else
               cooxt(j)=cooxt(2)
               cooyt(j)=0.
               coozt(j)=bd
            endif
c            stop
         else if (j.eq.3.and.idummy(j).eq.1) then
cc it is assumed a dummy atom is on the z axis with y=0 perp to atom2 or 1 
            open (unit=99,status='unknown')
            rewind (99)
            write (99,*) bname(j)
            write (99,*) anname(j)
            rewind (99)
            read(99,*)bd
            read(99,*)ang
            close(99)
c            if(iabs.eq.1.and.ireact.gt.3) write (99,*) dname(j)
            iangindex=0
            do ik=4,natomt
               if(idummy(ik).ne.1)then
                  iangname2=0
                  do ij=1,nint
                     if(intcoor(ij).eq.dname(ik))then
                        iangname2=1
                     endif
                  enddo
                  if (iangname2.eq.0)then
                     iangindex=ik
                  endif
               endif
            enddo
c            write(*,*)'iangindex is',iangindex
c            stop
c
            if (iangindex.ne.0) then
               adist2=bd**2
               iat1=ibconn(j)
               iat2=iaconn(j)
               da=bd
               db=sqrt((coox(iat1)-coox(iat2))**2+
     $                    (cooy(iat1)-cooy(iat2))**2+
     $                    (cooz(iat1)-cooz(iat2))**2)
               alfa=ang
               dc=0.
               call distang(da,db,dc,alfa)
c               write(*,*)'dc is',dc
c               stop
               bdist2=dc**2
               iat3=iangindex
               itestind=0
               if(bconnt(iat3).eq.atname(j))itestind=1
c               write(*,*)'itestind is ',itestind
c               write(*,*)'atname is ',atname(j)
c               write(*,*)'bconn is ',bconnt(iat3)
               if(aconnt(iat3).eq.atname(j))itestind=1
c               write(*,*)'itestind is ',itestind
c               write(*,*)'anname is ',aconnt(iat3)
               if(itestind.eq.1)then
                  write(*,*)'found reference atom for dummy atom3'
               endif

               iatemp1=iangindex
               icorr=0
               do jk=1,iatemp1-1
                  if(idummy(jk).eq.1)icorr=icorr+1
               enddo
               iatemp1=iatemp1-icorr
               iai=iatemp1
               iatemp2=ibconn(iangindex)
               icorr=0
               do jk=1,iatemp2-1
                  if(idummy(jk).eq.1)icorr=icorr+1
               enddo
               iatemp2=iatemp2-icorr
               iab=iatemp2
               iatemp3=idconn(iangindex)
               icorr=0
               do jk=1,iatemp3-1
                  if(idummy(jk).eq.1)icorr=icorr+1
               enddo
               iatemp3=iatemp3-icorr
               iac=iatemp3

               write(*,*)'ibconn',iatemp2
               da=sqrt((coox(iai)-coox(iab))**2+
     $                    (cooy(iai)-cooy(iab))**2+
     $                    (cooz(iai)-cooz(iab))**2)

               db=sqrt((coox(iac)-coox(iatemp2))**2+
     $                    (cooy(iac)-cooy(iatemp2))**2+
     $                    (cooz(iac)-cooz(iatemp2))**2)
               write(*,*)'da is',da
               write(*,*)'db is',db

c               db=bd

               open (unit=99,status='unknown')
               rewind (99)
               write (99,*) dname(iangindex)
               rewind (99)
               read(99,*)dih_ia
               rewind (99)
               write (99,*) anname(j)
               rewind (99)
               read(99,*) ang_abc
               close(99)

               write(*,*)'dihed is',dih_ia
               write(*,*)'ang_abc is',ang_abc
               if(dih_ia.eq.180.)then
                  write(*,*)'recognized opposite planar configuration'
               else
                  write(*,*)'failed z-mat conv. in routine upfate_zmat' 
               endif
               dac=sqrt((coox(iai)-coox(iac))**2+
     $                    (cooy(iai)-cooy(iac))**2+
     $                    (cooz(iai)-cooz(iac))**2)

c               dac=dac-0.01
c               write(*,*)'dac is',dac
c               write(*,*)'da+db is',da+db
c               if(abs(dac-da-db).gt.0.5)then
c                  alfa=180.-ang_abc
c               else
               alfa2=0.
               call distang2(da,db,dac,alfa2)
               alfa=alfa2-ang_abc
c               endif
               
c               stop
               db=bd
c              write(*,*)'alfa is ',alfa
c              stop
               call distang(da,db,dc,alfa)
c               write(*,*)'dc ',dc
c               stop
c               adist2=da**2
c               bdist2=db**2
               cdist2=dc**2
               icorr=0
               do jk=1,iat3-1
                  if(idummy(jk).eq.1)icorr=icorr+1
               enddo
               iat3=iat3-icorr
c               write(*,*)'adist2 ', adist2
c               write(*,*)'bdist2 ', bdist2
c               write(*,*)'cdist2 ', cdist2
c               write(*,*)'iat1 ',iat1 
c               write(*,*)'iat2 ',iat2 
c               write(*,*)'iat3 ',iat3 
c               write(*,*)'cooxa1 ',coox(iat1)
c               write(*,*)'cooya1 ',cooy(iat1)
c               write(*,*)'cooza1 ',cooz(iat1)
c               write(*,*)'cooxa2 ',coox(iat2)
c               write(*,*)'cooya2 ',cooy(iat2)
c               write(*,*)'cooza2 ',cooz(iat2)
c               write(*,*)'cooxa3 ',coox(iat3)
c               write(*,*)'cooya3 ',cooy(iat3)
c               write(*,*)'cooza3 ',cooz(iat3)
c               check if colinear on x, y or z axis
c               ichangex=0
c               ichangey=0
c               ichangez=0
c               if(abs(coox(iat1)-coox(iat2)).lt.0.2.and.
c     +            abs(coox(iat2)-coox(iat3)).lt.0.2.and.
c     +            abs(coox(iat1)-coox(iat3)).lt.0.2)then
c                  ichangex=1
c                  coox(iat3)= coox(iat3)+0.01
c               endif
c               if(abs(cooy(iat1)-cooy(iat2)).lt.0.2.and.
c     +            abs(cooy(iat2)-cooy(iat3)).lt.0.2.and.
c     +            abs(cooy(iat1)-cooy(iat3)).lt.0.2)then
c                  ichangey=1
c                  cooy(iat3)= cooy(iat3)+0.01
c               endif
c               if(abs(cooz(iat1)-cooz(iat2)).lt.0.2.and.
c     +            abs(cooz(iat2)-cooz(iat3)).lt.0.2.and.
c     +            abs(cooz(iat1)-cooz(iat3)).lt.0.2)then
c                  ichangez=1
c                  cooz(iat3)= cooz(iat3)+0.01
c               endif

               call twoan_to_xyz(adist2,bdist2,cdist2,coox(iat1)
     $  ,cooy(iat1),cooz(iat1),coox(iat2),cooy(iat2),cooz(iat2),
     $  coox(iat3),cooy(iat3),cooz(iat3),xa,ya,za)

c              if(ichangex.eq.1)then
c                 coox(iat1)= coox(iat1)-0.01
c              endif
c              if(ichangex.eq.1)then
c                  cooy(iat1)= cooy(iat1)-0.01
c               endif
c               if(ichangex.eq.1)then
c                  cooz(iat1)= cooz(iat1)-0.01
c               endif

c               write(*,*)'xa is ',xa
c               write(*,*)'ya is ',ya
c              write(*,*)'za is ',za

               cooxt(j)=xa
               cooyt(j)=ya
               coozt(j)=za
c               stop
            else
cc it is assumed a dummy atom is on the z axis with y=0 perp to atom2 or 1 
               open (unit=99,status='unknown')
               rewind (99)
               write (99,*) bname(j)
               write (99,*) anname(j)
               rewind (99)
               read(99,*)bd
               read(99,*)ang
               close(99)
               if(ibconn(j).eq.1)then
                  cooxt(j)=0.
                  cooyt(j)=0.
                  coozt(j)=bd
               else
                  cooxt(j)=cooxt(2)
                  cooyt(j)=0.
                  coozt(j)=bd
               endif
c               write(*,*) 'working on that'
            endif
         else
c            write(*,*)'j is ', j
            open (unit=99,status='unknown')
            rewind (99)
            write (99,*) bname(j)
            write (99,*) anname(j)
c            if(iabs.eq.1.and.ireact.gt.3) write (99,*) dname(j)
            if(ireact.gt.3) then
               do ij=1,nint
                  if(intcoor(ij).eq.dname(j))then
                     iangname=1
c                     write(99,*)'iname is ',iname
                     write (99,*) xintt(ij)
c                     stop
                 endif
               enddo
               if(iangname.eq.0)then
                  write (99,*) dname(j)
                  iangname=0
               endif
            endif
            rewind (99)
            read(99,*)bd
            read(99,*)ang
c            if(iabs.eq.1.and.ireact.gt.3)read(99,*)dihed
            if(ireact.gt.3)read(99,*)dihed
            close(99)
c            write(bd,"(f8.4)")bname(j)
c            write(ang,"(f8.4)")anname(j)
c            write(dihed,"(f8.4)")dname(j)
            
c            write(*,*)'iaconnj is ', iaconn(j)
c            dihed=dihed+0.1
c            ang=ang+0.1
c            write(*,*)'bd is ',bd
c            write(*,*)'ang is ', ang
c            write(*,*)'dihed is ', dihed

c            call zmat_to_xyz(xa,ya,za,coox(ibconn(j)),cooy(ibconn(j)),
c     $ cooz(ibconn(j)),coox(iaconn(j)),cooy(iaconn(j)),cooz(iaconn(j)),
c     $ coox(idconn(j)),cooy(idconn(j)),cooz(idconn(j)),
c     $ bd,ang,dihed)
            if(iangname.eq.0)then
            call zmat_to_xyz(xa,ya,za,cooxt(ibconn(j)),cooyt(ibconn(j)),
     $ coozt(ibconn(j)),cooxt(iaconn(j)),cooyt(iaconn(j)),
     $ coozt(iaconn(j)),
     $ cooxt(idconn(j)),cooyt(idconn(j)),coozt(idconn(j)),
     $ bd,ang,dihed)
            else if (iangname.eq.1) then
               dist1=bd
               iat2=ibconn(j)
               iat1=iaconn(j)
               ang1=ang
c               write(*,*) 'working on that'
               write(*,*)'dist1 ', dist1
               write(*,*)'iat1 ', iat1
               icorr=0
               do jk=1,iat1-1
                  if(idummy(jk).eq.1)icorr=icorr+1
               enddo
               iat1=iat1-icorr
               write(*,*)'iat1 ', iat1
               icorr=0
               do jk=1,iat2-1
                  if(idummy(jk).eq.1)icorr=icorr+1
               enddo
               iat2=iat2-icorr
               write(*,*)'iat2 ', iat2
               write(*,*)'ang1 is ',ang
               write(*,*)'j is ',j
               if(ibconn(j+1).eq.ibconn(j).and.iaconn(j+1).eq.j)then
                  write(*,*)'recognized 2 successive ang dihedral def'
                  iat3=j+1
                  write(*,*)'iat3 ', iat3
                  icorr=0
                  do jk=1,iat3-1
                     if(idummy(jk).eq.1)icorr=icorr+1
                  enddo
                  iat3=iat3-icorr
                  write(*,*)'iat3 ', iat3
                  do ij=1,nint
                     if(intcoor(ij).eq.anname(j+1))then
                        write (*,*)'failed to recog dummy 2 ang conf'
                        write (*,*)'stopping now'
                        stop
                     endif
                  enddo
                  open (unit=99,status='unknown')
                  write (99,*) anname(j+1)
                  rewind (99) 
                  read(99,*) ang2
                  close(99)
c                  write(*,*)'ang2 is ',ang2
               else
                  write (*,*)'failed2 to recog dummy 2 ang conf'
                  write (*,*)'stopping now'
                  stop
               endif

c               write(*,*)'iat1 ', iat1
c               write(*,*)'iat2 ', iat2
c               write(*,*)'iat3 ', iat3

               dist2=sqrt((coox(iat1)-coox(iat2))**2+
     $                    (cooy(iat1)-cooy(iat2))**2+
     $                    (cooz(iat1)-cooz(iat2))**2)
               dist3=sqrt((coox(iat2)-coox(iat3))**2+
     $                    (cooy(iat2)-cooy(iat3))**2+
     $                    (cooz(iat2)-cooz(iat3))**2)
c               write(*,*)'dist2 ', dist2
c               write(*,*)'dist3 ', dist3
               adist2=dist2**2+(dist1+dist2*
     $                cos(3.1415/180*(180-ang1)))**2
               bdist2=dist1*dist1
               cdist2=dist3**2+(dist1+dist3*
     $                cos(3.1415/180*(180-ang1)))**2
c               write(*,*)'adist2 ', adist2
c               write(*,*)'bdist2 ', bdist2
c               write(*,*)'cdist2 ', cdist2
               call twoan_to_xyz(adist2,bdist2,cdist2,coox(iat1)
     $  ,cooy(iat1),cooz(iat1),coox(iat2),cooy(iat2),cooz(iat2),
     $  coox(iat3),cooy(iat3),cooz(iat3),xa,ya,za)
c               stop
            endif

            cooxt(j)=xa
            cooyt(j)=ya
            coozt(j)=za
c            write(*,*)'xa is ,',xa
c            write(*,*)'ya is ,',ya
c            write(*,*)'za is ,',za
c            write(*,*)'x1 is ,',cooxt(ibconn(j))
c            write(*,*)'y1 is ,',cooyt(ibconn(j))
c            write(*,*)'z1 is ,',coozt(ibconn(j))
c            write(*,*)'x2 is ,',cooxt(iaconn(j))
c            write(*,*)'y2 is ,',cooyt(iaconn(j))
c            write(*,*)'z2 is ,',coozt(iaconn(j))
c            write(*,*)'x3 is ,',cooxt(idconn(j))
c            write(*,*)'y3 is ,',cooyt(idconn(j))
c            write(*,*)'z3 is ,',coozt(idconn(j))
c            stop

         endif
c        write(*,*)'cooxj is ,',cooxt(j)
c        write(*,*)'cooyj is ,',cooyt(j)
c        write(*,*)'coozj is ,',coozt(j)
      enddo

c      open (unit=99,status='unknown')
      do j=1,natomt
         write(*,*)atname(j),cooxt(j),cooyt(j),coozt(j)
      enddo
c      close(99)
c      stop



cc  update  nint intcoor variables

      do j=1,nint
         ind=0
         do i=1,natomt
            if(intcoor(j).eq.bname(i)) then
               xint(j)=sqrt((cooxt(i)-cooxt(ibconn(i)))**2
     $        +(cooyt(i)-cooyt(ibconn(i)))**2
     $        +(coozt(i)-coozt(ibconn(i)))**2)      
            endif
            if(intcoor(j).eq.anname(i)) then
               open (unit=99,file='dihed.dat',status='unknown')
               write(99,*)cooxt(i),cooyt(i),coozt(i)
               ianum=ibconn(i)
               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
               ianum=iaconn(i)
               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
               ianum=idconn(i)
               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
               close(99)
               call dihedral
               open(unit=23,file='dihed.res',status='unknown')
               read(23,*)xint(j)
               close(23)
            endif
            if(intcoor(j).eq.dname(i)) then
               open (unit=99,file='dihed.dat',status='unknown')
               write(99,*)cooxt(i),cooyt(i),coozt(i)
               ianum=ibconn(i)
               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
               ianum=iaconn(i)
               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
               ianum=idconn(i)
               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
               close(99)
               call dihedral
               open(unit=23,file='dihed.res',status='unknown')
               read(23,*)cjunk
               read(23,*)xint(j)
               close(23)
            endif
         enddo
         write(*,*)'xint is',j,xint(j)
      enddo

      do j=1,ntau
         nprog=natom*3-6-ntau+j
         do i=1,natomt
            if(bislab(j).eq.bname(i)) then
               tauopt(j)=sqrt((cooxt(i)-cooxt(ibconn(i)))**2+
     $       (cooyt(i)-cooyt(ibconn(i)))**2+(coozt(i)-coozt(ibconn(i))))      
            endif
            if(intcoor(j).eq.anname(i)) then
               open (unit=99,file='dihed.dat',status='unknown')
               write(99,*)cooxt(i),cooyt(i),coozt(i)
               ianum=ibconn(i)
               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
               ianum=iaconn(i)
               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
               ianum=idconn(i)
               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
               close(99)
               call dihedral
               open(unit=23,file='dihed.res',status='unknown')
               read(23,*)tauopt(j)
               close(23)
            endif
            if(intcoor(j).eq.dname(i)) then
               open (unit=99,file='dihed.dat',status='unknown')
               write(99,*)cooxt(i),cooyt(i),coozt(i)
               ianum=ibconn(i)
               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
               ianum=iaconn(i)
               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
               ianum=idconn(i)
               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
               close(99)
               call dihedral
               open(unit=23,file='dihed.res',status='unknown')
               read(23,*)cjunk
               read(23,*)tauopt(j)
               close(23)
            endif
         enddo
         if (tauopt(itau).gt.360.0d0) tauopt(itau) = 
     $        tauopt(itau) - 360.0d0
         if (tauopt(itau).lt.0.0d0) tauopt(itau) = 
     $        tauopt(itau) + 360.0d0
         write(*,*)'tauopt is',j,tauopt(j)
      enddo
cc now update atomlabel if iangindex not equal to zero (dummy atom at position
cc 3)
      
c      if(iangindex.ne.0)then
c      write(*,*)'iangindex is',iangindex
c      write(*,*)'atomlabel ',atomlabel(iangindex)
c      write(*,*)'atname(iangindex) ',atname(iangindex)
c      write(*,*)'bconnt ',bconnt(iangindex)
c      write(*,*)'bname ',bname(iangindex)
c      write(*,*)'aconnt ',aconnt(iangindex)
c      write(*,*)'aname ',anname(iangindex)
c      write(*,*)'dconnt ',dconnt(iangindex)
c      write(*,*)'dname ',dname(iangindex)
c         ia=iangindex
c         open (unit=99,file='dihed.dat',status='unknown')
c         write(99,*)cooxt(ia),cooyt(ia),coozt(ia)
c         ianum=ibconn(ia)
c         write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c         ianum=iaconn(ia)
c         write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c         ianum=idconn(ia)
c         write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c         close(99)
c         call dihedral
c         open(unit=23,file='dihed.res',status='unknown')
c         read(23,*)cjunk
c         read(23,*)cname
c         close(23)
c         write(*,*)'cname is ',cname
c      stop
c
c         open (unit=99,status='unknown')
c         write (99,2604) atname(ia),bconnt(ia),bname(ia),
c     $        aconnt(ia),anname(ia),dconnt(ia),cname
c         rewind(99)
c         read (99,'(A60)') atomlabel(ia)         
c         close(99)
c         write(*,*)'atomlabel ',atomlabel(iangindex)
c      endif

c 2604 format (1x,3a6,1x,a6,1x,a6,1x,a6,1x,a6)

c      stop

      return
      end

c************************************************************
      
      subroutine zmat_to_xyz(xa,ya,za,xb,yb,zb,xc,yc,zc,
     $ xd,yd,zd,bd,ang,dihed)

cc returns coordinates of a given coords of b,c,d and ab distance, abc angle 
cc and abcd dihedral angle

      implicit double precision (a-h,o-z)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
c     
      LOGICAL leof,lsec,ltit

      CHARACTER*160 line,sename,string,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'

c this is a geometric approach for calculating the xyz coordinates of atom A
c when the xyz coordinates of the B C and D are known and A position is defined with respect 
c to B C and D in Z-matrix notation.
c the adopted approach consists in translating to B the frame and then rotating 
c so that B C and D are on the xy plane, with C on the y axis
c A coordinates are immediately deterined in this frame, which is then
c rototranlated back to the BCD original frame of reference   
c _rt variables are in the rototranslated reference system
c _t variables are in the translated ref system 
c this is similar to the NERF algorithm (see for example 
c Parsons et al. J.Comp.Chem. 26(10) 1063, 2005.)

c first determines coordinates of ABCD in the RT frame of reference
      
      ang=ang*3.14159/180.
      dihed=dihed*3.14159/180.

      bdist=bd*cos(ang)
      ddist=bd*sin(ang)

      xb_rt=0.
      yb_rt=0.
      zb_tr=0.

      xa_rt=ddist*cos(dihed)
      ya_rt=bdist
      za_rt=ddist*sin(dihed)


      bc_dist=sqrt((xc-xb)**2+(yc-yb)**2+(zc-zb)**2)
      bd_dist=sqrt((xb-xd)**2+(yb-yd)**2+(zb-zd)**2)
      cd_dist=sqrt((xc-xd)**2+(yc-yd)**2+(zc-zd)**2)

      xc_rt=0.
      yc_rt=bc_dist
      zc_rt=0.

      yd_rt=(bc_dist**2+bd_dist**2-cd_dist**2)/2./bc_dist
      xd_rt=sqrt(bd_dist**2-yd_rt**2)
      zd_rt=0.

cc translate original frame of ref coords so that B is at (0,0,0)

      xb_t=0.
      yb_t=0.
      zb_t=0.
      xc_t=xc-xb
      yc_t=yc-yb
      zc_t=zc-zb
      xd_t=xd-xb
      yd_t=yd-yb
      zd_t=zd-zb

c      if(yd_t.eq.0)yd_t=0.01
c      if(xd_t.eq.0)xd_t=0.015
cc now determine the rotation matrix to rotate back to the original ref system
      write(*,*)'yc_rt ',yc_rt
      write(*,*)'xd_rt ',xd_rt

      r12=(xc-xb)/yc_rt
      r22=(yc-yb)/yc_rt
      r32=(zc-zb)/yc_rt

      r11=(xd-xb-yd_rt*r12)/xd_rt
      r21=(yd-yb-yd_rt*r22)/xd_rt
      r31=(zd-zb-yd_rt*r32)/xd_rt

      aconst=(yc_t-yd_t/xd_t*xc_t)/(zc_t-zd_t/xd_t*xc_t)
      den1=(yd_t/xd_t-aconst*zd_t/xd_t)
      if(den1.lt.1.0e-6)den1=1.0e-6

c      bconst=1./(yd_t/xd_t-aconst*zd_t/xd_t)
      bconst=1./den1
c      bconst=999999.
      xe_t=-1/sqrt(1+(bconst**2)*(1+aconst**2))
      ye_t=-xe_t*bconst
      ze_t=-ye_t*aconst

c      write(*,*)'den ',(yd_t/xd_t-aconst*zd_t/xd_t)
c      write(*,*)'den1 ',(yd_t/xd_t)
c      write(*,*)'den2 ',(aconst*zd_t/xd_t)
c      write(*,*)'den21 ',(aconst)
c      write(*,*)'den22 ',(zd_t/xd_t)

c      write(*,*)'aconst ',aconst
c      write(*,*)'bconst ',bconst
c      write(*,*)'xd ',xd_t
c      write(*,*)'yd ',yd_t
c      write(*,*)'zd ',zd_t

      r13=xe_t
      r23=ye_t
      r33=ze_t

      xe_rt=0.
      ye_rt=0.
      ze_rt=1.

cc now rotate and translate back

      xa=xb+r11*xa_rt+r12*ya_rt+r13*za_rt
      ya=yb+r21*xa_rt+r22*ya_rt+r23*za_rt
      za=zb+r31*xa_rt+r32*ya_rt+r33*za_rt

      return
      end

c************************************************************
      subroutine rototrasl(natom,coox,cooy,cooz,ilin)

      implicit double precision (a-h,o-z)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
 
c     
      dimension coox(natommx),cooy(natommx),cooz(natommx)
      dimension cooxt(natommx),cooyt(natommx),coozt(natommx)
      dimension cooxt1(natommx),cooyt1(natommx),coozt1(natommx)

      LOGICAL leof,lsec,ltit

      CHARACTER*160 line,sename,string,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'


cc performs traslation to center of axis
cc and zy rotation to get second atom on x axis
      if(ilin.eq.0)then
         do j=1,natom
            cooxt(j)=coox(j)-coox(1)
            cooyt(j)=cooy(j)-cooy(1)
            coozt(j)=cooz(j)-cooz(1)
         enddo

c         do j=1,natom
c            write(*,*)'coox ',cooxt(j)
c            write(*,*)'cooy ',cooyt(j)
c            write(*,*)'cooz ',coozt(j)
c         enddo
c         stop

cc first rotate around z
         ialfajump=0
         dist=sqrt(cooxt(2)*cooxt(2)+cooyt(2)*cooyt(2))
         if(dist.gt.1.0e-2) then
            alfa=acos((cooxt(2))/dist)
         else
            alfa=0.
            ialfajump=1
         endif
c         if(cooyt(2).gt.0)then
c            if(alfa.lt.0)alfa=3.1415-alfa
c         endif
c         if(cooyt(2).lt.0)then
c            if(alfa.gt.0)alfa=3.1415*2-alfa
c            if(alfa.lt.0)alfa=3.1415-alfa
c         endif

         write(*,*)'dist is = ',dist
         write(*,*)'alfa is = ',alfa

         if(cooyt(2).gt.0)then
            alfa=3.14159*2-alfa
c            if(cooxt(2).gt.0)alfa=3.1415-alfa+3.1415
c            if(cooxt(2).lt.0)alfa=3.1415-alfa+3.1415
c            if(alfa.lt.0)alfa=3.1415-abs(alfa)
         endif

c      sinpar=(cooyt(2))/dist

c         write(*,*)'alfa z = ',alfa
c      write(*,*)'alfa y = ',asin(sinpar)
         if(ialfajump.eq.0)then
            do j=1,natom
               cooxt1(j)=cooxt(j)*cos(alfa)-cooyt(j)*sin(alfa)
               cooyt1(j)=cooxt(j)*sin(alfa)+cooyt(j)*cos(alfa)
               coozt1(j)=coozt(j)
            enddo
         else
            do j=1,natom
               cooxt1(j)=cooxt(j)
               cooyt1(j)=cooyt(j)
               coozt1(j)=coozt(j)
            enddo
         endif
         do j=1,natom
            write(*,*)'coox1 ',cooxt1(j)
            write(*,*)'cooy1 ',cooyt1(j)
            write(*,*)'cooz1 ',coozt1(j)
         enddo


cc now rotate around y
         ialfajump=0
         dist=sqrt(cooxt1(2)*cooxt1(2)+coozt1(2)*coozt1(2))
         if(dist.gt.1.0e-2) then
            alfa=acos((cooxt1(2))/dist)
         else
            ialfajump=1
            alfa=0
         endif
      write(*,*)'alfa y = ',alfa
      write(*,*)'ialfa y = ',ialfajump
c      do j=1,natom
c         write(*,*)cooxt1(j),cooyt1(j),coozt1(j)
c      enddo

         if(coozt1(2).lt.0)then
            alfa=3.14159*2-alfa
c            if(cooxt(2).gt.0)alfa=3.1415-alfa+3.1415
c            if(cooxt(2).lt.0)alfa=3.1415-alfa+3.1415
c            if(alfa.lt.0)alfa=3.1415-abs(alfa)
         endif

         if(alfajump.eq.0)then
            do j=1,natom
               coox(j)=cooxt1(j)*cos(alfa)+coozt1(j)*sin(alfa)
               cooy(j)=cooyt1(j)
               cooz(j)=-cooxt1(j)*sin(alfa)+coozt1(j)*cos(alfa)
            enddo
         else
            do j=1,natom
               coox(j)=cooxt1(j)
               cooy(j)=cooyt1(j)
               cooz(j)=coozt1(j)
            enddo
         endif
         do j=1,natom
            write(*,*)'coox2 ',coox(j)
            write(*,*)'cooy2 ',cooy(j)
            write(*,*)'cooz2 ',cooz(j)
         enddo


c      do j=1,natom
c         write(*,*)coox(j),cooy(j),cooz(j)
c      enddo
c      write(*,*)'ok 2'
      else
cc check molecule axis: x (1), y(2) or z(3)
c         write(*,*)'arrived here'
c         stop

         cooxtot=0.
         cooytot=0.
         cooztot=0.
         do j=1,natom
            cooxtot=cooxtot+abs(coox(j))
            cooytot=cooytot+abs(cooy(j))
            cooztot=cooztot+abs(cooz(j))
         enddo
         idir=0
         if(cooxtot.gt.cooytot.and.cooxtot.gt.cooztot)idir=1
         if(cooytot.gt.cooxtot.and.cooytot.gt.cooztot)idir=2
         if(cooztot.gt.cooxtot.and.cooztot.gt.cooytot)idir=3
         if(idir.eq.2)then
            do j=1,natom
               coox(j)=cooy(j)
               cooy(j)=0.
            enddo
            do j=1,natom
               coox(j)=coox(j)-coox(1)
            enddo
         endif
         if(idir.eq.3)then
            do j=1,natom
               coox(j)=cooz(j)
               cooz(j)=0.
            enddo
            do j=1,natom
               coox(j)=coox(j)-coox(1)
            enddo
         endif
      endif

      return
      end

c************************************************************
         subroutine twoan_to_xyz(da2,db2,dc2,x1,y1,z1,x2,y2,z2,x3,y3,z3
     $  ,xa,ya,za)

cc returns coordinates xa ya za for an atom defined with respect to
cc three atoms whose coordinates are known and from which the distances 
cc are known
cc usee Newton Raphson to solve non linear 3 eq 3 unkown problem
cc using analytic Jacobian

      implicit double precision (a-h,o-z)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
c     
      LOGICAL leof,lsec,ltit

      CHARACTER*160 line,sename,string,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'

cc assume guess at 0.
      
c            if(ibconn(j).eq.1)then
c               cooxt(j)=0.
c               cooyt(j)=0.
c               coozt(j)=bd
c            else
c               cooxt(j)=cooxt(2)
c               cooyt(j)=0.
c               coozt(j)=bd

      write(*,*)'x1 is ',x1
      write(*,*)'y1 is ',y1
      write(*,*)'z1 is ',z1
      write(*,*)'x2 is ',x2
      write(*,*)'y2 is ',y2
      write(*,*)'z2 is ',z2
      write(*,*)'x3 is ',x3
      write(*,*)'y3 is ',y3
      write(*,*)'z3 is ',z3
      write(*,*)'da2 is ',da2
      write(*,*)'db2 is ',db2
      write(*,*)'dc2 is ',dc2

c      stop

      xa=(da2-db2+x2*x2)/2/x2
      Ac=-z3/y3
      Bc=-(dc2-da2+2*xa*x3-x3*x3-y3*y3-z3*z3)/2/y3
      za=(-Ac*Bc+sqrt(abs(Ac*Ac*Bc*Bc-(1+Ac*Ac)*(xa*xa-da2+Bc*Bc))))
     +   /(1+Ac*Ac)
      ya=Ac*za+Bc

      write(*,*)'xa is ',xa
      write(*,*)'ya is ',ya
      write(*,*)'za is ',za
      write(*,*)'sqrt is ',Ac*Ac*Bc*Bc-(1+Ac*Ac)*(xa*xa-da2+Bc*Bc)

c      xa=0.
c      ya=0.
c      za=0.
c      if(abs(xa-x1).lt.0.1)xa=0.5
c      if(abs(ya-y1).lt.0.1)ya=0.5
c      if(abs(za-z1).lt.0.1)za=0.5
c      ya=1.5
c      za=1.5

c      do j=1,100
c
c         d1=-((xa-x1)**2+(ya-y1)**2+(za-z1)**2-da2)/2
c         d2=-((xa-x2)**2+(ya-y2)**2+(za-z2)**2-db2)/2
c         d3=-((xa-x3)**2+(ya-y3)**2+(za-z3)**2-dc2)/2
c         a1=xa-x1
c         a2=xa-x2
c         a3=xa-x3
c         b1=ya-y1
c         b2=ya-y2
c         b3=ya-y3
c         c1=za-z1
c         c2=za-z2
c         c3=za-z3

c         write(*,*)'d1 is',d1
c         write(*,*)'d2 is',d2
c         write(*,*)'d3 is',d3
      
c         gr_i=(d2-d1/a1*a2)/(b2-b1/a1*a2)
c         gr_h=(c1/a1*a2-c2)/(b2-b1/a1*a2)
c         gr_m=-gr_h*b1/a1*a3+b3*gr_h-c1/a1*a3+c3
c         gr_n=d3-d1/a1*a3-gr_i*(-b1/a1*a3+b3)
c         delta3=gr_n/gr_m
c         delta2=delta3*gr_h+gr_i
c         delta1=d1/a1-b1/a1*delta2-c1/a1*delta3

c         write(*,*)'delta1 is',delta1
c         write(*,*)'delta2 is',delta2
c         write(*,*)'delta3 is',delta3
         
c         xa=xa+delta1
c         ya=ya+delta2
c         za=za+delta3
c         error_tot=abs(d1)+abs(d2)+abs(d3)
c         if (error_tot.lt.1.0e-10) goto 100
         
c         write(*,*)'xa is',xa
c         write(*,*)'ya is',ya
c         write(*,*)'za is',za
c         write(*,*)' error 1 is ',d1
c         write(*,*)' error 2 is ',d2
c         write(*,*)' error 2 is ',d3
c      enddo
c 100  continue
c      if (error_tot.gt.1.0) then
c         write(*,*)'failed xyz translation for dummy atoms '
c         write(*,*)'error in zmat-update'
c         stop
c      endif


c      write(*,*)'convergence reached in dummy atom xyz coo search'
c      stop

      return
      end
c************************************************************
      subroutine distang(da,db,dc,alfa)
c
c returns dc given da,db and alfa
c
      implicit double precision (a-h,o-z)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
c     
      LOGICAL leof,lsec,ltit

      CHARACTER*160 line,sename,string,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'
      pi=3.14159265359
      alfa=pi*alfa/180.

      da1=db*cos(pi-alfa)
      db1=db*sin(pi-alfa)

      dc=sqrt((da+da1)**2+db1**2)

      return
      end
c************************************************************
      subroutine distang2(da,db,dc,alfa)
c
c returns alfa given da, db, dc
c
      implicit double precision (a-h,o-z)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
c     
      LOGICAL leof,lsec,ltit

      CHARACTER*160 line,sename,string,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'
      write(*,*)'da is',da
      write(*,*)'db is',db
      write(*,*)'dc is',dc
      pi=3.14159265359
      beta=acos((dc*dc-da*da-db*db)/2/da/db)
c     write(*,*)'acos is ',(dc*dc-da*da-db*db)/2/da/db
c     write(*,*)'beta is ',beta
c      stop
      alfa=(pi-beta)*180./pi

      return
      end
c************************************************************

