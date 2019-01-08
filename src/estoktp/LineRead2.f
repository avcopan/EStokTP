

      SUBROUTINE LineRead2(IUnit)

C     Subroutine reads the next non-blank non-comment line from IUNIT.
C     WORD contains the first keyword (uppercase).
C     WORD2 contains the second keyword, if any.
C     WORD3 contains the third keyword, if any.
C     WORD4 contains the fourth keyword, if any.
C     WORD5 contains the fifth keyword, if any.
C     NWORD is 1 if one keyword, 2 if two keywords, 3 if three keywords.
C     If the first non-blank character is *, LSEC is set to
C     .TRUE. and WORD is the section name.
C     If the first non-blank character is &, the line is printed. 
C     Comments begin with ! and may begin anywhere in a line.
C     Comments following the keywords are detected and ignored.
C     If end of file is detected LEOF is set to .TRUE.

      implicit double precision(a-h,o-z)
      implicit integer (i-n)

c     keyword common block

      logical leof,lsec,ltit
      character*160 line,sename,string,word,word2,word3,title,title1
     #,word4,word5

      include 'filcomm.f'

 5    Continue
      LSEC = .FALSE.
      LEOF = .FALSE.
      LTIT=.FALSE.

 10   READ(iunit,'(A160)',END=9000) STRING
c     write (6,*) 'string2 test',string

      IP = 0
      NWORD=1

C     Find the first nonblank character in the line

 20   IP = IP + 1
      IF (((STRING(IP:IP).EQ.' ').OR.(STRING(IP:IP).EQ.'	'))
     &     .AND.(IP .Lt. 80)) GO TO 20           

C     If it  is all blank or is a comment line then read the next line

      IF (IP .ge. 80)  GO TO 10

C     If it is a title line write it and read next line.

      IF (STRING(IP:IP) .EQ. '&') THEN
         IPP=IP+1 
         TITLE1=STRING(IPP:)
         Write (6,27) TITLE1
 27      Format (a)
         Go to 5
      ENDIF

C     Convert to uppercase

      call upcase

C     Find the position of a comment character, if any,
C     in a nonblank noncomment line

      ICOM = IP 
 50   ICOM = ICOM + 1
      IF (ICOM .LT. 80)  GOTO 50      

 55   IF (IP .LT. ICOM) THEN
         IBEG = IP
 60      IP = IP + 1
         IF (LINE(IP:IP) .NE. ' ' .AND. IP .LT. ICOM) GOTO 60 
         WORD = LINE(IBEG:IP-1)
      ENDIF

c     Find second keyword

 70   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 70
      NWORD=2
      IBEG2=IP

 80   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 80 
      WORD2 = LINE(IBEG2:IP-1)

C     Find third keyword

 75   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 75
      NWORD=3
      IBEG3=IP

 85   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 85 
      WORD3 = LINE(IBEG3:IP-1)

C     Find fourth keyword

175   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 175
      NWORD=4
      IBEG4=IP

185   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 185 
      WORD4 = LINE(IBEG4:IP-1)
      IF (WORD4 .EQ. '=') GOTO 175

C     Find fifth keyword

275   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 275
      NWORD=5
      IBEG5=IP

285   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 285 
      WORD5 = LINE(IBEG5:IP-1)
      IF (WORD5 .EQ. '=') GOTO 275

C     See if WORD is a section header

 90   IF (WORD(1:1) .EQ. '*') THEN
         LSEC=.TRUE.
         SENAME = WORD(2:)
         WORD = SENAME
         If (WORD.eq.'END') Then
            LEOF = .TRUE.
            Return
         Endif
         LSEC = .TRUE.
      ENDIF

      RETURN

 9000 Continue
      Write (6,*) 'Unexpected end of input file upon calling LineRead'
      Write (6,*) 'Last keyword read was ',WORD
      go to 9900

 9900 continue

      Stop
      end


