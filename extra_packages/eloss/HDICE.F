c
      SUBROUTINE HDICE(vertex, cdir, dist, icell_ptg)

c     $Id: HDICE.F,v 1.3 2014/04/11 01:15:21 fklein Exp $
      IMPLICIT NONE
      SAVE
c
c_begin_var
c  input/output variables:
c  -----------------------
c
      REAL vertex(3)
      REAL cdir(3)
      REAL dist(20)
      INTEGER icell_ptg            ! =1 for g9a, =2 for g9b =3 for g14

c-- DIST array contains calculated pathlengths for all pieces of material downstream of upstream end of HDIce
c--   DIST(1)  target cell wall (KelF)
c--   DIST(2)  IVC wall (KelF)
c--   DIST(3)  part 4 of NMR coil (KelF)
c--   DIST(4)  part 5 of NMR coil (KelF)
c--   DIST(5)  part 7 of NMR coil (kelF)
c--   DIST(6)  part 8 of NMR coil (Al)
c--   DIST(7)  LHe magnet can inner wall (Al)
c--   DIST(8)  SC wire 4 layers 
c--   DIST(9)  Saddle coil 
c--   DIST(10) LHe jacket (L4He)
c--   DIST(11) He magnet can outer wall (Al)
c--   DIST(12) LHe magnet can end plate 
c--   DIST(13) 80K shield shell (Al)
c--   DIST(14) OVC inner wall region A+B (Al)
c--   DIST(15) OVC water jacket (H2O)
c--   DIST(16) OVC outer wall region A+B (Al)
c--   DIST(17) OVC wall region C (Al)
c--   DIST(18) backup coil (region A) (Al)
C--   DIST(19) OVC inner wall region D (Al)



c     Local pre-defined variables:
c     ---------------------------
c
c
c  Local pre-defined variables:
c  ---------------------------
c  RCS information: 
      CHARACTER*132  CFILE, CREVIS, CSTATE, CDATE, CAUTHO, CRCSID
      PARAMETER (CFILE = '$RCSfile: HDICE.F,v $')
      PARAMETER (CREVIS = '$Revision: 1.3 $')
      PARAMETER (CSTATE = '$State: Exp $')
      PARAMETER (CDATE = '$Date: 2014/04/11 01:15:21 $')
      PARAMETER (CAUTHO = '$Author: fklein $')
      DATA CRCSID/   
     1'$Id: HDICE.F,v 1.3 2014/04/11 01:15:21 fklein Exp $'   
     2/   
c  Module information:
      CHARACTER*(*)  CRNAME, CRAUTH
      CHARACTER*100  CRMESS
      PARAMETER (CRNAME='FROST')
      PARAMETER (CRAUTH='Brian Morrison')
c
c
c
      REAL TWOPI
      PARAMETER (TWOPI = 2.0 * 3.14159265)
      REAL dir(3)

      REAL RC(35)               ! radius of cylinders
      REAL ZPOS(26)             ! z-positions of cylinders


c fk: Data look like all Z-positions 8mm upstream 
c     (correction after data statement)
c
c--   GEOMETRY TAKEN FROM HDICE technical note 22
c--
c--   ALL DISTANCES ARE IN CENTIMETERS
c--
c--   RC(1):  Target Cell inner radius  = 1.6/2 - 0.05 = 0.75
c--   RC(2):  Target Cell outer radius = 1.6/2 = 0.8       
c--   RC(3):  IVC inner radius = 1.9/2 - 0.05 = 0.9
c--   RC(4):  IVC outer radius = 1.9/2 = 0.95 
c---      NMR coil see page 5 of the HDICE technical note 22
c--   RC(5):  part 4 inner radius = 4.308/2 - 0.122 = 2.032
c--   RC(6):  part 4 outer radius = 4.308/2 = 2.154
c--   RC(7):  part 5 inner radius = 4.674/2 - 0.305 = 2.032
c--   RC(8):  part 5 outer radius = 4.674/2 = 2.337 
c--   RC(9)   part 7 inner radius = 6.96/2 - 0.178 = 3.302
c--   RC(10)  part 7 outer radius = 6.96/2 = 3.48
c--   RC(11)  part 8 inner radius = 6.96/2 - 0.051 = 3.429
c---  RC(12)  part 8 outer radius = 6.96/2 = 3.48
c--
c--   RC(13): LHe magnet can inner wall inner radius = 3.49    
c--   RC(14): LHe magnet can inner wall outer radius = 3.566        
c--   RC(15): SC wire 4 layers inner radius = 3.566      
c--   RC(16): SC wire 4 layers outer radius = 3.566 + 0.0254*4*0.65 = 3.632
c--   RC(17)  Saddle coil outer radius  = 3.632 + 0.025 = 3.657
c--   RC(18): LHE jacket outer inner readius = 3.668
c--   RC(19): LHe jacket outer radius = 4.445  
c--   RC(20): LHe magnet can outer wall inner radius = 4.445       
c--   RC(21): LHe magnet can outer wall outer radius = 4.445 + 0.0813 = 4.5263
c--
c--   RC(22): 80K shield shell inner radius = 5.398
c--   RC(23): 80K shield shell inner radius = 5.398 + 0.0508 = 5.4488
c--
c--   RC(24): OVC inner wall region A+B inner radius = 6.668 
c--   RC(25): OVC inner wall region A+B outer radius = 6.668 + 0.1 = 6.768
c--   RC(26): OVC water jacket inner radius = 6.768
c--   RC(27): OVC water jacket inner radius = 6.768 + 0.1 = 6.868
c--   RC(28): OVC outer wall region A+B inner radius = 6.868
c--   RC(29): OVC outer wall region A+B outer radius = 6.868 + 0.1 = 6.968
c--
c--   RC(30): backup coil inner radius = 6.968
c--   RC(31): backup coil region outer radius = 6.968 + 0.1024*0.747 = 7.0445
c--
c--   RC(32): OVC wall region C inner radius = 6.668
c--   RC(33): OVC wall region C outer radius = 6.868 + 0.38 = 7.048
c--   RC(34): OVC wall region D inner radius = 6.668
c--   RC(35): OVC wall region D outer radius = 6.668 + 0.1 = 6.768
c--



      DATA RC/
     &     0.75,               ! RC1
     &     0.8,                 ! RC2
     &     0.9,                 ! RC3
     &     0.95,                ! RC4
     &     2.032,               ! RC5
     &     2.154,               ! RC6
     &     2.032,               ! RC7
     &     2.337,               ! RC8
     &     3.302,               ! RC9
     &     3.48,                ! RC10
     &     3.429,               ! RC11
     &     3.48,                ! RC12
     &     3.49,                ! RC13
     &     3.566,               ! RC14
     &     3.566,               ! RC15
     &     3.632,               ! RC16
     &     3.657,               ! RC17
     &     3.668,               ! RC18
     &     4.445,               ! RC19
     &     4.445,               ! RC20
     &     4.5623,              ! RC21
     &     5.398,               ! RC22
     &     5.4488,              ! RC23
     &     6.668,               ! RC24
     &     6.768,               ! RC25
     &     6.768,               ! RC26
     &     6.868,               ! RC27
     &     6.868,               ! RC28
     &     6.968,               ! RC29
     &     6.968,               ! RC30
     &     7.0445,              ! RC31
     &     6.668,               ! RC32
     &     7.048,               ! RC33
     &     6.668,               ! RC34
     &     6.768               ! RC35
     &     /            
      


c--   IMPORTANT NOTE: THIS CODE ONLY LOOKS AT THE GEOMETRY THAT BEGINS AT THE upstream end of HDICE  z = -4.129 - 5 = -9.129 cm)
c--   HDICE center is at -6.629 

c--   ZPOS(1):  Target cell end cap start  = -4.129
c--   ZPOS(2):  Target cell end cup finish = -4.129 + 0.05 = -4.079
c--
c--   ZPOS(3):  IVC end cap start          = -6.629 + 8.5  = 1.871
c--   ZPOS(4):  IVC end cap finish         = 1.871 + 0.05  = 1.921
c--
c--   ZPOS(5):  NMR part 4 end    = -4.129 + 0.66          = -3.469
c--   ZPOS(6):  NMR part 5 end    = -3.469 + 0.927         = -2.542
c--   ZPOS(7):  NMR part 7a start = -4.129 + 2.629         = -1.5
c--   ZPOS(8):  NMR part 7b start = -1.5   + 0.356         = -1.144
c--   ZPOS(9):  NMR part 7c start = -1.144 + 0.152         = -0.992
c--   ZPOS(10): NMR part 7c end   = -0.992 + 0.610         = -0.382
c--   ZPOS(11): NMR part 8 end    = -0.382 + 16.281        = 15.899 -> trim to 15.2
c--
c--   ZPOS(12): LHe magnet can end = -6.629 + 21.88        = 15.251
c--   ZPOS(13): SC coil end        = -6.629 + 20           = 13.371
c--
c--   ZPOS(14): Saddle coil start = -6.629+ + 12.93/2      = -0.164
c--   ZPOS(15): Saddle coil end   = -0.164+ 1.62           = 1.456
  
c--   ZPOS(16): 80K shield end     = -6.629 + 26.14        = 19.511
c--
c--   ZPOS(17): OVC region A (backup coil) end = -4.129 + 12.019 = 7.59
c--   ZPOS(18): OVC region A+B end             = -4.129 + 13.528 = 9.099
c--   ZPOS(19): OVC region C end               = -4.129 + 14.163 = 9.734
c--   ZPOS(20): OVC region D end               = -4.129 + 26.306 = 21.877
c--
c--
      DATA ZPOS/
     &     -4.129,              !ZPOS1
     &     -4.079,              !ZPOS2
     &     1.871,               !ZPOS3
     &     1.921,               !ZPOS4
     &     -3.469,              !ZPOS5
     &     -2.542,              !ZPOS6
     &     -1.5,                !ZPOS7
     &     -1.144,              !ZPOS8
     &     -0.992,              !ZPOS9
     &     -0.382,              !ZPOS10
     &     15.2,                !ZPOS11
     &     15.251,              !ZPOS12
     &     13.371,              !ZPOS13
     &     -0.164,              !ZPOS14
     &     1.456,               !ZPOS15
     &     19.511,              !ZPOS16
     &     7.59,                !ZPOS17
     &     9.099,               !ZPOS18
     &     9.734,               !ZPOS19
     &     21.877,              !ZPOS20
     &     0.,                  !ZPOS21
     &     0.,                  !ZPOS22
     &     0.,                  !ZPOS23
     &     0.,                  !ZPOS24
     &     0.,                  !ZPOS25
     &     0.                   !ZPOS26
     &     /

      DATA dir /0., 0., 1./     !target axis direction cosines
      
      REAL pathLength1, pathLength2, pathLength,intercept(3)
      INTEGER ierrD1,ierrD2,ierrOuter,ierrInner,exitWindowR
      INTEGER i,hit1,hit2,ierr
      REAL  rInner, rOuter, rDisk, zLow, zHigh,z0,z0Tmp
      REAL scYwidth
      REAL  vc1(3),vc2(3),pos1(3),pos2(3),pos1w(3),pos2w(3)
      REAL a1(3),a2(3),center(3),d,dSq,czAbs
      DATA center/0.,0.,0./
      DATA pos1/0.,0.,0./
      DATA pos2/0.,0.,0./
      DATA pos1w/0.,0.,0./
      DATA pos2w/0.,0.,0./
      DATA vc1/0.,0.,21.427754295/
      DATA vc2/0.,0.,22.062002433/

      REAL myZ, myPhi
      REAL phiCen, rCylIn, rCylOut
      REAL zCen, zTrack, zCirIn, zCirOut
      REAL dS, dZ
      LOGICAL Lfirstcall
      DATA Lfirstcall /.TRUE./

c     _end_var

      if(Lfirstcall) then
         Lfirstcall=.FALSE.
         do i=1,20
            ZPOS(i) = ZPOS(i)-0.8
         enddo
      endif

      do i=1,20
         dist(i) = 0.0
      enddo
      z0 = -9.129 -0.8

c--  Target cell Wall     

      rInner = RC(1)
      rOuter = RC(2)
      zLow = ZPOS(1)
      zHigh = ZPOS(2)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zLow,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(1) = pathLength

c--   IVC 

      rInner = RC(3)
      rOuter = RC(4)
      zLow = ZPOS(3)
      zHigh = ZPOS(4)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zLow,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(2) = pathLength

c--   Part 4 of NMR coil

      rInner = RC(5)
      rOuter = RC(6)
      zHigh = ZPOS(5)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(3) = pathLength

c--   Part 5 of NMR 

      rInner = RC(7)
      rOuter = RC(8)
      zLow   = ZPOS(5)
      zHigh  = ZPOS(6)
      
      call cylPath(vertex,cdir,rOuter,zLow,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(4) = pathLength
         
c--   Part 7 of NMR 

      rInner = RC(9)
      rOuter = RC(10)
      zLow   = ZPOS(7)
      zHigh  = ZPOS(10)
      
      call cylPath(vertex,cdir,rOuter,zLow,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif
      
      dist(5) = pathLength
         
c--   Part 8 of NMR 
      
      rInner = RC(11)
      rOuter = RC(12)
      zLow = ZPOS(10)
      zHigh = ZPOS(11)
      
      call cylPath(vertex,cdir,rOuter,zLow,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif
      
      dist(6) = pathLength

c--   LHe magnet can inner wall

      rInner = RC(13)
      rOuter = RC(14)
      zHigh = ZPOS(12)
      
      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif
      
      dist(7) = pathLength

c--   Superconducting wire

      rInner = RC(15)
      rOuter = RC(16)
      zHigh = ZPOS(13)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(8) = pathLength

C--  Saddle coil
      
      rInner = RC(16)
      rOuter = RC(17)
      zHigh = ZPOS(15)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)

      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
         call CylInt(vertex, cdir, center, dir, rInner, a1, a2, ierr)        
         
         if(ierr .ne. 0) then
            if(a1(3) .lt. ZPOS(14) .and. abs(a1(2)) .gt. 1.62) then ! goes through cutout  
               pathLength = 0.0
            else
               pathLength = 0.0
            endif
         else
            pathLength = 0.0
         endif
      else 
         pathLength = 0.0         
      endif

      dist(9) = pathLength


c--   LHe jacket

      rInner = RC(18)
      rOuter = RC(19)
      zHigh = ZPOS(12)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(10) = pathLength

c--   LHe magnet can outer wall

      rInner = RC(20)
      rOuter = RC(21)
      zHigh = ZPOS(12)
      
      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif
      
      dist(11) = pathLength

c--   LHe magnet can endplate

      rInner = RC(13)
      rOuter = RC(20)
      zLow = ZPOS(12)
      zHigh = zLow + 0.076

      call cylPath(vertex,cdir,rOuter,zLow,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif
      
      dist(12) = pathLength

c--   80K shield

      rInner = RC(22)
      rOuter = RC(23)
      zHigh = ZPOS(16)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(12) = pathLength

c--   OVC inner wall region A+B

      rInner = RC(24)
      rOuter = RC(25)
      zHigh = ZPOS(18)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(14) = pathLength

c--   OVC water jacket region A+B

      rInner = RC(26)
      rOuter = RC(27)
      zHigh = ZPOS(18)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(15) = pathLength

c--   OVC outer wall region A+B

      rInner = RC(28)
      rOuter = RC(29)
      zHigh = ZPOS(18)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(16) = pathLength

c--   OVC wall region C

      rInner = RC(24)
      rOuter = RC(29)
      zLow   = ZPOS(18)
      zHigh  = ZPOS(19)

      call cylPath(vertex,cdir,rOuter,zLow,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(17) = pathLength

c--   OVC inner wall region D

      rInner = RC(24)
      rOuter = RC(25)
      zLow   = ZPOS(19)
      zHigh  = ZPOS(20)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(19) = pathLength

c--   backup coil

      rInner = RC(30)
      rOuter = RC(31)
      zHigh = ZPOS(17)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(18) = pathLength


      RETURN
      end
