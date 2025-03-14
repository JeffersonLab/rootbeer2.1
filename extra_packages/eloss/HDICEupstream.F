      SUBROUTINE HDICEupstream(vertex, cdir, ups)

c  $Id: HDICEupstream.F,v 1.3 2014/04/11 01:15:21 fklein Exp $
      IMPLICIT NONE
      SAVE
c     
c     _begin_var
c     input/output variables:
c     -----------------------
c     
      REAL vertex(3)
      REAL cdir(3)
      REAL ups(20)


c--   UPS array contains calculated pathlengths for all pieces of material upstream end of HDIce
c--   UPS(1)  target cell wall (KelF)
c--   UPS(2)  IVC wall (KelF)
c--   UPS(3)  part 1 of NMR coil (KelF)
c--   UPS(4)  part 3 of NMR coil (KelF)
c--   UPS(5)  part 4 of NMR coil (kelF)
c--   UPS(6)  LHe magnet can inner wall (Al)
c--   UPS(7)  SC wire 4 layers 
c--   UPS(8)  Saddle coil 
c--   UPS(9)  LHe jacket (L4He)
c--   UPS(10) He magnet can outer wall (Al)
c--   UPS(11) 80K shield shell (Al)
c--   UPS(12) OVC inner wall region A+B (Al)
c--   UPS(13) OVC water jacket (H2O)
c--   UPS(14) OVC outer wall region A+B (Al)
c--   UPS(15) backup coil (region A) (Al)


c     Local pre-defined variables:
c     ---------------------------
c     

      REAL RC(35)               ! radius of cylinders
      REAL ZPOS(15)             ! z-positions of cylinders


c fk: Data look like all Z-positions 8mm upstream 
c     (correction after data statement)
c
c--   GEOMETRY TAKEN FROM HDICE technical note 22
c--   
c--   IMPORTANT NOTE: THIS CODE ONLY LOOKS AT THE GEOMETRY THAT ENDS AT THE upstream end HDICE (z0 = -9.129cm)
c--   
c--   ALL DISTANCES ARE IN CENTIMETERS
c--   
c--   RC(1):  Target Cell inner radius  = 1.6/2 - 0.05 = 0.795
c--   RC(2):  Target Cell outer radius = 1.6/2 = 0.8       
c--   RC(3):  IVC inner radius = 1.9/2 - 0.05 = 0.9
c--   RC(4):  IVC outer radius = 1.9/2 = 0.95 
c--   RC(5):  IVC upstream part inner radius = 
c--   RC(6):  IVC upstream part outer radius = 
c---  NMR coil see page 5 of the HDICE techical note 22
c--   RC(7):  part 4 inner radius = 4.308/2 - 0.122 = 2.032
c--   RC(8):  part 4 outer radius = 4.308/2 = 2.154
c--   RC(9):  part 3 inner radius = 4.674/2 - 0.305 = 2.032
c--   RC(10): part 3 outer radius = 4.674/2 = 2.337 
c--   RC(11): part 1 inner radius = 6.96/2 - 0.178 = 3.302
c--   RC(12): part 1 outer radius = 6.96/2 = 3.48
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

      DATA RC/
     &     0.75,                ! RC1
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
     &     6.768                ! RC35
     &     /            

      

c--   IMPORTANT NOTE: THIS CODE ONLY LOOKS AT THE GEOMETRY THAT ENDS AT THE upstream end of HDICE  z0 = -4.129 - 5 = -9.129 cm)
c--   HDICE center is at -6.629 

c--   ZPOS(1):  Target cell start = -9.129 - 1.            = -10.129
c--   
c--   ZPOS(2):  IVC end cap start = -9.129 - 1.            = -10.129
c--   
c--   ZPOS(3):  NMR part 4 start  = -4.129 - 5.662         = -9.791
c--   ZPOS(4):  NMR part 3 start  = -4.129 - 6.589         = -10.781
c--   ZPOS(5):  NMR part 1 start  = -4.129 - 8.138         = -12.267
c--   ZPOS(6):  NMR part 1 end    = -12.267 + 0.508         = -11.759
c--   
c--   ZPOS(7):  LHe magnet can start = -6.629 - 20.        = -26.629
c--   ZPOS(8): SC coil start        = -6.629 - 20.        = -26.629
c--   
c--   ZPOS(9): Saddle coil start = -6.629 - 12.93/2 -1.62   = -14.714
c--   ZPOS(10): Saddle coil end   = -6.629 - 12.93/2           = -13.094
      
c--   ZPOS(11): 80K shield start     =                           = -20.
c--   
c--   ZPOS(12): OVC region A (backup coil) estart = -6.629 - 13.97 = -20.599






      DATA ZPOS/
     &     -10.129,             !ZPOS(1)
     &     -10.129,             !ZPOS(2) 
     &     -9.791,              !ZPOS(3) 
     &     -10.781,             !ZPOS(4) 
     &     -12.267,             !ZPOS(5)
     &     -11.759,             !ZPOS(6) 
     &     -26.629,             !ZPOS(7) 
     &     -26.629,             !ZPOS(8)
     &     -14.714,             !ZPOS(9)
     &     -13.094,             !ZPOS(10) 
     &     -20.0,               !ZPOS(11)
     &     -20.599,             !ZPOS(12)
     &     0.0,                 !ZPOS(13)
     &     0.0,                 !ZPOS(14)
     &     0.0                  !ZPOS(15)
     &     /


c     RCS information: 

c     CHARACTER*132  CFILE, CREVIS, CSTATE, CDATE, CAUTHO, CRCSID
c     PARAMETER (CFILE = '$RCSfile: HDICEupstream.F,v $')
c     PARAMETER (CREVIS = '$Revision: 1.3 $')
c     PARAMETER (CSTATE = '$State: Exp $')
c     PARAMETER (CDATE = '$Date: 2014/04/11 01:15:21 $')
c     PARAMETER (CAUTHO = '$Author: fklein $')
c     DATA CRCSID/   
c     1     '$Id: HDICEupstream.F,v 1.3 2014/04/11 01:15:21 fklein Exp $'/
      
c     Module information:

c     CHARACTER*(*)  CRNAME, CRAUTH
c     CHARACTER*100  CRMESS
c     PARAMETER (CRNAME='FROST')
c     PARAMETER (CRAUTH='Brian Morrison')
c     
c     Local User defined variables:
c     -----------------------------
c     


      REAL pathLength1, pathLength2, pathLength,intercept(3)
      INTEGER ierrD1,ierrD2,ierrOuter,ierrInner,i,ierr
      REAL  rInner, rOuter, rDisk, zLow, zHigh,z0,z0Tmp
      REAL a1(3),a2(3),center(3),dir(3)

      DATA dir /0., 0., 1./     !target axis direction cosines
      DATA center/0.,0.,0./
      LOGICAL lfirstcall
      DATA lfirstcall /.TRUE./

c     _end_var

      if(lfirstcall) then
         lfirstcall=.FALSE.
         do i=1,12
            ZPOS(i) = ZPOS(i) -0.8
         enddo
      endif

      z0 = -9.129 -0.8
      do i=1,20
         ups(i) = 0.0
      enddo

c--   Target cell Wall     

      rInner = RC(1)
      rOuter = RC(2)
      zLow = ZPOS(1)
      
      call cylPath(vertex,cdir,rOuter,zLow,z0,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,z0,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      ups(1) = pathLength
      
c--   IVC 

      rInner = RC(3)
      rOuter = RC(4)
      zLow = ZPOS(2)


      call cylPath(vertex,cdir,rOuter,zLow,z0,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,z0,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      ups(2) = pathLength

c--   Part 4 of NMR coil

      rInner = RC(5)
      rOuter = RC(6)
      zLow = ZPOS(3)

      call cylPath(vertex,cdir,rOuter,zLow,z0,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,z0,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      ups(3) = pathLength

c--   Part 3 of NMR 

      rInner = RC(7)
      rOuter = RC(8)
      zLow   = ZPOS(4)
      zHigh  = ZPOS(3)
      
      call cylPath(vertex,cdir,rOuter,zLow,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      ups(4) = pathLength
      
c--   Part 1 of NMR 

      rInner = RC(9)
      rOuter = RC(10)
      zLow   = ZPOS(5)
      zHigh  = ZPOS(6)
      
      call cylPath(vertex,cdir,rOuter,zLow,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif
      
      ups(5) = pathLength

c--   LHe magnet can inner wall

      rInner = RC(13)
      rOuter = RC(14)
      zLow = ZPOS(7)
      
      call cylPath(vertex,cdir,rOuter,zLow,z0,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,z0,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif
      
      ups(6) = pathLength

c--   Superconducting wire

      rInner = RC(15)
      rOuter = RC(16)
      zLow = ZPOS(8)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      ups(7) = pathLength

C--   Saddle coil
      
      rInner = RC(16)
      rOuter = RC(17)
      zLow = ZPOS(9)

      call cylPath(vertex,cdir,rOuter,zLow,z0,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,z0,pathLength1,ierrInner)

      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
         call CylInt(vertex, cdir, center, dir, rInner, a1, a2, ierr)        
         if(ierr .ne. 0) then
            if(a1(3) .gt. ZPOS(10) .and. abs(a1(2)) .gt. 1.62) then ! goes through cutout  
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
      

      ups(8) = pathLength


c--   LHe jacket

      rInner = RC(18)
      rOuter = RC(19)
      zLow = ZPOS(7)

      call cylPath(vertex,cdir,rOuter,zLow,z0,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,z0,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      ups(9) = pathLength

c--   LHe magnet can outer wall

      rInner = RC(20)
      rOuter = RC(21)
      zLow = ZPOS(7)
      
      call cylPath(vertex,cdir,rOuter,zLow,z0,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,z0,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif
      
      ups(10) = pathLength

c--   80K shield

      rInner = RC(22)
      rOuter = RC(23)
      zLow = ZPOS(11)

      call cylPath(vertex,cdir,rOuter,zLow,z0,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,z0,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      ups(11) = pathLength

c--   OVC inner wall region A

      rInner = RC(24)
      rOuter = RC(25)
      zLow = ZPOS(12)

      call cylPath(vertex,cdir,rOuter,zLow,z0,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,z0,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      ups(12) = pathLength

c--   OVC water jacket region A

      rInner = RC(26)
      rOuter = RC(27)
      zLow = ZPOS(12)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      ups(13) = pathLength

c--   OVC outer wall region A

      rInner = RC(28)
      rOuter = RC(29)
      zLow = ZPOS(12)

      call cylPath(vertex,cdir,rOuter,zLow,z0,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,z0,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      ups(14) = pathLength

c--   backup coil

      rInner = RC(30)
      rOuter = RC(31)
      zLow = ZPOS(12)

      call cylPath(vertex,cdir,rOuter,zLow,z0,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,zLow,z0,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      ups(15) = pathLength

      RETURN
      end

