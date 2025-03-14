c
      SUBROUTINE FROST(vertex, cdir, dist, icell_ptg)


      IMPLICIT NONE
      SAVE
c
c_begin_var
c  input/output variables:
c  -----------------------
c
      REAL vertex(3)
      REAL cdir(3)
      REAL dist(17)
      INTEGER icell_ptg            ! =1 for g9a, =2 for g9b

c     Local pre-defined variables:
c     ---------------------------
c
c
c  Local pre-defined variables:
c  ---------------------------
c  RCS information: 
      CHARACTER*132  CFILE, CREVIS, CSTATE, CDATE, CAUTHO, CRCSID
      PARAMETER (CFILE = '$RCSfile: FROST.F,v $')
      PARAMETER (CREVIS = '$Revision: 1.10 $')
      PARAMETER (CSTATE = '$State: Exp $')
      PARAMETER (CDATE = '$Date: 2012/04/13 11:13:27 $')
      PARAMETER (CAUTHO = '$Author: fklein $')
      DATA CRCSID/   
     1'$Id: FROST.F,v 1.10 2012/04/13 11:13:27 fklein Exp $'   
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

      REAL RC(22)               ! radius of cylinders
      REAL ZPOS(26)             ! z-positions of cylinders
      REAL ZPOS1K_G9B(9:12)     ! for 1k pot (different because of transv. holding coils)
      REAL PHIPOS(1)            ! angular data for holding coil


c--   GEOMETRY TAKEN FROM http://www.jlab.org/~zarecky/Frozen%20Spin/66840-04674-04.pdf
c--   AND FROM http://www.jlab.org/~zarecky/Frozen%20Spin/Scattering%20Chamber/66850-04635.pdf
c--   Transverse Holding Coil (THC) geometry taken from email (INCLUDE LINK!!!!!!!!!!!!!!!!!!)
c--
c--   IMPORTANT NOTE: THIS CODE ONLY LOOKS AT THE GEOMETRY THAT BEGINS AT THE BUTANOL TSH (z = -2.64cm)
c--
c--   ALL DISTANCES ARE IN CENTIMETERS
c--
c--   RC(1):  Target Sample Holder (TSH) inner radius = (0.650-2*0.025)*2.54/2 = 0.762
c--   RC(2):  TSH outer radius = 0.650*2.54/2 = 0.8255        
c--   RC(3):  Mixing Chamber (MC) inner radius = (0.9-2*0.025)*2.54/2 = 1.0795
c--   RC(4):  Mixing Chamber (MC) outer radius = 0.9*2.54/2 = 1.143 
c-g9a  RC(5):  1K inner radius1 = 1.855*2.54/2 = 2.35585
c-g9a  RC(6):  1K outer radius1 = (1.969-2*0.017)*2.54/2 = 2.45745
c-g9a  RC(7):  1K inner radius2 = 2.35585 
c-g9a  RC(8):  1K outer radius2 = 2.44475
c--   RC(9)   1K inner radius3 = 2.35585
c--   RC(10)  1K outer radius3 = 2.50063
c-g9a  RC(11)  1K dummy var1 = 2.35585
c-g9a  RC(12)  1K dummy var2 = 2.50063
c-g9a  RC(13): Super conducting wire inner radius1 (3-layer) = RC(6) = 2.45745         
c-g9a  RC(14): Super conducting wire outer radius1 = 1.969*2.54/2 = 2.50063        
c-g9a  RC(15): Super conducting wire inner radius2 (4-layer) = RC(8) = 2.44475         
c-g9a  RC(16): Super conducting wire outer radius2 = 1.969*2.54/2 = 2.50063
c--   RC(17): 20K inner radius = (2.517-0.040)*2.54/2 - 10*0.000631957 = 3.1395.(SI layers=0.000631957)
c--   RC(18): 20K outer radius = 2.517*2.54/2 = 3.19659         
c--   RC(19): Scattering chamber inner radius = 3.510*2.54/2 = 4.4577        
c--   RC(20): Scattering chamber outer radius = 4.375*2.54/2 = 5.55625      
c--   RC(21): THC shell inner radius = 1.969*2.54/2
c--   RC(22): THC shell outer radius = (1.969+0.070)*2.54/2
c--
      DATA RC/0.762, 0.8255, 1.0795, 1.143, 2.35585, 2.45745, 2.35585,
     &     2.44475, 2.35585, 2.50063, 2.35585, 2.50063, 2.45745, 2.50063, 
     &     2.44475, 2.50063, 3.1395, 3.19659, 4.4577, 5.55625,
     &     2.50063, 2.50763/

c--   ZPOS(1):  TSH end cap start = 2.64
c--   ZPOS(2):  TSH end cap finish = 2.70
c--   ZPOS(3):  Mixing Chamber (MC) end cap start = 3.2
c--   ZPOS(4):  Mixing Chamber (MC) end cap finish = 3.2 + 0.064 = 3.264
c-g9a  ZPOS(5):  1k start1 = -2.64  
c-g9a  ZPOS(6):  1k finish1 = 3.02707
c-g9a  ZPOS(7):  1k start2 = 3.02707  
c-g9a  ZPOS(8):  1k finish2 = 5.29529
c-g9a  ZPOS(9):  1k start3 = 5.29529  
c-g9b  ZPOS1K_G9B(9):  1k start3 = -2.64 (changed)
c-g9a  ZPOS(10): 1k finish3 = 6.15 
c-g9b  ZPOS1k_G9B(10): 1k finish3 = 12.5 (changed)
c-g9a  ZPOS(11): 1K end cap start = 6.15
c-g9b  ZPOS1K_G9B(11): 1K end cap start = 12.5 (changed)
c-g9a  ZPOS(12): 1K end cap finish = 6.15 + 0.02 = 6.17
c-g9b  ZPOS1K_G9B(12): 1K end cap finish = 12.5 + 0.02 = 12.52 (changed)
c-g9a  ZPOS(13): 3-layer wire cylinder start = -2.64  
c-g9a  ZPOS(14): 3-layer wire cylinder finish = 3.02707
c-g9a  ZPOS(15): 4-layer wire cylinder start = 3.02707
c-g9a  ZPOS(16): 4-layer wire cylinder finish = 5.29529
c--   ZPOS(17): 20K end cap start = 15.95
c--   ZPOS(18): 20K end cap finish = 15.95 + 0.02 = 15.97
c--   ZPOS(19): Scattering chamber cylinder start = -2.64
c--   ZPOS(20): Scattering chamber cylinder finish = 18.9934
c--   ZPOS(21): exit window start = 21.61
c--   ZPOS(22): exit window finish = 21.61 + 0.025 = 21.635
c--   ZPOS(23): THC z center = 5" behind end of 1k shield = 12.5 - 2.54*5 = -0.2
c--   ZPOS(24): THC half length of racetrack = 7.874*2.54-pi*2.50063/2 = 6.0720
c--   ZPOS(25): THC radius of curvature of inner coil = 2.54 * 0.715 / 2
c--   ZPOS(26): THC radius of curvature of outer coil = half the cylinder = pi R
c--
      DATA ZPOS/2.64, 2.7, 3.2, 3.264, -2.64, 3.02707, 3.02707, 5.29529,
     &     5.29529, 6.15, 6.15, 6.17, -2.64, 3.02707, 3.02707, 5.29529,
     &     15.95, 15.97, -2.64, 18.9934, 21.61, 21.635, 
     &     -0.200, 6.0720, 0.0719, 0.9055/
      DATA ZPOS1K_G9B/ -2.64, 12.5, 12.5, 12.52 /

c--   PHIPOS: Radian position of center of one coil
      DATA PHIPOS / 0.5236 /

      DATA dir /0., 0., 1./     !target axis direction cosines
      
      REAL pathLength1, pathLength2, pathLength,intercept(3)
      INTEGER ierrD1,ierrD2,ierrOuter,ierrInner,exitWindowR
      INTEGER i,hit1,hit2,ierr
      REAL  rInner, rOuter, rDisk, zLow, zHigh,z0,z0Tmp
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
c     _end_var

      do i=1,17
         dist(i) = 0.0
      enddo

      if((icell_ptg.lt.1) .or. (icell_ptg.gt.2)) return

      if(icell_ptg.eq.2) then   ! g9b
         do i=9,12
            ZPOS(i) = ZPOS1K_G9B(i)
         enddo
      endif
      
      z0 = -2.64
      exitWindowR = 2.54/2

c--   TSH (Target Sample Holder)     

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

c--   MC (Mixing Chamber interior)

      rInner = RC(2)
      rOuter = RC(3)
      zLow = ZPOS(2)
      zHigh = ZPOS(3)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zLow,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(2) = pathLength

c--   MC (Mixing Chamber cylinder and end caps)

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

      dist(3) = pathLength


c--   1K stage 1 
c--
c--   NOTE: The 1K containment surface also acts as the support structure for the
c--   superconducting holding coils.
c-g9a Stage 1 is defined as the part of the 1K cylinder that has z>= -2.64 cm and 
c-g9a supports the 3-layers of wire. Stage 2 supports the 4-layers of wire. 
c--   Stage 3 does not provide support for any wires.

      if(icell_ptg.eq.1) then
         rInner = RC(5)
         rOuter = RC(6)
         zLow = ZPOS(5)
         zHigh = ZPOS(6)
         
         call cylPath(vertex,cdir,rOuter,zLow,zHigh,pathLength2,ierrOuter)
         call cylPath(vertex,cdir,rInner,zLow,zHigh,pathLength1,ierrInner)
         if (ierrOuter.ne.0) then
            pathLength = pathLength2 - pathLength1
         else
            pathLength = 0.0
         endif

         dist(4) = pathLength
         
c--   1K stage 2

         rInner = RC(7)
         rOuter = RC(8)
         zLow = ZPOS(7)
         zHigh = ZPOS(8)
         
         call cylPath(vertex,cdir,rOuter,zLow,zHigh,pathLength2,ierrOuter)
         call cylPath(vertex,cdir,rInner,zLow,zHigh,pathLength1,ierrInner)
         if (ierrOuter.ne.0) then
            pathLength = pathLength2 - pathLength1
         else
            pathLength = 0.0
         endif

         dist(5) = pathLength

      endif    ! g9a

c--   1K stage 3

      rInner = RC(9)
      rOuter = RC(10)
      zLow = ZPOS(11)
      zHigh = ZPOS(12)
      z0Tmp = ZPOS(9)

      call cylPath(vertex,cdir,rOuter,z0Tmp,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0Tmp,zLow,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(6) = pathLength

      if(icell_ptg.eq.1) then

c-g9a  Super conducting wire (3 layer)

         rInner = RC(13)
         rOuter = RC(14)
         zLow = ZPOS(13)
         zHigh = ZPOS(14)

         call cylPath(vertex,cdir,rOuter,zLow,zHigh,pathLength2,ierrOuter)
         call cylPath(vertex,cdir,rInner,zLow,zHigh,pathLength1,ierrInner)
         if (ierrOuter.ne.0) then
            pathLength = pathLength2 - pathLength1
         else
            pathLength = 0.0
         endif

         dist(7) = pathLength

c-g9a  Super conducting wire (4 layer)

         rInner = RC(15)
         rOuter = RC(16)
         zLow = ZPOS(15)
         zHigh = ZPOS(16)

         call cylPath(vertex,cdir,rOuter,zLow,zHigh,pathLength2,ierrOuter)
         call cylPath(vertex,cdir,rInner,zLow,zHigh,pathLength1,ierrInner)
         if (ierrOuter.ne.0) then
            pathLength = pathLength2 - pathLength1
         else
            pathLength = 0.0
         endif

         dist(8) = pathLength

      endif           ! g9a

c--   20k

      rInner = RC(17)
      rOuter = RC(18)
      zLow = ZPOS(17)
      zHigh = ZPOS(18)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zLow,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(9) = pathLength


c--   SC (Scattering Chamber cylinder wall)

      rInner = RC(19)
      rOuter = RC(20)
      zHigh = ZPOS(20)

      call cylPath(vertex,cdir,rOuter,z0,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rInner,z0,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif

      dist(10) = pathLength

c-- Exit window

      rOuter = exitWindowR
      zLow = ZPOS(21)
      zHigh = ZPOS(22)
      call diskIntercept(vertex,cdir,rOuter,zLow,intercept,ierrD1)
      czAbs = sqrt(cdir(3)*cdir(3))
      if (ierrD1.eq.1) then
         dist(11) = (zHigh - zLow)/czAbs
      endif

c--   Nose end cap
c--   GEOMETRY: 
c--   Conic structure with apex1 = (0.0, 0.0, 21.4278), apex2 = (0.0, 0.0, 22.0620)
c--   Angle that conic structure makes with the z-axis is 120 degrees


c--   FIRST CHECK IF TRACK GOES THROUGH NOSE END CAP

      rOuter = RC(19)
      zLow = ZPOS(20)
      call diskIntercept(vertex,cdir,rOuter,zLow,intercept,ierrD2)

c--   MAKE SURE THAT TRACK GOES THROUGH NOSE END CAP BUT NOT EXIT WINDOW

      if (ierrD2.eq.1.and.ierrD1.eq.0) then
c--   Need to get another point on track to send to IntersectionCone subroutine
         rInner = RC(19)
         call CylInt(vertex, cdir, center, dir, rInner, a1, a2, ierr)
         
c--   NOTE: Conic surface is at 60 degrees with respect to z-axis 
         
         call IntersectionCone(vertex, a1, vc1, dir, 60.0,
     &        pos1, pos2, hit1)
         call IntersectionCone(vertex, a1, vc2, dir, 60.0,
     &        pos1w, pos2w, hit2)
         
         dSq = 0
         do i=1,3
            dSq = dSq + (pos1w(i) - pos1(i))**2
         enddo
         d = sqrt(dSq)
         dist(12) = d
      endif

      if(icell_ptg.eq.1) return      ! g9a

c-g9b   Transverse holding coil
c-g9b   Strategy: create cylinder and 'cut' non-coil parts out of it

      rCylIn  = RC(21)
      rCylOut = RC(22)
      phiCen  = PHIPOS(1)
      zCen    = ZPOS(23)
      zTrack  = ZPOS(24)
      zCirIn  = ZPOS(25)
      zCirOut = ZPOS(26)
      
c-g9b   Check z0 and zHigh, zLow here!
      zLow    = zCen - zTrack - zCirOut
      zHigh   = zCen + zTrack + zCirOut
      call cylPath(vertex,cdir,rCylOut,zLow,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rCylIn, zLow,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif
      
c-g9b  Find position of hit on inner cylinder: myZ, myPhi
      call CylInt(vertex, cdir, center, dir, rCylIn, a1, a2, ierr)
      if (ierr .eq. 0) then
         pathLength = 0.0
         RETURN
      end if
      myZ   = a1(3)
      myPhi = atan(a1(2)/a1(1)) ! image included because atan returns acute angle!

c-g9b   Snip anything beyond outside azimuthal range
      if (rCylIn*abs(myPhi-phiCen) .gt. zCirOut)
     &     pathLength = 0.0
      
c-g9b   Racetracks
      if ( (abs(myZ-zCen) .lt. zTrack) .and.
     &     (rCylIn*abs(myPhi-phiCen) .lt. zCirIn) ) then
         pathLength = 0.0
c-g9b   Arcs
      elseif (zCen-myZ .gt. zTrack) then
         dZ = -myZ + zCen - zTrack
         dS = rCylIn * abs(myPhi-phiCen)
         if (dZ**2 + dS**2 .gt. zCirOut**2)   pathLength = 0.0
         if (dZ**2 + dS**2 .lt. zCirIn**2)    pathLength = 0.0
      elseif (myZ-zCen .gt. zTrack) then
         dZ = myZ - zCen - zTrack
         dS = rCylIn * abs(myPhi-phiCen)
         if (dZ**2 + dS**2 .gt. zCirOut**2)   pathLength = 0.0
         if (dZ**2 + dS**2 .lt. zCirIn**2)    pathLength = 0.0
      endif

c-g9b   Symmetry implies pathlength through cylinder is symmetric
c-g9b   MUST COUNT IF PIERCES COIL TWICE, THOUGH
c-g9b   if (ierr .eq. 2)
c-g9b      dist(13) = pathLength * 2.0

      RETURN
      end
