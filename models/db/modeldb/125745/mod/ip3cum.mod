COMMENT
IP3 accumulation with breakdown and radial and longitudinal diffusion
patterned after cadifus
ENDCOMMENT

NEURON {
  SUFFIX ip3cum
  USEION ip3 READ ip3i, iip3 WRITE ip3i VALENCE 1
  GLOBAL vrat  : vrat must be GLOBAL--see INITIAL block
}

DEFINE Nannuli 4

UNITS {
  (molar) = (1/liter)
  (mM)    = (millimolar)
  (uM)    = (micromolar)
  (um)    = (micron)
  (mA)    = (milliamp)
  FARADAY = (faraday)  (coulomb)
  PI      = (pi)       (1)
}

PARAMETER {
  DIP3 = 0.283 (um2/ms)
  kdegr = 0.14e-3 (/ms)  : degredation rate
}

CONSTANT {
  ip3i0 = 0.16e-3 (mM)  : [IP3]0  initial and resting ip3i conc
}

ASSIGNED {
  diam      (um)
  iip3      (mA/cm2)
  ip3i      (mM)
  vrat[Nannuli]  : numeric value of vrat[i] equals the volume 
                 : of annulus i of a 1um diameter cylinder
                 : multiply by diam^2 to get volume per um length
}

STATE {
  : ip3[0] is equivalent to ip3i
  : ip3[] are very small, so specify absolute tolerance
  ip3[Nannuli]       (mM) <1e-6>
}

BREAKPOINT { SOLVE state METHOD sparse }

LOCAL factors_done

INITIAL {
   if (factors_done == 0) {  : flag becomes 1 in the first segment
      factors_done = 1       :   all subsequent segments will have
      factors()              :   vrat = 0 unless vrat is GLOBAL
   }

  ip3i = ip3i0
  FROM i=0 TO Nannuli-1 {
    ip3[i] = ip3i
  }
}

LOCAL frat[Nannuli]  : scales the rate constants for model geometry

PROCEDURE factors() {
  LOCAL r, dr2
  r = 1/2                : starts at edge (half diam)
  dr2 = r/(Nannuli-1)/2  : full thickness of outermost annulus,
                         : half thickness of all other annuli
  vrat[0] = 0
  frat[0] = 2*r
  FROM i=0 TO Nannuli-2 {
    vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2  : interior half
    r = r - dr2
    frat[i+1] = 2*PI*r/(2*dr2)  : outer radius of annulus
                                : div by distance between centers
    r = r - dr2
    vrat[i+1] = PI*(r+dr2/2)*2*dr2  : outer half of annulus
  }
}

LOCAL dsq, dsqvol  : can't define local variable in KINETIC block
                   :   or use in COMPARTMENT statement

KINETIC state {
  COMPARTMENT i, diam*diam*vrat[i] {ip3 ip3i0}
  LONGITUDINAL_DIFFUSION i, DIP3*diam*diam*vrat[i] {ip3}
  ~ ip3[0] << (-iip3*PI*diam*(1e4)/FARADAY)  : -iip3 is IP3 production at the cell membrane
  FROM i=0 TO Nannuli-2 {
    ~ ip3[i] <-> ip3[i+1]  (DIP3*frat[i+1], DIP3*frat[i+1])
  }
  dsq = diam*diam
  FROM i=0 TO Nannuli-1 {
    dsqvol = dsq*vrat[i]
    ~ ip3[i] <-> ip3i0  (kdegr*dsqvol, kdegr*dsqvol)
  }
  ip3i = ip3[0]
}
