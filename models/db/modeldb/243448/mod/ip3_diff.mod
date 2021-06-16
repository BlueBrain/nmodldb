TITLE ip3 dyanmics for bladder small DRG neuron soma model

: Adapted from Fink et al., 2000
: Author: Darshan Mandge (darshanmandge@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

:For details refer: 
:A biophysically detailed computational model of bladder small DRG neuron soma 
:Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

NEURON {
	 SUFFIX ip3dif
 	 USEION ip3 READ ip3i WRITE ip3i VALENCE 1
  	 GLOBAL vol, DIP3, ip3i0, kdegr
     RANGE ip3i
	 THREADSAFE
}

DEFINE NANN 12:  :This needs to be changed if NANN in cadyn.mod is changed.

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
  	kdegr = 0.14e-3 (/ms)  : degredation rate Fink et al.,2000
  	DIP3 = 0.283(um2/ms)
  	ip3i0 = 0.16e-3 (mM)   : [IP3]0  initial and resting ip3i conc
}


ASSIGNED {
  	diam      (um)
 	ip3i      (mM)
  	vol[NANN]  		: numeric value of vol[i] equals the volume 
					: of annulus i of a 1um diameter cylinder
					: multiply by diam^2 to get volume per um length
}

STATE {
  	ip3[NANN]       (mM) <1e-6>
}

LOCAL factors_done

BREAKPOINT { SOLVE state METHOD sparse}


INITIAL {
   	if (factors_done == 0) {   : flag becomes 1 in the first segment
     	factors_done = 1       :   all subsequent segments will have
      	factors()              :   vol = 0 unless vol is GLOBAL
       }

  	ip3i = ip3i0
  	FROM i=0 TO NANN-1 {
    	ip3[i] = ip3i
	}
}


LOCAL frat[NANN]  : scales the rate constants for model geometry

PROCEDURE factors() {
  	LOCAL r, dr2
  	r = 1/2                : starts at edge (half diam)
  	dr2 = r/(NANN-1)/2     : full thickness of outermost annulus,
						   : half thickness of all other annuli
  	vol[0] = 0
  	frat[0] = 2*r
 	 FROM i=0 TO NANN-2 {
    		vol[i] = vol[i] + PI*(r-dr2/2)*2*dr2  : interior half
    		r = r - dr2
    		frat[i+1] = 2*PI*r/(2*dr2)  : outer radius of annulus
										: div by distance between centers
   		 r = r - dr2
    		vol[i+1] = PI*(r+dr2/2)*2*dr2  : outer half of annulus
  	}
}

LOCAL dsq, dsqvol  : can't define local variable in KINETIC block

KINETIC state {
  	COMPARTMENT i, diam*diam*vol[i]*0.81 {ip3 ip3i0}  : cytoplasmic volume is 0.81 of total cytoplasmic volume in each shell
	
	dsq = diam*diam
	
  	FROM i=0 TO NANN-2 {
   		 ~ ip3[i] <-> ip3[i+1]  (DIP3*frat[i+1], DIP3*frat[i+1])
  	}

  	
  	FROM i=0 TO NANN-1 {
    		dsqvol = dsq*vol[i]*0.81 						 : cytoplasmic volume is 0.81 of the total cell volume
    		~ ip3[i] <-> ip3i0  (kdegr*dsqvol, kdegr*dsqvol) : ip3i0 is inactive ip3
  	}

  	ip3i = ip3[0]
}