: Biophysical model of STDP
: Author: Mathilde Badoual, 2004

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
        POINT_PROCESS KINASE
	USEION ca READ cai
	
	RANGE an,bn,trans,am,bm,ah,bh

}
UNITS {
    (celsius) = (degC)        
	(nA) = (nanoamp)
    (mV) = (millivolt)
    (umho) = (micromho)
    (mM) = (milli/liter)

}

PARAMETER {

	trans
	cai		(mM)		: Ca concentration inside
	an	(1/mM/mM/mM/mM/ms)
	bn	(1/ms)
	am	(1/mM/ms)
	bm	(1/ms)
	ah	(1/mM/ms)
	bh	(1/ms)
 }


STATE {
        n nbegin m h mbegin hbegin:fraction of open sites by Ca
}




INITIAL {
	n=0.0
	m=0.0
	h=0.0
	nbegin=1.0
	mbegin=1.0
	hbegin=1.0
}

KINETIC scheme {
	~ nbegin + 4 cai <-> n (an, bn)
	~ hbegin + trans <-> h (ah, bh)
	~ mbegin + cai <-> m (am, bm)
}

BREAKPOINT {
	SOLVE scheme METHOD sparse
	if (n<0) {
 	  n=0.0
	}

	if (m<0) {
 	  m=0.0
	}
	if (h<0) {
 	  h=0.0
	}

	if (mbegin<0) {
 	  mbegin=0.0
	}
	if (nbegin<0) {
 	  nbegin=0.0
	}
	if (hbegin<0) {
 	  hbegin=0.0
	}
}




