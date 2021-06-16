COMMENT
: FORREST MD (2014) Two Compartment Model of the Cerebellar Purkinje Neuron
Longitudinal diffusion of sodium (no buffering)
(equivalent modified euler with standard method and
equivalent to diagonalized linear solver with CVODE )
ENDCOMMENT

NEURON {
	SUFFIX nadifl
	USEION na READ ina WRITE nai, ena
:        USEION Na READ iNa WRITE Nai VALENCE 1
	RANGE D, Nai, Total, neo, tau
}

UNITS {
	(mM) = (milli/liter)
	(um) = (micron)
	FARADAY = (faraday) (coulomb)
	PI = (pi) (1)
}

PARAMETER {
	D = .6 (um2/ms)
         k1buf = 0.01          : 1000                    : 
        k2buf =  0          : 1000                   :   
: the smaller these numbers the less equilibriated these are. 
: 0.001 gets a seperation of 3. though note that the seperation grows over time.  
: 0.0001 see a seperation of 10 after a while
       tau = 5000 (ms)
}

ASSIGNED {
	ina (milliamp/cm2)
	diam (um)
        iNa (milliamp/cm2)
        Total (mM)
        neo (milliamp/cm2)
        ena (mV)
}

STATE {
	nai (mM)
        Nai (mM)
}

INITIAL {
lates()
:	nai = 6
nai = 10
Nai = 10
Total = 20
neo = ina	
}

BREAKPOINT {
	SOLVE conc METHOD sparse
}

KINETIC conc {
        lates()
	COMPARTMENT PI*diam*diam/4 {nai}
	LONGITUDINAL_DIFFUSION D {nai}
~ nai << (-neo/(FARADAY)*PI*diam*(1e4))

:	~ nai << (-ina/(FARADAY)*PI*diam*(1e4))
:        ~ Nai << (-iNa/(FARADAY)*PI*diam*(1e4))
:        ~ nai <-> Nai (k1buf,k2buf)
: Total = nai + Nai + ina + iNa
: ~ nai << ((-ina/(FARADAY)*PI*diam*(1e4))/6)
: ~ nai << ((-ina/(FARADAY)*PI*diam*(1e4))+(10-nai)/0.3)
}


PROCEDURE lates() {
LAG ina BY tau
  neo = lag_ina_tau
if (ena < 70) {ena = 70}
if (nai < 10) {nai = 10}
}