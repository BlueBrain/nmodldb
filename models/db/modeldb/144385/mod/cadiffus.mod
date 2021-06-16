TITLE Calcium ion accumulation with longitudinal and radial diffusion

COMMENT
PROCEDURE factors_cadiffus() sets up the scale factors 
needed to model radial diffusion.
These scale factors do not have to be recomputed
when diam is changed.
The amount of calcium in an annulus is ca[i]*diam^2*vol[i] 
with ca[0] being the 2nd order correct concentration at the exact edge
and ca[NANN-1] being the concentration at the exact center.

ENDCOMMENT

NEURON {
	SUFFIX cadiffus
	USEION ca READ cai, ica WRITE cai
	GLOBAL vrat
}

DEFINE Nannuli  4

UNITS {
	(molar) =	(1/liter)
	(mM) =	(millimolar)
	(um) =	(micron)
	(mA) =	(milliamp)
	FARADAY =	(faraday)	(10000 coulomb)
	PI = (pi)	(1)
}

PARAMETER {
	DCa = 	0.23		(um2/ms) 
}

ASSIGNED {
	diam	(um)
	ica		(mA/cm2)
	cai		(mM)
	vrat[Nannuli]		: numeric value of vrat[i] equals the volume
                        : of annulis i of a 1um diameter cylinder
                        : multiply by diam^2 to get volume per um length
	B0		(mM)
}

CONSTANT { volo = 1e10 (um2)}

STATE {
	ca[Nannuli]		(mM) <1e-6>	: ca[0] is equivalent to cai
}

BREAKPOINT {
	SOLVE state METHOD sparse
}

LOCAL factors_done

INITIAL {
	if (factors_done == 0) {
		factors_done = 1
		factors()
	}

	FROM i=0 TO Nannuli-1 {
		ca[i] = cai
	}
}


LOCAL frat[Nannuli]

PROCEDURE factors() {
	LOCAL r, dr2
	r = 1/2			:starts at edge (half diam)
	dr2 = r/(Nannuli-1)/2	:half thickness of annulus
	vrat[0] = 0
	frat[0] = 2*r
	FROM i=0 TO Nannuli-2 {
		vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2	:interior half
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)	:exterior edge of annulus
					: divided by distance between centers
		r = r - dr2
		vrat[i+1] = PI*(r+dr2/2)*2*dr2	:outer half of annulus
	}
}

LOCAL dsq, dsqvol	: can't define local variable in KINETIC block 
			: or use in COMPARTMENT

KINETIC state {
	COMPARTMENT i, diam*diam*vrat[i] {ca CaBuffer Buffer}
	LONGITUDINAL_DIFFUSION i, DCa*diam*diam*vrat[i] {ca}
	~ ca[0] << (-ica*PI*diam/(2*FARADAY))
	FROM i=0 TO Nannuli-2 {
		~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
	}
cai = ca[0]
}
