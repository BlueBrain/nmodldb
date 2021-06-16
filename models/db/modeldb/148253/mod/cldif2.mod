COMMENT

Chloride accumulation and diffusion with chloride pump (Lineweaver-Burke equation) and chloride leak

Diffusion model is modified from Ca diffusion model in Hines & Carnevale: 
Expanding NEURON with NMODL, Neural Computation 12: 839-851, 2000 (Example 8)

ENDCOMMENT

NEURON {
	SUFFIX cldifus2
	USEION cl READ icl WRITE cli VALENCE -1
	USEION hco3 READ hco3i, hco3o VALENCE -1
	GLOBAL vrat		:vrat must be GLOBAL
	RANGE cli0, vmax, leak, Kd
}

DEFINE Nannuli 4

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
	(mA) = (milliamp)
	FARADAY = (faraday) (10000 coulomb)
	PI = (pi) (1)
}

PARAMETER {
	DCl = 2 (um2/ms) : Kuner & Augustine, Neuron 27: 447
	leak = .000588235 (mM/ms) 	:0.3 (mM/sec) 
	vmax = .005 (mM/ms)	:5 (mM/sec)
	Kd = 15 (mM)
	cli0 = 2 (mM)
	: Requires explicit use in INITIAL block in order 
	: to take precedence over the global cli0_cl_ion
	: Do not forget to initialize in hoc if different
	: from this default.
}

ASSIGNED {
	diam 	(um)
	icl 	(mA/cm2)
	cli 	(mM)
	hco3i	(mM)
	hco3o	(mM)
	vrat[Nannuli]	: numeric value of vrat[i] equals the volume
			: of annulus i of a 1um diameter cylinder
			: multiply by diam^2 to get volume per um length
}

STATE {
	: cl[0] is equivalent to cli
	: cl[] are very small, so specify absolute tolerance
	cl[Nannuli]	(mM) <1e-10>
}


BREAKPOINT { SOLVE state METHOD sparse}

LOCAL factors_done

INITIAL {
	if (factors_done == 0) {  	: flag becomes 1 in the first segment	
		factors_done = 1	: all subsequent segments will have
		factors()		: vrat = 0 unless vrat is GLOBAL
	}
	cli = cli0
	FROM i=0 TO Nannuli-1 {
		cl[i] = cli
	}
}

LOCAL frat[Nannuli]	: scales the rate constants for model geometry

PROCEDURE factors() {
	LOCAL r, dr2
	r = 1/2			: starts at edge (half diam), diam = 1, length = 1
	dr2 = r/(Nannuli-1)/2	: full thickness of outermost annulus,
				: half thickness of all other annuli
	vrat[0] = 0
	frat[0] = 2*r		: = diam
	FROM i=0 TO Nannuli-2 {
		vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2	: interior half
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)	: outer radius of annulus Ai+1/delta_r=2PI*r*1/delta_r
						: div by distance between centers 
		r = r - dr2
		vrat[i+1] = PI*(r+dr2/2)*2*dr2	: outer half of annulus
	}
}

KINETIC state {
	COMPARTMENT i, diam*diam*vrat[i] {cl}
	LONGITUDINAL_DIFFUSION i, DCl*diam*diam*vrat[i] {cl}
	~ cl[0] << ((icl*PI*diam/FARADAY) + (leak - vmax*(cl[0]/(Kd + cl[0])))*diam*diam*vrat[0]) : icl is Cl- influx 
	FROM i=0 TO Nannuli-2 {
		~ cl[i] <-> cl[i+1]	(DCl*frat[i+1], DCl*frat[i+1])
	}
	cli = cl[0]
}
