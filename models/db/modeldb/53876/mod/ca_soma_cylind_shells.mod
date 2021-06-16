COMMENT
This file, ca_soma_cylind_shells.mod, for Quadroni and Knopfel 1994, was modified from
cadifus.mod from Chapter 9 Hines and Carnevale NEURON
This file is identical to the other ca_*.mod files except has
the number of annulli, Nannuli (number of compartments) set to 25.
ENDCOMMENT
: Calcium ion accumulation with radial and longitudinal diffusion
NEURON {
SUFFIX ca_soma_cyl
USEION ca READ cai, ica WRITE cai
GLOBAL vrat : vrat must be GLOBAL --see INITIAL block
: however B which in cadifus.mod was called TotalBuffer may be and is here RANGE
RANGE K2f_ex, K2f_ATPase, B
RANGE i_Na_Ca_ex, i_ATPase, I	: Na Ca exchanger can pump one Ca out for 3 Na's 
				: that enter the cell.
				: the ATPase current pumps Ca out of the cell by
				: using up ATP
}
DEFINE Nannuli 50 : must be >=2 (i.e. at least shell and core)
UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
	(mA) = (milliamp)
	(mV) = (millivolt)
	FARADAY = (faraday) (10000 coulomb)
	PI = (pi) (1)
}

PARAMETER {
	DCa = 0.6 (um2/ms)
	: Ca buffer reaction rates :
	k1buf = 30  (/mM ms)  : these rates from Sala and Hernandez-Cruz 1990
	k2buf = 0.03 (/ms) : and are labeled f and b in 
		   : Quadroni and Knopfel 1994 p. 1916
	B = 0.025 (mM) : this is [B] in Quadroni and Knopfel 94 Table 4
	: compare above with these rates from
	: k1buf = 100 (/mM ms) : Yamada et al. 1989
	: k2buf = 0.1 (/ms)
	: B = 0.003 (mM)

	K2f_ex = 0 : should be 2*FARADAY*3.3e-13 (mA/cm2mM2): Quadroni and Knopfel 94 
	: K2f_ex = 5e-6 (mA/cm2mM2): a value from Lytton and Sejnowski's (LS '91) 5e-6
			: (Soma value only - dendrites are different)
	cao = 2 (mM) : [Ca]_outside is set constant because it
		   :	doesn't change in this simulation.
	E_1  = 0.01315 (/mV) : Quadroni and Knopfel 94
	E_2 = 0.0255 (/mV)   :  "	
	nai = 7.6 (mM)	     :  "	: (LS's '91 = 10  mM)
	nao = 152 (mM)	     :  "	: (LS's '91 = 140 mM)
	K2f_ATPase = 0 : 2*FARADAY*2.65e-9 (mA/cm2)  : type A soma only : dendrites and type B cell different
	f_ATPase = 100 (/mM ms)	: simply called f for forward rate in Quadroni Knopfel
	b_ATPase = 0.005 (/ms)	: 1994 - this one is just called b for backward
}

ASSIGNED {
	v (mV)
	diam (um)
	ica (mA/cm2)
	I (mA/cm2)	: same as ica but used for diagnostic purposes
	i_Na_Ca_ex (mA/cm2)
	i_ATPase (mA/cm2)
	cai (mM)
	vrat[Nannuli] (1) : dimensionless
	: numeric value of vrat[i] equals the volume
	: of annulus iof a 1um diameter cylinder
	: multiply by diam^2 to get volume per um length
	Kd (/mM)
	B0 (mM)
}
STATE {
	: ca[0] is equivalent to cai
	: ca[] are very small, so specify absolute tolerance
	ca[Nannuli] (mM) <1e10>
	CaBuffer[Nannuli] (mM)
	Buffer[Nannuli] (mM)
	n (1)
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	i_Na_Ca_ex = -K2f_ex * (nai^3 * cao * exp(E_1 * v) - nao^3 * cai * exp(-E_2*v))
	i_ATPase = K2f_ATPase * n
	I= ica	: diagnostic purposes only
	SOLVE state METHOD sparse
}

DERIVATIVE states {
	: compute state variable n at present v and t
	n' = f_ATPase * cai * (1 - n) - b_ATPase * n
}

LOCAL factors_done
INITIAL {
	if (factors_done == 0) { : flag becomes 1 in the first segment
		factors_done = 1 : all subsequent segments will have
		factors() : vrat = 0 unless vrat is GLOBAL
	}

	n = f_ATPase * cai / (f_ATPase * cai + b_ATPase)

	Kd = k1buf/k2buf
	B0 = B/(1 + Kd*cai)
	FROM i=0 TO Nannuli-1 {
		cai = 5e-5	: arbitrary initialization 
		ca[i] = cai
		Buffer[i] = B0
		CaBuffer[i] = B - B0
	}
}
LOCAL frat[Nannuli] : scales the rate constants for model geometry
PROCEDURE factors() {
	LOCAL r, dr2
	r = 1/2 : starts at edge (half diam)
	dr2 = r/(Nannuli-1)/2 : full thickness of outermost annulus,
	: half thickness of all other annuli
	vrat[0] = 0
	frat[0] = 2*r
	FROM i=0 TO Nannuli-2 {
		vrat[i] = vrat[i] + PI*(r-dr2/2)*4*dr2 : interior half
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2) : outer radius of annulus
		: div by distance between centers
		r = r - dr2
		vrat[i+1] = PI*(r+dr2/2)*2*dr2 : outer half of annulus
	}
}
LOCAL dsq, dsqvol : can't define local variable in KINETIC block
: or use in COMPARTMENT statement
KINETIC state {
	COMPARTMENT i, diam*diam*vrat[i] {ca CaBuffer Buffer}
	: LONGITUDINAL_DIFFUSION i, DCa*diam*diam*vrat[i] {ca}
	~ ca[0] << ( (-ica - i_Na_Ca_ex - i_ATPase)*PI*diam/(2*FARADAY)) 
			: ica is Ca current from hva or lva or both
			: i_Na_Ca_ex is the current from Na-Ca exhanger
			: i_ATPase is the current from Ca-ATPase
	FROM i=0 TO Nannuli-2 {
		~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
	}
	dsq = diam*diam
	FROM i=0 TO Nannuli-1 {
		dsqvol = dsq*vrat[i]
		~ ca[i] + Buffer[i] <-> CaBuffer[i] (k1buf*dsqvol, k2buf*dsqvol)
	}
	cai = ca[0]
}
