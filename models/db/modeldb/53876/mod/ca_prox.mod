COMMENT
This file, ca_dist.mod or ca_prox.mod, for Quadroni and Knopfel 1994, was modified
from cadifus.mod from Chapter 9 Hines and Carnevale NEURON
These files ca_dist.mod, ca_prox.mod, only have differences in number
of annuli, Nannuli, for distal dendrites (5 prox vs 2 dist compartments)
and of course in the suffix names ( ca_dist and ca_prox).
Here ca diffusion in a cylinder is modeled along with an ATPase pump
(the real pump consumes ATP as it pumps Ca out of the cell) and a
Na/Ca exchanger (each real exchanger pumps 1 Ca ions out of the cell for
3 Na ions let into the cell)
ENDCOMMENT

: Calcium ion accumulation with radial (and uncommentable longitudinal) diffusion
NEURON {
SUFFIX ca_prox		:	 			     ***************************************
USEION ca READ cai, ica WRITE cai	: note that the membrane ica is approximated as
			: not being changed by the little ATPase pump and the Na/Ca exchangers
			: and therefore ica is not written to with WRITE by this mechanism
NONSPECIFIC_CURRENT ifake : used to move the calc. of the breakpt. to same t as other currents
GLOBAL vrat : vrat must be GLOBAL --see INITIAL block
: however B which in cadifus.mod was called TotalBuffer may be and is here RANGE
RANGE K2f_ex, K2f_ATPase, B
}

DEFINE Nannuli 5 : must be >=2 (i.e. at least shell and core) ***************************************

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
	k1buf = 30 (/mM ms)  : these rates from Sala and Hernandez-Cruz 1990
	k2buf = 0.03 (/ms) : and are labeled f and b in 
		   : Quadroni and Knopfel 1994 p. 1916
	B = 0.025 (mM) : this is [B] in Quadroni and Knopfel 94 Table 4
	: compare above with these rates from
	: k1buf = 100 (/mMms) : Yamada et al. 1989
	: k2buf = 0.1 (/ms)
	: B = 0.003 (mM)

	K2f_ex = 4.85e-13 (mA/cm2mM2): Quadroni and Knopfel 94 (distal dendrite A type cells)
		: this above value will be overwritten to give each section it's proper value.
		: int the Q+R 94 files this is done through bothcells4.ses
	cao = 2 (mM) : 2e-3 (M) : [Ca]_outside is set constant to 2mM because it
		   :	doesn't change in this simulation (thought to be approx true in real).
	E_1  = 0.01315 (/mV) : Quadroni and Knopfel 94
	E_2 = 0.0255 (/mV)   :  "	
	nai = 7.6 (mM)	     :  "
	nao = 152 (mM)	     :  "
	K2f_ATPase = 9.2625897e-06 (mA/cm2 mM4)  : 9.6e-11 (/umol ms cm2)  : type A distal dendrites only
	f_ATPase = 100 (/mM ms)	: simply called f for forward rate in Quadroni Knopfel
	b_ATPase = 0.005 (/ms)	: 1994 - this one is just called b for backward
	mM2M = 1e-3 (1)	: mM to M conversion for cai concentration in rate eq
}

ASSIGNED {
	v (mV)
	diam (um)
	ica (mA/cm2)
	i_Na_Ca_ex (mA/cm2)
	i_ATPase (mA/cm2)
	cai (mM)
	vrat[Nannuli] (1) : dimensionless
	: numeric value of vrat[i] equals the volume
	: of annulus iof a 1um diameter cylinder
	: multiply by diam^2 to get volume per um length
	Kd (/mM)
	B0 (mM)
	ifake (mA/cm2) : fake current used to make solver move
		       : calculation of breakpt. to same time as currents.
}
STATE {
	: ca[0] is equivalent to cai
	: ca[] are very small, so specify absolute tolerance
	ca[Nannuli] (mM) <1e-6>
	CaBuffer[Nannuli] (mM) <1e-6>
	Buffer[Nannuli] (mM) <1e-6>
	n (1)
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	i_Na_Ca_ex = -K2f_ex * (nai^3 * cao * exp(E_1 * v) - nao^3 * cai * exp(-E_2*v))
	i_ATPase = K2f_ATPase * n
	SOLVE state METHOD sparse
	ifake=0 : causes solver to execute the breakpoint block at the
:	same time as the calculation of the other currents
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
	B0 = B/(1 + Kd * cai )
	FROM i=0 TO Nannuli-1 {
		cai = 5e-5 :	initialization value of 50 uM
		ca[i] = cai	: keep stored values of cai in millimolar units
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
		vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2 : interior half
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
			: ica is Ca efflux from lva and hva mechanisms
			: i_Na_Ca_ex is the flux from Na-Ca exhanger
			: i_ATPase is the flux from Ca-ATPase
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
