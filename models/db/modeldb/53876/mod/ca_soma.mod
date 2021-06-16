COMMENT
This file, ca_soma.mod, for Quadroni and Knopfel 1994, was modified from
the cylindrical diffusion model cadifus.mod from Chapter 9 Hines and Carnevale NEURON
This file is different than the other ca_[dist|prox].mod files in that it models spherical
shell diffusion rather than cylindrical shells.
ENDCOMMENT
: Calcium ion accumulation with radial diffusion in shells
NEURON {
SUFFIX ca_soma
USEION ca READ cai, ica WRITE cai
NONSPECIFIC_CURRENT ifake : used to move the calc. of the breakpt. to same t as other currents
GLOBAL vrat : vrat must be GLOBAL --see INITIAL block
: however B which in cadifus.mod was called TotalBuffer may be and is here RANGE
GLOBAL factor	: diagnostic variable for concentration dependent studies
RANGE K2f_ex, K2f_ATPase, B	: Note capital K2f_'s are intensive qty's, i.e. (mA/cm2 mM4)
RANGE i_Na_Ca_ex, i_ATPase, I	: I is used to display ica_ values for this mechanism
				: Na Ca exchanger can pump one Ca out for 3 Na's
				: that enter the cell.
				: the ATPase current pumps Ca out of the cell by
				: using up ATP
}
DEFINE Nthin	5	: Nthin = # of thin spherical shells that are all within and take up
			: the outermost thick shell.  The thick shells are evenly distributed
			: throughout the sphere in this first model.  A future version closer
			: to the paper will model the distribution with each thick shell having
			: constant volume
DEFINE Nthick	21	: # of thick spherical shells including the last one which is divided by above

DEFINE Nshells	25	: ((Nthick-1) + Nthin)	: total number of shells

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

	K2f_ex =  6.3680304e-08 (mA/cm2 mM4) : 3.3e-13 (/umol ms mM4 cm2)*2*FARADAY : Q and K 94 
	: K2f_ex = 5e-6 (mA/cm2mM2): a value from Lytton and Sejnowski's (LS '91) 5e-6
			: (Soma value only - dendrites are different)
	cao = 2 (mM) : [Ca]_outside is set constant to 2 mM because it
		   :	is to a good approximation (assumed not to change much in real world).
	E_1  = 0.01315 (/mV) : Quadroni and Knopfel 94
	E_2 = 0.0255 (/mV)   :  "	
	nai = 7.6 (mM)	     :  "	: (LS's '91 = 10  mM)
	nao = 152 (mM)	     :  "	: (LS's '91 = 140 mM)
	K2f_ATPase =0.00051137214  (mA/cm2) :  2.65e-9   (/umol ms cm2)*2*FARADAY : type A soma only
							 : dendrites and type B cell different
	f_ATPase = 100 (/mM ms)	: simply called f for forward rate in Quadroni Knopfel
	b_ATPase = 0.005 (/ms)	: 1994 - this one is just called b for backward
	factor = 1 (1)	: diagnostic for cai concentration in rate eq
	mM2M = 1e-3 (1)	: mM to M conversion for cai concentration in rate eq
}

ASSIGNED {
	v (mV)
	diam (um)
	ica (mA/cm2)
	I (mA/cm2)	: same as ica but used for diagnostic purposes
	i_Na_Ca_ex (mA/cm2)
	i_ATPase (mA/cm2)
	cai (mM)
	vrat[Nshells] (um3)	: vrat[i] equals the volume of a spherical shell in cubic microns
	Kd (/mM)
	B0 (mM)
	ifake (mA/cm2) : fake current used to make solver move
		       : calculation of breakpt. to same time as currents.
}
STATE {
	: ca[0] is equivalent to cai
	: ca[] are very small, so specify absolute tolerance
	ca[Nshells] (mM) <1e-6>
	CaBuffer[Nshells] (mM) <1e-6>
	Buffer[Nshells] (mM) <1e-6>
	n (1)
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	SOLVE state METHOD sparse
	i_Na_Ca_ex = -K2f_ex * (nai^3 * cao * exp(E_1 * v) - nao^3 * cai * exp(-E_2*v))
	i_ATPase = K2f_ATPase * n
	I= ica	: diagnostic purposes only
	ifake=0 : causes solver to execute the breakpoint block at the
:	same time as the calculation of the other currents
}

DERIVATIVE states {
	: compute state variable n at present v and t
	n' = f_ATPase * cai * (1 - n) - b_ATPase * n
}

INITIAL {
	factors() : vrat = 0 unless vrat is GLOBAL
	n = f_ATPase * cai / (f_ATPase * cai + b_ATPase)

	Kd = k1buf/k2buf
	B0 = B/(1 + Kd*cai)
	FROM i=0 TO Nshells-1 {
		cai = 5e-5 :	initialization value of 50 uM in mM
		ca[i] = cai
		Buffer[i] = B0
		CaBuffer[i] = B - B0
	}
}
LOCAL frat[Nshells] : scales the rate constants for model geometry
LOCAL radii[Nthick] : these are the constant volume thick shells radii
LOCAL drthick[Nthick] : these are simply the delta between the above
	: note that r_0 = r, (1/21)(4PI/3)r^3=(4PI/3)(r_i^3-r_(i+1)^3)
	: from which we derive r_(i+1) = ( r_i^3 - (1/21) * r_(i+1)^3 )^(1/3)
LOCAL drthick2[Nthick]	: in our weird notation will be just drthick[]/2
LOCAL i_	: an index for drthick[]'s

PROCEDURE factors() {
	LOCAL r, drthin, drthin2	: the X2 variables are divided by 2 (unusual notation)
	r = diam/2 : starts at edge (half diam)
	radii[0] = r
	FROM i=1 TO Nthick-1 {
		radii[i] = ( radii[i-1]^3 - (1/21)*r^3 )^(1/3)
		drthick[i-1] = radii[i-1] - radii[i]
		drthick2[i-1]=drthick[i-1]/2
	}
: if drthick was constant:	drthick=r/Nthick	: thick shells divided up evenly
: if drthick was constant:	drthick2 = drthick/2
:
: to figure out the below note that 3/2 drthin's fit if Nthin=2, 5/2 drthin's fit if Nthin=3...
	drthin = 2 * (drthick[0]) / ( 2*Nthin - 1 )	: the outter thick shell divided up
	drthin2 = drthin/2 : full thickness of outermost shell
	: half thickness of all other shells
	vrat[0] = 0
	frat[0] = 2*r	: where does this get used?
	: first do outer thin shells
	FROM i=0 TO Nthin-2 {
		: whole outermost shell, otherwise outer half of inner shells
		vrat[i] = vrat[i] + (4*PI/3) * (3*r^2*drthin2-3*r*drthin2^2+drthin2^3)
		r = r - drthin2
		frat[i+1] = 4*PI*r*r/drthin
		r = r-drthin2
		vrat[i+1] = (4*PI/3) * (3*r^2*drthin2 + 3*r*drthin2^2 + drthin2^3)
	}
	: next do inner thick shells
	FROM i=Nthin-1 TO Nshells-2 {
		i_ = 1	: special index for drthick starts at second index, one, since
			: the zeroth index is where the thin shells are.
		if (i==Nthin-1) {
			: special case - outer volumne is half a thin shell and inner volumner
			: is half a thick shell
			vrat[i] = vrat[i] + (4*PI/3) * (3*r^2*drthin2-3*r*drthin2^2+drthin2^3)
			r = r - drthin2
			frat[i+1] = 4*PI*r*r/(drthin2+drthick2[i_])
			r = r - drthick2[i_]
			vrat[i+1] = (4*PI/3) * (3*r^2*drthick2[i_] + 3*r*drthick2[i_]^2 + drthick2[i_]^3)
		} else {
			vrat[i] = vrat[i] + (4*PI/3) * (3*r^2*drthick2[i_]-3*r*drthick2[i_]^2+drthick2[i_]^3)
			r = r - drthick2[i_]
			frat[i+1] = 4*PI*r*r/drthick[i_]
			r = r-drthick2[i_]
			vrat[i+1] = (4*PI/3) * (3*r^2*drthick2[i_] + 3*r*drthick2[i_]^2 + drthick2[i_]^3)
		}
		i_ = i_ + 1
	}
}

: LOCAL dsq, dsqvol : can't define local variable in KINETIC block
: or use in COMPARTMENT statement
KINETIC state {
	COMPARTMENT i, vrat[i] {ca CaBuffer Buffer}
	: note that LONGITUDINAL_DIFFUSION doesn't make sense for spherical shells
	~ ca[0] << ((-ica - i_Na_Ca_ex - i_ATPase)*PI*diam*diam /(2*FARADAY)) 
			: ica is Ca efflux from lva and hva mechanisms
			: i_Na_Ca_ex is the current from Na-Ca exhanger
			: i_ATPase is the current from Ca-ATPase
	FROM i=0 TO Nshells-2 {
		~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
	}
:	dsq = diam*diam
	FROM i=0 TO Nshells-1 {
:		dsqvol = dsq*vrat[i]
		~ ca[i] + Buffer[i] <-> CaBuffer[i] (k1buf*vrat[i], k2buf*vrat[i])
	}
	cai = ca[0]
}
