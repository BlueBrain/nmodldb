COMMENT
This file, ca_soma_const_delta_thick_shells.mod, for Quadroni and Knopfel 1994, was modified from
the cylindrical diffusion model cadifus.mod from Chapter 9 Hines and Carnevale NEURON
This file is different than the other ca_[dist|prox].mod files in that it models spherical
shell diffusion rather than cylindrical shells.
ENDCOMMENT
: Calcium ion accumulation with radial diffusion in shells
NEURON {
SUFFIX ca_soma_const_drthick
USEION ca READ cai, ica WRITE cai
GLOBAL vrat : vrat must be GLOBAL --see INITIAL block
: however B which in cadifus.mod was called TotalBuffer may be and is here RANGE
GLOBAL factor	: diagnostic variable for concentration dependent studies
RANGE K2f_ex, K2f_ATPase, B
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
	k1buf = 30  (/mM ms)  : 30 (/mM ms) these rates from Sala and Hernandez-Cruz 1990
	k2buf = 0.03 (/ms) : and are labeled f and b in 
		   : Quadroni and Knopfel 1994 p. 1916
	B = 0.025 (mM) : this is [B] in Quadroni and Knopfel 94 Table 4
	: compare above with these rates from
	: k1buf = 100 (/mM ms) : Yamada et al. 1989
	: k2buf = 0.1 (/ms)
	: B = 0.003 (mM)

	K2f_ex = 0	: this should be 2*faraday*3.3e-13 (mA/cm2mM2): Quadroni and Knopfel 94 
	: K2f_ex = 5e-6 (mA/cm2mM2): a value from Lytton and Sejnowski's (LS '91) 5e-6
			: (Soma value only - dendrites are different)
	cao = 2 (mM) : [Ca]_outside is set constant to 2mM because it
		   :	doesn't change in this simulation.
	E_1  = 0.01315 (/mV) : Quadroni and Knopfel 94
	E_2 = 0.0255 (/mV)   :  "	
	nai = 7.6 (mM)	     :  "	: (LS's '91 = 10  mM)
	nao = 152 (mM)	     :  "	: (LS's '91 = 140 mM)
	K2f_ATPase = 0	: should be 2*faraday*2.65e-9 (mA/cm2)  : type A soma only : dendrites and type B cell different
	f_ATPase = 100 (/mM ms): 100 (/mM ms)	: simply called f for forward rate in 
						: Quadroni Knopfel
	b_ATPase = 0.005 (/ms)	: 1994 - this one is just called b for backward
	factor = 1 (1)	: diagnostic for cai concentration in rate eq
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
}
STATE {
	: ca[0] is equivalent to cai
	: ca[] are very small, so specify absolute tolerance
	ca[Nshells] (mM) <1e10>
	CaBuffer[Nshells] (mM)
	Buffer[Nshells] (mM)
	n (1)
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	i_Na_Ca_ex = -K2f_ex * (nai^3 * cao * exp(E_1 * v) - nao^3 * cai * exp(-E_2*v))
	i_ATPase = K2f_ATPase * n
	I= ica	: diagnostic purposes only
	: for (index=0;index<Nshells-1; index= index+1) { ca_[index] = ca]index] } : see if works
	SOLVE state METHOD sparse
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
		cai = 1e-5	: initialization
		ca[i] = cai
		Buffer[i] = B0
		CaBuffer[i] = B - B0
	}
}
LOCAL frat[Nshells] : scales the rate constants for model geometry
PROCEDURE factors() {
	LOCAL r, drthin, drthin2, drthick, drthick2	: the X2 variables are divided by 2 (unusual notation)
	r = diam/2 : starts at edge (half diam)
	drthick=r/Nthick	: thick shells divided up evenly
	drthick2 = drthick/2
	drthin = 2 * drthick / ( 2*Nthin - 1 )
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
		if (i==Nthin-1) {
			: special case - outer volumne is half a thin shell and inner volumner
			: is half a thick shell
			vrat[i] = vrat[i] + (4*PI/3) * (3*r^2*drthin2-3*r*drthin2^2+drthin2^3)
			r = r - drthin2
			frat[i+1] = 4*PI*r*r/(drthin2+drthick2)
			r = r - drthick2
			vrat[i+1] = (4*PI/3) * (3*r^2*drthick2 + 3*r*drthick2^2 + drthick2^3)
		} else {
			vrat[i] = vrat[i] + (4*PI/3) * (3*r^2*drthick2-3*r*drthick2^2+drthick2^3)
			r = r - drthick2
			frat[i+1] = 4*PI*r*r/drthick
			r = r-drthick2
			vrat[i+1] = (4*PI/3) * (3*r^2*drthick2 + 3*r*drthick2^2 + drthick2^3)
		}
	}
}

: LOCAL dsq, dsqvol : can't define local variable in KINETIC block
: or use in COMPARTMENT statement
KINETIC state {
	COMPARTMENT i, vrat[i] {ca CaBuffer Buffer}
	: note that LONGITUDINAL_DIFFUSION doesn't make sense for spherical shells
	~ ca[0] << ( (-ica - i_Na_Ca_ex - i_ATPase)*PI*diam*diam /(2*FARADAY)) 
			: ica is Ca efflux
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
