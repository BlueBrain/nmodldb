NEURON {
SUFFIX dipole
: POINT_PROCESS Dipole
RANGE ri, ia, Q, ztan
POINTER pv

POINTER Qsum : for density. sums into Dipole at section position 1
: RANGE Qsum : for POINT_PROCESS. Gets additions from dipole
}

UNITS {
	(nA) = (nanoamp)
	(mV) =(millivolt)
	(Mohm) = (megaohm)
	(um) = (micrometer)
	(Am)= (amp meter)
	(fAm) = (femto amp meter)
}

ASSIGNED {
	ia (nA)
	ri (Mohm)
	pv (mV)
	v (mV)
	ztan (um)
	Q  (fAm)
	Qsum (fAm) :human dipole order of 10nAm
}

AFTER SOLVE {     	: solve for v's first then use them
	ia=(pv-v)/ri
	Q=ia*ztan
	Qsum = Qsum + Q
}
	
AFTER INITIAL {
	ia=(pv-v)/ri
	Q=ia*ztan
	Qsum = Qsum + Q
}

: following needed for POINT_PROCESS only but will work if also in SUFFIX
: BEFORE INITIAL {
:	Qsum = 0
: }
: BEFORE BREAKPOINT {
:	Qsum = 0
: }

