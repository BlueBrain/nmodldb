TITLE Calcium ion accumulation and diffusion with linear pump and pulse input.

: This model is designed solely for the Stockbridge and Moore (1984)
: J. Neurosci. 4:803-811 simulations.
: This model should not be used with other calcium channels since it
: does not use those as contributions to ion accumulation nor does it change
: use the system wide calcium ion.

: The internal coordinate system is set up in PROCEDURE coord_cadifus()
: and must be executed before computing the concentrations.
: The scale factors set up in this procedure do not have to be recomputed
: when diam or DFree are changed.
: The amount of calcium in an annulus is ca[i]*diam^2*vol[i] with
: ca[0] being the second order correct concentration at the exact edge
: and ca[NANN-1] being the concentration at the exact center

? interface
NEURON {
	SUFFIX cab
	NONSPECIFIC_CURRENT idummy : force use of breakpoint block
	GLOBAL vol, frat
	RANGE in
}

: The cylindrical coordinate system is set up so that second order correct calculations
: for ca[i] are at positions (i)*dx for i=0,20 and ca[21] is at the center
: of the cell.

DEFINE NANN  99

UNITS {
	(mV)	= (millivolt)
	(molar) = (1/liter)
	(uM)	= (micromolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(mol)	= (1)
	(pmol)	= (picomol)
	FARADAY = (faraday)	 (coulomb)
	PI	= (pi)		(1)
	R 	= (k-mole)	(joule/degC)
}

PARAMETER {
	DFree = .6	(um2/ms) <0,1e9>
	beta = 50		<0, 1e9>

	ke = 1e-3	(cm/s)	: pump rate

	dur = 1		(ms) : input pulse duration
	t2 = 1e9	(ms)
	amp = 1000	(pmol/cm2-s) : input pulse amplitude

	dx = 100	(angstrom) <1, 1e9> : note that 20*dx must be < diam/2
	ca0 = 0		(um)
}

ASSIGNED {
	celsius		(degC)
	diam		(um)
	vol[NANN]	(um2)	: volume per unit length
	frat[NANN] 	()	: area per unit length divided by distance between centers
	in		(pmol/cm2-s)
	idummy (mA/cm2)
}


STATE {
	ca[NANN]	(uM) <1e-6>
}

INITIAL {LOCAL total
	parms()
	FROM i=0 TO NANN-1 {
		ca[i] = ca0
	}
}

BREAKPOINT {
	SOLVE state METHOD sparse
	if (at_time(0)) {
		in = amp
	}
	if (at_time(dur)) {
		in = 0
	}
	if (at_time(t2)) {
		in = amp
	}
	if (at_time(t2+dur)) {
		in = 0
	}
	idummy = 0 : breakpoint blocks only used if there is a current
}


PROCEDURE coord() {
	LOCAL r, dr2
	: cylindrical coordinate system  with dx annuli thickness to
	: center of cell. Note however that the first annulus is half thickness
	: so that the concentration is second order correct spatially at
	: the membrane or exact edge of the cell.
	: note ca[0] is at edge of cell and location of ca[i] is i*dx
	: except ca[NANN-1] is at center of cell
	r = diam/2			:starts at edge (half diam)
	dr2 = dx/2*(1e-4)		:relative half thickness of annulus
	vol[0] = 0
	frat[0] = 0 : unused
	FROM i=0 TO NANN-2 {
		vol[i] = vol[i] + PI*(r-dr2/2)*2*dr2	: interior half
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)	:exterior edge of annulus
					: divided by distance between centers
		r = r - dr2
		vol[i+1] = PI*(r+dr2/2)*2*dr2	:outer half of annulus
	}
	: adjust for interior bulk compartment
	vol[NANN-1] = vol[NANN-1] + PI*r*r
}

KINETIC state {
	COMPARTMENT i, (1+beta)*vol[i]*1(um) {ca}
? kinetics
	~ ca[0] -> (ke*PI*diam*(10)*1(um)) : pump
	~ ca[0] << (in*PI*diam*(1e-2)*1(um)) : input pulse
	:diffusion
	FROM i=0 TO NANN-2 {
		~ ca[i] <-> ca[i+1] (DFree*frat[i+1]*1(um), DFree*frat[i+1]*1(um))
	}
}
	
PROCEDURE parms() {
	coord()
}
