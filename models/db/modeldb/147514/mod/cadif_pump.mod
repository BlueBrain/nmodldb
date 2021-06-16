
TITLE Calcium ion accumulation with longitudinal and radial diffusion

COMMENT
PROCEDURE factors_cadifus() sets up the scale factors 
needed to model radial diffusion.
These scale factors do not have to be recomputed
when diam or DFree is changed.
The amount of calcium in an annulus is ca[i]*diam^2*vol[i] 
with ca[0] being the 2nd order correct concentration at the exact edge
and ca[NANN-1] being the concentration at the exact center.
Buffer concentration and rates are based on Yamada et al. 1989
model of bullfrog sympathetic ganglion cell.
ENDCOMMENT

NEURON {
	SUFFIX cadp
	USEION ca READ cao, cai, ica WRITE cai, ica
	RANGE ica_pump
	GLOBAL vol, TotalBuffer, TotalPump
	RANGE cai0
	THREADSAFE
}

DEFINE NANN  4

UNITS {
	(molar) =	(1/liter)
	(mol) = (1)
	(mM) =	(millimolar)
	(um) =	(micron)
	(mA) =	(milliamp)
	FARADAY =	(faraday)	(10000 coulomb)
	PI = (pi)	(1)
}

PARAMETER {
	DCa = 0.6		(um2/ms)
	: to change rate of buffering without disturbing equilibrium
	: multiply the following two by the same factor
	k1buf	= 100			(/mM-ms)
	k2buf	= 0.1			(/ms)
	TotalBuffer = 0.003	(mM)
	cai0 = 50e-6 (mM)	: Requires explicit use in INITIAL block
	
	k1 = 1 (/mM-ms)
	k2 = 0.005(/ms)
	k3 = 1 (/ms)
	k4 = 0.005 (/mM-ms)
	: to eliminate pump, set Totalpump to 0 in hoc 
	TotalPump = 1e-11 (mol/cm2)
	
}

ASSIGNED {
	diam		(um)
	ica		(mA/cm2)
	cai		(mM)
	vol[NANN]	(1)	: gets extra um2 when multiplied by diam^2
	Kd		(/mM)
	B0		(mM)
	
	cao (mM)
	ica_pmp (mA/cm2)
	parea (um)
}

CONSTANT{ volo = 1e10 (um2) }

STATE {
	ca[NANN]		(mM) <1e-6>	: ca[0] is equivalent to cai
	CaBuffer[NANN]	(mM)
	Buffer[NANN]	(mM)
	
	pump (mol/cm2)
	pumpca (mol/cm2)
}

BREAKPOINT {
	SOLVE state METHOD sparse
	ica = ica_pmp
}

LOCAL factors_done

INITIAL {
	MUTEXLOCK
	if (factors_done == 0) {
		factors_done = 1
		factors()
	}
	MUTEXUNLOCK

	cai = cai0
	Kd = k1buf/k2buf
	B0 = TotalBuffer/(1 + Kd*cai)

	FROM i=0 TO NANN-1 {
		ca[i] = cai
		Buffer[i] = B0
		CaBuffer[i] = TotalBuffer - B0
	}
	
	parea = PI*diam
	pump = TotalPump/(1+(cai*k1/k2))
	pumpca = TotalPump - pump
	
}

COMMENT
factors() sets up factors needed for radial diffusion 
modeled by NANN concentric compartments.
The outermost shell is half as thick as the other shells 
so the concentration is spatially second order correct 
at the surface of the cell.
The radius of the cylindrical core 
equals the thickness of the outermost shell.
The intervening NANN-2 shells each have thickness = r/(NANN-1)
(NANN must be >= 2).

ca[0] is at the edge of the cell, 
ca[NANN-1] is at the center of the cell, 
and ca[i] for 0 < i < NANN-1 is 
midway through the thickness of each annulus.
ENDCOMMENT

LOCAL frat[NANN]

PROCEDURE factors() {
	LOCAL r, dr2
	r = 1/2			:starts at edge (half diam)
	dr2 = r/(NANN-1)/2	:half thickness of annulus
	vol[0] = 0
	frat[0] = 2*r
	FROM i=0 TO NANN-2 {
		vol[i] = vol[i] + PI*(r-dr2/2)*2*dr2	:interior half
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)	:exterior edge of annulus
					: divided by distance between centers
		r = r - dr2
		vol[i+1] = PI*(r+dr2/2)*2*dr2	:outer half of annulus
	}
}

LOCAL dsq, dsqvol	: can't define local variable in KINETIC block 
			: or use in COMPARTMENT

KINETIC state {
	COMPARTMENT i, diam*diam*vol[i] {ca CaBuffer Buffer}
	
	COMPARTMENT (1e10)*parea {pump pumpca}
	COMPARTMENT volo {cao}	
	LONGITUDINAL_DIFFUSION i, DCa*diam*diam*vol[i] {ca}
	
	:pump
	~ ca[0] + pump <-> pumpca (k1*parea*(1e10), k2*parea*(1e10))
	~ pumpca <-> pump + cao  (k3*parea*(1e10), k4*parea*(1e10))
	CONSERVE pump + pumpca = TotalPump *parea * (1e10)
	ica_pmp = 2*FARADAY*(f_flux - b_flux)/parea
	
	:all currents except pump
	
	~ ca[0] << (-(ica - ica_pmp)*PI*diam/(2*FARADAY))
	
	FROM i=0 TO NANN-2 {
		~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
	}
	dsq = diam*diam
	FROM i=0 TO NANN-1 {
		dsqvol = dsq*vol[i]
		~ ca[i] + Buffer[i] <-> CaBuffer[i] (k1buf*dsqvol, k2buf*dsqvol)
	}
	cai = ca[0]
}

