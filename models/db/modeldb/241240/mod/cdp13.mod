
TITLE Calcium ion accumulation with longitudinal and radial diffusion buffering and pumping

COMMENT
Modified from examples/nrniv/nmodl/cabpump.mod to match The         
  NEURON Book write-up and to include more comments
See examples/nrniv/nmodl/capump.mod for the pump by itself
PROCEDURE factors_cdp() sets up the scale factors 
  needed to model radial diffusion.  These scale factors do not 
  have to be recomputed when diam or DCa is changed
The amount of calcium in an annulus is ca[i]*diam^2*vrat[i] 
  with ca[0] being the 2nd order correct concentration at the  
  exact edge and ca[Nannuli-1] being the concentration at the  
  exact center.
Buffer concentration and rates are based on Yamada et al. 1989
  model of bullfrog sympathetic ganglion cell.
ENDCOMMENT

NEURON {
	SUFFIX cdp13
	USEION ca READ cao, cai, ica WRITE cai, ica
		: need to read cai to initialize buffer and
		: calcium in the compartments
		: cao and ica added for pump.  Writes to ica
	GLOBAL vrat
      		: vrat must be GLOBAL--see 
      		: INITIAL block, but TotalBuffer and TotalPump may be RANGE
		: if it varies among compartments
      RANGE cai0, TotalBuffer, TotalPump, ica_pmp   
		: pump current is ica_pmp
		: TotalPump is density of pump sites
		: Make it RANGE if it varies among sections	
      	: pump current can be referenced as ica_pmp_cdp
}

DEFINE Nannuli  13	: number of concentric shells  
			: must be >= 2, at least shell and core
UNITS {
	(molar) =	(1/liter)
	(mM) =	(millimolar)
	(um) =	(micron)
	(mA) =	(milliamp)
	(mol) =	(1)  : pump density is in mol/cm2
	FARADAY =	(faraday)	(10000 coulomb)
	PI = (pi)	(1)
}

PARAMETER {
	DCa = 0.223	(um2/ms)
	: to change rate of buffering without disturbing equilibrium
	: multiply the following two by the same factor
	k1buf	= 100			(/mM-ms)
	k2buf	= 0.1			(/ms)
	TotalBuffer = 0.03	(mM)
	cai0 = 50e-6 (mM)  : Requires explicit use in INITIAL block
	k1	=	1		(/mM-ms)
	k2	=	0.005	(/ms)
	k3	=	1		(/ms)
	k4	=	0.005	(/mM-ms)
	: rates mean 50 nM is equilibrium
	: to eliminate the pump, set TotalPump to 0 in hoc
	TotalPump	=	3e-13	(mol/cm2)
}

ASSIGNED {
	diam		(um)
	ica		(mA/cm2)
	cai		(mM)
	cao		(mM)		: here cao is assumed constant
	ica_pmp	(mA/cm2)	
	ica_pmp_last	(mA/cm2) : used so pump not counted twice
	parea	(um)		: surface area for pump

	vrat[Nannuli]	(1)	: dimensionless.  Value of vrat[i] 
		:  is the volume of annulus i of a 1 um diameter cyl
		:  per unit length
		:  gets extra um2 when multiplied by diam^2 to give vol 
	Kd		(/mM)
	B0		(mM)	: initial free buffer
}

CONSTANT {volo = 1e+10 (um2) }  : need an extracellular volume
	: 1 liter per um of length.  Actual value not relevant
	: since cao is constant

STATE {
	: ca[0] is equivalent to cai
	: ca[] are very small, so specify absolute tolerance
	ca[Nannuli]		(mM) <1e-7>	
	CaBuffer[Nannuli]	(mM) <1e-5>
	Buffer[Nannuli]	(mM)  <1e-5>
	pump				(mol/cm2) <1e-15>	: free pump sites
	pumpca			(mol/cm2) <1e-15>	: pump sites with ca bound
}

BREAKPOINT {
	SOLVE state METHOD sparse
	ica_pmp_last = ica_pmp	: used so pump not counted twice
	ica = ica_pmp			: pump current needed for cai
						: calculation
}

LOCAL factors_done	: like a static variable in C

INITIAL {
	if (factors_done == 0) {   : flag becomes 1 in the first seg
		factors_done = 1	   : all subsequent segs will have
		factors()		   : vrat =0 unless vrat is GLOBAL
					   : vrat and frat have to be calculated
					   : only once
	}

	cai = cai0	: cai is set to cai0_ca_ion by default.  
      : This overrides
	Kd = k1buf/k2buf			
	B0 = TotalBuffer/(1 + Kd*cai)   :  various initializations

	FROM i=0 TO Nannuli-1 {
		ca[i] = cai
		Buffer[i] = B0
		CaBuffer[i] = TotalBuffer - B0
	}
	parea = PI*diam	: defines pump area per unit length
	pump = TotalPump/(1 + (cai*k1/k2))	: initial free sites
	pumpca = TotalPump -  pump	: initial ca bound sites
}

COMMENT
Note, the above initializations may not give a true "rest"
  If Ca currents are included, as these may affect cai.  It 
  may be necessary to do an "initialization run" and do 
  SaveState (see The NEURON Book Chap 8.4) 
factors() sets up factors needed for radial diffusion 
  modeled by Nannuli concentric compartments.
The outermost shell is half as thick as the other shells 
  so the concentration is spatially second order correct 
  at the surface of the cell.
The radius of the cylindrical core 
  equals the thickness of the outermost shell.
The intervening Nannuli-2 shells each have thickness = 
  r/(Nannuli-1).  (Nannuli must be >= 2).

ca[0] is at the edge of the cell, ca[Nannuli-1] is at the 
  center of the cell,  and ca[i] for 0 < i < Nannuli-1 is 
  midway through the thickness of each annulus.
ENDCOMMENT

LOCAL frat[Nannuli]	: scales the rate constants for model
				: geometry.  Local since applies to all
				: segments that have a caldif, but not
				: of general interest to be global

PROCEDURE factors() {
	LOCAL r, dr2
	r = 1/2			:radius starts at edge (half diam)
	dr2 = r/(Nannuli-1)/2	:full thickness of outermost annulus
					:half thickness of all others
	vrat[0] = 0
	frat[0] = 2*r
	FROM i=0 TO Nannuli-2 {
		vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2	:interior half
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)	:exterior edge of annulus
					: divided by distance between centers
					: frat[i+1] equals A(i+1)/(delta-r)
					: where A is the shell surface area
		r = r - dr2
		vrat[i+1] = PI*(r+dr2/2)*2*dr2	
      :outer half of annulus
	}
}
LOCAL dsq, dsqvol	: can't define local variable in KINETIC
      : block or use in COMPARTMENT

KINETIC state {
	COMPARTMENT i, diam*diam*vrat[i] {ca CaBuffer Buffer}
		: specifies shell volumes for state variables
		: for all compartments.  i is the index
		: converts states to mass from concentration
	COMPARTMENT (1e10)*parea {pump pumpca}
	COMPARTMENT volo {cao}
		: these COMPARTMENT statements and 1e10 needed
		: for dimensional consistency
	:LONGITUDINAL_DIFFUSION i, DCa*diam*diam*vrat[i] {ca}
		: diffusion into neighboring sections
		: i is the index, flux expression (DCa * cross-sec A),
		: variable (ca)
	: pump
	~ ca[0] + pump <-> pumpca (k1*parea*(1e10), k2*parea*(1e10))
	~ pumpca <-> pump + cao (k3*parea*(1e10), k4*parea*(1e10))
	CONSERVE pump + pumpca = TotalPump * parea * (1e10)
		: ensures conservation.  Helpful for accuracy
	ica_pmp = 2*FARADAY*(f_flux - b_flux)/parea
		: after each reaction the forward and backward fluxes
		: are assigned to f_flux and b_flux automatically
		: f_flux is ca that is pumped out

      ~ ca[0] << (-(ica - ica_pmp_last)*PI*diam/(2*FARADAY))
		: consider all currents except pump since pump
      : effect on ca[0] was considered a few statements
: above.  ica contains pump current ica_pmp from the
		: previous time step, so subtract ica_pmp_last here. 
		: (ica_pmp is the "new" value, but ca[0] must be
		: computed using the "old" value, i.e. ica_pmp_last)
      
		: uses kinetic syntax in place of an ode where ica
		: comes into the outer compartment from outside.
		: Note, mass/time here--not concentration

	FROM i=0 TO Nannuli-2 {	: radial diffusion
		~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
	}

	dsq = diam*diam
	FROM i=0 TO Nannuli-1 {
		dsqvol = dsq*vrat[i]
		~ ca[i] + Buffer[i] <-> CaBuffer[i] (k1buf*dsqvol, k2buf*dsqvol)
	}
		: the dsqvol is needed to convert from conc to mass
	cai = ca[0]
		: updates cai from ca[0] so eca can be computed 
      : by NEURON
}

COMMENT
	Be sure to use modlunit to do units checking
	When inserted into a section, 
	  ca_cadifus[], Buffer_cadifus[], and CaBuffer-cadifus[]
	  are available to the .hoc file and for plotting
	if total buffer or DCa are non-uniform in the cell, make 
        them RANGE variables
ENDCOMMENT
