COMMENT
%W%                  %G%
ENDCOMMENT

TITLE Calcium ion accumulation and diffusion
: The internal coordinate system is set up in PROCEDURE coord_cadifus()
: and must be executed before computing the concentrations.
: The scale factors set up in this procedure do not have to be recomputed
: when diam or DFree are changed.
: The amount of calcium in an annulus is ca[i]*diam^2*vol[i] with
: ca[0] being the second order correct concentration at the exact edge
: and ca[NANN-1] being the concentration at the exact center

DEFINE NANN  12

INDEPENDENT {t FROM 0 TO 1 WITH 10 (ms)}

NEURON {
	SUFFIX cagrcdifus
	USEION ca READ cao, cai, ica WRITE cai
	RANGE vol, farea, pcdc, Buffer0
	GLOBAL CaRest
}


UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	FARADAY = (faraday)	 (coulomb)
	PI	= (pi) (1)
}

PARAMETER {
	DFree = .6	(um2/ms)
	cao		(mM)
	ica		(mA/cm2)
	k1buf		(/mM-ms)
	k2buf		(/ms)
	CaRest = 6e-5	(mM)
}

ASSIGNED {
	cai		(mM)
	vol[NANN]	(cm3)
        farea[NANN]     (cm2)
        pcdc[NANN]      (cm/ms)

}

STATE {
	ca[NANN]	(mM) : ca[0] is equivalent to cai
	CaBuffer[NANN]  (mM)
	Buffer[NANN]    (mM)
}

BREAKPOINT {
	SOLVE state METHOD sparse
}

LOCAL coord_done

 INITIAL {
	if (coord_done == 0) {
		coord_done = 1
		coord()
	}
	: note Buffer gets set to Buffer0 automatically
	: and CaBuffer gets set to 0 (Default value of CaBuffer0) as well
	cai = CaRest
	FROM i=0 TO NANN-1 {
		ca[i] = cai
                CaBuffer[i] = Buffer0/(1+k2buf/k1buf/ca[i])
                Buffer[i] = Buffer0 - CaBuffer[i]
	}
}


PROCEDURE coord() {
                  vol[0] = 3.0793e-11
                  vol[1] = 2.956e-11
                  vol[2] = 2.8354e-11
                  vol[3] = 2.7173e-11 
                  vol[4] = 2.6017e-11
                  vol[5] = 4.8665e-11
                  vol[6] = 6.4955e-11
                  vol[7] = 8.8489e-11
                  vol[8] = 6.6497e-11
                  vol[9] = 7.9587e-11
                  vol[10] = 2.9322e-11
                  vol[11] = 4.1888e-12
                  farea[0] = 3.094e-06
                  farea[1] = 3.0172e-06
                  farea[2] = 2.8953e-06
                  farea[3] = 2.7759e-06
                  farea[4] = 2.659e-06
                  farea[5] = 2.5447e-06
                  farea[6] = 2.3235e-06
                  farea[7] = 2.0106e-06
                  farea[8] = 1.5394e-06
                  farea[9] = 1.131e-06
                  farea[10] = 5.0265e-07
                  farea[11] = 1.2566e-07
                  pcdc[0] = 0
                  pcdc[1] = 1.8103e-09
                  pcdc[2] = 1.7372e-09
                  pcdc[3] = 1.6656e-09
                  pcdc[4] = 1.5954e-09
                  pcdc[5] = 1.0179e-09
                  pcdc[6] = 5.5764e-10
                  pcdc[7] = 3.0159e-10
                  pcdc[8] = 1.8473e-10
                  pcdc[9] = 9.0478e-11
                  pcdc[10] = 3.0159e-11
                  pcdc[11] = 5.0266e-12
}

LOCAL dsq, dsqvol : can't define local variable in KINETIC block or use
		:  in COMPARTMENT
KINETIC state {
	COMPARTMENT i, vol[i] {ca }
	~ ca[0] << (-ica*farea[0]/(2*FARADAY))
	FROM i=1 TO NANN-1 {
		~ ca[i-1] <-> ca[i] (pcdc[i],pcdc[i])
	}
	FROM i=0 TO NANN-1 {
		~ ca[i] + Buffer[i] <-> CaBuffer[i] (k1buf*vol[i],k2buf*vol[i])
	}
	cai = ca[0]
}

	
COMMENT
At this time, conductances (and channel states and currents are
calculated at the midpoint of a dt interval.  Membrane potential and
concentrations are calculated at the edges of a dt interval.  With
secondorder=2 everything turns out to be second order correct.
ENDCOMMENT

