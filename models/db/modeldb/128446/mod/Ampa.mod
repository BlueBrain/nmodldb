TITLE AMPA

COMMENT
ENDCOMMENT

NEURON {
	POINT_PROCESS Ampa

	NONSPECIFIC_CURRENT i
	RANGE Cdur, Erev, g, gmax, kB
	RANGE r1FIX, r6FIX, r1, r2, r5, r6
	RANGE tau_1, tau_rec, tau_facil, U	 
	RANGE T, Tmax, Trelease		
	RANGE M, Diff, R, lamd
	RANGE tspike, PRE 
	RANGE NTdiffusion 	
	RANGE xview,yview,zview,Pview
}

UNITS {
	(nA) 	= (nanoamp)
	(mV) 	= (millivolt)
	(umho) 	= (micromho)
	(mM) 	= (milli/liter)
	(pS)	= (picosiemens)
	(nS) 	= (nanosiemens)
	(um) 	= (micrometer)
	PI 	= (pi)		(1)
}

PARAMETER {
	: postsynaptic parameters
	gmax		= 1200  (pS)		
	Cdur		= 0.3	(ms)	
	Erev		= 0	(mV)
	kB		= 0.44	(mM)
		 
	r1FIX		= 5.4	(/ms/mM) 
	r6FIX		= 1.12	(/ms/mM)		
	r2		= 0.82	(/ms)		
	r5		= 0.013	(/ms)		 
	
	: presynaptic parameters
	tau_1 		= 3 (ms) 	< 1e-9, 1e9 >
	tau_rec 	= 35.1 (ms) 	< 1e-9, 1e9 > 	
	tau_facil 	= 10.8 (ms) 	< 0, 1e9 > 	

	U 		= 0.416 (1) 	< 0, 1 >
	u0 		= 0 (1) 	< 0, 1 >	 
	Tmax		= 1  (mM)
	
	: Diffusion			
	M		= 21500				 
	R		= 1.033 (um)
	Diff		= 0.223 (um2/ms)
	lamd		= 20 (nm)
}


ASSIGNED {
	v		(mV)		 
	i 		(nA)		 
	g 		(pS)		 
	T		(mM)
	
	r1		(/ms)
	r6		(/ms)

	Trelease	(mM)
	tspike[50]	(ms)
	x 
	tsyn		(ms)
	PRE[50]
	
	Mres		(mM)	
	NTdiffusion	(mM)
	numpulses
	
	xview   
	yview  
	zview
	Pview   
}

STATE {	
	C
	O
	D
}

INITIAL {
	C=1
	O=0
	D=0
	
	T=0 (mM)
	Trelease=0 (mM)
	tspike[0]=1e12	(ms)

	Mres = ( 1e3 * 1e15 / 6.022e23 * M )   : (M) to (mM) so 1e3, 1um^3=1dm^3*1e-15 so 1e15   
	numpulses = 0
	
	xview = 1 
	yview = 0 
	zview = 0 
	Pview = 0	
}

FUNCTION NTdiffWave(){
	LOCAL ijk,t0
	: sums up diffusion contributes
	NTdiffusion=0
	FROM ijk=1 TO numpulses{
		t0=tspike[ijk-1]
		if(t>t0){		
			NTdiffusion=NTdiffusion+PRE[ijk-1]*Mres*exp(-R*R/(4*Diff*(t-t0)))/(4*PI*Diff*((1e-3)*lamd)*(t-t0))	
		}
	}					
	NTdiffWave=NTdiffusion
}

BREAKPOINT {	
	Trelease = T + NTdiffWave()  
	SOLVE kstates METHOD sparse
	g = gmax * O
	i = (1e-6) * g * (v-Erev) 
}

KINETIC kstates {
	r1 = r1FIX * Trelease^2 / (Trelease + kB)^2
	r6 = r6FIX * Trelease^2 / (Trelease + kB)^2
	~ C  <-> O (r1,r2)
	~ D  <-> C (r5,r6)
	CONSERVE C+O+D = 1
}

NET_RECEIVE(weight, on, nspike, t0 (ms),y, z, u, tsyn (ms)) {
	INITIAL {
		y = 0
		z = 0
		u = u0
		tsyn = t
		nspike = 1
	}
  	if (flag == 0) { 
		: presynaptic modulation
		nspike = nspike + 1
		if (!on) {
			t0 = t
			on = 1		
					
			z = z*exp( - (t - tsyn) / tau_rec )	
			z = z + ( y*(exp(-(t - tsyn)/tau_1) - exp(-(t - tsyn)/tau_rec))/((tau_1/tau_rec)-1) )  
			y = y*exp(-(t - tsyn)/tau_1)			
			x = 1-y-z
			
			if (tau_facil > 0) { 
				u = u*exp(-(t - tsyn)/tau_facil)
				u = u + U * ( 1 - u )							
			} else { u = U }
			y = y + x * u

			xview = x	 
			yview = y  
			Pview = u

			T = Tmax * y			
			PRE[numpulses] = y	 
			tspike[numpulses] = t
			numpulses = numpulses + 1
			tsyn = t
						
		}
		net_send(Cdur, nspike)	 
   	}
	if (flag == nspike) { 
			t0 = t
			T = 0
			on = 0
	}
}

