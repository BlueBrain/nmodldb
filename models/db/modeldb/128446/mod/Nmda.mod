TITLE NMDA  

COMMENT
ENDCOMMENT

NEURON {
	POINT_PROCESS Nmda
	
	NONSPECIFIC_CURRENT i
	RANGE Rb, Ru, Rd, Rr, Ro, Rc, rb
	RANGE g, gmax, Erev 	 			
	RANGE MgBlock,v0_block,k_block
	RANGE Trelease  
	RANGE tau_1, tau_rec, tau_facil, U	 
	RANGE R, M, Diff
}

UNITS {
	(nA) 	= (nanoamp)
	(mV) 	= (millivolt)
	(pS) 	= (picosiemens)
	(umho) 	= (micromho)
	(mM) 	= (milli/liter)
	(uM) 	= (micro/liter)
	PI	= (pi)		(1)
}

PARAMETER {
	: postsynaptic parameters
	gmax		= 1.88e4  	(pS)
	Erev		= 0		(mV)	 
	v0_block 	= -20 		(mV)	 
	k_block 	= 13		(mV)
	
	Rb		=  5		(/ms/mM)  	: binding  
	Ru		=  0.1		(/ms)		: unbinding
	Rd		=  12e-4  	(/ms)		: desensitization
	Rr		=  9e-3		(/ms)		: resensitization 
	Ro		=  3e-2 	(/ms)		: opening
	Rc		=  0.966	(/ms)		: closing	
 
	: presynaptic parameters
	tau_1 		= 3 (ms) 	< 1e-9, 1e9 >
	tau_rec 	= 35.1 (ms) 	< 1e-9, 1e9 > 	
	tau_facil 	= 10.8 (ms) 	< 0, 1e9 > 	

	U 		= 0.416 (1) 	< 0, 1 >
	u0 		= 0 (1) 	< 0, 1 >	
	
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
	rb		(/ms)    
	MgBlock
	
	Trelease	(mM)
	tspike[50]	(ms)
	x 
	tsyn		(ms)
	PRE[50]
	
	Mres		(mM)	
	NTdiffusion	(mM)
	numpulses
}

STATE {
	C0		: unbound
	C1		: single bound
	C2		: double bound
	D		: desensitized
	O		: open
}

INITIAL {
	rates(v)
	C0=1
	C1=0
	C2=0
	O=0
	D=0
		 
	Trelease=0 (mM)
	tspike[0]=1e12	(ms)

	Mres = ( 1e3 * 1e15 / 6.022e23 * M )   : (M) to (mM) so 1e3, 1um^3=1dm^3*1e-15 so 1e15   
	numpulses = 0
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
	Trelease = NTdiffWave()  
	SOLVE kstates METHOD sparse
	g = gmax * O
	i = (1e-6) * g * (v-Erev) * MgBlock 
}

KINETIC kstates {	
	rb = Rb * Trelease 
	~ C0 <-> C1	(rb,Ru) 	 
	~ C1 <-> C2	(rb,Ru)		 
	~ C2 <-> D	(Rd,Rr)
	~ C2 <-> O	(Ro,Rc)
	CONSERVE C0+C1+C2+D+O = 1
}

PROCEDURE rates(v(mV)) {
	: update the tables 
	TABLE MgBlock DEPEND v0_block,k_block FROM -120 TO 30 WITH 150
	MgBlock = 1 / ( 1 + exp ( - ( v - v0_block ) / k_block ) )
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
			
			PRE[numpulses] = y	 
			tspike[numpulses] = t
			numpulses = numpulses + 1
			tsyn = t	
		}
		on = 0
   	}
}
