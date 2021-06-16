COMMENT

Multiple presynaptic spike generator with correlated outputs
------------------------------------------------------------

ALGORITHM

 This mechanism has been written to be able to use synapses in a single
 neuron receiving various types of presynaptic trains.  Several randomly
 spiking outputs are generated here, but correlated with each other 
 according to a given correlation coefficient (correl). 

 "Distributed generator"
 The algorithm generates N output variables with correlation according 
 to the following algorithm.  N2 independent Poisson-distributed variables 
 are generated and distributed among the N output variables.  If N2<N, there
 will be some redundancy in the N outputs and the correlation is a complex
 function of N and N2:
 If N2 = 1, correl=1 (all outputs are identical to the same rnd variable)
 If N2 = N, correl=0 (all outputs are independent rnd variables)
 If N2>1 and N2<N, the correlation is intermediate
 The program calculates N2 from the linear interpolation:
	N2 = N + (1-N) * sqrt(correl)

 The algorithm starts by generating N2 independent variables (R[i]), then 
 distribute these N2 variables among the N outputs using a coin-tossing 
 procedure.

INPUT PARAMETERS

 N	: number of random channels generated
 freq	: mean frequency of firing of each channel (Hz, nb spikes per second)
 correl	: value of the desired correlation
 refract: minimal period between spikes (ms)
 min_val: min value of presynaptic variable (mV)
 max_val: max value of presynaptic variable (mV)
 on	: on=1 the generator is on, on=0 it is interrupted (default=1)
	  on=2, a spike is forced for all outputs, then on is reset to 1
	  on=3, a spike is forced for outputs selected by the vector "sync"
 sync   : array of flags; when sinc[i] is set to 1, then channel i will
	  fire when on=3
 latency: latency at which spikes begin (ms; initialized to 0)
 shutoff: time at which the generator is shutoff (ms; initialized to -10000)


OUTPUT PARAMETERS

 x	: array of N elements, contain the values of the outputs
 ns	: array of N elements, spike count for each channel (reset by init)
 ls	: array of N elements, time of last spike (reset by init)
 spont	: spontaneous probability of firing, calculated at each dt
 ext	: external trigger, calculated at each dt
 prob	: probability of firing, calculated at each dt

PROCEDURES

 new_seed : sets the seed to the value passed as argument
 printvec : prints the vectors


EXAMPLE OF HOW TO USE

 access PRE		// presynaptic compartment is fake here

 objectvar pg
 pg = new corrGen(0.5)  // create random spike generator

 pg.N = 10		// number of output channels
 pg.freq	= 40	// mean frequency of firing of each channel (Hz)
 pg.correl = 0		// correlation between outputs
 pg.refract = 1		// refractory period for spikes (ms)
 pg.min_val = -70  	// min value of presynaptic variable (mV)
 pg.max_val = 40 	// max value of presynaptic variable (mV)
 pg.latency = 50	// time at which generator begins (ms)


Alain Destexhe, Laval University, 1995

Modif: June 98: added special case for correl=0 to accelerate

-------------------------------------------------------------------
ENDCOMMENT

DEFINE MAXCHANNELS 25000		: maximum number of outputs

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON	{ 
	POINT_PROCESS corrGen8
	RANGE N, freq, correl, refract, min_val, max_val
	RANGE on, latency, shutoff
	RANGE x, ns, ls, sync
	RANGE spont, prob, N2
	RANGE new_seed, printvec
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	dt	(ms)

	N	= 100		: number of outputs
	N2	= 1		: number of independent variables
	freq	= 40 (/s)	: mean frequency of firing of each channel
	correl	= 0		: correlation between outputs
	refract = 1 (ms)	: refractory period for spikes
	min_val	= -70 (mV)	: min value of presynaptic variable
	max_val	= 40 (mV)	: max value of presynaptic variable
	on	= 1		: logical on/off variable	
	latency	= 0 (ms)	: time at which spikes begin
	shutoff	= 1e6 (ms)	: shutoff time

}

ASSIGNED {
	x[MAXCHANNELS]	(mV)	: outputs
	R[MAXCHANNELS]		: independent variables
	ns[MAXCHANNELS]		: spike counters for each channel
	ls[MAXCHANNELS]		: time of last spike for each channel
	sync[MAXCHANNELS]	: flag for sync firing
	spont			: probability of spontaneous firing
	prob			: probability of firing
}
	
INITIAL { LOCAL i
	spont = (0.001) * freq * dt	: initiate spontaneous probability
	if(spont > 0.5) {
	   VERBATIM
	   printf("\n\nERROR in correlated random generator\n");
	   printf("firing probability is too high: spont=%g\n",(float)spont);
	   printf("decrease integration step or mean firing frequency\n");
	   exit(-1);
	   ENDVERBATIM
	}
	FROM i=0 TO N-1 {
	  ns[i] = 0
	  ls[i] = -10000
	}
	N2 = N + sqrt(correl) * (1-N)
	go()
}

BREAKPOINT {
	SOLVE go
}


UNITSOFF

PROCEDURE go() { LOCAL i, j, sum

: reset all channels

   FROM i=0 TO N-1 {
	x[i] = min_val
   }

   if( (on==1) && (t>=latency) && (t<=shutoff) )  {

: Determine how to distribute random variables among the N outputs:

     if(correl==0) {		: If correl=0, create N random variables

	FROM i=0 TO N-1 {
	   if(get_random(1) <= spont) { 	: toss coin...
		x[i] = max_val			: spike!
	   } else {
		x[i] = min_val
	   }
	}

     } else {			: if correl>0, distribute N2 rnd in N outputs

	FROM i=0 TO N2-1 {		: first generate N2 random variables
	   if(get_random(1) <= spont) { 	: toss coin...
		R[i] = max_val			: spike!
	   } else {
		R[i] = min_val
	   }
	}

	FROM i=0 TO N-1 {		: scan each output variable
	   j = get_random(N2)			: chose one of the rnd variable
	   x[i] = R[j]				: assign the rnd variable
	}

     }

: update the counters and check for refractory period

     FROM i=0 TO N-1 {
	   if((t-ls[i]) <= refract) {		: force to zero if refractory
		x[i] = min_val
	   }
	   if(x[i] == max_val) {		: if variable has fired
		   ns[i] = ns[i] + 1			: increase spike count
		   ls[i] = t				: memorize last spike
	   }
     }


   } else if(on==2) {
	FROM i=0 TO N-1 {		: fire all channels
	   x[i] = max_val				: spike !
	   ns[i] = ns[i] + 1				: increase spike count
	   ls[i] = t					: memorize last spike
	}
	on = 1				: reset to normal

   } else if(on==3) {
	FROM i=0 TO N-1 {		: fire channels selected by sync
	   if(sync[i]) {
		x[i] = max_val				: spike !
		ns[i] = ns[i] + 1			: increase spike count
		ls[i] = t				: memorize last spike
	   }
	}
	on = 1				: reset to normal
   }
}
UNITSON


FUNCTION get_random(maxval) {			: simple random
	get_random = maxval * random() / (2^31)
}


PROCEDURE new_seed(seed) {		: procedure to set the seed
	srandom(seed)
	VERBATIM
	  printf("Setting random generator with seed = %g\n", _lseed);
	ENDVERBATIM
}

PROCEDURE printvec() { LOCAL i		: procedure to print the vectors
   VERBATIM 
   {
	int i;
	printf("i\tx\tns\tls\n");
	for(i=0; i<N; i++) {
	  printf("%d\t%g\t%g\t%g\n",i,(float)x[i],(float)ns[i],(float)ls[i]);
	}
   } 
   ENDVERBATIM
}

