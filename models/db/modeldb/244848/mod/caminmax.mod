NEURON { 
        SUFFIX caminmax     
		USEION ca READ cai
        RANGE camax, camin, tmax, tmin
} 

UNITS 
{
	(mM) = (milli/liter)
}



ASSIGNED { 
        cai (mM) 
        camax (mM) 
		camin (mM)
        tmax (ms) 
		tmin (ms)
} 

INITIAL { 
        camax = cai
        camin = cai
        tmax = t 
        tmin = t 
} 

BREAKPOINT { 

VERBATIM 
	if (cai > camax) { 
		camax = cai; 
		tmax = t;               
	}           

	if (cai < camin) { 
		camin = cai; 
		tmin = t;               
	}           
   
ENDVERBATIM 
}               
