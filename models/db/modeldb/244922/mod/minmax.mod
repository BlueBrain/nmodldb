NEURON { 
        SUFFIX minmax     
        RANGE vmax, vmin, tmax, tmin
} 

ASSIGNED { 
        v (millivolt) 
        vmax (millivolt) 
		vmin (millivolt)
        tmax (ms) 
		tmin (ms)
} 

INITIAL { 
        vmax = v 
        vmin = v 
        tmax = t 
        tmin = t 
} 

BREAKPOINT { 

VERBATIM 
	if (v > vmax) { 
		vmax = v; 
		tmax = t;               
	}           

	if (v < vmin) { 
		vmin = v; 
		tmin = t;               
	}           
   
ENDVERBATIM 
}               
