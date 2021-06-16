TITLE Internal noisy channel

:
: Include internal noisy channels with density form
: With cut off frequency 

: in=-norm(mean,std)  

: by Yi Zhou for MSO use


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS  current_gauss
	
	NONSPECIFIC_CURRENT  in  : negative current depolarizes the membrane
	RANGE del,dur
        
	RANGE in
	RANGE rand
	RANGE mean, std0, std
	RANGE f0  : the sampling frequency
	RANGE N_smooth : =1000/f0/dt
	RANGE count
	RANGE noise_seed
	
}


UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
}

PARAMETER {
	del=0 (ms)
	dur=100 (ms) <0,1e9>
	
	mean=0         (nA)  : default
	std0=1e-4       (nA)  : default
	std=1e-4
	noise_seed=1
    
	f0=4000   (1/s)   : default 4000Hz 
	
}



ASSIGNED {
	v      (mV)
	dt	(ms) 

	in	(nA)
	
	rand
	count
        N_smooth
	
}

PROCEDURE seed(x) {
	set_seed(x)
}


BREAKPOINT {

	if (t<del+dur && t>del) {

		if(count/N_smooth==1){  :error with the first N_smooth-1 zeros paddings, i.e., 1/f0 ms
			rand = normrand(mean, std)
			in  = -rand   : depolaring current
			count=0
		}else{
			count=count+1
		} 
	} else {
	  in=0 
	}

:printf("count=%g\n",count)

}



INITIAL {
	
	N_smooth=floor(1000/f0/dt)*2  : break point called twice,here 1000 is a scalar

		:### equalize the power of the sampled white noise
		:### std_a/std_b=sqrt(f_a/f_b) with f_a=4000 and std0=std_a : see notebook p186
	std=std0/sqrt(4000/f0)  
        	:###note that for Wiener process should be std_a/std_b=sqrt(f_b/f_a)
	
	rand = normrand(mean, std)
	in  = -rand

	count=0
	
	}

UNITSON



