
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS IVNMDA5
	POINTER C
	RANGE C0, C1, C2, D, O, B
	RANGE g, gmax, rb, gluc, w, timer, llabel
	GLOBAL Erev, mg, Rb, Ru, Rd, Rr, Ro, Rc
	GLOBAL vmin, vmax

      NONSPECIFIC_CURRENT inmda
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(uM) = (micro/liter)
}

PARAMETER {

	Erev	= 0    (mV)	: reversal potential
	gmax	= 500  (pS)	: maximal conductance
	mg	= 0    (mM)	: external magnesium concentration
	vmin = -120	(mV)
	vmax = 100	(mV)
: Rates
	gluc = 0.00003 (mM)
	timer = 110 (ms)
	llabel = 0
	
	: Destexhe, Mainen & Sejnowski, 1996
	Rb	= 5e-3    (/uM /ms)	: binding 		
	Ru	= 12.9e-3  (/ms)	: unbinding		
	Rd	= 8.4e-3   (/ms)	: desensitization
	Rr	= 6.8e-3   (/ms)	: resensitization 
	Ro	= 46.5e-3   (/ms)	: opening
	Rc	= 73.8e-3   (/ms)	: closing
}

COMMENT
	: Clements et al. 1992
	Rb	= 5e-3    (/uM /ms)	: binding 		
	Ru	= 9.5e-3  (/ms)	: unbinding		
	Rd	= 16e-3   (/ms)	: desensitization
	Rr	= 13e-3   (/ms)	: resensitization 
	Ro	= 25e-3   (/ms)	: opening
	Rc	= 59e-3   (/ms)	: closing

	: Hessler Shirke & Malinow 1993
	Rb	= 5e-3    (/uM /ms)	: binding 		
	Ru	= 9.5e-3  (/ms)	: unbinding		
	Rd	= 16e-3   (/ms)	: desensitization
	Rr	= 13e-3   (/ms)	: resensitization 
	Ro	= 25e-3   (/ms)	: opening
	Rc	= 59e-3   (/ms)	: closing

	: Clements & Westbrook 1991
	Rb	=  5    (uM /s)	: binding 		
	Ru	=  5	(/s)	: unbinding -> gives Kd = Rb/Ru = 1 uM
	Rd	=  4.0  (/s)	: desensitization
	Rr	=  0.3  (/s)	: resensitization 
	Ro	= 10  (/s)	: opening
	Rc	= 322   (/s)	: closing

	: Edmonds & Colquhoun 1992
	Rb	=  5    (uM /s)	: binding 		
	Ru	=  4.7  (/s)	: unbinding		
	Rd	=  8.4  (/s)	: desensitization
	Rr	=  1.8  (/s)	: resensitization 
	Ro	= 46.5  (/s)	: opening
	Rc	= 91.6  (/s)	: closing

	: Lester & Jahr 1992
	Rb	= 5    (uM /s)	: binding 		
	Ru	= 6.7   (/s)	: unbinding		
	Rd	= 15.2  (/s)	: desensitization
	Rr	= 9.4   (/s)	: resensitization 
	Ro	= 83.8  (/s)	: opening
	Rc	= 83.8  (/s)	: closing

ENDCOMMENT


ASSIGNED {
	v		(mV)		: postsynaptic voltage
        inmda               (nA)            : current = g*(v - Erev)
	g 		(pS)		: conductance
	C 		(mM)		: pointer to glutamate concentration

	rb		(/ms)    : binding
	w
}

STATE {
	: Channel states (all fractions)
	C0		: unbound
	C1		: single bound
	C2		: double bound
	D		: desensitized
	O		: open

	B		: fraction free of Mg2+ block
}

INITIAL {
	rates(v)
	C0 = 1
	llabel = 0.0
}

BREAKPOINT {
	rates(v)
	SOLVE kstates METHOD sparse

	g = gmax * O * B
      inmda = (1e-6) * g * (v - Erev)

:	 VERBATIM
                      
       llabel = writeIV(inmda, gluc, timer, llabel)
	
:          return 0;
:       ENDVERBATIM

	
}

KINETIC kstates {
	
	rb = Rb * (1e3) * C 

        ~ C0 <-> C1     (2*rb,Ru)
        ~ C1 <-> C2     (rb,2*Ru)
	~ C2 <-> D	(Rd,Rr)
	~ C2 <-> O	(Ro,Rc)

	CONSERVE C0+C1+C2+D+O = 1
}

PROCEDURE rates(v(mV)) {
	TABLE B
	DEPEND mg
	FROM vmin TO vmax WITH 200

	: from Jahr & Stevens

	B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}

VERBATIM
double peak = 0.0;
double writeIV(double(in),double(cc), double(timerr), double(label)) {

	/*stops at peak not used now*/		 
                      FILE *fp;

/*                 		printf("%f ", cc);*/
			    if (peak*1e10 >= in*1e10) {
					peak = in;
/*					printf("%e	%e\n", peak, in);*/
					if(label == 1.0) {
						label = 0.0;
					}
					if(label == 2.0) {
						label = 0.0;
					}
					if(label == 3.0) {
						peak = 0;
					}
								
			    }

			    if ((peak*1e10 < in*1e10) && (label != 3.0)) {
  			         
				   label = label + 1.0;
				   if(label == 1.0){
					peak = peak;
				   } else if (label == 3.0) {
    					   fp = fopen("IV.txt", "a");

	                      	   if(fp == NULL)  {
#if 0
						      printf("Can't open IV.txt file\n");
						      exit(1);
#else
							hoc_execerror("Can't open IV.txt file", (char*)0);
#endif
	                           }
					
					   printf("peak value is %e\n", peak);	
					   fprintf(fp, "%e	%e \n", cc, peak);
					   fclose(fp);
					   peak = 0.0;
				   }				
			   }
			   if ((t+dt/2) > timerr ) {
				printf("baby one");
		            if (label == 0.0) {
  			         fp = fopen("IV.txt", "a");
                      	   if(fp == NULL)  {
#if 0
					      printf("Can't open IV.txt file\n");
					      exit(1);
#else
							hoc_execerror("Can't open IV.txt file", (char*)0);
#endif
	                     }
				   fprintf(fp, "%e %e \n", cc, in);
					   
                           fclose(fp);
				   label = 3.0;
				   peak = 0.0;		
				}	
					  
			   }	
			                			
				return(label);		
}
ENDVERBATIM


               


