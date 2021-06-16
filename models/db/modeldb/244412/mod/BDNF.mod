COMMENT
BDNF kinetics

Includes extracellular and intracellular mechanisms.

Based on data provided by V. Lessmann and Kurt ???

Developers: Solinas & Migliore 2017

BDNF release 
The proBDNF-containing vesicles are released with 100 s delay 
after an increase of cai. To model this we must explicity model BDNF vesicles
as events triggered by a cai thereshold that execute a net_send with a probability 
that is proportional to the inverse of the time step.

ENDCOMMENT

NEURON {
    POINT_PROCESS BDNF
    USEION ca READ cai	: Weight update requires cai 
    : USEION bdnf READ bdnfi WRITE bdnfi VALENCE 0
    
    RANGE max_BDNF_rel_delay, theta_cai_BDNF, max_cai_BDNF, BDNF_prel, fused_vesicles, duration_BDNF_release
    RANGE proBDNF_uptake,PC_uptake,mBDNF_uptake, TrkB, proBDNF_fraction
    RANGE tau_LTP14, theta_gAMPA, sigma_gAMPA, alpha_gAMPA, shift_gAMPA, v_BDNF, is, intracell_signaling
    POINTER randObjPtr
    POINTER gAMPA
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (mM) = (millimolar)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

PARAMETER {
    cai (mM)
    
    dt (ms)
       
    :BDNF
    : [Ca]i threshold for BDNF vesicle release
    theta_cai_BDNF = 0.11 (mM)
    max_cai_BDNF = 0.13 (mM)
    
    max_BDNF_rel_delay = 400e3 (ms)
    duration_BDNF_release = 1800e3 :1800e3 (ms) :30*60*1e3 = 30 min
    proBDNF_uptake = 0.00001 (mM/ms)    
    mBDNF_uptake = 0.00001 (mM/ms)    
    PC_uptake = 0.00001 (mM/ms)    
    
    : Time step for BDNF vesicle release
    cai_integration_time_step = 1 (ms)
    : Normalisation factor to convert [Ca]i * cai_integration_time_step to a unitless number
    alpha_dt_cai = 1.2 (1/mM 1/ms)
    
    : Vesicles
    n_vesicles = 300 (1)
    v_BDNF = 0.002 (mM)
    proBDNF_fraction = 0.7
    v_PC = 0.002 (mM)
    
    : Cleavage
    tau_cleave = 10000.0 (ms/mM) : alpha_cleave = 1e-4 (mM/ms)
    :rb_BDNF = 0 : uncleaving is not allowd ??
    
    :TrkB
    tau_LTP14 = 180e3 :1800e3 (ms) : 1 min = 60e3 ms, 30 min = 1800e3 ms
    theta_TrkB = 0.0002 (mM)
    sigma_TrkB = 0.00001 (mM)
    
    gbar_AMPA = 1 (nS)
    scale_AMPA = 100 (1)
    alpha_gAMPA = 0.5
    theta_gAMPA = 0.5
    sigma_gAMPA = 0.1
    shift_gAMPA = 0    
    
    : intracell signaling is udes to have a decay of [Ca]i efficacy on BDNF release
    : so that only stimuli given at 0.5 Hz are effective on triggering BDNF release 
    : while stimuli given at 0.05 Hz are not effective
    is_decay = 0.1e-3 (/ms) : tau = 10 sec, decay rate of intracell_signaling (is)
    
    : max_fused_vesicles = 0
}

ASSIGNED {
    BDNF_prel (1)
    randObjPtr
    cai_th_crossed (1)
    gAMPA (nS)
    gAMPA_g
    alpha_LTP14 (mM/ms)
    mBDNF_fraction
    fusion_delay (ms)
    v_prel_norm (1)
    n_avail_vesicles (1)
    cai_factor (1)
    : alpha_cleave = 1 / tau_cleave
}

STATE {
    : v_proBDNF (mM)
    : v_PC (mM)
    fused_vesicles (1)
    proBDNF (mM)
    mBDNF (mM)
    ppBDNF (mM) 
    PC (mM)  : tPA in paper, PC in Lessman sketch
    proBDNF_PC (mM)
    TrkB (mM)
    proBDNF_removed (mM)
    mBDNF_removed (mM)
    PC_removed (mM)
    is (1)
    intracell_decay (1)
    intracell_signaling (mM)
}

INITIAL {    
    fused_vesicles = 0
    cai_th_crossed = 0
        
    proBDNF = 0 
    mBDNF = 0 
    ppBDNF = 0
    PC = 0
    proBDNF_removed = 0 
    mBDNF_removed = 0 
    ppBDNF = 0
    PC_removed = 0
    is = 0    
    intracell_decay = 0
    intracell_signaling = 0
    mBDNF_fraction = 1 - proBDNF_fraction

    alpha_LTP14 = 1/tau_LTP14
    
    TrkB = 0 (mM)
    : v_prel_norm = 1/(n_vesicles * alpha_dt_cai * cai_integration_time_step * theta_cai_BDNF)
    n_avail_vesicles = n_vesicles
}

BREAKPOINT {
    SOLVE kstates METHOD sparse
    : gAMPA = gbar_AMPA * TrkB
    : gAMPA = gAMPA_g
}

KINETIC kstates {
    
    
    : if (t/60e3 > 0.1) {
    : 	if (t/60e3 < 4) {
    : 	    fused_vesicles = max_fused_vesicles
    : 	    : printf("fused %g maxfu %g\n",fused_vesicles,max_fused_vesicles)
    : 	}
    : }
    
    
    ~ is <-> intracell_decay (is_decay,0)
    : printf("isk %g\n",is)
    : Release the proBDNF, mBDNF, and PC
    : Release is mantianed for a time of x
    : Protracted vesicular release 
    : printf("BDNF curr: %g\n",v_BDNF * proBDNF_fraction * fused_vesicles / duration_BDNF_release)
    ~ proBDNF << (v_BDNF * proBDNF_fraction * fused_vesicles / duration_BDNF_release)
    ~ mBDNF << (v_BDNF * mBDNF_fraction * fused_vesicles / duration_BDNF_release)
    ~ PC << (v_PC * fused_vesicles / duration_BDNF_release)
    
    : proBDNF cleavage
    ~ proBDNF + PC <-> proBDNF_PC (1/tau_cleave,0)  
    ~ proBDNF_PC <-> ppBDNF + mBDNF + PC (1/tau_cleave,0)  
    : Uptake or diffusion from synaptic cleft
    ~ proBDNF <-> proBDNF_removed (proBDNF_uptake, 0)
    ~ mBDNF <-> mBDNF_removed (mBDNF_uptake, 0)
    ~ PC <-> PC_removed (PC_uptake, 0)
    ~ fused_vesicles <-> fused_vesicles (0,0)

    : TrkB activation
    TrkB = mBDNF * sigh(mBDNF, theta_TrkB, sigma_TrkB)
    ~ TrkB <-> intracell_signaling (alpha_LTP14, 0)
    
    
    :gAMPA_g = gbar_AMPA * TrkB/scale_AMPA
    gAMPA = 1 + alpha_gAMPA * sigh(intracell_signaling, theta_gAMPA, sigma_gAMPA) + shift_gAMPA
    
}


FUNCTION sigh(x (mM), theta (mM), sigma (mM)) {
    : LOCAL e
    : e = (x - theta) / sigma
    : if ( -e > 700 ) {
    : 	printf("%f\t",(-(x - theta) / sigma))
    : }
    sigh = 1 / (1 + exp((theta - x) / sigma))
}

VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

FUNCTION randGen() {
VERBATIM
   if (_p_randObjPtr) {
      /*
      :Supports separate independent but reproducible streams for
      : each instance. However, the corresponding hoc Random
      : distribution MUST be set to Random.uniform(0,1)
      */
      _lrandGen = nrn_random_pick(_p_randObjPtr);
   }else{
      hoc_execerror("Random object ref not set correctly for randObjPtr"," only via hoc Random");
   }
ENDVERBATIM
}
 
PROCEDURE setRandObjRef() {
VERBATIM
   void** pv4 = (void**)(&_p_randObjPtr);
   if (ifarg(1)) {
      *pv4 = nrn_random_arg(1);
   }else{
      *pv4 = (void*)0;
   }
ENDVERBATIM
}

FUNCTION min(x,y) { if (x<=y){ min = x }else{ min = y } }

NET_RECEIVE (weight (1)) { LOCAL prel, is_effect
    if ((flag == 0) && (cai_th_crossed == 0) ) { : the netcon is indicating that the theta_cai_BDNF was crossed
	cai_th_crossed = 1
	: step increase of intracell signaling caused by [Ca]i LTP14 threshold crossing
	is = is + 0.1
	
	: printf("is %g\n",is)
	: printf("Cai crossed th at t %g\n",t)
	: Calculate a probability of release proportional to cai e to dt only if there are vesicles for fusion
	if (n_avail_vesicles > 0) {
	    net_send(cai_integration_time_step,2) : keep on watching cai
	}
    }
    
    if (flag == 2 && n_avail_vesicles > 0) {
    	: Calculate a probability of release proportional to cai, cai_dt, :and number of available vesicles
	cai_factor = (cai - theta_cai_BDNF) / (max_cai_BDNF - theta_cai_BDNF)
	if (is > 0.15) {
	    is_effect = 1
	} else {
	    is_effect = 0
	}
	prel = alpha_dt_cai * cai_integration_time_step * min(1,cai_factor) * is_effect: * n_avail_vesicles/n_vesicles
	: prel = (1-((n_vesicles-n_avail_vesicles)/n_vesicles)^2) * alpha_dt_cai * cai_integration_time_step * cai
	: printf("Prel: vesicles, constant, cai_ratio %g\t%g\t%g\n", n_avail_vesicles/n_vesicles, alpha_dt_cai * cai_integration_time_step, (cai - theta_cai_BDNF) / (max_cai_BDNF - theta_cai_BDNF))
    	if ( randGen() < prel ) {
	    BDNF_prel = prel
    	    : Relase one BDNF vesicle in the future
    	    : printf("Release bdnf t %g\n",t)
	    : Calc the delay of vesicle fusion, :when many vesicles are available the delay is shorter
	    fusion_delay = max_BDNF_rel_delay * (1 - min(1,cai_factor)) * randGen()
    	    net_send(fusion_delay, 3)
	    n_avail_vesicles = n_avail_vesicles - 1
    	    : printf("future: dt %g\t%g\t%g\t%g\n", t+fusion_delay, max_BDNF_rel_delay, max_BDNF_rel_delay  * (1 - min(1,cai_factor)), (cai - theta_cai_BDNF) / (max_cai_BDNF - theta_cai_BDNF))
    	    : printf("Trigger bdnf release in the future: dt %g\t%g\t%g\t%g\n",BDNF_prel, t+fusion_delay, max_BDNF_rel_delay, n_avail_vesicles)
    	    :printf("Trigger bdnf release in the future: dt %g\t%g\t%g\n",BDNF_prel, cai, dt)
    	}
	: Continue releasing each ms till cai is above threshold only if there are vesicles available for fusion
    	if (cai > theta_cai_BDNF) {
	    if (n_avail_vesicles > 0) {
    		net_send(cai_integration_time_step,2) : keep on watching cai    
	    }
	} else {
	    cai_th_crossed = 0
	}
    }
    
    if (flag == 3) { 
	: Increase counter of fused vesicles
	fused_vesicles = fused_vesicles + 1
	: printf("Fusing a vesicle %g\t%g\n",fused_vesicles,t)
	net_send(duration_BDNF_release,4)
    }
    
    if (flag == 4) { 
	: Decrease counter of fused vesicles
	fused_vesicles = fused_vesicles - 1
    }
}
