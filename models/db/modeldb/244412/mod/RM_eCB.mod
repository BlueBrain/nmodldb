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
    POINT_PROCESS RM_eCB
    USEION ca READ cai	: Weight update requires cai 
    
    :RM
    RANGE tau_RM
    :cai->RM
    RANGE alpha_cai_RM,theta_cai_RM, sigma_cai_RM
    : RM->RMr
    RANGE tau_RMLTP11, RMr
    : RMr->U_SE
    RANGE alpha_RMru, theta_RMru, sigma_RMru
    :cai->RMBLK
    RANGE alpha_cai_RMBLK, theta_cai_RMBLK, sigma_cai_RMBLK
    :RMBLK
    RANGE tau_RMBLK
    :RM_RMr
    RANGE tau_RM_RMr, alpha_RM_RMr, theta_RM_RMr, sigma_RM_RMr

    : delta_U->U_SE
    POINTER delta_U
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
    
    : RM
    tau_RM = 1800e3 (ms)
    RM_inf = 0 (mM)
    alpha_cai_RM = 1e-3 (mM/ms)
    theta_cai_RM = 1.5e-3 (mM)
    sigma_cai_RM = 0.01e-3 (mM)
    
    tau_RMLTP11 = 1800e3 (ms) :1800e3 (ms) : 1 min = 60e3 ms, 30 min = 1800e3 ms
    RMr_0 = 0 (mM)
    
    alpha_RMru = 0.3 (1)
    theta_RMru = 0.5 (mM)
    sigma_RMru = 0.1 (mM)
    
    :RMBLK
    tau_RMBLK = 1800e3 (ms)
    RMBLK_inf = 0 (mM)
    alpha_cai_RMBLK = 1e-3 (mM/ms)
    theta_cai_RMBLK = 1.8e-3 (mM)
    sigma_cai_RMBLK = 0.01e-3 (mM)
    
    theta_RMBLK = 0.015 (mM)
    sigma_RMBLK = 0.001 (mM)
    
    :RM
    tau_RM_RMr = 1e3 (ms)
    alpha_RM_RMr = 1
    theta_RM_RMr = 0.02 (mM)
    sigma_RM_RMr = 0.001 (mM)
}

ASSIGNED {
    delta_U (1)
}

STATE {
    RM (mM)
    RMr (mM)
    RMBLK (mM)
    post_intra (mM)
}

INITIAL {
    RM = RM_inf
    RMr = RMr_0
    : printf("sigcai%f\t", sigh(RMr,theta_RMr,sigma_RMr))
    RMBLK = RMBLK_inf
    : RMBLKe = RMBLKe_0
    post_intra = 0
    
    delta_U = alpha_RMru * sigh(RMr,theta_RMru,sigma_RMru): * (1 - sigh(RMBLKe,theta_RMBLKe,sigma_RMBLKe))
}

BREAKPOINT {
    SOLVE state METHOD cnexp
}

DERIVATIVE state {
        
    RMBLK' = (RMBLK_inf - RMBLK)/tau_RMBLK  + alpha_cai_RMBLK * sigh(cai,theta_cai_RMBLK,sigma_cai_RMBLK): - alpha_rmep * RMBLK: + (RMBLK_inf - RMBLK) / tau_RMBLKLTP11
    : RMBLKe' = (RMBLK - RMBLK_inf) / tau_RMBLKLTP11 : thereshould be a decay here otherwize the LTP11 accumulates, this is a memory eraser
    
    : We have to use a negative rate to decrease RM if during induction the RMBLK is activated to prevent RM accumulation.
    : This way is preferreble to the product by a sigmoid of delta_u, see below, as it does not depotentiate the effects induced by an LTP11.
    : This should be further clarified as depotentiation could be important here.
    : Note that for this to work must be alpha_rmep > alpha_cai_RM
    RM' = (RM_inf - RM)/tau_RM  + alpha_cai_RM * sigh(cai,theta_cai_RM,sigma_cai_RM) * (1 - sigh(RMBLK,theta_RMBLK,sigma_RMBLK))  + (RM_inf - RM) / tau_RM_RMr * alpha_RM_RMr * sigh(RM,theta_RM_RMr,sigma_RM_RMr)
    RMr' = (RM - RM_inf) / tau_RM_RMr * alpha_RM_RMr * sigh(RM,theta_RM_RMr,sigma_RM_RMr) - RMr / tau_RMLTP11
    : there should be a decay here otherwize the effects of LTP11 accumulates, this would a memory eraser
    
    : post_intra is a intracellular buffer for RMr states that will be converted to post_intra with time const tau_RMLTP11 
    post_intra' = RMr / tau_RMLTP11
    delta_U = alpha_RMru * sigh(post_intra,theta_RMru,sigma_RMru)
}

FUNCTION sigh(x (mM), theta (mM), sigma (mM)) {
    : LOCAL e
    : e = (x - theta) / sigma
    : if ( -e > 700 ) {
    : 	printf("%f\t",(-(x - theta) / sigma))
    : }
    sigh = 1 / (1 + exp((theta - x) / sigma))
}
