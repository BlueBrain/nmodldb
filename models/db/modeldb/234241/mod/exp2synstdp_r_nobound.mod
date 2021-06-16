: STDP by Hines, modified by Michiel to implement additive STDP as in Song et al 2000, 
: modified by Katha to inhibitory synapse with Exp2Syn conductance change and 
: STDP rule is reversed and can be shifted
: use net_send to flag3 for positive shift, and to flag 4 for negative shift

NEURON {
        POINT_PROCESS Exp2Syn_Inh_STDP_nobound
        RANGE tau1, tau2, e, i, dd, dp, dtau, ptau, thresh, wmax, wmin, D, shift, g, mean, std, learning_rate
        NONSPECIFIC_CURRENT i
}

UNITS {
        (nA) = (nanoamp)
        (mV) = (millivolt)
        (uS) = (microsiemens)
}

PARAMETER {
        tau1     = 0.1 (ms) <1e-9,1e9>
        tau2    = 10 (ms) <1e-9,1e9>
        e       = -73     (mV)
        dd      = 0.001 <0,1>   : depression factor (relative!)
        dp      = 0.001       : potentiation factor (relative!)
        dtau    = 20 (ms)   : depression effectiveness time constant
        ptau    = 20 (ms)   : Bi & Poo (1998, 2001)
        thresh  = -20 (mV)      : postsynaptic voltage threshold
        wmax    = 0.001 (uS)
        wmin    = 0 (uS)
        shift   = 0 (ms) 	: shift 
        mean    = 0 (ms)
        std     = 0 (ms)
        learning_rate = 1

}

ASSIGNED {
        v (mV)
        i (nA)
        D
        tpost (ms)
        g        
        factor
        jitter
}

STATE {
        A (uS)
        B (uS)
}

INITIAL {
        LOCAL tp
        A = 0
        B = 0
        D = 0
        tpost = -1e9
        jitter = normrand(mean, std)
        net_send(0, 1)
        if (tau1/tau2 > .9999){
            tau1 = .9999*tau2
        }
        tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
        factor = -exp(-tp/tau1) + exp(-tp/tau2)
        factor = 1/factor
}

BREAKPOINT {
        SOLVE state METHOD cnexp
        g = B-A
        i = g*(v - e)
}

DERIVATIVE state {
        A' = -A/tau1
        B' = -B/tau2
}

NET_RECEIVE(w (uS), P, tpre (ms), wsyn) {
        INITIAL { P = 0 tpre = -1e9 wsyn = w}

        if (flag == 0) { : presynaptic spike  (after last post so potentiate)      
                :printf("entry flag=%g t=%g w=%g A=%g tpre=%g tpost=%g g=%g\n", flag, t, w, P, tpre, tpost,g)
                :printf("wsyn=%g\n", wsyn)
                :printf("shift=%g\n",shift)
        net_send(jitter,5)
        net_send(jitter+shift, 3)        
        :P = P*exp((tpre-t)/ptau) - dp
        :tpre = t+shift
        :wsyn = wsyn + wmax * D * exp((tpost-t)/dtau) : interval is negative
        :if (wsyn > wmax) { wsyn = wmax }
        
        }else if (flag == 2) { : postsynaptic spike                             
                :printf("entry flag=%g t=%g tpost=%g\n", flag, t, tpost)
        FOR_NETCONS(w1, P1, tpre1, wsyn1) {
            wsyn1 = wsyn1 + learning_rate * wmax * P1 * exp(-(t-tpre1)/ptau) : interval is positive
            if (wsyn1 < wmin) { wsyn1 = wmin }
        }
        :net_send(shift, 4)     
        D = D*exp((tpost-t)/dtau) + dd
        tpost = t

        } else if (flag == 3) { 
                :printf("entry flag=%g \n", flag)
        P = P*exp((tpre-t)/ptau) - dp
        tpre = t
        wsyn = wsyn + learning_rate * wmax * D * exp((tpost-t)/dtau) : interval is negative
        :if (wsyn > wmax) { wsyn = wmax }

        :} else if (flag == 4) { 
        :        :printf("entry flag=%g \n", flag)
        :D = D*exp((tpost-t)/dtau) + dd
        :tpost = t
        :FOR_NETCONS(w1, P1, tpre1, wsyn1) {
        :    wsyn1 = wsyn1 + learning_rate * wmax * P1 * exp(-(t-tpre1)/ptau) : interval is positive
        :    if (wsyn1 < wmin) { wsyn1 = wmin }
        :}
        } else if (flag == 5) {
        A = A + wsyn*factor
        B = B + wsyn*factor
        } else { : flag == 1 from INITIAL block                                 
        :printf("entry flag=%g t=%g\n", flag, t)
                WATCH (v > thresh) 2
        }
}