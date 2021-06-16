: STDP by Hines, modified by Michiel to implement additive STDP as in Song et al 2000

NEURON {
        POINT_PROCESS ExpSynCaSTDP
        RANGE tau, e, i, dd, dp, dtau, ptau, thresh, ca_thresh, wmax, wmin, D, cai
        NONSPECIFIC_CURRENT i
        USEION ca READ cai
}

UNITS {
        (nA) = (nanoamp)
        (mV) = (millivolt)
        (uS) = (microsiemens)
}

PARAMETER {
        tau     = 3 (ms) <1e-9,1e9>
        e       = 0     (mV)
        dd      = 0.001 <0,1>   : depression factor (relative!)
        dp      = 0.00106       : potentiation factor (relative!)
        dtau    = 20 (ms)   : depression effectiveness time constant
        ptau    = 20 (ms)   : Bi & Poo (1998, 2001)
        thresh  = -20 (mV)      : postsynaptic voltage threshold
        ca_thresh  = 0.5 (mM)      : postsynaptic calcium threshold
        wmax    = 0.001 (uS)
        wmin    = 0 (uS)
}

ASSIGNED {
        v (mV)
        cai (mM)
        i (nA)
    D
        tpost (ms)
}

STATE {
        g (uS)
}

INITIAL {
        g = 0
        D = 0
        tpost = -1e9
        net_send(0, 1)
}

BREAKPOINT {
        SOLVE state METHOD cnexp
        i = g*(v - e)
}

DERIVATIVE state {
        g' = -g/tau
}

NET_RECEIVE(w (uS), P, tpre (ms), wsyn) {
        INITIAL { P = 0 tpre = -1e9 wsyn = w}
        if (flag == 0) { : presynaptic spike  (after last post so depress)      
	printf("entryca flag=%g t=%g w=%g A=%g tpre=%g tpost=%g g=%g ca=%g ca_thresh=%g\n", flag, t, w, P, tpre, tpost,g,cai,ca_thresh)
        printf("wsynca=%g\n", wsyn)
        
        g = g + wsyn
        P = P*exp((tpre-t)/ptau) + dp
        tpre = t
        wsyn = wsyn + wmax * D * exp((tpost-t)/dtau) : interval is negative
        if (wsyn < wmin) { wsyn = wmin }
        }else if (flag == 2) { : postsynaptic spike                             
        printf("entryca flag=%g t=%g tpost=%g ca=%g\n", flag, t, tpost,cai)
        FOR_NETCONS(w1, P1, tpre1, wsyn1) {
            :D1 = D1*exp((tpost-t)/dtau) - dd
            wsyn1 = wsyn1 + wmax * P1 * exp(-(t - tpre1)/ptau) : interval is positive
            if (wsyn1 > wmax) { wsyn1 = wmax }
        }
        D = D*exp((tpost-t)/dtau) - dd
        :wsyn =  0
        :wsyn = wsyn + wmax * P * exp(-(t - tpre)/ptau)
        :if (wsyn > wmax) { wsyn = wmax }
        tpost = t
        } else { : flag == 1 from INITIAL block                                 
        printf("entryca flag=%g t=%g ca=%g\n", flag, t,cai)
                WATCH (cai > ca_thresh) 2
        }
}