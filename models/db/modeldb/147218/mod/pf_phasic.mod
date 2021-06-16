: phasic synaptic current of parallel fibers

NEURON {
    SUFFIX pf_phasic
    RANGE nb_epsc,freq,del,tauOpf,tauCpf
    NONSPECIFIC_CURRENT  i
    RANGE i, e, g
}

PARAMETER {
    nb_epsc = 10 : number of EPSCs 
    freq = 50 (/s) : frequency of EPSCs
    g = 0 (siemens/cm2)  < 0, 1e9 >
    e = 0    (millivolts)
    del = 100 (ms)
    tauOpf = 2.4 (ms)
    tauCpf = 6.3 (ms)
    i1 =0 (millivolts*siemens/cm2) 
    i2 =0 (millivolts*siemens/cm2)
    i3 =0 (millivolts*siemens/cm2)
    i4 =0 (millivolts*siemens/cm2)
    i5 =0 (millivolts*siemens/cm2)
    i6 =0 (millivolts*siemens/cm2)
    i7 =0 (millivolts*siemens/cm2)
    i8 =0 (millivolts*siemens/cm2)
    i9 =0 (millivolts*siemens/cm2)
    i10 =0 (millivolts*siemens/cm2)
}

ASSIGNED {
    i   (milliamp/cm2)
    v   (millivolt) 
}

INITIAL  { i = 0   }

BREAKPOINT {
    at_time(del)            if (t<del)            {i1 =0}  else {if (1>nb_epsc)  {i1 = 0}  else  { i1  = g*(1-exp(-(t-(del))/tauOpf))*exp(-(t-(del))/tauCpf)*(v - e) }}
    at_time(del+1e3/freq)   if (t<del+1e3/freq)   {i2 =0}  else {if (2>nb_epsc)  {i2 = 0}  else  { i2  = g*(1-exp(-(t-(del+1e3/freq))/tauOpf))  *exp(-(t-(del+1e3/freq))/tauCpf)*(v - e) }}
    at_time(del+2*1e3/freq) if (t<del+2*1e3/freq) {i3 =0}  else {if (3>nb_epsc)  {i3 = 0}  else { i3  = g*(1-exp(-(t-(del+2*1e3/freq))/tauOpf))*exp(-(t-(del+2*1e3/freq))/tauCpf)*(v - e) }}
    at_time(del+3*1e3/freq) if (t<del+3*1e3/freq) {i4 =0}  else {if (4>nb_epsc)  {i4 = 0}  else { i4  = g*(1-exp(-(t-(del+3*1e3/freq))/tauOpf))*exp(-(t-(del+3*1e3/freq))/tauCpf)*(v - e) }}
    at_time(del+4*1e3/freq) if (t<del+4*1e3/freq) {i5 =0}  else {if (5>nb_epsc)  {i5 = 0}  else { i5  = g*(1-exp(-(t-(del+4*1e3/freq))/tauOpf))*exp(-(t-(del+4*1e3/freq))/tauCpf)*(v - e) }}
    at_time(del+5*1e3/freq) if (t<del+5*1e3/freq) {i6 =0}  else {if (6>nb_epsc)  {i6 = 0}  else { i6  = g*(1-exp(-(t-(del+5*1e3/freq))/tauOpf))*exp(-(t-(del+5*1e3/freq))/tauCpf)*(v - e) }}
    at_time(del+6*1e3/freq) if (t<del+6*1e3/freq) {i7 =0}  else {if (7>nb_epsc)  {i7 = 0}  else { i7  = g*(1-exp(-(t-(del+6*1e3/freq))/tauOpf))*exp(-(t-(del+6*1e3/freq))/tauCpf)*(v - e) }}
    at_time(del+7*1e3/freq) if (t<del+7*1e3/freq) {i8 =0}  else {if (8>nb_epsc)  {i8 = 0}  else { i8  = g*(1-exp(-(t-(del+7*1e3/freq))/tauOpf))*exp(-(t-(del+7*1e3/freq))/tauCpf)*(v - e) }}
    at_time(del+8*1e3/freq) if (t<del+8*1e3/freq) {i9 =0}  else {if (9>nb_epsc)  {i9 = 0}  else { i9  = g*(1-exp(-(t-(del+8*1e3/freq))/tauOpf))*exp(-(t-(del+8*1e3/freq))/tauCpf)*(v - e) }}
    at_time(del+9*1e3/freq) if (t<del+9*1e3/freq) {i10 =0} else {if (10>nb_epsc) {i10 = 0} else { i10 = g*(1-exp(-(t-(del+9*1e3/freq))/tauOpf))*exp(-(t-(del+9*1e3/freq))/tauCpf)*(v - e) }}
    i=i1+i2+i3+i4+i5+i6+i7+i8+i9+i10
}