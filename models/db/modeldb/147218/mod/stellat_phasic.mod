: phasic synaptic current of stellate cells

NEURON {
    SUFFIX stellate_phasic
    RANGE nb_ipsc,freq,del,tauOsc,tauCsc
    NONSPECIFIC_CURRENT  i
    RANGE i, e, g
}

PARAMETER {
    nb_ipsc = 6 : number of IPSCs 
    freq = 50 (/s) : frequency of IPSCs
    g = 0 (siemens/cm2)  < 0, 1e9 >
    e = -80    (millivolts)
    del = 100 (ms)
    tauOsc = 0.9 (ms)
    tauCsc = 9.0 (ms)
    i1 =0 (millivolts*siemens/cm2) 
    i2 =0 (millivolts*siemens/cm2)
    i3 =0 (millivolts*siemens/cm2)
    i4 =0 (millivolts*siemens/cm2)
    i5 =0 (millivolts*siemens/cm2)
    i6 =0 (millivolts*siemens/cm2)
}

ASSIGNED {
    i   (milliamp/cm2)
    v   (millivolt) 
}

INITIAL  { i = 0   }

BREAKPOINT {
    at_time(del)            if (t<del)            {i1 =0}  else {if (1>nb_ipsc)  {i1 = 0}  else  { i1  = g*(1-exp(-(t-(del))/tauOsc))*exp(-(t-(del))/tauCsc)*(v - e) }}
    at_time(del+1e3/freq)   if (t<del+1e3/freq)   {i2 =0}  else {if (2>nb_ipsc)  {i2 = 0}  else  { i2  = g*(1-exp(-(t-(del+1e3/freq))/tauOsc))  *exp(-(t-(del+1e3/freq))/tauCsc)*(v - e) }}
    at_time(del+2*1e3/freq) if (t<del+2*1e3/freq) {i3 =0}  else {if (3>nb_ipsc)  {i3 = 0}  else { i3  = g*(1-exp(-(t-(del+2*1e3/freq))/tauOsc))*exp(-(t-(del+2*1e3/freq))/tauCsc)*(v - e) }}
    at_time(del+3*1e3/freq) if (t<del+3*1e3/freq) {i4 =0}  else {if (4>nb_ipsc)  {i4 = 0}  else { i4  = g*(1-exp(-(t-(del+3*1e3/freq))/tauOsc))*exp(-(t-(del+3*1e3/freq))/tauCsc)*(v - e) }}
    at_time(del+4*1e3/freq) if (t<del+4*1e3/freq) {i5 =0}  else {if (5>nb_ipsc)  {i5 = 0}  else { i5  = g*(1-exp(-(t-(del+4*1e3/freq))/tauOsc))*exp(-(t-(del+4*1e3/freq))/tauCsc)*(v - e) }}
    at_time(del+5*1e3/freq) if (t<del+5*1e3/freq) {i6 =0}  else {if (6>nb_ipsc)  {i6 = 0}  else { i6  = g*(1-exp(-(t-(del+5*1e3/freq))/tauOsc))*exp(-(t-(del+5*1e3/freq))/tauCsc)*(v - e) }}
    i=i1+i2+i3+i4+i5+i6
}