TITLE ZAP current

COMMENT
-----------------------------------------------------------------------------

    ZAP current model for membrane impedance analysis
    ==================================================

 IMPLEMENTATION

  This mechanism is implemented as a nonspecific current defined as a
  point process, mimicking a current-clamp stimulation protocol, injecting
  a sinusoidally oscillating waveform I(t), with instantaneous frequency
  changing in time (i.e. chirp or ZAP waveform).
  
  I(t) = A * sin (2 pi (f(t) - Fstart) ( t - ttstart) / 2 )
  
  f(t) = Fstart + (Fstop - Fstart) * ( t - ttstart) / ( ttstop - ttstart)
  A(t) = Astart + (Astop - Astart) * ( t - ttstart) / ( ttstop - ttstart)
 
  Note: Although counterintuitive at a first glance, the above expression of I(t)
        indeed correspond to a sinusoid starting with the initial frequency Fstart
        that is linearly increasing up to Fstop while the time goes from ttstart to
        ttstop.

  Note: 
  Since this is an electrode current, positive values of i depolarize the cell and in the
  presence of the extracellular mechanism there will be a change in vext since i is not a
  transmembrane current but a current injected directly to the inside of the cell.
 
  Refer to: Cali' et al. (2007)

 PARAMETERS

  This mechanism takes the following parameters:

  Ioff      = 0. (nA) : initial current offset.
  Astart    = 0. (nA) : initial value of the (linearly changing) amplitude of the ZAP current.
  Astop     = 0. (nA) : final value of the (linearly changing) amplitude of the ZAP current.
  ttstart   = 0. (ms) : starting time of the stimulation.
  ttstop    = 0. (ms) : final time of the stimulation.
  Fstart    = 0. (Hz) : initial value of the (linearly changing) frequency of the ZAP current.
  Fstop     = 0. (Hz) : final value of the (linearly changing) frequency of the ZAP current.

 Written by M. Giugliano and C. Cali', Brain Mind Institute, EPFL, March 2006

-----------------------------------------------------------------------------
ENDCOMMENT


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS Izap
    RANGE Astart, Astop, ttstart, ttstop, Fstart, Fstop, Ioff
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp) 
    (mV) = (millivolt)
}

PARAMETER {
  Ioff      = 0. (nA) : initial current offset
  Astart    = 0. (nA) : initial value of the (linearly changing) amplitude of the ZAP current
  Astop     = 0. (nA) : final value of the (linearly changing) amplitude of the ZAP current
  ttstart   = 0. (ms) : starting time of the stimulation..
  ttstop    = 0. (ms) : final time of the stimulation..
  Fstart    = 0. (Hz) : initial value of the (linearly changing) frequency of the ZAP current
  Fstop     = 0. (Hz) : final value of the (linearly changing) frequency of the ZAP current
}

ASSIGNED {
    i     (nA)        : fluctuating current
}


BREAKPOINT {
    if ((t < ttstart) || (t > ttstop)) {  i = - Ioff }
    else { 
    i    = - (Ioff + (  ((Astart*(ttstop-t)/(ttstop-ttstart)) + ((t-ttstart)/(ttstop-ttstart)) * Astop)) * sin(  (0.0062831853071795866 * Fstart + ((t-ttstart)/(ttstop-ttstart)) * 0.5 * 0.0062831853071795866 * (Fstop - Fstart) ) * (t-ttstart)))        
    }
}


