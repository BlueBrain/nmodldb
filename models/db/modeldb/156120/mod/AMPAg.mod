: ampa.mod
: saturating synapse model using discrete events

NEURON {
  POINT_PROCESS AMPA_S
  RANGE g, g_eff,g_specif,cellu_area
  RANGE Cdur, Alpha, Beta, Erev, Rinf, Rtau
  NONSPECIFIC_CURRENT i
}

UNITS {
  (nA)   = (nanoamp)
  (mV)   = (millivolt)
  (uS) 	= (microsiemens)
  (um)  = (microns)
}

PARAMETER {
  Cdur  = 1.0   (ms)  : transmitter duration (rising phase)
  Alpha = 1.1   (/ms) : forward (binding) rate
  Beta  = 0.19  (/ms) : backward (dissociation) rate
  Erev  = 0     (mV)  : equilibrium potential
  g_specif = 1e-5	(uS/um2)		: conductance
  cellu_area= 1   (um2)
  g                 (uS)
}

ASSIGNED {
  v     (mV)   : postsynaptic voltage
  i     (nA)   : current = g*(v - Erev)
  g_eff (uS)   : conductance
  Rtau  (ms)   : time constant of channel binding
  Rinf  : fraction of open channels if xmtr was present "forever"
  synon : sum of weights of all synapses that are in the "onset" state
}

STATE { Ron Roff }  : initialized to 0 by default
: Ron and Roff are the total conductances of all synapses
:   that are in the "onset" (transmitter pulse ON)
:   and "offset" (transmitter pulse OFF) states, respectively

INITIAL {
  Rinf = Alpha / (Alpha + Beta)
  Rtau = 1 / (Alpha + Beta)
  synon = 0
  g=g_specif*cellu_area
}

BREAKPOINT {
  SOLVE release METHOD cnexp
  g_eff = g*(Ron + Roff)
  i = g_eff*(v - Erev)
}

DERIVATIVE release {
  Ron' = (synon*Rinf - Ron)/Rtau
  Roff' = -Beta*Roff
}

NET_RECEIVE(weight, on, r0, t0 (ms)) {
  : on == 1 if transmitter is present ("onset" state), otherwise 0
  : flag is an implicit argument of NET_RECEIVE, normally 0
  if (flag == 0) {
    : a spike happened, so start onset state if not already in onset state
    if (!on) {
      : this synapse joins the set of synapses in the onset state
      synon = synon + weight
      r0 = r0*exp(-Beta*(t - t0)) : r0 at start of onset state
      :printf("r0aaa is %g\n", r0)
      : r0 joins the "onset" conductance pool,
      :   which grows according to Ron' = ...
      :   and leaves the "offset" conductance pool,
      :   which decays according to Roff' = ...
      Ron = Ron + r0
      Roff = Roff - r0
      t0 = t
      on = 1
      net_send(Cdur, 1)
    } else {
      : already in onset state, so move offset time
      net_move(t+Cdur)
    }
  }
  if (flag == 1) {
    : "turn off transmitter"
    : i.e. this synapse joins the set of synapses in the offset state
    synon = synon - weight
    : r0 at start of offset state
    :printf("Ron and Roff are %g\t%g\n", Ron,Roff)
    r0 = weight*Rinf + (r0 - weight*Rinf)*exp(-(t - t0)/Rtau)
    :printf("r0 is %g\n", r0)
    : r0 leaves the "onset" conductance pool,
    :   and joins the "offset" conductance pool
    Ron = Ron - r0
    Roff = Roff + r0
    t0 = t
    on = 0
  }
  :printf("the conductance is %g\n",g_eff )
}
