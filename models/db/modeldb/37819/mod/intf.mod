: $Id: intf.mod,v 1.40 2004/05/11 22:45:50 billl Exp $

NEURON {
  ARTIFICIAL_CELL INTF
  RANGE tau1, tau2, tau3, tau4, refrac, m1, m2, m3, m4, thresh, refractory, fflag
  RANGE adap,adapwt,tauadap
}

PARAMETER {
  tau1 = 10 (ms)
  tau2 = 10 (ms)
  tau3 = 10 (ms)
  tau4 = 10 (ms)
  adapwt = 0
  tauadap= 10 (ms)
  refrac = 5 (ms)
  thresh = 1
  fflag           = 1             : don't change
}

ASSIGNED {
  m1
  m2
  m3
  m4
  adap
  t0(ms)
  refractory
}

INITIAL {
  adap = 0
  m1 = 0
  m2 = 0
  m3 = 0
  m4 = 0
  t0 = t
  refractory = 0 : 0-integrates input, 1-refractory
}

NET_RECEIVE (w1,w2,w3,w4) {
  INITIAL { w2=w2 w3=w3 w4=w4 }
  : always update the state vars
  m1 = m1*exp(-(t - t0)/tau1)
  m2 = m2*exp(-(t - t0)/tau2)
  m3 = m3*exp(-(t - t0)/tau3)
  m4 = m4*exp(-(t - t0)/tau4)
  adap = adap*exp(-(t - t0)/tauadap)
  t0 = t
  if (flag==0) { : only add weights if an external excitation
    m1 = m1 + w1
    m2 = m2 + w2
    m3 = m3 + w3
    m4 = m4 + w4
   }
  if (flag==1) { refractory = 0 }
  if ((refractory==0 || flag==1) && (m1+m2+m3+m4>thresh)) {
    refractory = 1
    adap = adap + adapwt
    net_send(refrac+adap*adap, refractory)
    net_event(t)
  }
}

FUNCTION M1() { M1 = m1*exp(-(t - t0)/tau1) }
FUNCTION M2() { M2 = m2*exp(-(t - t0)/tau2) }
FUNCTION M3() { M3 = m3*exp(-(t - t0)/tau3) }
FUNCTION M4() { M4 = m4*exp(-(t - t0)/tau4) }
FUNCTION AD() { AD = adap*exp(-(t - t0)/tauadap) }
FUNCTION MT() { 
 MT = m1*exp(-(t - t0)/tau1)+m2*exp(-(t - t0)/tau2)+m3*exp(-(t - t0)/tau3)+m4*exp(-(t - t0)/tau4) 
}
