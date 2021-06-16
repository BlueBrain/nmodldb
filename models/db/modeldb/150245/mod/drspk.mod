: $Id: drspk.mod,v 1.27 2011/11/22 04:23:40 samn Exp $

UNITS {
    (mV) = (millivolt)
    (nA) = (nanoamp)
}

NEURON {
  POINT_PROCESS DRSPK
  GLOBAL refrac, vrefrac, rdmthresh
  RANGE drive,rand,inrefrac,fflag,qq,thresh
  NONSPECIFIC_CURRENT i
}

PARAMETER {
  refrac = 5 (ms)
  vrefrac = 0 (mV)
  drive = 0
  rand = 1
  fflag=1
  i = 0 (nA)
  rdmthresh = 0
}

ASSIGNED {
  : i (nA)
  v (mV)
  inrefrac
  qq
  thresh
}

CONSTRUCTOR {
  VERBATIM 
  ENDVERBATIM
}

INITIAL {
  net_send(0, 3)
  i=0
  drive=0
  qq=0
  if (rdmthresh) {
    rand=1
    thresh=1
  }
  inrefrac=0
}

BREAKPOINT {
  if (inrefrac) {
    qq = 0
  } else {
    qq = drive
    thresh = rand
  }
  i = -qq
}

NET_RECEIVE(w) {
  if (flag == 1 && !inrefrac) {
    net_event(t)
    net_send(refrac, 2)
    v = vrefrac
    inrefrac=1
    qq = 0
  } else if (flag == 2) {
    inrefrac=0
    if(v > thresh) { net_send(0,1) }
  } else if (flag == 3) {
    WATCH (v>thresh) 1
  }	
}

:** vers gives version
PROCEDURE vers () {
  printf("$Id: drspk.mod,v 1.27 2011/11/22 04:23:40 samn Exp $\n")
}
