COMMENT
Mod file to find maximum voltage and dV/dt estimate
ENDCOMMENT

NEURON {
    SUFFIX vms
    RANGE vm, vs, vstemp, vlast
}

PARAMETER {
       dt (ms)
       v (millivolt)
}

ASSIGNED {
       vlast (millivolt)
       vm (millivolt)
       vs (millivolt/ms)
       vstemp (millivolt/ms)
}

INITIAL {
    vm = v
    vlast = v
    vs = 0 (millivolt/ms)
    vstemp = 0 (millivolt/ms)
}

BREAKPOINT { 
   if (v>vm) { vm=v }
   if (dt>0 (ms)) {
      vstemp=(v-vlast)/dt
      }
      if (vstemp>vs) { vs = vstemp }
      vlast=v
}




