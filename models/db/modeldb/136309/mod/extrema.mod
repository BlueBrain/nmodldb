NEURON {
        SUFFIX extr     : for "extrema"
        RANGE vmax, vmin, tmax, tmin, r2, vmax2, tmin2
}

ASSIGNED {
        v (millivolt)
        vmin (millivolt)
        tmin (ms)
        vmax (millivolt)
        tmax (ms)
        vmax2 (millivolt)
        tmax2 (ms)
        r2 (1)
}

INITIAL {
        vmin = v
        tmin = t
        vmax = v
        tmax = t
        vmax2 = v
        tmax2 = t
        r2 = 0
}

AFTER SOLVE {
        if (v < vmin) {
                vmin = v
                tmin = t
        }
        if (v > vmax) {
                vmax = v
                tmax = t
        }
        if(r2 == 0 && vmax > -70 && v < -70){
          r2 = 1
        }
        if(r2 && v>vmax2){
          vmax2 = v
          tmax2 = t
        }
}
