NEURON {
	POINT_PROCESS STEP_REL
	RANGE GLU, release_time, amplitude, duration
}

PARAMETER {
    release_time = 30 (ms)
    duration = 1     (ms)
    amplitude = 1     (mM)
}

ASSIGNED{
        GLU (mM)
}

INITIAL {
	GLU = 0 (mM)
}

BREAKPOINT {
    if (t>release_time){
        if (t<release_time+duration){
            GLU=amplitude
        }
        else{
            GLU=0 (mM)
        }
    }

}
