TITLE Counter.mod  simple counter

COMMENT
 counter which can be run by a NetCon
 written by Jonathan Z Simon, jzsimon@isr.umd.edu
ENDCOMMENT

NEURON {
	POINT_PROCESS Counter
	RANGE count 
}

PARAMETER {
	count = 0
}

NET_RECEIVE(change) {
    count = count + change
}