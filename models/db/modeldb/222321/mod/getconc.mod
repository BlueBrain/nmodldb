NEURON {
	SUFFIX getconc
	USEION na WRITE nai, nao
	USEION  k WRITE  ki,  ko
	USEION ca WRITE cai, cao
	USEION cl WRITE cli, clo
	POINTER naip, naop, kip, kop, caip, caop, clip, clop
}

ASSIGNED {
	naip naop kip kop caip caop clip clop
}

STATE { nai nao ki ko cai cao cli clo }

BREAKPOINT {
	  nai=naip nao=naop
	  ki=kip   ko=kop
	  cai=caip cao=caop
          cli=clip clo=clop
}

INITIAL {
	at_time(t) {
	  nai=naip nao=naop
          ki=kip   ko=kop
	  cai=caip cao=caop
	  cli=clip clo=clop
	}
}