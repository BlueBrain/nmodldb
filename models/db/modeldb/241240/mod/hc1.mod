COMMENT

Hair cell version 1

Fires as Poisson process with rate constant tau

ENDCOMMENT

VERBATIM
extern double du_dev0( );
extern double dexp_dev( );
ENDVERBATIM

NEURON {
  ARTIFICIAL_CELL 	hc1
  RANGE			tau, xp, yp, zp, del
}

PARAMETER {
  idebug	= 0
  tau		= 200 (ms)
  xp		= 0	: x,y location
  yp		= 0
  zp		= 0
  del		= 0
}

ASSIGNED {

}

INITIAL {
  if( del > 0 ) {
     net_send( del,3)
  }
  net_send( 0, 1 )
}

NET_RECEIVE( weight ){
  LOCAL wait
  if( idebug ) {
    printf( "Net_receive hc1: t %g flag %g", t, flag )
    if( flag == 0 ){ printf( " weight %g", weight ) }
    printf( "\n" )
  }
  if( flag == 1 && tau > 0 ){
    wait = dexp_dev( 1/tau )
    net_send( wait, 2 )
  }
  if( flag == 2 ){
    net_event( t )
    wait = dexp_dev( 1/tau )
    net_send( wait, 2 )
  }
  if( flag == 3 ){	: self-event from del > 0
    net_event(t)
:   wait=del
:    net_send(wait,2)
  }
}
