COMMENT

NetStimG() is based on NetStim() but has a Gaussian stochastic interval. The Poisson noise capability is 

retained. New parametrs are MeanInterval and SD (standard deviation) which control the interval

Mike Sikora July,2006

Added-The Start time is also Gaussian randomized vars MeanStart + StartSD

ENDCOMMENT



NEURON	{ 

  ARTIFICIAL_CELL NetStimG

  RANGE interval, number, MeanStart, StartSD

  RANGE noise,MeanInterval,SD

}



PARAMETER {

        MeanInterval = 10 (ms)

        SD = 0.5 (ms)

        MeanStart = 50 (ms) : Start of 1st spike

        StartSD = 1 : Standard deviation of 1st spike

	number	= 10 <0,1e9>	: number of spikes (independent of noise)

	noise		= 0 <0,1>	: amount of randomness (0.0 - 1.0)

}



ASSIGNED {

	event (ms)

	on

	ispike

	interval

}



PROCEDURE seed(x) {

	set_seed(x)

}



INITIAL {

        interval = normrand(MeanInterval,SD)

	on = 0

	ispike = 0

	if (noise < 0) {

		noise = 0

	}

	if (noise > 1) {

		noise = 1

	}

	if (number > 0) {

		: allowing MeanStart to be <0

		event = normrand(MeanStart,StartSD)

		: 

		if (event < 0) {

			event = 0

		}

		net_send(event, 3)

	}

}	



PROCEDURE init_sequence(t(ms)) {

	if (number > 0) {

		on = 1

		event = 0

		ispike = 0

	}

}



FUNCTION invl(mean (ms)) (ms) {

	if (mean <= 0.) {

		mean = .01 (ms) : I would worry if it were 0.

	}

	if (noise == 0) {

		invl = mean

	}else{

		invl = (1. - noise)*mean + noise*mean*exprand(1)

	}

}



PROCEDURE next_invl() {

	if (number > 0) {

                interval=normrand(MeanInterval,SD)

		event = invl(interval)

	}

	if (ispike >= number) {

		on = 0

	}

}



NET_RECEIVE (w) {

	if (flag == 0) { : external event

		if (w > 0 && on == 0) { : turn on spike sequence

			init_sequence(t)

			: randomize the first spike so on average it occurs at

			: noise*interval (most likely interval is always 0)

			next_invl()

			event = event - interval*(1. - noise)

			net_send(event, 1)

		}else if (w < 0 && on == 1) { : turn off spiking

			on = 0

		}

	}

	if (flag == 3) { : from INITIAL

		if (on == 0) {

			init_sequence(t)

			net_send(0, 1)

		}

	}

	if (flag == 1 && on == 1) {

		ispike = ispike + 1

		net_event(t)

		next_invl()

		if (on == 1) {

			net_send(event, 1)

		}

	}

}



