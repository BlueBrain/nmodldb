: width = yvec.width(xvec, yval)  sensitive to noise in yvec
NEURON {
	SUFFIX nothing
}

VERBATIM

/*********************************************************************

Modified the linfit() function to be compatible with NEURON 7.1.
Only minor changes were necessary, in particular to change the 
linfit() variable 't' (used universally throughout NEURON) to 'tmp'.  

March 17, 2010
Christina Weaver
new contact email:  christina.weaver@fandm.edu

**********************************************************************/


extern double hoc_epsilon;
#define EPS hoc_epsilon
static double lin_interp();

/***********************************************************

Same as yfitness below, except that this function reports an
error proportional to the size of the experimental window.  
This is most useful in APs near the start or the end of a 
model run, where data in the width of the entire experimental 
window may not be available.  

    This function is called by the Multiple Run Fitter error functions that use AP shape.  The function is called by, e.g.:
	    eval = $o1.ywnscl_fitness($o2, peak, ytmp, xtmp,ytmp.firstpeak,tmp_modind,mnind,mxind)

    so that the variables defined below correspond to the following:

    	y		the y-values of model output
        x		the t-values of model output
        xpeak		time of the specified model AP peak
	yval		the y-values of the target output
	xval		the t-values of the target output
	val_pkind	index of the AP peak of the target data
        pkind		index of peak of the model AP whose shape error will be calculated
        mnind		first index of model data within time bounds
        mxind		last  index of model data within time bounds

***********************************************************/
static double ywnscl_fitness(void* vv) {
	int nx, ny, nyval, nxval, i, j;
	double sum, d, xpeak, *y, *x, *yval, *xval;
	int val_pkind, pkind, mod_start, mod_end, mnind, mxind, n_pts;
	double ytmp, val_winsz, mod_winsz;

	//
	// get parameters from NEURON function
	ny = vector_instance_px(vv, &y);		// y is the vector
	nx = vector_arg_px(1, &x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	xpeak = *getarg(2);
	nyval = vector_arg_px(3, &yval);
	nxval = vector_arg_px(4, &xval);
	val_pkind = *getarg(5);
	pkind = *getarg(6);
	mnind = *getarg(7);
	mxind = *getarg(8);

	j = 0;
	sum = 0.;
	n_pts = 0;

	//
	// Calculate the error directly at, before, and after the
	// peak, for each of the experimental data points included in (xval,yval).
	// Interpolate when data points do not occur at the same x-values.

        // start at the beginning of the model window, and the beginning 
        // of the experimental window. 
        mod_start = mxind+1;
        mod_end = mnind-1;
        i = mnind;
        j = 0; 

	while( i < mxind && j < nyval ) {

	    // model & experiment may have different (even variable) values of dt;
	    // find the appropriate x value, and interpolate accordingly.

	    // find next model value which is greater than the current experiment value
	    if( xval[j] > (x[i] - xpeak) ) {
	        while( (i<ny) && xval[j] > x[i] - xpeak && (fabs(xval[j] - (x[i] - xpeak)) > EPS) ) {
	            i++;
		}
	    }


	    if( fabs(xval[j] - (x[i]-xpeak)) < EPS ) {
	        // both model & experiment were evaluated at the same t value
	        d = y[i] - yval[j];
		sum += d*d;
                if( n_pts == 0 ) { mod_start = i; }
		n_pts++;
	    } else if((i == ny-1) && (j==nyval-1) && (xval[j] > x[i]-xpeak)) {
		// end of both experiment & model; interpolate experiment backward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
		d = y[i] - ytmp;
		sum += d*d;
                if( n_pts == 0 ) { mod_start = i; }
		n_pts++;
	    } else if((i==ny-1) && (j < nyval) ) {
		// j is not at end of list; interpolate experiment forward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
		d = y[i] - ytmp;
		sum += d*d;
                if( n_pts == 0 ) { mod_start = i; }
		n_pts++;
	    } else {
	        //interpolate model backward
	        ytmp = lin_interp(xval[j],x[i]-xpeak,y[i],x[i-1]-xpeak,y[i-1]);
	        d = ytmp - yval[j];
		sum += d*d;
                if( n_pts == 0 ) { mod_start = i; }
		n_pts++;
	    }

	    j++;
	}

	// end of model AP window
	mod_end = i;
	if( mod_end == nx ) { 
	     mod_end--; 
	}

	// root mean squared error
        if( n_pts == 0 ) { 
	    sum = -1;
        } else {
  	    sum = sqrt(sum/(double)n_pts);
	}

        // total size (in ms) of the experimental & model windows;
	// scale the calculated error proportional to the experimental error size
	val_winsz = xval[nxval-1] - xval[0];
        if( n_pts == 0 ) { 
	    mod_winsz = 0;
	} else {
	    mod_winsz = x[mod_end] - x[mod_start];
	}
        if( mod_winsz > 0 ) {
	    sum = mod_winsz*sum/val_winsz; 
	} else {
	    sum = -1.0;
	}
	return sum;

}


// Modified, Aug 2004 by Christina Weaver
static double yfitness_weaver(void* vv) {
	int nx, ny, nyval, nxval, i, j;
	double sum, d, xpeak, *y, *x, *yval, *xval;
	int val_pkind, pkind;
	double ytmp;
	ny = vector_instance_px(vv, &y);
	nx = vector_arg_px(1, &x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	xpeak = *getarg(2);
	nyval = vector_arg_px(3, &yval);
	nxval = vector_arg_px(4, &xval);
	val_pkind = *getarg(5);
	pkind = *getarg(6);
	j = 0;
	sum = 0.;

	// DON'T PENALIZE if the model AP window is not complete at start of the run,
	// or at end of the run.
	//
	// Calculate the error directly at, before, and after the
	// peak, for each of the experimental data points included in (xval,yval).
	// Interpolate when data points do not occur at the same x-values.

	// first, error at the peak
	d = y[pkind] - yval[val_pkind];
	sum += d*d;

	// now, error before the peak
        i = pkind; // -1;
	j = val_pkind -1;
	while( i > 0 && j >= 0 ) {

	    // model & experiment may have different (even variable) values of dt;
	    // find the appropriate x value, and interpolate accordingly.

	    // find next model value which is less than the current experiment value
	    if( xval[j] < x[i] - xpeak ) {
	        while( (i>0) && xval[j] < x[i] - xpeak ) {

	            i--;
		}
	    }


	    if( fabs(xval[j] - (x[i]-xpeak)) < EPS ) {
	        // both model & experiment were evaluated at the same t value
	        d = y[i] - yval[j];
		sum += d*d;
	    } else if( j==0 && i==0 && (xval[j] < x[i]-xpeak) ) {
		// interpolate experiment forward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
		d = y[i] - ytmp;
		sum += d*d;
	    } else if( i==0 && j>0 ) {
		// j is nonzero; interpolate experiment backward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
		d = y[i] - ytmp;
		sum += d*d;
	    } else {
	        // interpolate model forward
	        ytmp = lin_interp(xval[j],x[i]-xpeak,y[i],x[i+1]-xpeak,y[i+1]);
	        d = ytmp - yval[j];
		sum += d*d;
	    }

	    j--;
	}

	// error after the peak
        i = pkind; // + 1;
	j = val_pkind + 1;
	while( i < ny-1 && j < nyval ) {

	    // model & experiment may have different (even variable) values of dt;
	    // find the appropriate x value, and interpolate accordingly.

	    // find next model value which is greater than the current experiment value
	    if( xval[j] > x[i] - xpeak ) {
	        while( (i<ny) && xval[j] > x[i] - xpeak ) {
	            i++;
		}
	    }



	    if( fabs(xval[j] - (x[i]-xpeak)) < EPS ) {
	        // both model & experiment were evaluated at the same t value
	        d = y[i] - yval[j];
		sum += d*d;
	    } else if((i == ny-1) && (j==nyval-1) && (xval[j] > x[i]-xpeak)) {
		// end of both experiment & model; interpolate experiment backward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
		d = y[i] - ytmp;
		sum += d*d;
	    } else if((i==ny-1) && (j < nyval) ) {
		// j is not at end of list; interpolate experiment forward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
		d = y[i] - ytmp;
		sum += d*d;
	    } else {
	        //interpolate model backward
	        ytmp = lin_interp(xval[j],x[i]-xpeak,y[i],x[i-1]-xpeak,y[i-1]);
	        d = ytmp - yval[j];
		sum += d*d;
	    }

	    j++;
	}

	return sum;

}




static double lin_interp(xstar, x1, y1, x2, y2) double xstar, x1, y1, x2, y2; {
        double ystar;

	if( fabs(x2-x1) < EPS ) return 0.5*(y1+y2);

	ystar = y1 + ((y2-y1)/(x2-x1))*(xstar-x1);
	return ystar;
}

/************************************

    firstmax

    this function finds the first local max of the vector y.

************************************/
static double firstmax(void* vv) { 
	int ny, i;
	double *y;
	ny = vector_instance_px(vv, &y) - 1;
	i = 0;
	while (i < ny) {
		if (y[i] > y[i+1]) {
			return (double) i;
		}
		i = i + 1;
	}
	return 0.;
}


static double nextpeak(void* vv) {
	int ny, i;
	double *y;
	ny = vector_instance_px(vv, &y) - 1;
	i = *getarg(1);
	while (i < ny) {
		if (y[i] >= -20) {
			if (y[i] > y[i+1]) {
				return (double) i;
			}
			i = i + 1;
		} else {
			i = i + 2;
		}
	}
	return 0.;
}
 

static	double	sqrarg;

#define	SQR(a) (sqrarg=(a),sqrarg*sqrarg)

/*
*       From NUMERICAL RECIPES IN C.
*   
*	Given a set of points x[strt ... end], y[strt ... end], with standard
*	deviations sig[strt ... end], fit them to
*
*			y = a + bx		(!!!!!!!!!)
*
*	by minimizing chi_sq.
*	Returned are a,b, siga, sigb, chi_sq (chi2), and the goodness of fit
*	probability q (that the fit would have chi_sq this large or larger).
*	If mwt = 0 *	on input, then the standard deviations are assumed
*	unavailable, q is returned as 1.0, and the normalization of chi2
*	is to unit standard deviation on all points.
*/

void	linfit(void* vv) {

  /** mwt = 0, no 'sig' array is given , no q returned.  
      references to these values have been deleted. **/
	double	*x;
	double 	*y;
	double	*a; 
	double	*b;
	double 	*siga;
	double	*sigb;
	double 	*chi2;

	int	strt=0;
	int	end=0;

	int	i=0;
	int	nx=0;
	int	ny=0;

	double	wt=0.;
	double	tmp=0.;
	double 	sxoss=0.;
	double	ss=0.;
	double	sigdat=0.;
	double	sx  = 0.0;
	double	sy  = 0.0;
	double	st2 = 0.0;

	ny = vector_instance_px(vv,&y);
	nx = vector_instance_px(1,&x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	strt = *getarg(2);
	end  = *getarg(3);
	*a    = *getarg(4);
	*b    = *getarg(5);
	*siga = *getarg(6);
	*sigb = *getarg(7);
	*chi2 = *getarg(8);

	*b = 0.0;
	for( i = strt;  i <= end;  i++ )
	{
		sx += x[i];	sy += y[i];
	}
	ss = end+1-strt;

	sxoss = sx/ss;

	for( i = strt;  i <= end;  i++ )
	{
		tmp = x[i] - sxoss;
		st2 += tmp*tmp;
		*b += tmp*y[i];
	}

	*b /= st2;
	*a = (sy-sx*(*b))/ss;
	*siga = sqrt((1.0 + sx*sx/(ss*st2))/ss);
	*sigb = sqrt(1.0/st2);
	*chi2 = 0.0;

	for( i = strt;  i <= end;  i++ ) {
		*chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
	}
	sigdat = sqrt((*chi2)/(end-2));
	*siga *= sigdat;
	*sigb *= sigdat;

}
ENDVERBATIM


PROCEDURE install_weaver_fitness() {
VERBATIM
  {static int once; if (!once) { once = 1;
	install_vector_method("ywnscl_fitness", ywnscl_fitness);
	install_vector_method("yfitness_weaver", yfitness_weaver);
	install_vector_method("firstmax", firstmax);
	install_vector_method("nextpeak", nextpeak);
	install_vector_method("linfit", linfit);
  }}
ENDVERBATIM
}

