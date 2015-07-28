/*
 * fitspec.h
 *
 *  Created on: Mar 16, 2015
 *      Author: lerma
 */

#ifndef FITSPEC_H_
#define FITSPEC_H_

#define M_PI           3.14159265358979323846  /* pi */



/*
 ********* Parameters for the fitting function
 */
struct Parameters {
	valarray<double> x;
	valarray<double> f;
	valarray<double> y;
	valarray<double> yErr;
};

/*
 ********* Parameters for threads
 */
struct ParamsThread {
	size_t index;		// Index of thread
	size_t jGlobal;	// Global index of the cube, from where to start computations and where to write the result
	size_t jobSize;		// Size of the computations for thread
};


#endif /* FITSPEC_H_ */




