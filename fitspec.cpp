/*
 * fitspec.cpp
 *
 *  Created on: Mar 16, 2015
 *      Author: lerma
 */
#include <iostream>
#include <math.h>
#include <cstring>
#include <valarray>
#include <omp.h>
#include <time.h>
#include <fstream>
#include <vector>
#include <iterator>
#include <set>
#include <unistd.h>

#include <limits>

#include <CCfits/CCfits>

#include "fitspec.h"
#include "mpfit.h"
//#define THIRD_DIMENTION 240

using namespace std;
using namespace CCfits;



// Global variables********************
double piSqr = sqrt( 2.0 * M_PI );
size_t sliceSize,
	   aSize,
	   bSize;

valarray<double> finalCube,
				 sky,
				 wl,
				 var,
				 sig2d,
				 alast,
				 brlf,
				 aResult,
				 bResult;

vector<unsigned int> finalCubeDimentions,
					 skyDimentions,
					 wlDimentions;

//*************************************


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
	size_t iGlobal;	// Global index of the cube, from where to start computations and where to write the result
	size_t jGlobal;		//
	size_t jobSize;		// Size of the computations for thread
};


void* runThreads (void*);
int  triplet	(int, int, double*, double*, double**, void*);
int  tripletbr	(int, int, double*, double*, double**, void*);
void readImage	(valarray<double>&, vector<unsigned int>&, string);
int  writeImage	(string, long, vector<long>, valarray<double>);







/*
********* Function to read the input file 
*/
void readInputFile(string fileName, valarray<double> &alast, valarray<double> &brlf, size_t *start, size_t *sigStart, size_t *sliceSize, int *nbThreads) {
	ifstream myFile;
	myFile.open( fileName.c_str() );

	if ( myFile.is_open() ) {
		size_t alastSize, brlfSize;
		myFile >> (*nbThreads);
		myFile >> (*start);
		myFile >> (*sigStart);
		myFile >> (*sliceSize);
		myFile >> alastSize;
		alast.resize( alastSize );
		for(size_t i = 0; i < alastSize; ++i) {
			myFile >> alast[i];
		}
		myFile >> brlfSize;
		brlf.resize( brlfSize );
		for(size_t i = 0; i < brlfSize; ++i) {
			myFile >> brlf[i];
		}
//		cout << "sliceSize is: " << (*sliceSize) << endl;
	}
	else
		cout << "Can't open file: " << fileName << endl;
	myFile.close();
}


/*
********** Compute average value of the valarray
*/
double average(valarray<double> va) {
	return va.sum() / va.size();
}

/*
********** Compute standard deviation of the valarray
*/
double stdev(valarray<double> va)
{
    double E = 0;
    double ave = average( va );
    double inverse = 1.0 / static_cast<double>(va.size());
    for(size_t i = 0; i < va.size(); ++i) {
        E += pow( static_cast<double>( va[i] ) - ave, 2 );
    }
    return sqrt( inverse * E );
}



/*
********* Helper mfunction for quicksort
*/
int partition(valarray<double> &A, int lo, int hi) {
  int pivotIndex = ( hi + lo ) / 2;
  double pivotValue = A[pivotIndex];
  int storeIndex = lo;

  swap( A[pivotIndex], A[hi] );

  for(int i = lo; i < hi; ++i) {
    if ( A[i] < pivotValue ) {
      swap( A[i], A[storeIndex] );
      ++storeIndex;
    }
  }

  swap( A[storeIndex], A[hi] );

  return storeIndex;
}

/*
********** Sort the part of the valarray from *lo* to *hi* using quicksort algorithm
*/
void quickSort(valarray<double> &A, int lo, int hi) {
  if ( lo < hi ) {
    int p = partition( A, lo, hi );
    if ( ( p == (A.size() / 2) ) && ( (A.size() % 2) == 1 ) ) {
      cout << "returning" << endl;
      return;
    }
//#pragma omp sections
//    {
//#pragma omp section
  //  	{
    		quickSort( A, lo, p - 1 );
//    	}
//#pragma omp section
 //   	{
    		quickSort( A, p + 1, hi );
 //   	}
    }
//  }
}

/*
********* Sort the entire valarray using quicksort algorithm
*/
void quickSort(valarray<double> &A) {
  quickSort( A, 0, A.size() - 1 );
}




/*
********* Compute the median for valarray.
********* If number of elements is odd, median is a central element of the sorted values.
********* If it's even, than median is the average of two central values.
*/
double median(valarray<double> va) {

	size_t size = va.size();

	quickSort( va );

	return size % 2 ? va[size / 2] : ( va[size / 2] + va[size / 2 + 1] ) / 2;

}


/*
 ******** compute the average of 8 neighbors and the element itself in 3 dimensional array
 */
double avOfNeighbors(valarray<double> va, vector<unsigned int> dimentions,  size_t i, size_t j, size_t k) {

	double sum = 0.0;

	sum += va[i * dimentions[2] * dimentions[1] + dimentions[2] * j + k]
		 + va[i * dimentions[2] * dimentions[1] + dimentions[2] * ( j - 1 ) + k]
		 + va[i * dimentions[2] * dimentions[1] + dimentions[2] * ( j + 1 ) + k]
		 + va[( i - 1 ) * dimentions[2] * dimentions[1] + dimentions[2] * j + k]
		 + va[( i - 1 ) * dimentions[2] * dimentions[1] + dimentions[2] * ( j - 1 ) + k]
		 + va[( i - 1 ) * dimentions[2] * dimentions[1] + dimentions[2] * ( j + 1 ) + k]
		 + va[( i + 1 ) * dimentions[2] * dimentions[1] + dimentions[2] * j + k]
		 + va[( i + 1 ) * dimentions[2] * dimentions[1] + dimentions[2] * ( j - 1 ) + k]
		 + va[( i + 1 ) * dimentions[2] * dimentions[1] + dimentions[2] * ( j + 1 ) + k];

	return sum / 9.0;
}


/*
 ****** Resize one dimensional array, so that *startIndex* will be first element and (*startIndex*+*count*-1) will be last
 */
void resizeOneDimentional(valarray<double> &va, size_t startIndex, size_t count){
	valarray<double> tmp( count );
	tmp = va[slice( startIndex, count, 1 )];
	va.resize( count );
	va = tmp;
}


int main(int argc, char* argv[]) {
	clock_t ts = clock();

	cout << numeric_limits<double>::max() << endl;
//	return 0;


	//size_t startIndex = 0, sigStartIndex = 0;
	string inputFile = "inputfile.txt";

	if ( argc == 1 ) {
		cout << "You should provide input file name. \nIf not, the default will be used: inputfile.txt" << '\n';
	} else if ( argc == 2 ) {
			inputFile = argv[1];
		   } else
			   cout << "Only one parameter should be passed. \nThe rest will be ignored" << '\n';

/*
	valarray<double> finalCube,
					 sky,
					 wl,
					 var;

	vector<unsigned int> finalCubeDimentions,
						 skyDimentions,
						 wlDimentions;
	valarray<double> alast,
				   brlf;
*/
	size_t startIndex,
		   sigStartIndex;
//		   sliceSize;
	int nbThreads;

	readInputFile( "inputfile.txt", alast, brlf, &startIndex, &sigStartIndex, &sliceSize, &nbThreads );
	readImage( finalCube, finalCubeDimentions, "coadd-HaNII.fits" );
	readImage( sky, skyDimentions, "sky_HaNII.fits" );
	readImage( wl, wlDimentions, "wl_HaNII.fits" );

	//* Remove all elements except the range [startIndex, sliceSize-1]

/*
	for(size_t i = 0; i < finalCubeDimentions[0]; ++i)
		for(size_t j = 0; j < finalCubeDimentions[1]; ++j)
			for(size_t k = 373; k < 613; ++k) {
				size_t index = k * finalCubeDimentions[0] * finalCubeDimentions[1] + j * finalCubeDimentions[0] + i;
				cout << "finalCube[" << i << "][" << j << "][" << k << "]= " << finalCube[index] << endl;
			}
		*/
		//cout << finalCube[i + finalCubeDimentions[0]*finalCubeDimentions[1]] << "\t\t" << sky[i] << "\t\t" << wl[i] << endl;
		//cout << "finalCube[" << 0 << "][" << 0 << "][" << i << "]= " << finalCube[i * finalCubeDimentions[1] * finalCubeDimentions[2]] << endl;

	cout << "startIndex is: " << startIndex << endl;
	cout << "sigStartIndex is: " << sigStartIndex << endl;
	cout << "sky size is: " << sky.size() << endl;

	//cout << "finalCube: " << finalCube[slice( startIndex, sliceSize, 1 )] << endl;

	//valarray<double> tmpFinalCube( finalCubeDimentions[0] * finalCubeDimentions[1] * sliceSize );
/*	for(size_t i = 0; i < finalCubeDimentions[0]; ++i)
		for( size_t j = 0, k = 0; j < finalCubeDimentions[1]; ++j, k += sliceSize) {
			tmpFinalCube[slice( 0, sliceSize, 1 )] = finalCube[slice( startIndex, sliceSize, 1 )];
		}*/
//	tmpFinalCube[slice( 0, sliceSize, 1 )] = finalCube[slice( startIndex, sliceSize, 1 )];
/*
	for(size_t i = 0; i < sliceSize; ++i) {
		cout << "finalCube[" << i << "]= " << tmpFinalCube[i] << endl;
	}
*/


	resizeOneDimentional( sky, startIndex, sliceSize );
	resizeOneDimentional( wl, startIndex, sliceSize );
	//sky /= 10000;
	var.resize( sliceSize );
	var = 1. / sky;

	valarray<double> fitcube( finalCubeDimentions[0] * finalCubeDimentions[2] * sliceSize, 0.0 ),
					 acube( finalCubeDimentions[0] * finalCubeDimentions[2] * 5, 0.0 ),
					 bcube( finalCubeDimentions[0] * finalCubeDimentions[2] * 3, 0.0 ),
					 aerrorcube( finalCubeDimentions[0] * finalCubeDimentions[2] * 8, 0.0 ),
					 akeep( alast[0], alast.size() ), bkeep( brlf[0], brlf.size() );


	 sig2d.resize( finalCubeDimentions[0] * finalCubeDimentions[1] * sliceSize );
/*
	for(size_t outterIndex = 3; outterIndex < finalCubeDimentions[1] - 3; ++outterIndex)
		for(size_t innerIndex = 3; innerIndex < finalCubeDimentions[0] - 3; ++innerIndex) {


		}
*/
/*
	for(size_t i = 0; i < sky.size(); ++i)
		cout << "var[" << i << "]= " << var[i] << endl;
	for(size_t i = 0; i < alast.size(); ++i)
		cout << "alast[" << i << "]= " << alast[i] << endl;
	for(size_t i = 0; i < brlf.size(); ++i)
			cout << "brlf[" << i << "]= " << brlf[i] << endl;
*/
	cout << "The start index: " << startIndex << " The sig start index: " << sigStartIndex <<endl;

	cout << "The dimension of final cube is " << finalCubeDimentions.size() << endl;

	for(size_t i = 0; i < finalCubeDimentions.size(); ++i)
		cout << "The size of dimension " << i << " is " << finalCubeDimentions[i] << endl;
	//readTable();
	cout << "The size of final cube is " << finalCube.size() << endl;


	cout << "The size of sig2d is " << sig2d.size() << endl;
//	omp_set_num_threads( 4 );
//#pragma omp parallel for schedule(static) collapse(3)
	size_t first = startIndex,
		   last = sigStartIndex + sliceSize;
	if( sigStartIndex < startIndex ) {
		first = sigStartIndex;
		last = startIndex + sliceSize;
	}

	valarray<double> tmpFinalCube( finalCubeDimentions[0] * finalCubeDimentions[1] * sliceSize );

	for(size_t i = 0, l = 0, m = 0; i < finalCubeDimentions[0]; ++i)
		for(size_t j = 0; j < finalCubeDimentions[1]; ++j) {
			size_t  ij = j * finalCubeDimentions[0] + i;

			for(size_t k = first; k < last; ++k) {

				size_t ijk = k * finalCubeDimentions[0] * finalCubeDimentions[1] + ij;

				if ( k >= sigStartIndex && k < ( sigStartIndex + sliceSize ) ) {
					sig2d[l++] = finalCube[ijk];// / 10000;
					//slice( startIndex, sliceSize, 1 )
				}
				if ( k >= startIndex && k < ( startIndex + sliceSize ) ) {
					tmpFinalCube[m++] = finalCube[ijk];// / 10000;
					//cout << "finalCube[" << i << "][" << j << "][" << k << "]= " << finalCube[ijk] << '\t';
					//cout << "tmpFinalCube[" << i << "][" << j << "][" << k << "]= " << tmpFinalCube[m-1] << endl;
				}
				//cout << sig2d[index] << " ";
			}
//			cout << "l= " << l << " m= " << m << endl;
		}


	finalCube.resize( tmpFinalCube.size() );
	finalCubeDimentions[2] = sliceSize;
	finalCube = tmpFinalCube;
	tmpFinalCube.resize( 0 );
	cout << "The size of sig2d is " << sig2d.size() << endl;
	cout << "The size of tmpFinalCube is " << tmpFinalCube.size() << endl;
	cout << "The size of finalCube is " << finalCube.size() << endl;



/*
	for(size_t i = 0; i < finalCubeDimentions[0]; ++i)
		for(size_t j = 0; j < finalCubeDimentions[1]; ++j)
			for(size_t k = 0; k < finalCubeDimentions[2]; ++k) {
				//size_t index = k * finalCubeDimentions[0] * finalCubeDimentions[1] + j * finalCubeDimentions[0] + i;
				size_t index = i * finalCubeDimentions[2] * finalCubeDimentions[1] + finalCubeDimentions[2] * j + k;
				cout << "finalCube[" << i << "][" << j << "][" << k << "]= " << finalCube[index] << endl;
			}
*/
	//valarray<double> akeep ( alast ), bkeep( brlf );

//	double a[5] = { 0.0543693, 10.1, 2.0, 0.0, 10.1 };

	aSize = alast.size();
	bSize = brlf.size();
//	size_t aCount = 0,
//		   bCount = 0;

	vector<pthread_t> thread( nbThreads );
	vector<ParamsThread> paramsThreads( nbThreads ); //= new ParamsThread[nbThreads];

	for (size_t i = 0; i < nbThreads - 1; ++i) {
		paramsThreads[i].index = i;
		paramsThreads[i].jobSize = finalCubeDimentions[0] / nbThreads;
		paramsThreads[i].iGlobal = i * ( finalCubeDimentions[0] / nbThreads );
	}

	for (size_t i = 0; i < nbThreads - 1; ++i) {
		pthread_create( &thread[i], NULL, runThreads, (void*)&paramsThreads[i] );
	}

	for (size_t i = 0; i < nbThreads - 1; ++i) {
		pthread_join( thread[i], NULL );
	}



//	omp_set_num_threads( 2 );
//#pragma omp parallel
//	{

	aResult.resize( finalCubeDimentions[0] * finalCubeDimentions[1] * 5 );
	bResult.resize( finalCubeDimentions[0] * finalCubeDimentions[1] * 3 );


//	omp_set_num_threads( 2 );

//	omp_set_nested( 1 );
//#pragma omp parallel for schedule(static) collapse(2)


//	} End of parallel section
/*
	vector<long> aCubeDimentions( 3 );
	aCubeDimentions[0] = finalCubeDimentions[0] - 6;
	aCubeDimentions[1] = finalCubeDimentions[1] - 6;
	aCubeDimentions[2] = aSize;
*/
	//writeImage( "acube.fits", 3, aCubeDimentions, acube );


	cout << "Execution time is " << (double)(clock() - ts)/CLOCKS_PER_SEC << endl;
	return 0;
}


/*
 ******* Function to be executed at threads' launch
 */
void* runThreads(void* param) {

	ParamsThread *paramTh = (ParamsThread *) param;


	for(size_t j = 178; j < paramTh->jobSize/*finalCubeDimentions[1] - 3*/; ++j)
		for(size_t i = 206; i < 217/*finalCubeDimentions[0] - 3*/; ++i) {

//			cout << "Thread# " << omp_get_thread_num() << " entered in for" << endl;

//			cout << "nesed is active? " << omp_get_nested() << " nested level= " << omp_get_active_level() << endl;

			//valarray<double> a( alast.data(), alast.size() ), b( brlf.data(), brlf.size() );
			//double *a = alast.data(), *b = brlf.data();
			valarray<double> a( alast ), b( brlf );

			//doule a[5] = { 0.0543693, 10.1, 2.0, 0.0 10.1 }:
		//	size_t i = 35;
		//	j = 45;
			valarray<double> x( wl ), y( sliceSize ), sig( sliceSize );
			size_t  ij =  i * finalCubeDimentions[2] * finalCubeDimentions[1] + finalCubeDimentions[2] * j;

			for(size_t k = 0; k < sliceSize; ++k) {

				size_t ijk = k + ij;

				//cout << "finalCube[" << i << "][" << j << "][" << k << "]= " << finalCube[ijk] << endl;
				//cout << "finalCube[" << k << "]= " << finalCube[k] << endl;
				y[k] = finalCube[ijk];
				sig[k] = sig2d[ijk];
			}
			//exit(0);

			double sig2 = stdev( sig );
			sig2 *= sig2;

			valarray<double> w( sliceSize );
			w = var / median( var );

			double av = average( y );

/*			valarray<double> justFor( sliceSize );
			justFor = y - av;
			for(size_t ii = 0; ii < y.size(); ++ii)
				cout << "y[" << ii << "]= " << y[ii] << endl;
*/
			double chi0 = ( ( w * ( y - av ) * ( y - av ) ) / sig2 ).sum();
			//double chi0r = chi0 / y.size();
//			cout << "y average= " << average( y ) << endl;
//			cout << "av= " << av << endl;
/*
			for(size_t ind = 0; ind < w.size(); ++ind)
				cout << "w[" << ind << "]= " << w[ind] << endl;
*/
			struct Parameters p;
			p.x.resize( sliceSize );
			p.y.resize( sliceSize );
			p.yErr.resize( sliceSize );
			p.f.resize( sliceSize );
			p.x = x;
			p.y = y;
			p.yErr = w;
		//	p.par = alast.data();
			mp_result result;
			memset( &result, 0, sizeof( result ) );


//			for(size_t i = 0; i < aSize; ++i)
//				cout << " before a[" << i << "]= " << a[i] << endl;

			mpfit( tripletbr, sliceSize, aSize, &a[0], 0, 0, (void *) &p, &result );

			if( result.status < 0 ) {
				cout << "Error while fitting on index i= " << i << " j= " << j << " !" << endl;
				cout << "Error status: " << result.status << endl;
		//		continue;
			}

			double chi1 = ( w * ( y - p.f ) * ( y - p.f ) / sig2 ).sum();
			//double chi1r = chi1 / ( sliceSize - aSize + 1 );  // y.size() is replaced  with sliceSize as tey'r equal

//			double chi2lim = 49.0;

			if( chi0 - chi1 > 49.0 ) {

				mpfit( triplet, sliceSize, bSize, &b[0], 0, 0, (void *) &p, &result );

				double chi2 = ( w * ( y - p.f ) * ( y - p.f ) / sig2 ).sum();  // p.f is different from p.f used in chi1
				//double chi2r = chi1 / ( sliceSize - bSize + 1 );

				if ( chi1 - chi2 >= 64 ) {
					aResult[slice( paramTh->iGlobal, aSize, 1 )] = b[slice( 0, 5, 1 )]; // first 5 elements of b
					bResult[slice( paramTh->iGlobal, bSize, 1 )] = b[slice( 5, 3, 1 )]; // last 3 elements of b
	//				++bCount;
				} else
					aResult[slice( paramTh->iGlobal, aSize, 1 )] = a;

//				++aCount;
			} else {

				for(size_t k = 0; k < sliceSize; ++k) {
					y[k] = avOfNeighbors( finalCube, finalCubeDimentions, i, j, k );        // finalCube[ijk];
					sig[k] = avOfNeighbors( sig2d, finalCubeDimentions, i, j, k );              //sig2d[ijk];
				}

				sig2 = stdev( sig );
				sig2 *= sig2;
				av = average( y );
				chi0 = ( ( w * ( y - av ) * ( y - av ) ) / sig2 ).sum();
				//chi0r = chi0 / y.size();

				p.y = y;

				mpfit( tripletbr, sliceSize, aSize, &a[0], 0, 0, (void *) &p, &result );

				if( result.status < 0 ) {
					cout << "Error while fitting on index i= " << i << " j= " << j << " !" << endl;
					cout << "Error status: " << result.status << endl;
				//	continue;
				}

				chi1 = ( w * ( y - p.f ) * ( y - p.f ) / sig2 ).sum();
				//chi1r = chi1 / ( sliceSize - aSize + 1 );

				if( chi0 - chi1 > 49.0 ) {

					mpfit( triplet, sliceSize, bSize, &b[0], 0, 0, (void *) &p, &result );

					double chi2 = ( w * ( y - p.f ) * ( y - p.f ) / sig2 ).sum();  // p.f is different from p.f used in chi1
					//double chi2r = chi1 / ( sliceSize - bSize + 1 );

					if ( chi1 - chi2 >= 64 ) {
						aResult[slice( paramTh->iGlobal, aSize, 1 )] = b[slice( 0, 5, 1 )]; // first 5 elements of b
						bResult[slice( paramTh->iGlobal, bSize, 1 )] = b[slice( 5, 3, 1 )]; // last 3 elements of b
//						++bCount;
					} else
						aResult[slice( paramTh->iGlobal, aSize, 1 )] = a;

//					++aCount;
				} else
					cout << "No line detected at i= " << i << ", j= " << j << endl;
			}


/*
//			cout << "mpfit status: " << status << " or " << result.status << endl;
//			cout << "mpfit chi: " << result.bestnorm << endl;

			for(size_t i = 0; i < p.f.size(); ++i)
				cout << "F[" << i << "]= " << p.f[i] << endl;

			for(size_t i = 0; i < alast.size(); ++i)
				cout << "a[" << i << "]= " << a[i] << endl;
			exit( 0 );
*/
	} // End of double for

	return 0;
}



int triplet(int m, int n, double *p, double *deviates, double **derivs, void *data) {
	//double *par = ;
	struct Parameters *pr = (struct Parameters *) data;
	valarray<double> x( pr->x ),
					 y( pr->y ),
					 yErr( pr->yErr );
/*
	for(size_t i = 0; i < yErr.size(); ++i)
		cout << "yErr[" << i << "]= " << yErr[i] << endl;
*/
	double z = p[0],
		  inten = p[1],
		  ww = p[2],
		  contin = p[3],
		  NIIinten = p[4],
		  ha = 6562.8 * ( 1 + z ),
		  NII_r = 6583 * ( 1 + z ),
		  dNII_b = 6548.1 * ( 1 + z );

	//cout << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4] << endl;


	ww = sqrt( 0.6 * 0.6 + ww * ww );//         ; default width to match sky line.

	//cout << "ww= " << ww << endl;



	valarray<double> Ha = inten / ( ww * piSqr ) * exp( -0.5 * ( x - ha ) * ( x - ha ) / ( ww * ww ) );
	valarray<double> NII_a =  NIIinten / ( ww * piSqr ) * exp( -0.5 * ( x - NII_r ) * ( x - NII_r ) / ( ww * ww ) );
	valarray<double> NII_b = NIIinten / ( 3. * ww * piSqr ) * exp( -0.5 * ( x - dNII_b ) * ( x - dNII_b ) / ( ww * ww ) );

	pr->f = NII_b + Ha + NII_a  + contin;
	pr->f *= 10000;


	for( size_t i = 0; i < y.size(); ++i ) {
		deviates[i] = ( y[i] - pr->f[i] ) / yErr[i];
	}

	//pr->pder.resize( x.size() * par.size() );
/*
	for(size_t i = 0; i < pr->f.size(); ++i)
		cout << "pr->f[" << i << "]= " << pr->f[i] << endl;
	for(size_t i = 0; i < pr->pder.size(); ++i)
		cout << "pr->par[" << i << "]= " << pr->pder[i] << endl;
*/
	return 0;
	//pder = fltarr(n_elements(x),n_elements(par));//   ; no value returned.
}


void readImage(valarray<double> &contents, vector<unsigned int> &dimentions, string fileName) {

	  FITS::setVerboseMode(true);

	  try {

		   auto_ptr<FITS> pInfile( new FITS( fileName, Read, true ) );

		   PHDU& image = pInfile->pHDU();

		  // valarray<double>  contents;


		   // read all user-specifed, coordinate, and checksum keys in the image
		   image.readAllKeys();

		   image.read( contents );

		   int size =  image.axes();
		   dimentions.resize( size );
		  // dimentions[0] = size;
		   for(int i = 0; i < size; ++i )
			   dimentions[i] = image.axis( i );

		   cout << "size of the array is: " << contents.size() << '\n';

		   // this doesn't print the data, just header info.
	 //      std::cout << image << std::endl;
	/*       if ( size > 2 ) {
			   cout << fileName << ":\n";
			   for(size_t i = 0; i < dimentions[0] * 4; ++i)
				   cout << contents[i] << ' ';
			   cout << endl;
		   }
	*/
	/*
		   long ax1(image.axis(0));
		   long ax2(image.axis(1));


		   for (long j = 0; j < ax2; j+=10)
		   {
				   std::ostream_iterator<double> c(std::cout,"\t");
				   std::copy(&contents[j*ax1],&contents[(j+1)*ax1-1],c);
				   std::cout << '\n';
		   }
	*/
	 //      std::cout << std::endl;
	  //     return ;


	    } catch (FitsException&) {

	      cerr << " Fits Exception Thrown by FitsImage class \n";
	      cerr << " Fits file name : " << fileName << endl;
	     // rv=1; // problem
	    }
}


int tripletbr(int m, int n, double *p, double *deviates, double **derivs, void *data) {
	//double *par = ;
	struct Parameters *pr = (struct Parameters *) data;
	valarray<double> x( pr->x ),
					 y( pr->y ),
					 yErr( pr->yErr );

	double z = p[0],
		  inten = p[1],
		  ww = p[2],
		  contin = p[3],
		  NIIinten = p[4],
		  brint = p[5],
		  brww = p[6],
		  zbr = p[7],
		  ha = 6562.8 * ( 1 + z ),
		  dNII_r = 6583 * ( 1 + z ),
		  dNII_b = 6548.1 * ( 1 + z ),
		  dbr = 6562.8 * ( 1 + zbr );

	ww = sqrt( 0.6 * 0.6 + ww * ww ); //        ; default width to match sky line.
	brww = sqrt( 0.6 * 0.6 + brww * brww );


	//;if ww ge 10 then ww = 10

	valarray<double> Ha = inten / ( ww * piSqr ) *  exp( -0.5 * ( x - ha ) * ( x - ha ) / ( ww * ww ) );
	valarray<double> NII_r = NIIinten / ( ww * piSqr ) * exp( -0.5 * ( x - dNII_r ) * ( x - dNII_r ) / ( ww * ww ) );
	valarray<double> NII_b = NIIinten / ( 3. * ww * piSqr ) * exp( -0.5 * ( x - dNII_b ) * ( x - dNII_b ) / ( ww * ww ) );
	valarray<double> br = brint / ( brww * piSqr ) * exp( -0.5 * ( x - dbr ) * ( x - dbr ) / ( brww * brww ) );

	pr->f = Ha + NII_r + NII_b + br + contin;
	pr->f *= 10000;

	for( size_t i = 0; i < y.size(); ++i ) {
		deviates[i] = ( y[i] - pr->f[i] ) / yErr[i];
	}
	//pder = fltarr( n_elements( x ), n_elements( par ) );
	return 0;
}



int writeImage( string fileName, long naxis, vector<long> naxes, valarray<double> va )
{

    // Create a FITS primary array containing a 2-D image
    // declare axis arrays.
 //   long naxis    =   2;
  //  long naxes[2] = { 300, 200 };

    // declare auto-pointer to FITS at function scope. Ensures no resources
    // leaked if something fails in dynamic allocation.
    auto_ptr<FITS> pFits( 0 );

    try
    {
        // overwrite existing file if the file already exists.

//        const std::string fileName("!atestfil.fit");

        // Create a new FITS object, specifying the data type and axes for the primary
        // image. Simultaneously create the corresponding file.

        // this image is unsigned short data, demonstrating the cfitsio extension
        // to the FITS standard.

        pFits.reset( new FITS( fileName , static_cast<int>(DOUBLE_IMG) , naxis, naxes.data() ) );
    }
    catch (FITS::CantCreate*)
    {
          // ... or not, as the case may be.
          return -1;
    }

    long nelements = 1;

    for (long i = 0; i < naxis; ++i)
    	nelements *= naxes[i];


    // references for clarity.

 //   long& vectorLength = naxes[0];
 //   long& numberOfRows = naxes[1];
 //   long nelements(1);


    // Find the total size of the array.
    // this is a little fancier than necessary ( It's only
    // calculating naxes[0]*naxes[1]) but it demonstrates  use of the
    // C++ standard library accumulate algorithm.

 //   nelements = std::accumulate(&naxes[0],&naxes[naxis],1,std::multiplies<long>());

    // create a new image extension with a 300x300 array containing float data.
/*
    std::vector<long> extAx(2,300);
    string newName ("NEW-EXTENSION");
    ExtHDU* imageExt = pFits->addImage(newName,FLOAT_IMG,extAx);

    // create a dummy row with a ramp. Create an array and copy the row to
    // row-sized slices. [also demonstrates the use of valarray slices].
    // also demonstrate implicit type conversion when writing to the image:
    // input array will be of type float.

    std::valarray<int> row(vectorLength);
    for (long j = 0; j < vectorLength; ++j) row[j] = j;
    std::valarray<int> array(nelements);
    for (int i = 0; i < numberOfRows; ++i)
    {
        array[std::slice(vectorLength*static_cast<int>(i),vectorLength,1)] = row + i;
    }

    // create some data for the image extension.
    long extElements = std::accumulate(extAx.begin(),extAx.end(),1,std::multiplies<long>());
    std::valarray<float> ranData(extElements);
    const float PIBY (M_PI/150.);
    for ( int jj = 0 ; jj < extElements ; ++jj)
    {
            float arg = PIBY*jj;
            ranData[jj] = std::cos(arg);
    }

    long  fpixel(1);

    // write the image extension data: also demonstrates switching between
    // HDUs.
    imageExt->write(fpixel,extElements,ranData);

    //add two keys to the primary header, one long, one complex.

    long exposure(1500);
    std::complex<float> omega(std::cos(2*M_PI/3.),std::sin(2*M_PI/3));
    pFits->pHDU().addKey("EXPOSURE", exposure,"Total Exposure Time");
    pFits->pHDU().addKey("OMEGA",omega," Complex cube root of 1 ");

*/
    // The function PHDU& FITS::pHDU() returns a reference to the object representing
    // the primary HDU; PHDU::write( <args> ) is then used to write the data.

    pFits->pHDU().write( 1, nelements-1, va );


    // PHDU's friend ostream operator. Doesn't print the entire array, just the
    // required & user keywords, and is provided largely for testing purposes [see
    // readImage() for an example of how to output the image array to a stream].

  //  std::cout << pFits->pHDU() << std::endl;

    return 0;
}










