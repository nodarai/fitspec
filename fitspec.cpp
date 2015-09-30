/*
 * fitspec.cpp
 *
 *  Created on: Mar 16, 2015
 *      Author: lerma
 */
#include <iostream>
#include <math.h>
#include <cmath>
#include <cstring>
#include <valarray>
#include <time.h>
#include <fstream>
#include <vector>
#include <iterator>
#include <set>
#include <unistd.h>
#include <cstdio>
#include <assert.h>
#include <map>
#include <new>

#include <CCfits/CCfits>
#include <boost/math/distributions/fisher_f.hpp>

#include "functions.h"
#include "mpfit.h"
#include "cJSON/cJSON.h"

//#define INFO()
typedef valarray<double> vad;

using namespace std;
using namespace CCfits;
using boost::math::fisher_f;





// Global variables********************

map<string, int (*)(int, int, double*, double*, double**, void*)> fittingFunctions;

size_t sliceSize,
	   aSize,
	   bSize;

int info,
	alwaysFg,
	compare_limits,
	compare_guess;

double alpha,
	   chi2lim,
	   chi2lim2;

string fittingFunction;

valarray<double> finalCube,
				 sky,
				 wl,
				 var,
				 sig2d,
				 alast,
				 brlf,
				 aResult,
				 bResult,
				 cResult,
				 limits;

vector<unsigned int> finalCubeDimentions,
					 skyDimentions,
					 wlDimentions;

//*************************************




void*  runThreads   (void*);
void   readImage	(valarray<double>&, vector<unsigned int>&, string);
int    writeImage	(string, long, vector<long>, valarray<double>);
void   printArray	(size_t, valarray<double>, string);
void   printArray	(size_t, vector<double>, string);
double chiSquared   (valarray<double>, valarray<double>);
double fTest	    (size_t, size_t, double, double);
int 	oneEqual	(valarray<double> a, valarray<double> b);
vad vecToValar ( vector<double> );




cJSON* checkParam( cJSON* p, int type ) {
	for(cJSON* iter = p->child; iter; iter = iter->next) {
		if( !strcmp( iter->string, "value" ) ) {
			if( iter->type == type ) {
				return iter;
			}
		}
	}
	return 0;
}


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
	}
	else
		cout << "Can't open file: " << fileName << endl;
	myFile.close();
}


vector<double> readArrayJSON(cJSON **iter) {

	cJSON* tmpIter = *iter;
	vector<double> v;

	while( strcmp( tmpIter->string, "value" ) ) {    //find the "value" parameter
		tmpIter = tmpIter->next;
	}

	if( tmpIter->type != cJSON_Array )
		cout << "Error in the input file. Wrong type of value parameter in pads: " << tmpIter->type << endl;

	for( cJSON* it = tmpIter->child; it; it = it->next ) {
		v.push_back( it->valuedouble );
	}
	return v;
}


int readInputFileJSON(string fileName, string &resultfile, vector<string> &inputFits,
					  	 size_t *start, size_t *sigStart, size_t *nbThreads) {

	ifstream in( fileName.c_str() );

	if ( !in.is_open() ) {
		cout << "Unable to open the file: " << fileName << endl;
		return 1;
	}

	string contents( ( istreambuf_iterator<char>( in ) ), istreambuf_iterator<char>() );
	char* inputString = (char*)contents.c_str();

	cJSON* inputParams;
	inputParams = cJSON_Parse( inputString );

	assert( inputParams );

	cJSON* iter = inputParams->child;
	size_t i;

	for( i = 0; iter; iter = iter->next, ++i ) {

		if( !strcmp( iter->string, "number_of_threads" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_Number );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				(*nbThreads) = tmpIter->valueint;
			}
		}
		else
		if( !strcmp( iter->string, "start_index" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_Number );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				(*start) = tmpIter->valueint;
			}
		}
		else
		if( !strcmp( iter->string, "sig_start_index" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_Number );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				(*sigStart) = tmpIter->valueint;
			}
		}
		else
		if( !strcmp( iter->string, "inputs_limits" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_Array );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
					vector<double> v =  readArrayJSON( &iter->child );
					valarray<double> vaTmp;
					vaTmp.resize(v.size());
					vaTmp = vecToValar(v);
					limits.resize( vaTmp.size() );
					limits = vaTmp;
				}
		}
		else
		if( !strcmp( iter->string, "slice_size" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_Number );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				sliceSize = tmpIter->valueint;
			}
		}
		else
		if( !strcmp( iter->string, "alpha" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_Number );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				alpha = tmpIter->valuedouble;
			}
		}
		else
		if( !strcmp( iter->string, "always_fg" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_Number );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				alwaysFg = tmpIter->valueint;
			}
		}
		else
		if( !strcmp( iter->string, "compare_guess" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_Number );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				compare_guess = tmpIter->valueint;
			}
		}
		else
		if( !strcmp( iter->string, "compare_limits" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_Number );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				compare_limits = tmpIter->valueint;
			}
		}
		else
		if( !strcmp( iter->string, "chi2lim" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_Number );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				chi2lim = tmpIter->valuedouble;
			}
		}
		else
		if( !strcmp( iter->string, "chi2lim2" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_Number );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				chi2lim2 = tmpIter->valuedouble;
			}
		}
		else
		if( !strcmp( iter->string, "finalcube" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_String );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				inputFits[0] = tmpIter->valuestring;
			}
		}
		else
		if( !strcmp( iter->string, "fitting_function" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_String );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				fittingFunction = tmpIter->valuestring;
				cout << "Fitting function is:    " << fittingFunction << endl;
			}
		}
		else
		if( !strcmp( iter->string, "sky" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_String );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				inputFits[1] = tmpIter->valuestring;
			}
		}
		else
		if( !strcmp( iter->string, "wl" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_String );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				inputFits[2] = tmpIter->valuestring;
			}
		}
		else
		if( !strcmp( iter->string, "resultfile" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_String );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
				resultfile = tmpIter->valuestring;
			}
		}
		else
		if( !strcmp( iter->string, "a_first_guess" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_Array );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
					vector<double> v =  readArrayJSON( &iter->child );
					valarray<double> vaTmp;
					vaTmp.resize(v.size());
					vaTmp = vecToValar(v);
					alast.resize( vaTmp.size() );
					alast = vaTmp;
				}
		}
		else
		if( !strcmp( iter->string, "b_first_guess" ) ) {
			cJSON* tmpIter = checkParam( iter, cJSON_Array );
			if( !tmpIter ) {
				printf("Wrong parameter type for %s\n",  iter->string);
				return 1;
			}
			else {
					vector<double> v =  readArrayJSON( &iter->child );
					valarray<double> vaTmp;
					vaTmp.resize(v.size());
					vaTmp = vecToValar(v);
					brlf.resize( vaTmp.size() );
					brlf = vaTmp;
				}
		}
	}

	in.close();

	return 0;
}


valarray<double> vecToValar( vector<double> v ) {
	valarray<double> va;
	va.resize( v.size() );
	for(size_t i=0; i< v.size(); i++) {
		va[i] = v[i];
	}
	return va;
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
double stdev(valarray<double> va) {
    double E = 0;
    double ave = average( va );
    double inverse = 1.0 / static_cast<double>(va.size());

    E = ( ( va - ave ) * ( va - ave ) ).sum();
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
    size_t p = partition( A, lo, hi );
    if ( ( p == ( A.size() / 2 ) ) && ( ( A.size() % 2 ) == 1 ) ) {
      return;
    }

	quickSort( A, lo, p - 1 );
	quickSort( A, p + 1, hi );
    }
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

void avOfNeighbors (valarray<double> &y, bool fc, size_t i, size_t j, size_t neibSize) {

	double* va = static_cast<double*>(  fc ? &finalCube[0] : &sig2d[0] );

	double sum = 0.0;

	size_t dim0 = finalCubeDimentions[0];
	size_t dim01 = finalCubeDimentions[0] * finalCubeDimentions[1];

	for(size_t k = 0; k < sliceSize; ++k) {
		for(size_t l = j - neibSize / 2; l <= j + neibSize / 2; ++l) {
			size_t ij = l * dim0 + k * dim01;
			for(size_t m = i - neibSize / 2; m <= i + neibSize / 2; ++m) {
				size_t ijk = m + ij;
				sum+= *(va + ijk);

			}
		}
		y[k] = sum / 9.0;
		sum = 0.0;
	}
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

void printArray(size_t index, valarray<double> a, string name) {

	cout << "Thread #" << index << " " << name << " = { ";
	for(size_t i = 0; i < a.size(); ++i)
		cout << a[i] << ", ";
	cout << "}" << endl;
}



void printArray(size_t index, vector<double> a, string name) {

	cout << "Thread #" << index << " " << name << " = { ";
	for(size_t i = 0; i < a.size(); ++i)
		cout << a[i] << ", ";
	cout << "}" << endl;
}



int main(int argc, char* argv[]) {

	time_t ts = time( NULL );
	time( &ts );

	cout.rdbuf()->pubsetbuf( 0, 0 ); //Turn off buffered output



	//Mapping the strings to functions
	fittingFunctions["tripletHb"] 	 = tripletHb;
	fittingFunctions["tripletHb_br"] = tripletHb_br;
	fittingFunctions["triplet"] 	 = triplet;
	fittingFunctions["triplet_br"] 	 = triplet_br;
	fittingFunctions["SII"] 		 = SII;
	fittingFunctions["SII_br"] 		 = SII_br;
	fittingFunctions["Hb"] 			 = Hb;
	fittingFunctions["Hb_br"] 		 = Hb_br;
	fittingFunctions["Ha"] 			 = Ha;
	fittingFunctions["Ha_br"] 		 = Ha_br;
	fittingFunctions["OI"] 			 = OI;
	fittingFunctions["OI_br"] 		 = OI_br;
	fittingFunctions["Cont"] 		 = Cont;


	string inputFile = "inputfile.json",
		   resultFile;
	/*This vector contains filenames of FITS input data
	  0: finalcube, 1: sky, 2: wl */
	vector<string> inputFits( 3 );
	info = 0;

	switch( argc ) {
	case 1:
		cout << "You should provide input file name. \nIf not, the default will be used: inputfile.json" << endl;
		break;
	case 2:
		if( !strcmp( argv[1], "-info" ) ) {
			cout << "You should provide input file name. \nIf not, the default will be used: inputfile.json" << endl;
			info = 1;
		} else
			inputFile = argv[1];
		break;
	case 3:
		if( !strcmp( argv[1], "-info" ) ) {
			inputFile = argv[2];
			info = 1;
		} else
			if( !strcmp( argv[2], "-info" ) ) {
				inputFile = argv[1];
				info = 1;
			}
		break;
	default:
		 cout << "At most two parameters should be passed. \n1)Input file name \n2) -info if you want informative output " << endl;
	}


	size_t startIndex,
		   sigStartIndex,
		   nbThreads;

	readInputFileJSON( inputFile, resultFile, inputFits, &startIndex, &sigStartIndex, &nbThreads );

	readImage( finalCube, finalCubeDimentions, inputFits[0] );
	readImage( sky, skyDimentions, inputFits[1] );
	readImage( wl, wlDimentions, inputFits[2] );

	cout << "resultfile= " << resultFile << endl;
	cout << "input fits= " << inputFits[0] << endl;
	cout << "startIndex is: " << startIndex << endl;
	cout << "sigStartIndex is: " << sigStartIndex << endl;
	cout << "sky size is: " << sky.size() << endl;


	resizeOneDimentional( sky, startIndex, sliceSize );
	resizeOneDimentional( wl, startIndex, sliceSize );
	var.resize( sliceSize );
	var = 1. / sky;

	valarray<double> fitcube( finalCubeDimentions[0] * finalCubeDimentions[2] * sliceSize, 0.0 ),
					 acube( finalCubeDimentions[0] * finalCubeDimentions[2] * 5, 0.0 ),
					 bcube( finalCubeDimentions[0] * finalCubeDimentions[2] * 3, 0.0 ),
					 aerrorcube( finalCubeDimentions[0] * finalCubeDimentions[2] * 8, 0.0 ),
					 akeep( alast[0], alast.size() ), bkeep( brlf[0], brlf.size() );


	sig2d.resize( finalCubeDimentions[0] * finalCubeDimentions[1] * sliceSize );

	cout << "The start index: " << startIndex << " The sig start index: " << sigStartIndex <<endl;

	cout << "The dimension of final cube is " << finalCubeDimentions.size() << endl;

	for(size_t i = 0; i < finalCubeDimentions.size(); ++i)
		cout << "The size of dimension " << i << " is " << finalCubeDimentions[i] << endl;

	cout << "The size of final cube is " << finalCube.size() << endl;
	cout << "The size of sig2d is " << sig2d.size() << endl;



	valarray<double> tmpFinalCube( finalCubeDimentions[0] * finalCubeDimentions[1] * sliceSize );

	size_t fcd01 = finalCubeDimentions[0] * finalCubeDimentions[1];

	cout << "first index of new final cube = " << fcd01 * startIndex << endl;
	cout << "last index of new final cube = " << fcd01 * startIndex + fcd01 * sliceSize << endl;


	tmpFinalCube = finalCube[slice( fcd01 * startIndex, fcd01 * sliceSize, 1 )];
	sig2d = finalCube[slice( fcd01 * sigStartIndex, fcd01 * sliceSize, 1 )];

	finalCube.resize( tmpFinalCube.size() );
	finalCubeDimentions[2] = sliceSize;
	finalCube = tmpFinalCube;
	tmpFinalCube.resize( 0 );


	aSize = alast.size();
	bSize = brlf.size();

	vector<long> aResDims( 3 );
	aResDims[0] = finalCubeDimentions[0];
	aResDims[1] = finalCubeDimentions[1];
	aResDims[2] = aSize;

	vector<long> bResDims( 3 );
	bResDims[0] = finalCubeDimentions[0];
	bResDims[1] = finalCubeDimentions[1];
	bResDims[2] = bSize - aSize;

	vector<long> cResDims( 3 );
	cResDims[0] = finalCubeDimentions[0];
	cResDims[1] = finalCubeDimentions[1];
	cResDims[2] = 2;

	aResult.resize( aResDims[0] * aResDims[1] * aResDims[2] );
	bResult.resize( bResDims[0] * bResDims[1] * bResDims[2] );
	cResult.resize( cResDims[0] * cResDims[1] * cResDims[2] );

	vector<pthread_t> thread( nbThreads );
	vector<ParamsThread> paramsThreads( nbThreads );

	size_t jobSize = ( finalCubeDimentions[1] - 6 ) / nbThreads;
	size_t remainder = ( finalCubeDimentions[1] - 6 ) % nbThreads;
	paramsThreads[0].index = 0;
	paramsThreads[0].jobSize = 0 < remainder ? jobSize + 1 : jobSize;            //Calculations are made from index 3 to  *finalCubeDimentions[0]*-3
	paramsThreads[0].jGlobal = 3;

	for (size_t i = 1; i < nbThreads; ++i) {
		paramsThreads[i].index = i;
		paramsThreads[i].jobSize = i < remainder ? jobSize + 1 : jobSize;            //Calculations are made from index 3 to  *finalCubeDimentions[0]*-3
		paramsThreads[i].jGlobal = paramsThreads[i - 1].jGlobal + paramsThreads[i - 1].jobSize;	  //Calculations are made from index 3 to  *finalCubeDimentions[0]*-3
	}


	for (size_t i = 0; i < remainder; ++i) {
		++paramsThreads[i].jobSize;
		if( i < ( paramsThreads.size() - 1 ) )
			++paramsThreads[i + 1].jGlobal;
	}

	for (size_t i = 0; i < nbThreads; ++i) {
		pthread_create( &thread[i], NULL, runThreads, (void*)&paramsThreads[i] );
	}

	for (size_t i = 0; i < nbThreads; ++i) {
		pthread_join( thread[i], NULL );
	}



	writeImage( "a" + resultFile, 3, aResDims, aResult );
	writeImage( "b" + resultFile, 3, bResDims, bResult );
	writeImage( "c" + resultFile, 3, cResDims, cResult );

	cout << "Execution time is " <<  time( NULL ) - ts << " seconds." << endl;
	return 0;
}


/*
 ******* Function to be executed at threads' launch
 */
void* runThreads(void* param) {


	ParamsThread *paramTh = (ParamsThread *) param;

	cout << "Thread #" << paramTh->index << " started working." << endl;


	size_t fcDim01 = finalCubeDimentions[0] * finalCubeDimentions[1];
	size_t resDim01 = ( finalCubeDimentions[0] ) * ( finalCubeDimentions[1]  );

	valarray<double> x( wl ), y( sliceSize ), sig( sliceSize );
	valarray<double> a( alast ), b( brlf );
	valarray<double> a_old( alast ), b_old( brlf );
	valarray<double> w( sliceSize );
	w = var / median( var );

	if( info ) {
		printArray( paramTh->index, alast, "alast" );
		printArray( paramTh->index, a, "first_a" );
	}

	mp_par pars[b.size()];
	memset( pars, 0, sizeof( pars ) );

	for(size_t ind = 0; ind < limits.size(); ind += 3) {
		pars[(int)limits[ind]].limited[0] = true;
		pars[(int)limits[ind]].limited[1] = true;
		pars[(int)limits[ind]].limits[0] = limits[ind + 1];
		pars[(int)limits[ind]].limits[1] = limits[ind + 2];
		pars[(int)limits[ind]].fixed = 0;
	}


	for(size_t j = 0; ( j < paramTh->jobSize ) && ( j + paramTh->jGlobal < finalCubeDimentions[1] - 3 ); ++j) {

		size_t realInd = j + paramTh->jGlobal; // j index of the finalcube
		a_old = alast;
		b_old = brlf;


		for(size_t i = 3; i < finalCubeDimentions[0] - 3; ++i) {

			size_t ij = realInd * finalCubeDimentions[0] + i;
			size_t resIj = realInd * ( finalCubeDimentions[0] ) + i;

			for(size_t k = 0; k < sliceSize; ++k) {

				size_t ijk = k * fcDim01 + ij;
				y[k] = finalCube[ijk];
				sig[k] = sig2d[ijk];
			}

			double sig2 = stdev( sig );
			sig2 *= sig2;

			if ( sig2 == 0.0 )
				continue;

			double av = average( y );
			double chi0 = ( ( w * ( y - av ) * ( y - av ) ) / sig2 ).sum();

			struct Parameters p;
			p.x.resize( sliceSize );
			p.y.resize( sliceSize );
			p.yErr.resize( sliceSize );
			p.f.resize( sliceSize );
			p.x = x;
			p.y = y;
			p.yErr = 1. / w;

			mp_result result;
			memset( &result, 0, sizeof( result ) );


			a = a_old;
			b = b_old;

			if( info ) {
				cout << "At i=" << i << ", j=" << realInd << " ";
				printArray( paramTh->index, a, "aBefore" );
			}

			mpfit( fittingFunctions[fittingFunction], sliceSize, aSize, &a[0], pars, 0, (void *) &p, &result );

			if( result.status < 0 ) {
				cout << "Error while fitting on index i= " << realInd << " j= " << j << " !" << endl;
				cout << "Error status: " << result.status << endl;
			}

			if( info ) {
				cout << "At i=" << i << ", j=" << realInd << " ";
				printArray( paramTh->index, a, "a" );
			}

			double chi1 = ( w * ( y - p.f ) * ( y - p.f ) / sig2 ).sum();
			double bestNorm = chiSquared( y, p.f );
			size_t N0 = y.size() - 1,
				   N1 = y.size() - a.size(),
				   N2 = y.size() - b.size();

			double pFTest = fTest( N0, N1, chiSquared( y, sig ), bestNorm );
			double pFTest_old;

			if( ( chi0 - chi1 > chi2lim ) && pFTest < alpha && !oneEqual( a, alast) ) { // p < alpha
				pFTest_old = pFTest;

				mpfit( fittingFunctions[fittingFunction + "_br"], sliceSize, bSize, &b[0], pars, 0, (void *) &p, &result );

				if( info ) {
					cout << "At i=" << i << ", j=" << realInd << " ";
					printArray( paramTh->index, b, "b" );
				}

				double chi2 = ( w * ( y - p.f ) * ( y - p.f ) / sig2 ).sum();  // p.f is different from p.f used in chi1

				pFTest = fTest( N1, N2, bestNorm, chiSquared( y, p.f ) );

				if ( chi1 - chi2 >= chi2lim2 && pFTest < alpha && !oneEqual( b, brlf) ) {

					if( !alwaysFg ) {
						a_old = b[slice( 0, 5, 1 )];
						b_old = b;
					}

					aResult[slice( resIj, aSize, resDim01 )] = b[slice( 0, aSize, 1 )]; // first 5 elements of b
					bResult[slice( resIj, bSize - aSize, resDim01 )] = b[slice( aSize, bSize - aSize, 1 )]; // last 3 elements of b
					cResult[resIj] =  sqrt( chi1 - chi2 );
					cResult[resIj + resDim01] = pFTest;

					if( info ) {
						cout << "At i=" << i << ", j=" << realInd << " pFTest= " << pFTest << endl;
						cout << "i= " << i << " j= " << realInd << endl;
						printArray( paramTh->index, b, "b" );
					}

				} else {

					if( !alwaysFg )
						a_old = a;

					aResult[slice( resIj, aSize, resDim01 )] = a;
					cResult[resIj] =  sqrt( chi0 - chi1 );
					cResult[resIj + resDim01] = pFTest_old;

					if( info ) {
						cout << "At i=" << i << ", j=" << realInd << " pFTest_old= " << pFTest_old << endl;
						cout << "i= " << i << " j= " << realInd << endl;
						printArray( paramTh->index, a, "a" );
					}
				}

			} else {

				avOfNeighbors( y, true, i, realInd, 3 );
				avOfNeighbors( sig, false, i, realInd, 3 );

				sig2 = stdev( sig );
				sig2 *= sig2;
				av = average( y );
				chi0 = ( ( w * ( y - av ) * ( y - av ) ) / sig2 ).sum();
				p.y = y;


				mpfit( fittingFunctions[fittingFunction], sliceSize, aSize, &a[0], pars, 0, (void *) &p, &result );

				if( result.status < 0 ) {
					cout << "Error while fitting on index i= " << i << " j= " << realInd << " !" << endl;
					cout << "Error status: " << result.status << endl;
				}

				if( info ) {
					cout << "At i=" << i << ", j=" << realInd << " ";
					printArray( paramTh->index, a, "aN" );
				}

				chi1 = ( w * ( y - p.f ) * ( y - p.f ) / sig2 ).sum();
				bestNorm = chiSquared( y, p.f );
				pFTest = fTest( N0, N1, chiSquared( y, sig ), bestNorm );

				if( chi0 - chi1 > chi2lim && pFTest < alpha && !oneEqual( a, alast)  ) {

					pFTest_old = pFTest;

					mpfit( fittingFunctions[fittingFunction + "_br"], sliceSize, bSize, &b[0], pars, 0, (void *) &p, &result );

					if( result.status < 0 ) {
						cout << "Error while fitting on index i= " << i << " j= " << realInd << " !" << endl;
						cout << "Error status: " << result.status << endl;
					}

					if( info ) {
						cout << "At i=" << i << ", j=" << realInd << " ";
						printArray( paramTh->index, b, "bN" );
					}

					double chi2 = ( w * ( y - p.f ) * ( y - p.f ) / sig2 ).sum();  // p.f is different from p.f used in chi1

					pFTest = fTest( N1, N2, bestNorm, chiSquared( y, p.f ) );

					if ( chi1 - chi2 >= chi2lim2 && pFTest < alpha  && !oneEqual( b, brlf ) ) { // p < alpha
						if( !alwaysFg ) {
							a_old = b[slice( 0, 5, 1 )];
							b_old = b;
						}

						aResult[slice( resIj, aSize, resDim01 )] = b[slice( 0, aSize, 1 )]; // first 5 elements of b
						bResult[slice( resIj, bSize - aSize, resDim01 )] = b[slice( aSize, bSize - aSize, 1 )]; // last 3 elements of b
						cResult[resIj] =  sqrt( chi1 - chi2 );
						cResult[resIj + resDim01] = pFTest;

						if( info ) {
							cout << "At i=" << i << ", j=" << realInd << " pFTest= " << pFTest << endl;
							cout << "i= " << i << " j= " << realInd << endl;
							printArray( paramTh->index, b, "b9" );
						}

					} else {

						if( !alwaysFg )
							a_old = a;

						aResult[slice( resIj, aSize, resDim01 )] = a;
						cResult[resIj] =  sqrt( chi0 - chi1 );
						cResult[resIj + resDim01] = pFTest_old;

						if( info ) {
							cout << "At i=" << i << ", j=" << realInd << " pFTest_old= " << pFTest_old << endl;
							cout << "i= " << i << " j= " << realInd << endl;
							printArray( paramTh->index, a, "a9" );
						}
					}
				}  else {
					if( info )
						cout << "Thread #" << paramTh->index << ": No line detected at i= " << i << ", j= " << realInd << endl;
				}
			}

		} // End of first for
	}  // End of second for

	return 0;
}



void readImage(valarray<double> &contents, vector<unsigned int> &dimentions, string fileName) {

	if( info )
	  FITS::setVerboseMode(true);

  try {
	   auto_ptr<FITS> pInfile( new FITS( fileName, Read, true ) );
	   PHDU& image = pInfile->pHDU();

	   // read all user-specifed, coordinate, and checksum keys in the image
	   image.readAllKeys();
	   image.read( contents );

	   int size =  image.axes();
	   dimentions.resize( size );

	   for(int i = 0; i < size; ++i )
		   dimentions[i] = image.axis( i );

	} catch (FitsException&) {

	  cerr << " Fits Exception Thrown by FitsImage class \n";
	  cerr << " Fits file name : " << fileName << endl;
	}
	catch (std::bad_alloc& ba)
	 {
	   cerr << "bad_alloc caught: " << ba.what() << endl;
	   cerr << "Not enough space in RAM to load FITS file" << endl;
	   exit(1);
	 }

}



int writeImage( string fileName, long naxis, vector<long> naxes, valarray<double> va ) {

	if( info )
	  FITS::setVerboseMode(true);

    // declare auto-pointer to FITS at function scope. Ensures no resources
    // leaked if something fails in dynamic allocation.
    auto_ptr<FITS> pFits( 0 );

    try
    {
        // overwrite existing file if the file already exists.
    	fileName = "!" + fileName;

        // Create a new FITS object, specifying the data type and axes for the primary
        // image. Simultaneously create the corresponding file.

        pFits.reset( new FITS( fileName , static_cast<int>(FLOAT_IMG) , naxis, naxes.data() ) );
    }
    catch (FITS::CantCreate*)
    {
          return -1;
    }

    long nelements = 1;

    for (long i = 0; i < naxis; ++i)
    	nelements *= naxes[i];

    // The function PHDU& FITS::pHDU() returns a reference to the object representing
    // the primary HDU; PHDU::write( <args> ) is then used to write the data.

    pFits->pHDU().write( 1, nelements, va );

    return 0;
}


double chiSquared (valarray<double> observed, valarray<double> expected) {
	observed -= expected;
	return ( observed * observed / expected ).sum();
}

double fTest( size_t n0, size_t n1, double bN0, double bN1 ) {

	bN0 -= bN1;
	double F = bN0 / ( ( n0 - n1 ) * ( bN1 / n1 ) );

	if ( F < 0 )
		return 1.0;

	if ( !isnan( F ) && !isinf( F )  ) {
		fisher_f dist( n0 - 1, n1 - 1 );
		return cdf( complement( dist, F ) );
	}
	return 0.0;
}

int oneEqual(valarray<double> a, valarray<double> b) {

	if( compare_guess ) {
		for(size_t i=0; i < a.size(); ++i) {
			if(a[i] == b[i])
				return 1;
		}
	}

	if( compare_limits ) {
		for(size_t i = 0; i < limits.size(); i+=3) {
			if( limits[i] >= a.size() )
				continue;
			if(a[limits[i]] == limits[i+1] || a[limits[i]] == limits[i+2])
				return 1;
		}
	}

	return 0;
}
