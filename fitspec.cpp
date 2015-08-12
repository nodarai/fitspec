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
#include <cstdio>
#include <limits>
#include <assert.h>

#include <mpi.h>

#include <CCfits/CCfits>

//#include <boost/math/distributions/fisher_f.hpp>

//#include "fitspec.h"
#include "functions.h"
#include "mpfit.h"
#include "cJSON/cJSON.h"
//#define THIRD_DIMENTION 240

using namespace std;
using namespace CCfits;
//using namespace boost::math;



// Global variables********************

size_t sliceSize,
	   aSize,
	   bSize;

double chi2lim,
	   chi2lim2;

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




void* runThreads   (void*);
void readImage	   (valarray<double>&, vector<unsigned int>&, string);
int  writeImage	   (string, long, vector<long>, valarray<double>);
void printArray	   (size_t, valarray<double>, string);
void printMpiError (int);



/*

int writeToAscii(size_t i, size_t j, valarray<double> va, ofstream& myFile) {

	myFile << i << " " << j << " ";
	for(size_t ind = 0; ind < va.size(); ++ind) {
		myFile << va[ind] << " ";
	}
	myFile << "\n";

	return 0;
}
*/






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
//		cout << "sliceSize is: " << (*sliceSize) << endl;
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

	for( i = 0/*cJSON* iter = inputParams->child->next*/; iter; iter = iter->next, ++i ) {

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
					valarray<double> vaTmp( v.data()[0], v.size() );
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
					valarray<double> vaTmp ( v.data()[0], v.size() );
					brlf.resize( vaTmp.size() );
					brlf = vaTmp;
				}
		}
	}

	in.close();

	return 0;
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
 /*
    for(size_t i = 0; i < va.size(); ++i) {
        E += pow( va[i] - ave, 2 );
    }
*/
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
    if ( ( p == ( A.size() / 2 ) ) && ( ( A.size() % 2 ) == 1 ) ) {
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
 *
void avOfNeighbors (valarray<double> &y,/* valarray<double> va, * size_t i, size_t j, size_t neibSize) {

	double sum = 0.0;
	unsigned int dim12 = finalCubeDimentions[2] * finalCubeDimentions[1];
	unsigned int dim2 = finalCubeDimentions[2];
//	unsigned int jdim2 = finalCubeDimentions[2] * j;
//	size_t k = 0;
/*
	sum += va[i 		* dim12 + dim2 * j 		   + k]
		 + va[i 		* dim12 + dim2 * ( j - 1 ) + k]
		 + va[i 		* dim12 + dim2 * ( j + 1 ) + k]
		 + va[( i - 1 ) * dim12 + dim2 * j 		   + k]
		 + va[( i - 1 ) * dim12 + dim2 * ( j - 1 ) + k]
		 + va[( i - 1 ) * dim12 + dim2 * ( j + 1 ) + k]
		 + va[( i + 1 ) * dim12 + dim2 * j 		   + k]
		 + va[( i + 1 ) * dim12 + dim2 * ( j - 1 ) + k]
		 + va[( i + 1 ) * dim12 + dim2 * ( j + 1 ) + k];

	cout << sum / 9.0 << " ";
	sum = 0.0;

	for(k = 1; k < sliceSize; ++k) {
		for(size_t l = i - neibSize / 2; l <= i + neibSize / 2; ++l)
			for(size_t m = j - neibSize / 2; m <= j + neibSize / 2; ++m)
				sum+= va[l * dim12 + m * dim2 + k];
		cout << sum / 9.0 << " ";
		sum = 0.0;
	}
	cout << endl;
*
	//size_t ar[9];// =
	vector<size_t> v ( neibSize * neibSize );
	for(size_t n = 0, l = i - neibSize / 2; l <= i + neibSize / 2; ++l)
		for(size_t m = j - neibSize / 2; m <= j + neibSize / 2; ++m)
			v[n++] = l * dim12 + m * dim2;

	for(size_t k = 0; k < sliceSize; ++k, sum = 0.0) {
		for(size_t l = 0; l < v.size(); ++l)
			sum += finalCube[v[l] + k];
		y[k] =  sum / (double) v.size();
//		cout << y[k] << " ";
//		sum = 0.0;
	}
//	cout << endl;
/*

	for(size_t l = i - neibSize / 2; l <= i + neibSize / 2; ++l)
		for(size_t m = j - neibSize / 2; m <= j + neibSize / 2; ++m)
			sum+= va[l * dim12 + m * dim2 + k];

	y[k] = sum/ 9.0;


	for(k = 1; k < sliceSize; ++k) {
		for(size_t l = i - neibSize / 2, m = j - neibSize / 2 - 1; l < i + neibSize / 2; ++l)
			sum -= va[l * dim12 + j * dim2 + k];
		for(size_t l = i - neibSize / 2, m = j + neibSize / 2; l < i + neibSize / 2; ++l)
			sum -= va[l * dim12 + j * dim2 + k];
	}
*

//	return sum / 9.0;
}
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



int main(int argc, char* argv[]) {


/*

	// standard normal distribution object:
	normal norm;
	// print survival function for x=2.0:
	cout << cdf(norm, 2.0) << endl;

	return 0;
*/



	time_t ts = time( NULL );
	time( &ts );

	cout.rdbuf()->pubsetbuf( 0, 0 ); //Turn off buffered output

//	cout << numeric_limits<double>::max() << endl;
//	return 0;


	//size_t startIndex = 0, sigStartIndex = 0;
	string inputFile = "inputfile.json",
		   resultFile;
	/*This vector contains filenames of FITS input data
	  0: finalcube, 1: sky, 2: wl */
	vector<string> inputFits( 3 );

	if ( argc == 1 ) {
		cout << "You should provide input file name. \nIf not, the default will be used: inputfile.json" << '\n';
	} else if ( argc == 2 ) {
			inputFile = argv[1];
		   } else
			   cout << "Only one parameter should be passed. \nThe rest will be ignored" << '\n';

	size_t startIndex,
		   sigStartIndex,
		   nbThreads;

	//readInputFile( "inputfile.txt", alast, brlf, &startIndex, &sigStartIndex, &sliceSize, &nbThreads );
	readInputFileJSON( inputFile, resultFile, inputFits, &startIndex, &sigStartIndex, &nbThreads );

	readImage( finalCube, finalCubeDimentions, inputFits[0] );
	readImage( sky, skyDimentions, inputFits[1] );
	readImage( wl, wlDimentions, inputFits[2] );


	cout << "resultfile= " << resultFile << endl;

	cout << "input fits= " << inputFits[0] << endl;

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

/*
	finalCubeDimentions[0] = finalCubeDimentions[1] = finalCubeDimentions[2] = 5;
	finalCube.resize( finalCubeDimentions[0] * finalCubeDimentions[1] * finalCubeDimentions[2] );

	for (size_t i = 0; i < finalCube.size(); ++i)
		finalCube[i] = i + 1;
*/
	//sliceSize = 3;



/*
	valarray<double> tmpFinalCube( finalCubeDimentions[0] * finalCubeDimentions[1] * sliceSize );

	for (size_t k = startIndex; k < startIndex + sliceSize; ++k)
		for (size_t j = 0; j < finalCubeDimentions[1]; ++j) {
			size_t ij =
			for(size_t i = 0; i < finalCubeDimentions[0]; ++i)
				tmpFinalCube[k - startIndex] = finalCube[];
		}



	size_t first = startIndex,
		   last = sigStartIndex + sliceSize;
	if( sigStartIndex < startIndex ) {
		first = sigStartIndex;
		last = startIndex + sliceSize;
	}
*/




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

/*
	cout << endl;

	for(size_t i = 0, l = 0, m = 0; i < finalCubeDimentions[0]; ++i)
		for(size_t j = 0; j < finalCubeDimentions[1]; ++j) {

			size_t  ij = i * finalCubeDimentions[1] + j;

			for(size_t k = first; k < last; ++k) {

				size_t ijk = k * finalCubeDimentions[1] * finalCubeDimentions[2] + ij;

				cout << ijk << " ";

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

	cout << endl;

	*/

/*
	vector<long> aResDims2( 3 );
	aResDims2[0] = finalCubeDimentions[0];
	aResDims2[1] = finalCubeDimentions[1];
	aResDims2[2] = finalCubeDimentions[2];



	writeImage( "FCreadcut" + resultFile, 3, aResDims2, finalCube );

	return 0;
*/


	/*
	cout << "The size of sig2d is " << sig2d.size() << endl;
	cout << "The size of tmpFinalCube is " << tmpFinalCube.size() << endl;
	cout << "The size of finalCube is " << finalCube.size() << endl;
*/


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


	vector<long> aResDims( 3 );
	aResDims[0] = finalCubeDimentions[0];// - 6;
	aResDims[1] = finalCubeDimentions[1];// - 6;
	aResDims[2] = aSize;

	vector<long> bResDims( 3 );
	bResDims[0] = finalCubeDimentions[0];// - 6;
	bResDims[1] = finalCubeDimentions[1];// - 6;
	bResDims[2] = bSize - aSize;

	aResult.resize( aResDims[0] * aResDims[1] * aResDims[2] );
	bResult.resize( bResDims[0] * bResDims[1] * bResDims[2] );
	//aResult.resize( finalCube.size() );
	//bResult.resize( finalCube.size() );





/*
	vector<long> aResDims2( 3 );
	aResDims2[0] = finalCubeDimentions[0];
	aResDims2[1] = finalCubeDimentions[1];
	aResDims2[2] = finalCubeDimentions[2];



	writeImage( "FCread" + resultFile, 3, aResDims2, finalCube );

	return 0;
*/



	int rank,
	    worldSize,
		provided;

	MPI_Init_thread ( &argc, &argv, MPI_THREAD_MULTIPLE, &provided );      /* starts MPI */
	if ( provided != MPI_THREAD_MULTIPLE )
	{
		printf( "Sorry, this MPI implementation does not support multiple threads\n" );
		MPI_Abort( MPI_COMM_WORLD, 1 );
		return 1;
	}

	MPI_Comm_rank ( MPI_COMM_WORLD, &rank );        /* get current process id */
	MPI_Comm_size ( MPI_COMM_WORLD, &worldSize );        /* get number of processes */

	vector<pthread_t> thread( nbThreads );
	vector<ParamsThread> paramsThreads( nbThreads ); //= new ParamsThread[nbThreads];

	size_t jobSize = ( finalCubeDimentions[1] - 6 ) / ( nbThreads * worldSize);
	size_t remainder = ( finalCubeDimentions[1] - 6 ) % (nbThreads * worldSize );
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

	//paramsThreads[0].jGlobal = 3; //We are starting calculations from index 3

	//paramsThreads[nbThreads - 1].index = nbThreads - 1;
	//paramsThreads[nbThreads - 1].jobSize = finalCubeDimentions[0] / nbThreads + finalCubeDimentions[0] % nbThreads;
	//paramsThreads[nbThreads - 1].iGlobal = ( nbThreads - 1 ) * ( finalCubeDimentions[0] / nbThreads );



	for (size_t i = 0; i < nbThreads; ++i) {
		pthread_create( &thread[i], NULL, runThreads, (void*)&paramsThreads[i] );
	}

	for (size_t i = 0; i < nbThreads; ++i) {
		pthread_join( thread[i], NULL );
	}
/*
	if ( rank == 0 ) {
		//receive
		int error;
		valarray<double> aOthers( aResult.size() );
		valarray<double> bOthers( bResult.size() );

		for(size_t i = 1; i < worldSize; ++i) {
			error = MPI_Recv( &aOthers[0], aResult.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
			printMpiError( error );
			aResult += aOthers;

			error = MPI_Recv( &bOthers[0], bResult.size(), MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
			printMpiError( error );
			bResult += bOthers;
		}

	} else {
		//send
		int error;
		error = MPI_Send( &aResult, aResult.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		printMpiError( error );

		error = MPI_Send( &bResult, bResult.size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
		printMpiError( error );

	}
*/
	MPI_Finalize();

//	omp_set_num_threads( 2 );
//#pragma omp parallel
//	{




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




	writeImage( "a" + resultFile, 3, aResDims, aResult );
	writeImage( "b" + resultFile, 3, bResDims, bResult );


	//Here ends the process-parallel routines


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

/*
	ofstream aMyFile, bMyFile;
	aMyFile.open( "aASCIIRes.txt" );
	bMyFile.open( "bASCIIRes.txt" );
*/

//	cout << "jobsize= " << paramTh->jobSize << endl;
//	cout << "glovbal I" << paramTh->iGlobal << endl;

	for(size_t j = 0; ( j < paramTh->jobSize ) && ( j + paramTh->jGlobal < finalCubeDimentions[1] - 3 ); ++j) {

		size_t realInd = j + paramTh->jGlobal; // j index of the finalcube
		a_old = alast;
		b_old = brlf;


		//cout << "Global J= " << paramTh->jGlobal << endl;


		for(size_t i = 3; i < finalCubeDimentions[0] - 3; ++i) {



			size_t ij = realInd * finalCubeDimentions[0] + i;
			size_t resIj = realInd * ( finalCubeDimentions[0] ) + i;


		//	size_t ij = realInd * finalCubeDimentions[1] + j;
		//	size_t resIj = realInd * ( finalCubeDimentions[1] - 6 ) + j;

		//	size_t  ij = j * finalCubeDimentions[0] + realInd;
		//	size_t resIj = j * ( finalCubeDimentions[0] - 6 ) + realInd;

	//			cout << "Thread# " << omp_get_thread_num() << " entered in for" << endl;

	//			cout << "nesed is active? " << omp_get_nested() << " nested level= " << omp_get_active_level() << endl;

			//valarray<double> a( alast.data(), alast.size() ), b( brlf.data(), brlf.size() );
			//double *a = alast.data(), *b = brlf.data();


			//doule a[5] = { 0.0543693, 10.1, 2.0, 0.0 10.1 }:
		//	size_t i = 35;
		//	j = 45;


			//old version of indexing
			//	size_t  ij =  i * finalCubeDimentions[2] * finalCubeDimentions[1] + finalCubeDimentions[2] * j;



			/*
			size_t ijk = k * finalCubeDimentions[0] * finalCubeDimentions[1] + ij;
			*/

			for(size_t k = 0; k < sliceSize; ++k) {

				//old version of indexing
				//size_t ijk = k + ij;
				size_t ijk = k * fcDim01 + ij;

				//cout << "finalCube[" << i << "][" << j << "][" << k << "]= " << finalCube[ijk] << endl;
				//cout << "finalCube[" << k << "]= " << finalCube[k] << endl;
				y[k] = finalCube[ijk];
				sig[k] = sig2d[ijk];
			}
			//exit(0);

			double sig2 = stdev( sig );
			sig2 *= sig2;

			double av = average( y );

	/*			valarray<double> justFor( sliceSize );
			justFor = y - av;
			for(size_t ii = 0; ii < y.size(); ++ii)
				cout << "y[" << ii << "]= " << y[ii] << endl;
	*/
			double chi0 = ( ( w * ( y - av ) * ( y - av ) ) / sig2 ).sum();
			//printArray(y.size(), y, "y");
			cout << "At i=" << i << ", j=" << realInd << " chi0= " << chi0 << endl;

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



	//		printArray( 6, a, "old_a" )

			a = a_old;
			b = b_old;

			mpfit( triplet, sliceSize, aSize, &a[0], 0, 0, (void *) &p, &result );

			if( result.status < 0 ) {
				cout << "Error while fitting on index i= " << realInd << " j= " << j << " !" << endl;
				cout << "Error status: " << result.status << endl;
		//		continue;
			}
			cout << "At i=" << i << ", j=" << realInd << " ";
			printArray( paramTh->index, a, "a" );



			double chi1 = ( w * ( y - p.f ) * ( y - p.f ) / sig2 ).sum();
			//double chi1r = chi1 / ( sliceSize - aSize + 1 );  // y.size() is replaced  with sliceSize as tey'r equal

			cout << "At i=" << i << ", j=" << realInd << " chi1= " << chi1 << endl;
/*			if ( ( realInd == 166 && j == 69 ) ||
				 ( realInd == 168 && j == 61 ) ||
				 ( realInd == 158 && j == 214 ) ||
				 ( realInd == 143 && j == 221 ) ||
				 ( realInd == 110 && j == 235 )
				 ) {
				cout << "At i=" << realInd << ", j=" << " chi1= " << chi1 << endl;
			}
*/

			if( chi0 - chi1 > chi2lim ) {

				mpfit( tripletbr, sliceSize, bSize, &b[0], 0, 0, (void *) &p, &result );

				cout << "At i=" << i << ", j=" << realInd << " ";
				printArray( paramTh->index, b, "b" );

				double chi2 = ( w * ( y - p.f ) * ( y - p.f ) / sig2 ).sum();  // p.f is different from p.f used in chi1
				//double chi2r = chi1 / ( sliceSize - bSize + 1 );

//				cout << "At i=" << i << ", j=" << realInd << " chi2= " << chi2 << endl;
/*
				if ( ( realInd == 166 && j == 69 ) ||
					 ( realInd == 168 && j == 61 ) ||
					 ( realInd == 158 && j == 214 ) ||
					 ( realInd == 143 && j == 221 ) ||
					 ( realInd == 110 && j == 235 )
					 ) {
					cout << "At i=" << realInd << ", j=" << " chi2= " << chi2 << endl;
				}
*/

				if ( chi1 - chi2 >= chi2lim2 ) {
			//		a_old = b[slice( 0, 5, 1 )];
			//		b_old = b;

			//		writeToAscii( realInd, j, b[slice( 0, 5, 1 )], aMyFile );
			//		writeToAscii( realInd, j, b[slice( 5, 3, 1 )], bMyFile );
					aResult[slice( resIj, aSize, resDim01 )] = b[slice( 0, aSize, 1 )]; // first 5 elements of b
					bResult[slice( resIj, bSize - aSize, resDim01 )] = b[slice( aSize, bSize - aSize, 1 )]; // last 3 elements of b
					cout << "i= " << i << " j= " << realInd << endl;
					printArray( paramTh->index, b, "b" );
	//				++bCount;
				//	exit(0);
				} else {
			//		a_old = a;
			//		writeToAscii( realInd, j, a, aMyFile );
					aResult[slice( resIj, aSize, resDim01 )] = a;
					cout << "i= " << i << " j= " << realInd << endl;
			//		cout << "i= " << paramTh->iGlobal + i << " j= " << j << endl;
					printArray( paramTh->index, a, "a" );
				}

	//				++aCount;
			} else {

	//				clock_t neighb = clock();
	//			cout << "i= " << paramTh->iGlobal + i << " j= " << j << endl;
				avOfNeighbors( y, true,/*finalCube,*/ i, realInd, 3 );
				avOfNeighbors( sig, false,/*sig2d,*/ i, realInd, 3 );
	//				cout << "Neighbor count time is " << (double)(clock() - neighb)/CLOCKS_PER_SEC << endl;

				sig2 = stdev( sig );
				sig2 *= sig2;
				av = average( y );
				chi0 = ( ( w * ( y - av ) * ( y - av ) ) / sig2 ).sum();
				//printArray(y.size(), y, "y after");

				cout << "At i=" << i << ", j=" << realInd << " after binfing chi0= " << chi0 << endl;

				//chi0r = chi0 / y.size();

				p.y = y;


				mpfit( triplet, sliceSize, aSize, &a[0], 0, 0, (void *) &p, &result );

				if( result.status < 0 ) {
					cout << "Error while fitting on index i= " << i << " j= " << realInd << " !" << endl;
					cout << "Error status: " << result.status << endl;
				//	continue;
				}

				cout << "At i=" << i << ", j=" << realInd << " ";
				printArray( paramTh->index, a, "aN" );

				chi1 = ( w * ( y - p.f ) * ( y - p.f ) / sig2 ).sum();
				cout << "At i=" << i << ", j=" << realInd << " after binding chi1= " << chi1 << endl;

				//chi1r = chi1 / ( sliceSize - aSize + 1 );

				if( chi0 - chi1 > chi2lim ) {

					mpfit( tripletbr, sliceSize, bSize, &b[0], 0, 0, (void *) &p, &result );

					cout << "At i=" << i << ", j=" << realInd << " ";
					printArray( paramTh->index, b, "bN" );

					double chi2 = ( w * ( y - p.f ) * ( y - p.f ) / sig2 ).sum();  // p.f is different from p.f used in chi1
					//double chi2r = chi1 / ( sliceSize - bSize + 1 );

					if ( chi1 - chi2 >= chi2lim2 ) {
			//			a_old = b[slice( 0, 5, 1 )];
			//			b_old = b;
		//				writeToAscii( realInd, j, b[slice( 0, 5, 1 )], aMyFile );
		//				writeToAscii( realInd, j, b[slice( 5, 3, 1 )], bMyFile );
						aResult[slice( resIj, aSize, resDim01 )] = b[slice( 0, aSize, 1 )]; // first 5 elements of b
						bResult[slice( resIj, bSize - aSize, resDim01 )] = b[slice( aSize, bSize - aSize, 1 )]; // last 3 elements of b
						cout << "i= " << i << " j= " << realInd << endl;
						printArray( paramTh->index, b, "b9" );

						//						++bCount;
					} else {
				//		a_old = a;
				//		writeToAscii( realInd, j, a, aMyFile );
						aResult[slice( resIj, aSize, resDim01 )] = a;
						cout << "i= " << i << " j= " << realInd << endl;
						printArray( paramTh->index, a, "a9" );
					}

	//					++aCount;
				}  else
					cout << "Thread #" << paramTh->index << ": No line detected at i= " << i << ", j= " << realInd << endl;
			}
		} // End of first for
	}  // End of second for



//	aMyFile.close();
//	bMyFile.close();

	return 0;
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

	//	   cout << "size of the array is: " << contents.size() << '\n';

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



int writeImage( string fileName, long naxis, vector<long> naxes, valarray<double> va ) {

    // declare auto-pointer to FITS at function scope. Ensures no resources
    // leaked if something fails in dynamic allocation.
    auto_ptr<FITS> pFits( 0 );

    try
    {
        // overwrite existing file if the file already exists.

    	fileName = "!" + fileName;
//        const std::string fileName("!atestfil.fit");

        // Create a new FITS object, specifying the data type and axes for the primary
        // image. Simultaneously create the corresponding file.

        // this image is unsigned short data, demonstrating the cfitsio extension
        // to the FITS standard.

        pFits.reset( new FITS( fileName , static_cast<int>(FLOAT_IMG) , naxis, naxes.data() ) );
    }
    catch (FITS::CantCreate*)
    {
          // ... or not, as the case may be.
          return -1;
    }

    long nelements = 1;

    for (long i = 0; i < naxis; ++i)
    	nelements *= naxes[i];

    // The function PHDU& FITS::pHDU() returns a reference to the object representing
    // the primary HDU; PHDU::write( <args> ) is then used to write the data.

    pFits->pHDU().write( 1, nelements, va );


    // PHDU's friend ostream operator. Doesn't print the entire array, just the
    // required & user keywords, and is provided largely for testing purposes [see
    // readImage() for an example of how to output the image array to a stream].

  //  std::cout << pFits->pHDU() << std::endl;

    return 0;
}

void printMpiError( int error ) {
    switch( error ) {
    case MPI_SUCCESS:
  	  //cout << "MPI succeed to start procedure\n";
  	  break;
    case MPI_ERR_COMM:
  	  cout << "MPI failed at invalid communicator" << endl; break;
    case MPI_ERR_COUNT:
  	  cout << "MPI failed at invalid count" << endl; break;
    case MPI_ERR_TYPE:
  	  cout << "MPI failed at invalid data type argument" << endl; break;
    case MPI_ERR_TAG:
  	  cout << "MPI failed at invalid tag" << endl; break;
    case MPI_ERR_RANK:
  	  cout << "MPI failed at invalid rank" << endl; break;
    case MPI_ERR_INTERN:
  	  cout << "MPI failed at invalid MPICH implementation" << endl; break;
    }
}





