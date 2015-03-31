/*
 * fitspec.cpp
 *
 *  Created on: Mar 16, 2015
 *      Author: lerma
 */
#include <iostream>
#include <math.h>
#include <valarray>
#include <omp.h>
#include <time.h>
#include <fstream>
#include <vector>
#include <iterator>
#include <set>

#include <CCfits/CCfits>

#define M_PI           3.14159265358979323846  /* pi */
#define THIRD_DIMENTION 240

using namespace std;
using namespace CCfits;


void readImage(valarray<double> &contents, vector<unsigned int> &dimentions, string fileName)
{
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

     //  cout << "size of the array is: " << contents.size() << '\n';

       // this doesn't print the data, just header info.
       std::cout << image << std::endl;
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
       std::cout << std::endl;
  //     return ;
}


void readData(valarray<double> &content, vector<unsigned int> &dimentions, string filename) {
	ifstream myFile;
	myFile.open( filename.c_str() );

	if( myFile.is_open() ) {
		unsigned int size;
		myFile >> size;
		dimentions.push_back( size );
		cout << "Size of the wl is: " << dimentions[0] << endl;
		content.resize( dimentions[0] );
		double d;
		for(size_t i = 0; i < dimentions[0]; ++i) {
			myFile >> d;
			content[i] = d;
			cout << content[i] << ' ';
		}
		cout << endl;
	} else
		cout << "Can't open file: " << filename << endl;

	myFile.close();
}

void readInputFile(string fileName, vector<double> &alast, vector<double> &brlf, size_t *start, size_t *sigStart) {
	ifstream myFile;
	myFile.open( fileName.c_str() );

	if ( myFile.is_open() ) {
		size_t alastSize, brlfSize;
		myFile >> (*start);
		myFile >> (*sigStart);
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
}


double stdev(valarray<double> va)
{
    double E = 0;
    double ave = va.sum() / va.size();
    double inverse = 1.0 / static_cast<double>(va.size());
    for(size_t i = 0; i < va.size(); ++i) {
        E += pow( static_cast<double>( va[i] ) - ave, 2 );
    }
    return sqrt( inverse * E );
}
/*
double median(valarray<double> va) {
	//multiset<double> ms( va );
	size_t size = ms.size();
	if ( size % 2 )
//		return ms[size / 2];
	return 0.0;
}
*/
/*
void triplet(double x, double *par, double *f, double *pder) {
	double z = par[0],
		  inten = par[1],
		  ww = par[2],
		  contin = par[3],
		  NIIinten = par[4],
		  Ha = 6562.8 * ( 1 + z ),
		  NII_r = 6583 * ( 1 + z ),
		  NII_b = 6548.1 * ( 1 + z );

	ww = sqrt( 0.6 ^ 2 + ww ^ 2 );//         ; default width to match sky line.

	Ha = ( inten / ww / sqrt( 2.0 * M_PI ) ) * ( exp( -0.5 * ( x - Ha ) ^ 2 / ww ^ 2 ) );
	NII_r = ( NIIinten / ww / sqrt( 2.0 * M_PI ) ) * ( exp( -0.5 * ( x - NII_r ) ^ 2 / ww ^ 2 ) );
	NII_b = ( NIIinten / 3. / ww / sqrt( 2.0 * M_PI ) ) * ( exp( -0.5 * ( x - NII_b ) ^ 2 / ww ^ 2 ) );

	f = contin + Ha + NII_r + NII_b;

	//pder = fltarr(n_elements(x),n_elements(par));//   ; no value returned.
}


void tripletbr(double x, double *par, double *f, double *pder) {

	double z = par[0],
		  inten = par[1],
		  ww = par[2],
		  contin = par[3],
		  NIIinten = par[4],
		  brint = par[5],
		  brww = par[6],
		  zbr = par[7],
		  Ha = 6562.8 * ( 1 + z ),
		  NII_r = 6583 * ( 1 + z ),
		  NII_b = 6548.1 * ( 1 + z ),
		  br = 6562.8 * ( 1 + zbr );

	ww = sqrt( 0.6 ^ 2 + ww ^ 2 ); //        ; default width to match sky line.
	brww = sqrt( 0.6 ^ 2 + brww ^ 2 );


	//;if ww ge 10 then ww = 10

	Ha = ( inten / ww / sqrt( 2.0 * M_PI ) ) * ( exp( -0.5 * ( x - Ha ) ^ 2 / ww ^ 2 ) );
	NII_r = ( NIIinten / ww / sqrt( 2.0 * M_PI ) ) * ( exp( -0.5 * ( x - NII_r ) ^ 2 / ww ^ 2 ) );
	NII_b = ( NIIinten / 3. / ww / sqrt( 2.0 * M_PI ) ) * ( exp( -0.5 * ( x - NII_b ) ^ 2 / ww ^ 2 ) );
	br = ( brint / brww / sqrt( 2.0 * M_PI ) ) * ( exp( -0.5 * ( x - br ) ^ 2 / brww ^ 2 ) );

	f = contin + Ha + NII_r + NII_b + br;

	//pder = fltarr( n_elements( x ), n_elements( par ) );
}
*/


int main(int argc, char* argv[]) {

	//size_t startIndex = 0, sigStartIndex = 0;
	string inputFile = "inputfile.txt";

	if ( argc == 1 ) {
		cout << "You should provide input file name. \nIf not, the default will be used: inputfile.txt" << '\n';
	} else if ( argc == 2 ) {
			inputFile = argv[1];
		} else
			cout << "Only one parameter should be passed. \nThe rest will be ignored" << '\n';


/*
	sig2d = finalcube[*,*,2200:2440];
	finalcube = finalcube[*,*,2605:2845];
	wl = wl[2605:2845];
	sz = size(finalcube);
	sky = sky[2605:2845];
	var = 1/sky;
*/
	clock_t ts = clock();
	valarray<double> finalCube, sky, wl, var;
	vector<unsigned int> finalCubeDimentions, skyDimentions, wlDimentions;
	vector<double> alast, brlf;
	size_t startIndex, sigStartIndex;

	readInputFile( "inputfile.txt", alast, brlf, &startIndex, &sigStartIndex );
	readImage( finalCube, finalCubeDimentions, "coadd-HaNII.fits" );
	readImage( sky, skyDimentions, "sky_HaNII.fits" );
	readImage( wl, wlDimentions, "wl_HaNII.fits" );





	sky.resize( THIRD_DIMENTION );
	var.resize( THIRD_DIMENTION );
	var = 1. / sky;

	valarray<double> fitcube( finalCubeDimentions[0] * finalCubeDimentions[2] * THIRD_DIMENTION, 0.0 ),
					 acube( finalCubeDimentions[0] * finalCubeDimentions[2] * 5, 0.0 ),
					 bcube( finalCubeDimentions[0] * finalCubeDimentions[2] * 3, 0.0 ),
					 aerrorcube( finalCubeDimentions[0] * finalCubeDimentions[2] * 8, 0.0 );

	valarray<double> sig2d( finalCubeDimentions[0] * finalCubeDimentions[1] * THIRD_DIMENTION );
	//double alast[] = { 0.034354, 0.1, 2.0, 0, 0.1, 5 };
	//double brlf[] = { 0.034354, 0.5, 2.0, 0, 0.3, 1, 10, 0.034354 };


	valarray<double> akeep( alast[0], alast.size() ), bkeep( brlf[0], brlf.size() );

	for(size_t outterIndex = 3; outterIndex < finalCubeDimentions[1] - 3; ++outterIndex)
		for(size_t innerIndex = 3; innerIndex < finalCubeDimentions[0] - 3; ++innerIndex) {


		}


	for(size_t i = 0; i < alast.size(); ++i)
		cout << "alast[" << i << "]= " << alast[i] << endl;
	for(size_t i = 0; i < brlf.size(); ++i)
			cout << "brlf[" << i << "]= " << brlf[i] << endl;

	cout << "The dimension of final cube is " << finalCubeDimentions.size() << endl;

	for(size_t i = 0; i < finalCubeDimentions.size(); ++i)
		cout << "The size of dimension " << i << " is " << finalCubeDimentions[i] << endl;
	//readTable();
	cout << "The size of final cube is " << finalCube.size() << endl;


	cout << "The size of sig2d is " << sig2d.size() << endl;
//	omp_set_num_threads( 4 );
//#pragma omp parallel for schedule(static) collapse(3)
	for(size_t i = 0; i < finalCubeDimentions[0]; ++i)
		for(size_t j = 0; j < finalCubeDimentions[1]; ++j)
			for(size_t k = 0; k < THIRD_DIMENTION; ++k) {
				size_t index = i * finalCubeDimentions[1] * THIRD_DIMENTION + j * THIRD_DIMENTION + k;
				sig2d[index] = finalCube[index];
				//cout << sig2d[index] << " ";
			}
	cout << "The size of sig2d is " << sig2d.size() << endl;

	//valarray<double> akeep ( alast ), bkeep( brlf );

	for(size_t j = 3; j < finalCubeDimentions[1] - 3; ++j)
		for(size_t i = 3; i < finalCubeDimentions[0]; ++i) {
			valarray<double> a( alast[0], alast.size() ), b( brlf[0], brlf.size() );
		//	size_t i = 35;
		//	j = 45;
			valarray<double> x( wl ), y( 0.0, THIRD_DIMENTION ), sig( 0.0, THIRD_DIMENTION );
			for(size_t k = 0; finalCube[2] - 1; ++k) {
				size_t index = i * finalCubeDimentions[1] * THIRD_DIMENTION + j * THIRD_DIMENTION + k;
				y[k] = finalCube[index];
				sig[k] = sig2d[index];
			}

			double sig2 = stdev( sig );
			sig2 *= sig2;

			valarray<double> w( 1.0, THIRD_DIMENTION );
		//	w = var / median( var );

	}

	cout << "Execution time is " << (double)(clock() - ts)/CLOCKS_PER_SEC << endl;
	return 0;
}



