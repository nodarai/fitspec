/*
 * fitspec.cpp
 *
 *  Created on: Mar 16, 2015
 *      Author: lerma
 */
#include <iostream>

#include <CCfits/CCfits>

//using namespace std;
using namespace CCfits;


int readImage()
{
       std::auto_ptr<FITS> pInfile(new FITS("coadd-C2.fits",Read,true));

       PHDU& image = pInfile->pHDU();

       std::valarray<unsigned long>  contents;

       // read all user-specifed, coordinate, and checksum keys in the image
       image.readAllKeys();

       image.read(contents);

       // this doesn't print the data, just header info.
       std::cout << image << std::endl;

       long ax1(image.axis(0));
       long ax2(image.axis(1));


       for (long j = 0; j < ax2; j+=10)
       {
               std::ostream_iterator<short> c(std::cout,"\t");
               std::copy(&contents[j*ax1],&contents[(j+1)*ax1-1],c);
               std::cout << '\n';
       }

       std::cout << std::endl;
       return 0;
}


int readTable()
{
       // read a table and explicitly read selected columns. To read instead all the
       // data on construction, set the last argument of the FITS constructor
       // call to 'true'. This functionality was tested in the last release.
       std::vector<string> hdus(2);
       hdus[0] = "PLANETS_ASCII";
   //    hdus[1] = "TABLE_BINARY";


       std::auto_ptr<FITS> pInfile(new FITS("example.fits",Read,hdus,false));

       ExtHDU& table = pInfile->extension(hdus[1]);



       std::vector < valarray <int > > pp;
       table.column("powerSeq").readArrays( pp, 1,3 );

       std::vector < valarray <std::complex<double> > > cc;
       table.column("dcomplex-roots").readArrays( cc, 1,3 );

       std::valarray < std::complex<float> > ff;
       table.column("fcomplex-roots").read( ff, 4 );

       std::cout << pInfile->extension(hdus[0]) << std::endl;

       std::cout << pInfile->extension(hdus[1]) << std::endl;

       return 0;
}



int main() {

	readImage();

	//readTable();

	return 0;
}



