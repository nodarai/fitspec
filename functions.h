#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "fitspec.h"


const double piSqr = sqrt( 2.0 * M_PI );




int tripletHb(int m, int n, double *p, double *deviates, double **derivs, void *data) {
	//std::cout << "in function tripletHb" << std::endl;
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
		  ha = 4861.32 * ( 1 + z ),
		  NII_r = 5007.00 * ( 1 + z ),
		  dNII_b = 4959.00 * ( 1 + z );
		 

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

	return 0;
}


int tripletHb_br(int m, int n, double *p, double *deviates, double **derivs, void *data) {
	//std::cout << "in function tripletHb_br" << std::endl;
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
		  ha = 4861.32 * ( 1 + z ),
		  dNII_r = 5007.00 * ( 1 + z ),
		  dNII_b = 4959.00 * ( 1 + z ),
		  dbr = 4861.32 * ( 1 + zbr );

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



int Hb(int m, int n, double *p, double *deviates, double **derivs, void *data) {
	//std::cout << "in function Hb" << std::endl;
	struct Parameters *pr = (struct Parameters *) data;
	valarray<double> x( pr->x ),
					 y( pr->y ),
					 yErr( pr->yErr );

	double z = p[0],
		   ww = p[1],
		   Hb = 4861.32 * ( 1 + z ),
		   OIII_b = 4959.00 * ( 1 + z ),
		   OIII_r = 5007.00 * ( 1 + z ),
		   Hbinten = p[2],
		   OIIIinten = p[3],
		   cont = p[4];

	valarray<double> Hb_ln = ( Hbinten / ww / piSqr ) * ( exp( -0.5 * ( x - Hb ) * ( x - Hb ) / ( ww * ww ) ) );
	valarray<double> OIII_a_ln = ( OIIIinten / ww / piSqr ) * ( exp( -0.5 * ( x - OIII_r ) * ( x - OIII_r ) / ( ww * ww ) ) );
	valarray<double> OIII_b_ln = ( OIIIinten /3. / ww / piSqr ) * ( exp( -0.5 * ( x - OIII_b ) * ( x - OIII_b ) / ( ww * ww ) ) );

	pr->f = cont + Hb_ln + OIII_a_ln + OIII_b_ln;
	pr->f *= 10000;

	for( size_t i = 0; i < y.size(); ++i ) {
		deviates[i] = ( y[i] - pr->f[i] ) / yErr[i];
	}

	return 0;
}


int Hb_br(int m, int n, double *p, double *deviates, double **derivs, void *data) {
//	std::cout << "in function tripletHb_br" << std::endl;
	struct Parameters *pr = (struct Parameters *) data;
	valarray<double> x( pr->x ),
					 y( pr->y ),
					 yErr( pr->yErr );

	double z = p[0],
		   ww = p[1],
		   Hb = 4861.32 * ( 1 + z ),
		   OIII_b = 4959.00 * ( 1 + z ),
		   OIII_r = 5007.00 * ( 1 + z ),
		   Hbinten = p[2],
		   OIIIinten = p[3],
		   cont = p[4],
		   z_br = p[5],
		   ww_br = p[6],
		   Hb_br = 4861.32 * ( 1 + z_br ),
		   OIII_b_br = 4959.00 * ( 1 + z_br ),
		   OIII_r_br = 5007.00 * ( 1 + z_br ),
		   Hbinten_br = p[7],
		   OIIIinten_br = p[8];

	valarray<double> Hb_ln = ( Hbinten / ww / piSqr ) * ( exp( -0.5 * ( x - Hb ) * ( x - Hb ) / ( ww * ww ) ) );
	valarray<double> OIII_a_ln = ( OIIIinten / ww / piSqr ) * ( exp( -0.5 * ( x - OIII_r ) * ( x - OIII_r ) / ( ww * ww ) ) );
	valarray<double> OIII_b_ln = ( OIIIinten /3. / ww / piSqr ) * ( exp( -0.5 * ( x - OIII_b ) * ( x - OIII_b ) / ( ww * ww ) ) );

	valarray<double> Hb_ln_br = ( Hbinten / ww_br / piSqr ) * ( exp( -0.5 * ( x - Hb_br ) * ( x - Hb_br ) / ( ww_br * ww_br ) ) );
	valarray<double> OIII_a_ln_br = ( OIIIinten /ww_br / piSqr ) * ( exp( -0.5 * ( x - OIII_r_br ) * ( x - OIII_r_br ) / ( ww_br * ww_br ) ) );
	valarray<double> OIII_b_ln_br = ( OIIIinten / 3. / ww_br / piSqr ) * ( exp( -0.5 * ( x - OIII_b_br ) * ( x - OIII_b_br ) / ( ww_br * ww_br ) ) );

	pr->f  = cont + Hb_ln + OIII_a_ln + OIII_b_ln + Hb_ln_br + OIII_a_ln_br + OIII_b_ln_br;
	pr->f *= 10000;

	for( size_t i = 0; i < y.size(); ++i ) {
		deviates[i] = ( y[i] - pr->f[i] ) / yErr[i];
	}

	return 0;
}


int OI(int m, int n, double *p, double *deviates, double **derivs, void *data) {
//	std::cout << "in function OI" << std::endl;
	struct Parameters *pr = (struct Parameters *) data;
	valarray<double> x( pr->x ),
					 y( pr->y ),
					 yErr( pr->yErr );

	double z = p[0],
		   ww = p[1],
		   OI_b = 6300.00 * ( 1 + z ),
		   OI_r = 6366.00 * ( 1 + z ),
		   OI_b_inten = p[2],
		   OI_r_inten = p[3],
		   cont = p[4];

	valarray<double> OI_b_ln = ( OI_b_inten / ww / piSqr ) * ( exp( -0.5 * ( x - OI_b ) * ( x - OI_b ) / ( ww * ww ) ) );
	valarray<double> OI_r_ln = ( OI_r_inten / ww / piSqr ) * ( exp( -0.5 * ( x - OI_r ) * ( x - OI_r ) / ( ww * ww ) ) );

	pr->f = cont + OI_b_ln + OI_r_ln;
	pr->f *= 10000;

	for( size_t i = 0; i < y.size(); ++i ) {
		deviates[i] = ( y[i] - pr->f[i] ) / yErr[i];
	}

	return 0;
}


int OI_br(int m, int n, double *p, double *deviates, double **derivs, void *data) {
//	std::cout << "in function OI_br" << std::endl;
	struct Parameters *pr = (struct Parameters *) data;
	valarray<double> x( pr->x ),
					 y( pr->y ),
					 yErr( pr->yErr );

	double z = p[0],
		   ww = p[1],
		   OI_b = 6300.00 * ( 1 + z ),
		   OI_r = 6366.00 * ( 1 + z ),
		   OI_b_inten = p[2],
		   OI_r_inten = p[3],
		   cont = p[4],
		   z_br = p[5],
		   ww_br = p[6],
		   OI_b_br = 6300.00 * ( 1 + z_br ),
		   OI_r_br = 6366.00 * ( 1 + z_br ),
		   OI_b_inten_br = p[2],
		   OI_r_inten_br = p[3];

	valarray<double> OI_b_ln = ( OI_b_inten / ww / piSqr ) * ( exp( -0.5 * ( x - OI_b ) * ( x - OI_b ) / ( ww * ww ) ) );
	valarray<double> OI_r_ln = ( OI_r_inten / ww / piSqr ) * ( exp( -0.5 * ( x - OI_r ) * ( x - OI_r ) / ( ww * ww ) ) );

	valarray<double> OI_b_br_ln = ( OI_b_inten_br / ww_br / piSqr ) * ( exp( -0.5 * ( x - OI_b_br ) * ( x - OI_b_br ) / ( ww_br * ww_br ) ) );
	valarray<double> OI_r_br_ln = ( OI_r_inten_br / ww_br / piSqr ) * ( exp( -0.5 * ( x - OI_r_br ) * ( x - OI_r_br ) / ( ww_br * ww_br ) ) );

	pr->f = cont + OI_b_ln + OI_r_ln + OI_b_br_ln + OI_r_br_ln;
	pr->f *= 10000;

	for( size_t i = 0; i < y.size(); ++i ) {
		deviates[i] = ( y[i] - pr->f[i] ) / yErr[i];
	}

	return 0;
}


int SII(int m, int n, double *p, double *deviates, double **derivs, void *data) {
	//std::cout << "in function SII" << std::endl;
	struct Parameters *pr = (struct Parameters *) data;
	valarray<double> x( pr->x ),
					 y( pr->y ),
					 yErr( pr->yErr );

	double z = p[0],
		   ww = p[1],
		   SII_b = 6716.00 * ( 1 + z ),
		   SII_r = 6731.00 * ( 1 + z ),
		   SII_b_inten = p[2],
		   SII_r_inten = p[3],
		   cont = p[4];

	valarray<double> SII_b_ln = ( SII_b_inten / ww / piSqr ) * ( exp( -0.5 * ( x - SII_b ) * ( x - SII_b ) / ( ww * ww ) ) );
	valarray<double> SII_r_ln = ( SII_r_inten / ww / piSqr ) * ( exp( -0.5 * ( x - SII_r ) * ( x - SII_r ) / ( ww * ww ) ) );

	pr->f = cont + SII_b_ln + SII_r_ln;
	pr->f *= 10000;

	for( size_t i = 0; i < y.size(); ++i ) {
		deviates[i] = ( y[i] - pr->f[i] ) / yErr[i];
	}

	return 0;
}



int SII_br(int m, int n, double *p, double *deviates, double **derivs, void *data) {
	//std::cout << "in function SII_br" << std::endl;
	struct Parameters *pr = (struct Parameters *) data;
	valarray<double> x( pr->x ),
					 y( pr->y ),
					 yErr( pr->yErr );

	double z = p[0],
		   ww = p[1],
		   SII_b = 6716.00 * ( 1 + z ),
		   SII_r = 6731.00 * ( 1 + z ),
		   SII_b_inten = p[2],
		   SII_r_inten = p[3],
		   cont = p[4],
		   z_br = p[5],
		   ww_br = p[6],
		   SII_b_br = 6716.00 * ( 1 + z_br ),
		   SII_r_br = 6731.00 * ( 1 + z_br ),
		   SII_b_inten_br = p[7],
		   SII_r_inten_br = p[8];


	valarray<double> SII_b_ln = ( SII_b_inten / ww / piSqr ) * ( exp( -0.5 * ( x - SII_b ) * ( x - SII_b ) / ( ww * ww ) ) );
	valarray<double> SII_r_ln = ( SII_r_inten / ww / piSqr ) * ( exp( -0.5 * ( x - SII_r ) * ( x - SII_r ) / ( ww * ww ) ) );

	valarray<double> SII_b_br_ln = ( SII_b_inten_br / ww_br / piSqr ) * ( exp( -0.5 * ( x - SII_b_br ) * ( x - SII_b_br ) / ( ww_br * ww_br ) ) );
	valarray<double> SII_r_br_ln = ( SII_r_inten_br / ww_br / piSqr ) * ( exp( -0.5 * ( x - SII_r_br ) * ( x - SII_r_br ) / ( ww_br * ww_br ) ) );

	pr->f = cont + SII_b_ln + SII_r_ln + SII_b_br_ln + SII_r_br_ln;
	pr->f *= 10000;

	for( size_t i = 0; i < y.size(); ++i ) {
		deviates[i] = ( y[i] - pr->f[i] ) / yErr[i];
	}

	return 0;
}




int triplet(int m, int n, double *p, double *deviates, double **derivs, void *data) {
	//std::cout << "in function triplet" << std::endl;
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

int triplet_br(int m, int n, double *p, double *deviates, double **derivs, void *data) {
	//std::cout << "in function triplet_br" << std::endl;
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

int Ha(int m, int n, double *p, double *deviates, double **derivs, void *data) {
	//std::cout << "in function Ha" << std::endl;
	struct Parameters *pr = (struct Parameters *) data;
	valarray<double> x( pr->x ),
					 y( pr->y ),
					 yErr( pr->yErr );

	double z = p[0],
		   ww = p[1],
		   Ha = 6562.8 * ( 1 + z ),
		   NII_b = 6548.1 * ( 1 + z ),
		   NII_r = 6583 * ( 1 + z ),
		   Hainten = p[2],
		   NIIinten = p[3],
		   cont = p[4];

	valarray<double> Ha_ln = ( Hainten / ww / piSqr ) * ( exp( -0.5 * ( x - Ha ) * ( x - Ha ) / ( ww * ww ) ) );
	valarray<double> NII_a_ln = ( NIIinten / ww / piSqr ) * ( exp( -0.5 * ( x - NII_r ) * ( x - NII_r ) / ( ww * ww ) ) );
	valarray<double> NII_b_ln = ( NIIinten / 3. / ww / piSqr ) * ( exp( -0.5 * ( x - NII_b ) * ( x - NII_b ) / ( ww * ww ) ) );

	pr->f = cont + Ha_ln + NII_a_ln + NII_b_ln;
	pr->f *= 10000;

	for( size_t i = 0; i < y.size(); ++i ) {
		deviates[i] = ( y[i] - pr->f[i] ) / yErr[i];
	}

	return 0;
}


int Ha_br(int m, int n, double *p, double *deviates, double **derivs, void *data) {
	//std::cout << "in function Ha_br" << std::endl;
	struct Parameters *pr = (struct Parameters *) data;
	valarray<double> x( pr->x ),
					 y( pr->y ),
					 yErr( pr->yErr );

	double z = p[0],
		   ww = p[1],
		   Ha = 6562.8 * ( 1 + z ),
		   NII_b = 6548.1 * ( 1 + z ),
		   NII_r = 6583 * ( 1 + z ),
		   Hainten = p[2],
		   NIIinten = p[3],
		   cont = p[4],
		   z_br = p[5],
		   ww_br = p[6],
		   Ha_br = 6562.8 * ( 1 + z_br ),
		   NII_b_br = 6548.1 * ( 1 + z_br ),
		   NII_r_br = 6583 * ( 1 + z_br ),
		   Hainten_br = p[7],
		   NIIinten_br = p[8];


	valarray<double> Ha_ln = ( Hainten / ww / piSqr ) * ( exp( -0.5 * ( x - Ha ) * ( x - Ha ) / ( ww * ww ) ) );
	valarray<double> NII_a_ln = ( NIIinten / ww / piSqr ) * ( exp( -0.5 * ( x - NII_r ) * ( x - NII_r ) / ( ww * ww ) ) );
	valarray<double> NII_b_ln = ( NIIinten / 3. / ww / piSqr ) * ( exp( -0.5 * ( x - NII_b ) * ( x - NII_b ) / ( ww * ww ) ) );

	valarray<double> Ha_br_ln = ( Hainten_br / ww_br / piSqr ) * ( exp( -0.5 * ( x - Ha_br ) * ( x - Ha_br ) / ( ww * ww ) ) );
	valarray<double> NII_a_br_ln = ( NIIinten_br / ww_br / piSqr ) * ( exp( -0.5 * ( x - NII_r_br ) * ( x - NII_r_br ) / ( ww * ww ) ) );
	valarray<double> NII_b_br_ln = ( NIIinten_br / 3. / ww_br / piSqr ) * ( exp( -0.5 * ( x - NII_b_br ) * ( x - NII_b_br ) / ( ww * ww ) ) );

	pr->f = cont + Ha_ln + NII_a_ln + NII_b_ln + Ha_br_ln + NII_a_br_ln + NII_b_br_ln;
	pr->f *= 10000;

	for( size_t i = 0; i < y.size(); ++i ) {
		deviates[i] = ( y[i] - pr->f[i] ) / yErr[i];
	}

	return 0;
}



int Cont(int m, int n, double *p, double *deviates, double **derivs, void *data) {
	struct Parameters *pr = (struct Parameters *) data;
	valarray<double> x( pr->x ),
					 y( pr->y ),
					 yErr( pr->yErr );
	double cont = p[0];

	pr->f = cont;
	pr->f *= 10000;

	for( size_t i = 0; i < y.size(); ++i ) {
		deviates[i] = ( y[i] - pr->f[i] ) / yErr[i];
	}

	return 0;

}

#endif /* FUNCTIONS_H_ */
