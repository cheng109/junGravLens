/*
 * Image.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: cheng109
 */

#include "Image.h"
#include "fitsio.h"
#include <string>
#include <vector>
#include <iostream>
#include "commons.h"
#include <iomanip>
#include <cmath>
#include <random>
//using namespace arma;
using namespace std;
#define buffsize 100000


Image::Image() {
	// TODO Auto-generated constructor stub

}


Image::Image(vector<double> xpos, vector<double> ypos, vector<double> *briList, long naxis1, long naxis2, int bitpix):
			naxis(2), naxis1(naxis1), naxis2(naxis2), bitpix(bitpix), data(naxis1*naxis2, 0) {
	npixels = naxis1*naxis2;

	long iList=0, x, y;
	for(int i=0; i< briList->size(); ++i) {
		x = nearbyint(xpos[i]);
		y = nearbyint(ypos[i]);

		if(x>0 && x< naxis1 && y>0 && y<naxis2) {
			iList = naxis1*y+x;
			data[iList] += briList->at(i);

		}

	}
}


/***************************
Function:   	Image Constructor
Description:    Construct an image object from given vectors of position and brightness;
Arguments:		(1): Pixel position along X direction;
				(2): Pixel position along Y direction;
				(3): Pixel brightness;
				(4)-(6): naxis1, naxis2, bitpix;
Returns:
****************************/
Image::Image(vector<int> xpos, vector<int> ypos, vector<double> *briList, long naxis1, long naxis2, int bitpix):
			naxis(2), naxis1(naxis1), naxis2(naxis2), bitpix(bitpix), data(naxis1*naxis2, 0) {

	npixels = naxis1*naxis2;
	long iList=0, x, y;
	for(int i=0; i< briList->size(); ++i) {
		x = xpos[i];
		y = ypos[i];

		if(x>=0 && x< naxis1 && y>=0 && y<naxis2) {
			iList = naxis1*y+x;
			data[iList] += briList->at(i);

		}
	}
}

/***************************
Function:   	Image Constructor
Description:    Construct an image object from a given FITS image file;
Arguments:		(1): Fits image name;
Returns:
Notes: 			This constructor is used to read the given lensed image,
				variance image, or PSF fits image;
****************************/
Image::Image(string imgFileName) {

	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status,  nfound, anynull;
	long naxes[2], fpixel, nbuffer;


	float nullval, buffer[buffsize];
	//char filename[]  = "atestfil.fit";     /* name of existing FITS file   */

	status = 0;

	if ( fits_open_file(&fptr, imgFileName.c_str(), READONLY, &status) )
		printerror( status );
	if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
		printerror( status );
	if(fits_get_img_type(fptr, &bitpix, &status))
		printerror( status );
	if (fits_get_img_dim(fptr, &naxis, &status))
		printerror( status );
	naxis1 = naxes[0] ;
	naxis2 = naxes[1] ;
	npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
	fpixel   = 1;
	nullval  = 0;                /* don't check for null values in the image */

	while (npixels > 0)
	{
		nbuffer = npixels;
		if (npixels > buffsize)
			nbuffer = buffsize;     /* read as many pixels as will fit in buffer */

		if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,
				buffer, &anynull, &status) )
			printerror( status );
		for(long i=0; i<nbuffer; ++i) {
			data.push_back(buffer[i]);
		}
		npixels -= nbuffer;    /* increment remaining number of pixels */
		fpixel  += nbuffer;    /* next pixel to be read in image */
	}

	npixels  = naxes[0] * naxes[1];
	if ( fits_close_file(fptr, &status) )
		printerror( status );
}



void Image::getConstants(long *len, long* naxis1, long* naxis2, double *res, int *bit){
	//cout << "length " << length << endl;
	*len=this->length;
	*naxis1 = this->naxis1;
	*naxis2 = this->naxis2;
	*res = this->res;
	*bit = this->bitpix;
}


/***************************
Function:   	printImageInfo
Description:    Print out some basic information of the FITS image, and value of the image;
Arguments:		(1): left down x1;
				(2): left down y1;
				(3): right up x2;
				(4): right up y2;
Returns:		void
Notes: 			Basic info: filename, naxis1, naxis2, bitpix, npixels and so on;
				Image info: left-down(x1,y1);  right-up(x2, y2);
****************************/
void Image::printImageInfo(int x1, int y1, int x2, int y2) {
	cout << "************************" << endl;
	cout << "File:   " << fileName << endl;
	cout << "NAXIS:  " << naxis << endl;
	cout << "NAXIS1: " << naxis1 << endl;
	cout << "NAXIS2: " << naxis2 << endl;
	cout << "Total:  " << npixels << endl;
	cout << "TYPE:   " << bitpix << endl;
	cout << "************************" << endl;

	//cout << counter << endl;
	x1 = (x1>0)?x1:0;
	y1 = (y1>0)?y1:0;
	x2 = (x2<naxis1)?x2:naxis1;
	y2 = (y2<naxis2)?y2:naxis2;
	// in order of DS9 show.
	for (int y=y2; y>y1; --y) {
		for(int x=x1; x<x2; ++x) {
			cout <<setprecision(4) << this->data[y*naxis1+x] << "\t";
		}
		cout << endl;

	}
	cout << endl;
	cout << data.size() << endl;
}


void Image::writeFilterImage(string imgFileName) {
	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status, ii;
	long  fpixel;
	long naxes[2] = { naxis1, naxis2 };   /* image is 300 pixels wide by 200 rows */
	remove(imgFileName.c_str());               /* Delete old file if it already exists */
	status = 0;         /* initialize status before calling fitsio routines */
	if (fits_create_file(&fptr, imgFileName.c_str(), &status)) /* create new FITS file */
		printerror( status );           /* call printerror if error occurs */
	if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
		printerror( status );
	double **array = (double **)malloc(naxis2* sizeof(double**));
	array[0] = (double *)malloc( naxis1 * naxis2* sizeof( double ) );
	for( ii=1; ii<naxis2; ii++ )
		array[ii] = array[ii-1] + naxis1;
	for (int i = 0; i < naxis2; i++) {
		for (int j = 0; j < naxis1; j++)
			array[i][j] = 0.0;
	}
	for( ii=1; ii<naxis2; ii++ )
		array[ii] = array[ii-1] + naxis1;
	for (int i=0; i<dataList.size();  ++i) {
		array[yList[i]][xList[i]] = dataList[i];
	}
	//cout << "After filter2:" <<  dataList.size() << endl;

	fpixel = 1;                               /* first pixel to write      */
	if (fits_write_img(fptr, TDOUBLE, fpixel, npixels, array[0], &status))
		printerror( status );
	free( array[0] );  /* free previously allocated memory */
	free( array);
	if ( fits_close_file(fptr, &status) )                /* close the file */
		printerror( status );


}



/***************************
Function:   	writeToFile
Description:    Write an image object to a fits file;
Arguments:		(1): Name of the fits file you want to output;
Notes: 			This function only write 'data' member of the image object;
				If you want to specify other data, please use the other 'writeFilterImage' function;
Returns:		void;
****************************/
void Image::writeToFile(string imgFileName) {
	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status, ii;
	long  fpixel;

	long naxes[2] = { naxis1, naxis2 };   /* image is 300 pixels wide by 200 rows */
	remove(imgFileName.c_str());               /* Delete old file if it already exists */
	status = 0;         /* initialize status before calling fitsio routines */
	if (fits_create_file(&fptr, imgFileName.c_str(), &status)) /* create new FITS file */
		printerror( status );           /* call printerror if error occurs */
	if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
		printerror( status );
	double **array = (double **)malloc(naxis2* sizeof(double**));
	array[0] = (double *)malloc( naxis1 * naxis2* sizeof( double ) );
	for( ii=1; ii<naxis2; ii++ )
	      array[ii] = array[ii-1] + naxis1;
	for (int i = 0; i < naxis2; i++) {
		for (int j = 0; j < naxis1; j++){
			array[i][j] = data[i*naxis1+j];
		}
	}
	fpixel = 1;                               /* first pixel to write      */
	/* write the array of unsigned integers to the FITS file */
	if (fits_write_img(fptr, TDOUBLE, fpixel, npixels, array[0], &status))
		printerror( status );

	free( array[0] );  /* free previously allocated memory */
	free( array);
	if ( fits_close_file(fptr, &status) )                /* close the file */
		printerror( status );
}


/***************************
Function:   	writeToFile
Description:    Write an image object to a fits file;
Arguments:		(1): Name of the fits file you want to output;
				(2): Mean of background; 
				(3): STD of background; 
Notes: 			This function only write 'data' member of the image object;
				If you want to specify other data, please use the other 'writeFilterImage' function;
Returns:		void;
****************************/

void Image::writeToFile(string imgFileName, double back_mean, double back_std) {
	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status, ii;
	long  fpixel;
	default_random_engine generator; 
	normal_distribution<double> distribution(back_mean, back_std); 


	long naxes[2] = { naxis1, naxis2 };   /* image is 300 pixels wide by 200 rows */
	remove(imgFileName.c_str());               /* Delete old file if it already exists */
	status = 0;         /* initialize status before calling fitsio routines */
	if (fits_create_file(&fptr, imgFileName.c_str(), &status)) /* create new FITS file */
		printerror( status );           /* call printerror if error occurs */
	if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
		printerror( status );
	double **array = (double **)malloc(naxis2* sizeof(double**));
	array[0] = (double *)malloc( naxis1 * naxis2* sizeof( double ) );
	for( ii=1; ii<naxis2; ii++ )
	      array[ii] = array[ii-1] + naxis1;
	for (int i = 0; i < naxis2; i++) {
		for (int j = 0; j < naxis1; j++){
			double number = distribution(generator); 
			array[i][j] = data[i*naxis1+j] + number ;
		}
	}
	fpixel = 1;                               /* first pixel to write      */
	/* write the array of unsigned integers to the FITS file */
	if (fits_write_img(fptr, TDOUBLE, fpixel, npixels, array[0], &status))
		printerror( status );

	free( array[0] );  /* free previously allocated memory */
	free( array);
	if ( fits_close_file(fptr, &status) )                /* close the file */
		printerror( status );
}

/***************************
Function:   	updateBackSubtract
Description:    subtract background from the regioin filtered image. 
Arguments:		(1): Mean of background; 
				(2): STD of background; 
Returns:		void
Notes: 			
****************************/

void Image::updateBackSubtract(double back_mean, double back_std) {

	std::default_random_engine generator; 
	std::normal_distribution<double> distribution(back_mean, back_std); 
	for (int i=0; i<dataList.size(); ++i ) {
		double number = distribution(generator); 
		dataList[i] -= number ; 
		d[i] -= number; 
	}

}

/***************************
Function:   	updateFilterImage
Description:    Filter image using DS9 region file.
Arguments:		(1): DS9 region file;
				(2): Flag = 0, 1
Returns:		void
Notes: 			Flag = 0: 'WITHOUT' filter;
				Flag = 1: 'WITH' filter (DS9 region file);
				regionType = 0:   Not supported yet; 
				regionType = 1:   Polygon region;
				regionType = 2:   Point region; 
****************************/
void Image::updateFilterImage(string regionFileName, int flag) {
	
	// Flag = 1, use region file to filter the image;
	// Flag = 0, no region filter.

	
	vector<double> xpos, ypos;
	if(flag==1) {  // with region; 
		string regionType = parseReagionFile(regionFileName, &xpos, &ypos);
		if(regionType == "polygon") {  // region type = polygon; 

			for(int i=0; i<naxis1*naxis2; ++i) {
				int y=i/naxis1;
				int x=i%naxis1;
				if(pnpoly(xpos.size(), &xpos, &ypos,  double(x+0.5) ,  double(y+0.5))) {
					dataList.push_back(data[i]);
					xList.push_back(x);
					yList.push_back(y);
				}
			}	
		}
		
		else if(regionType == "point" ) {
			for(int i=0; i< xpos.size(); ++i) {
				xList.push_back(int(xpos[i])); 
				yList.push_back(int(ypos[i])); 
				int index = int(ypos[i]) * naxis1 + int(xpos[i]);  
				dataList.push_back(data[index]); 
				
			}
		}		
		else if(regionType == "box") {

			// xpos and ypos contain the same values:  [centerx, centery, sizex, sizey, rotation] // rotation will be always '0'
			double boxCenterX 	= xpos[0]; 
			double boxCenterY 	= xpos[1]; 
			double boxSizeX 	= xpos[2]; 
			double boxSizeY 	= xpos[3]; 
			double boxRotation 	= xpos[4]; 
			for(int i=0; i<naxis1; ++i) {
				for (int j=0; j<naxis2; ++j) {
					int dy= j - boxCenterY;
					int dx= i - boxCenterX;
					if(     dx <  0.5*boxSizeX and dy <  0.5*boxSizeY 
						and dx > -0.5*boxSizeX and dy > -0.5*boxSizeY ) {
						dataList.push_back(data[j*naxis1 + i]);
						xList.push_back(i);
						yList.push_back(j);
					}
				}
			}	

		}
		else if(regionType == "circle") {

			// xpos and ypos contain the same values:  [centerx, centery, sizex, sizey, rotation] // rotation will be always '0'
			double circleCenterX 	= xpos[0]; 
			double circleCenterY 	= xpos[1]; 
			double circleR 			= xpos[2]; 

			for(int i=0; i<naxis1; ++i) {
				for (int j=0; j<naxis2; ++j) {
					int dy= j - circleCenterY;
					int dx= i - circleCenterX;
					double dr2 = dx*dx + dy*dy; 
					if( dr2 < circleR*circleR) {
						dataList.push_back(data[j*naxis1 + i]);
						xList.push_back(i);
						yList.push_back(j);
					}
				}
			}	

		}

		else if(regionType == "ellipse") {

			// xpos and ypos contain the same values:  [centerx, centery, sizex, sizey, rotation] // rotation will be always '0'
			double ellipseCenterX 	= xpos[0]; 
			double ellipseCenterY 	= xpos[1]; 
			double a				= xpos[2]; 
			double b				= xpos[3]; 
			double A 				= xpos[4];   // rotation angle; 


			for (int i=0; i<5; ++i) {


				cout << xpos[i] << endl; 
			}

			for(int i=0; i<naxis1; ++i) {
				for (int j=0; j<naxis2; ++j) {
					int dx= i - ellipseCenterX;
					int dy= j - ellipseCenterY;

					double dr2 = ((dx*cos(A) + dy*sin(A))/a)*((dx*cos(A) + dy*sin(A))/a) 
									+ ((dx*sin(A) - dy*cos(A))/b)*((dx*sin(A) - dy*cos(A))/b); 
					if( dr2 < 1) {
						dataList.push_back(data[j*naxis1 + i]);
						xList.push_back(i);
						yList.push_back(j);
					}
				}
			}	

		}

		else {
			cout << "Region type is not supported! " << endl; 
			//return 1; 
		}
	}
	else if(flag==0) {  // without region; 
		for(int i=0; i<naxis1*naxis2; ++i) {
			int y=i/naxis1;
			int x=i%naxis1;
			dataList.push_back(data[i]);
			xList.push_back(x);
			yList.push_back(y);

		}

	}

	else 
		cout << "Wrong region flag is used:  " << flag << endl; 
	length = dataList.size() ;
	//cout << "length = " << length << endl;
	for(int i=0; i<dataList.size(); ++i) {
		iList.push_back(i);
	}


	d = cV_to_eigenV (&dataList); 
}


vec Image::getMatrixD() {
	vec d(length);
	for(int i=0; i<length; ++i) {
		d[i]= dataList[i];
	}
	return d;
}


/***************************
Function:   	updateGridPointType
Description:    Classify the grid points as node or interpolate.
Arguments:		void
Notes: 			'1' indicates node point.
				'0' indicates interpolate point.
Returns:		void;
****************************/
void Image::updateGridPointType() {
	for(int i=0; i<length; ++i) {
		if((xList[i]+yList[i])%2==0)
			type.push_back(1); // should be zero; 
		else
			type.push_back(1);
	}

}

sp_mat Image::getVarMatrix()  {

	//this->updateFilterImage(regionFileName);
	long n= length;
	//n = 100;
	sp_mat invC(n, n);
	invC.reserve(n);
	for(int i=0; i<n; ++i) {
		invC.insert(i,i)= 1.0/varList[i];
		//cout << dataList[i] << endl;
	}
	return invC;
}

void Image::updateVarList(double threshold, double back_mean, double back_std) {
	double back_var = back_std*back_std ; 
	for(int i=0; i<length; ++i) {
		if(dataList[i]<threshold) {
			varList.push_back(back_var);
		}
		else
			varList.push_back(back_var);
			//varList.push_back(dataList[i]);
	}
}

void Image::updateVarList(string varFileName, string regionFileName) {
	Image varImage(varFileName);
	varImage.updateFilterImage(regionFileName, 0);
	for(int i=0; i<length; ++i) {
		varList.push_back(dataList[i]);
	}
}

void Image::normalizeData() {
	double sum=0;
	for(int i=0; i<length; ++i)
		sum += data[i];
	cout << "normSum: " << sum << endl;
	for(int i=0; i<length; ++i)
		data[i] = data[i]/sum;
}





Image::~Image() {
	// TODO Auto-generated destructor stub
}
