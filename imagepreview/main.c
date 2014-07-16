//
//  main.c
//  imagepreview
//
//  Created by Eduardo on 06/02/2013.
//  Copyright (c) 2013 Institute of Astronomy. All rights reserved.
//
// gcc main.c torben.c coords.c fopen.c -o preview -I/usr/local/include/wcslib -lwcs -lcurl -lcfitsio -lcpgplot

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "fitsio.h"
#include "cpgplot.h"
#include "wcs.h"

#ifndef NOCURL
#include "wcs.h"
#include <curl/curl.h>
#endif

#define TWOMASS_URL "http://casu.ast.cam.ac.uk/vistasp/conesearch/twomass?ra=%f&dec=%f&rad=%f"
#define SDSS_URL "http://casu.ast.cam.ac.uk/vistasp/conesearch/sdss?ra=%f&dec=%f&rad=%f"

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )

float torben(float a[], int n) ;
float mad(float a[], int n) ;
float *zscale(float m[], int n);
double str2ra (const char *in);
double str2dec (const char *in);

char *replace_str(char *str, char *orig, char *rep)
{
    static char buffer[4096];
    char *p;
    
    if(!(p = strstr(str, orig)))  // Is 'orig' even in 'str'?
        return str;
    
    strncpy(buffer, str, p-str); // Copy characters from 'str' start to 'orig' st$
    buffer[p-str] = '\0';
    
    sprintf(buffer+(p-str), "%s%s", rep, p+strlen(orig));
    
    return buffer;
}

char *strip_str(char *str)
{
	static char buffer[4096];
	char *p;
    
	if(!(p = strstr(str, "[")))  // Is 'orig' even in 'str'?
        return str;
	
	strncpy(buffer, str, p-str);
	buffer[p-str] = '\0';
	
	return buffer;
}

float *get_section(char *str)
{
	static char buffer[80];
	static char x1a[80], x1b[80];
	static float xx[2] = {0.0, 0.0};
	
	char *p, *p1, *p2;
    
	if(!(p = strstr(str, "[")))
		return xx;
    
	strcpy(buffer, p);
	buffer[*p] = '\0';
	
	if(!(p1 = strstr(buffer, ",")))
		return xx;
    
	strncpy(x1a, buffer+1, p1-buffer);
	x1a[p1-buffer-1] = '\0';
	
	if(!(p2 = strstr(x1a, ":")))
		return xx;
	
	strncpy(x1b, x1a, p2-x1a);
	x1b[p2-x1a] = '\0';
    
	xx[0]=atof(x1b)-1;
	
	if(!(p1 = strstr(buffer, ",")))
		return xx;
    
	strcpy(x1a, p1+1);
	x1a[*p1] = '\0';
	
	if(!(p2 = strstr(x1a, ":")))
		return xx;
    
	strncpy(x1b, x1a, p2-x1a);
	x1b[p2-x1a] = '\0';
	
	xx[1]=atof(x1b)-1;
    
	return xx;
}



char *trimwhitespace(char *str)
{
    char *end;
    
    // Trim leading space
    while(isspace(*str)) str++;
    
    // Trim trailing space
    end = str + strlen(str) - 1;
    while(end > str && isspace(*end)) end--;
    
    // Write new null terminator
    *(end+1) = 0;
    
    return str;
    
}

double *get_radec(char *str)
{
	static char buffer[256];
	static char x1a[256], x1[256];
	static double xx[2] = {0.0, 0.0};
	
	char *p, *p1, *p2;
	
	while(isspace(*str)) str++;
	p = strstr(str, " ");
	strncpy(buffer, str, p-str);
	buffer[p-str] = '\0';
	xx[0] =  (double) atof(buffer);
    
	if (!(p = strstr(str, " ")))
		return xx;
    
	while(isspace(*p)) p++;
	strcpy(x1a, p);
	x1a[*p] = '\0';
    
	if ((p = strstr(x1a, " ")))
	{
		strncpy(buffer, x1a, p-x1a);
		buffer[p-x1a] = '\0';
		xx[1] =  (double) atof(buffer);
	} else {
		xx[1] =  (double) atof(x1a);
	}
	
	return xx;
}

#ifndef NOCURL

enum fcurl_type_e { CFTYPE_NONE=0, CFTYPE_FILE=1, CFTYPE_CURL=2 };

struct fcurl_data
{
    enum fcurl_type_e type;     /* type of handle */
    union {
        CURL *curl;
        FILE *file;
    } handle;                   /* handle */
    
    char *buffer;               /* buffer to store cached data*/
    int buffer_len;             /* currently allocated buffers length */
    int buffer_pos;             /* end of data in buffer*/
    int still_running;          /* Is background url fetch still in progress */
};

typedef struct fcurl_data URL_FILE;

URL_FILE *url_fopen(const char *url,const char *operation);

#endif

int main (int argc, char *argv[]) {
	int c, j;
	long i, n;
	int dozscale=0, catalogue=0, interactive=0, pawprint=0, hdunum, isclassified=1;
	unsigned int iseed = (unsigned int)time(NULL);
	char *p;
	int pawnum[16] = {13, 14, 15, 16, 9, 10, 11, 12, 5, 6, 7, 8, 1, 2, 3, 4};
	int twomass=0, sdss=0;
	char rastr[32], decstr[32];
	
	/* 2MASS */
#ifndef NOCURL
	
	URL_FILE *handle;
	FILE *outf;
    
#endif
    
	int nread;
	char buffer[256];
	char url[512];
	char *urlpath;
	
	/* WCSLIB */
	double cd[2][2], cd11, cd22, cd12, cd21, pv21, pv22, pv23, pv24, pv25;
	float x, y, distance = 0.0, pixscl;
	float nxpix, nypix, crpix1, crpix2;
#ifndef NOCURL
	struct wcsprm *wcs;
	struct pvcard pv[5];
#endif
	char ctype[2][9] = {"RA---ZPN", "DEC--ZPN"};
	char cunit[2][9] = {"deg", "deg"};
	double xy[2],std[2],phi,theta,radec[2],crval1, crval2;
	double *dox;
	int stat;
    
	/* CFITSIO */
	fitsfile *infptr, *catfptr;
	int status = 0, ii = 1, iteration = 0, single = 0, hdupos, thdupos, colnum;
    int hdutype, bitpix, bytepix, naxis = 0, nkeys, datatype = 0, anynul;
	
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
	long nrows;
	
    long first, totpix = 0, npix;
    float *array, *ranarray, bscale = 1.0, bzero = 0.0, nulval = 0.;
	float *xcoord, *ycoord, *classification, *ellipticity, *posang, *gaussian;
	float skylevel, skynoise;
	float *zs;
    char comment[81];
	char instrument[20];
	
	/* PGPLOT */
	int symbol=4;
	float tr[6] = {0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
	char *device = "/xserve", chout[10], section=0;
	float z1, z2, width=9.0, sigma=10.0, cheight=2.0, radius=20.0;
	float gl[2] = {0.0, 1.0};
	float gr[2] = {0.0, 1.0};
	float gg[2] = {0.0, 1.0};
	float gb[2] = {0.0, 1.0};
	int x1, x2, y1, y2, dx, dy;
	float *ox, xout, yout;
	float xe[60], ye[60];
	
	if (argc < 2) {
		printf("Usage:\n");
		printf("\n");
		printf("    preview image.fit+5\n");
		printf("\n");
		printf("Options:\n\n");
		printf("  -a 2mass/sdss : query 2mass or sdss archive [experimental]\n");
		printf("  -c            : plots sources from catalogue if present\n");
		printf("  -d /xserve    : output graphics device\n");
		printf("  -h 1          : sets the symbol size when overlaying a catalogue\n");
		printf("  -i            : returns (x,y) on cursor input\n");
		printf("  -p            : plots all 16 chips from VISTA\n");
		printf("  -s 4          : sets the symbol type when overlaying a catalogue\n");
	    printf("  -t 10         : adjust the contrast (sigmas around background)\n");
		printf("  -w 9          : sets the size of the output image\n");
		printf("  -x bl         : displays only a section [bl, tl, tr, br, cc]\n");
		printf("  -z            : autoscale the image [default is to read SKYLEVEL, \n");
		printf("                                               SKYNOISE from header]\n");
		
		printf("\n");
		printf("Examples:\n");
		printf("\n");
		printf("    preview -h 1 -c -w 6 v20091103_00368_st.fit+12\n");
		printf("\n");
		return(0);
	}
	
	
	while ((c = getopt (argc, argv, "a:ch:izpd:s:t:w:x:")) != -1)
        switch(c) {
            case 'a':
                if (strstr(optarg, "2mass")) twomass=1;
                if (strstr(optarg, "sdss")) sdss=1;
                break;
            case 'i':
                interactive=1;
                break;
            case 'c':
                catalogue=1;
                break;
            case 'h':
                cheight=atof(optarg);
                break;
            case 'd':
                device = optarg;
                break;
            case 'p':
                pawprint=1;
                break;
            case 's':
                symbol=atoi(optarg);
                break;
            case 't':
                sigma = atof(optarg);
                break;
            case 'w':
                width = atof(optarg);
                break;
            case 'x':
                if (p = strstr(optarg, "bl")) section=1;
                if (p = strstr(optarg, "tl")) section=2;
                if (p = strstr(optarg, "tr")) section=3;
                if (p = strstr(optarg, "br")) section=4;
                if (p = strstr(optarg, "cc")) section=5;
                break;
            case 'z':
                dozscale=1;
                break;
            default:
                abort();
        }
	
	if (optind == argc) {
		printf("Image file not set. Type 'preview' for usage instructions.\n");
		return(1);
	}
	
	/* Open input file for read */
    if (fits_open_image(&infptr, argv[optind], READONLY, &status)) {
        fits_report_error(stderr, status);
        return(status);
    };
    
	fits_get_hdu_num(infptr, &hdupos);  /* Get the current HDU position */
	fits_get_num_hdus(infptr, &hdunum, &status);
	
	
	
	cpgopen(device);
	cpgpap(width, 1.0);
	cpgvsiz(0,width,0,width);
	
	if (pawprint) {
		fits_movabs_hdu(infptr, 1, &hdutype, &status);
		fits_read_key(infptr, TSTRING, "INSTRUME", &instrument, comment, &status);
		if (strstr(instrument, "VIRCAM")) {
			cpgsubp(4,4);
		} else if (strstr(instrument, "WFCAM")) {
			pawnum[0]=1;
			pawnum[1]=2;
			pawnum[2]=3;
			pawnum[3]=4;
			cpgsubp(2,2);
		} else if (strstr(instrument, "WFC")) {
			pawnum[0]=1;
			pawnum[1]=2;
			pawnum[2]=3;
			pawnum[3]=4;
			cpgsubp(2,2);
		} else if (strstr(instrument, "MOSAIC")) {
			pawnum[0]=1;
			pawnum[1]=2;
			pawnum[2]=3;
			pawnum[3]=4;
			pawnum[4]=5;
			pawnum[5]=6;
			pawnum[6]=7;
			pawnum[7]=8;
			cpgsubp(4,2);
		} else if (strstr(instrument, "SuprimeCam")) {
			pawnum[0]=1;
			pawnum[1]=2;
			pawnum[2]=3;
			pawnum[3]=4;
			pawnum[4]=5;
			pawnum[5]=6;
			pawnum[6]=7;
			pawnum[7]=8;
			pawnum[8]=9;
			pawnum[9]=10;
			cpgsubp(5,2);
		} else {
			printf("Instrument: %s %d\n", instrument, hdunum);
		}
	}
    
	
	if (hdunum==1) hdunum++;
	for (hdupos=0; hdupos<hdunum-1; hdupos++) {
		if (pawprint) {
			cpgpage();
			fits_movabs_hdu(infptr, pawnum[hdupos]+1, &hdutype, &status);
		}
        
		fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
        
		totpix = naxes[0] * naxes[1];
		npix = totpix;
		datatype=TFLOAT;
		first=1;
		//array = (float *) calloc(npix, datatype);
		array = (float *) malloc(totpix * sizeof(float));
		fits_read_img(infptr, datatype, first, npix, &nulval, array, &anynul, &status);
		fits_read_key(infptr, TFLOAT, "SKYLEVEL", &skylevel, comment, &status);
		fits_read_key(infptr, TFLOAT, "SKYNOISE", &skynoise, comment, &status);
		
		if (status || dozscale) {
			srand (iseed);
			ranarray = (float *) calloc(1000, datatype);
			zs = (float *) calloc(2, datatype);
			for (i=0; i<1000; i++) {
				j = rand()*1.0*npix/RAND_MAX;
				ranarray[i]=array[j];
			}
			zs=zscale(ranarray, 1000);
			skylevel = zs[0];
			skynoise = zs[1];
			//skylevel = torben(ranarray, 1000);
			//skynoise = mad(ranarray, 1000);
			z1 = skylevel - sigma * skynoise / 1.2;
			z2 = skylevel + sigma * skynoise;
			status=0;
			fits_get_hdu_num(infptr, &thdupos);
			printf("HDU %d - Median: %f Mad: %f\n", thdupos-1, zs[0], zs[1]);
		} else {
			z1 = skylevel - sigma * skynoise / 1.2;
			z2 = skylevel + sigma * skynoise;
		}
		
		if (catalogue) {
			if(fits_open_table(&catfptr, strip_str(replace_str(argv[optind], ".fit", "_cat.fits")), READONLY, &status)) {
				//fits_report_error(stderr, status);
		        //return(status);
                status=0;
                continue;
			}
			fits_get_hdu_num(infptr, &j);
			if (j==1) { j=2; }
			fits_movabs_hdu(catfptr, j, &hdutype, &status);
			fits_get_num_rows(catfptr, &nrows, &status);
			xcoord = (float *) calloc(nrows, TFLOAT);
			ycoord = (float *) calloc(nrows, TFLOAT);
			classification = (float *) calloc(nrows, TFLOAT);
			gaussian = (float *) calloc(nrows, TFLOAT);
			ellipticity = (float *) calloc(nrows, TFLOAT);
			posang = (float *) calloc(nrows, TFLOAT);
			fits_get_colnum(catfptr, 0, "x_coordinate", &colnum, &status);
			fits_read_col(catfptr, TFLOAT, colnum, 1, 1, nrows, &nulval, xcoord, &anynul, &status);
			fits_get_colnum(catfptr, 0, "y_coordinate", &colnum, &status);
			fits_read_col(catfptr, TFLOAT, colnum, 1, 1, nrows, &nulval, ycoord, &anynul, &status);
			fits_get_colnum(catfptr, 0, "classification", &colnum, &status);
			if (status) {
				status=0;
				isclassified=0;
			} else {
				fits_read_col(catfptr, TFLOAT, colnum, 1, 1, nrows, &nulval, classification, &anynul, &status);
			}
			fits_get_colnum(catfptr, 0, "gaussian_sigma", &colnum, &status);
			fits_read_col(catfptr, TFLOAT, colnum, 1, 1, nrows, &nulval, gaussian, &anynul, &status);
			fits_get_colnum(catfptr, 0, "ellipticity", &colnum, &status);
			fits_read_col(catfptr, TFLOAT, colnum, 1, 1, nrows, &nulval, ellipticity, &anynul, &status);
			fits_get_colnum(catfptr, 0, "position_angle", &colnum, &status);
			fits_read_col(catfptr, TFLOAT, colnum, 1, 1, nrows, &nulval, posang, &anynul, &status);
			fits_close_file(catfptr, &status);
		}
        
		x1 = 1;
		x2 = naxes[0];
		y1 = 1;
		y2 = naxes[1];
		dx = x2-x1+1;
		dy = y2-y1+1;
		
		switch (section) {
			case(1):
                x1=1;
                x2=600;
                y1=1;
                y2=600;
                break;
			case(2):
                x1=1.0;
                x2=600.;
                y1=naxes[1]-600.0;
                y2=naxes[1];
                break;
			case(3):
                x1=naxes[0]-600.0;
                x2=naxes[0];
                y1=naxes[1]-600.0;
                y2=naxes[1];
                break;
			case(4):
                x1=naxes[0]-600.0;
                x2=naxes[0]-1;
                y1=1.0;
                y2=600.0;
                break;
			case(5):
                x1=naxes[0]/2-300.0;
                x2=naxes[0]/2+300.0;
                y1=naxes[1]/2-300.0;
                y2=naxes[1]/2+300.0;
                break;
		}
        
		
		cpgwnad(x1,x2,y1,y2);
		cpgctab(gl, gr, gg, gb, 2, 1.5, 0.5);
		cpgimag(array, dx, dy, x1, x2, y1, y2, z1, z2, tr);
        
		if (catalogue) {
			cpgbbuf();
			ox=get_section(argv[optind]);
			cpgsci(2);
			cpgsch(cheight);
			for (i=0; i<nrows; i++) {
				if (isclassified) {
					if (classification[i]==-1 || classification[i]==-2){
						cpgsci(4);
					} else if (classification[i]==0) {
						cpgsci(2);
					} else if (classification[i]==1 || classification[i]==2) {
						cpgsci(3);
					}
				} else {
					if (ellipticity[i]<=0.2) {
						cpgsci(4);
					} else if (ellipticity[i]>0.4) {
						cpgsci(2);
					}
				}
				//cpgpt1(xcoord[i]-ox[0], ycoord[i]-ox[1], symbol);
				for (j=0; j<60; j++) {
					xe[j] = xcoord[i]-ox[0] + 2.4*cheight*gaussian[i]*cos(6*j/180.0*3.141592)*cos(posang[i]/180.*3.141592) - 2.4*cheight*gaussian[i]*(1-ellipticity[i])*sin(6*j/180.0*3.141592)*sin(posang[i]/180.*3.141592);
					ye[j] = ycoord[i]-ox[1] + 2.4*cheight*gaussian[i]*cos(6*j/180.0*3.141592)*sin(posang[i]/180.*3.141592) + 2.4*cheight*gaussian[i]*(1-ellipticity[i])*sin(6*j/180.0*3.141592)*cos(posang[i]/180.*3.141592);
				}
				cpgline(60, xe, ye);
				
			}
			cpgebuf();
			//cpgpt(nrows, xcoord, ycoord, symbol);
			
		}
        
#ifndef NOCURL
        
		if (twomass || sdss) {
			status=0;
			fits_read_key(infptr, TDOUBLE, "CRVAL1", &crval1, comment, &status);
			fits_read_key(infptr, TDOUBLE, "CRVAL2", &crval2, comment, &status);
			fits_read_key(infptr, TFLOAT, "CRPIX1", &crpix1, comment, &status);
			fits_read_key(infptr, TFLOAT, "CRPIX2", &crpix2, comment, &status);
			fits_read_key(infptr, TDOUBLE, "CD1_1", &cd11, comment, &status);
			fits_read_key(infptr, TDOUBLE, "CD1_2", &cd12, comment, &status);
			fits_read_key(infptr, TDOUBLE, "CD2_2", &cd22, comment, &status);
			fits_read_key(infptr, TDOUBLE, "CD2_1", &cd21, comment, &status);
			status=0;
            
			pixscl=cd11*cd22-cd12*cd21;
			if (pixscl<0) pixscl=-pixscl;
			pixscl=sqrt(pixscl)*3600.0;
            
			radius = max(x2, y2) * pixscl /60.0/1.5;
			
			wcs = malloc(sizeof(struct wcsprm));
			wcs->flag=-1;
			wcsini(1, 2, wcs);
			wcsnpv(5);
            
			strcpy(wcs->ctype[0], &ctype[0][0]);
			strcpy(wcs->ctype[1], &ctype[1][0]);
			strcpy(wcs->cunit[0], &cunit[0][0]);
			strcpy(wcs->cunit[1], &cunit[1][0]);
            
			wcs->crval[0] = crval1;
			wcs->crval[1] = crval2;
			wcs->crpix[0] = ((double) crpix1);
			wcs->crpix[1] = ((double) crpix2);
            
			cd[0][0] = ((double) cd11);
			cd[0][1] = ((double) cd12);
			cd[1][0] = ((double) cd21);
			cd[1][1] = ((double) cd22);
            
			wcs->altlin|=2;
			wcs->cd = *cd;
            
            
			pv[0].i=2;
			pv[0].m=1;
			pv[0].value=1.0;
			fits_read_key(infptr, TDOUBLE, "PV2_1", &pv21, comment, &status);
			if (status==0) {
				pv[0].value=pv21;
			} else {
				status=0;
			}
            
			pv[1].i=2;
			pv[1].m=2;
			pv[1].value=0.0;
			fits_read_key(infptr, TDOUBLE, "PV2_2", &pv22, comment, &status);
			if (status==0) {
				pv[1].value=pv22;
			} else {
				status=0;
			}
            
			pv[2].i=2;
			pv[2].m=3;
			pv[2].value = (double) 0.0;
			fits_read_key(infptr, TDOUBLE, "PV2_3", &pv23, comment, &status);
			if (status==0) {
				pv[2].value=pv23;
			} else {
				status=0;
			}
            
			pv[3].i=2;
			pv[3].m=4;
			pv[3].value=0.0;
			fits_read_key(infptr, TDOUBLE, "PV2_4", &pv24, comment, &status);
			if (status==0) {
				pv[3].value=pv24;
			} else {
				status=0;
			}
			
			pv[4].i=2;
			pv[4].m=5;
			pv[4].value = (double) 0.0;
			fits_read_key(infptr, TDOUBLE, "PV2_5", &pv25, comment, &status);
			if (status==0) {
				pv[4].value=pv25;
			} else {
				status=0;
			}
            
			wcs->npv=5;
			wcs->pv[0] = pv[0];
			wcs->pv[1] = pv[1];
			wcs->pv[2] = pv[2];
			wcs->pv[3] = pv[3];
			wcs->pv[4] = pv[4];
            
			(void)wcsset(wcs);
            
			//printf("Querying 2MASS....\n");
			xy[0]=x2/2.0;
			xy[1]=y2/2.0;
			
			(void)wcsp2s(wcs,1,2,xy,std,&phi,&theta,radec,&stat);
			
			if (twomass) {
				if(!(urlpath = getenv("TWOMASS_URL"))) urlpath=TWOMASS_URL;
				
			} else if (sdss) {
				if(!(urlpath = getenv("SDSS_URL"))) urlpath=SDSS_URL;
			}
			
			sprintf(url, urlpath, radec[0], radec[1], radius);
			handle = url_fopen(url, "r");
			
			cpgbbuf();
			cpgsci(3);
			cpgsch(cheight);
			while(!url_feof(handle))
            {
                url_fgets(buffer,sizeof(buffer),handle);
                
                if( strstr(buffer, "#")) continue;
                if( strstr(buffer, "RAJ")) continue;
                if( strstr(buffer, "---")) continue;
                if( strstr(buffer, "   ")) continue;
                if( strstr(buffer, "2MASS")) continue;
                if( strstr(buffer, "deg")) continue;
                if ( strlen(buffer)<10 ) continue;
                if( strstr(buffer, "xmlns")) break;
                dox=get_radec(buffer);
                radec[0]=dox[0];
                radec[1]=dox[1];
                (void)wcss2p(wcs,1,2,radec,&phi,&theta,std,xy,&stat);
                
                cpgpt1(xy[0],xy[1], symbol);
                
                //printf("%f %f\n", x, y);
            }
            
			cpgebuf();
            
			wcsfree(wcs);
			free(wcs);
		}
        
#endif
		
		free(array);
		
		if (!pawprint) break;
	}
	
	fits_close_file(infptr, &status);
	
	if (interactive & !pawprint) {
		ox=get_section(argv[optind]);
		while(!(p = strstr(chout, "q"))) {
			cpgband(7, 1, dx/2.0, dy/2.0, &xout, &yout, chout);
			printf("x = %d, y = %d\n", (int )(xout+ox[0]), (int )(yout+ox[1]) );
		}
	}
	
	cpgend();
	//cpgclos();
	
	/* if error occurred, print out error message */
	if (status) fits_report_error(stderr, status);
	
	return(0);
}
