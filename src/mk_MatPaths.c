//Matthieu Landes, IPGP avril 2009

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FILE_MATG "matG.bin"

#define ABS(a) (((a)>0) ? (a):-(a))
enum {
	RIEN=0, BAS, DROITE, HAUT, GAUCHE 
};

double *xGrid;
int nX;
double *yGrid;
int nY;
double *G;
int nG;

typedef struct {
	double x;
	double y;
} POINT;

typedef struct {
	double a;
	double b;
	double c;
} LINE;

double G_xmin, G_xmax, G_ymin, G_ymax, G_dx, G_dy;


int isInGrid(POINT );
int isInRect(POINT pt, POINT , POINT );
void distForOnePath(POINT , POINT ); 
void calcDroite(LINE *, POINT, POINT);
void writeG(FILE*);
int CoordInter(double *, double *, LINE, LINE);
int SIAM (POINT, POINT, POINT );
int isIntersection ( POINT, POINT, POINT, POINT );

int main(int argc, char *argv[]) {
	FILE *flux, *flin;
	int i;
	POINT p1, p2;

	if( argc != 3 ) {
		printf("usage : <%s> filePath filegrid.\n", argv[0]);
		exit(0);
	}


	if( (flin=fopen(argv[1],"r"))==NULL) {
		printf("ERROR : cannot open %s.\n",argv[1]);
		exit(0);
	}

    if( (flux=fopen(argv[2],"r"))==NULL) {
		printf("ERROR : cannot open %s.\n",argv[2]);
		exit(0);
	}


	fscanf(flux, "%lf%lf", &G_xmin, &G_xmax);
	fscanf(flux, "%lf%lf", &G_ymin, &G_ymax);
	fscanf(flux, "%lf%lf", &G_dx, &G_dy);
	
	fclose(flux);

	printf("Lecture : x=(%f,%f), y=(%f,%f) and delta=(%f,%f)\n",
		G_xmin,G_xmax,G_ymin,G_ymax,G_dx,G_dy);

	//allocation memoire
	nX=(int)floor((G_xmax-G_xmin)/G_dx)+1;
	printf("Dim xGrid=%d\n",nX);
	if( (xGrid = (double*) malloc(sizeof(double)*(nX+1))) == NULL ) {
		printf("ERROR : not enough memory\n");
		exit(0);
	}

	nY=(int)floor((G_ymax-G_ymin)/G_dy)+1;
	printf("Dim yGrid=%d\n",nY);
	if( (yGrid = (double*) malloc(sizeof(double)*(nY+1))) == NULL ) {
		printf("ERROR : not enough memory\n");
		exit(0);
	}
	
	nG=nX*nY;
	if( (G=(double*)malloc(sizeof(double)*nG)) == NULL) {
		printf("ERROR : not enough memory\n");
		exit(0);
	}
//	printf("xGrid\n");
	for(i=0; i<=nX; i++) {xGrid[i]=G_xmin+i*G_dx;}// printf("%2.2lf ", xGrid[i]);}
//	printf("yGrid\n");
	for(i=0; i<=nY; i++) {yGrid[i]=G_ymin+i*G_dy;} //printf("%2.2lf ", yGrid[i]);}

	if( (flux=fopen(FILE_MATG, "wb"))==NULL) {
		printf("ERROR : cannot open %s.\n", FILE_MATG);
		exit(0);
	}


	while( fscanf(flin,"%lg %lg %lg %lg", &p1.x, &p1.y, &p2.x, &p2.y)!=EOF) {
		//printf("%lf %lf %lf %lf\n", p1.x, p1.y, p2.x, p2.y);
		//p1.x=0.5; p1.y=0.5;
		//p2.x=7.5; p2.y=3.0;
		for(i=0; i<nG; i++) G[i]=0.0;
		if( isInGrid(p1) && isInGrid(p2) ) {
			//calculer dt pour ce trajet
			distForOnePath(p1, p2);
		}
		else
			printf("Attention: (%.5f,%.5f) or (%.5f,%.5f) out of the grid\n", p1.x, p1.y, p2.x, p2.y);
		fwrite(G, sizeof(double), nG, flux);
//		writeG(flux);

	}
	
	fclose(flin);
	fclose(flux);
	free(xGrid);
	free(yGrid);
	free(G);
    
    return 1;

}



void distForOnePath(POINT p1, POINT p2) {
	int iX, iY;
	int incG, etat;
	double *ptG;
	double x,y;
	LINE dpath, dtmp;
	POINT bd,bg,hg,hd, ptCur, *p3, *p4;

	iX = (int)floor((p1.x-G_xmin)/G_dx);
	iY = (int)floor((p1.y-G_ymin)/G_dy);
	//printf("iX=%d, iY=%d : ",iX, iY);

	calcDroite(&dpath, p1, p2); 	

	ptCur.x=p1.x;
	ptCur.y=p1.y;


	ptG = &(G[iY*nX+iX]);

	bg.x=xGrid[iX]; bg.y=yGrid[iY];
	bd.x=xGrid[iX+1]; bd.y=yGrid[iY];
	hd.x=xGrid[iX+1]; hd.y=yGrid[iY+1];
	hg.x=xGrid[iX]; hg.y=yGrid[iY+1];
	
	etat=RIEN;

	while( !isInRect(p2, bg, hd) ) {
		//teste intersection avec segment du bas
		if( isIntersection(p1, p2, bg, bd) && (etat!=HAUT)) {
			iY--;
			incG=-nX;
			p3=&bg;
			p4=&bd;
			etat=BAS;
		}
		// avec le segment de droite
		else if( isIntersection(p1, p2, bd, hd) && (etat!=GAUCHE) ) {
			iX++;
			incG=1;
			p3=&bd;
			p4=&hd;
			etat=DROITE;
		}
		// avec segment du haut
		else if( isIntersection(p1, p2, hg, hd) && (etat!=BAS)) {
			iY++;
			incG=nX;
			p3=&hg;
			p4=&hd;
			etat=HAUT;
		}
		//avec celui de gauche
		else {
			iX--;
			incG=-1;
			p3=&bg;
			p4=&hg;
			etat=GAUCHE;
		}

		calcDroite(&dtmp, *p3, *p4);
		CoordInter(&x, &y, dpath, dtmp);
		//printf("(%d,%d) ", iX, iY);fflush(stdout);

		*ptG=sqrt( (x-ptCur.x)*(x-ptCur.x)+(y-ptCur.y)*(y-ptCur.y));

		//G[iX+iY*nX]=sqrt( (x-ptCur.x)*(x-ptCur.x)+(y-ptCur.y)*(y-ptCur.y));

		ptCur.x=x;
		ptCur.y=y;

		ptG+=incG;

		bg.x=xGrid[iX]; bg.y=yGrid[iY];
		bd.x=xGrid[iX+1]; bd.y=yGrid[iY];
		hd.x=xGrid[iX+1]; hd.y=yGrid[iY+1];
		hg.x=xGrid[iX]; hg.y=yGrid[iY+1];
	}
//	putchar('\n');
	*ptG=sqrt( (p2.x-ptCur.x)*(p2.x-ptCur.x)+(p2.y-ptCur.y)*(p2.y-ptCur.y));
}

void writeG(FILE *flux) {
	int i;
	for(i=0; i<nG; i++) fprintf(flux, "%f ", G[i]);
	fputc('\n', flux);
}
	
int CoordInter(double *x, double *y, LINE d1, LINE d2) {
	double det;

	det=d1.b*d2.a-d1.a*d2.b;
	if( ABS(det) < 1e-10) return 0;

	*x=(d1.c*d2.b-d1.b*d2.c)/det;
	*y=-(d1.c*d2.a-d1.a*d2.c)/det;
    
    return 1;
}


void calcDroite(LINE *d, POINT p1, POINT p2) {
	if( ABS(p1.y-p2.y) < 1e-10 ) {
		d->a=0.0;
		d->b=1.0;
		d->c=-p1.y;
	}
	else if( ABS(p1.x-p2.x) < 1e-10 ) {
		d->a=1.0;
		d->b=0.0;
		d->c=-p1.x;
	}
	else {
		d->a=1.0;
		d->b=-(p1.x-p2.x)/(p1.y-p2.y);
		d->c=(p1.x*p2.y-p2.x*p1.y)/(p1.y-p2.y);
	}
}



int isInRect(POINT pt, POINT basGauche, POINT hautDroit) {
	return ((pt.x>=basGauche.x) && (pt.x<=hautDroit.x) && 
		(pt.y>=basGauche.y) && (pt.y<=hautDroit.y));
}

/*
int isInGrid(POINT pt) {
	return ((pt.x <= G_xmax) && (pt.x >= G_xmin) && (pt.y <= G_ymax) && (pt.y >= G_ymin));
}*/

int isInGrid(POINT pt) {
	return ((pt.x <= xGrid[nX]) && (pt.x >= G_xmin) && (pt.y <= yGrid[nY]) && (pt.y >= G_ymin));
}


int SIAM (POINT p0, POINT p1, POINT p2 ) {
	double dx1,dx2,dy1,dy2;
	dx1=p1.x-p0.x;dy1=p1.y-p0.y;
	dx2=p2.x-p0.x;dy2=p2.y-p0.y;
	if (dx1*dy2>dy1*dx2) return 1;
	if (dx1*dy2<dy1*dx2) return -1;
	if((dx1*dx2<0)||(dy1*dy2<0))return -1;
	if ((dx1*dx1+dy1*dy1)<(dx2*dx2+dy2*dy2))return 1;
	return 0;
}


int isIntersection( POINT p1, POINT p2, POINT p3, POINT p4) {
	return (( SIAM(p1, p2, p3)*SIAM(p1, p2, p4))<=0)
		&&((SIAM(p3, p4, p1)*SIAM(p3, p4, p2))<=0);
}
	
