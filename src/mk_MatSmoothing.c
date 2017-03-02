// Matthieu Landes, IPGP, avril 2009
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define FILE_MATF "matF.bin"

#define ABS(a) (((a)>0) ? (a):-(a))


#define LAMBDA 2.


double *xGrid_center;
int nX;
double *yGrid_center;
int nY;

double *dist;
int n;

float *F;

double G_xmin, G_xmax, G_ymin, G_ymax, G_dx, G_dy;

int main(int argc, char *argv[]) {
	FILE *flux;
	int i, j, k, l, m, iX, iY;
	double *ptDouble;
	double norm, Lcorr2=LAMBDA, Sk, tmp;


	if( argc < 3) {
		printf("usage : <%s> Lcorr Gridfile\n", argv[0]);
		exit(0);
	}
    
    if( (flux=fopen(argv[2],"r"))==NULL) {
		printf("ERREUR impossible d'ouvrir le fichier de grille %s.\n",argv[2]);
		exit(0);
	}

	sscanf(argv[1],"%lf", &Lcorr2);
	printf("Lcorr=%f\n", Lcorr2);


	fscanf(flux, "%lf%lf", &G_xmin, &G_xmax);
    printf("lonmin=%lf lonmax=%lf\n",G_xmin, G_xmax);
	fscanf(flux, "%lf%lf", &G_ymin, &G_ymax);
    printf("latmin=%lf latmax=%lf\n",G_ymin, G_ymax);
	fscanf(flux, "%lf%lf", &G_dx, &G_dy);
	printf("binsize_lon=%lf binsizelat=%lf\n",G_dx, G_dy);
	fclose(flux);


	//allocation memoire
	nX=(int)floor((G_xmax-G_xmin)/G_dx)+1;
	printf("Dim xGrid=%d\n",nX);
	if( (xGrid_center = (double*) malloc(sizeof(double)*(nX+1))) == NULL ) {
		printf("ERREUR : pas assez dde memoire\n");
		exit(0);
	}

	nY=(int)floor((G_ymax-G_ymin)/G_dy)+1;
	printf("Dim yGrid=%d\n",nY);
	if( (yGrid_center = (double*) malloc(sizeof(double)*(nY+1))) == NULL ) {
		printf("ERREUR : pas assez dde memoire\n");
		exit(0);
	}

	n=nX*nY;
	if( ( dist = (double*) malloc(sizeof(double)*n)) == NULL ) {
		printf("ERREUR : pas assez dde memoire\n");
		exit(0);
	}

	if( ( F = (float*) malloc(sizeof(float)*n*n)) == NULL ) {
		printf("ERREUR : pas assez dde memoire\n");
		exit(0);
	}

	for(i=0; i<nX; i++) xGrid_center[i]=G_xmin+(i+1/2.)*G_dx;
	for(i=0; i<nY; i++) yGrid_center[i]=G_ymin+(i+1/2.)*G_dy;

	ptDouble=dist;

	for(i=0; i<nY; i++) for(j=0; j<nX; j++) {
		*ptDouble++=sqrt((xGrid_center[j]-xGrid_center[0])*(xGrid_center[j]-xGrid_center[0])
			+ (yGrid_center[i]-yGrid_center[0])* (yGrid_center[i]-yGrid_center[0]));
		//dist[i*nY+j]
	}

	printf("LCORR2 = %f",Lcorr2);
	for(i=0,iY=0; i<nY; i++) for(j=0; j<nX; j++, iY+=n) {
	
			norm=0.0;
			for(k=0,iX=0; k<nY; k++) for(l=0; l<nX; l++, iX++) {
				tmp=(float)ABS(dist[ABS(i-k)*nX+ABS(j-l)]);
				Sk=exp(-0.5*tmp/(Lcorr2));
				norm+=Sk;
				F[iX+iY]=(float)Sk;
			}
			//if( norm > 1e-10) for(m=0; m<n; m++) F[iY+m]/=norm;
	}

	for(m=0; m<n*n; m+=n+1) F[m]=1.0;


	if( (flux=fopen(FILE_MATF, "wb"))==NULL) {
		printf("ERREUR: impossible d'ouvrir %s.\n", FILE_MATF);
		exit(0);
	}

	fwrite(F, sizeof(float), n*n, flux);
    
	fclose(flux);
    printf("%i\n", sizeof(float));
    printf("n*n=%d\n",(n*n));
	free(xGrid_center); free(yGrid_center);
	free(dist);
	free(F);
    
    return 1;
}
