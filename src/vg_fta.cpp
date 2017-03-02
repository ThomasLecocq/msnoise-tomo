/*********************************************
 * auteur : Matthieu LANDES
 * mail   : landes@ipgp.fr
 * date :   12/09/2013
 * Institut de Physique du Globe de Paris
 * Equipe de sismologie
 ********************************************/

#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>

#include "libfta.h"
#include "readsac.h"
#include "configparser.h"
#include "fta_param.h"

using namespace std;


void usage() {
	cerr << "---------" << endl;
	cerr << "Options:" << endl;
	cerr << "fmin|tmax =" << endl;
	cerr << "fmax|tmin =" << endl;
cerr << "bmin=, bmax= | b= (0.1)" << endl;
	cerr << "nfreq= (100)" << endl;
	cerr << "diag=[FP][TV] (PT)" << endl;
	cerr << "disp=no|max|all|cont|pt[:tp_i:vgtg_i]*i  (cont)" << endl;
	cerr << "tinit|finit = (at maximum)" << endl;
	cerr << "vginit|tginit= (all signal)" << endl;
	cerr << "ampMin= (no)" << endl;
//	cerr << "verbose=[01]" << endl;
	cerr << "out=none|xyz|mat (mat)" << endl;

}


int main(int argc, char *argv[]) {

	if( argc < 2 ) {
		cerr << "usage : " << argv[0] << " SacFile [options]" << endl;
		usage();
		return 1;
	}

	LANDES::ReadSac onesac( argv[1]);
	if ( !onesac.isOk() ) {
		cerr << "ERREUR : impossible d'ouvrir le fichier " << argv[1] << endl;
		return 1;
	}
	ConfigParser parseArg;
	for(int i=2; i<argc; i++) {
		parseArg.parse(argv[i]);
	}

	LANDES::FTA::Parameter p;
	p.init(parseArg);

//	if( p.verbose )
	p.print();


	LANDES::FTA::Signal sac_sig;
	sac_sig.initFromSac(onesac, p.vgMin, p.vgMax);


//	cout << "Read Arguments:" << endl;
//	parseArg.write(cout);



	vector<double> freqs = LANDES::FTA::powspace(p.fmin, p.fmax, p.nfreqs);
	vector<double> b = LANDES::FTA::linspace(p.bmin, p.bmax, p.nfreqs);

	LANDES::FTA::TimeFreq_Diagram ftaDiag(sac_sig);

	for(size_t i=0; i<freqs.size(); ++i) {
		ftaDiag.addFreq(freqs[i], b[i]);
	}

	ftaDiag.normalizeMax();

	cout << "max: " << freqs[ftaDiag.findFreq_forMax()] << endl;
	if( p.finit > 0.0)
		cout << "freq near " << p.finit << " = " << freqs[ftaDiag.findFreq_near(p.finit)] << endl;

	//ftaDiag.normalizePerFreq();

	LANDES::FTA::DispersionCurve disp;

	if( p.dispOpt != LANDES::FTA::NO_DISP ) {
		vector<int> ifreq, itg;

		int iTginit=-1;
		switch(p.dispOpt) {
			case LANDES::FTA::MAX_DISP:
				ftaDiag.calcMaxDisp(disp);
				break;
			case LANDES::FTA::ALL_DISP:
				ftaDiag.calcAllPtDisp(disp);
				break;
			case LANDES::FTA::CONT_DISP:
				size_t iFinit;
				if( p.finit > 0.0)
					iFinit = ftaDiag.findFreq_near(p.finit);
				else
					iFinit = ftaDiag.findFreq_forMax();

				if( p.vginit > 0.0 ) 
					p.tginit = sac_sig.vToTime(p.vginit);

				if( p.tginit > 0.0)
					iTginit = sac_sig.timeToIndice(p.tginit);

				ftaDiag.calcContinuousDisp(iFinit, disp, p.ampMin, iTginit);
				break;
			case LANDES::FTA::CUSTOM_PT:
				for(size_t i=2; i<=p.customPt.size(); i+=2) {
					double val = p.customPt[i-2];
					if( !p.fta_freq ) val = 1.0/val;
					
					ifreq.push_back(ftaDiag.findFreq_near(val));
					//cout << ftaDiag.findFreq_near(val);
					val = p.customPt[i-1];
					if( !p.fta_time ) val = sac_sig.vToTime(val);
					itg.push_back(sac_sig.timeToIndice(val));
					//cout << " " << sac_sig.timeToIndice(val) << endl;

				}
				ftaDiag.calclCustomDisp(disp, ifreq, itg);
				break;
			case LANDES::FTA::NO_DISP:
				break;

			default:
				break;
		}
		if( !p.fta_time) convertDispTimeToVg(disp, sac_sig);

		writeDisp("write", disp, p.fta_freq);
	}

	if( p.outputFta ) {
		if( p.outputFta_matrixFormat )
			ftaDiag.writeMat("write", p.fta_freq, p.fta_time);
		else
			ftaDiag.writeXYZ("write", p.fta_freq, p.fta_time);

	}

	return 0;
}

