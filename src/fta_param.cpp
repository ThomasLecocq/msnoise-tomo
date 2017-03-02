/*********************************************
 * auteur : Matthieu LANDES
 * mail   : landes@ipgp.fr
 * date :   12/09/2013
 * Institut de Physique du Globe de Paris
 * Equipe de sismologie
 ********************************************/


#include <string>
#include <iterator>

#include "fta_param.h"
#include "configparser.h"
#include "stringutils.h"

using namespace std;



void LANDES::FTA::Parameter::print() const {
	cout << "fmin=" << fmin << ", fmax=" << fmax << ", nfreqs=" << nfreqs << endl;
	cout << "bmin=" << bmin << ", bmax=" << bmax << endl;

	if( ampMin > 0.0)
		cout << "ampMin=" << ampMin << endl;


	if( vgMin>=0.0)
		cout << "vgMin=" << vgMin << " ";
	if( vgMax>=0.0)
		cout << "vgMax=" << vgMax;

	cout << endl;

	cout << "disp: ";
	switch(dispOpt) {
		case MAX_DISP:
			cout << "max";
			break;
		case CONT_DISP:
			cout << "cont, ";
			if( finit > 0.0) cout << " finit=" << finit;
			else cout << "finit=max";
			cout << ' ';
			if( tginit > 0.0 )
				cout << "tginit=" << tginit;
			else if( vginit > 0.0)
				cout << "vginit=" << vginit;
			break;
		case ALL_DISP:
			cout << "all";
			break;
		case CUSTOM_PT:
			cout << "custom pt : ";
			copy(customPt.begin(), customPt.end(), ostream_iterator<double>(cout," "));
			break;
		default:
			cout << "no";
			break;
	}

	cout << endl << "diag: ";
	if( fta_freq ) cout << "Frequency";
	else cout << "Period";
	cout << " - ";
	if( fta_time ) cout << "Time";
	else cout << "Vgroup";
	cout << endl;

	cout << "out: "; 
	if( outputFta ) {
		cout << "enable";
		if( outputFta_matrixFormat ) 
			cout << " (matrix)";
		else cout << " (xyz)";

	}
	else
		cout << "disable";

	cout << endl;

}

void LANDES::FTA::Parameter::init( const ConfigParser& config) {
	askFmin = true, askFmax=true;

	outputFta = true;
	if( config.hasKey("out") ) {
		outputFta_matrixFormat = true;
		string val = config.get<string>("out");
		if( val.compare("none") == 0) {
			outputFta = false;
		}
		else if( val.compare("xyz") == 0) {
			outputFta_matrixFormat = false;
		}
	}

	dispOpt = CONT_DISP;
	if( config.hasKey("disp") ) {
		string val = config.get<string>("disp");
		if( val.compare("max") == 0) {
			dispOpt = MAX_DISP;
		}
		else if( val.compare("all") == 0) {
			dispOpt = ALL_DISP;
		}
		else if( val.compare("cont") == 0) {
			dispOpt = CONT_DISP;
		}
		else if( val.compare(0,3,"pt:") == 0) {
			dispOpt = CUSTOM_PT;
			customPt.clear();
			itemize<double>(val.substr(3), back_inserter(customPt),':');
		} 
		else
			dispOpt = NO_DISP;
	}

	ampMin = -1.0;
	if( config.hasKey("ampMin") ) {
		ampMin = config.get<double>("ampMin");
	}

	if( config.hasKey("fmin") ) {
		fmin = config.get<double>("fmin");
		askFmin = false;
	}

	if( config.hasKey("fmax") ) {
		fmax = config.get<double>("fmax");
		askFmax = false;
	}

	if( config.hasKey("tmin") ) {
		fmax = 1.0/config.get<double>("tmin");
		askFmax = false;
	}

	if( config.hasKey("tmax") ) {
		fmin = 1.0/config.get<double>("tmax");
		askFmin = false;
	}

	if( config.hasKey("tinit") ) {
		finit = 1.0/config.get<double>("tinit");
	}
	else if( config.hasKey("finit") ) {
		finit = config.get<double>("finit");
	}
	else finit = -1.0;

	if( config.hasKey("tginit") ) {
		tginit = config.get<double>("tginit");
	}
	else if(  config.hasKey("vginit") ) {
		vginit = config.get<double>("vginit");
	}
	else {
		tginit=-1.0;
		vginit=-1.0;
	}



	if( config.hasKey("nfreq") ) {
		nfreqs = config.get<size_t>("nfreq");
	}
	else
		nfreqs = 100;

	if( config.hasKey("bmin") ) {
		bmin = config.get<double>("bmin");
	}
	else bmin = 0.1;

	if( config.hasKey("bmax") ) {
		bmax = config.get<double>("bmax");
	}
	else bmax = 0.1;

	if( config.hasKey("b") ) {
		bmin = bmax = config.get<double>("b");
	}

	fta_freq = false;
	fta_time = true;
	if( config.hasKey("diag") ) {
		string val = config.get<string>("diag");
		if( val[0] == 'F' ) fta_freq = true;
		if( val[1] == 'V' ) fta_time = false;
	}

	if( config.hasKey("verbose") ) {
		verbose = (config.get<int>("verbose") == 1);
	}

	vgMin = -1.0;
	vgMax = -1.0;
	if( config.hasKey("vgMin") ) {
		vgMin = config.get<double>("vgMin");
	}
	if( config.hasKey("vgMax") ) {
		vgMax = config.get<double>("vgMax");
	}

}

