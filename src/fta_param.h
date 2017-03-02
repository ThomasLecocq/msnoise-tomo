/*********************************************
 * auteur : Matthieu LANDES
 * mail   : landes@ipgp.fr
 * date :   12/09/2013
 * Institut de Physique du Globe de Paris
 * Equipe de sismologie
 ********************************************/

#ifndef FTA_PARAM_HEADER
#define FTA_PARAM_HEADER

class ConfigParser;

#include <vector>

namespace LANDES {
namespace FTA {

enum DispOptions {
	NO_DISP=0,
	MAX_DISP=1,
	CONT_DISP=2,
	ALL_DISP=3,
	CUSTOM_PT=4
};

struct Parameter {
	double fmin, fmax, bmin, bmax, finit, tginit, ampMin;
	double vgMin, vgMax, vginit;
	size_t nfreqs;

	std::vector<double> customPt;

	bool fta_freq, fta_time;
	DispOptions dispOpt;
	bool outputFta, outputFta_matrixFormat;
	bool verbose;
	bool askFmin, askFmax;

	void init(const ConfigParser& config);
	void print() const;
};



}
}

#endif
