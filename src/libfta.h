/*********************************************
 * auteur : Matthieu LANDES
 * mail   : landes@ipgp.fr
 * date :   12/09/2013
 * Institut de Physique du Globe de Paris
 * Equipe de sismologie
 ********************************************/

#ifndef LANDES_LIBFTA_HEADER
#define LANDES_LIBFTA_HEADER

#include <vector>

namespace LANDES {

class ReadSac;

namespace FTA {

struct Signal {
	std::vector<double> data;
	size_t log_2N, npts;
	double dt, t0, bt, dist, iMin, iMax;


	void initFromSac(const LANDES::ReadSac& sacheader, double vgMin, double vgMax);

	size_t vToIndice(double v) const;
	size_t timeToIndice(double t) const;

	double indiceToV(size_t indice) const;
	double indiceToTime(size_t indice) const;
	double timeToV(double t) const;
	double vToTime(double v) const;

};

struct OneMeasure {
	double freq, grp_mes, snr, amp;
	OneMeasure(double f, double t) : freq(f), grp_mes(t), snr(1.0),  amp(1.0) {}
	OneMeasure(double f, double t, double a) : freq(f), grp_mes(t), snr(1.0),  amp(a) {}
};

typedef std::vector<OneMeasure> DispersionCurve;

void convertDispTimeToVg(DispersionCurve& disp, const Signal& sig);

void writeDisp(const char * prefix, const DispersionCurve& disp, bool fpQ);

struct TimeFreq_Diagram {
	std::vector<double> amplitude, phase;
	std::vector<double> freq, time;


	TimeFreq_Diagram(Signal& sig);

	void addFreq(double freq, double beta);

	size_t findFreq_forMax() const;
	size_t findFreq_near(double f) const;
	void normalizeMax();
	void normalizePerFreq();

	void calcMaxDisp(DispersionCurve& curve);
	void calcContinuousDisp(size_t iFreqRef, DispersionCurve& curve, double ampMin, int iTgRef);
	void calcAllPtDisp(DispersionCurve& curve);
	void calclCustomDisp(DispersionCurve& curve, const std::vector<int>& ifreq, const std::vector<int>& itg);

	bool writeMat(const char* prefix, bool fpQ, bool tVgQ);
	bool writeXYZ(const char* prefix, bool fpQ, bool tVgQ);

	private:
		std::vector<double> buffer;
		Signal& m_sig;


		template<class Iterator>
		void writeOneTable(std::ostream& flux,  Iterator start, int incFreq, int incDt) {
			int imin = static_cast<int>(m_sig.iMin), imax = static_cast<int>(m_sig.iMax);
			if( incDt < 0.0 ) {
				int tmp = imin;
				imin = imax;
				imax = tmp;
			}

			Iterator its, ite;
			for(size_t i=0; i<freq.size(); ++i) {
				its = start + imin;
				ite = start + imax;
				while(its != ite) {
					flux << *its << " ";
					its += incDt;
				}
				start += incFreq;
				flux << std::endl;
			}
		}

};


typedef std::pair<size_t, double> PtInfo;
typedef std::vector<LANDES::FTA::PtInfo> VPtInfo;



template<class T>
T mySquare(T a) {
	return a*a;
}


std::vector<double> powspace(double mi, double ma, size_t n);
std::vector<double> linspace(double a, double b, size_t n);

}
}


#endif
