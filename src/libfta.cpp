/*********************************************
 * auteur : Matthieu LANDES
 * mail   : landes@ipgp.fr
 * date :   12/09/2013
 * Institut de Physique du Globe de Paris
 * Equipe de sismologie
 ********************************************/


#include <iostream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <limits>
#include <functional>


#include "libfta.h"
#include "readsac.h"
#include "fft_NR.h"
#include "mathfunc.h"

using namespace std;

const double Pi = 3.1415926535897;
const double Pi2 = 2*Pi;



double filtre_Gaussien(double f, double fh, double beta) 
{
  //printf("freq : %f filtrage : %f beta : %f\n",f, fh, beta);
  return 1.0/beta*exp(-0.5*(f-fh)*(f-fh)/(beta*beta));
}





std::vector<double> LANDES::FTA::powspace(double a, double b, size_t n) {
	vector<double> res(n);
	for(size_t i=0; i<n; ++i)
		res[i] = a*pow(b/a, static_cast<double>(i)/(n-1.));

	return res;
}


std::vector<double> LANDES::FTA::linspace(double a, double b, size_t n) {
	vector<double> res(n);
	for(size_t i=0; i<n; ++i)
		res[i] = a + (b-a)*static_cast<double>(i)/(n-1.);

	return res;
}

template< class Iterator>
void unwrap(Iterator inCurrent, Iterator inEnd, Iterator outCurent) {
	typename Iterator::value_type cp, dp, pm1 = *inCurrent++, cadran(0),
		seuil(Pi-std::numeric_limits<typename Iterator::value_type>::epsilon());
	*outCurent++ = pm1;

	while( inCurrent!= inEnd) {
		cp = *inCurrent + cadran;
		dp = cp - pm1;
		pm1 = cp;
		if( dp>seuil) {
			while( dp>seuil) {
				cadran -= Pi2;
				dp -= Pi2;
			}
		}
		if( dp<-seuil) {
			while(dp<-seuil) {
				cadran += Pi2;
				dp += Pi2;
			}
		}
		cp = *inCurrent++ + cadran;
		pm1 = cp;
		*outCurent++ = pm1;
	}
}


template<class Iterator>
void normalizeAray(Iterator beg, Iterator end) {
	typename Iterator::value_type m = *max_element(beg, end);
	if( m != 0)
		m = static_cast<typename Iterator::value_type>(1)/m;

	while(beg != end) {
		*beg++ *= m;
	}
}

template<class Iterator>
void findAllMaxArray(Iterator beg, size_t n, vector<LANDES::FTA::PtInfo>& listMax) {
	listMax.clear();
	
	typename Iterator::value_type d1, d2;
	beg++;
	for(size_t i=1; i<n; ++i, ++beg) {
		d1 = *beg - *(beg-1);
		d2 = *(beg+1) - *beg;

		if( (d1 >= 0.0) && (d2 <= 0.0) ) {
			listMax.push_back( make_pair(i,*beg) );
		}
	}
}


/*
double diff(const vector<double>& v, double dx, size_t icenter) {
	if( icenter == 0) 
		return 1.0/Pi2*(v[icenter+1]-v[icenter])/dx;
	else if( icenter+1 == v.size())
		return 1.0/Pi2*(v[icenter]-v[icenter-1])/dx;
	else
		return 0.5/Pi2*(v[icenter+1]-v[icenter-1])/dx;
}*/

template<class Iterator>
double diff(Iterator beg, Iterator end, double dx, Iterator center) {
	if( center == beg )
		return 1.0/Pi2*( *(center+1)-*center)/dx;
	else if( center+1 == end) 
		return 1.0/Pi2*( *center - *(center-1) )/dx;
	else
		return 0.5/Pi2*( *(center+1) - *(center-1))/dx;

}


void LANDES::FTA::Signal::initFromSac(const LANDES::ReadSac& sacheader, double vgMin, double vgMax) {
	data.assign(sacheader.trace.begin(), sacheader.trace.end());
	npts = data.size();
	dt = sacheader.header.delta;
	bt = sacheader.header.b;
	dist = sacheader.header.dist;

	iMin=0.0, iMax=npts;
	if( vgMin > 0.0 ) {
		double tmax = dist/vgMin;
		iMax = min( timeToIndice(tmax), npts);
	}

	if( vgMax > 0.0) {
		double tmin = max(dist/vgMax, bt);
		 iMin = timeToIndice(tmin);
	}

	Tukey<double> w(0.01);
//	Hamming<double> w;
	for(size_t i=0; i<npts; ++i)
		data[i] *= w( (static_cast<double>(i)-iMin)/(iMax-iMin));



	log_2N = static_cast<size_t>(ceil(log(npts)/log(2.0)));

	data.resize( 1 << log_2N, 0.0);

	vector<double> save = data;

	realft(data,1);

	data.resize(data.size()*2,0.0);

}

double LANDES::FTA::Signal::indiceToV(size_t i) const {
	return dist/indiceToTime(i);
}

double LANDES::FTA::Signal::indiceToTime(size_t i) const {
	return bt + static_cast<double>(i) * dt;
}

double  LANDES::FTA::Signal::timeToV(double t) const {
	return dist/t;
}

double  LANDES::FTA::Signal::vToTime(double v) const {
	return dist/v;
}

size_t LANDES::FTA::Signal::vToIndice(double v) const {
	return timeToIndice(dist/v);
}

size_t LANDES::FTA::Signal::timeToIndice(double t) const {
	double it = max((t - bt)/dt,0.0);
	return min(static_cast<size_t>(it), npts);
}


LANDES::FTA::TimeFreq_Diagram::TimeFreq_Diagram(Signal& sig) : 
	m_sig(sig) 
{

	time.reserve(m_sig.npts);

	for(size_t i=0; i<m_sig.npts; ++i) {
		time.push_back(m_sig.bt + static_cast<double>(i) * m_sig.dt);
	}
	buffer.resize(m_sig.data.size());
	
}


void LANDES::FTA::TimeFreq_Diagram::addFreq(double f, double beta) {
	fill(buffer.begin(), buffer.end(), 0.0);
	double w = 1.0/(static_cast<double>(m_sig.data.size())*m_sig.dt);

	for(size_t i=0; i<m_sig.data.size(); i+=2) {
		double t = filtre_Gaussien( w*static_cast<double>(i), f, beta);
		buffer[i] = m_sig.data[i] * t;
		buffer[i+1] = m_sig.data[i+1] * t;
	}

	buffer[0] /= 2.0;
	buffer[1] /= 2.0;

	four1(buffer, -1);

	double a,p;


	for(size_t i=0; i<m_sig.npts; i++) {
		a = sqrt(LANDES::FTA::mySquare(buffer[i*2]) + LANDES::FTA::mySquare(buffer[2*i+1]));

		if( a == 0.0) {
			a = 0.0;
			p = 0.0;
		}
		else {
			p = atan2(buffer[2*i+1], buffer[2*i]);
		}

		amplitude.push_back(a);
		phase.push_back(p);

	}
	freq.push_back(f);

	unwrap(phase.end() - m_sig.npts, phase.end(), phase.end() - m_sig.npts);

}


size_t LANDES::FTA::TimeFreq_Diagram::findFreq_forMax() const {

	size_t d = distance(amplitude.begin(), 
			max_element(amplitude.begin(), amplitude.end()));

	return d/m_sig.npts;
}

size_t LANDES::FTA::TimeFreq_Diagram::findFreq_near(double f) const {
	return distance( freq.begin(), lower_bound(freq.begin(), freq.end()-1, f));
}

void LANDES::FTA::TimeFreq_Diagram::normalizeMax() {
	normalizeAray(amplitude.begin(), amplitude.end());
}

void LANDES::FTA::TimeFreq_Diagram::normalizePerFreq() {
	vector<double>::iterator it = amplitude.begin();
	for(size_t i=0; i<freq.size(); ++i) {
		normalizeAray(it, it + m_sig.npts);
		it += m_sig.npts;
	}
}

void LANDES::FTA::TimeFreq_Diagram::calcMaxDisp(DispersionCurve& curve) {
	curve.clear();

	vector<double>::iterator it = amplitude.begin(), itmax, itphase = phase.begin();
	for(size_t i=0; i<freq.size(); ++i) {
		itmax = max_element(it, it+m_sig.npts);
		size_t d = distance(it, itmax);

		double f = -diff(itphase, itphase+m_sig.npts, m_sig.dt, itphase+d);
		curve.push_back(OneMeasure(f, m_sig.indiceToTime(d), *itmax));
		it += m_sig.npts;
		itphase += m_sig.npts;
	}
}

bool comp_amp_PtInfo( LANDES::FTA::PtInfo a, LANDES::FTA::PtInfo b) {
	return a.second < b.second;
}


void LANDES::FTA::TimeFreq_Diagram::calcAllPtDisp(DispersionCurve& curve) {

	LANDES::FTA::VPtInfo listMax;
	vector<double>::iterator itphase = phase.begin(), itamplitude = amplitude.begin();
	double f;

	curve.clear();

	for(size_t i=0; i<freq.size(); ++i, itamplitude += m_sig.npts, itphase += m_sig.npts) {
		findAllMaxArray(itamplitude, m_sig.npts, listMax);
		LANDES::FTA::VPtInfo::iterator it = listMax.begin();
		for(; it != listMax.end(); ++it) {
			f = -diff(itphase, itphase+m_sig.npts, m_sig.dt, itphase+it->first);
			curve.push_back(OneMeasure(f, m_sig.indiceToTime(it->first), it->second));
		}

	}
}

void LANDES::FTA::TimeFreq_Diagram::calclCustomDisp(DispersionCurve& curve, const std::vector<int>& ifreq, const std::vector<int>& itg) {
	vector<double>::iterator itphase = phase.begin(), itamplitude = amplitude.begin();
	double f;
	size_t n = min(ifreq.size(), itg.size());

	for(size_t i=0; i<n; ++i) {
		int deltaF = ifreq[i]*m_sig.npts;
		f = -diff(itphase+deltaF, itphase+deltaF+m_sig.npts, m_sig.dt, itphase+deltaF+itg[i]);
		curve.push_back(OneMeasure(f, m_sig.indiceToTime(itg[i]), *(itamplitude+deltaF+itg[i])));
	}


}

struct DistPtInfo {
	private:
		double i,a;

	public:
		DistPtInfo(LANDES::FTA::PtInfo ptRef) : i(static_cast<double>(ptRef.first)), a(ptRef.second) {}

		double operator() (LANDES::FTA::PtInfo pt) {
			return abs(static_cast<double>(pt.first) - i);
		}


};

struct NearestPt {
	double indiceRef;

	NearestPt(int i) : indiceRef(static_cast<double>(i)) {}

	double operator()(LANDES::FTA::PtInfo pt) {
		return abs(static_cast<double>(pt.first) - indiceRef);
	}

};




void LANDES::FTA::TimeFreq_Diagram::calcContinuousDisp(size_t iFreqRef, DispersionCurve& curve, double ampMin, int iTgRef=-1) {
	LANDES::FTA::VPtInfo listMax;
	LANDES::FTA::PtInfo firstPt, refPt;
	vector<double>::iterator itphase, itamplitude;
	curve.clear();

	itamplitude = amplitude.begin() + iFreqRef*m_sig.npts;
	itphase = phase.begin() + iFreqRef*m_sig.npts;

	findAllMaxArray(itamplitude, m_sig.npts, listMax);
	LANDES::FTA::VPtInfo::iterator itmax;
	if( iTgRef < 0 ) {
		itmax = max_element(listMax.begin(), listMax.end(), comp_amp_PtInfo);
		refPt = *itmax;
	}
	else {
		//vector<double> l;
		//l.resize(listMax.size());
		//transform(listMax.begin(), listMax.end(), l.begin(), NearestPt(iTgRef));
		//refPt = listMax[distance(l.begin(), min_element(l.begin(), l.end()))];

		refPt.first = iTgRef;
		refPt.second = *(itamplitude + iTgRef);
	}


	firstPt = refPt;

	double f = -diff(itphase, itphase+m_sig.npts, m_sig.dt, itphase+refPt.first);
	curve.push_back(OneMeasure(f, m_sig.indiceToTime(refPt.first), refPt.second));

	vector<double> dist;
	for(long int i=iFreqRef+1; i<static_cast<long int>(freq.size()); ++i) {
		itamplitude = amplitude.begin() + i*m_sig.npts;
		itphase = phase.begin() + i*m_sig.npts;

		findAllMaxArray(itamplitude, m_sig.npts, listMax);
		dist.resize(listMax.size());
		transform(listMax.begin(), listMax.end(), dist.begin(), DistPtInfo(refPt));

		size_t d = distance( dist.begin(), min_element(dist.begin(), dist.end()));
		refPt = *(listMax.begin()+d);

		if( refPt.second <= ampMin ) break;
		f = -diff(itphase, itphase+m_sig.npts, m_sig.dt, itphase+refPt.first);
//		curve.push_back(OneMeasure(freq[i], m_sig.indiceToTime(refPt.first), refPt.second));
		curve.push_back(OneMeasure(f, m_sig.indiceToTime(refPt.first), refPt.second));
	}

	refPt = firstPt;
	for(long int i=iFreqRef-1; i>0; --i) {
		itamplitude = amplitude.begin() + i*m_sig.npts;
		itphase = phase.begin() + i*m_sig.npts;

		findAllMaxArray(itamplitude, m_sig.npts, listMax);
		dist.resize(listMax.size());
		transform(listMax.begin(), listMax.end(), dist.begin(), DistPtInfo(refPt));

		size_t d = distance( dist.begin(), min_element(dist.begin(), dist.end()));
		refPt = *(listMax.begin()+d);

		if( refPt.second <= ampMin ) break;
		f = -diff(itphase, itphase+m_sig.npts, m_sig.dt, itphase+refPt.first);
//		curve.push_back(OneMeasure(freq[i], m_sig.indiceToTime(refPt.first), refPt.second));
		curve.push_back(OneMeasure(f, m_sig.indiceToTime(refPt.first), refPt.second));
	}



}

bool LANDES::FTA::TimeFreq_Diagram::writeXYZ(const char* prefix, bool fpQ, bool tVgQ) {


	vector<double>::iterator itwrite;

	string filename(prefix);
	filename += "_amp.txt";
	ofstream flux(filename.c_str());
	if( flux.bad() )
		return false;

	double ft, tv;
	itwrite = amplitude.begin();
	for(size_t i=0; i<freq.size(); ++i, itwrite += freq.size()) {
		ft = fpQ?freq[i]:1.0/freq[i];
		for(size_t itime=m_sig.iMin; itime<m_sig.iMax; ++itime) {
			tv = tVgQ?time[itime]:m_sig.timeToV(time[itime]);
			flux << ft << " " << tv << " " << *(itwrite+itime) << endl;
		}
	}
	flux.close();


	return true;
}

bool LANDES::FTA::TimeFreq_Diagram::writeMat(const char* prefix, bool fpQ, bool tVgQ) {

	int difreq;
	vector<double>::iterator ita, itp;
	int imin = static_cast<int>(m_sig.iMin), imax = static_cast<int>(m_sig.iMax);
	if( fpQ ) {
		ita =  amplitude.begin();
		itp = phase.begin();
		difreq = time.size();
	}
	else {
		itp = phase.end() - time.size();
		ita = amplitude.end() - time.size();
		difreq = -time.size();
	}

	int dit = (tVgQ?1:-1);
	string filename(prefix);
	filename += "_amp.txt";
	ofstream flux(filename.c_str());
	if( flux.bad() )
		return false;

	cout << "write amplitude" << endl;
	writeOneTable(flux, ita, difreq, dit);
	flux.close();

	filename = prefix;
	filename += "_ph.txt";
	flux.open(filename.c_str());
	if( flux.bad() )
		return false;

	cout << "write phase" << endl;
	writeOneTable(flux, itp, difreq, dit);
	flux.close();


	cout << "write FT" << endl;

	filename = prefix;
	filename += "_FP.txt";
	flux.open(filename.c_str());
	if( flux.bad() )
		return false;
	if( fpQ ) {
		copy(freq.begin(), freq.end(),  ostream_iterator<double>(flux, " "));
	}
	else {
		vector<double>::reverse_iterator itout = freq.rbegin();
		for(; itout != freq.rend(); ++itout) flux << 1.0/ *itout << " ";
	}
	flux.close();


	cout << "write TV" << endl;
	filename = prefix;
	filename += "_TV.txt";
	flux.open(filename.c_str());
	if( flux.bad() )
		return false;
	if( tVgQ ) {
		copy(time.begin()+imin, time.begin()+imax,  ostream_iterator<double>(flux, " "));
	}
	else {
		vector<double>::reverse_iterator itout = time.rbegin()+(m_sig.npts-imax);
		for(; itout != time.rend()-imin; ++itout) flux << m_sig.timeToV(*itout) << " ";
	}
	flux.close();


	return true;

}



void LANDES::FTA::writeDisp(const char * prefix, const DispersionCurve& disp, bool fpQ) {

	string filename(prefix);
	filename += "_disp.txt";

	ofstream flux(filename.c_str());

	for(size_t i=0; i<disp.size(); ++i) {
		if( fpQ) {
			flux << disp[i].freq;
		}
		else 
			flux << 1.0/disp[i].freq;

		flux << " " << disp[i].grp_mes << " " << disp[i].amp << endl;
	}
	flux.close();

}

void LANDES::FTA::convertDispTimeToVg(DispersionCurve& disp, const Signal& sig) {
	DispersionCurve::iterator it = disp.begin();
	for(;it != disp.end(); ++it) {
		it->grp_mes = sig.timeToV(it->grp_mes);
	}

}

