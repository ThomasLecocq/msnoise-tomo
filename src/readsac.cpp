/*********************************************
 * auteur : Matthieu LANDES
 * mail   : landes@ipgp.fr
 * date :   12/09/2013
 * Institut de Physique du Globe de Paris
 * Equipe de sismologie
 ********************************************/

#include "readsac.h"

#include <fstream>
#include <iostream>

using namespace std;
using namespace LANDES;



LANDES::ReadSac::ReadSac(const char* filename): m_isReadOk(false) {

	ifstream flux(filename, ios::binary);

	if( flux.bad() ) {
		cerr << "Error: unable to open [" << filename << "] !" << endl;
		return;
	}


	//read header
	flux.read(reinterpret_cast<char*>(&header), sizeof(SacHeader));
	if( !flux ) {
		cerr << "Error: read " << flux.gcount() <<
		   	" of expected header size(" << sizeof(SacHeader) << ")" << endl;
		return;
	}

	//read trace
	int32_t npts = header.npts;
	trace.resize(npts);
	flux.read(reinterpret_cast<char*>(trace.data()), npts*sizeof(float));
	if( !flux ) {
		cerr << "Error: read " << flux.gcount() <<
		   	" of expected trace size(" << npts*sizeof(float) << ")" << endl;
		return;
	}


	m_isReadOk = true;

}
