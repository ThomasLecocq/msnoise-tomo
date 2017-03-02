/*********************************************
 * auteur : Matthieu LANDES
 * mail   : landes@ipgp.fr
 * date :   12/09/2013
 * Institut de Physique du Globe de Paris
 * Equipe de sismologie
 ********************************************/

#ifndef READSAC_HEADER
#define READSAC_HEADER

#include <vector>
//#include <cstdint>
#include <stdint.h>

namespace LANDES {

	struct SacHeader {
		float	delta,     depmin,    depmax,    scale,     odelta;    
		float	b,         e,         o,         a,         internal1; 
		float	t0,        t1,        t2,        t3,        t4;        
		float	t5,        t6,        t7,        t8,        t9;        
		float	f,         resp0,     resp1,     resp2,     resp3;     
		float	resp4,     resp5,     resp6,     resp7,     resp8;     
		float	resp9,     stla,      stlo,      stel,      stdp;      
		float	evla,      evlo,      evel,      evdp,      unused1;   
		float	user0,     user1,     user2,     user3,     user4;     
		float	user5,     user6,     user7,     user8,     user9;     
		float	dist,      az,        baz,       gcarc,     internal2; 
		float	internal3, depmen,    cmpaz,     cmpinc,    unused2;   
		float	unused3,   unused4,   unused5,   unused6,   unused7;   
		float	unused8,   unused9,   unused10,  unused11,  unused12;  
		int32_t	nzyear,    nzjday,    nzhour,    nzmin,     nzsec;     
		int32_t	nzmsec,    internal4, internal5, internal6, npts;      
		int32_t	internal7, internal8, unused13,  unused14,  unused15;  
		int32_t	iftype,    idep,      iztype,    unused16,  iinst;     
		int32_t	istreg,    ievreg,    ievtyp,    iqual,     isynth;    
		int32_t	unused17,  unused18,  unused19,  unused20,  unused21;  
		int32_t	unused22,  unused23,  unused24,  unused25,  unused26;  
		int32_t	leven,     lpspol,    lovrok,    lcalda,    unused27;  
		char	kstnm[8],  kevnm[16];           
		char	khole[8],  ko[8],     ka[8];               
		char	kt0[8],    kt1[8],    kt2[8];              
		char	kt3[8],    kt4[8],    kt5[8];              
		char	kt6[8],    kt7[8],    kt8[8];              
		char	kt9[8],    kf[8],     kuser0[8];           
		char	kuser1[8], kuser2[8], kcmpnm[8];           
		char	knetwk[8], kdatrd[8], kinst[8];  
	};





	class ReadSac {
		private:
			bool m_isReadOk;
		public:
			SacHeader header;
			std::vector<float> trace;

			ReadSac(const char* filename);

			bool isOk() const { return m_isReadOk; }


	};



}


#endif
