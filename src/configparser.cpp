/*********************************************
 * auteur : Matthieu LANDES
 * mail   : landes@ipgp.fr
 * date :   12/09/2013
 * Institut de Physique du Globe de Paris
 * Equipe de sismologie
 ********************************************/



#include "configparser.h"


using namespace std;

string whitespaces(" \t\f\v\n\r");

string rlstrip(const string& txt) {
    size_t first = txt.find_first_not_of(whitespaces);
    size_t last = txt.find_last_not_of(whitespaces);

    if( first == string::npos)
        return "";

    return txt.substr(first, last-first+1);

}


bool isAssignation(const string& txt, string& key, string& val) {
    size_t indice = txt.find_first_of('=');;
    bool res = false;

    if( indice != string::npos ) {
        res = true;
        key = rlstrip(txt.substr(0,indice));
        val = rlstrip(txt.substr(indice+1, string::npos));

        if( key.empty() || val.empty() )
            res = false;
    }

    return res;
}


bool ConfigParser::hasKey(const string& key) const {
	return m_listeEntry.find(key) != m_listeEntry.end();
}

void ConfigParser::parse(const string& txt) {
	string key, val;

	if( isAssignation(txt, key, val) ) {
		if( !key.empty())
			m_listeEntry[key] = val;
	}

}

void ConfigParser::write(ostream& flux) const {
	MapEntry::const_iterator it = m_listeEntry.begin();

	for(; it != m_listeEntry.end(); ++it) {
		flux << it->first << " = " << it->second << endl;
	}

}
