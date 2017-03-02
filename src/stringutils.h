/*********************************************
 * auteur : Matthieu LANDES
 * mail   : landes@ipgp.fr
 * date :   12/09/2013
 * Institut de Physique du Globe de Paris
 * Equipe de sismologie
 ********************************************/

#ifndef STRINGUTILS_HEADER
#define STRINGUTILS_HEADER

#include <string>
#include <sstream>
#include <vector>
#include <iostream>


template<class T>
T from_string( const std::string & Str)
{
	T Dest;
    std::istringstream iss( Str );
    iss >> Dest;
	return Dest;
}

template<typename T>
std::string to_string( const T & Value )
{
    std::ostringstream oss;
    oss << Value;
    return oss.str();
}

/*
template<>
std::string from_string(const std::string & Str) {
	return Str;
}
*/

template<class T, class IteratorOut>
IteratorOut itemize(const std::string& txt, IteratorOut itemIt, const char delim = ' ') {
//	typedef typename IteratorOut::value_type T;
	std::string dummy;

	std::stringstream istr( txt);
//	std::cout << "Read : ";
	while(! istr.eof()) {
		getline(istr, dummy, delim);
//		std::cout << '[' << dummy << ']';
		if( !dummy.empty()) {
//			std::cout << '*';
			*itemIt++ = from_string<T>(dummy);
		}
//		std::cout << ' ';

	}

//	std::cout << std::endl;
	return itemIt;

}





#endif
