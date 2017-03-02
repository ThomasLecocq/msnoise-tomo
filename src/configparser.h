/*********************************************
 * auteur : Matthieu LANDES
 * mail   : landes@ipgp.fr
 * date :   12/09/2013
 * Institut de Physique du Globe de Paris
 * Equipe de sismologie
 ********************************************/

#ifndef CONFIG_PARSER_HEADER
#define CONFIG_PARSER_HEADER

#include <string>
#include <map>
#include <iostream>

#include "stringutils.h"

class ConfigParser {
	typedef std::map<std::string, std::string> MapEntry;
	MapEntry m_listeEntry;


	public:
		void parse(const std::string& txt);

		bool hasKey(const std::string& key) const;

		template<class T>
		void add(const std::string& key, const T& value) {
			m_listeEntry[key] = toString(value);
		}

		template<class T>
		T get(const std::string& key) const {
			MapEntry::const_iterator it = m_listeEntry.find(key);
			if( it == m_listeEntry.end()) {
				return T();
			}
			return from_string<T>(it->second);
		}

		void write(std::ostream& flux) const;


};

#endif
