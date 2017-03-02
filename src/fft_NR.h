/*********************************************
 * auteur : Matthieu LANDES
 * mail   : landes@ipgp.fr
 * date :   12/09/2013
 * Institut de Physique du Globe de Paris
 * Equipe de sismologie
 ********************************************/

#ifndef FFT_NR_HEADER
#define FFT_NR_HEADER

#include <vector>

void four1(double *data, const int n, const int isign);
void four1(std::vector<double>& data, const int isign);
void realft(std::vector<double>& data, const int isign);


#endif
