/*********************************************
 * auteur : Matthieu LANDES
 * mail   : landes@ipgp.fr
 * date :   12/09/2013
 * Institut de Physique du Globe de Paris
 * Equipe de sismologie
 ********************************************/

#ifndef MATHFUNC_HEADER
#define MATHFUNC_HEADER

#include <cmath>
#include <limits>

const double pi = 3.14159265358979323846264338327950288419716939937510;

template< typename T>
struct Hamming {
	T operator() (T x) const {
		T const two = static_cast<T>(2);

		return static_cast<T>(0.54) - static_cast<T>(0.46)*std::cos(two * static_cast<T>(pi)*x);
	}
};

template< typename T>
class Blackman {
	private:
		T alpha_over2_;
	public:
		Blackman(const T a = 0.16) : alpha_over2_(a/static_cast<T>(2)) {}
		T operator() (const T x) const {
			T const one = static_cast<T>(1);
			T const two = static_cast<T>(2);
			T const demi = one/static_cast<T>(2);


			return demi - alpha_over2_ - demi * std::cos(two * pi * x) 
				+  alpha_over2_ * std::cos(two*two * pi *x);
		}
};

template<typename T>
class Tukey {
	private:
		T alpha_over2_;

	public:
		Tukey(const T a) : alpha_over2_(a/static_cast<T>(2)) {}

		T operator() (const T x) const {

			T const one = static_cast<T>(1);
			T const zero = static_cast<T>(0);
			T const demi = one/static_cast<T>(2);

			if( x < zero)
				return zero;
			else if( x < alpha_over2_)
				return demi*(one+std::cos(static_cast<T>(pi)*(x/alpha_over2_ - one)));

			else if( x < one-alpha_over2_)
				return one;
			else if( x < one)
				return demi*(one+std::cos(static_cast<T>(pi)*((x-one)/alpha_over2_ + one)));
			else
				return zero;

		}

};

template <typename T>
T sinc(const T x) {

	static T const    taylor_0_bound = std::numeric_limits<T>::epsilon();
	static T const    taylor_2_bound = std::sqrt(taylor_0_bound);
	static T const    taylor_n_bound = std::sqrt(taylor_2_bound);


	if (std::abs(x) >= taylor_n_bound) {
		return(std::sin(x)/x);
	}
	else {
		// approximation by taylor series in x at 0 up to order 0
		T    result = static_cast<T>(1);
		if    (std::abs(x) >= taylor_0_bound) {
			T    x2 = x*x;
			// approximation by taylor series in x at 0 up to order 2
			result -= x2/static_cast<T>(6);
			if    (std::abs(x) >= taylor_2_bound) {
				// approximation by taylor series in x at 0 up to order 4
				result += (x2*x2)/static_cast<T>(120);
			}
		}
		return(result);
	}

}



template <typename T>
T sinc_pi ( const T x) {
	return sinc<T>(x * static_cast<T>(pi));
}



#endif
