//TODO: add '/n' to end of error messages, maybe add color & font stuff?

#include "poly.h"


template<typename T>
Polynomial<T>::Polynomial() {
	t = l = 1;
	c = new T[1];
}

template<typename T>
Polynomial<T>::Polynomial(uint32_t terms, T coeffs[]) : t(terms), l((uint64_t)pow(2, ceil(log2(terms)))) {
	c = new T[l];
	memcpy(c, coeffs, terms*sizeof(T));
}
//TODO: add constructor that takes auto array
template<typename T>
Polynomial<T>::Polynomial(const Polynomial<T>& p) : Polynomial(p.t, p.c) {}


template<typename T>
Polynomial<T>::~Polynomial() {
	delete[] c;
}



template<typename T>
Polynomial<T>& Polynomial<T>::pad(uint64_t len) {
	uint64_t padded_len = (uint64_t)pow(2, ceil(log2(len)));
	if(padded_len > l) {
		c = (T*)realloc(c, padded_len*sizeof(T));
		memset(c+t*sizeof(T), 0, (l-t)*sizeof(T));
		l = padded_len;
	}

	return *this;
}
template<typename T>
Polynomial<T>& Polynomial<T>::strip() {
	uint32_t highest = 0;
	for(uint32_t i=0; i<t; ++i) {
		if(!(c[i] == (T)0)) {
			highest = i;
		}
	}
	t = highest+1;
	l = (uint64_t)pow(2, ceil(log2(t)));
	c = (T*)realloc(c, l*sizeof(T));

	return *this;
}//has 0 -> T cast


template<typename T>
Polynomial<T>& Polynomial<T>::setCoeffs(uint64_t len, std::complex<double> coeffs[]) {
	pad(len);
	t = (uint32_t)len;
	//to type T coefficients
	if(std::is_arithmetic<T>::value || std::is_same<T, std::complex<float> >::value || std::is_same<T, std::complex<double> >::value || std::is_same<T, std::complex<long double> >::value) {
		for(uint32_t i=0; i<t; ++i) {
			c[i] = cast(coeffs[i]);
		}
		memset(c+t*sizeof(T), 0, (l-t)*sizeof(T));
	}
	else {
		printf("Unable to cast to correct type coefficients.");
	}
	strip();
	return *this;
}//has complex(precision) -> T cast

//TODO: second setCoeffs function for type T array, maybe make nicer with realloc/memset, not using pad and strip?



template<typename T>
std::complex<double>* Polynomial<T>::recursive_fft(uint64_t n, std::complex<double> p[], std::complex<double>* q, bool inverse) {
	if(n == 1) {
		q[0] = p[0];
	}
	else {
		//w_n is an nth root of unity
		std::complex<double> w[n];
		for(uint64_t i=0; i<n; ++i) {
			double theta = inverse ? -2*M_PI*i/n : 2*M_PI*i/n;
			w[i] = std::complex<double>(cos(theta), sin(theta));
		}

		
		//decompose p_x = p_e + x*p_o
		std::complex<double> p_e[n/2]; //even
		std::complex<double> p_o[n/2]; //odd
		//get coefficents out of p_x
		for(uint64_t i=0; i<n/2; ++i) {
			p_e[i] = p[2*i];
			p_o[i] = p[2*i+1];
		}

		//recursive call
		std::complex<double> y_e[n/2];
		recursive_fft(n/2, p_e, y_e, inverse);

		std::complex<double> y_o[n/2];
		recursive_fft(n/2, p_o, y_o, inverse);

		//combine results: y[k] = y_e[k] + w_n^+-k y_o[k]
		//				   y[k+t/2] = y_e[k] - w_n^+-k y_o[k]
		std::complex<double> y[n];
		for(uint64_t k = 0; k < n/2; ++k) {
			y[k] = y_e[k] + w[k]*y_o[k];
			y[k+n/2] = y_e[k] - w[k]*y_o[k];
		}

		memcpy(q, y, n*sizeof(std::complex<double>));
	}
	return q;
}//has precision element
template<typename T>
Polynomial<T>& Polynomial<T>::fft(bool inverse) {
	//convert coefficients to complex numbers
	std::complex<double> p[l];
	if(std::is_arithmetic<T>::value || std::is_same<T, std::complex<float> >::value || std::is_same<T, std::complex<double> >::value || std::is_same<T, std::complex<long double> >::value) {
		for(uint64_t i=0; i<l; ++i) {
			p[i] = c[i];
		}
		std::complex<double> q[l];
		recursive_fft(l, p, q, inverse);
		if(inverse) {
			//divide coefficients by t
			for(uint64_t i=0; i<l; ++i) {
				q[i] /= l;
			}
		}
		setCoeffs(l, q);
	}
	else {
		printf("Cannot run FFT on polynomial because the coefficients couldn't be casted to std::complex<double> type.");
	}

	return *this;
}//has T -> complex cast, precision element

template<typename T>
Polynomial<T>& Polynomial<T>::dft() {
	return fft(false); //use counterclockwise roots of unity
}

template<typename T>
Polynomial<T>& Polynomial<T>::idft() {
	return fft(true); //use clockwise roots of unity & divide by l at end
}




template<typename T>
T Polynomial<T>::operator()(T x) {
	//Horner's method
	T y = c[t];
	for(uint32_t i=t; i>0; --i) {
		y *= x;
		y = y+c[i-1];
	}
	return y;
}

template<typename T>
T Polynomial<T>::operator[](uint32_t i) {
	return c[i];
}



template<typename T>
Polynomial<T>& Polynomial<T>::operator+=(const Polynomial<T>& q) {
	if(q.t > t) {
		if((uint64_t)q.t > l) {
			pad((uint64_t)q.t);
		}
		t = q.t;
	}
	for(uint32_t i=0; i<q.t; ++i) {
		c[i] += q.c[i];
	}
	return *this;
}
template<typename T>
Polynomial<T>& Polynomial<T>::operator+=(T a) {
	c[0] += a;
	return *this;
}

template<typename T>
Polynomial<T>& Polynomial<T>::operator-=(const Polynomial<T>& q) {
	if(q.t > t) {
		if((uint64_t)q.t > l) {
			pad((uint64_t)q.t);
		}
		t = q.t;
	}
	for(uint32_t i=0; i<q.t; ++i) {
		c[i] -= q.c[i];
	}
	return *this;
}
template<typename T>
Polynomial<T>& Polynomial<T>::operator-=(T a) {
	c[0] -= a;
	return *this;
}

template<typename T>
Polynomial<T>& Polynomial<T>::operator*=(const Polynomial<T>& q) {
	//pad input polynomials so we can get enough info to compute pq
	//make coefficients complex
	uint16_t deg = (t-1)+(q.t-1);
	uint64_t len = (uint64_t)pow(2, ceil(log2(deg+1)));
	std::complex<double> p_padded[len];
	std::complex<double> q_padded[len];
	if(std::is_arithmetic<T>::value || std::is_same<T, std::complex<float> >::value || std::is_same<T, std::complex<double> >::value || std::is_same<T, std::complex<long double> >::value) {
		for(uint32_t i=0; i<t; ++i) {
			p_padded[i] = c[i];
		}
		for(uint32_t i=0; i<q.t; ++i) {
			q_padded[i] = q.c[i];
		}

		//DFT
		std::complex<double> p_hat[len];
		std::complex<double> q_hat[len];
		recursive_fft(len, p_padded, p_hat);
		recursive_fft(len, q_padded, q_hat);

		//multiply
		std::complex<double> y_hat[len];
		for(uint64_t i=0; i<len; ++i) {
			y_hat[i] = p_hat[i]*q_hat[i];
		}

		//Inverse-DFT
		std::complex<double> y[len];
		recursive_fft(len, y_hat, y, true);
		for(uint64_t i=0; i<len; ++i) {
			y[i] /= len;
		}

		setCoeffs(len, y);
	}
	else {
		printf("Cannot multiply polynomials because the coefficients couldn't be casted to std::complex<double> type.");
	}

	return *this;
}//has precision
template<typename T>
Polynomial<T>& Polynomial<T>::operator*=(T alpha) {
	for(uint32_t i=0; i<t; ++i) {
		c[i] *= alpha;
	}
	return *this;
}

template<typename T>
Polynomial<T>& Polynomial<T>::operator/=(const Polynomial<T>& q) {
	//TODO: implement
	return *this;
}//has precision
template<typename T>
Polynomial<T>& Polynomial<T>::operator/=(T alpha) {
	for(uint32_t i=0; i<t; ++i) {
		c[i] /= alpha;
	}
	return *this;
}

template<typename T>
Polynomial<T>& Polynomial<T>::operator%=(const Polynomial<T>& q) {
	//TODO: implement
	return *this;
}
template<typename T>
Polynomial<T>& Polynomial<T>::operator%=(T alpha) {
	if(std::is_integral<T>::value) {
		for(uint32_t i=0; i<t; ++i) {
			c[i] %= alpha;
		}
	}
	else {
		printf("Unable to compute modulus on coefficients;");
	}
	
	return *this;
}
template<>
Polynomial<float>& Polynomial<float>::operator%=(float alpha) {
	for(uint32_t i=0; i<t; ++i) {
		float cq = trunc(c[i]/alpha);
		c[i] -= cq*alpha;
	}
	
	return *this;
}
template<>
Polynomial<double>& Polynomial<double>::operator%=(double alpha) {
	for(uint32_t i=0; i<t; ++i) {
		double cq = trunc(c[i]/alpha);
		c[i] -= cq*alpha;
	}
	
	return *this;
}
template<>
Polynomial<long double>& Polynomial<long double>::operator%=(long double alpha) {
	for(uint32_t i=0; i<t; ++i) {
		long double cq = trunc(c[i]/alpha);
		c[i] -= cq*alpha;
	}
	
	return *this;
}
template<>
Polynomial<std::complex<float> >& Polynomial<std::complex<float> >::operator%=(std::complex<float> alpha) {
	for(uint32_t i=0; i<t; ++i) {
		std::complex<float> cq = std::complex<float>(trunc(real(c[i]/alpha)), trunc(imag(c[i]/alpha)));
		c[i] -= cq*alpha;
	}
	
	return *this;
}
template<>
Polynomial<std::complex<double> >& Polynomial<std::complex<double> >::operator%=(std::complex<double> alpha) {
	for(uint32_t i=0; i<t; ++i) {
		std::complex<double> cq = std::complex<double>(trunc(real(c[i]/alpha)), trunc(imag(c[i]/alpha)));
		c[i] -= cq*alpha;
	}
	
	return *this;
}
template<>
Polynomial<std::complex<long double> >& Polynomial<std::complex<long double> >::operator%=(std::complex<long double> alpha) {
	for(uint32_t i=0; i<t; ++i) {
		std::complex<long double> cq = std::complex<long double>(trunc(real(c[i]/alpha)), trunc(imag(c[i]/alpha)));
		c[i] -= cq*alpha;
	}
	
	return *this;
}









//has precision element
//cast
template<>
bool Polynomial<bool>::cast(std::complex<double> z) {
	return (bool)((int)round(real(z))%2);
}
template<>
char Polynomial<char>::cast(std::complex<double> z) {
	return (char)round(real(z));
}
template<>
char16_t Polynomial<char16_t>::cast(std::complex<double> z) {
	return (char16_t)round(real(z));
}
template<>
char32_t Polynomial<char32_t>::cast(std::complex<double> z) {
	return (char32_t)round(real(z));
}
template<>
wchar_t Polynomial<wchar_t>::cast(std::complex<double> z) {
	return (wchar_t)round(real(z));
}
template<>
signed char Polynomial<signed char>::cast(std::complex<double> z) {
	return (signed char)round(real(z));
}
template<>
uint8_t Polynomial<uint8_t>::cast(std::complex<double> z) {
	return (uint8_t)abs(round(real(z)));
}
template<>
uint16_t Polynomial<uint16_t>::cast(std::complex<double> z) {
	return (uint16_t)abs(round(real(z)));
}
template<>
uint32_t Polynomial<uint32_t>::cast(std::complex<double> z) {
	return (uint32_t)abs(round(real(z)));
}
template<>
uint64_t Polynomial<uint64_t>::cast(std::complex<double> z) {
	return (uint64_t)abs(round(real(z)));
}
template<>
short int Polynomial<short int>::cast(std::complex<double> z) {
	return (short int)round(real(z));
}
template<>
int Polynomial<int>::cast(std::complex<double> z) {
	return (int)round(real(z));
}
template<>
long int Polynomial<long int>::cast(std::complex<double> z) {
	return (long int)round(real(z));
}
template<>
long long int Polynomial<long long int>::cast(std::complex<double> z) {
	return (long long int)round(real(z));
}
template<>
float Polynomial<float>::cast(std::complex<double> z) {
	return (float)real(z);
}
template<>
double Polynomial<double>::cast(std::complex<double> z) {
	return (double)real(z);
}
template<>
long double Polynomial<long double>::cast(std::complex<double> z) {
	return (long double)real(z);
}
template<>
std::complex<float> Polynomial<std::complex<float> >::cast(std::complex<double> z) {
	return std::complex<float>(real(z), imag(z));
}
template<>
std::complex<double> Polynomial<std::complex<double> >::cast(std::complex<double> z) {
	return z;
}
template<>
std::complex<long double> Polynomial<std::complex<long double> >::cast(std::complex<double> z) {
	return std::complex<long double>(real(z), imag(z));
}





//print
template<>
Polynomial<bool>& Polynomial<bool>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<char>& Polynomial<char>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<char16_t>& Polynomial<char16_t>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<char32_t>& Polynomial<char32_t>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<wchar_t>& Polynomial<wchar_t>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<signed char>& Polynomial<signed char>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %d\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %d\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<uint8_t>& Polynomial<uint8_t>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<uint16_t>& Polynomial<uint16_t>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<uint32_t>& Polynomial<uint32_t>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %u\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<uint64_t>& Polynomial<uint64_t>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %llu\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %llu\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<short int>& Polynomial<short int>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %d\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %d\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<int>& Polynomial<int>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %d\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %d\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<long int>& Polynomial<long int>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %ld\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %ld\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<long long int>& Polynomial<long long int>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %lld\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %lld\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<float>& Polynomial<float>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %f\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %f\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<double>& Polynomial<double>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %f\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %f\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<long double>& Polynomial<long double>::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, %Lf\n", i, c[i]);
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, %Lf\n", i, c[i]);
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<std::complex<float> >& Polynomial<std::complex<float> >::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, (%f, %f)\n", i, real(c[i]), imag(c[i]));
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, (%f, %f)\n", i, real(c[i]), imag(c[i]));
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<std::complex<double> >& Polynomial<std::complex<double> >::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, (%f, %f)\n", i, real(c[i]), imag(c[i]));
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, (%f, %f)\n", i, real(c[i]), imag(c[i]));
	}
	printf("\n");
	return *this;
}
template<>
Polynomial<std::complex<long double> >& Polynomial<std::complex<long double> >::print() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n", l);
	for(uint64_t i=0; i<t; ++i) {
		printf("%llu, (%Lf, %Lf)\n", i, real(c[i]), imag(c[i]));
	}
	printf("(padding)\n");
	for(uint64_t i=t; i<l; ++i) {
		printf("%llu, (%Lf, %Lf)\n", i, real(c[i]), imag(c[i]));
	}
	printf("\n");
	return *this;
}


//fprint
template<>
Polynomial<bool>& Polynomial<bool>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%u", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%u)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<char>& Polynomial<char>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%u", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%u)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<char16_t>& Polynomial<char16_t>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%u", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%u)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<char32_t>& Polynomial<char32_t>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%u", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%u)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<wchar_t>& Polynomial<wchar_t>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%u", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%u)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<signed char>& Polynomial<signed char>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%d", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%d)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<uint8_t>& Polynomial<uint8_t>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%u", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%u)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<uint16_t>& Polynomial<uint16_t>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%u", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%u)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<uint32_t>& Polynomial<uint32_t>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%u", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%u)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<uint64_t>& Polynomial<uint64_t>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%llu", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%llu)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<short int>& Polynomial<short int>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%d", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%d)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<int>& Polynomial<int>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%d", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%d)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<long int>& Polynomial<long int>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%ld", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%ld)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<long long int>& Polynomial<long long int>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%lld", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%lld)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<float>& Polynomial<float>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%f", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%f)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<double>& Polynomial<double>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%f", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%f)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<long double>& Polynomial<long double>::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[%Lf", l, c[0]);
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%Lf)x^%llu", c[i], i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<std::complex<float> >& Polynomial<std::complex<float> >::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[(%f + i%f)", l, real(c[0]), imag(c[0]));
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%f + i%f)x^%llu", real(c[0]), imag(c[0]), i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<std::complex<double> >& Polynomial<std::complex<double> >::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[(%f + i%f)", l, real(c[0]), imag(c[0]));
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%f + i%f)x^%llu", real(c[0]), imag(c[0]), i);
	}
	printf("]\n");
	return *this;
}
template<>
Polynomial<std::complex<long double> >& Polynomial<std::complex<long double> >::fprint() {
	printf("\nPolynomial %p:\n", this);
	printf("t=%u:\n", t);
	printf("l=%llu:\n[(%Lf + i%Lf)", l, real(c[0]), imag(c[0]));
	for(uint64_t i=1; i<t; ++i) {
		printf(" + (%Lf + i%Lf)x^%llu", real(c[0]), imag(c[0]), i);
	}
	printf("]\n");
	return *this;
}


















template class Polynomial<bool>;

template class Polynomial<char>;
template class Polynomial<char16_t>;
template class Polynomial<char32_t>;
template class Polynomial<wchar_t>;
template class Polynomial<signed char>;

template class Polynomial<uint8_t>;
template class Polynomial<uint16_t>;
template class Polynomial<uint32_t>;
template class Polynomial<uint64_t>;

template class Polynomial<short int>;
template class Polynomial<int>;
template class Polynomial<long int>;
template class Polynomial<long long int>;

template class Polynomial<float>;
template class Polynomial<double>;
template class Polynomial<long double>;



template class Polynomial<std::complex<float> >;
template class Polynomial<std::complex<double> >;
template class Polynomial<std::complex<long double> >;











