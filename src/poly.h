#include <cstdio>
#include <type_traits>
#include <complex>
#include <cmath>
#define _USE_MATH_DEFINES //legacy feature of C


//TODO: increase precision to long double/parameterize precision?
//TODO: figure out problem with type based control statements


template<typename T>
struct Polynomial {
	
	uint32_t t; //# of terms

	uint64_t l; //next highest power of 2 to t
	T* c; //pointer to array of coefficients

	
	//constructs the 0 polynomial
	Polynomial();
	
	//constructs a polynomial
	/*
		@ terms: number of terms to be in the polynomial
				 (length of coeffs)
		@ coeffs: coefficients of polynomial,
				  size of this array should be >= terms
	*/
	Polynomial(uint32_t terms, T coeffs[]);

	//copies polynomial
	Polynomial(const Polynomial<T>& p);


	~Polynomial();


	//pads polynomial with 0s
	/*
		@ len: pads coefficient array to next highest power of 2 of len

		c will be reallocated
	*/
	Polynomial<T>& pad(uint64_t len);

	//strips excess 0s
	/*
		t will be reset to the highest degree of the polynomial currently +1
		l will be reset to the next highest power of 2 from new t
		c will be reallocated
	*/
	Polynomial<T>& strip();


	
	//sets coefficients
	/*
		@ terms: number of terms
		@ coeffs: array coeffecients that polynomial coefficients are set to
	*/
	Polynomial<T>& setCoeffs(uint64_t len, std::complex<double> coeffs[]);

	//TODO: second setCoeffs funciton for T array?

	

	//FFT recursive helper function
	/*
		@ p: polynomial with complex coefficients, as coeffiecient array
		@ q: pointer to destination of fft output,
			 polynomial with complex coefficients, as coefficient array,
			 size of this should be the same as p
		@ t: # of terms,
			 should be a power of 2 (pad polynomial with 0s if necessary)
		@ inverse: is this fft being used as an inverse fourier transform (wether to use counterclockwise or clockwise roots of unity)
	*/
	static std::complex<double>* recursive_fft(uint64_t n, std::complex<double> p[], std::complex<double>* q, bool inverse = false);

	//FFT algorithm for (Inverse)DFT on polynomial with complex coefficients
	/*
		@ inverse: is this fft being used as an inverse fft (wether to use counterclockwise or clockwise roots of unity)
	*/
	Polynomial<T>& fft(bool inverse = false);
	
	
	//DFT of polynomial with complex coefficients
	/*
		sets this polynomial to its DFT
	*/
	Polynomial<T>& dft();
	
	//Inverse-DFT of polynomial with complex coefficients
	/*
		sets this polynomial to its inverse DFT
	*/
	Polynomial<T>& idft();


	//evaluates polynomial for given value
	/*
		@ x: input value to polynomial
	*/
	T operator()(T x);

	//gets coefficient of polynomial
	/*
		@ i: get coefficient of x^i
	*/
	T operator[](uint32_t i);


	//adds given polynomial
	/*
		@ q: polynomial being added
	*/
	Polynomial<T>& operator+=(const Polynomial<T>& q);
	//adds to constant term
	/*
		@ a: value to add to constant term
	*/
	Polynomial<T>& operator+=(T a);

	//subtracts given polynomial
	/*
		@ q: polynomial being subtracted
	*/
	Polynomial<T>& operator-=(const Polynomial<T>& q);
	//subtracts from constant term
	/*
		@ a: value to subtract from constant term
	*/
	Polynomial<T>& operator-=(T a);

	//multiplies by given polynomial
	/*
		@ q: polynomial to multiply by

		undefined behavior could oocur if this polynomial and q have greater than 2^32 terms
	*/
	Polynomial<T>& operator*=(const Polynomial<T>& q);
	//scales polynomial
	/*
		@ alpha: value to scale the polynomial coefficients by
	*/
	Polynomial<T>& operator*=(T alpha);

	//TODO: implement
	//divides polynomial by the given polynomial
	/*
		@ q: divisor polynomial
	*/
	Polynomial<T>& operator/=(const Polynomial<T>& q);
	//divides all coefficients by given value
	/*
		@ alpha: divisor of all coefficients
	*/
	Polynomial<T>& operator/=(T alpha);


	//TODO: implement
	//takes the polynomial modulo given polynomial
	/*
		@ q: modulus polynomial
	*/
	Polynomial<T>& operator%=(const Polynomial<T>& q);
	//TODO: how to have different functions for different types
	//takes all coefficients modulo given value
	/*
		@ alpha: modulus to apply to all coefficients
	*/
	Polynomial<T>& operator%=(T alpha);


	//casts complex double to T
	/*
		@ z: complex double to be casted to type T
	*/
	T cast(std::complex<double> z);


	//prints all information about polynomial
	Polynomial<T>& print();
	//prints polynomial in format (can be pasted into wolfram alpha)
	Polynomial<T>& fprint();
};










