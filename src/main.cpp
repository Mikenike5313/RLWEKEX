#include <random>



#include "poly.h"



//take a polynomial modulus a polynomial Phi
/*
	@ p: polynomial
	@ n: exponent of Phi, a polynomial of form x^n + 1,
		 p will be taken modulo this Phi
	@ r: polynomial to store remainder
*/
template <typename T>
Polynomial<T>& poly_mod_Phi(Polynomial<T>& p, uint16_t n, Polynomial<T>& r) {
	//initialize r_x to 0
	std::complex<double> r_x[p.l];

	//use rule that x^n + 1 = 0 => x^n = -1
	for(uint32_t i=0; i<p.l; ++i) {
		r_x[i%n] += pow(-1, i/n)*p[i];
	}
	r.setCoeffs(p.l, r_x);
	return r;
}//has precision




//prints samples from a discrete normal distribution in a finite ring
/*
	@ rng: random engine, pre-seeded
	@ d: any normal distribution
	@ m: modulus
	@ num_samples: # of samples to take
	@ scl: scale printed histogram (# samples equating to a "*" in histogram)
*/
void print_discrete_normal_distribution_sample(std::default_random_engine rng, std::normal_distribution<double> d, uint32_t m, uint64_t num_samples = 10000, uint16_t scl = 200) {
	uint64_t sample[m];
	memset(sample, 0, m*sizeof(uint32_t));
	for(uint64_t i = 0; i < num_samples; ++i) {
		int c = (int)round(d(rng))%(int)m;
		while(c < 0) {
			c+=m;
		}
		c = c%(int)m;
		++sample[c];
	}

	printf("\n");
	for(uint32_t i = 0; i < m; ++i) {
		if(sample[i] != 0) {
			printf("\n % 6u| ", i);
			for(uint64_t j = 1; j < sample[i]; ++j) {
				if(j % scl == 0) {
					printf("*");
				}
			}
		}
	}
	printf("\n");
}









int main(int argc, char* argv[]) {


	
	//for 128 bits of security (Singh)
	/*
	uint16_t n = 512; //power of 2
	uint64_t q = 25601; // 1 mod 4, 1 mod 2n
	//*/
	
	//for 256 bits of security (Singh)
	/*
	uint16_t n = 1024; //power of 2
	uint64_t q = 40961; // 1 mod 4, 1 mod 2n
	//*/

	//for 200 bits of security (NewHope)
	///*
	uint16_t n = 1024; //power of 2
	uint32_t q = 12289; // 1 mod 4, 1 mod 2n
	//*/


	//small example
	/*
	uint16_t n = 8; //power of 2
	uint32_t q = 17; // 1 mod 4, 1 mod 2n
	//*/






	//rng setup
	
	std::random_device rd;
	std::default_random_engine generator(rd());


	double sigma = 4*M_SQRT1_2*M_2_SQRTPI;
	std::normal_distribution<double> chi_a(0, sigma); //gaussian sampling

	uint8_t b = 5; // bound used if uniform sampling is used
				   // 5 reccomended by Singh
	std::uniform_int_distribution<int> uniform_b(-b, b); //uniform sampling for small polynomials


	std::uniform_int_distribution<int> uniform_a(-(int)q/2, (int)q/2); //uniform sampling for public polynomial

	

	//WORKING RLWE-KEX IMPLEMENTATION (Inspiration from - GitHub: https://github.com/vscrypto/ringlwe)

	//polynomial of size n where every coefficient is q
	long long add_q_x[n];
	memset(add_q_x, 0, n*sizeof(long long));
	for(uint16_t i=0; i<n; ++i) {
		add_q_x[i] = q;
	}
	Polynomial<long long> add_q(n, add_q_x);




	//public polynomial
	long long a_x[n];
	memset(a_x, 0, n*sizeof(long long));
	for(uint16_t i=0; i<n; ++i) {
		a_x[i] = (long long)uniform_a(generator)%(long long)q;
	}
	Polynomial<long long> a(n, a_x);
	//a.fprint();
	

	//GENERATE
	
	//Alice
	//alice's secret
	long long s_A_x[n];
	memset(s_A_x, 0, n*sizeof(long long));
	for(uint16_t i=0; i<n; ++i) {
		s_A_x[i] = (long long)round(chi_a(generator))%(long long)q;
	}
	Polynomial<long long> s_A(n, s_A_x);
	//s_A.fprint();

	//alice's error
	long long e_A_x[n];
	memset(e_A_x, 0, n*sizeof(long long));
	for(uint16_t i=0; i<n; ++i) {
		e_A_x[i] = (long long)round(chi_a(generator))%(long long)q;
	}
	Polynomial<long long> e_A(n, e_A_x);
	//e_A.fprint();

	//alice's public key
	Polynomial<long long> b_A = a;
	b_A *= s_A;
	poly_mod_Phi(b_A, n, b_A);
	b_A += e_A;
	b_A %= q;
	//b_A.fprint();


	//ENCAPSULATE

	//Bob
	//bob's secret
	long long s_B_x[n];
	memset(s_B_x, 0, n*sizeof(long long));
	for(uint16_t i=0; i<n; ++i) {
		s_B_x[i] = (long long)round(chi_a(generator))%(long long)q;
	}
	Polynomial<long long> s_B(n, s_B_x);
	//s_B.fprint();

	//bob's error
	long long e_B_x[n];
	memset(e_B_x, 0, n*sizeof(long long));
	for(uint16_t i=0; i<n; ++i) {
		e_B_x[i] = (long long)round(chi_a(generator))%(long long)q;
	}
	Polynomial<long long> e_B(n, e_B_x);
	//e_B.fprint();

	//bob's public key
	Polynomial<long long> b_B = a;
	b_B *= s_B;
	poly_mod_Phi(b_B, n, b_B);
	b_B += e_A;
	b_B %= q;
	//b_B.fprint();

	//bob's 2nd error
	long long ep_B_x[n];
	memset(ep_B_x, 0, n*sizeof(long long));
	for(uint16_t i=0; i<n; ++i) {
		ep_B_x[i] = (long long)round(chi_a(generator))%(long long)q;
	}
	Polynomial<long long> ep_B(n, ep_B_x);
	//ep_B.fprint();

	//bob's initial key
	Polynomial<long long> k_B = b_A;
	k_B *= s_B;
	poly_mod_Phi(k_B, n, k_B);
	k_B += ep_B;
	k_B %= q;
	k_B += add_q; //make sure all coeffs are positive
	k_B %= q;

	//find bob's shared key with modular rounding
	bool sk_B[n];
	memset(sk_B, 0, n);
	for(uint16_t i=0; i<n; ++i) {
		sk_B[i] = (k_B[i] > (long long)q/4 && k_B[i] < 3*(long long)q/4);
	}

	//find reconciliation info using cross rounding
	bool w[n];
	memset(w, 0, n);
	for(uint16_t i=0; i<n; ++i) {
		w[i] = ((k_B[i] > (long long)q/4 && k_B[i] <= (long long)q/2) || k_B[i] >= 3*(long long)q/4);
	}


	//DECAPSULATE

	//alice's initial key
	Polynomial<long long> k_A = b_B;
	k_A *= s_A;
	poly_mod_Phi(k_A, n, k_A);
	k_A %= q;
	k_A += add_q; //make sure all coeffs are positive
	k_A %= q;

	//find alice's shared key with reconciliation info
	bool sk_A[n];
	memset(sk_A, 0, n);
	for(uint16_t i=0; i<n; ++i) {
		if(w[i]) {
			//printf("%u: %lld <= %llu <= %lld --> %u\n", i, (long long)q/8, k_A[i], 5*(long long)q/8, (k_A[i] >= (long long)q/8 && k_A[i] <= 5*(long long)q/8));
			sk_A[i] = (k_A[i] >= (long long)q/8 && k_A[i] <= 5*(long long)q/8);
		}
		else {
			//printf("%u: %lld < %llu <= %lld --> %u\n", i, 3*(long long)q/8, k_A[i], 7*(long long)q/8, (k_A[i] >= 3*(long long)q/8 && k_A[i] <= 7*(long long)q/8));
			sk_A[i] = (k_A[i] > 3*(long long)q/8 && k_A[i] <= 7*(long long)q/8);
		}
	}




	//check that the keys are in fact established
	printf("\nKeys:\n");
	printf("  Bit| _A_|_B_\n");
	for(uint16_t i=0; i<n; ++i) {
		printf("% 5u|  %u | %u\n", i, sk_A[i], sk_B[i]);
	}
	printf("\n(end keys)\n\n");

	//check if the keys are the same
	uint16_t errors = 0;
	for(uint16_t i=0; i<n; ++i) {
		if(sk_A[i] != sk_B[i]) {
			++errors;
		}
	}
	printf("\n%u errors in keys, at bits:\n", errors);
	for(uint16_t i=0; i<n; ++i) {
		if(sk_A[i] != sk_B[i]) {
			printf("% 5d: % 1u != % 1u\n", i, sk_A[i], sk_B[i]);
		}
	}
	printf("\n(end errors)\n\n");

	

	










	printf("\n\nfinished main!");
	return 0;
}























