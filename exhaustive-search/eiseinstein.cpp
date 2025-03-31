#include "eiseinstein.h"

Eiseinstein::Eiseinstein(long long aval, long long bval) : a(aval), b(bval) {}

long long Eiseinstein::norm() const {
	return a*a - a*b + b*b;
};


Eiseinstein Eiseinstein::conjugate() const {
	return Eiseinstein(a-b, -b);
};

Eiseinstein Eiseinstein::operator +(const Eiseinstein& other){
	return Eiseinstein(a + other.a, b + other.b);
};

Eiseinstein Eiseinstein::operator -(const Eiseinstein& other){
	return Eiseinstein(a - other.a, b - other.b);
};

Eiseinstein Eiseinstein::operator *(const Eiseinstein& other){

	long long newa = a * other.a - b * other.b;
	long long newb = a * other.b + b * other.a - b * other.b;
	return Eiseinstein(newa, newb);
};

Eiseinstein Eiseinstein::operator-() const {
        return Eiseinstein(-a, -b);
};

Eiseinstein Eiseinstein::pow(int k){
	Eiseinstein result(1,0);
	for (int i = 1; i <= k; i++){

		result = (*this) * result;
	}
	return result;
};

pair<Eiseinstein, Eiseinstein> Eiseinstein::operator /(const Eiseinstein& other){

	long long const norm = other.norm();

	if (norm != 0)
	{

		Eiseinstein numerator = (*this) * other.conjugate();
	        double u = static_cast<double>(numerator.a)/norm;
	        double v = static_cast<double>(numerator.b)/norm;
	        // rounding off u and v
	        long long m = static_cast<long long>(ceil(u - 0.5));
	        long long n = static_cast<long long>(ceil(v - 0.5));
	        Eiseinstein quotient(m,n);
	        Eiseinstein remainder = (*this) - quotient * other;

	        return make_pair(quotient, remainder);

	} 

	else{
		throw std::invalid_argument("division by zero");
	}

};

pair<Eiseinstein, Eiseinstein> Eiseinstein::operator /(const long long& other){
	Eiseinstein new_other(other, 0);
	return (*this)/new_other;
}


int Eiseinstein::get_sde(int exponent){
	if (exponent == 0) return 0;

	Eiseinstein chi(1,2);
	Eiseinstein factor(1,0);

	for (int i = exponent; i > -1; i--){

		auto division_data = (*this) / factor;
		if (division_data.second.norm() != 0) return i + 1; 


		factor = factor * chi;

	}

	return exponent; // we should never get here
}


void Eiseinstein::print() const {
        cout << a << " + " << b << "*ω" << endl;
};