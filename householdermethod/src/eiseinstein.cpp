#include "eiseinstein.h"

Eiseinstein::Eiseinstein(){
	a = 0;
	b = 0;
}

Eiseinstein::Eiseinstein(mpreal aval, mpreal bval){
	a = aval;
	b = bval;
}


mpreal Eiseinstein::norm() const {
	return a*a - a * b + b * b;
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

Eiseinstein Eiseinstein::operator-() const {
        return Eiseinstein(-a, -b);
};

Eiseinstein Eiseinstein::operator *(const Eiseinstein& other){

	mpreal newa = a * other.a - b * other.b;
	mpreal newb = a * other.b + b * other.a - b * other.b;
	return Eiseinstein(newa, newb);
};


Eiseinstein Eiseinstein::pow(int k){
	Eiseinstein result(1,0);
	for (int i = 1; i <= k; i++){

		result = (*this) * result;
	}
	return result;
};

pair<Eiseinstein, Eiseinstein> Eiseinstein::operator /(const Eiseinstein& other){

	mpreal const norm = other.norm();

	if (norm != 0)
	{

		Eiseinstein numerator = (*this) * other.conjugate();
	        mpreal u = static_cast<mpreal>(numerator.a)/norm;
	        mpreal v = static_cast<mpreal>(numerator.b)/norm;
	        // rounding off u and v
	        mpreal m = static_cast<mpreal>(ceil(u - 0.5));
	        mpreal n = static_cast<mpreal>(ceil(v - 0.5));
	        Eiseinstein quotient(m,n);
	        Eiseinstein remainder = (*this) - quotient * other;

	        return make_pair(quotient, remainder);

	} 

	else{
		throw std::invalid_argument("division by zero");
	}

};

pair<Eiseinstein, Eiseinstein> Eiseinstein::operator /(const mpreal& other){
	Eiseinstein new_other(other, 0);
	return (*this)/new_other;
}


complex<mpreal> Eiseinstein::value(){
	mpreal three(3.0);
	return complex<mpreal>(a - 0.5 * b, sqrt(three) * 0.5 * b);
}


void Eiseinstein::print() const {
        cout << a << " + " << b << "*ω" << endl;
};