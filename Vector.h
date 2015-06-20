#ifndef _Vector_h__
#define _Vector_h__

#include <iostream>
#include <cmath>

struct Vector
{
	explicit Vector(double x = 0, double y = 0, double z = 0, double w = 1)
	{
		a[0] = x;
		a[1] = y;
		a[2] = z;
		a[3] = w;
	}
	double a[4];

	double& operator[](int x)	{return a[x];}
	const double& operator[](int x) const {return a[x];}
	
	// Member Binary Operator
	Vector& operator+=(const Vector& v);
	Vector& operator-=(const Vector& v);

	// Non-member Binary Operator
	friend const Vector operator+(const Vector& lhs, const Vector& rhs);
	friend const Vector operator-(const Vector& lhs, const Vector& rhs);
	friend const double operator*(const Vector& lhs, const Vector& rhs);
	friend const Vector operator*(const Vector& v, const double k);
	friend const Vector operator*(const double k, const Vector& v);
	friend std::ostream& operator<<(std::ostream& os, const Vector& v);

	static const Vector cross(const Vector& lhs, const Vector& rhs)
	{
		return Vector(lhs[1]*rhs[2] - lhs[2]*rhs[1], lhs[2]*rhs[0] - lhs[0]*rhs[2], lhs[0]*rhs[1] - lhs[1]*rhs[0]);
	}
	static const Vector normalize(const Vector& v)
	{
		double sq = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		return Vector(v[0] / sq, v[1] / sq, v[2] / sq);
	}
	static const double distance(const Vector& lhs, const Vector& rhs)
	{
		Vector v = lhs - rhs;
		return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	}
};

inline const Vector operator+(const Vector& lhs, const Vector& rhs)
{
	Vector ret;
	for(int i = 0; i < 4; ++i)
		ret[i] = lhs[i] + rhs[i];
	ret[3] = 1;
	return ret;
}

inline const Vector operator-(const Vector& lhs, const Vector& rhs)
{
	Vector ret;
	for(int i = 0; i < 4; ++i)
		ret[i] = lhs[i] - rhs[i];
	ret[3] = 1;
	return ret;
}

inline const double operator*(const Vector& lhs, const Vector& rhs)
{
	double ret = 0;
	for(int i = 0; i < 4; ++i)
		ret += lhs[i] * rhs[i];
	return ret;
}

inline const Vector operator*(const Vector& v, const double k)
{
	return Vector(v.a[0]*k, v.a[1]*k, v.a[2]*k);
}

inline const Vector operator*(const double k, const Vector& v)
{
	return Vector(v.a[0]*k, v.a[1]*k, v.a[2]*k);
}

inline std::ostream& operator<<(std::ostream& os, const Vector& v)
{
	os << "(" << v.a[0] << ", " << v.a[1] << ", " << v.a[2] << ")";
	return os;
}

inline Vector& Vector::operator+=(const Vector& v)
{
	for (int i = 0; i < 4; ++i)
	{
		a[i] += v[i];
	}
	return *this;
}

inline Vector& Vector::operator-=(const Vector& v)
{
	for (int i = 0; i < 4; ++i)
	{
		a[i] -= v[i];
	}
	return *this;
}


#endif // Vector_h__
