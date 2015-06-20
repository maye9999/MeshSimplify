#ifndef _Matrix_h__
#define _Matrix_h__

#include "Vector.h"

struct Matrix
{
	double a[4][4];

	Matrix()
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				a[i][j] = 0;
			}
		}
	}

	void clear() {for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) a[i][j] = 0;}

	double& operator()(int i, int j) {return a[i][j];}
	const double& operator()(int i, int j) const {return a[i][j]; }

	Matrix& operator+=(const Matrix& b)
	{
		for(int i = 0; i < 4; ++i)
			for(int j = 0; j < 4; ++j)
				a[i][j] += b(i, j);
		return *this;
	}
	Matrix& operator-=(const Matrix& b)
	{
		for(int i = 0; i < 4; ++i)
			for(int j = 0; j < 4; ++j)
				a[i][j] -= b(i, j);
		return *this;
	}

	friend const Matrix operator+(const Matrix& lhs, const Matrix& rhs);
	friend const Matrix operator-(const Matrix& lhs, const Matrix& rhs);

	static const Vector mul(const Matrix& m, const Vector& v)
	{
		Vector ret(0, 0, 0, 0);
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				ret[i] += m(i, j) * v[j];
			}
		}
		return ret;
	}

	static const Vector mul(const Vector& v_t, const Matrix& m)
	{
		Vector ret(0, 0, 0, 0);
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				ret[i] += m(j, i) * v_t[j];
			}
		}
		return ret;
	}

	static void mul_fast(const Vector& v, const Vector& v_t, Matrix& m)
	{
		for(int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; j++)
			{
				m(i, j) += v[i] * v_t[j];
			}
	}

	static const Matrix invert(const Matrix& mm)
	{
		Matrix m = mm;
		int i, j, k;

		for (i = 1; i < 4; i++)
		{
			m.a[0][i] /= m.a[0][0];
		}

		for (i = 1; i < 4; i++)
		{
			for (j = i; j < 4; j++)
			{
				double sum = 0.0;
				for (k = 0; k < i; k++)
				{
					sum += m.a[j][k] * m.a[k][i];
				}
				m.a[j][i] -= sum;
			}
			if (i == 3) continue;
			for (j=i+1; j < 4; j++)
			{
				double sum = 0.0;
				for (int k = 0; k < i; k++)
					sum += m.a[i][k]*m.a[k][j];
				m.a[i][j] = 
					(m.a[i][j]-sum) / m.a[i][i];
			}
		}

		//
		// Invert L
		//
		for ( i = 0; i < 4; i++ )
		{
			for ( int j = i; j < 4; j++ )
			{
				double x = 1.0;
				if ( i != j )
				{
					x = 0.0;
					for ( int k = i; k < j; k++ ) 
						x -= m.a[j][k]*m.a[k][i];
				}
				m.a[j][i] = x / m.a[j][j];
			}
		}

		//
		// Invert U
		//
		for ( i = 0; i < 4; i++ )
		{
			for ( j = i; j < 4; j++ )
			{
				if ( i == j ) continue;
				float sum = 0.0;
				for ( int k = i; k < j; k++ )
					sum += m.a[k][j]*( (i==k) ? 1.0 : m.a[i][k] );
				m.a[i][j] = -sum;
			}
		}

		//
		// Final Inversion
		//
		for ( i = 0; i < 4; i++ )
		{
			for ( int j = 0; j < 4; j++ )
			{
				float sum = 0.0;
				for ( int k = ((i>j)?i:j); k < 4; k++ )  
					sum += ((j==k)?1.0:m.a[j][k])*m.a[k][i];
				m.a[j][i] = sum;
			}
		}
		return m;
	}

};

inline const Matrix operator+(const Matrix& lhs, const Matrix& rhs)
{
	Matrix ret;
	for(int i = 0; i < 4; ++i)
		for(int j = 0; j < 4; ++j)
			ret.a[i][j] = lhs(i, j) + rhs(i, j);
	return ret;
}

inline const Matrix operator-(const Matrix& lhs, const Matrix& rhs)
{
	Matrix ret;
	for(int i = 0; i < 4; ++i)
		for(int j = 0; j < 4; ++j)
			ret.a[i][j] = lhs(i, j) - rhs(i, j);
	return ret;
}
#endif // Matrix_h__
