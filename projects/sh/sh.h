#pragma once

#include <math.h>

// The spherical harmonic basis functions are formed through the combination 
// of the Legrende polynomials (parameterized by cos(theta)) and the Fourier basis functions (based on sine and cosine).

int factorial(int x)
{
   int f = 1;

   while (x > 0)
   {
      f *= x;
      x -= 1;
   }

   return f;
}

int doubleFactorial(int x)
{
   int f = 1;

   while (x > 0)
   {
      f *= x;
      x -= 2;
   }

   return f;
}

int pow(int x, int e)
{
   int a = 1;

   for (int i=0; i < e; ++i)
      a *= x;

   return a;
}
   

double shK(int l, int m)
{
	assert(l >= 0);
	assert(m >= 0);

    return sqrt((2.0*l+1.0)*factorial(l-m)/(4.0*kPi*factorial(l+m)));
}

double shL(int l, int m, double x)
{
	assert(x >= -1.0 && x <= 1.0);
	assert(m >= 0 && m <= l);
	assert(l >= 0);

    if (m == l)
    {
		return pow(-1, m)*doubleFactorial(2*m-1)*pow(1.0-x*x, m*0.5); 
    }
    else if (m == l-1)
    {
		return x*(2*m+1)*shL(m, m, x); 
    }
    else
    {
		return 1.0/(l-m)*(x*(2*l-1)*shL(l-1, m, x)-(l+m-1)*shL(l-2, m, x));
    }
}

double shY(int l, int m, double theta, double phi)
{
	assert(l >= 0);
	assert(m <= l && m >= -l);

   if (m > 0)
   {
      return sqrt(2.0)*shK(l, m)*cos(m*phi)*shL(l, m, cos(theta));
   }
   else if (m < 0)
   {
      return sqrt(2.0)*shK(l, -m)*sin(-m*phi)*shL(l, -m, cos(theta));
   }
   else
   {
      return shK(l, 0)*shL(l, 0, cos(theta));
   }
}

template <int l, int m>
double shY(double theta, double phi) { return shY(l, m, theta, phi); }

// numerically integrate a 1 variable function over a range
template <typename func>
double integrate1(func f, double xstart, double xend, int steps=8192)
{
	double dx = (xend-xstart)/steps;
	double x = xstart;
	double r = 0.0;

	for (int i=0; i < steps; ++i)
	{
		r += f(x)*dx;
		x += dx;
	}

	return r;
}

// numerically integrate a 2 variable function over a range
template<typename func>
double integrate2(func f, double xstart, double xend, double ystart, double yend, int steps=8192)
{
	const double dx = (xend-xstart)/steps;
	const double dy = (yend-ystart)/steps;

	double y = ystart;
	double a = 0.0;

	for (int i=0; i < steps; ++i)
	{
		double x = xstart;

		for (int j=0; j < steps; ++j)
		{
			a += f(x, y)*dx*dy;
			x += dx;
		}

		y += dy;
	}

	return a;
}

int shIndex(int l, int m)
{
	return (l*l)+m+l;
}

template<typename Func, typename Value>
void shProject(const Func& f, int lmax, Value* coefficients)
{
	const int kSteps = 1024;

	const double dTheta = kPi/kSteps;
	const double dPhi = 2.0*kPi/(kSteps*2);

	for (int l=0; l < lmax; ++l)
	{
		for (int m=-l; m <=l; ++m)
		{
			int c = shIndex(l, m);

			double theta = 0.0;

			for (int i=0; i < kSteps; ++i)
			{
				double phi = 0.0;

				for (int j=0; j < kSteps*2; ++j)
				{
					coefficients[c] += f(theta, phi)*float(shY(l, m, theta, phi)*sin(theta)*dTheta*dPhi);
					phi += dPhi;
				}

				theta += dTheta;
			}
		}
	}
}

template <typename Value>
void shReduceRinging(Value coefficient[], int lmax, double lambda)
{
	for (int l=0; l < lmax; ++l)
	{
		for (int m=-l; m <= l; ++m)
		{
			coefficient[shIndex(l, m)] *= 1.0 / (1.0 + lambda*l*l*(l+1.0)*(l+1.0));
		}
	}
}

template <typename Value>
Value shExpand(Value coefficients[], int lmax, double theta, double phi)
{
	Value x = Value();

	for (int l=0; l < lmax; ++l)
	{
		for (int m=-l; m <= l; ++m)
		{
			x += shY(l, m, theta, phi)*coefficients[shIndex(l, m)];
		}
	}

	return x;
}

template <typename Value>
void shConvolve(Value out[], const Value f[], const double h[], int lmax)
{
	for (int l=0; l < lmax; ++l)
	{
		double k = sqrt(4.0*kPi/(2.0*l+1.0));

		for (int m=-l; m <= l; ++m)
		{
			out[shIndex(l, m)] = k*f[shIndex(l, m)]*h[shIndex(l, 0)];
		}
	}
}	

