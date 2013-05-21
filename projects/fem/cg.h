#include "core/types.h"
#include "core/maths.h"
#include "mat22.h"

#include <vector>


typedef std::vector<Matrix22> NodeRow;

std::vector<Vec2> CgMul(const std::vector<NodeRow>& A, const std::vector<Vec2>& x)
{
	assert(A[0].size() == x.size());

	std::vector<Vec2> b(A.size());

	for (size_t i=0; i < A.size(); ++i)
	{
		for (size_t j=0; j < x.size(); ++j)
		{
			b[i] += A[i][j]*x[j];
		}
	}	

	return b;
}

std::vector<Vec2> CgMul(float c, const std::vector<Vec2>& x)
{
	std::vector<Vec2> r(x.size());

	for (size_t i=0; i < r.size(); ++i)
		r[i] = c*x[i];

	return r;
}

std::vector<Vec2> CgAdd(const std::vector<Vec2>& a, const std::vector<Vec2>& b)
{
	assert(a.size() == b.size());

	std::vector<Vec2> c(a.size());

	for (size_t i=0; i < a.size(); ++i)
		c[i] = a[i] + b[i];

	return c;
}

std::vector<Vec2> CgSub(const std::vector<Vec2>& a, const std::vector<Vec2>& b)
{
	assert(a.size() == b.size());

	std::vector<Vec2> c(a.size());

	for (size_t i=0; i < a.size(); ++i)
		c[i] = a[i] - b[i];

	return c;
}

float CgDot(const std::vector<Vec2>& a, const std::vector<Vec2>& b)
{
	assert(a.size() == b.size());

	float d = 0.0f;

	for (size_t i=0; i < a.size(); ++i)
	{
		d += Dot(a[i], b[i]);
	}

	return d;
}

void CgDebug(const Matrix22& m)
{
	printf("{ %f, %f }\n", m(0, 0), m(0, 1));
	printf("{ %f, %f }\n", m(1, 0), m(1, 1));
}

void CgDebug(const char* s, const std::vector<Vec2>& x)
{
	for (size_t i=0; i < x.size(); ++i)
		printf("%s[%u] = {%f, %f}\n", s, uint32_t(i), x[i].x, x[i].y);
}

void CgDebug(const char* s, const std::vector<NodeRow>& A)
{
	for (size_t i=0; i < A.size(); ++i)
	{

		for (int r=0; r < 2; ++r)
		{
			for (size_t j=0; j < A[i].size(); ++j)
			{
				printf("%f, %f, ", A[i][j](r,0), A[i][j](r, 1));
			}
			printf("\n");
		}
	}
}


// returns x
std::vector<Vec2> CgSolve(const std::vector<NodeRow>& A, const std::vector<Vec2>& b, uint32_t imax, float e)
{
	assert(A.size() == A[0].size());

	uint32_t i=0;
	std::vector<Vec2> x(A[0].size());
	std::vector<Vec2> r = CgSub(b, CgMul(A, x));
	std::vector<Vec2> d = r;

	float sigmaNew = CgDot(r, r); 
	//const float sigma0 = sigmaNew;

	assert(A.size() == A[0].size());
	assert(x.size() == A.size());
	assert(r.size() == A.size());
	assert(d.size() == A.size());
	
	//CgDebug("b", b);
	std::vector<Vec2> q;

	while (i < imax)// && sigmaNew > e*e*sigma0)
	{
		q = CgMul(A, d);
		
	//	CgDebug("q", q);

		float a = sigmaNew / CgDot(d, q);
		x = CgAdd(x, CgMul(a,d));

		// could reseed residual here to account for numerical inaccuracy
		r = CgSub(r, CgMul(a, q));

		float sigmaOld = sigmaNew;
		if (sigmaOld <= 0.0001f)
			return x;
		sigmaNew = CgDot(r, r);

		//printf("%d %f\n", i, sigmaNew);
	
		float beta = sigmaNew / sigmaOld;

		d = CgAdd(r, CgMul(beta, d));

		++i;
	}

	return x;	
}


