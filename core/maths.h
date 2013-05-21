#pragma once

#include <cmath>
#include <float.h>
#include <cassert>

#include "core.h"
#include "types.h"

#ifdef __CUDACC__
#define CUDA_CALLABLE __host__ __device__
#else
#define CUDA_CALLABLE
#endif

const float kPi = 3.141592653589f;
const float k2Pi = 2.0f*kPi;
const float kInvPi = 1.0f/kPi;
const float kInv2Pi = 0.5f/kPi;
const float kDegToRad = kPi/180.0f;
const float kRadToDeg = 180.0f/kPi;

CUDA_CALLABLE inline float DegToRad(float t)
{
	return t * kDegToRad;
}

CUDA_CALLABLE inline float RadToDeg(float t)
{
	return t * kRadToDeg;
}

CUDA_CALLABLE inline float Sin(float theta)
{
	return sinf(theta);
}

CUDA_CALLABLE inline float Cos(float theta)
{
	return cosf(theta);
}

CUDA_CALLABLE inline void SinCos(float theta, float& s, float& c)
{
	// no optimizations yet
	s = sinf(theta);
	c = cosf(theta);
}

CUDA_CALLABLE inline float Tan(float theta)
{
	return tanf(theta);
}

CUDA_CALLABLE inline float Sqrt(float x)
{
	return sqrtf(x);
}

CUDA_CALLABLE inline double Sqrt(double x)
{
	return sqrt(x);
}

CUDA_CALLABLE inline float ASin(float theta)
{
	return asinf(theta);
}

CUDA_CALLABLE inline float ACos(float theta)
{
	return acosf(theta);
}

CUDA_CALLABLE inline float ATan(float theta)
{
	return atanf(theta);
}

CUDA_CALLABLE inline float ATan2(float x, float y)
{
	return atan2f(x, y);
}

CUDA_CALLABLE inline float Abs(float x)
{
	return fabsf(x);
}

CUDA_CALLABLE inline float Pow(float b, float e)
{
	return powf(b, e);
}

CUDA_CALLABLE inline float Sgn(float x)
{
	return (x < 0.0f ? -1.0f : 1.0f);
}

CUDA_CALLABLE inline float Sign(float x)
{
	return x < 0.0f ? -1.0f : 1.0f;
}

CUDA_CALLABLE inline double Sign(double x)
{
	return x < 0.0f ? -1.0f : 1.0f;
}

CUDA_CALLABLE inline float Mod(float x, float y)
{
	return fmod(x, y);
}

template <typename T>
CUDA_CALLABLE inline T Min(T a, T b)
{
	return a < b ? a : b;
}

template <typename T>
CUDA_CALLABLE inline T Max(T a, T b)
{
	return a > b ? a : b;
}

template <typename T>
CUDA_CALLABLE inline T Clamp(T a, T low, T high)
{
	if (low > high)
		Swap(low, high);
	
	return Max(low, Min(a, high));
}

template <typename V, typename T> 
CUDA_CALLABLE inline V Lerp(const V& start, const V& end, const T& t)
{
	return start + (end-start)*t;
}

template <typename T> 
CUDA_CALLABLE inline void Swap(T& a, T& b)
{
	T tmp = a;
	a = b;
	b = tmp;
}

CUDA_CALLABLE inline float InvSqrt(float x)
{
	return 1.0f / sqrtf(x);
}

// round towards +infinity
CUDA_CALLABLE inline int Round(float f)
{
	return int(f+0.5f);
}

template <typename T>
CUDA_CALLABLE T Normalize(const T& v)
{
	T a(v);
	a /= Length(v);
	return a;
}

template <typename T>
CUDA_CALLABLE inline typename T::value_type LengthSq(const T v)
{
	return Dot(v,v);
}

template <typename T>
CUDA_CALLABLE inline typename T::value_type Length(const T& v)
{
	typename T::value_type lSq = LengthSq(v);
	if (lSq)
		return Sqrt(LengthSq(v));
	else
		return 0.0f;
}

// this is mainly a helper function used by script
template <typename T>
CUDA_CALLABLE inline typename T::value_type Distance(const T& v1, const T& v2)
{
	return Length(v1-v2);
}

template <typename T>
CUDA_CALLABLE inline T SafeNormalize(const T& v, const T& fallback=T())
{
	float l = LengthSq(v);
	if (l > 0.0f)
	{
		return v * InvSqrt(l);
	}
	else
		return fallback;
}

#include "vec2.h"
#include "vec3.h"
#include "vec4.h"
#include "point3.h"
#include "mat22.h"
#include "mat33.h"
#include "mat44.h"

// represents a plane in the form ax + by + cz - d = 0
class Plane : public Vec4
{
public:

	Plane() {}
	Plane(float x, float y, float z, float w) : Vec4(x, y, z, w) {}
	Plane(const Point3& p, const Vector3& n)
	{
		x = n.x;
		y = n.y;
		z = n.z;
		w = -Dot3(p, n);
	}

	Vec3 GetNormal() const { return Vec3(x, y, z); }
	Point3 GetPoint() const { return Point3(x*-w, y*-w, z*-w); }

	Plane(const Vec3& v) : Vec4(v.x, v.y, v.z, 1.0f) {}
	Plane(const Vec4& v) : Vec4(v) {}
};

template <typename T>
inline T Dot(const XVector4<T>& v1, const XVector4<T>& v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z + v1.w*v2.w;
}

// helper function that assumes a w of 0
inline float Dot(const Plane& p, const Vector3& v)
{
	return p.x*v.x + p.y*v.y + p.z*v.z;
}

// helper function that assumes a w of 1
inline float Dot(const Plane& p, const Point3& v)
{
	return p.x*v.x + p.y*v.y + p.z*v.z + p.w;
}

//----------------------------------------------------------------------------
inline float RandomUnit()
{
	float r = (float)rand();
	r /= RAND_MAX;
	return r;
}

// Random number in range [-1,1]
inline float RandomSignedUnit()
{
	float r = (float)rand();
	r /= RAND_MAX;
	r = 2.0f * r - 1.0f;
	return r;
}

inline float Random(float lo, float hi)
{
	float r = (float)rand();
	r /= RAND_MAX;
	r = (hi - lo) * r + lo;
	return r;
}

extern uint32_t seed1;
extern uint32_t seed2;

void RandInit();

// random number generator
inline uint32_t Rand()
{
	seed1 = ( seed2 ^ ( ( seed1 << 5 ) | ( seed1 >> 27 ) ) ) ^ ( seed1*seed2 );
	seed2 = seed1 ^ ( ( seed2 << 12 ) | ( seed2 >> 20 ) );

	return seed1;
}

// returns a random number in the range [min, max)
inline uint32_t Rand(uint32_t min, uint32_t max)
{
	return min + Rand()%(max-min);
}

// returns random number between 0-1
inline float Randf()
{
	uint32_t value = Rand();
	uint32_t limit = 0xffffffff;

	return ( float )value*( 1.0f/( float )limit );
}

// returns random number between min and max
inline float Randf(float min, float max)
{
	//	return Lerp(min, max, ParticleRandf());
	float t = Randf();
	return (1.0f-t)*min + t*(max);
}

// returns random number between 0-max
inline float Randf(float max)
{
	return Randf()*max;
}

// returns a random unit vector (also can add an offset to generate around an off axis vector)
inline Vec3 RandomUnitVector()
{
	float phi = Randf(kPi*2.0f);
	float theta = Randf(kPi*2.0f);

	float cosTheta = Cos(theta);
	float sinTheta = Sin(theta);

	float cosPhi = Cos(phi);	
	float sinPhi = Sin(phi);

	return Vec3(cosTheta*sinPhi,cosPhi,sinTheta*sinPhi);
}

inline Vec3 RandVec3() { return Vec3(Randf(-1.0f, 1.0f), Randf(-1.0f, 1.0f), Randf(-1.0f, 1.0f)); }

inline Vec3 UniformSampleSphere()
{
	float u1 = Randf(0.0f, 1.0f);
	float u2 = Randf(0.0f, 1.0f);

	float z = 1.f - 2.f * u1;
	float r = sqrtf(Max(0.f, 1.f - z*z));
	float phi = 2.f * kPi * u2;
	float x = r * cosf(phi);
	float y = r * sinf(phi);

	return Vector3(x, y, z);
}

inline Vec3 UniformSampleHemisphere()
{
	// generate a random z value
	float z = Randf(0.0f, 1.0f);
	float w = Sqrt(1.0f-z*z);

	float phi = k2Pi*Randf(0.0f, 1.0f);
	float x = Cos(phi)*w;
	float y = Sin(phi)*w;

	return Vec3(x, y, z);
}

inline Vec2 UniformSampleDisc()
{
	float r = Sqrt(Randf(0.0f, 1.0f));
	float theta = k2Pi*Randf(0.0f, 1.0f);

	return Vec2(r * Cos(theta), r * Sin(theta));
}

inline void UniformSampleTriangle(float& u, float& v)
{
	float r = Sqrt(Randf());
	u = 1.0f - r;
	v = Randf() * r;
}

inline Vec3 CosineSampleHemisphere()
{
	Vec2 s = UniformSampleDisc();
	float z = Sqrt(Max(0.0f, 1.0f - s.x*s.x - s.y*s.y));

	return Vec3(s.x, s.y, z);
}

inline Vec3 SphericalToXYZ(float theta, float phi)
{
	float cosTheta = cos(theta);
	float sinTheta = sin(theta);

	return Vec3(sin(phi)*sinTheta, cosTheta, cos(phi)*sinTheta);
}

// returns random vector between -range and range
inline Vec4 Randf(const Vec4 &range)
{
	return Vec4(Randf(-range.x, range.x), 
			    Randf(-range.y, range.y),
			    Randf(-range.z, range.z),
			  	Randf(-range.w, range.w));
}

// generates a transform matrix with v as the z axis, taken from PBRT
inline void BasisFromVector(const Vec3& w, Vec3* u, Vec3* v)
{
	if (fabsf(w.x) > fabsf(w.y))
	{
		float invLen = 1.0f / sqrtf(w.x*w.x + w.z*w.z);
		*u = Vec3(-w.z*invLen, 0.0f, w.x*invLen);
	}
	else
	{
		float invLen = 1.0f / sqrtf(w.y*w.y + w.z*w.z);
		*u = Vec3(0.0f, w.z*invLen, -w.y*invLen);
	}

	*v = Cross(w, *u);	

	assert(fabsf(Length(*u)-1.0f) < 0.01f);
	assert(fabsf(Length(*v)-1.0f) < 0.01f);
}

// same as above but returns a matrix
inline Mat44 TransformFromVector(const Vec3& w, const Point3& t=Point3(0.0f, 0.0f, 0.0f))
{
	Mat44 m = Mat44::kIdentity;
	m.SetCol(2, Vec4(w.x, w.y, w.z, 0.0));
	m.SetCol(3, Vec4(t.x, t.y, t.z, 1.0f));

	BasisFromVector(w, (Vec3*)m.columns[0], (Vec3*)m.columns[1]);

	return m;
}

// todo: sort out rotations
inline Mat44 ViewMatrix(const Point3& pos) 
{
	
	float view[4][4] = { { 1.0f, 0.0f, 0.0f, 0.0f },
						{ 0.0f, 1.0f, 0.0f, 0.0f },
						{ 0.0f, 0.0f, 1.0f, 0.0f },
						{ -pos.x, -pos.y, -pos.z, 1.0f } };
	
	return Mat44(&view[0][0]);
}


inline Mat44 LookAtMatrix(const Point3& viewer, const Point3& target)
{
	// create a basis from viewer to target (OpenGL convention looking down -z)
	Vec3 forward = -Normalize(target-viewer);
	Vec3 up(0.0f, 1.0f, 0.0f);
	Vec3 left = Cross(up, forward);
	up = Cross(forward, left);

	float xform[4][4] = {  { left.x, left.y, left.z, 0.0f },
	  					     { up.x, up.y, up.z, 0.0f},
						     { forward.x, forward.y, forward.z, 0.0f},
							     { viewer.x, viewer.y, viewer.z, 1.0f} };

	return AffineInverse(Mat44(&xform[0][0]));		
}

// generate a rotation matrix around an axis, from PBRT p74
inline Mat44 RotationMatrix(float angle, const Vec3& axis)
{
	Vec3 a = Normalize(axis);
	float s = sinf(angle);
	float c = cosf(angle);

	float m[4][4];

	m[0][0] = a.x * a.x + (1.0f - a.x * a.x) * c;
	m[0][1] = a.x * a.y * (1.0f - c) + a.z * s;
	m[0][2] = a.x * a.z * (1.0f - c) - a.y * s;
	m[0][3] = 0.0f;

	m[1][0] = a.x * a.y * (1.0f - c) - a.z * s;
	m[1][1] = a.y * a.y + (1.0f - a.y * a.y) * c;
	m[1][2] = a.y * a.z * (1.0f - c) + a.x * s;
	m[1][3] = 0.0f;

	m[2][0] = a.x * a.z * (1.0f - c) + a.y * s;
	m[2][1] = a.y * a.z * (1.0f - c) - a.x * s;
	m[2][2] = a.z * a.z + (1.0f - a.z * a.z) * c;
	m[2][3] = 0.0f;

	m[3][0] = 0.0f;
	m[3][1] = 0.0f;
	m[3][2] = 0.0f;
	m[3][3] = 1.0f;

	return Mat44(&m[0][0]);
}

inline Mat44 TranslationMatrix(const Point3& t)
{
	Mat44 m(Mat44::kIdentity);
	m.SetTranslation(t);
	return m;
}

inline Mat44 OrthographicMatrix(float left, float right, float bottom, float top, float n, float f)
{
	
	float m[4][4] = { { 2.0f/(right-left), 0.0f, 0.0f, 0.0f },
					  { 0.0f, 2.0f/(top-bottom), 0.0f, 0.0f },			
					  { 0.0f, 0.0f, -2.0f/(f-n), 0.0f },
					  { -(right+left)/(right-left), -(top+bottom)/(top-bottom), -(f+n)/(f-n), 1.0f } };
	

	return Mat44(&m[0][0]);
}

// this is designed as a drop in replacement for gluPerspective
inline Mat44 ProjectionMatrix(float fov, float aspect, float znear, float zfar) 
{
	float f = 1.0f / tanf(DegToRad(fov*0.5f));
	float zd = znear-zfar;

	float view[4][4] = { { f/aspect, 0.0f, 0.0f, 0.0f },
						 { 0.0f, f, 0.0f, 0.0f },
						 { 0.0f, 0.0f, (zfar+znear)/zd, -1.0f },
						 { 0.0f, 0.0f, (2.0f*znear*zfar)/zd, 0.0f } };
 
	return Mat44(&view[0][0]);
}

// encapsulates an orientation encoded in Euler angles, not the sexiest 
// representation but it is convenient when manipulating objects from script

class Rotation
{
public:

	Rotation() : yaw(0), pitch(0), roll(0) {}
	Rotation(float inYaw, float inPitch, float inRoll) : yaw(inYaw), pitch(inPitch), roll(inRoll) {}

	Rotation& operator +=(const Rotation& rhs) {yaw += rhs.yaw; pitch += rhs.pitch; roll += rhs.roll; return *this;}
	Rotation& operator -=(const Rotation& rhs) {yaw -= rhs.yaw; pitch -= rhs.pitch; roll -= rhs.roll; return *this;}

	Rotation operator + (const Rotation& rhs) const { Rotation lhs(*this); lhs += rhs; return lhs; }
	Rotation operator - (const Rotation& rhs) const { Rotation lhs(*this); lhs -= rhs; return lhs; }

	// all members are in degrees (easy editing)
	float yaw;
	float pitch;
	float roll;
};

inline Mat44 ScaleMatrix(const Vector3& s)
{
	float m[4][4] = { {s.x, 0.0f, 0.0f, 0.0f },
					  { 0.0f, s.y, 0.0f, 0.0f},
					  { 0.0f, 0.0f, s.z, 0.0f},
					  { 0.0f, 0.0f, 0.0f, 1.0f} };

	return Mat44(&m[0][0]);
}

// assumes yaw on y, then pitch on z, then roll on x
inline Mat44 TransformMatrix(const Rotation& r, const Point3& p)
{
	const float yaw = DegToRad(r.yaw);
	const float pitch = DegToRad(r.pitch);
	const float roll = DegToRad(r.roll);

	const float s1 = Sin(roll);
	const float c1 = Cos(roll);
	const float s2 = Sin(pitch);
	const float c2 = Cos(pitch);
	const float s3 = Sin(yaw);
	const float c3 = Cos(yaw);

	// interprets the angles as yaw around world-y, pitch around new z, roll around new x
	// i feel this is the most natural mapping for the majority of games i'll be making
	// technically these are "Tait-Bryan rotations" (subsequent rotations depend on the previous ones)

	// you can derive this matrix by constructing the three rotation matrices according to yaw, pitch, roll
	// and multiplying them in in reverse order 
	float mr[4][4] = {	{ c2*c3, s2, -c2*s3, 0.0f},
						{ s1*s3-c1*c3*s2, c1*c2, c3*s1+c1*s2*s3, 0.0f},
						{ c3*s1*s2+c1*s3, -c2*s1, c1*c3-s1*s2*s3, 0.0f},	
						{ p.x, p.y, p.z, 1.0f} };

	Mat44 m1(&mr[0][0]);

	return m1;//m2 * m1;
}

class Transform
{
public:

	Transform() {};
	Transform(const Point3& p, const Rotation& r) : position(p), rotation(r) {}

	Mat44 ToMatrix() const
	{
		return TransformMatrix(rotation, position);
	}

	// helper function to translate object
	void Translate(const Vec3& delta)
	{
		position += delta;
	}

	// helper function to rotate an object
	void Rotate(const Rotation& delta)
	{
		rotation += delta;
	}

	void RotateToLookAt(const Point3& target)
	{
		// get vector to point
		Vec3 delta = target-position;

		// yaw
		rotation.yaw = atan2f(delta.z, delta.x);
		rotation.pitch = atan2f(delta.y, sqrt(delta.x*delta.x+delta.z*delta.z));
		rotation.roll = 0.0f;
	}

	Vec3 GetXAxis() const
	{
		return ToMatrix().columns[0];
	}

	Vec3 GetYAxis() const
	{
		return ToMatrix().columns[1];
	}

	Vec3 GetZAxis() const
	{
		return ToMatrix().columns[2];
	}

	Point3 position;
	Rotation rotation;
};

// aligns the z axis along the vector
inline Rotation AlignToVector(const Vec3& vector)
{
	// todo: fix, see spherical->cartesian coordinates wikipedia
	return Rotation(0.0f, RadToDeg(atan2(vector.y, vector.x)), 0.0f);
}

// creates a vector given an angle measured clockwise from horizontal (1,0)
inline Vec2 AngleToVector(float a)
{	
	return Vec2(Cos(a), Sin(a));
}

inline float VectorToAngle(const Vec2& v)
{
	return atan2f(v.y, v.x);
}

CUDA_CALLABLE inline float SmoothStep(float a, float b, float t)
{
	t = Clamp(t-a / (b-a), 0.0f, 1.0f);
	return t*t*(3.0f-2.0f*t);
}

// hermite spline interpolation
template <typename T>
T HermiteInterpolate(const T& a, const T& b, const T& t1, const T& t2, float t)
{
	// blending weights
	const float w1 = 1.0f - 3*t*t + 2*t*t*t;
	const float w2 = t*t*(3.0f-2.0f*t);
	const float w3 = t*t*t - 2*t*t + t;
	const float w4 = t*t*(t-1.0f);

	// return weighted combination
	return a*w1 + b*w2 + t1*w3 + t2*w4;

}
	
template <typename T>
T HermiteTangent(const T& a, const T& b, const T& t1, const T& t2, float t)
{
	// first derivative blend weights
	const float w1 = 6.0f*t*t-6*t;
	const float w2 = -6.0f*t*t + 6*t;
	const float w3 = 3*t*t - 4*t + 1;
	const float w4 = 3*t*t - 2*t;

	// weighted combination
	return a*w1 + b*w2 + t1*w3 + t2*w4;
}

template <typename T>
T HermiteSecondDerivative(const T& a, const T& b, const T& t1, const T& t2, float t)
{
	// first derivative blend weights
	const float w1 = 12*t - 6.0f;
	const float w2 = -12.0f*t + 6;
	const float w3 = 6*t - 4.0f;
	const float w4 = 6*t - 2.0f;

	// weighted combination
	return a*w1 + b*w2 + t1*w3 + t2*w4;
}

inline float Log(float base, float x)
{
	// calculate the log of a value for an arbitary base, only use if you can't use the standard bases (10, e)
	return logf(x) / logf(base);
}

// function which maps a value to a range
template <typename T>
T RangeMap(T value, T lower, T upper)
{
	assert(upper >= lower);
	return (value-lower)/(upper-lower);
}

// simple colour class
class Colour 
{
public:

	enum Preset
	{
		kRed,
		kGreen,
		kBlue,
		kWhite,
		kBlack
	};
	
	Colour(float r_=0.0f, float g_=0.0f, float b_=0.0f, float a_=1.0f) : r(r_), g(g_), b(b_), a(a_) {}
	Colour(float* p) : r(p[0]), g(p[1]), b(p[2]), a(p[3]) {}
	Colour(uint32_t rgba)
	{
		a = ((rgba)&0xff)/255.0f;
		r = ((rgba>>24)&0xff)/255.0f;
		g = ((rgba>>16)&0xff)/255.0f;
		b = ((rgba>>8)&0xff)/255.0f;
	}
	Colour(Preset p);

	// cast operator
	operator const float*() const { return &r; }
	operator float*() { return &r; }

	Colour operator * (float scale) const { Colour r(*this); r *= scale; return r; }
	Colour operator / (float scale) const { Colour r(*this); r /= scale; return r; }
	Colour operator + (const Colour& v) const { Colour r(*this); r += v; return r; }
	Colour operator - (const Colour& v) const { Colour r(*this); r -= v; return r; }
	Colour operator * (const Colour& scale) const { Colour r(*this); r *= scale; return r;}

	Colour& operator *=(float scale) {r *= scale; g *= scale; b*= scale; a*= scale;  return *this;}
	Colour& operator /=(float scale) {float s(1.0f/scale); r *= s; g *= s; b *= s; a *=s; return *this;}
	Colour& operator +=(const Colour& v) {r += v.r; g += v.g; b += v.b; a += v.a; return *this;}
	Colour& operator -=(const Colour& v) {r -= v.r; g -= v.g; b -= v.b; a -= v.a; return *this;}
	Colour& operator *=(const Colour& v) {r *= v.r; g *= v.g; b *= v.b; a *= v.a; return *this;}

	float r, g, b, a;

};

inline bool operator == (const Colour& lhs, const Colour& rhs)
{
	return lhs.r == rhs.r && lhs.g == rhs.g && lhs.b == rhs.b && lhs.a == rhs.a;
}

inline bool operator != (const Colour& lhs, const Colour& rhs)
{
	return !(lhs == rhs);
}

inline Colour ToneMap(const Colour& s)
{
	//return Colour(s.r / (s.r+1.0f),	s.g / (s.g+1.0f), s.b / (s.b+1.0f), 1.0f);
	float Y = 0.3333f*(s.r + s.g + s.b);
	return s / (1.0f + Y);
}

// lhs scalar scale
inline Colour operator * (float lhs, const Colour& rhs)
{
	Colour r(rhs);
	r *= lhs;
	return r;
}

inline Colour YxyToXYZ(float Y, float x, float y)
{
	float X = x * (Y / y);
	float Z = (1.0f - x - y) * Y / y;

	return Colour(X, Y, Z, 1.0f);
}

inline Colour HSVToRGB( float h, float s, float v )
{
	float r, g, b;

	int i;
	float f, p, q, t;
	if( s == 0 ) {
		// achromatic (grey)
		r = g = b = v;
	}
	else
	{
		h *= 6.0f;			// sector 0 to 5
		i = int(floor( h ));
		f = h - i;			// factorial part of h
		p = v * ( 1 - s );
		q = v * ( 1 - s * f );
		t = v * ( 1 - s * ( 1 - f ) );
		switch( i ) {
			case 0:
				r = v;
				g = t;
				b = p;
				break;
			case 1:
				r = q;
				g = v;
				b = p;
				break;
			case 2:
				r = p;
				g = v;
				b = t;
				break;
			case 3:
				r = p;
				g = q;
				b = v;
				break;
			case 4:
				r = t;
				g = p;
				b = v;
				break;
			default:		// case 5:
				r = v;
				g = p;
				b = q;
				break;
		};
	}

	return Colour(r, g, b);
}

inline Colour XYZToLinear(float x, float y, float z)
{
	float c[4];
	c[0] =  3.240479f * x + -1.537150f * y + -0.498535f * z;
	c[1] = -0.969256f * x +  1.875991f * y +  0.041556f * z;
	c[2] =  0.055648f * x + -0.204043f * y +  1.057311f * z;
	c[3] = 1.0f;

	return Colour(c);
}

inline uint32_t ColourToRGBA8(const Colour& c)
{
	union SmallColor
	{
		uint8_t u8[4];
		uint32_t u32;
	};

	SmallColor s;
	s.u8[0] = (uint8_t)(Clamp(c.r, 0.0f, 1.0f) * 255);
	s.u8[1] = (uint8_t)(Clamp(c.g, 0.0f, 1.0f) * 255);
	s.u8[2] = (uint8_t)(Clamp(c.b, 0.0f, 1.0f) * 255);
	s.u8[3] = (uint8_t)(Clamp(c.a, 0.0f, 1.0f) * 255);

	return s.u32;
}

inline Colour LinearToSrgb(const Colour& c)
{
	const float kInvGamma = 1.0f/2.2f;
	return Colour(powf(c.r, kInvGamma), powf(c.g, kInvGamma), powf(c.b, kInvGamma), c.a); 
}

inline Colour SrgbToLinear(const Colour& c)
{
	const float kInvGamma = 2.4f;
	return Colour(powf(c.r, kInvGamma), powf(c.g, kInvGamma), powf(c.b, kInvGamma), c.a); 
}

// intersection routines
inline bool IntersectRaySphere(const Point3& sphereOrigin, float sphereRadius, const Point3& rayOrigin, const Vec3& rayDir, float& t, Vec3* hitNormal=NULL)
{
	Vec3 d(sphereOrigin-rayOrigin);
	float deltaSq = LengthSq(d);
	float radiusSq = sphereRadius*sphereRadius;

	// if the origin is inside the sphere return no intersection
	if (deltaSq > radiusSq)
	{
		float dprojr = Dot(d, rayDir);
		
		// if ray pointing away from sphere no intersection
		if (dprojr < 0.0f)
			return false;

		// bit of Pythagoras to get closest point on ray
		float dSq = deltaSq-dprojr*dprojr;
		
		if (dSq > radiusSq)
			return false;
		else
		{
			// length of the half cord
			float thc = sqrt(radiusSq-dSq);
			
			// closest intersection
			t = dprojr - thc;

			// calculate normal if requested
			if (hitNormal)
				*hitNormal = Normalize((rayOrigin+rayDir*t)-sphereOrigin);

			return true;
		}
	}
	else
	{
		return false;
	}
}

template <typename T>
CUDA_CALLABLE inline bool SolveQuadratic(T a, T b, T c, T& minT, T& maxT)
{
	if (a == 0.0f && b == 0.0f)
	{
		minT = maxT = 0.0f;
		return true;
	}

	T discriminant = b*b - T(4.0)*a*c;

	if (discriminant < 0.0f)
	{
		return false;
	}

	// numerical receipes 5.6 (this method ensures numerical accuracy is preserved)
	T t = T(-0.5) * (b + Sign(b)*Sqrt(discriminant));
	minT = t / a;
	maxT = c / t;

	if (minT > maxT)
	{
		Swap(minT, maxT);
	}

	return true;
}

// alternative ray sphere intersect, returns closest and furthest t values
inline bool IntersectRaySphere(const Point3& sphereOrigin, float sphereRadius, const Point3& rayOrigin, const Vector3& rayDir, float& minT, float &maxT, Vec3* hitNormal=NULL)
{
	Vector3 q = rayOrigin-sphereOrigin;

	float a = 1.0f;
	float b = 2.0f*Dot(q, rayDir);
	float c = Dot(q, q)-(sphereRadius*sphereRadius);

	bool r = SolveQuadratic(a, b, c, minT, maxT);

	if (minT < 0.0)
		minT = 0.0f;

	// calculate the normal of the closest hit
	if (hitNormal && r)
	{
		*hitNormal = Normalize((rayOrigin+rayDir*minT)-sphereOrigin);
	}

	return r;
}

inline bool IntersectRayPlane(const Point3& p, const Vector3& dir, const Plane& plane, float& t)
{
    float d = Dot(plane, dir);
    
    if (d == 0.0f)
    {
        return false;
    }
	else
    {
        t = -Dot(plane, p) / d;
    }

	return (t > 0.0f);	
}

inline bool IntersectLineSegmentPlane(const Vec3& start, const Vec3& end, const Plane& plane, Vec3& out)
{
	Vec3 u(end - start);

	float dist = -Dot(plane, start) / Dot(plane, u);

	if (dist > 0.0f && dist < 1.0f)
	{
		out = (start + u * dist);
		return true;
	}
	else
		return false;
}

// Moller and Trumbore's method
inline bool IntersectRayTriTwoSided(const Point3& p, const Vec3& dir, const Point3& a, const Point3& b, const Point3& c, float& t, float& u, float& v, float& w, float& sign)//Vec3* normal)
{
    Vector3 ab = b - a;
    Vector3 ac = c - a;
    Vector3 n = Cross(ab, ac);

    float d = Dot(-dir, n);
    float ood = 1.0f / d; // No need to check for division by zero here as infinity aritmetic will save us...
    Vector3 ap = p - a;

    t = Dot(ap, n) * ood;
    if (t < 0.0f)
        return false;

    Vector3 e = Cross(-dir, ap);
    v = Dot(ac, e) * ood;
    if (v < 0.0f || v > 1.0f) // ...here...
        return false;
    w = -Dot(ab, e) * ood;
    if (w < 0.0f || v + w > 1.0f) // ...and here
        return false;

    u = 1.0f - v - w;
    //if (normal)
        //*normal = n;
	sign = d;

    return true;
}



// mostly taken from Real Time Collision Detection - p192
inline bool IntersectRayTri(const Point3& p, const Vec3& dir, const Point3& a, const Point3& b, const Point3& c,  float& t, float& u, float& v, float& w, Vec3* normal)
{
	const Vec3 ab = b-a;
	const Vec3 ac = c-a;

	// calculate normal
	Vec3 n = Cross(ab, ac);

	// need to solve a system of three equations to give t, u, v
	float d = Dot(-dir, n);

	// if dir is parallel to triangle plane or points away from triangle 
	if (d <= 0.0f)
        return false;

	Vec3 ap = p-a;
	t = Dot(ap, n);

	// ignores tris behind 
	if (t < 0.0f)
		return false;

	// compute barycentric coordinates
	Vec3 e = Cross(-dir, ap);
	v = Dot(ac, e);
	if (v < 0.0f || v > d) return false;

	w = -Dot(ab, e);
	if (w < 0.0f || v + w > d) return false;

	float ood = 1.0f / d;
	t *= ood;
	v *= ood;
	w *= ood;
	u = 1.0f-v-w;

	// optionally write out normal (todo: this branch is a performance concern, should probably remove)
	if (normal)
		*normal = n;

	return true;
}

inline float minf(const float a, const float b) { return a < b ? a : b; }
inline float maxf(const float a, const float b) { return a > b ? a : b; }

inline bool IntersectRayAABBOmpf(const Point3& pos, const Vector3& rcp_dir, const Vector3& min, const Vector3& max, float& t) {
       
    float
        l1	= (min.x - pos.x) * rcp_dir.x,
        l2	= (max.x - pos.x) * rcp_dir.x,
        lmin	= minf(l1,l2),
        lmax	= maxf(l1,l2);

    l1	= (min.y - pos.y) * rcp_dir.y;
    l2	= (max.y - pos.y) * rcp_dir.y;
    lmin	= maxf(minf(l1,l2), lmin);
    lmax	= minf(maxf(l1,l2), lmax);

    l1	= (min.z - pos.z) * rcp_dir.z;
    l2	= (max.z - pos.z) * rcp_dir.z;
    lmin	= maxf(minf(l1,l2), lmin);
    lmax	= minf(maxf(l1,l2), lmax);

    //return ((lmax > 0.f) & (lmax >= lmin));
    //return ((lmax > 0.f) & (lmax > lmin));
    bool hit = ((lmax >= 0.f) & (lmax >= lmin));
    if (hit)
        t = lmin;
    return hit;
}


inline bool IntersectRayAABB(const Point3& start, const Vector3& dir, const Vector3& min, const Vector3& max, float& t, Vector3* normal)
{
	//! calculate candidate plane on each axis
	float tx = -1.0f, ty = -1.0f, tz = -1.0f;
	bool inside = true;
			
	//! use unrolled loops

	//! x
	if (start.x < min.x)
	{
		if (dir.x != 0.0f)
			tx = (min.x-start.x)/dir.x;
		inside = false;
	}
	else if (start.x > max.x)
	{
		if (dir.x != 0.0f)
			tx = (max.x-start.x)/dir.x;
		inside = false;
	}

	//! y
	if (start.y < min.y)
	{
		if (dir.y != 0.0f)
			ty = (min.y-start.y)/dir.y;
		inside = false;
	}
	else if (start.y > max.y)
	{
		if (dir.y != 0.0f)
			ty = (max.y-start.y)/dir.y;
		inside = false;
	}

	//! z
	if (start.z < min.z)
	{
		if (dir.z != 0.0f)
			tz = (min.z-start.z)/dir.z;
		inside = false;
	}
	else if (start.z > max.z)
	{
		if (dir.z != 0.0f)
			tz = (max.z-start.z)/dir.z;
		inside = false;
	}

	//! if point inside all planes
	if (inside)
    {
        t = 0.0f;
		return true;
    }

	//! we now have t values for each of possible intersection planes
	//! find the maximum to get the intersection point
	float tmax = tx;
	int taxis = 0;

	if (ty > tmax)
	{
		tmax = ty;
		taxis = 1;
	}
	if (tz > tmax)
	{
		tmax = tz;
		taxis = 2;
	}

	if (tmax < 0.0f)
		return false;

	//! check that the intersection point lies on the plane we picked
	//! we don't test the axis of closest intersection for precision reasons

	//! no eps for now
	float eps = 0.0f;

	Point3 hit = start + dir*tmax;

	if ((hit.x < min.x-eps || hit.x > max.x+eps) && taxis != 0)
		return false;
	if ((hit.y < min.y-eps || hit.y > max.y+eps) && taxis != 1)
		return false;
	if ((hit.z < min.z-eps || hit.z > max.z+eps) && taxis != 2)
		return false;

	//! output results
	t = tmax;
			
	return true;
}
	
// 2d rectangle class
class Rect
{
public:

	Rect() : m_left(0), m_right(0), m_top(0), m_bottom(0) {}

	Rect(uint32_t left, uint32_t right, uint32_t top, uint32_t bottom) : m_left(left), m_right(right), m_top(top), m_bottom(bottom)
	{
		assert(left <= right);
		assert(top <= bottom);	
	}

	uint32_t Width() const { return m_right - m_left; }
	uint32_t Height() const { return m_bottom - m_top; }

	// expand rect x units in each direction
	void Expand(uint32_t x)
	{
		m_left -= x;
		m_right += x;
		m_top -= x;
		m_bottom += x;
	}

	uint32_t Left() const { return m_left; }
	uint32_t Right() const { return m_right; }
	uint32_t Top() const { return m_top; }
	uint32_t Bottom() const { return m_bottom; }

	bool Contains(uint32_t x, uint32_t y) const
	{
		return (x >= m_left) && (x <= m_right) && (y >= m_top) && (y <= m_bottom);
	}

	uint32_t m_left;
	uint32_t m_right;
	uint32_t m_top;
	uint32_t m_bottom;
};

// doesn't really belong here but efficient (and I believe correct) in place random shuffle based on the Fisher-Yates / Knuth algorithm
template <typename T>
void RandomShuffle(T begin, T end)
{
	assert(end > begin);
	uint32_t n = distance(begin, end);

	for (uint32_t i=0; i < n; ++i)
	{
		// pick a random number between 0 and n-1
		uint32_t r = Rand() % (n-i);

		// swap that location with the current randomly selected position
		swap(*(begin+i), *(begin+(i+r)));
	}
}


CUDA_CALLABLE inline Vec4 QuatFromAxisAngle(const Vec3& axis, float angle)
{
	Vec3 v = Normalize(axis);

	float half = angle*0.5f;
	float w = cosf(half);

	const float sin_theta_over_two = sinf(half);
	v *= sin_theta_over_two;

	return Vec4(v.x, v.y, v.z, w);
}


// rotate by quaternion (q, w)
CUDA_CALLABLE inline Vec3 rotate(const Vec3& q, float w, const Vec3& x)
{
	return x*(2.0f*w*w-1.0f) + Cross(q, x)*w*2.0f + q*Dot(q, x)*2.0f;
}

// rotate x by inverse transform in (q, w)
CUDA_CALLABLE inline Vec3 rotateInv(const Vec3& q, float w, const Vec3& x)
{
	return x*(2.0f*w*w-1.0f) - Cross(q, x)*w*2.0f + q*Dot(q, x)*2.0f;
}