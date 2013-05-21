#pragma once

#include <cassert>

#if 0 //defined(_DEBUG) && defined(_WIN32)
#define VEC4_VALIDATE() {	\
	assert(_finite(x));\
	assert(!_isnan(x));\
	\
	assert(_finite(y));\
	assert(!_isnan(y));\
	\
	assert(_finite(z));\
	assert(!_isnan(z));\
	\
	assert(_finite(w));\
	assert(!_isnan(w));\
}
#else
#define VEC4_VALIDATE()
#endif

template <typename T>
class XVector4
{
public:

	typedef T value_type;

	CUDA_CALLABLE XVector4() : x(0), y(0), z(0), w(0) {}
	CUDA_CALLABLE XVector4(T a) : x(a), y(a), z(a), w(a) {}
	CUDA_CALLABLE XVector4(const T* p) : x(p[0]), y(p[1]), z(p[2]), w(p[3]) {}
	CUDA_CALLABLE XVector4(T x_, T y_, T z_, T w_) : x(x_), y(y_), z(z_), w(w_)
	{
		VEC4_VALIDATE();
	}

	CUDA_CALLABLE operator T* () { return &x; }
	CUDA_CALLABLE operator const T* () const { return &x; };

	CUDA_CALLABLE void Set(T x_, T y_, T z_, T w_) { VEC4_VALIDATE(); x = x_; y = y_; z = z_; w = w_; }

	CUDA_CALLABLE XVector4<T> operator * (T scale) const { XVector4<T> r(*this); r *= scale; VEC4_VALIDATE(); return r;}
	CUDA_CALLABLE XVector4<T> operator / (T scale) const { XVector4<T> r(*this); r /= scale; VEC4_VALIDATE(); return r; }
	CUDA_CALLABLE XVector4<T> operator + (const XVector4<T>& v) const { XVector4<T> r(*this); r += v; VEC4_VALIDATE(); return r; }
	CUDA_CALLABLE XVector4<T> operator - (const XVector4<T>& v) const { XVector4<T> r(*this); r -= v; VEC4_VALIDATE(); return r; }
	CUDA_CALLABLE XVector4<T> operator * (XVector4<T> scale) const { XVector4<T> r(*this); r *= scale; VEC4_VALIDATE(); return r; }

	CUDA_CALLABLE XVector4<T>& operator *=(T scale) {x *= scale; y *= scale; z*= scale; w*= scale; VEC4_VALIDATE(); return *this;}
	CUDA_CALLABLE XVector4<T>& operator /=(T scale) {T s(1.0f/scale); x *= s; y *= s; z *= s; w *=s; VEC4_VALIDATE(); return *this;}
	CUDA_CALLABLE XVector4<T>& operator +=(const XVector4<T>& v) {x += v.x; y += v.y; z += v.z; w += v.w; VEC4_VALIDATE(); return *this;}
	CUDA_CALLABLE XVector4<T>& operator -=(const XVector4<T>& v) {x -= v.x; y -= v.y; z -= v.z; w -= v.w; VEC4_VALIDATE(); return *this;}
	CUDA_CALLABLE XVector4<T>& operator *=(const XVector4<T>& v) {x *= v.x; y *= v.y; z *= v.z; w *= v.w; VEC4_VALIDATE(); return *this;}

	CUDA_CALLABLE bool operator != (const XVector4<T>& v) const { return (x != v.x || y != v.y || z != v.z || w != v.w); }

	// negate
	CUDA_CALLABLE XVector4<T> operator -() const { VEC4_VALIDATE(); return XVector4<T>(-x, -y, -z, -w); }

	T x,y,z,w;
};

typedef XVector4<float> Vector4;
typedef XVector4<float> Vec4;

// lhs scalar scale
template <typename T>
CUDA_CALLABLE XVector4<T> operator *(T lhs, const XVector4<T>& rhs)
{
	XVector4<T> r(rhs);
	r *= lhs;
	return r;
}

template <typename T>
CUDA_CALLABLE bool operator==(const XVector4<T>& lhs, const XVector4<T>& rhs)
{
	return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w);
}

