#pragma once

// stores column vectors in column major order
template <typename T>
class XMatrix44
{
public:

	XMatrix44() { memset(columns, 0, sizeof(columns)); }
	XMatrix44(const T* d) { assert(d); memcpy(columns, d, sizeof(*this)); }
	XMatrix44(T c11, T c21, T c31, T c41,
				 T c12, T c22, T c32, T c42,
			    T c13, T c23, T c33, T c43,
				 T c14, T c24, T c34, T c44)
	{
		columns[0][0] = c11;
		columns[0][1] = c21;
		columns[0][2] = c31;
		columns[0][3] = c41;

		columns[1][0] = c12;
		columns[1][1] = c22;
		columns[1][2] = c32;
		columns[1][3] = c42;

		columns[2][0] = c13;
		columns[2][1] = c23;
		columns[2][2] = c33;
		columns[2][3] = c43;

		columns[3][0] = c14;
		columns[3][1] = c24;
		columns[3][2] = c34;
		columns[3][3] = c44;
	}

	operator T* () { return &columns[0][0]; }
	operator const T* () const { return &columns[0][0]; }

	// right multiply
	XMatrix44<T> operator * (const XMatrix44<T>& rhs) const
	{
		XMatrix44<T> r;
		MatrixMultiply(*this, rhs, r);
		return r;
	}

	// right multiply
	XMatrix44<T>& operator *= (const XMatrix44<T>& rhs)
	{
		XMatrix44<T> r;
		MatrixMultiply(*this, rhs, r);
		*this = r;

		return *this;
	}

	// scalar multiplication
	XMatrix44<T>& operator *= (const T& s)
	{
		for (int c=0; c < 4; ++c)
		{
			for (int r=0; r < 4; ++r)
			{
				columns[c][r] *= s;
			}
		}

		return *this;
	}

	void MatrixMultiply(const T* __restrict lhs, const T* __restrict rhs, T* __restrict result) const
	{
		assert(lhs != rhs);
		assert(lhs != result);
		assert(rhs != result);
		
		for (int i=0; i < 4; ++i)
		{
			for (int j=0; j < 4; ++j)
			{
				result[j*4+i]  = rhs[j*4+0]*lhs[i+0]; 
				result[j*4+i] += rhs[j*4+1]*lhs[i+4];
				result[j*4+i] += rhs[j*4+2]*lhs[i+8];
				result[j*4+i] += rhs[j*4+3]*lhs[i+12];
			}
		}
	}	

	void SetCol(int index, const Vec4& c)
	{
		columns[index][0] = c.x;
		columns[index][1] = c.y;
		columns[index][2] = c.z;
		columns[index][3] = c.w;
	}

	// convenience overloads
	void SetAxis(uint32_t index, const XVector3<T>& a)
	{
		columns[index][0] = a.x;
		columns[index][1] = a.y;
		columns[index][2] = a.z;
		columns[index][3] = 0.0f;
	}

	void SetTranslation(const Point3& p)
	{
		columns[3][0] = p.x;	
		columns[3][1] = p.y;
		columns[3][2] = p.z;
		columns[3][3] = 1.0f;
	}

	const Vec3& GetAxis(int i) const { return *reinterpret_cast<const Vec3*>(&columns[i]); }
	const Vec4& GetCol(int i) const { return *reinterpret_cast<const Vec4*>(&columns[i]); }
	const Point3& GetTranslation() const { return *reinterpret_cast<const Point3*>(&columns[3]); }

	Vec4 GetRow(int i) const { return Vec4(columns[0][i], columns[1][i], columns[2][i], columns[3][i]); }

	float columns[4][4];

	static XMatrix44<T> kIdentity;

};

// right multiply a point assumes w of 1
template <typename T>
Point3 Multiply(const XMatrix44<T>& mat, const Point3& v)
{
	Point3 r;
	r.x = v.x*mat[0] + v.y*mat[4] + v.z*mat[8] + mat[12];
	r.y = v.x*mat[1] + v.y*mat[5] + v.z*mat[9] + mat[13];
	r.z = v.x*mat[2] + v.y*mat[6] + v.z*mat[10] + mat[14];

	return r;
}

// right multiply a vector3 assumes a w of 0
template <typename T>
XVector3<T> Multiply(const XMatrix44<T>& mat, const XVector3<T>& v)
{
	XVector3<T> r;
	r.x = v.x*mat[0] + v.y*mat[4] + v.z*mat[8];
	r.y = v.x*mat[1] + v.y*mat[5] + v.z*mat[9];
	r.z = v.x*mat[2] + v.y*mat[6] + v.z*mat[10];

	return r;
}

// right multiply a vector4
template <typename T>
XVector4<T> Multiply(const XMatrix44<T>& mat, const XVector4<T>& v)
{
	XVector4<T> r;
	r.x = v.x*mat[0] + v.y*mat[4] + v.z*mat[8] + v.w*mat[12];
	r.y = v.x*mat[1] + v.y*mat[5] + v.z*mat[9] + v.w*mat[13];
	r.z = v.x*mat[2] + v.y*mat[6] + v.z*mat[10] + v.w*mat[14];
	r.w = v.x*mat[3] + v.y*mat[7] + v.z*mat[11] + v.w*mat[15];

	return r;
}

template <typename T>
Point3 operator*(const XMatrix44<T>& mat, const Point3& v)
{
	return Multiply(mat, v);
}

template <typename T>
XVector4<T> operator*(const XMatrix44<T>& mat, const XVector4<T>& v)
{
	return Multiply(mat, v);
}

template <typename T>
XVector3<T> operator*(const XMatrix44<T>& mat, const XVector3<T>& v)
{
	return Multiply(mat, v);
}

template<typename T>
inline XMatrix44<T> Transpose(const XMatrix44<T>& m)
{
	XMatrix44<float> inv;

	// transpose
	for (uint32_t c=0; c < 4; ++c)
	{
		for (uint32_t r=0; r < 4; ++r)
		{
			inv.columns[c][r] = m.columns[r][c];
		}
	}

	return inv;
}

template <typename T>
XMatrix44<T> AffineInverse(const XMatrix44<T>& m)
{
	XMatrix44<T> inv;
	
	// transpose upper 3x3
	for (int c=0; c < 3; ++c)
	{
		for (int r=0; r < 3; ++r)
		{
			inv.columns[c][r] = m.columns[r][c];
		}
	}
	
	// multiply -translation by upper 3x3 transpose
	inv.columns[3][0] = -Dot3(m.columns[3], m.columns[0]);
	inv.columns[3][1] = -Dot3(m.columns[3], m.columns[1]);
	inv.columns[3][2] = -Dot3(m.columns[3], m.columns[2]);
	inv.columns[3][3] = 1.0f;

	return inv;	
}

// convenience
typedef XMatrix44<float> Mat44;
typedef XMatrix44<float> Matrix44;

