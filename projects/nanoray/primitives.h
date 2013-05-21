#pragma once

#include "core/core.h"
#include "core/maths.h"
#include "core/mesh.h"
#include "core/aabbtree.h"

#include "materials.h"

class Primitive
{
public:

    Primitive(Material* m) : m_material(m), m_emission(0.0f) {}

	virtual bool Intersect(const Point3& rayOrigin, const Vector3& rayDir, float& t, Vector3* normal) const = 0;
	virtual Point3 Sample(Vector3* normal) const { assert(0); return Point3(0.0f); };
	virtual float Area() const { assert(0); return 0.0f; }

	Material* m_material;

	Colour m_emission;		// flux / sq ?
};

class SpherePrimitive : public Primitive
{
public:

    SpherePrimitive(Material* mat, const Point3& center, float radius)
        : Primitive(mat)
        , m_position(center)
        , m_radius(radius)
	{		
	}

	virtual bool Intersect(const Point3& rayOrigin, const Vector3& rayDir, float& t, Vector3* normal) const
	{
		return IntersectRaySphere(m_position, m_radius, rayOrigin, rayDir, t, normal);
	}

	// return a random point on the primitive
	virtual Point3 Sample(Vector3* normal) const
	{
		// generate a random point on the sphere
		assert(0);
		return Point3(0.0f);
	}

	virtual float Area() const
	{
		return 4.0f*kPi*m_radius*m_radius;
	}

private:

	Point3 m_position;
	float m_radius;
};

typedef AABBTree SpatialTree;

class MeshPrimitive : public Primitive
{
public:

    MeshPrimitive(Material* mat, Mesh* mesh)
        : Primitive(mat)
        , m_mesh(mesh)
		, m_totalArea(0.0f)
		, m_cdf(mesh->GetNumFaces())
    {
        // create AABB tree 
        m_aabbTree = new SpatialTree(&mesh->m_positions[0], mesh->GetNumVertices(), &mesh->m_indices[0], mesh->GetNumFaces());

		float cdf = 0.0f;

		for (uint32_t i=0; i < m_mesh->GetNumFaces(); ++i)
		{
			Point3& a = m_mesh->m_positions[m_mesh->m_indices[i*3+0]];
			Point3& b = m_mesh->m_positions[m_mesh->m_indices[i*3+1]];
			Point3& c = m_mesh->m_positions[m_mesh->m_indices[i*3+2]];

			cdf += 0.5f*Length(Cross(b-a, c-a));

			m_cdf[i] = cdf;
		}

		// store total area
		m_totalArea = cdf;

		// convert cumulative areas to cdfs
		const float m_rcpCdf = 1.0f / cdf;
		
		for (uint32_t i=0; i < m_mesh->GetNumFaces(); ++i)
		{
			m_cdf[i] *= m_rcpCdf;
		}

		assert(abs(1.0f-m_cdf.back()) < 0.0001f);

	}

    virtual bool Intersect(const Point3& rayOrigin, const Vector3& rayDir, float& t, Vector3* normal) const
    {
        uint32_t triIndex;
        float u, v, w, s;
        
        if (m_aabbTree->TraceRay(rayOrigin, rayDir, t, u, v, w, s, triIndex))
        {
            // interpolate vertex normals
            uint32_t i0 = m_mesh->m_indices[triIndex*3+0];
            uint32_t i1 = m_mesh->m_indices[triIndex*3+1];
            uint32_t i2 = m_mesh->m_indices[triIndex*3+2];
            
            *normal = Sign(s)*Normalize(u*m_mesh->m_normals[i0] + v*m_mesh->m_normals[i1] + w*m_mesh->m_normals[i2]);

            return true;
        }

        return false;
    }

	inline uint32_t PickTri(float random) const 
	{
		float cdf = 0.0f;
		uint32_t face=0;

		for(; face < m_mesh->GetNumFaces(); ++face)
		{
			if (random <= m_cdf[face])
				break;			
		}

		return std::min(face, m_mesh->GetNumFaces()-1);
	}

	inline uint32_t PickTriFast(float random) const
	{
		std::vector<float>::const_iterator i = std::lower_bound(m_cdf.begin(), m_cdf.end(), random);
		return (i-m_cdf.begin());

		/*
		uint32_t low = 0;
		uint32_t high = m_cdf.size()-1;

		// binary search
		while (low < high)
		{
			uint32_t mid = (low+high)/2;

			if (random <= m_cdf[mid])
			{
				high = mid;
			}
			else if (random > m_cdf[mid+1])
			{
				low = mid+1;
			}
			else
			{				
				return mid;
			}
		}

		return low;
		*/
	}

	inline uint32_t PickTriRef(float random) const
	{
		float cdf = 0.0f;
		uint32_t face=0;

		for(; face < m_mesh->GetNumFaces(); ++face)
		{
			Point3& a = m_mesh->m_positions[m_mesh->m_indices[face*3+0]];
			Point3& b = m_mesh->m_positions[m_mesh->m_indices[face*3+1]];
			Point3& c = m_mesh->m_positions[m_mesh->m_indices[face*3+2]];
			
			float pdf = 0.5f*Length(Cross(b-a, c-a))/m_totalArea;

			cdf += pdf;
			
			if (random <= cdf)
				break;			
		}

		return std::min(face, m_mesh->GetNumFaces()-1);
	}

	// return a random point on the primitive
    virtual Point3 Sample(Vector3* normal) const
    {
		const uint32_t face = PickTriFast(Randf());

		float u, v;
		UniformSampleTriangle(u, v);

		Vector3 a = m_mesh->m_positions[m_mesh->m_indices[face*3+0]];
		Vector3 b = m_mesh->m_positions[m_mesh->m_indices[face*3+1]];
		Vector3 c = m_mesh->m_positions[m_mesh->m_indices[face*3+2]];

		return Point3(u*a + v*b + (1.0f-u-v)*c);
    }

    virtual float Area() const
    {
        return m_totalArea;
    }

    Mesh* m_mesh;
    SpatialTree* m_aabbTree;

	float m_totalArea;

	// used for sampling
	std::vector<float> m_cdf;
};

class PlanePrimitive : public Primitive
{
public:
	
	// construct with a plane equation
    PlanePrimitive(Material* mat, const Point3& center, const Vector3& normal, const Vector3& extents=FLT_MAX) 
        : Primitive(mat)
        , m_center(center)
        , m_plane(center, normal)
        , m_extents(extents)
	{
	}

	virtual bool Intersect(const Point3& rayOrigin, const Vector3& rayDir, float& t, Vector3* normal) const
	{
		if (normal)
			*normal = m_plane.GetNormal();

		return IntersectRayPlane(rayOrigin, rayDir, m_plane, t);
	}

private:

	Point3 m_center;
	Vector3 m_extents;

	Plane m_plane;	
};

class DiscPrimitive : public Primitive
{
public:

	DiscPrimitive(Material* mat, const Point3& pos, const Vector3& normal, float radius) 
        : Primitive(mat)
        , m_position(pos)
        , m_plane(pos, normal)
		, m_radius(radius)
	{
		m_objToWorld = TransformFromVector(normal, pos);
	}

	virtual bool Intersect(const Point3& rayOrigin, const Vector3& rayDir, float& t, Vector3* normal) const
	{
		Point3 hit;

		// intersect with plane
		if (IntersectRayPlane(rayOrigin, rayDir, m_plane, t))
		{
			hit = rayOrigin + rayDir*t;
			
			// check point lies inside disc radius
			if (LengthSq(hit-m_position) < m_radius*m_radius)
			{
				if (normal)
					*normal = m_plane.GetNormal();
				return true;
			}
		}

		return false;
	}

	// generate a random uniformly distributed point on the disc
	virtual Point3 Sample(Vector3* n) const
	{		
		if (n)
			*n = m_plane.GetNormal();

		Vec2 sample = UniformSampleDisc() * m_radius;		
		return m_objToWorld * Point3(sample.x, sample.y, 0.0f);
	}

	virtual float Area() const
	{
		return kPi*m_radius*m_radius;
	}

//private:

	//TODO: lots of data duplication here, tidy up
	Matrix44 m_objToWorld;

	Point3 m_position;
	Plane m_plane;

	float m_radius;

};

class ConePrimitive : public Primitive
{
public:

	ConePrimitive(Material* mat, float radius, float height) 
        : Primitive(mat)
        , m_radius(radius)
        , m_height(height)
	{
	}

	virtual bool Intersect(const Point3& rayOrigin, const Vector3& rayDir, float& t, Vector3* normal) const
	{

		float s = (m_radius*m_radius) / (m_height*m_height);
		float q = rayOrigin.z-m_height;

		float a = rayDir.x*rayDir.x + rayDir.y*rayDir.y - s*rayDir.z*rayDir.z;
		float b = 2.0f*(rayDir.x*rayOrigin.x + rayDir.y*rayOrigin.y - s*rayDir.z*q);
		float c = rayOrigin.x*rayOrigin.x + rayOrigin.y*rayOrigin.y - s*q*q;

		float maxT;
		
		if (SolveQuadratic(a, b, c, t, maxT))
		{
			// temp, write something to visualise
			*normal = (rayOrigin + rayDir * t) - Point3(0.0f, 0.0f, 0.0f);
			return true;
		}
		else
		{
			return false;
		}
	}

private:

	float m_radius;
	float m_height;
};

class DistanceFieldPrimitive : public Primitive
{
public:

	// function sig for distance field evaluation
	typedef float (*EvalFunc)(const Point3& p);

	DistanceFieldPrimitive(Material* mat, EvalFunc f, float l=1.0f) 
        : Primitive(mat)
        , m_func(f)
        , m_LipschitzConstant(l)
	{
	}

	virtual bool Intersect(const Point3& rayOrigin, const Vector3& rayDir, float& t, Vector3* normal) const
	{
		// distance to terminate search
		const float kMinDistance = 0.001f;
		const float kMaxDistance = 1000.0f;

		// need to offset slightly otherwise won't detect secondary intersection with self
		Point3 p = rayOrigin + rayDir*kMinDistance;
		t = 0.0f;
		
		float lastd = FLT_MAX;

		while (t < kMaxDistance)
		{
			// calculate distance, step along ray			
			float d = m_func(p) * m_LipschitzConstant;

			if (d < kMinDistance)
				break;

			p += rayDir*d;
			t += d;
		} 

		if (t >= kMaxDistance)
			return false;

		// calculate normal using finite differences
		const float kDelta = 0.001f;

		normal->x = m_func(p+Vector3(kDelta, 0.0f, 0.0f))-m_func(p-Vector3(kDelta, 0.0f, 0.0f));
		normal->y = m_func(p+Vector3(0.0f, kDelta, 0.0f))-m_func(p-Vector3(0.0f, kDelta, 0.0f));
		normal->z = m_func(p+Vector3(0.0f, 0.0f, kDelta))-m_func(p-Vector3(0.0f, 0.0f, kDelta));

		*normal = SafeNormalize(*normal);

		return true;
	}

	// must be set so that |f(x)-f(y)| <= L*|x-y|
	float m_LipschitzConstant;
	
	// field evaluator
	EvalFunc m_func;

};

