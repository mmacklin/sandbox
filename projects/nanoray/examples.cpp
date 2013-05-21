#include "scene.h"
#include "camera.h"


// scene descs
void InitScene0(Scene& scene, Camera& camera)
{
	DiscPrimitive d(Point3(0.0f, 25.0f, 20.0f), Vector3(0.0f, -1.0f, 0.0f), 7.0f);
	d.m_emission = Colour(20.0f, 20.0f, 20.0f);
	d.m_material = new MatteMaterial(new ConstantTexture(Colour(0.0f, 0.0f, 0.0f)));

	uint32_t di = scene.AddPrimitive(d);

	// attach area light to the disc
	Light l(Light::kArea, 1, di);
	scene.AddLight(l);

	
	ConePrimitive s2(10.0f, 20.0);
	s2.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.4f, 0.3f)));
	scene.AddPrimitive(s2);	
	
	SpherePrimitive s2;
	s2.m_position = Point3(3.5f, 3.0f, 20.0f);
	s2.m_radius = 3.0f;
	s2.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.4f, 0.3f)));
	scene.AddPrimitive(s2);


	
	// floor
	PlanePrimitive p;
	p.m_plane = Plane(0.0f, 1.0f, 0.0f, 0.001f);
	p.m_material = new MatteMaterial(new CheckerboardTexture(2.0f, Colour(0.6f, 0.6f, 0.6f), Colour(0.2f, 0.2f, 0.2f)));
	scene.AddPrimitive(p);

	// back
	PlanePrimitive p4;
	p4.m_plane = Plane(0.0f, 0.0f, -1.0f, 100.0f);
	p4.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.6f, 0.6f)));
	scene.AddPrimitive(p4);

	// roof
	PlanePrimitive p1;
	p1.m_plane = Plane(0.0f, -1.0f, 0.0f, 30.0f);
	p1.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.6f, 0.6f)));
	scene.AddPrimitive(p1);
	
}

float DistanceFunc0(const Point3& p)
{
	const Point3 s(-5.0f, 4.0f, 25.0f);

	return Length(s-p) - 5.0f + Abs(Perlin3D(p.x, p.y, p.z, 1, 1.0f));
}


// scene descs
void InitScene1(Scene& scene, Camera& camera)
{
	DiscPrimitive d(Point3(0.0f, 25.0f, 20.0f), Vector3(0.0f, -1.0f, 0.0f), 7.0f);
	d.m_emission = Colour(20.0f, 20.0f, 20.0f);
	d.m_material = new MatteMaterial(new ConstantTexture(Colour(0.0f, 0.0f, 0.0f)));

	uint32_t di = scene.AddPrimitive(d);

	// attach area light to the disc
	Light l(Light::kArea, 1, di);
	scene.AddLight(l);

	/*
	SpherePrimitive s1;
	s1.m_position = Point3(-4.0f, 5.0f, 25.0f);
	s1.m_radius = 5.0f;
	s1.m_material = new PlasticMaterial(new ConstantTexture(Colour(0.3f, 0.4f, 0.6f)));
	scene.AddPrimitive(s1);
	*/

	DistanceFieldPrimitive f(DistanceFunc0, 0.5f);
	//f.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.4f, 0.3f))); 
	f.m_material = new PlasticMaterial(new ConstantTexture(Colour(0.3f, 0.4f, 0.6f)));
	scene.AddPrimitive(f);

	SpherePrimitive s2;
	s2.m_position = Point3(3.5f, 3.0f, 20.0f);
	s2.m_radius = 3.0f;
	s2.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.4f, 0.3f)));
	scene.AddPrimitive(s2);

	// floor
	PlanePrimitive p;
	p.m_plane = Plane(0.0f, 1.0f, 0.0f, 0.001f);
	p.m_material = new MatteMaterial(new CheckerboardTexture(2.0f, Colour(0.6f, 0.6f, 0.6f), Colour(0.2f, 0.2f, 0.2f)));
	scene.AddPrimitive(p);

	// roof
	PlanePrimitive p1;
	p1.m_plane = Plane(0.0f, -1.0f, 0.0f, 30.0f);
	p1.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.6f, 0.6f)));
	scene.AddPrimitive(p1);

	// left
	PlanePrimitive p2;
	p2.m_plane = Plane(1.0f, 0.0f, 0.0f, 30.0f);
	p2.m_material = new MatteMaterial(new ConstantTexture(Colour(0.3f, 0.6f, 0.4f)));
	scene.AddPrimitive(p2);

	// right
	PlanePrimitive p3;
	p3.m_plane = Plane(-1.0f, 0.0f, 0.0f, 30.0f);
	p3.m_material = new MatteMaterial(new ConstantTexture(Colour(0.7f, 0.3f, 0.3f)));
	scene.AddPrimitive(p3);

	// back
	PlanePrimitive p4;
	p4.m_plane = Plane(0.0f, 0.0f, -1.0f, 50.0f);
	p4.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.6f, 0.6f)));
	scene.AddPrimitive(p4);

}


float DistanceFunc1(const Point3& p)
{
	const Point3 s(-5.0f, 4.0f, 25.0f);
	const Point3 t(4.0f, 4.0f, 25.0f);

	return Length(s-p) - 4.0f - Perlin3D(p.x, p.y, p.z, 1, 1.0f);
}

void InitScene2(Scene& scene, Camera& camera)
{
	DiscPrimitive d(Point3(0.0f, 25.0f, 20.0f), Vector3(0.0f, -1.0f, 0.0f), 7.0f);
	d.m_emission = Colour(20.0f, 20.0f, 20.0f);
	d.m_material = new MatteMaterial(new ConstantTexture(Colour(0.0f, 0.0f, 0.0f)));

	uint32_t di = scene.AddPrimitive(d);

	// attach area light to the disc
	Light l(Light::kArea, 4, di);
	scene.AddLight(l);

	SpherePrimitive s;
	s.m_position = Point3(5.0f, 4.0f, 25.0f);
	s.m_radius = 4.0f;
	s.m_material = new PlasticMaterial(new ConstantTexture(Colour(0.6f, 0.4f, 0.3f)));
	scene.AddPrimitive(s);

	DistanceFieldPrimitive f(DistanceFunc1, 0.5f);
	//f.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.4f, 0.3f))); 
	f.m_material = new MatteMaterial(new ConstantTexture(Colour(0.4f, 0.5f, 0.5f)));
	scene.AddPrimitive(f);

	// floor
	PlanePrimitive p;
	p.m_plane = Plane(0.0f, 1.0f, 0.0f, 0.01f);
	//p.m_material = new PlasticMaterial(new ConstantTexture(Colour(0.6f, 0.6f, 0.6f)));
	p.m_material = new MatteMaterial(new CheckerboardTexture(2.0f, Colour(0.6f, 0.6f, 0.6f), Colour(0.2f, 0.2f, 0.2f)));
	scene.AddPrimitive(p);
}

void InitScene3(Scene& scene, Camera& camera, uint32_t width, uint32_t height, float time)
{
	camera = Camera(TransformMatrix(Rotation(0.0f, 0.0f, DegToRad(-5.0f)), Point3(0.0f, 8.0f, 0.0f)), DegToRad(55.0f), 1.0f, 10000.0f, width, height);	

	const float kSunThetaStart = kPi/2.0f;
	const float kSunThetaEnd = -kPi/2.0f;
	const float kSunPhiStart = k2Pi;
	const float kSunPhiEnd = kPi/1.5f;

	const float t = time / 10.0f;

	scene.SetSkyParams(Abs(Lerp(kSunThetaStart, kSunThetaEnd, SmoothStep(0.0f, 1.0f, t))), Lerp(kSunPhiStart, kSunPhiEnd, t), 2.0f);

	SpherePrimitive s;
	s.m_position = Point3(8.0f, 6.0f, 25.0f);
	s.m_radius = 6.0f;
	s.m_material = new PlasticMaterial(new ConstantTexture(Colour(0.6f, 0.4f, 0.3f)));
	//s.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.6f, 0.6f)));
	scene.AddPrimitive(s);

	// floor
	PlanePrimitive p;
	p.m_plane = Plane(0.0f, 1.0f, 0.0f, 0.1f);
	p.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.6f, 0.6f)));
	scene.AddPrimitive(p);
}

// a simple terrain
float DistanceFunc2(const Point3& p)
{	
	const float kScale = 0.005f;

	return p.y - 50.0f*Perlin2D(kScale*p.x, kScale*p.z, 10, 0.5f);
}

void InitScene4(Scene& scene, Camera& camera)
{


	DistanceFieldPrimitive f(DistanceFunc2, 0.1f);
	//f.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.4f, 0.3f))); 
	f.m_material = new MatteMaterial(new ConstantTexture(Colour(0.4f, 0.5f, 0.5f)));
	scene.AddPrimitive(f);
}

extern float gTime;

float SuperEllipticBlend(float d1, float r1, float d2, float r2)
{
	return Min(d1,d2)*(1.0f-powf(Max(0.0f, (1.0f - d1/r1)), 2.0f)-powf(Max(0.0f, (1.0f - d2/r2)), 2.0f));
}

float EuclidBlend(float x, float n)
{
	return (Sqrt(x*x+n)-x)*0.5f;
}

float SonyBlend(float d1, float r1, float d2, float r2)
{
	if (d1 < r1 && d2 < r2)
	{
		float m1 = Clamp((r1-d1)/r1, 0.0f, 1.0f);
		float m2 = Clamp((r2-d2)/r2, 0.0f, 1.0f);

		return powf(m1*m2, 4.0f)*1.2f;
	}
	
	return 0.0f;
}


// metaballs
float DistanceFunc3(const Point3& p)
{	
	const Point3 c1(0.0f, 2.0f + gTime, 25.0f);
	const Point3 c2(0.0f, 6.0f, 25.0f);

	float d1 = (Length(p-c1)-2.0f);
	float d2 = p.y;// - Abs(Perlin2D(p.x*0.1, p.z*0.1, 1, 1.0f));

	float r = Perlin2D(p.x,p.y,1, 1.0f)*0.5f;

	float rd = Sqrt((p.x-c2.x)*(p.x-c2.x)+(p.y-c2.y)*(p.y-c2.y));
	float rq = Max(rd-4.0f,0.0f) - Min(rd-2.0f+r, 0.0f);
	float dd = Max(0.0f, Abs(p.z-c2.z)-1.0f);
	float d3 = Sqrt(dd*dd + rq*rq);
	
	return Min(d1, d3) - EuclidBlend(Abs(d1-d3), 2.0f);
}

void InitScene5(Scene& scene, Camera& camera)
{
	DiscPrimitive d(Point3(0.0f, 25.0f, 20.0f), Vector3(0.0f, -1.0f, 0.0f), 7.0f);
	d.m_emission = Colour(20.0f, 20.0f, 20.0f);
	d.m_material = new MatteMaterial(new ConstantTexture(Colour(0.0f, 0.0f, 0.0f)));

	uint32_t di = scene.AddPrimitive(d);

	// attach area light to the disc
	Light l(Light::kArea, 4, di);
	scene.AddLight(l);

	//scene.SetSkyParams(kPi/3.3f, kPi*2.0f, 2.0f);

	/*
	SpherePrimitive s;
	s.m_position = Point3(0.0f, 6.0f, 25.0f);
	s.m_radius = 6.0f;
	s.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.4f, 0.3f)));
	//s.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.6f, 0.6f)));
	scene.AddPrimitive(s);
	*/
	// floor
	PlanePrimitive p;
	p.m_plane = Plane(0.0f, 1.0f, 0.0f, 0.1f);
	p.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.6f, 0.6f)));
	scene.AddPrimitive(p);

	DistanceFieldPrimitive f(DistanceFunc3, 0.5f);
	//f.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.4f, 0.3f))); 
	f.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.4f, 0.3f)));
	scene.AddPrimitive(f);
}


class FlameVolume : public Volume
{
public:

	FlameVolume(const Point3& p) : m_position(p) {}

	Point3 m_position;

	virtual void Eval(const Point3& p, const Vector3& dir,  float maxT, Colour& transmittance, Colour& emission, Colour& inScatter) const
	{
		const float kStepSize = 0.1f;
		const float kRadius = 5.0f;

		// attenuate emissive term by distance to center
		float closestT, furthestT;
		float visibility = 1.0f;

		if (IntersectRaySphere(m_position, kRadius, p, dir, closestT, furthestT))
		{
			float t = Min(closestT, maxT);
			maxT = Min(furthestT, maxT);

			while (t < maxT)
			{
				Point3 s = p + t*dir;

				float distance = LengthSq(s-m_position);
				float nd = distance / (kRadius*kRadius);
				float d = Max(0.0f, 1.0f - nd);

				// calculate emission 
				float heat = powf(Max(0.0f, Perlin3D(s.x, s.y, s.z, 1, 2.0f)), 6.0f);
				emission += Lerp(Colour(0.8, 0.05, 0.0f, 1.0f), Colour(8.9f, 8.7f, 5.3f, 1.0f), heat)*0.3f*d;
				
				visibility -= Max(0.0f, (nd-0.8f)*5.0f)*Abs(Perlin3D(s.x, s.y, s.z, 1, 2.0f)*5.0f)*0.05f;

				t += kStepSize;
			}
		}

		visibility = Max(0.0f, visibility);
		transmittance = Colour(visibility, visibility, visibility);
	}
};

void InitScene6(Scene& scene, Camera& camera)
{
	
	DiscPrimitive d(Point3(0.0f, 25.0f, 20.0f), Vector3(0.0f, -1.0f, 0.0f), 7.0f);
	d.m_emission = Colour(20.0f, 20.0f, 20.0f);
	d.m_material = new MatteMaterial(new ConstantTexture(Colour(0.0f, 0.0f, 0.0f)));

	uint32_t di = scene.AddPrimitive(d);

	// attach area light to the disc
	Light l(Light::kArea, 1, di);
	scene.AddLight(l);

	SpherePrimitive s2;
	s2.m_position = Point3(3.5f, 3.0f, 20.0f);
	s2.m_radius = 3.0f;
	s2.m_material = new PlasticMaterial(new ConstantTexture(Colour(0.6f, 0.4f, 0.3f)));
	scene.AddPrimitive(s2);

	FlameVolume* f = new FlameVolume(Point3(-3.5, 6.0f, 20.0f));
	scene.AddVolume(f);

	// floor
	PlanePrimitive p;
	p.m_plane = Plane(0.0f, 1.0f, 0.0f, 0.001f);
	p.m_material = new MatteMaterial(new CheckerboardTexture(2.0f, Colour(0.6f, 0.6f, 0.6f), Colour(0.2f, 0.2f, 0.2f)));
	scene.AddPrimitive(p);

	// back
	PlanePrimitive p4;
	p4.m_plane = Plane(0.0f, 0.0f, -1.0f, 50.0f);
	p4.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.6f, 0.6f)));
	scene.AddPrimitive(p4);

	// roof
	PlanePrimitive p1;
	p1.m_plane = Plane(0.0f, -1.0f, 0.0f, 30.0f);
	p1.m_material = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.6f, 0.6f)));
	scene.AddPrimitive(p1);

}