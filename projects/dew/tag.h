#pragma once

#include "core/maths.h"
#include <vector>

struct Vertex
{
	Vertex() {}
	Vertex(Point3 p, Vec3 n, Vec3 c=Vec3(1.0f)) : position(p), normal(n), color(c) {}

	Point3 position;
	Vec3 normal;
	Vec3 color;
};

struct Brush
{
	virtual void Eval(float t, std::vector<Vertex>& verts, float w, float h)
	{

	}
};

struct SquareBrush: public Brush
{
	virtual void Eval(float t, std::vector<Vertex>& verts, float w, float h)
	{
		const Vertex shape[] = 
		{
			Vertex(Point3(-w,  h, 0.0f), Vec3(-1.0f, 0.0f, 0.0f)),
			Vertex(Point3(-w,  -h, 0.0f), Vec3(-1.0f, 0.0f, 0.0f)),
	
			Vertex(Point3( -w, -h, 0.0f), Vec3( 0.0f, -1.0f, 0.0f)),
			Vertex(Point3( w,  -h, 0.0f), Vec3( 0.0f, -1.0f, 0.0f)),
	
			Vertex(Point3( w, -h, 0.0f), Vec3( 1.0f, 0.0f, 0.0f)),
			Vertex(Point3( w,  h, 0.0f), Vec3( 1.0f, 0.0f, 0.0f)),
		
			Vertex(Point3( w,  h, 0.0f), Vec3( 0.0f, 1.0f, 0.0f)),
			Vertex(Point3( -w, h, 0.0f), Vec3( 0.0f, 1.0f, 0.0f))
		};

		verts.assign(shape, shape+8);
	}
};

static Vec3 FromHex(uint32_t rgba)
{
	float r = ((rgba>>24)&0xff)/255.0f;
	float g = ((rgba>>16)&0xff)/255.0f;
	float b = ((rgba>>8)&0xff)/255.0f;

	float gamma = 2.0f;
	float scale = 2.0f;

	return Vec3(powf(r, gamma), powf(g, gamma), powf(b, gamma))*scale;
}

static Vec3 colors1[4] =  {

	FromHex(0x93E46300),
	FromHex(0xCCE2BE00),
	FromHex(0x70A53E00),
	FromHex(0xFF701A00),
 };

static Vec3 colors2[4] =  {

	FromHex(0xFFF7A900),
	FromHex(0xF98D0400),
	FromHex(0xFE9B1E00),
	FromHex(0xF44E9000),
 };

static Vec3 colors3[4] =  {

	FromHex(0xDEE8FC00),
	FromHex(0x4C55FF00),
	FromHex(0x93E46300),
	FromHex(0x4856C300),
 };

struct SquareColorBrush: public Brush
{
	SquareColorBrush(Vec3* scheme) : colors(scheme) {}

	virtual void Eval(float t, std::vector<Vertex>& verts, float w, float h)
	{
		const Vertex shape[] = 
		{
			Vertex(Point3(-w,  h, 0.0f), Vec3(-1.0f, 0.0f, 0.0f), colors[0]),
			Vertex(Point3(-w,  -h, 0.0f), Vec3(-1.0f, 0.0f, 0.0f),colors[0]),
	
			Vertex(Point3( -w, -h, 0.0f), Vec3( 0.0f, -1.0f, 0.0f), colors[1]),
			Vertex(Point3( w,  -h, 0.0f), Vec3( 0.0f, -1.0f, 0.0f), colors[1]),
	
			Vertex(Point3( w, -h, 0.0f), Vec3( 1.0f, 0.0f, 0.0f), colors[2]),
			Vertex(Point3( w,  h, 0.0f), Vec3( 1.0f, 0.0f, 0.0f), colors[2]),
		
			Vertex(Point3( w,  h, 0.0f), Vec3( 0.0f, 1.0f, 0.0f), colors[3]),
			Vertex(Point3( -w, h, 0.0f), Vec3( 0.0f, 1.0f, 0.0f), colors[3])
		};

		verts.assign(shape, shape+8);
	}

	Vec3* colors;
};

struct TriangleBrush: public Brush
{
	virtual void Eval(float t, std::vector<Vertex>& verts, float w, float h)
	{
		/*
		Point3 s[] = { Point3(0.0f, -h*0.5f, 0.0f),
					   Point3(w, -h, 0.0f),
					   Point3(0.0f, h*0.5f, 0.0f),
					   Point3(-w, -h, 0.0f) };
					   */

		Point3 s[] = { Point3(0.0f, h*0.5f, 0.0f),
					   Point3(-w, -h, 0.0f),				   
					   Point3(w, -h, 0.0f) };


		verts.resize(0);

		const uint32_t numPoints = sizeof(s)/sizeof(Point3);

		// output edges
		for (uint32_t i=0; i < numPoints; ++i)
		{		
			uint32_t first = i;
			uint32_t next = (i+1)%numPoints;

			Vec3 n = Normalize(Cross(s[next]-s[first], Vec3(0.0f, 0.0f, 1.0f)));
			Vertex v0(s[first], n);
			Vertex v1(s[next], n);

			verts.push_back(v0);
			verts.push_back(v1);
		}
	}
};

struct CircleBrush : public Brush
{
	virtual void Eval(float t, std::vector<Vertex>& verts, float w, float h)
	{
		const uint32_t k = 8;
		const float kinc = -k2Pi/k;

		verts.resize(0);

		for (uint32_t i=0; i <= 8; ++i)
		{
			float x = sinf(i*kinc);
			float y = cosf(i*kinc);
			Vertex v(Point3(x*w, y*h, 0.0f), Vec3(x, y, 0.0f));
			verts.push_back(v);
		}
	}
};


struct Tag
{
	Tag(float smoothing, float width, float height, float velscale, Brush* b) : basis(Matrix44::kIdentity), smoothing(smoothing), width(width), height(height), velscale(velscale), velocity(0.0f), draw(false), shape(b)
	{
		samples.reserve(4096);
		vertices.reserve(100000);
		indices.reserve(100000);
	}

	void SetBrush(Brush* b) { shape = b; }

	void Start()
	{
		OutputCap(true);	

		draw = true;
	}

	void Stop()
	{
		OutputCap(false);

		draw = false;
	}

	void OutputCap(bool flip)
	{
		// draw cap
		shape->Eval(1.0f, brush, width, height);

		Point3 center(0.0f);

		// transform verts and create faces
		for (size_t i=0; i < brush.size(); ++i)		
		{
			center += Vec3(brush[i].position);
		}
		
		center /= brush.size();
	
		uint32_t i0 = vertices.size();
		
		float dir = flip?-1.0f:1.0f;
		Vec3 n = dir*Vec3(basis.GetCol(2));

		Vertex c(basis*center, n);
		vertices.push_back(c);

		// transform verts and create faces
		for (size_t i=0; i < brush.size(); ++i)
		{
			float s = GetCurrentSize();

			// transform position and normal to world space
			Vertex v(basis*(brush[i].position*s), n);

			vertices.push_back(v);
			
			if (i > 0)
			{
				if (!flip)
				{
					indices.push_back(i0);
					indices.push_back(i0+i);
					indices.push_back(i0+i+1);
				}
				else
				{
					indices.push_back(i0+i+1);
					indices.push_back(i0+i);
					indices.push_back(i0);
				}
			}	
		}
	}

	void PushSample(float t, Matrix44 m)
	{
		// evaluate brush
		shape->Eval(t, brush, width, height);

		size_t startIndex = vertices.size();

		Point3 prevPos = samples.empty()?m.GetTranslation():samples.back();
		Point3 curPos = m.GetTranslation();

		// low-pass filter position
		Point3 p = Lerp(prevPos, curPos, 1.0f-smoothing);
		m.SetTranslation(p);

		samples.push_back(m.GetTranslation());

		// need at least 4 points to construct valid tangents
		if (samples.size() < 4)
			return;

		// the point we are going to output
		size_t c = samples.size()-3;

		// calculate the tangents for the two samples using central differencing
		Vec3 tc = Normalize(samples[c+1]-samples[c-1]);
		Vec3 td = Normalize(samples[c+2]-samples[c]);
		float a = acosf(Dot(tc, td));

		// calculate moving average of velocity
		velocity = 0.8f*Length(samples[c+2]-samples[c])+0.2f*velocity;

		// save rotations
		rotations.push_back(a);

		if (fabsf(a) > 0.001f)
		{
			// use the parallel transport method to move the reference frame along the curve
			Vec3 n = Normalize(Cross(tc, td));

			if (samples.size() == 4)
				basis = TransformFromVector(Normalize(tc));
		
			// 'transport' the basis forward
			basis = RotationMatrix(a, n)*basis;
			
			m = basis;
			m.SetTranslation(samples[c]);
			basis = m;
		}
	
		if (!draw)
			return;
		
		// transform verts and create faces
		for (size_t i=0; i < brush.size(); ++i)
		{
			float s = GetCurrentSize();

			// transform position and normal to world space
			Vertex v(m*(brush[i].position*s), m*brush[i].normal, brush[i].color);

			vertices.push_back(v);
		}

		if (startIndex != 0)
		{
			size_t b = brush.size();

			for (size_t i=0; i < b; ++i)
			{
				size_t curIndex = startIndex + i;
				size_t nextIndex = startIndex + (i+1)%b; 

				indices.push_back(curIndex);
				indices.push_back(curIndex-b);
				indices.push_back(nextIndex-b);
				
				indices.push_back(nextIndex-b);
				indices.push_back(nextIndex);
				indices.push_back(curIndex);			
			}	
		}
	}

	float GetCurrentSize()
	{
		// scale size inversely proportionally to velocity
		//float s = Max(0.0f, 1.0f-velscale*avgSpeed);//Min(1.0f, 1.0f / (velscale*avgSpeed));
		float s = Lerp(1.0f, GetAvgSpeed(), velscale);//Min(1.0f, (1.0f-velscale)*GetAvgSpeed());		
		//float s = Max(1.0f, 1.0f-velscale*GetAvgAngularSpeed());		
		//float s = 1.0f;
		return s;
	}

	// returns brush size based on current velocity
	float GetAvgSpeed()
	{
		return velocity;

		/*
		if (velscale == 0.0f)
			return 1.0f;

		const int numSamples = 10;

		int start = Max(int(0), int(samples.size())-numSamples);
		int end = Min(size_t(start+numSamples), samples.size());

		float avgSpeed = 0.0f;

		for (int i=start+1; i < end; ++i)
		{
			avgSpeed += Length(samples[i] - samples[samples.size()-1]);
		}

		avgSpeed /= numSamples;
		return avgSpeed;
		*/
	}

	void Draw()
	{
		if (vertices.empty())
			return;

		// draw the tag
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, sizeof(Vertex), &vertices[0].position);
		glEnableClientState(GL_NORMAL_ARRAY);
		glNormalPointer(GL_FLOAT, sizeof(Vertex), &vertices[0].normal);
		glEnableClientState(GL_COLOR_ARRAY);
		glColorPointer(3, GL_FLOAT, sizeof(Vertex), &vertices[0].color);
		
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);

	}

	void ExportToObj(const char* path)
	{
		FILE* f = fopen(path, "w");

		if (!f)
			return;

		fprintf(f, "# %d positions\n", int(vertices.size()));

		for (uint32_t i=0; i < vertices.size(); ++i)
		{
			Point3 p = vertices[i].position;
			fprintf(f, "v %f %f %f\n", p.x, p.y, p.z);
		}

		fprintf(f, "# %d normals\n", int(vertices.size()));

		for (uint32_t i=0; i < vertices.size(); ++i)
		{
			Vec3 n = vertices[i].normal;
			fprintf(f, "vn %f %f %f\n", n.x, n.y, n.z);
		}

		fprintf(f, "# %d faces\n", int(indices.size()/3));

		for (uint32_t t=0; t < indices.size(); t+=3)
		{
			// obj is 1 based
			uint32_t i = indices[t+0]+1;
			uint32_t j = indices[t+1]+1;
			uint32_t k = indices[t+2]+1;

			fprintf(f, "f %d//%d %d//%d %d//%d\n", i, i, j, j, k, k); 
		}	

		fclose(f);
	}

	void GetBounds(Point3& lower, Point3& upper)
	{
		Point3 l(FLT_MAX);
		Point3 u(-FLT_MAX);

		for (uint32_t i=0; i < samples.size(); ++i)
		{
			l = Min(l, samples[i]);
			u = Max(u, samples[i]);	
		}	

		lower = l;
		upper = u;
	}

	void Clear()
	{
		samples.resize(0);
		brush.resize(0);
		vertices.resize(0);
		indices.resize(0);
	}

	Matrix44 basis;

	std::vector<float> rotations;
	std::vector<Point3> samples;
	std::vector<Vertex> brush;
	std::vector<Vertex> vertices;
	std::vector<uint32_t> indices;

	float smoothing;
	float width;
	float height;		
	float velscale;
	float velocity;
	bool draw;	
	
	Brush* shape;
};


