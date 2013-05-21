#include "sdf.h"

#include <vector>
#include <float.h>
#include <math.h>

using namespace std;

namespace
{
	inline float Sqr(float x) { return x*x; }
	inline int Clamp(int x, int lower, int upper) { return min(max(lower, x), upper); }

	uint32_t Sample(const uint32_t* image, uint32_t w, uint32_t h, int x, int y)
	{
		return image[Clamp(y, 0, h-1)*w + Clamp(x, 0, w-1)];
	}

	uint32_t Sample(const uint32_t* image, uint32_t w, uint32_t h, uint32_t d, int x, int y, int z)
	{
		return image[Clamp(z, 0, d-1)*w*h + Clamp(y, 0, h-1)*w + Clamp(x, 0, w-1)];
	}

	// returns true if point is on the surface
	bool EdgeDetect(const uint32_t* img, uint32_t w, uint32_t h, int x, int y)
	{
		bool center = Sample(img, w, h, x, y) != 0;

		for (int j=y-1; j <= y+1; ++j)
		{
			for (int i=x-1; i <= x+1; ++i)
			{
				if ((0 != Sample(img, w, h, i, j)) != center)
				{
					return true;
				}
			}
		}
		
		return false;
	}


	// returns true if point is on the surface
	bool EdgeDetect(const uint32_t* img, uint32_t w, uint32_t h, uint32_t d, int x, int y, int z)
	{
		bool center = Sample(img, w, h, d, x, y, z) != 0;

		for (int k=z-1; k <= z+1; ++k)
		{
			for (int j=y-1; j <= y+1; ++j)
			{
				for (int i=x-1; i <= x+1; ++i)
				{
					if ((0 != Sample(img, w, h, d, i, j, k)) != center)
					{
						return true;
					}
				}
			}
		}

		return false;
	}
}

// 2D fast marching method (FMM J. Sethian. A fast marching level set method for monotonically advancing fronts. Proc. Natl. Acad. Sci., 93:1591–1595, 1996.)
namespace 
{
	struct Coord2D
	{
		int i, j;
		
		float d;
		int si, sj;

		bool operator < (const Coord2D& c) const { return d > c.d; }
	};
}

void MakeSDF(const uint32_t* img, uint32_t w, uint32_t h, float* output)
{	
	const float scale = 1.0f / max(w, h);

	std::vector<Coord2D> queue;

	// find surface points
	for (uint32_t y=0; y < h; ++y)
	{
		for (uint32_t x=0; x < w; ++x)
		{
			if (EdgeDetect(img, w, h, x, y))
			{
				Coord2D c = {x, y, 0.0f, x, y};
				queue.push_back(c);
			}

			output[y*w + x] = FLT_MAX;
		}
	}

	std::make_heap(queue.begin(), queue.end());

	while (!queue.empty())
	{
		std::pop_heap(queue.begin(), queue.end());

		Coord2D c = queue.back();
		queue.pop_back();

		// freeze coord if not already frozen
		if (output[c.j*w + c.i] == FLT_MAX)
		{
			output[c.j*w + c.i] = c.d;

			// update neighbours
			int xmin = max(c.i-1, 0), xmax = min(c.i+1, int(w-1));
			int ymin = max(c.j-1, 0), ymax = min(c.j+1, int(h-1));

			for (int y=ymin; y <= ymax; ++y)
			{
				for (int x=xmin; x <= xmax; ++x)
				{
					if (c.i != x || c.j != y)
					{
						int dx = x-c.si;
						int dy = y-c.sj;

						// calculate distance to source coord
						float d = sqrtf(float(dx*dx + dy*dy));

						Coord2D newc = {x, y, d, c.si, c.sj};
						queue.push_back(newc);
						std::push_heap(queue.begin(), queue.end());
					}
				}
			}	
		}
	}

	for (uint32_t y=0; y < h; ++y)
	{
		for (uint32_t x=0; x < w; ++x)
		{
			assert(output[y*w + x] < FLT_MAX);

			// flip sign for interior
			output[y*w + x] *= (img[y*w + x]?-1.0f:1.0f)*scale;
		}
	}
}

// 3D fast marching method (FMM J. Sethian. A fast marching level set method for monotonically advancing fronts. Proc. Natl. Acad. Sci., 93:1591–1595, 1996.)
namespace 
{
	struct Coord3D
	{
		int i, j, k;
		
		float d;
		int si, sj, sk;

		bool operator < (const Coord3D& c) const { return d > c.d; }
	};
}

void MakeSDF(const uint32_t* img, uint32_t w, uint32_t h, uint32_t d, float* output)
{	
	const float scale = 1.0f / max(max(w, h), d);

	std::vector<Coord3D> queue;

	// find surface points
	for (uint32_t z=0; z < d; ++z)
	{
		for (uint32_t y=0; y < h; ++y)
		{
			for (uint32_t x=0; x < w; ++x)
			{
				if (EdgeDetect(img, w, h, d, x, y, z))
				{
					Coord3D c = {x, y, z, 0.0f, x, y, z};
					queue.push_back(c);
				}

				output[z*w*h + y*w + x] = FLT_MAX;
			}
		}
	}

	std::make_heap(queue.begin(), queue.end());

	while (!queue.empty())
	{
		std::pop_heap(queue.begin(), queue.end());

		Coord3D c = queue.back();
		queue.pop_back();

		// freeze coord if not already frozen
		if (output[c.k*w*h + c.j*w + c.i] == FLT_MAX)
		{
			output[c.k*w*h + c.j*w + c.i] = c.d;

			// update neighbours
			int xmin = max(c.i-1, 0), xmax = min(c.i+1, int(w-1));
			int ymin = max(c.j-1, 0), ymax = min(c.j+1, int(h-1));
			int zmin = max(c.k-1, 0), zmax = min(c.k+1, int(d-1));

			for (int z=zmin; z <= zmax; ++z)
			{
				for (int y=ymin; y <= ymax; ++y)
				{
					for (int x=xmin; x <= xmax; ++x)
					{
						if ((c.i != x || c.j != y || c.k != z) && output[z*w*h + y*w + x] == FLT_MAX)
						{
							int dx = x-c.si;
							int dy = y-c.sj;
							int dz = z-c.sk;

							// calculate distance to source coord
							float d = sqrtf(float(dx*dx + dy*dy + dz*dz));

							assert(d > 0.0f);

							Coord3D newc = {x, y, z, d, c.si, c.sj, c.sk};

							queue.push_back(newc);
							std::push_heap(queue.begin(), queue.end());
						}
					}
				}	
			}
		}
	}

	for (uint32_t z=0; z < d; ++z)
	{
		for (uint32_t y=0; y < h; ++y)
		{
			for (uint32_t x=0; x < w; ++x)
			{
				assert(output[z*w*h + y*w + x] < FLT_MAX);

				// flip sign for interior
				output[z*w*h + y*w + x] *= (img[z*w*h + y*w + x]?-1.0f:1.0f)*scale;
			}
		}
	}
}




/*
// Brute-force 2D SDF generation

void FindNeighbour(const uint32_t* image, uint32_t w, uint32_t h, uint32_t cx, uint32_t cy, uint32_t& i, uint32_t& j, float& d)
	{
		float minDistSq=FLT_MAX;

		float fx = float(cx);
		float fy = float(cy);

		for (uint32_t y=0; y < h; ++y)
		{
			for (uint32_t x=0; x < w; ++x)
			{
				if ((x != cx || y != cy) && image[y*w + x])
				{
					float dSq = Sqr(fx-float(x)) + Sqr(fy-float(y));
					if (dSq < minDistSq)
					{
						minDistSq = dSq;
						i = x;
						j = y;
					}
				}
			}	
		}

		d = sqrtf(minDistSq);
	}	


// brute force
void MakeSDF(const uint32_t* img, uint32_t w, uint32_t h, float* output)
{
	// find surface points
	vector<uint32_t> surface(w*h);

	for (uint32_t y=0; y < h; ++y)
	{
		for (uint32_t x=0; x < w; ++x)
		{
			if (EdgeDetect(img, w, h, x, y))
			{
				surface[y*w + x] = 1;
			}
		}
	}

	// brute force search
	for (uint32_t y=0; y < h; ++y)
	{
		for (uint32_t x=0; x < w; ++x)
		{
			uint32_t i, j;
			float d;
			FindNeighbour(&surface[0], w, h, x, y, i, j, d);
		
			// flip sign for pixels inside the shape	
			float sign = (img[y*w + x])?-1.0f:1.0f; 

			output[y*w + x] = d*sign/w;
		}	
	}
}
*/
