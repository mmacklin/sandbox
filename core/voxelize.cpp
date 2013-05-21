#include "core/aabbtree.h"
#include "core/mesh.h"

void Voxelize(const Mesh& mesh, uint32_t width, uint32_t height, uint32_t depth, uint32_t* volume, Vec3 minExtents, Vec3 maxExtents)
{
	memset(volume, 0, sizeof(uint32_t)*width*height*depth);

	// build an aabb tree of the mesh
	AABBTree tree(&mesh.m_positions[0], mesh.m_positions.size(), &mesh.m_indices[0], mesh.m_indices.size()/3); 

	// parity count method, single pass
	const Vec3 extents(maxExtents-minExtents);
	const Vec3 delta(extents.x/width, extents.y/height, extents.z/depth);
	const Vec3 offset(0.5f/width, 0.5f/height, 0.5f/depth);
	
	const float eps = 0.0000001f*extents.z;

	for (uint32_t x=0; x < width; ++x)
	{
		for (uint32_t y=0; y < height; ++y)
		{
			bool inside = false;

			//float sign = -1.0f;

			Vec3 rayDir = Vec3(0.0f, 0.0f, 1.0f);
			Vec3 rayStart = minExtents + Vec3(x*delta.x + offset.x, y*delta.y + offset.y, -0.0f*extents.z); // z-coord starts somewhat outside bounds 

			uint32_t lastTri = uint32_t(-1);
			//while (z < depth)
			for (;;)
			{
				// calculate ray start
				float t, u, v, w, s;
				uint32_t tri;

				if (tree.TraceRay(Point3(rayStart), rayDir, t, u, v, w, s, tri))
				//if (tree.TraceRaySlow(Point3(rayStart), rayDir, t, u, v, w, s, tri))
				{
					//uint32_t i = mesh.m_indices[tri*3+0];
					//uint32_t j = mesh.m_indices[tri*3+1];
					//uint32_t k = mesh.m_indices[tri*3+2];

					//Vec3 n = SafeNormalize(Cross(mesh.m_positions[j]-mesh.m_positions[i], mesh.m_positions[k]-mesh.m_positions[i]), Vec3(0.0f));

					//if (fabsf(n.z) > 0.01f)
					{

						// calculate cell in which intersection occurred
						const float zpos = rayStart.z + t*rayDir.z;
						const float zhit = (zpos-minExtents.z)/delta.z;
					
						uint32_t z = uint32_t((rayStart.z-minExtents.z)/delta.z);
						uint32_t zend = std::min(uint32_t(zhit), depth-1);

						// must be true for termination
						//assert(zend >= z);
						//if (zend == z && zend < depth-1)
	//						zend++;

						if (inside)
						{
							// march along column setting bits 
							for (uint32_t k=z; k <= zend; ++k)
								volume[k*width*height + y*width + x] = uint32_t(-1);
						}
					
						inside = !inside;
					}

					//if (Sign(s) != sign)
						//inside = !inside;
					//sign = Sign(s);

					/*
					if (Sign(s) == sign)
						printf("Mesh not closed tri: %d lasttri: %d\n", tri, lastTri);
					sign = Sign(s);
					*/

					if (tri == lastTri)
						printf("Error self-intersect\n");
					lastTri = tri;

					rayStart += rayDir*(t+eps);

				}
				else
					break;
			}
		}
	}	
}
