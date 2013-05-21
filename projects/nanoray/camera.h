#pragma once

class Camera
{
public:

	Camera() : m_cameraToWorld(Matrix44::kIdentity),
			   m_rasterToWorld(Matrix44::kIdentity)
	{
	}

	Camera(const Matrix44& cameraToWorld,
		 	 float fov,
			 float near,
			 float far,
			 uint32_t width,
			 uint32_t height) : m_cameraToWorld(cameraToWorld)

	{
		Matrix44 rasterToScreen( 2.0f / width, 0.0f, 0.0f, 0.0f,
								 0.0f, -2.0f / height, 0.0f, 0.0f,
								 0.0f,  0.0f, 1.0f, 0.0f,
								-1.0f,  1.0f, 1.0f, 1.0f);

		float f = tanf(fov*0.5f);
		float aspect = float(width) / height;

		Matrix44 screenToCamera(f*aspect, 0.0f, 0.0f, 0.0f,
										0.0f, f, 0.0f, 0.0f, 
										0.0f, 0.0f, -1.0f, 0.0f,
										0.0f, 0.0f, 0.0f, 1.0f);

		m_rasterToWorld = cameraToWorld*screenToCamera*rasterToScreen;
	}	

	void GenerateRay(uint32_t rasterX, uint32_t rasterY, Point3& origin, Vector3& dir) const
	{
		float xoff = Randf(-0.5f, 0.5f);
		float yoff = Randf(-0.5f, 0.5f);

		Point3 p = m_rasterToWorld * Point3(float(rasterX) + 0.5f + xoff, float(rasterY) + 0.5f + yoff, 0.0f);		

        origin = m_cameraToWorld.GetTranslation();
		dir = Normalize(p-m_cameraToWorld.GetTranslation());
	}

	void GenerateRayNoJitter(uint32_t rasterX, uint32_t rasterY, Point3& origin, Vector3& dir) const
	{
		Point3 p = m_rasterToWorld * Point3(float(rasterX) + 0.5f, float(rasterY) + 0.5f, 0.0f);		

        origin = m_cameraToWorld.GetTranslation();
		dir = Normalize(p-m_cameraToWorld.GetTranslation());
	}

	Matrix44 m_cameraToWorld;
	Matrix44 m_rasterToWorld;
};

