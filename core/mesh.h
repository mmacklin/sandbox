#pragma once

#include <vector>

#include "core.h"
#include "maths.h"

struct Mesh
{
    void AddMesh(Mesh& m);

    uint32_t GetNumVertices() const { return m_positions.size(); }
    uint32_t GetNumFaces() const { return m_indices.size() / 3; }

	void DuplicateVertex(uint32_t i);

    void CalculateNormals();
    void Transform(const Matrix44& m);

    void GetBounds(Vector3& minExtents, Vector3& maxExtents) const;

    std::vector<Point3> m_positions;
    std::vector<Vector3> m_normals;
    std::vector<Vector2> m_texcoords[2];
    std::vector<Colour> m_colours;

    std::vector<uint32_t> m_indices;    
};

// create mesh from file
Mesh* ImportMeshFromObj(const char* path);
Mesh* ImportMeshFromPly(const char* path);

// create procedural primitives
Mesh* CreateQuadMesh(float size, float y=0.0f);
Mesh* CreateDiscMesh(float radius, uint32_t segments);

// unwrap the mesh, 
void UvAtlas(Mesh* m);
