#pragma once

#include "core.h"
#include "mesh.h"

struct Surfel
{
	Surfel() : position(0.0f)
			 , normal(0.0f)
			 , area(0.0f)
			 , colour(0.0f)
			 , emissive(0.0f)
		 	 , irradiance(0.0f)
			 , next(0)
			 , numChildren(0)
			 , tri(0)
			 , depth(0)
	{
	}
	
	Vec3 position;
	Vec3 normal;
	float area;

	Colour colour;		// material colour
	Colour emissive;	// material emission
	Colour irradiance;	// current irradiance

	size_t next;		// for traversal
	size_t numChildren; // number of children in cluster
	size_t tri;			// the tri this surfel represents
	size_t depth;
};

typedef std::vector<Surfel> SurfelArray;

// helper class to manage tree nodes stored in a contiguous linear array
template <typename T>
class TreeNode : public T
{
public:

	TreeNode() : childPtr(NULL)
			   , firstChild(0)
			   , parent(0)
	{

	}

	// returns a surfels children
	const TreeNode& GetChild(size_t i) const
	{
		assert(i < numChildren);
		return (this-firstChild)[i];
	}

	TreeNode& GetChild(size_t i)
	{
		assert(i < numChildren);
		return (this-firstChild)[i];
	}

	TreeNode* childPtr;	// relative index to first child
	size_t firstChild;
	size_t parent;		// parent of the surfel
	//size_t depth;		// just for debug
	size_t subtreeSize; // just for debug
};

// some conveniences
typedef TreeNode<Surfel> SurfelNode;
typedef std::vector<SurfelNode> SurfelNodeArray;

class IlluminationSim
{
public:

	enum DrawMode
	{
		kDrawSurfels	= 1,
		kDrawIrradiance = 2,
		kDrawRadiance	= 3,
		kDrawEmissive	= 4,
		kDrawMesh		= 5
	};

	IlluminationSim(Mesh* m);
	~IlluminationSim();

	// steps the simulation
	void Step();

	// debug visualisation
	void DrawDebug(DrawMode mode, int treeDepth=0);

private:
	
	// adds a mesh into the sim
	void AddMesh(Mesh* m);

	// propogates radiance between surfels in the sim
	void UpdateSurfels();
	
	// hierarchical radiance update
	void UpdateSurfelsFast();
	void UpdateSurfelsFastPreOrder();
	void UpdateSurfelsGPU();

	// apply irradiance values back to the control mesh
	void UpdateMesh();

	// generate clusters
	size_t BuildLightTree(SurfelNode* inputs, size_t numInputs, SurfelNode* clusters, const size_t degree);

	// helper 
	void DrawSurfel(const Surfel& s, const Colour& c, bool solid);

	// surfels in our simulation
	SurfelNodeArray m_surfels;

	SurfelArray m_lightTreePreOrder;

	// source mesh
	Mesh* m_mesh;

	void InitGPUBuffers();
	void UpdateGPUBuffers();

	// shaders
	GLuint m_fragmentShader;
	GLuint m_vertexShader;
	GLuint m_glslProgram;
	
	// fragment shader params
	GLuint m_paramRcpScreen;
	GLuint m_paramTexPosition;
	GLuint m_paramTexNormal;
	GLuint m_paramTexIrradiance;
	GLuint m_paramTexAlbedo;
	GLuint m_paramTexIndices;

	// fbos
	GLuint m_fboIrradiance;

	// textures
	GLsizei m_gpuBufferWidth;
	GLsizei m_gpuBufferHeight;

	GLuint m_texPosition;
	GLuint m_texNormal;
	GLuint m_texIrradiance[2];
	GLuint m_texAlbedo;
	GLuint m_texIndices;
	
	Vec4* m_readBuffer;

	// stats
	float m_statAddMeshTime;
	float m_statStepTime;
	float m_statBakeTime;
	float m_statBuildTreeTime;
	size_t m_statMaxChildrenPerNode;
};