#include "illuminationsim.h"
#include "profiler.h"
#include "core.h"
#include "glsl.h"
#include "material.h"

#include <limits>

using namespace std;

// helper function to create a floating point texture
GLuint CreateGPUBuffer(size_t width, size_t height, GLuint format=GL_RGBA16F)
{
	// create a new texture name
	GLuint texID;
	glGenTextures (1, &texID);
	// bind the texture name to a texture target
	glBindTexture(GL_TEXTURE_2D,texID);
	// turn off filtering and set proper wrap mode 
	// (obligatory for float textures atm)
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	// set texenv to replace instead of the default modulate
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	
	// and allocate graphics memory
	GlVerify(glTexImage2D(GL_TEXTURE_2D, 0, format, GLsizei(width), GLsizei(height), 0, GL_RGBA, GL_FLOAT, 0));

	return texID;
}

IlluminationSim::IlluminationSim(Mesh* m)
{
	// initialize shaders
	m_fragmentShader = GlslCreateShader("illuminate.fp", GL_FRAGMENT_SHADER);
	m_vertexShader = GlslCreateShader("illuminate.vp", GL_VERTEX_SHADER);

	m_glslProgram = glCreateProgram();
	GlVerify(glAttachShader(m_glslProgram, m_fragmentShader));
	GlVerify(glAttachShader(m_glslProgram, m_vertexShader));
	GlVerify(glLinkProgram(m_glslProgram));

	// get fragment shader parameters
	m_paramRcpScreen = glGetUniformLocation(m_glslProgram, "rcpScreen");
	m_paramTexPosition = glGetUniformLocation(m_glslProgram, "texPosition");
	m_paramTexNormal = glGetUniformLocation(m_glslProgram, "texNormal");
	m_paramTexIrradiance = glGetUniformLocation(m_glslProgram, "texIrradiance");
	m_paramTexAlbedo = glGetUniformLocation(m_glslProgram, "texAlbedo");
	m_paramTexIndices = glGetUniformLocation(m_glslProgram, "texIndices");
	
	GlslPrintProgramLog(m_glslProgram);

	// initialize fbos
	glGenFramebuffersEXT(1, &m_fboIrradiance); 

	// create textures
	m_gpuBufferWidth = 1024;
	m_gpuBufferHeight = 64;

	// allocate buffer to read back into
	m_readBuffer = new Vec4[m_gpuBufferWidth*m_gpuBufferHeight];

	//create textures
	m_texPosition = CreateGPUBuffer(m_gpuBufferWidth, m_gpuBufferHeight, GL_RGBA32F);
	m_texNormal = CreateGPUBuffer(m_gpuBufferWidth, m_gpuBufferHeight);
	m_texIrradiance[0] = CreateGPUBuffer(m_gpuBufferWidth, m_gpuBufferHeight);
	m_texIrradiance[1] = CreateGPUBuffer(m_gpuBufferWidth, m_gpuBufferHeight);
	m_texIndices = CreateGPUBuffer(m_gpuBufferWidth, m_gpuBufferHeight);
	m_texAlbedo = CreateGPUBuffer(m_gpuBufferWidth, m_gpuBufferHeight);

	if (m)
		AddMesh(m);
}

IlluminationSim::~IlluminationSim()
{
	delete[] m_readBuffer;

	glDeleteProgram(m_glslProgram);
	glDeleteFramebuffersEXT(1, &m_fboIrradiance);

	glDeleteTextures(1, &m_texPosition);
	glDeleteTextures(1, &m_texNormal);
	glDeleteTextures(1, &m_texIndices);
	glDeleteTextures(2, m_texIrradiance);

	GlCheckErrors("Cleaning up");
}

float balance = 0.0f;

// recursive function to output a tree of surfels in pre-order
// and set next pointers accordingly for non-recursive traversal
size_t PreorderTree(SurfelNode& s, SurfelArray& output, size_t treeSize)
{
	size_t currentIndex = output.size();
	size_t subtreeSize = 0;

	// output the node
	output.push_back(s);
	output.back().colour.a = output.back().emissive.r;

	size_t smin = INT_MAX;
	size_t smax = 0;

	// calculate size of this nodes sub-tree
	for (size_t i=0; i < s.numChildren; ++i)
	{
		size_t st = PreorderTree(s.GetChild(i), output, treeSize);
		
		smin = Min(smin, st);
		smax = Max(smax, st);

		subtreeSize += st;
	}
	
	// calculate the range
	if (s.numChildren)
		balance += (smax-smin)/float(subtreeSize);

	size_t offset = subtreeSize+1;

	// if next index would be valid then store it
	if ((currentIndex + offset) < treeSize)
		output[currentIndex].next = offset;

	// sanity checking
	if (s.numChildren)
		assert(s.next == 0 || s.next > 1);
	
	if (s.numChildren == 0)
		assert(s.next == 0 || s.next == 1);

	return subtreeSize+1;	
}

void IlluminationSim::AddMesh(Mesh* m)
{
	{

	assert(m);
	m_mesh = m;

	Mesh::Vertex* vertices = m->GetVertices();
	
	size_t triIndex = 0;

	// loop over material subsets creating surfels
	for (size_t s=0; s < m->GetNumSubsets(); ++s)
	{
		Mesh::Subset* subset = m->GetSubset(int(s));

		uint16* indices = &m->GetIndices()[subset->m_firstIndex];
		size_t numtris = subset->m_numIndices/3;

		// loop over tris calculating surfel information
		for (size_t i=0; i < numtris; i++)
		{
			const Mesh::Vertex& v1 = vertices[indices[i*3+0]];
			const Mesh::Vertex& v2 = vertices[indices[i*3+1]];
			const Mesh::Vertex& v3 = vertices[indices[i*3+2]];

			// calculate tri normal
			Vec3 e1 = v2.pos-v1.pos;
			Vec3 e2 = v3.pos-v1.pos;
			
			Vec3 n = Cross(e1, e2);

			// throw away degenerates
			if (Length(n) > 0.000001f)
			{
				SurfelNode surfel;
				surfel.area = Length(n)*0.5f;						 
				surfel.position = (Vec3(v1.pos) + Vec3(v2.pos) + Vec3(v3.pos))/3.0f;	
				surfel.normal = Normalize(n);
				surfel.irradiance = 0.0f;		
				surfel.emissive = subset->m_material->m_emission*40;//*0.01f;
				surfel.emissive.a = 1.0f;
				surfel.colour = subset->m_material->m_diffuse*0.8f;
				
				//if (subset->m_material->m_diffuse.r == 1.0f && subset->m_material->m_diffuse.b == 0.0f)
					//surfel.colour = Colour(0.3f, 0.4f, 0.6f, 1.0f);
				surfel.colour += Colour(0.1f, 0.1f, 0.1f, 0.0f);

				if (surfel.position.y < 1.0f)
					surfel.emissive = Colour(0.0f, 0.0f, 0.0f, 0.0f);

				surfel.tri = triIndex;

				m_surfels.push_back(surfel);
			}

			++triIndex;
		}
	}

	const size_t kDegree = 5;
	const size_t maxNodes = (kDegree*m_surfels.size())/(kDegree-1);
	const size_t leafNodes = m_surfels.size();

	assert(maxNodes > m_surfels.size());

	// iteratively build the tree bottom up
	size_t currentCount = leafNodes;
	size_t totalCount = currentCount;

	m_surfels.resize(maxNodes);

	SurfelNode* cursor = &m_surfels[0];

	size_t depth = 0;

	while (currentCount > 1)
	{
		size_t c = BuildLightTree(cursor, currentCount, cursor+currentCount, kDegree);

		// DEBUG: set the depth of each surfel
		for (size_t i=0; i < currentCount; ++i)
		{
			cursor[i].depth = depth;

			//if (depth > 0)
			//	assert(cursor[i].numChildren > 0);
		}
		++depth;
		// DEBUG END

		// c should always be less than current count, if its not something went horribly wrong
		assert(c < currentCount);

		// move cursor forward
		cursor += currentCount;
		// update the total number of nodes
		totalCount += c;
		// update the number in the working set
		currentCount = c;
	}
	for (size_t i=0; i < currentCount; ++i)
	{
		cursor[i].depth = depth;
	}

	// final pass through the surfels to calculate relative child index as we copy and reallocate
	// the surfel array we need to use indices instead of pointers
	for (size_t i=0; i < totalCount; ++i)
	{
		SurfelNode& s = m_surfels[i];
		if (s.numChildren)
			s.firstChild = &s - s.childPtr;
	}

	m_surfels.resize(totalCount);
		
	// output in pre-order
	PreorderTree(m_surfels.back(), m_lightTreePreOrder, m_surfels.size());

	assert(m_lightTreePreOrder.size() == m_surfels.size());
	reverse(m_surfels.begin(), m_surfels.end());

	InitGPUBuffers();
	}

	cout << "Added " << m_surfels.size() << " surfels in " << m_statAddMeshTime*1000.0f << "ms" << endl;
	cout << "Built tree in " << m_statBuildTreeTime*1000.0f << "ms" << endl;

	// debug
//	size_t w = WalkTree(m_lightTree.front());
//	assert(w == totalCount);

}

// sorts surfels by parent
bool ParentSort(const SurfelNode& left, const SurfelNode& right)
{
	return left.parent < right.parent;
}

// take a set of inputs and build a set of clusters using Lloyd's algorithm
// it re-arranges the inputs so that parents can just store an index and child count
// returns the number of clusters created with appropriate child indices and counts
size_t IlluminationSim::BuildLightTree(SurfelNode* inputs, size_t numInputs, SurfelNode* clusters, const size_t degree)
{
	assert(inputs);
	assert(inputs < clusters);
	assert(inputs+numInputs <= clusters);	

	const size_t kIterations = 4;	// how many iterations of k-means clustering to perform
	const size_t kDegree = degree;  // aim for a n-ary tree, this is not a guarantee!

	// starting number of clusters
	size_t numClusters = max(size_t(1), numInputs/kDegree);

	// seeded using n/degree random surfels from the inputs array
	for (size_t i=0; i < numClusters; ++i)
	{
		clusters[i] = inputs[i*kDegree];
		
		// reset the hierarchy info
		clusters[i].firstChild = NULL;
		clusters[i].numChildren = 0;
		clusters[i].parent = 0;
		clusters[i].tri = -1;
	}
	
	for (int w=0; w < kIterations; ++w)
	{
		m_statMaxChildrenPerNode = 0;

		// assign each surfel to the  closest cluster
		for (size_t i=0; i < numInputs; ++i)
		{
			SurfelNode& s = inputs[i];

			size_t closestIndex=0;
			float closestError=FLT_MAX;

			// find the closest cluster
			for (size_t n=0; n < numClusters; ++n)
			{
				const SurfelNode& c = clusters[n];

				float distSq = LengthSq(c.position-s.position);
				float angle = 1.0f-Dot(c.normal, s.normal);

				// TODO: integrate UV distance into error metric to support deformable objects
				
				float r = abs(c.emissive.r-s.emissive.r);

				// error is computed as a combination of distance and angle delta				
				float error = distSq + (distSq*angle*2.0f) + (distSq*r*10000000.0);
				
				if (error < closestError)
				{
					closestIndex = n;
					closestError = error;
				}
			}

			assert(closestIndex < numClusters);
			s.parent = closestIndex;
		}

		// zero cluster positions
		for (size_t i=0; i < numClusters; ++i)
		{
			clusters[i].position = Vec3(0.0f);
			clusters[i].normal = Vec3(0.0f);
			clusters[i].area = 0.0f;
			clusters[i].colour = Colour(0.0f);
			clusters[i].emissive = Colour(0.0f);
			clusters[i].irradiance = Colour(0.0f);
			clusters[i].numChildren = 0;
		}

		// update cluster positions
		for (size_t i=0; i < numInputs; ++i)
		{
			const SurfelNode& s = inputs[i];

			assert(s.parent < numClusters);
			SurfelNode& p = clusters[s.parent];

			p.position += s.position;
			p.normal += s.normal*s.area;
			p.area += s.area;
			p.colour += s.colour;
			p.emissive += s.emissive*s.area;
		
			p.numChildren++;
		}
		
		// average
		for (size_t i=0; i < numClusters; ++i)
		{
			SurfelNode& c = clusters[i];
			
			if (c.numChildren == 0)
			{
				c.area = 1.0f;
				continue;
			}
		
/*
			// isolated cluster
			while (c.numChildren == 0)
			{
				cout << "removing isolated cluster" << endl;

				// swap one from the back and lower total number of clusters
				c = clusters[--numClusters];
			}
			assert(c.numChildren > 0);
*/	

			float rcpChildren = 1.0f / c.numChildren;
			float rcpArea = 1.0f / c.area;

			// average values, note we don't average area but sum it instead
			c.position *= rcpChildren;
			c.normal = SafeNormalize(c.normal);
			c.colour *= rcpChildren;
			c.emissive *= rcpArea;
		}
	}
	
	// when we're done sort the children and update the parent's child pointers
	sort(inputs, inputs+numInputs, ParentSort);

	SurfelNode* lastParent=NULL;
	
	// run over them updating the cluster indices
	for (size_t i=0; i < numInputs; ++i)
	{
		SurfelNode& s = inputs[i];
		SurfelNode& c = clusters[s.parent];
		
		if (lastParent == NULL || lastParent != &c)
		{
			// calculate relative index to the first child
			c.childPtr = &s;
				
			// sanity check
			assert(c.numChildren > 0);
			//assert(&c - &s >= c.numChildren);
		}

		lastParent = &c;
	}

	return numClusters;
}


// steps the simulation
void IlluminationSim::Step()
{
	//cout << "Added " << m_surfels.size() << " surfels in " << m_statAddMeshTime*1000.0f << "ms" << endl;
	//cout << "Built tree in " << m_statBuildTreeTime*1000.0f << "ms" << endl;

	// update sim
	//UpdateSurfels();
//	UpdateSurfelsFast();
//	cout << "UpdateSurfelsFast: " << m_statStepTime*1000.0f << endl;

	//UpdateSurfelsFastPreOrder();
	//cout << "UpdateSurfelsFastPreOrder: " << m_statStepTime*1000.0f << endl;

	UpdateGPUBuffers();
	UpdateSurfelsGPU();
	//cout << "UpdateSurfelsGPU() took: " << m_statStepTime*1000.0f << endl;

	// bake to mesh
	UpdateMesh();

	//cout << "BakeTime: " << m_statBakeTime*1000.0f << endl << endl;
}


// propagates radiance between surfels in the sim
void IlluminationSim::UpdateSurfels()
{
	/*
	// profiler
	ScopedTimer timer(m_statStepTime);

	SurfelArray newSurfels(m_surfels);

	// big old O(N^2)
	for (size_t i=0; i < newSurfels.size(); ++i)
	{
		Surfel& r = newSurfels[i];
		r.irradiance = Colour(0.0f, 0.0f, 0.0f, 0.0f);

		// process every other surfel as an emitter
		for (size_t j=0; j < m_surfels.size(); ++j)
		{
			const Surfel& e = m_surfels[j];

			// don't process ourself
			if (i != j)
			{
				Vec3 v(e.position-r.position);				
				Vec3 dir = Normalize(v);

				float ev = Dot(dir,e.normal);
				float rv = Clamp(Dot(dir,r.normal), 0.0f, 1.0f);

				float ff = fabsf((e.area*ev*rv)/(kPi*Dot(v,v) + e.area));

				// add glow sources
				if (e.emissive.r > 0.0f)
				{
					r.irradiance += ff*e.emissive;
					continue;
				}

				// add bounced light
				if (ev < 0.0f)
					r.irradiance += ff*e.irradiance*e.colour;				
				else
					r.irradiance -= ff*e.irradiance;
			}
		}
		r.irradiance.r = Max(0.0f, r.irradiance.r);
		r.irradiance.g = Max(0.0f, r.irradiance.g);
		r.irradiance.b = Max(0.0f, r.irradiance.b);
		r.irradiance.a = 1.0f;
	}

	// swap the new ones in
	m_surfels.swap(newSurfels);	
	*/
}

const float kDistanceThreshold = 4.0f;
/*
void GatherIrradiance(const Surfel& e, Surfel& r)
{
	Vec3 v(e.position-r.position);

	float distSq = LengthSq(v);

	// if we're far enough to the surfel or it's a leaf then use it
	if (e.numChildren == 0 || distSq > kDistanceThreshold*e.area/kPi)
	{
		// stops us processing ourselves
		if (distSq == 0.0f)
			return;

		Vec3 dir(v / sqrtf(distSq));

		float ev = Dot(dir,e.normal);
		float rv = Clamp(Dot(dir,r.normal), 0.0f, 1.0f);

		float ff = fabsf((e.area*ev*rv)/(kPi*distSq + e.area));

		// add glow sources
		if (e.emissive.r > 0.0f)
		{
			r.irradiance += ff*e.emissive;			
		}

		// add bounced light
		if (ev < 0.0f)
			r.irradiance += ff*e.irradiance*e.colour;				
		else
 			r.irradiance -= ff*e.irradiance;

		//assert(_finite(ff));
		//assert(!_isnan(ff));
	}
	else
	{
		for (size_t i=0; i < e.numChildren; ++i)
			GatherIrradiance(e.GetChild(i), r);
	}
}
*/

float CalculateClippedArea(Vec3 v, Vec3 nr, Vec3 ne, float a)
{
	float r2 = abs(a) / kPi;
	float r = sqrt(r2);
	float h = max(Dot(v, nr), 0.0f);
	float d = r;

	//float c = sqrt(r2-d*d);
	float theta = 2.0*acos(d/r);
	float stheta = sin(theta);

	float cseg = 0.5*r2*(theta-stheta);
	return a-cseg;	
}

void GatherIrradiancePreOrder(const Surfel& e, Surfel& r)
{
	const Surfel* current = &e;

	while (current)
	{
		Vec3 v(current->position-r.position);

		float distSq = LengthSq(v) + 1e-16f;

		// if we're far enough to the surfel or it's a leaf then use it
		if (current->numChildren == 0 || distSq > 4.0f*current->area/kPi)
		{
			// stops us processing ourselves
			if (distSq > 0.0f)
			{
				Vec3 dir(v / sqrtf(distSq));

				float ev = Dot(dir,current->normal);
				float rv = Clamp(Dot(dir,r.normal), 0.0f, 1.0f);

				float ff = fabsf((current->area*Clamp(-ev, 0.0f, 1.0f)*rv)/(kPi*distSq + current->area));

				// add glow sources
				//if (current->emissive.r > 0.0f)
				//{
				//	r.irradiance += ff*current->emissive;			
				//}
				float ca = CalculateClippedArea(v, r.normal, current->normal, current->area);

				r.irradiance += ff*current->colour*current->colour.a;			

				// add bounced light
				if (ev < 0.0f)
					r.irradiance += ff*current->irradiance*current->colour;				
				else
				{
					if (ff != 0.0f)
						cout << "yeah" << endl;

					float ffa = fabsf((current->area*ev*rv)/(kPi*distSq + current->area));
					r.irradiance -= ffa*current->irradiance*3.0f;
				}
			}
	
			//assert(_finite(ff));
			//assert(!_isnan(ff));

			if (current->next)
				current += current->next;
			else
				current = NULL;
		}
		else if (current->numChildren)
			++current;
		else
		{
			assert(0);
		}
	}
}

void IlluminationSim::UpdateSurfelsFast()
{
	/*
	// profiler
	ScopedTimer timer(m_statStepTime);

	SurfelArray newSurfels(m_lightTree);

	for (size_t i=0; i < newSurfels.size(); ++i)
	{
		Surfel& r = newSurfels[i];
		r.irradiance = Colour(0.0f, 0.0f, 0.0f, 0.0f);

		// if a parent node, then average children's irradiance, considerably 
		// faster than treating internal nodes the same as leaves
		if (r.numChildren)
		{
			for (size_t i=0; i < r.numChildren; ++i)
			{
				Surfel& e = (&r-r.firstChild)[i];
				r.irradiance += e.irradiance*e.area;
			}

			r.irradiance /= r.area;
		}
		// otherwise gather it from the light tree
		else
		{
			GatherIrradiance(m_lightTree.back(), r);
		}

		// process every other surfel as an emitter
		//r.irradiance.r = Max(0.0f, r.irradiance.r);
		//r.irradiance.g = Max(0.0f, r.irradiance.g);
		//r.irradiance.b = Max(0.0f, r.irradiance.b);
		//r.irradiance.a = 1.0f;
	}

	// swap the new ones in
	m_lightTree.swap(newSurfels);	
	*/
}

void IlluminationSim::UpdateSurfelsFastPreOrder()
{
	// profiler
	SurfelArray newSurfels(m_lightTreePreOrder);

	for (SurfelArray::iterator iter=newSurfels.begin(); iter != newSurfels.end(); ++iter)
	{
		Surfel& r = *iter;
		r.irradiance = Colour(0.0f, 0.0f, 0.0f, 0.0f);

		/*
		//GatherIrradiance(m_lightTree.back(), r);
		if (r.numChildren)
		{
			Surfel* c=&r+1;

			for (size_t i=0; i < r.numChildren; ++i)
			{
				r.irradiance += c->irradiance*c->area;
				c += c->next;
			}

			r.irradiance /= r.area;
		}
		else
		{*/
			GatherIrradiancePreOrder(m_lightTreePreOrder.front(), r);
		//}

		// process every other surfel as an emitter
		r.irradiance.r = Max(0.0f, r.irradiance.r);
		r.irradiance.g = Max(0.0f, r.irradiance.g);
		r.irradiance.b = Max(0.0f, r.irradiance.b);
		r.irradiance.a = Max(0.0f, r.irradiance.a);
		//r.irradiance.a = 1.0f;
	}

	// swap the new ones in
	m_lightTreePreOrder.swap(newSurfels);	
}

void IlluminationSim::InitGPUBuffers()
{
	Vec4* positions = new Vec4[m_gpuBufferWidth*m_gpuBufferHeight];
	Vec4* normals = new Vec4[m_gpuBufferWidth*m_gpuBufferHeight];
	Vec4* albedo = new Vec4[m_gpuBufferWidth*m_gpuBufferHeight];
	Vec4* indices = new Vec4[m_gpuBufferWidth*m_gpuBufferHeight];

	memset(positions, (uint)std::numeric_limits<float>::quiet_NaN(), sizeof(Vec4)*m_gpuBufferWidth*m_gpuBufferHeight);
	memset(normals, 0, sizeof(Vec4)*m_gpuBufferWidth*m_gpuBufferHeight);
	memset(albedo, 0, sizeof(Vec4)*m_gpuBufferWidth*m_gpuBufferHeight);
	memset(indices, 0, sizeof(Vec4)*m_gpuBufferWidth*m_gpuBufferHeight);

	GlVerify(glBindTexture(GL_TEXTURE_2D, m_texIrradiance[0]));
	GlVerify(glTexSubImage2D(GL_TEXTURE_2D,0,0,0,m_gpuBufferWidth,m_gpuBufferHeight, GL_RGBA, GL_FLOAT, indices));

	GlVerify(glBindTexture(GL_TEXTURE_2D, m_texIrradiance[1]));
	GlVerify(glTexSubImage2D(GL_TEXTURE_2D,0,0,0,m_gpuBufferWidth,m_gpuBufferHeight, GL_RGBA, GL_FLOAT, indices));

	for (size_t i=0; i < m_lightTreePreOrder.size(); ++i)
	{
		const Surfel& s = m_lightTreePreOrder[i];

		positions[i].x = s.position.x;
		positions[i].y = s.position.y;
		positions[i].z = s.position.z;	
		positions[i].w = s.area;

		// if this is a parent node then make the area negative
		if (s.numChildren)
			positions[i].w *= -1.0f;

		// pack normal and area in one buffer
		normals[i].x = s.normal.x;
		normals[i].y = s.normal.y;
		normals[i].z = s.normal.z;
		normals[i].w = s.numChildren;

		albedo[i].x = s.colour.r;
		albedo[i].y = s.colour.g;
		albedo[i].z = s.colour.b;
		albedo[i].w = s.colour.a;

		// next in pre-order traversal
		if (s.next)
		{
			indices[i].x = ((i + s.next) % m_gpuBufferWidth) / float(m_gpuBufferWidth);
			indices[i].y = ((i + s.next) / m_gpuBufferWidth) / float(m_gpuBufferHeight);
		}
		else
		{
			indices[i].x = 1.0f;
			indices[i].y = 1.0f;
		}

		// next child
		if (s.numChildren)
		{
			indices[i].z = ((i+1) % m_gpuBufferWidth) / float(m_gpuBufferWidth);
			indices[i].w = ((i+1) / m_gpuBufferWidth) / float(m_gpuBufferHeight);
		}
		else
		{
			indices[i].z = 1.0f;
			indices[i].w = 1.0f;
		}
	}

	GlVerify(glBindTexture(GL_TEXTURE_2D, m_texPosition));
	GlVerify(glTexSubImage2D(GL_TEXTURE_2D,0,0,0,m_gpuBufferWidth,m_gpuBufferHeight, GL_RGBA, GL_FLOAT, positions));

	GlVerify(glBindTexture(GL_TEXTURE_2D, m_texNormal));
	GlVerify(glTexSubImage2D(GL_TEXTURE_2D,0,0,0,m_gpuBufferWidth,m_gpuBufferHeight, GL_RGBA, GL_FLOAT, normals));

	GlVerify(glBindTexture(GL_TEXTURE_2D, m_texAlbedo));
	GlVerify(glTexSubImage2D(GL_TEXTURE_2D,0,0,0,m_gpuBufferWidth,m_gpuBufferHeight, GL_RGBA, GL_FLOAT, albedo));

	GlVerify(glBindTexture(GL_TEXTURE_2D, m_texIndices));
	GlVerify(glTexSubImage2D(GL_TEXTURE_2D,0,0,0,m_gpuBufferWidth,m_gpuBufferHeight, GL_RGBA, GL_FLOAT, indices));

	delete[] positions;
	delete[] albedo;
	delete[] normals;
	delete[] indices;

	GlCheckErrors("Creating source textures");
}

void IlluminationSim::UpdateGPUBuffers()
{
	/*
	Vec4* positions = new Vec4[m_gpuBufferWidth*m_gpuBufferHeight];

	for (size_t i=0; i < m_lightTreePreOrder.size(); ++i)
	{
		Surfel& s = m_lightTreePreOrder[i];
	
		float offset = 0.0f;

		if (s.emissive.r > 0.0f)
		{
			offset = sin(fmod(float(GetSeconds()), 2.0f*kPi));
			s.position.y -= offset;
		}

		positions[i].x = s.position.x;
		positions[i].y = s.position.y;
		positions[i].z = s.position.z;
		positions[i].w = s.area;

		// if this is a parent node then make the area negative
		if (s.numChildren)
			positions[i].w *= -1.0f;
	}

	glBindTexture(GL_TEXTURE_2D, m_texPosition);
	glTexSubImage2D(GL_TEXTURE_2D,0,0,0,m_gpuBufferWidth,m_gpuBufferHeight, GL_RGBA, GL_FLOAT, positions);

	GlCheckErrors("Creating source textures");

	delete[] positions;
	*/
}

void IlluminationSim::UpdateSurfelsGPU()
{
	// bind offscreen buffer 
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_fboIrradiance);

	// bind texture
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_texIrradiance[0], 0);	

	// debug
	GlFramebufferStatus();

	// setup projection
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, m_gpuBufferWidth, 0.0, m_gpuBufferHeight);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glViewport(0, 0, m_gpuBufferWidth, m_gpuBufferHeight);

	// draw full screen quad
	glPolygonMode(GL_FRONT,GL_FILL);
	glUseProgram(m_glslProgram);	

	// set params
	GlVerify(glActiveTexture(GL_TEXTURE0));
	GlVerify(glBindTexture(GL_TEXTURE_2D, m_texPosition));
	GlVerify(glUniform1i(m_paramTexPosition,0)); // texunit 0

	GlVerify(glActiveTexture(GL_TEXTURE1));
	GlVerify(glBindTexture(GL_TEXTURE_2D, m_texNormal));
	GlVerify(glUniform1i(m_paramTexNormal,1)); // texunit 1

	GlVerify(glActiveTexture(GL_TEXTURE2));
	GlVerify(glBindTexture(GL_TEXTURE_2D, m_texIrradiance[1]));
	GlVerify(glUniform1i(m_paramTexIrradiance,2)); // texunit 2

	GlVerify(glActiveTexture(GL_TEXTURE3));
	GlVerify(glBindTexture(GL_TEXTURE_2D, m_texAlbedo));
	GlVerify(glUniform1i(m_paramTexAlbedo,3)); // texunit 3

	GlVerify(glActiveTexture(GL_TEXTURE4));
	GlVerify(glBindTexture(GL_TEXTURE_2D, m_texIndices));
	GlVerify(glUniform1i(m_paramTexIndices,4)); // texunit 4

	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0); 
	glVertex2f(0.0, 0.0);
	glTexCoord2f(1.0, 0.0); 
	glVertex2f(m_gpuBufferWidth, 0.0);
	glTexCoord2f(1.0, 1.0);
	glVertex2f(m_gpuBufferWidth, m_gpuBufferHeight);
	glTexCoord2f(0.0, 1.0); 
	glVertex2f(0.0, m_gpuBufferHeight);
	glEnd();

	// read back results
	GlVerify(glReadBuffer(GL_COLOR_ATTACHMENT0_EXT));
	GlVerify(glReadPixels(0,0,m_gpuBufferWidth,m_gpuBufferHeight, GL_RGBA, GL_FLOAT, m_readBuffer));

	// unbind stuff
	GlVerify(glActiveTexture(GL_TEXTURE2));
	GlVerify(glBindTexture(GL_TEXTURE_2D, 0));	
	GlVerify(glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0));	

	swap(m_texIrradiance[0], m_texIrradiance[1]);
/*
	// assert results are correct
	for (size_t i=0; i < m_lightTreePreOrder.size(); ++i)
	{
		Surfel& s = m_lightTreePreOrder[i];
		
		if (Length((Vec3)s.irradiance - (Vec3)m_readBuffer[i]) > 0.001f)
		{
			cout << "wtf" << endl;
		}
	}
*/
}

// apply irradiance values to source mesh
void IlluminationSim::UpdateMesh()
{
	Mesh* m = m_mesh;

	Mesh::Vertex* vertices = m->GetVertices();
	uint16* indices = m->GetIndices();

	// zero vert colours
	for (size_t i=0; i < m->GetNumVertices(); ++i)
		vertices[i].col = Colour(0.0f);

	// loop over surfels adding irradiance contribution back to verts
	for (size_t i=0; i < m_lightTreePreOrder.size(); ++i)
	{	
		const Surfel& s = m_lightTreePreOrder[i];

		if (s.numChildren != 0 || s.tri == size_t(-1))
			continue;

		size_t t = s.tri*3;

		Mesh::Vertex& v1 = vertices[indices[t+0]];	
		Mesh::Vertex& v2 = vertices[indices[t+1]];
		Mesh::Vertex& v3 = vertices[indices[t+2]];
	
		const Colour& irradiance = ((Colour*)m_readBuffer)[i];
		//const Colour& irradiance = s.irradiance;
		
		Colour c(s.colour*irradiance*kInvPi + s.emissive);
		c.a = 1.0f;
		
		v1.col += c;
		v2.col += c;
		v3.col += c;
	}
	
	// average colours
	for (size_t i=0; i < m->GetNumVertices(); ++i)
	{
		vertices[i].col = vertices[i].col / vertices[i].col.a;
		vertices[i].col = LinearToSrgb(vertices[i].col);
	}
}

void IlluminationSim::DrawSurfel(const Surfel& s, const Colour& c, bool solid)
{
	// draw discs
	const Vec3 up(0.0f, 1.0, 0.0f);
	const Vec3 right(1.0f, 0.0f, 0.0f);

	glBegin(s.emissive.r > 0.0f?GL_TRIANGLE_FAN:GL_LINE_STRIP);

	float kNumSteps = 12;
	float kStepSize = 2.0f*kPi / kNumSteps;
	float kRadius = sqrtf(s.area/kPi);

	// calculate plane vectors
	Vec3 u = kRadius*Normalize(Cross(s.normal, fabsf(Dot(s.normal, up))>0.95f?right:up));
	Vec3 v = kRadius*Normalize(Cross(s.normal, u));

	glColor3fv(LinearToSrgb(c));

	if (solid)
		glVertex3fv(s.position);

	for (size_t j=0; j < kNumSteps+1; ++j)
	{
		float a = j*kStepSize;
		glVertex3fv(s.position+Cos(a)*u+Sin(a)*v);
	}

	glEnd();
}

void DrawText(int x, int y, const char *string);

// debug visualisation
void IlluminationSim::DrawDebug(DrawMode mode, int treeDepth)
{
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);	 
	glUseProgram(0);		

	const SurfelArray& surfels = m_lightTreePreOrder;
	//const SurfelNodeArray& surfels = m_surfels;

	glDisable(GL_LIGHTING);

	if (mode == kDrawSurfels)
	{
		// draw light tree
		for (size_t i=0; i < surfels.size(); ++i)
		{
			const Surfel& s = surfels[i];
			
			if (s.depth == treeDepth)
				DrawSurfel(s, s.colour, false);
		}

		const float kNormalScale=5.0f;

		// draw normals
		glBegin(GL_LINES);

		for (size_t i=0; i < surfels.size(); ++i)
		{
			const Surfel& s = surfels[i];

			if (s.depth == treeDepth)
			{
				glColor3fv(s.colour);
				glVertex3fv(s.position);
				glVertex3fv(s.position+s.normal*kNormalScale);
			}
		}

		glEnd();
	}
	else
	{	/*
		if (mode == kDrawSurfels)
		{
			for (size_t i=0; i < m_surfels.size(); ++i)
			{
				const Surfel& s = m_surfels[i];

				DrawSurfel(s, s.colour, false);
			}

			const float kNormalScale=5.0f;

			// draw normals
			glBegin(GL_LINES);

			for (size_t i=0; i < m_surfels.size(); ++i)
			{
				const Surfel& s = m_surfels[i];

				glColor3fv(s.colour);
				glVertex3fv(s.position);
				glVertex3fv(s.position+s.normal*kNormalScale);
			}

			glEnd();
		}
		*/
	}


	glEnable(GL_DEPTH_TEST);

	// draw surfels
	for (size_t i=0; i < surfels.size(); ++i)
	{
		const Surfel& s = surfels[i];

		if (s.numChildren != 0)
			continue;

		switch (mode)
		{
		case kDrawIrradiance:
			DrawSurfel(s, SrgbToLinear(s.irradiance), true);
			break;
		case kDrawEmissive:
			DrawSurfel(s, SrgbToLinear(s.emissive), true);
			break;
		case kDrawRadiance:
			DrawSurfel(s, SrgbToLinear(s.colour*s.irradiance), true);
			break;
		};
	}

	if (mode == kDrawMesh)
	{
		// set up lighting state
		glDisable(GL_LIGHTING);

		// Create light components
		glShadeModel (GL_SMOOTH);

		m_mesh->Render();

	}
}
