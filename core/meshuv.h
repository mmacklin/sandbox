#pragma once

/*
 *  MeshUv.h
 *  Core
 *
 *  Created by Miles Macklin on 25/04/11.
 *  Copyright 2011 None. All rights reserved.
 *
 */

#include "mesh.h"

#include <set>
#include <vector>
#include <algorithm>

struct Edge;

struct Face
{
	Face()
	{
		e[0] = e[1] = e[2] = 0;
	}
	
	Edge* e[3];
	uint32_t v[3];
	
	Vector3 n;	// normal
	Point3 c;	// centroid
};

struct Edge
{
	Edge(int v0, int v1)
	{
		assert(v0 != v1);
		
		v[0] = std::min(v0, v1);
		v[1] = std::max(v0, v1);
		
		f[0] = f[1] = NULL;
		
		//cout << "Inserting: " << v[0] << ", " << v[1] << endl;
	}
	
	int v[2];
	Face* f[2];
};

bool operator<(const Edge& lhs, const Edge& rhs)
{
	return (lhs.v[0] < rhs.v[0] || (lhs.v[0] == rhs.v[0] && lhs.v[1] < rhs.v[1]));
}

void BindEdgeToFace(Face& f, Edge& e, int edgeIndex)
{
	
	// set face to first edge slot available
	if (!e.f[0])
	{
		e.f[0] = &f;
	}
	else
	{
		// check mesh is 2-manifold
		assert(!e.f[1]);
		e.f[1] = &f;
	}
	
	// set edge to face
	assert(f.e[edgeIndex] == NULL);
	f.e[edgeIndex] = &e;
	
}

struct Chart
{
	Chart(Face* start, const Mesh* m) 
	: mesh(m) 
	, worldArea(0.0f)
	, scale(1.0f)
	, translation(0.0f, 0.0f)
	, minuv(FLT_MAX, FLT_MAX)
	, maxuv(-FLT_MAX, -FLT_MAX)
	{
		AddEdge(start->e[0]);
		AddEdge(start->e[1]);
		AddEdge(start->e[2]);
		AddFace(start);
	}
	
	void AddFace(Face* f)
	{
		const Point3& a = mesh->m_positions[f->v[0]];
		const Point3& b = mesh->m_positions[f->v[1]];
		const Point3& c = mesh->m_positions[f->v[2]];
		
		// update chart world area
		float lengthSqr = LengthSq(Cross(b-a, c-a));
		if (lengthSqr > 0.0f)
		{
			worldArea += sqrt(lengthSqr)*0.5f;
		}
		
		const Vector2& u = mesh->m_texcoords[0][f->v[0]];
		const Vector2& v = mesh->m_texcoords[0][f->v[1]];
		const Vector2& w = mesh->m_texcoords[0][f->v[2]];
		
		// update uv bounds
		minuv = Min(minuv, u);
		minuv = Min(minuv, v);
		minuv = Min(minuv, w);
		
		maxuv = Max(maxuv, u);
		maxuv = Max(maxuv, v);
		maxuv = Max(maxuv, w);
		
		faces.insert(f);
	}
	
	void AddEdge(Edge* e)
	{
		boundary.insert(e);
	}
	
	void RemoveEdge(Edge* e)
	{
		boundary.erase(e);
	}
	
	void RemapVertex(uint32_t oldIndex, uint32_t newIndex)
	{
		// find any references to old vertex and remap to new index
		for (std::set<Face*>::iterator it=faces.begin(), end=faces.end(); it != end; ++it)
		{
			Face* face = *it;
			
			// for each vertex in face
			for (int i=0; i < 3; ++i)
			{
				uint32_t v = face->v[i];
				
				if (v == oldIndex)
				{
					face->v[i] = newIndex;
				}
			}
		}
	}
	
	float worldArea;
	
	float scale;
	Vector2 translation;
	
	Vector2 minuv;
	Vector2 maxuv;
	
	std::set<Face*> faces;
	std::set<Edge*> boundary;
	
	const Mesh* mesh;
};

typedef std::set<Edge> EdgeSet;
typedef std::vector<Face> FaceSet;
typedef std::vector<Chart*> ChartSet;

// unwrap the mesh
void UvAtlas(Mesh* mesh, EdgeSet& edges, FaceSet& faces, ChartSet& charts)
{	
	uint32_t numFaces = mesh->GetNumFaces();
	
	// pre-allocate
	faces.resize(numFaces);
	
	for (uint32_t i=0; i < numFaces; ++i)
	{
		int v0 = mesh->m_indices[i*3+0];
		int v1 = mesh->m_indices[i*3+1];
		int v2 = mesh->m_indices[i*3+2];
		
		Face& f = faces[i];
		
		const Point3& a = mesh->m_positions[v0];
		const Point3& b = mesh->m_positions[v1];
		const Point3& c = mesh->m_positions[v2];
		
		f.n = Normalize(Cross(b-a, c-a));
		f.c = (a+b+c)/3.0f;		
		f.v[0] = v0;
		f.v[1] = v1;
		f.v[2] = v2;
		
		BindEdgeToFace(f, const_cast<Edge&>(*edges.insert(Edge(v0, v1)).first), 0);
		BindEdgeToFace(f, const_cast<Edge&>(*edges.insert(Edge(v1, v2)).first), 1);
		BindEdgeToFace(f, const_cast<Edge&>(*edges.insert(Edge(v2, v0)).first), 2);
	}
	
	// create clusters
	std::set<Face*> unchartedFaces;
	
	for (uint32_t i=0; i < numFaces; ++i)
		unchartedFaces.insert(&faces[i]);
	
	std::vector<Chart*> vertexToChart(mesh->m_positions.size());
	
	while (!unchartedFaces.empty())
	{
		// pick a random start face
		Face* seed = *unchartedFaces.begin();
		unchartedFaces.erase(seed);
		
		// create new chart
		Chart* chart = new Chart(seed, mesh);
		
		std::vector<Edge*> edgeQueue(chart->boundary.begin(), chart->boundary.end());
		
		// try and grow the chart by expanding along boundary edges
		while (!edgeQueue.empty())
		{
			Edge* e = edgeQueue.back();
			edgeQueue.pop_back();
			
			// find the first connected face that isn't in the chart set already
			Face* candidateFace = NULL;
			
			for (int i=0; i < 2; ++i)
			{
				if (chart->faces.find(e->f[i]) == chart->faces.end())
				{
					candidateFace = e->f[i];
					break;
				}
			}
			
			// if it is an open edge then we may find no candidate face
			if (!candidateFace)
				continue;
			
			const float kThreshold = -0.6f;
			
			// calculate error
			float error = Dot(e->f[0]->n, e->f[1]->n);
			
			if (true || error > kThreshold)
			{
				// insert new boundary edges (excluding currently processed edge)
				for (int i=0; i < 3; ++i)
				{
					if (candidateFace->e[i] != e)
					{
						chart->AddEdge(candidateFace->e[i]);
						
						edgeQueue.push_back(candidateFace->e[i]);
					}
				}
				
				// remove old edge from the boundary
				chart->RemoveEdge(e);
				
				// add face to chart
				chart->AddFace(candidateFace);
				
				// remove from uncharted faces
				unchartedFaces.erase(candidateFace);
			}
		}
		
		// duplicate any vertices already referenced from other charts
		for (std::set<Face*>::iterator it=chart->faces.begin(), end=chart->faces.end(); it != end; ++it)
		{
			const Face* face = *it;
			
			// for each vertex in face
			for (int i=0; i < 3; ++i)
			{
				uint32_t v = face->v[i];
				
				// if vertex already referenced then
				if (vertexToChart[v] == NULL)
				{
					vertexToChart[v] = chart;
				}
				else if (vertexToChart[v] != chart)
				{
					// copy vertex, fix up all references
					mesh->DuplicateVertex(v);
										
					chart->RemapVertex(v, mesh->GetNumVertices()-1);
					
					vertexToChart.push_back(chart);
				}
			}
		}
		
		// add chart
		charts.push_back(chart);
	}
}

struct ChartSizeSort
{
	bool operator()(Chart* const & lhs, Chart* const & rhs)
	{
		return lhs->worldArea > rhs->worldArea;
	}
};

template <typename ChartIterator>
void PackChartsStrip(ChartIterator begin, ChartIterator end)
{		
	float totalWorldArea = 0.0f;
	
	// calculate total world surface area
	for (ChartIterator cIt=begin; cIt != end; ++cIt)
	{		
		totalWorldArea += (*cIt)->worldArea;
	}
	
	std::sort(begin, end, ChartSizeSort());
	
	Vector2 packPos(0.0f, 0.0f);
	float bottom = 0.0f;
	
	for (ChartIterator it=begin; it != end; ++it)
	{
		Chart& chart = *(*it);
		
		// normalize chart bounds
		Vector2 minuv = chart.minuv;
		Vector2 maxuv = chart.maxuv;
		
		Vector2 edgeLength = (maxuv-minuv);
		
		float scale = sqrtf(chart.worldArea / (totalWorldArea*edgeLength.x*edgeLength.y));
		
		if (edgeLength.x == 0.0f || edgeLength.y == 0.0f)
			scale = 0.0f;
		
		edgeLength *= scale;
		
		chart.scale = scale;
		chart.translation = packPos;
		
		// move along strip
		packPos.x += edgeLength.x;
		
		bottom = std::max(bottom, packPos.y + edgeLength.y);
		
		if (packPos.x > 1.0f)
			packPos = Vector2(0.0f, bottom);
	}
}

struct PackGrid
{
	PackGrid(uint32_t gridWidth, uint32_t gridHeight) 
		: m_width(0)
		, m_height(0)
		, m_gridWidth(gridWidth)
		, m_gridHeight(gridHeight)
	{
		
		m_grid = new uint8_t[gridWidth*gridHeight];
		memset(m_grid, 0, sizeof(uint8_t)*gridWidth*gridHeight);
	}
	
	~PackGrid()
	{
		delete[] m_grid;
	}
	
	// tries to pack a rect of w,h into the current bounds, returns true on success
	inline bool TryPack(uint32_t w, uint32_t h, uint32_t& x, uint32_t& y) const
	{
		if (w > m_width || h > m_height)
			return false;
		
		uint32_t endx = m_width-w;
		uint32_t endy = m_height-h;
		
		for (y=0; y < endy; ++y)
		{
			for (x=0; x < endx; ++x)
			{					
				if (!Occupied(x, y, w, h))
				{
					return true;
				}
			}
		}
		
		return false;
	}
		
	// return true if rect is already occupied
	inline bool Occupied(uint32_t px, uint32_t py, uint32_t w, uint32_t h) const
	{			   
		uint32_t xend = px+w;
		uint32_t yend = py+h;

		assert(xend <= m_width);
		assert(yend <= m_height);
		
		for (uint32_t y=py; y < yend; ++y)
		{
			const uint8_t* row = &m_grid[y*m_gridWidth];
			
			for (uint32_t x=px; x < xend; ++x)
			{
				if (row[x] != 0)
					return true;
			}
		}
		
		return false;
	}
	
	inline void Fill(uint32_t px, uint32_t py, uint32_t w, uint32_t h)
	{		
		uint32_t xend = px+w;
		uint32_t yend = py+h;

		assert(xend <= m_width);
		assert(yend <= m_height);
		
		for (uint32_t y=py; y < yend; ++y)
		{
			memset(&m_grid[y*m_gridWidth + px], 1, w);
	//		for (uint32_t x=px; x < xend; ++x)
	//		{
	//			m_grid[y*m_gridWidth + x] = 1;
	//		}
		}
	}
	
	void SetSize(uint32_t w, uint32_t h)
	{
		assert(w <= m_gridWidth);
		assert(h <= m_gridHeight);
		
		m_width = w;
		m_height = h;
	}
	
private:
	
	// valid dimensions
	uint32_t m_width;
	uint32_t m_height;
	
	const uint32_t m_gridWidth;
	const uint32_t m_gridHeight;
	
	uint8_t* m_grid;
};

float PackCost(uint32_t mw, uint32_t mh, uint32_t nw, uint32_t nh)
{	
	// calculate the relative increase in map size
	uint32_t oldArea = mw*mh;
	uint32_t newArea = nw*nh;
	
	float aspect = float(nw)/nh;
	
	if (aspect < 1.0f)
		aspect = 1.0f / aspect;
	
	return float(newArea)/oldArea + aspect;
}

template <typename ChartIterator>
void PackChartsBrute(ChartIterator begin, ChartIterator end, uint32_t resolution, float efficiency)
{
	const uint32_t kPadding = 1;
	
	PackGrid* rasterGrid = new PackGrid(resolution, resolution);
	
	float totalWorldArea = 0.0f;
	
	// calculate total world surface area
	for (ChartIterator cIt=begin; cIt != end; ++cIt)
	{		
		totalWorldArea += (*cIt)->worldArea;
	}
	
	std::sort(begin, end, ChartSizeSort());
	
	uint32_t mw = 0;
	uint32_t mh = 0;
	
	for (ChartIterator it=begin; it != end; ++it)
	{
		Chart& chart = *(*it);
		
		// skip degenerate charts
		Vector2 chartBounds = (chart.maxuv-chart.minuv);
		if (chartBounds.x == 0.0f || chartBounds.y == 0.0f)
			continue;		
		
		chart.scale = sqrtf(efficiency*chart.worldArea / (totalWorldArea*chartBounds.x*chartBounds.y));
		
		// scale chart bounds
		chartBounds *= chart.scale;
		
		// search for position in current bounds to hold chart
		uint32_t cw = ceilf(resolution*chartBounds.x) + 2*kPadding;
		uint32_t ch = ceilf(resolution*chartBounds.y) + 2*kPadding;
		
		// pack pos
		uint32_t px, py;
		
		if (!rasterGrid->TryPack(cw, ch, px, py))
		{				
			// failed to pack see if we can expand the map bounds sufficiently, otherwise we will have to give up
			if ((mw + cw) > resolution && (mh + ch) > resolution)
			{
				std::cout << "Failed to fit with the efficiency: " << efficiency << std::endl;
				return;
			}

			// evalulate the cost of expanding horizontally and vertically
			float hcost = PackCost(mw, mh, mw+cw, std::max(mh, ch));
			float vcost = PackCost(mw, mh, std::max(mw, cw), mh+ch);
			
			// take lowest and update pack position and bounds
			if (hcost < vcost)
			{
				px = mw;
				py = 0;
				
				mw += cw;
				mh = std::max(mh, ch);
			}
			else
			{
				px = 0;
				py = mh;
				
				mw = std::max(mw, cw);
				mh += ch;
			}

			// expand bounds
			rasterGrid->SetSize(mw, mh);
		}

		// succeeded, fill bounds
		rasterGrid->Fill(px, py, cw, ch);
	   
		// set chart position
		chart.translation = Vector2((px+kPadding+0.5f)/float(resolution), (py+kPadding+0.5f)/float(resolution));
	}
	
	delete rasterGrid;
}

// generate new chart tex coords based on transforms
template <typename ChartIterator>
void BakeCharts(ChartIterator begin, ChartIterator end, const std::vector<Vector2>& texin, std::vector<Vector2>& texout)
{
	// create new texcoord set
	texout.resize(texin.size());
	
	for (ChartIterator it=begin; it != end; ++it)
	{
		const Chart& chart = *(*it);
		
		for (std::set<Face*>::iterator fIt=chart.faces.begin(), fEnd=chart.faces.end(); fIt != fEnd; ++fIt)
		{
			const Face* face = *fIt;
			
			texout[face->v[0]] = chart.scale*(texin[face->v[0]]-chart.minuv) + chart.translation;
			texout[face->v[1]] = chart.scale*(texin[face->v[1]]-chart.minuv) + chart.translation;
			texout[face->v[2]] = chart.scale*(texin[face->v[2]]-chart.minuv) + chart.translation;
		}
	}
}
