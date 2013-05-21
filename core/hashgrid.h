#pragma once

#include "maths.h"

// fixed size grid with wrapping, supports point queries only
template <typename T, int NumBuckets>
class HashGrid
{
public:

	HashGrid(float cellSize=1.0f) : m_cellSize(cellSize)
								  , m_rcpCellSize(1.0f / cellSize)
	{
		memset(m_data, 0, sizeof(m_data));
	}

	~HashGrid()
	{
		Clear();
	}

	void Insert(const Point3& pos, const T& val)
	{
		int idx = GetCell(pos);

		Entry* p = new Entry;
		p->pos = pos;
		p->val = val;
		p->next = m_data[idx];

		m_data[idx] = p;

		assert(GetValue(pos));
	}

	// will call the functor specified for every entry within the query bounds
	template <typename Func>
	void QuerySphere(const Point3& p, float r, Func& f)
	{
		// by requiring cell size to be greater than the largest query radius we only
		// can limit our search to the closest 8 neighboring cells
		assert(r < m_cellSize);

		const int x = int(floorf(p.x * m_rcpCellSize));
		const int y = int(floorf(p.y * m_rcpCellSize));
		const int z = int(floorf(p.z * m_rcpCellSize));

		const float radiusSq = r*r;

		// search the current cell and all neighbors
		for (int i=x-1; i <= x+1; ++i)
		{
			for (int j=y-1; j <= y+1; ++j)
			{
				for (int k=z-1; k <= z+1; ++k)
				{
					int idx = GetCell(i, j, k);

					Entry* head = m_data[idx];
					Entry** prev = &m_data[idx];	// will be updated in case an element is deleted

					for (Entry* iter=head; iter;)
					{
						Entry* next = iter->next;

						// finally, the sphere test
						if (LengthSq(iter->pos-p) <= radiusSq)
						{
							// invoke the functor, return value allows erasing element
							bool erase = f(iter->pos, iter->val);

							if (erase)
							{
								// update previous nodes pointer
								*prev = next;
								delete iter;
							}
							else
							{
								prev = &iter->next;
							}
						}
						else
						{
							prev = &iter->next;
						}
						
						iter = next;
					}
				}
			}
		}
	}

	// iterates over entries
	void QueryAll(T& func)
	{
		for (int i=0; i < NumBuckets; ++i)
		{
			for (Entry* iter=m_data[i]; iter; iter=iter->next)
			{
				func(iter->pos, iter->val);
			}
		}
	}

	void Clear()
	{
		for (int i=0; i < NumBuckets; ++i)
		{
			for (Entry* iter=m_data[i]; iter;)
			{
				Entry* next = iter->next;
				delete iter;
				iter = next;
			}
		}

		memset(m_data, 0, sizeof(m_data));
	}


	const T* GetValue(const Point3& p) const
	{
		const int idx = GetCell(p);

		// iterator over entries looking for given position
		for (Entry* iter=m_data[idx]; iter; iter=iter->next)
		{
			if (iter->pos == p)
			{
				return &iter->val;
			}
		}

		return NULL;
	}

	int GetCell(const Point3& p) const
	{
		return GetCell((int)floorf(p.x*m_rcpCellSize), (int)floorf(p.y*m_rcpCellSize), (int)floorf(p.z*m_rcpCellSize));
	}

	int GetCell(int x, int y, int z) const
	{
		// from RTCD p288
		const int32_t h1 = 0x8da6b343;
		const int32_t h2 = 0xd8163841;
		const int32_t h3 = 0xcb1ab31f;

		int32_t n = h1 * x + h2 * y + h3 * z;
		
		n %= NumBuckets;
		if (n < 0)
			n += NumBuckets;

		return n;
	}

	int GetSize() const
	{
		int s=0;
		for (int i=0; i < NumBuckets; ++i)
		{
			for (Entry* iter=m_data[i]; iter; iter=iter->next)
			{
				++s;
			}
		}

		return s;
	}

	struct Entry
	{
		Point3 pos;
		T val;

		Entry* next;
	};

	float m_cellSize;
	float m_rcpCellSize;	// 1.0f / cellSize;

	// array big enough for xdiv*ydiv*zdiv, each entry points to the head of each cell
	Entry* m_data[NumBuckets];
};

/*
// fixed size grid with wrapping, supports point queries only
template <typename T, int NumBuckets>
class HashGrid2D
{
public:
	
	HashGrid2D(float cellSize=1.0f) : m_cellSize(cellSize)
	, m_rcpCellSize(1.0f / cellSize)
	{
		memset(m_data, 0, sizeof(m_data));

		m_pool = malloc(sizeof(Entry)*MaxItems);
	}
	
	~HashGrid2D()
	{
		Clear();

		free(m_pool);
	}
	
	void Insert(const Vec2& pos, const T& val)
	{
		int idx = GetCell(pos);
		
		Entry* p = &m_pool[m_numItems++];
		
		p->pos = pos;
		p->val = val;
		p->next = m_data[idx];
		
		m_data[idx] = p;
		
		assert(GetValue(pos));
	}
	
	// will call the functor specified for every entry within the query bounds
	template <typename Func>
	void QuerySphere(const Vec2& p, float r, Func& f)
	{
		// by requiring cell size to be greater than the largest query radius we only
		// can limit our search to the closest 8 neighboring cells
		assert(r < m_cellSize);
		
		const int x = int(floorf(p.x * m_rcpCellSize));
		const int y = int(floorf(p.y * m_rcpCellSize));
		
		const float radiusSq = r*r;
		
		// search the current cell and all neighbors
		for (int i=x-1; i <= x+1; ++i)
		{
			for (int j=y-1; j <= y+1; ++j)
			{
					int idx = GetCell(i, j);
					
					Entry* head = m_data[idx];
					Entry** prev = &m_data[idx];	// will be updated in case an element is deleted
					
					for (Entry* iter=head; iter;)
					{
						Entry* next = iter->next;
						
						// finally, the sphere test
						if (LengthSq(iter->pos-p) <= radiusSq)
						{
							// invoke the functor, return value allows erasing element
							f(iter->pos, iter->val);
							
							prev = &iter->next;
						}
						else
						{
							prev = &iter->next;
						}
						
						iter = next;
					}
				}
			}
		}
	}
	
	// iterates over entries
	void QueryAll(T& func)
	{
		for (int i=0; i < NumBuckets; ++i)
		{
			for (Entry* iter=m_data[i]; iter; iter=iter->next)
			{
				func(iter->pos, iter->val);
			}
		}
	}
	
	void Clear()
	{
		m_numItems = 0;
		memset(m_data, 0, sizeof(m_data));
	}
	
	
	const T* GetValue(const Vec2& p) const
	{
		const int idx = GetCell(p);
		
		// iterator over entries looking for given position
		for (Entry* iter=m_data[idx]; iter; iter=iter->next)
		{
			if (iter->pos == p)
			{
				return &iter->val;
			}
		}
		
		return NULL;
	}
	
	int GetCell(const Vec2& p) const
	{
		return GetCell((int)floorf(p.x*m_rcpCellSize), (int)floorf(p.y*m_rcpCellSize));
	}
	
	int GetCell(int x, int y) const
	{
		// from RTCD p288
		const int32_t h1 = 0x8da6b343;
		const int32_t h2 = 0xd8163841;
		
		int32_t n = h1 * x + h2 * y;
		
		n %= NumBuckets;
		if (n < 0)
			n += NumBuckets;
		
		return n;
	}
	
	int GetSize() const
	{
		return m_numItems;
	}
	
	struct Entry
	{
		Vec2 pos;
		T val;
		
		Entry* next;
	};
	
	float m_cellSize;
	float m_rcpCellSize;	// 1.0f / cellSize;
	
	// array big enough for xdiv*ydiv*zdiv, each entry points to the head of each cell
	Entry* m_data[NumBuckets];
	Entry* m_pool[MaxItems];
	uint32_t m_numItems;
};*/