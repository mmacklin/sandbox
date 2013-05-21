#pragma once

struct Mesh;

// voxelizes a mesh using a single pass parity algorithm
void Voxelize(const Mesh& mesh, uint32_t width, uint32_t height, uint32_t depth, uint32_t* volume, Vec3 minExtents, Vec3 maxExtents);