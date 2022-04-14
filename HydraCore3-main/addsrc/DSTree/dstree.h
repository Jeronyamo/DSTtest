#pragma once

#include "dst.h"
#include "../../../RayTracing/CrossRT.h"
//#include "../AddHeaders/CrossRT.h"
#include <iostream>
#include <fstream>
#include <iomanip>


struct simpleAABB {
	LiteMath::float3 min;
	LiteMath::float3 max;

	simpleAABB() : min(), max() {}
	simpleAABB(const LiteMath::float3& new_min, const LiteMath::float3& new_max) : min(new_min), max(new_max) {}
};

struct simpleInstance {
	uint32_t geomID;
	uint32_t firstInstID;
	uint32_t lastInstID;
	LiteMath::float4x4 transf;
};

struct simpleMeshInfo {
	uint32_t firstVertID;
	uint32_t lastVertID;
	uint32_t firstIndID;
	uint32_t lastIndID;
};


class DSTree : public ISceneObject {
	std::vector <LiteMath::float4> vertices;
	std::vector <unsigned> indices;
	std::vector <unsigned> instances;
	std::vector <simpleInstance> instances_info;
	std::vector <simpleMeshInfo> meshes;
	DST dst_nodes;
	size_t instances_size;
	size_t instances_info_size;


/* ================  Functions  ================ */
	friend float getMinPos(const void* index);
	friend float getMaxPos(const void* index);
  /* ========  SAH calculation  ======== */
	const unsigned DSTREE_MAX_POLY = 8u;
	const float EMPTY_COST = 3e3f;
	const float SAH_MAX = 10e10f;
	unsigned current_scene_axis = 3u;

	//float sceneAABB_arr[6u];
	simpleAABB scene_AABB;
	simpleAABB getAABB(unsigned index);
	simpleAABB mergeAABB(unsigned index, simpleAABB aabb);

	uint32_t getTrIndex(uint32_t instTrIndex);
	uint32_t getInstID(uint32_t instTrIndex);
	LiteMath::float3& getVertex(unsigned instance_ind, unsigned* instances_ptr, unsigned geom_vert_ind);
	LiteMath::float3x3& getTrVertices(unsigned instTrIndex);
	LiteMath::float3& getTrVerticesAxis(unsigned instTrIndex);


	float calcSurf(const simpleAABB& aabb) { return 2.0f * ((aabb.max.x - aabb.min.x) +
															(aabb.max.y - aabb.min.y) +
															(aabb.max.z - aabb.min.z));
	}

	float calculateSAH(unsigned& new_axis, unsigned& new_index, unsigned *instances_ptr,
					   simpleAABB& parentAABB, simpleAABB& leftAABB, simpleAABB& rightAABB,
					   unsigned elem_count, float parentSAH);


  /* ========  Tree builder  ======== */
	unsigned buildHeader(bool is_leaf, bool is_carve, unsigned info);
	bool sortTempPlanes(std::vector <unsigned>& target_plane_vec, const LiteMath::float3& minDif, const LiteMath::float3& maxDif,
						float curr_difference, unsigned plane_info, bool dif_consition);
	unsigned buildCarvingNodes(simpleAABB& parent_aabb, simpleAABB& child_aabb, unsigned axis,
								bool is_left_child, bool is_leaf, unsigned triangle_index);
	void dstBuilderRecur(unsigned first, unsigned last, float parentSAH);


  /* ========  Ray traverse  ======== */
	//enum BoolParams { CHECK_LEAF_POLY, IS_CARVING_NODE, IS_DOUBLE_CARVE };
	/*
	struct bool2 {
		bool x = false, y = false;

		bool2() {}
		bool2(bool temp) : x(temp), y(temp) {}
		bool2(bool temp_x, bool temp_y) : x(temp_x), y(temp_y) {}
		bool2(const bool2& temp_b2) : x(temp_b2.x), y(temp_b2.y) {}

		bool& operator[](const unsigned& ind) { if (ind) return y; return x; }
	};

	struct bool3 {
		bool x = false, y = false, z = false;

		bool3() {}
		bool3(bool temp) : x(temp), y(temp), z(temp) {}
		bool3(bool temp_x, bool temp_y, bool temp_z) : x(temp_x), y(temp_y), z(temp_z) {}
		bool3(const bool3& temp_b3) : x(temp_b3.x), y(temp_b3.y), z(temp_b3.z) {}

		bool& operator[](const unsigned& ind) { if (ind > 1u) return z; if (ind) return y; return x; }
	};
*/
	bool traceAABB(LiteMath::float3 Position, LiteMath::float3 Direction, LiteMath::float3 InvDir);
	CRT_Hit traceTriangle(LiteMath::float3 Position, LiteMath::float3 Direction, unsigned trIndex);
	CRT_Hit findHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, bool findAny);

public:
	//void setVertices(const std::vector <LiteMath::float4>& new_vertices) { vertices = new_vertices; }
	//void setIndices(const std::vector <unsigned>& new_indices) { indices = new_indices; }
	//std::vector <unsigned>& getIndices() { return indices; }
	//float* getSceneAABB();
	//std::vector <DSNode>& getDSTNodes() { return dst_nodes.nodes; }


	unsigned AddGeom_Triangles4f(const LiteMath::float4* a_vpos4f, size_t a_vertNumber, const unsigned* a_triIndices, size_t a_indNumber);
	void UpdateGeom_Triangles4f(unsigned a_geomId, const LiteMath::float4* a_vpos4f, size_t a_vertNumber, const unsigned* a_triIndices, size_t a_indNumber);

	CRT_Hit RayQuery_NearestHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar);
	bool RayQuery_AnyHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar);

	uint32_t AddInstance(uint32_t a_geomId, const LiteMath::float4x4& a_matrix);
	void UpdateInstance(uint32_t a_instanceId, const LiteMath::float4x4& a_matrix);
	void ClearGeom();
	void ClearScene();
	void CommitScene();
};