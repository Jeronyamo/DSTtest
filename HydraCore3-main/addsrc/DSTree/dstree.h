#pragma once

#include "dst.h"
#include "../../../RayTracing/CrossRT.h"
//#include "../AddHeaders/CrossRT.h"
#include <iostream>
#include <fstream>
#include <iomanip>

struct simpleDSTinfo {
	bool is_leaf;
	bool is_single_carve;
	bool is_doulble_carve;
	bool is_n1_pos;
	bool is_n2_pos;
	int plane1;
	int plane2;
};
std::vector <simpleDSTinfo> carve_info;
struct simpleAABB {
	LiteMath::float3 min;
	LiteMath::float3 max;

	simpleAABB() : min(), max() {}
	simpleAABB(const LiteMath::float3& new_min, const LiteMath::float3& new_max) : min(new_min), max(new_max) {}
	simpleAABB(const simpleAABB &temp) : min(temp.min), max(temp.max) {}
};

struct simpleInstance {
	uint32_t geomID;
	//uint32_t firstInstID;
	//uint32_t lastInstID;
	LiteMath::float4x4 transf;
	simpleAABB instAABB;
};

struct simpleMeshInfo {
	uint32_t firstVertID;
	uint32_t lastVertID;
	uint32_t firstIndID;
	uint32_t lastIndID;
	uint32_t DSTreeOffset;
	simpleAABB meshAABB;
};


class DSTree : public ISceneObject {
	const float AABBeps = 1e-4f;
	std::vector <LiteMath::float4> vertices;
	std::vector <unsigned> indices;
	std::vector <unsigned> indices_sorted;
	std::vector <simpleMeshInfo> meshes;
	std::vector <simpleInstance> instances_info;

	std::vector <DSNode> upper_tree; //leaves contain instances' indices
	std::vector <DSNode> lower_tree; //stores arrays of meshes' DST
	const LiteMath::float4* tempVertices; //temp pointer of vertices added in AddGeom_Triangles4f
	const unsigned* tempIndices; //temp pointer of indices added in AddGeom_Triangles4f
	const unsigned* tempIndicesSorted; //temp pointer of sorted indices added in AddGeom_Triangles4f => vertices[ 3 * indices[sorted_indices] + 0/1/2]
	DSNode* dst_ptr; //temp pointer for traversal; points to current instance's DST in lower_tree
	//unsigned* instances_ptr;//temp pointer for traversal; points to current instance's instances array
	std::vector <unsigned> instances; //temp array for builder; stores indices while sorting
	std::vector <DSNode> dst_nodes; //temp array for builder; builder output



/* ================  Comparator functions  ================ */
	float getMinPos(const void* index);
	float getMaxPos(const void* index);
	float getMinAABBpos(const void* index);
	float getMaxAABBpos(const void* index);
	friend int cmpElemPos(const void* ind1, const void* ind2);
	friend int cmpElemPosAABB(const void* ind1, const void* ind2);


  /* ========  SAH calculation  ======== */
	const float EMPTY_COST = 3e3f;
	const float SAH_MAX = 10e10f;
	const unsigned DSTREE_MAX_POLY = 8u;
	unsigned current_scene_axis = 3u;
	unsigned current_instance = 0u;

	simpleAABB scene_AABB;
	simpleAABB getAABB(unsigned index);
	simpleAABB mergeAABB(unsigned index, simpleAABB aabb);
	simpleAABB getInstAABB(unsigned index);
	simpleAABB mergeInstAABB(unsigned index, simpleAABB aabb);

	//uint32_t getInstID(uint32_t instTrIndex);
	//uint32_t getTrIndex(uint32_t instTrIndex);
	//LiteMath::float3& getVertex(unsigned instance_ind, unsigned* instances_ptr, unsigned geom_vert_ind);
	//LiteMath::float3x3& getTrVertices(unsigned instTrIndex);
	//LiteMath::float3& getTrVerticesAxis(unsigned instTrIndex);


	float calcSurf(const simpleAABB& aabb) { return 2.0f * ((aabb.max.x - aabb.min.x) +
															(aabb.max.y - aabb.min.y) +
															(aabb.max.z - aabb.min.z));
	}

	float calculateSAH(unsigned& new_axis, unsigned& new_index, unsigned *instances_ptr,
					   simpleAABB& parentAABB, simpleAABB& leftAABB, simpleAABB& rightAABB,
					   unsigned elem_count, float parentSAH);
	void calculateUpperTreeSAH(unsigned& new_axis, unsigned& new_index, unsigned* instances_ptr,
						simpleAABB& parentAABB, simpleAABB& leftAABB, simpleAABB& rightAABB,
						unsigned elem_count);


  /* ========  Tree builder  ======== */
	unsigned buildHeader(bool is_leaf, bool is_carve, unsigned info);
	bool sortTempPlanes(std::vector <unsigned>& target_plane_vec, const LiteMath::float3& minDif, const LiteMath::float3& maxDif,
						float curr_difference, unsigned plane_info, bool dif_consition);
	unsigned buildCarvingNodes(simpleAABB& parent_aabb, simpleAABB& child_aabb, unsigned axis,
								bool is_left_child, bool is_leaf, unsigned triangle_index);
	void dstBuilderRecur(unsigned first, unsigned last, float parentSAH);
	void dstBuilderRecurUpper(unsigned first, unsigned last);


  /* ========  Ray traverse  ======== */
	std::vector <unsigned> DSTree::TreePath(unsigned init, unsigned find);	//finding the fastest way to a certain node
	bool traceAABB(const simpleAABB& tempAABB, LiteMath::float3 Position, LiteMath::float3 Direction, LiteMath::float3 InvDir, LiteMath::float2& tMinMax);
	CRT_Hit traceTriangle(LiteMath::float3 Position, LiteMath::float3 Direction, LiteMath::float4* tempVertices, unsigned* tempInsices);
	CRT_Hit findHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, bool findAny, uint32_t current_instance);
	void DSTree::findInstHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, std::vector <uint32_t> &insts);

public:
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