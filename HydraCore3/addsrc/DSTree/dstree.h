#pragma once

#include "../../../RayTracing/CrossRT.h"
#include "../Vis/vis.h"
//#include "../AddHeaders/CrossRT.h"
#include <iostream>
//#include <fstream>
//#include <iomanip>


//16 bytes
struct DSNode {
	unsigned leftChild = 0u;
	unsigned rightNode = 0u;
	float    planes[2] = { 0.f, 0.f };
};


struct simpleDSTinfo {
	float p1, p2;
	int p1axis, p2axis;
	bool is_1max, is_2max, is_leaf, is_carve;
	unsigned depth;
};

struct simpleAABB {
	LiteMath::float3 min;
	LiteMath::float3 max;

	simpleAABB() : min(), max() {}
	simpleAABB(const LiteMath::float3& new_min, const LiteMath::float3& new_max) : min(new_min), max(new_max) {}
	simpleAABB(const simpleAABB &temp) : min(temp.min), max(temp.max) {}
};

struct simpleInstance {
	uint32_t geomID;
	LiteMath::float4x4 transf;
	LiteMath::float4x4 transfInv;
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
	std::vector <LiteMath::float3> vertices;
	std::vector <unsigned> indices;
	std::vector <unsigned> indices_sorted;
	std::vector <simpleMeshInfo> meshes;
	std::vector <simpleInstance> instances_info;

	std::vector <DSNode> upper_tree;  //leaves contain instances' indices
	std::vector <DSNode> lower_tree;  //stores arrays of meshes' DST
	std::vector <unsigned> instances; //temp array for builder; stores indices while sorting
	std::vector <DSNode> dst_nodes;   //temp array for builder; builder output
	std::vector <std::vector <float>> tree_info; // DS trees info [mesh][tree layer]
	std::vector <std::vector <unsigned>> layers_info; // DS trees info [mesh][tree layer]


/*  ============  File i/o  ============  */
	void dst_input (char* path);
	void dst_output(char* path);


/* ================  Comparator functions  ================ */
	float getMaxElem(unsigned *indices_ptr, LiteMath::float3 *vertices_ptr, unsigned current_scene_axis);
	void qsortUpper(unsigned *inst_arr, unsigned count, unsigned current_scene_axis);
	void qsortLower(unsigned *inst_arr, unsigned count, unsigned current_scene_axis, LiteMath::float3* tempVertices, unsigned* tempIndices);


  /* ========  SAH calculation  ======== */
	const float EMPTY_COST = 3e3f;
	const float SAH_MAX = 10e10f;
	const unsigned DSTREE_MAX_POLY = 8u;

	simpleAABB scene_AABB;
	simpleAABB getAABB(unsigned index, const LiteMath::float3* tempVertices, const unsigned* tempIndices);
	simpleAABB mergeAABB(unsigned index, simpleAABB aabb, const LiteMath::float3* tempVertices, const unsigned* tempIndices);
	simpleAABB getInstAABB(unsigned index);
	simpleAABB mergeInstAABB(unsigned index, simpleAABB aabb);


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
	unsigned buildCarvingNodes(simpleAABB& parent_aabb, simpleAABB& child_aabb, unsigned axis,
								bool is_left_child, bool is_leaf, unsigned triangle_index);
	void dstBuilderRecur(unsigned first, unsigned last, float parentSAH);
	void dstBuilderRecurUpper(unsigned first, unsigned last);


  /* ========  Ray traverse  ======== */
	std::vector <unsigned> TreePath(unsigned init, unsigned find);	//finding the fastest way to a specific node
	bool traceAABB(const simpleAABB& tempAABB, LiteMath::float3 Position, LiteMath::float3 Direction, LiteMath::float3 InvDir, LiteMath::float2& tMinMax);
	CRT_Hit traceTriangle(LiteMath::float3 Position, LiteMath::float3 Direction, LiteMath::float3* tempVertices, unsigned* tempInsices);
	CRT_Hit findHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, bool findAny, uint32_t current_instance, float triang_t);
	void treeVisInfo();
	uint32_t findInstHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, uint32_t* insts);

public:
	unsigned AddGeom_Triangles3f(const float* a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel = BUILD_HIGH, size_t vByteStride = sizeof(float) * 3);
	void UpdateGeom_Triangles3f(uint32_t a_geomId, const float* a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel = BUILD_HIGH, size_t vByteStride = sizeof(float) * 3);

	CRT_Hit RayQuery_NearestHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar);
	bool RayQuery_AnyHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar);

	uint32_t AddInstance(uint32_t a_geomId, const LiteMath::float4x4& a_matrix);
	void UpdateInstance(uint32_t a_instanceId, const LiteMath::float4x4& a_matrix);
	void ClearGeom();
	void ClearScene();
	void CommitScene(BuildQuality a_qualityLevel = BUILD_MEDIUM);
};