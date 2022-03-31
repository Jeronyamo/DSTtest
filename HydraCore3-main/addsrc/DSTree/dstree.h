#pragma once

#include "dst.h"
#include "../../../RayTracing/CrossRT.h"

struct simpleAABB {
	LiteMath::float3 min;
	LiteMath::float3 max;

	simpleAABB() : min(), max() {}
	simpleAABB(const LiteMath::float3& new_min, const LiteMath::float3& new_max) : min(new_min), max(new_max) {}
}; 


class DSTree : public ISceneObject {
	std::vector <LiteMath::float4> vertices;
	std::vector <unsigned> indices;
	std::vector <unsigned> instances;
	std::vector <float> transform_mat;
	LiteMath::float4x4 invMVP;
	LiteMath::float3 camPos;

	DST dst_nodes;


/* ================  Functions  ================ */
	friend float getMinPos(const void* index);
	friend float getMaxPos(const void* index);
  /* ========  SAH calculation  ======== */
	const unsigned DSTREE_MAX_POLY = 16u;
	const float EMPTY_COST = 3e3f;
	const float SAH_MAX = 10e10f;
	unsigned current_scene_axis = 3u;

	float sceneAABB_arr[6u];
	simpleAABB scene_AABB;
	simpleAABB getAABB(unsigned index, unsigned * instances_ptr);
	simpleAABB mergeAABB(unsigned index, simpleAABB aabb, unsigned* instances_ptr);

	LiteMath::float4x4 getTransformMatrix(unsigned index) { return LiteMath::float4x4(&(transform_mat[index << 4])); }
	unsigned findSameMatrix(const LiteMath::float4x4& a_matrix);
	LiteMath::float3 getVertex(unsigned instance_ind, unsigned* instances_ptr, unsigned geom_vert_ind);


	float calcSurf(const simpleAABB& aabb) { return 2.0f * ((aabb.max.x - aabb.min.x) +
															(aabb.max.y - aabb.min.y) +
															(aabb.max.z - aabb.min.z));
	}

	//float getMinPos(const void* index);
	//float getMaxPos(const void* index);
	//int cmpMinPos(const void* ind1, const void* ind2); //Check if this works with qsort()
	//int cmpMaxPos(const void* ind1, const void* ind2); //Check if this works with qsort()
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
	enum BoolParams { CHECK_LEAF_POLY, IS_CARVING_NODE, IS_DOUBLE_CARVE };

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

	bool traceAABB(LiteMath::float3 Position, LiteMath::float3 Direction, LiteMath::float3 InvDir);
	CRT_Hit traceTriangle(LiteMath::float3 Position, LiteMath::float3 Direction, unsigned trIndex);
	CRT_Hit findHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, bool findAny);

public:
	void setVertices(const std::vector <LiteMath::float4>& new_vertices) { vertices = new_vertices; }
	void setIndices(const std::vector <unsigned>& new_indices) { indices = new_indices; }
	std::vector <unsigned>& getIndices() { return indices; }
	void setMVPinv(const LiteMath::float4x4& new_invMVP) { invMVP = new_invMVP; }
	float* getSceneAABB();
	std::vector <DSNode>& getDSTNodes() { return dst_nodes.nodes; }
	//void setCamPos(const LiteMath::float3& new_camPos) { camPos = new_camPos; }

	unsigned AddGeom_Triangles4f(const LiteMath::float4* a_vpos4f, size_t a_vertNumber, const unsigned* a_triIndices, size_t a_indNumber);
	void UpdateGeom_Triangles4f(unsigned a_geomId, const LiteMath::float4* a_vpos4f, size_t a_vertNumber, const unsigned* a_triIndices, size_t a_indNumber);

	CRT_Hit RayQuery_NearestHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar);
	bool RayQuery_AnyHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar);

	unsigned AddInstance(unsigned a_geomId, const LiteMath::float4x4& a_matrix);
	void UpdateInstance(unsigned a_instanceId, const LiteMath::float4x4& a_matrix);
	void ClearGeom();
	void ClearScene();
	void CommitScene();
};