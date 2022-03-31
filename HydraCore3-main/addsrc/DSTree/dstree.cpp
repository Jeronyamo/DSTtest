#include "dstree.h"
#include <stdlib.h>


unsigned DSTree::AddGeom_Triangles4f(const LiteMath::float4* a_vpos4f, size_t a_vertNumber, const unsigned* a_triIndices, size_t a_indNumber) {
	for (size_t j = 0; j < a_vertNumber; ++j)
		vertices.push_back(a_vpos4f[j]);

	for (size_t j = 0; j < a_indNumber; ++j)
		indices.push_back(a_triIndices[j]);

	unsigned index = indices.size() / 3;
	if (index) --index;
	return index;
}

void DSTree::UpdateGeom_Triangles4f(unsigned a_geomId, const LiteMath::float4* a_vpos4f, size_t a_vertNumber, const unsigned* a_triIndices, size_t a_indNumber) {
	for (size_t j = 0; j < (a_vertNumber > a_indNumber ? a_vertNumber : a_indNumber); ++j) {
		if (j < a_vertNumber)
			vertices[indices[a_vertNumber * a_geomId + j]] = a_vpos4f[j];
		if (j < a_indNumber)
			indices[a_indNumber * a_geomId + j] = a_triIndices[j];
	}
}

unsigned DSTree::AddInstance(unsigned a_geomId, const LiteMath::float4x4& a_matrix) {
	unsigned mat_index = findSameMatrix(a_matrix);
	instances.push_back(a_geomId);
	instances.push_back(mat_index);

	if (mat_index == (transform_mat.size() >> 4)) {
		for (unsigned j = 0u; j < 16u; ++j) {
			transform_mat.push_back(a_matrix(j >> 2u, j & 3u));
		}
	}
	return (instances.size() >> 1) - 1u;
}

void DSTree::UpdateInstance(unsigned a_instanceId, const LiteMath::float4x4& a_matrix) {
	unsigned new_mat_index = findSameMatrix(a_matrix);
	bool newmat_is_unique = (new_mat_index == (transform_mat.size() >> 4));

	unsigned old_mat_index = instances[unsigned(a_instanceId << 1) + 1u];
	bool oldmat_is_used_once = true;

	for (unsigned i = 1u; i < instances.size(); i += 2u) {
		if (instances[i] == old_mat_index && (i >> 1) != a_instanceId) {
			oldmat_is_used_once = false;
			break;
		}
	}

	if (oldmat_is_used_once && newmat_is_unique) {
		for (unsigned i = 0u; i < 16u; ++i) {
			transform_mat[(old_mat_index << 4) + i] = a_matrix(i >> 2u, i & 3u);
		}
		return;
	}

	if (newmat_is_unique) {
		for (unsigned i = 0u; i < 16u; ++i) {
			transform_mat.push_back(a_matrix(i >> 2u, i & 3u));
		}
	}
	instances[(a_instanceId << 1) + 1u] = new_mat_index;

	if (oldmat_is_used_once) {
		transform_mat.erase(transform_mat.begin() + (old_mat_index << 4), transform_mat.begin() + ((old_mat_index + 1u) << 4));

		for (unsigned i = 1u; i < instances.size(); i += 2u) {
			if (instances[i] > old_mat_index)
				--instances[i];
		}
	}
}

void DSTree::ClearGeom() {
	vertices.clear();
	indices.clear();
	instances.clear();
	transform_mat.clear();
}

void DSTree::ClearScene() {
	instances.clear();
	transform_mat.clear();
}


/* ========  SAH calculation  ======== */

DSTree *scene;

float max(const LiteMath::float2& params) { return params[params.x < params.y]; }
float min(const LiteMath::float2& params) { return params[params.x > params.y]; }


simpleAABB DSTree::getAABB(unsigned index, unsigned* instances_ptr) {
	simpleAABB tempAABB;
	LiteMath::float3 a = getVertex(index, instances_ptr, 0u),
					 b = getVertex(index, instances_ptr, 1u),
					 c = getVertex(index, instances_ptr, 2u);

	for (unsigned i = 0u; i < 3u; ++i) {
		tempAABB.min[i] = min(LiteMath::float2(a[i], min(LiteMath::float2(b[i], c[i]))));
		tempAABB.max[i] = max(LiteMath::float2(a[i], max(LiteMath::float2(b[i], c[i]))));
	}
	/*
	for (unsigned i = 0u; i < 3u; ++i) {
		LiteMath::float3 coords = LiteMath::float3(vertices[indices_ptr[3*index    ]][i],
												   vertices[indices_ptr[3*index + 1]][i],
												   vertices[indices_ptr[3*index + 2]][i]);
		tempAABB.min[i] = min(LiteMath::float2(coords.x, min(LiteMath::float2(coords.y, coords.z))));
		tempAABB.max[i] = max(LiteMath::float2(coords.x, max(LiteMath::float2(coords.y, coords.z))));
	}
	*/
	return tempAABB;
}

simpleAABB DSTree::mergeAABB(unsigned index, simpleAABB aabb, unsigned * instances_ptr) {
	LiteMath::float3 a = getVertex(index, instances_ptr, 0u),
					 b = getVertex(index, instances_ptr, 1u),
					 c = getVertex(index, instances_ptr, 2u);
	
	for (unsigned i = 0u; i < 3u; ++i) {
		aabb.min[i] = min(LiteMath::float2(min(LiteMath::float2(aabb.min[i], a[i])), min(LiteMath::float2(b[i], c[i]))));
		aabb.max[i] = max(LiteMath::float2(max(LiteMath::float2(aabb.max[i], a[i])), max(LiteMath::float2(b[i], c[i]))));
	}
	/*
	for (unsigned i = 0u; i < 3u; ++i) {
		LiteMath::float3 coords = LiteMath::float3(vertices[indices_ptr[3*index    ]][i],
												   vertices[indices_ptr[3*index + 1]][i],
												   vertices[indices_ptr[3*index + 2]][i]);

		aabb.min[i] = min(LiteMath::float2(min(LiteMath::float2(aabb.min[i], coords.x)), min(LiteMath::float2(coords.y, coords.z))));
		aabb.max[i] = max(LiteMath::float2(max(LiteMath::float2(aabb.max[i], coords.x)), max(LiteMath::float2(coords.y, coords.z))));
	}
	*/
	return aabb;
}

float *DSTree::getSceneAABB() {
	sceneAABB_arr[0u] = scene_AABB.min.x;
	sceneAABB_arr[1u] = scene_AABB.min.y;
	sceneAABB_arr[2u] = scene_AABB.min.z;
	sceneAABB_arr[3u] = scene_AABB.max.x;
	sceneAABB_arr[4u] = scene_AABB.max.y;
	sceneAABB_arr[5u] = scene_AABB.max.z;
	return sceneAABB_arr;
}

unsigned DSTree::findSameMatrix(const LiteMath::float4x4& a_matrix) {
	unsigned mat_index = (transform_mat.size() >> 4);

	for (unsigned i = 0u; i < mat_index; ++i) {
		LiteMath::float4x4 temp_mat = getTransformMatrix(i);
		bool mat_equal = true;

		for (unsigned j = 0u; (j < 16u) && mat_equal; ++j) {
			mat_equal = (temp_mat(j >> 2u, j & 3u) == a_matrix(j >> 2u, j & 3u));
		}

		if (mat_equal) {
			mat_index = i;
			break;
		}
	}

	return mat_index;
}

LiteMath::float3 DSTree::getVertex(unsigned instance_ind, unsigned* instances_ptr, unsigned geom_vert_ind) {
	LiteMath::float4 coords = LiteMath::float4x4(&(transform_mat[instances_ptr[(instance_ind << 1) + 1u] << 4])) *
													vertices[indices[3 * instances_ptr[instance_ind << 1] + geom_vert_ind]];

	return LiteMath::float3(coords.x, coords.y, coords.z);
}


float getMinPos(const void* index) {
	unsigned* temp_index = (unsigned*) index;
	unsigned mat_index = temp_index[1u];

	LiteMath::float3 coords = LiteMath::float3(scene->getVertex(0u, temp_index, 0u)[scene->current_scene_axis],
											   scene->getVertex(0u, temp_index, 1u)[scene->current_scene_axis],
											   scene->getVertex(0u, temp_index, 2u)[scene->current_scene_axis]);

	return min(LiteMath::float2(coords.x, min(LiteMath::float2(coords.y, coords.z))));
}

float getMaxPos(const void* index) {
	unsigned* temp_index = (unsigned*) index;
	unsigned mat_index = temp_index[1u];

	LiteMath::float3 coords = LiteMath::float3(scene->getVertex(0u, temp_index, 0u)[scene->current_scene_axis],
											   scene->getVertex(0u, temp_index, 1u)[scene->current_scene_axis],
											   scene->getVertex(0u, temp_index, 2u)[scene->current_scene_axis]);

	return max(LiteMath::float2(coords.x, max(LiteMath::float2(coords.y, coords.z))));
}

int cmpMinPos(const void* ind1, const void* ind2) {
	float min1 = getMinPos(ind1);
	float min2 = getMinPos(ind2);

	return int(min1 > min2) - int(min1 < min2);
}

int cmpMaxPos(const void* ind1, const void* ind2) {
	float max1 = getMaxPos(ind1);
	float max2 = getMaxPos(ind2);

	return int(max1 > max2) - int(max1 < max2);
}

float DSTree::calculateSAH(unsigned& new_axis, unsigned& new_index, unsigned *instances_ptr,
						   simpleAABB& parentAABB, simpleAABB& leftAABB, simpleAABB& rightAABB,
						   unsigned elem_count, float parentSAH) {
	unsigned index = 0u;
	new_axis = 3u;

	for (unsigned curr_axis = 0u; curr_axis < 3u; ++curr_axis) {
		std::vector <simpleAABB> right_surfs;

		current_scene_axis = curr_axis;
		qsort(instances_ptr, elem_count, 2u * sizeof(unsigned), cmpMaxPos);

		right_surfs.push_back(getAABB(elem_count - 1u, instances_ptr));
		for (unsigned j = elem_count - 2u; j > 0u; --j) {
			right_surfs.push_back(mergeAABB(j, right_surfs[right_surfs.size() - 1u], instances_ptr));
		}
		parentAABB = mergeAABB(0, right_surfs[right_surfs.size() - 1u], instances_ptr);

		simpleAABB temp_aabb(getAABB(0u, instances_ptr));
		for (unsigned k = 1u; k < elem_count - 1u; ++k) {
			temp_aabb = mergeAABB(k - 1u, temp_aabb, instances_ptr);

			float tempSAH = EMPTY_COST + k * calcSurf(temp_aabb) +
				(elem_count - k) * calcSurf(right_surfs[right_surfs.size() - k]);

			if (tempSAH < parentSAH) {
				parentSAH = tempSAH;
				rightAABB = right_surfs[right_surfs.size() - k];
				leftAABB = temp_aabb;
				new_index = k - 1u;
				new_axis = curr_axis;
			}
		}

		right_surfs.clear();
	}

	if (new_axis < 2u) {
		current_scene_axis = new_axis;
		qsort(instances_ptr, elem_count, 2u * sizeof(unsigned), cmpMaxPos);
	}

	return parentSAH;
}


/* ========  Tree builder  ======== */

unsigned DSTree::buildHeader(bool is_leaf, bool is_carve, unsigned info) {
	unsigned header = unsigned(is_leaf) + (unsigned(is_carve) << 1);
	unsigned field1 = (info >> 2) & 3u, field2 = info & 3u;

	return (field1 << 4u) + (field2 << 2u) + header;
}

bool DSTree::sortTempPlanes(std::vector <unsigned>& target_plane_vec, const LiteMath::float3& minDif, const LiteMath::float3& maxDif,
							float curr_difference, unsigned plane_info, bool dif_consition) {
	if (dif_consition) {
		for (std::vector <unsigned>::iterator planesInd = target_plane_vec.begin(); planesInd <= target_plane_vec.end(); ++planesInd) {

			if ((planesInd == target_plane_vec.end() ||
				(!(plane_info & 1u) && curr_difference <= minDif[(*planesInd) >> 1])  ||
				( (plane_info & 1u) && curr_difference <= maxDif[(*planesInd) >> 1]))) {
				target_plane_vec.insert(planesInd, plane_info);
				return true;
			}
		}
	}
	return false;
}

unsigned DSTree::buildCarvingNodes(simpleAABB& parent_aabb, simpleAABB& child_aabb, unsigned axis,
									bool is_left_child, bool is_leaf, unsigned triangle_index) {
	unsigned size1 = dst_nodes.nodes.size();
	std::vector <unsigned> new_possible_planes;
	std::vector <unsigned> new_planes;
	LiteMath::float3	minDif(child_aabb.min.x - parent_aabb.min.x,
							   child_aabb.min.y - parent_aabb.min.y,
							   child_aabb.min.z - parent_aabb.min.z),
						maxDif(parent_aabb.max.x - child_aabb.max.x,
							   parent_aabb.max.y - child_aabb.max.y,
							   parent_aabb.max.z - child_aabb.max.z);

	for (unsigned i = 0; i < 3u; ++i) {
		if (!(!is_left_child && i == axis)) {
			float curr_difference = child_aabb.min[i] - parent_aabb.min[i];
			if (!sortTempPlanes(new_planes, minDif, maxDif, curr_difference, i << 1, curr_difference >= 0.1f))
				sortTempPlanes(new_possible_planes, minDif, maxDif, curr_difference, i << 1,
								curr_difference > LiteMath::EPSILON && curr_difference < 0.1f);
		}

		if (!(is_left_child && i == axis)) {
			float curr_difference = parent_aabb.max[i] - child_aabb.max[i];
			if (!sortTempPlanes(new_planes, minDif, maxDif, curr_difference, (i << 1) + 1u, curr_difference >= 0.1f))
				sortTempPlanes(new_possible_planes, minDif, maxDif, curr_difference, (i << 1) + 1u,
								curr_difference > LiteMath::EPSILON && curr_difference < 0.1f);
		}
	}

	DSNode new_carving_node;

	for (unsigned ind = 0u; ind < (new_planes.size() >> 1u) + (new_planes.size() & 1u); ind += 2u) {
		unsigned plane1 = new_planes[ind], plane2 = new_planes[ind] ^ 1u;

		if (ind < new_planes.size() - 1u)
			plane2 = new_planes[ind + 1u];

		if (ind == new_planes.size() - 1u && !new_possible_planes.empty()) {
			plane2 = new_possible_planes[0u];
			new_possible_planes.erase(new_possible_planes.begin());
		}

		bool double_carve = (plane1 >> 1) != (plane2 >> 1);

		if (((plane1 >> 1) > (plane2 >> 1)) || (!double_carve && ((plane1 & 1u) > (plane2 & 1u)))) {
			unsigned temp = plane1;
			plane1 = plane2;
			plane2 = temp;
		}

		new_carving_node.planes[0u] = child_aabb.min[plane1 >> 1];
		if (plane1 & 1u)
			new_carving_node.planes[0u] = child_aabb.max[plane1 >> 1];

		new_carving_node.planes[1u] = child_aabb.min[plane2 >> 1];
		if (plane2 & 1u)
			new_carving_node.planes[1u] = child_aabb.max[plane2 >> 1];

		new_carving_node.rightNode = 0u;

		unsigned field1 = 2u, field2 = plane1 >> 1;
		bool last_carve = (ind + 2u) >= ((new_planes.size() >> 1u) + (new_planes.size() & 1u));
		unsigned left_child = dst_nodes.nodes.size() + 1u;
		if (is_leaf && last_carve) {
			left_child = triangle_index;
		}

		if (double_carve) {
			field1 = (unsigned((plane1 >> 1) == 1u) << 1) + unsigned((plane2 >> 1) == 2u);
			field2 = ((plane1 & 1u) << 1) + (plane2 & 1u);
		}

		new_carving_node.leftChild = (left_child << 6) +
			buildHeader(is_leaf && last_carve, true, (field1 << 2) + field2);

		dst_nodes.nodes.push_back(new_carving_node);
	}
	return dst_nodes.nodes.size() - size1;
}

void DSTree::dstBuilderRecur(unsigned first, unsigned last, float parentSAH) {
	unsigned curr_node_index = dst_nodes.nodes.size();
	unsigned new_index = last;

	//std::cout << "Node " << curr_node_index << std::endl;
	dst_nodes.nodes.emplace_back();
	dst_nodes.nodes[curr_node_index].leftChild = (first << 6u) + 1u;
	dst_nodes.nodes[curr_node_index].planes[0] = float(last - first + 1u);

	if (last - first >= DSTREE_MAX_POLY) {
		simpleAABB left_aabb, right_aabb, curr_aabb;
		unsigned new_axis = 3u;
		float SAH = calculateSAH(new_axis, new_index, &(instances[2u*first]),
								 curr_aabb, left_aabb, right_aabb, unsigned(last - first + 1u), parentSAH);

		if (new_axis < 3u && SAH <= parentSAH) {
			parentSAH = SAH;
			new_index += first;
			dst_nodes.nodes[curr_node_index].planes[0] = left_aabb.max[new_axis];
			dst_nodes.nodes[curr_node_index].planes[1] = right_aabb.min[new_axis];
			dst_nodes.nodes[curr_node_index].leftChild = (dst_nodes.nodes.size() << 6u) + buildHeader(false, false, new_axis);

			buildCarvingNodes(curr_aabb, left_aabb, new_axis, true, (new_index - first) < DSTREE_MAX_POLY, first);
			if (!((new_index - first) < DSTREE_MAX_POLY && (dst_nodes.nodes.size() - 1u != curr_node_index)))
				dstBuilderRecur(first, new_index, parentSAH);

			for (unsigned temp_index = curr_node_index + 1u;;) {
				dst_nodes.nodes[temp_index].rightNode = dst_nodes.nodes.size();
				if (dst_nodes.nodes[temp_index].leftChild & 1u)
					break;

				if (dst_nodes.nodes[temp_index].leftChild & 2u) {
					++temp_index;
				}
				else if (!(dst_nodes.nodes[temp_index].leftChild & 3u))
					temp_index = dst_nodes.nodes[dst_nodes.nodes[temp_index].leftChild >> 6].rightNode;
			}

			unsigned carve_nodes_num = buildCarvingNodes(curr_aabb, right_aabb, new_axis, false, (last - new_index - 1u) < DSTREE_MAX_POLY, new_index + 1u);

			if (!((last - new_index - 1u) < DSTREE_MAX_POLY && carve_nodes_num))
				dstBuilderRecur(new_index + 1u, last, parentSAH);
		}
	}
	return;
}

void DSTree::CommitScene() {
	if (indices.empty())
		return;

	scene = this;
	dstBuilderRecur(0u, (instances.size() >> 1) - 1u, SAH_MAX);

	scene_AABB = getAABB(0u, instances.data());
	for (unsigned i = 1u; i < (instances.size() >> 1); ++i) {
		scene_AABB = mergeAABB(i, scene_AABB, instances.data());
	}
	scene = nullptr;
	return;
}


/* ========  Ray traverse  ======== */

bool DSTree::traceAABB(LiteMath::float3 Position, LiteMath::float3 Direction, LiteMath::float3 InvDir) {
	float tMin = LiteMath::INF_NEGATIVE, tMax = LiteMath::INF_POSITIVE;
	LiteMath::float3 Min = scene_AABB.min, Max = scene_AABB.max;

	for (unsigned dim = 0u; dim < 3u; ++dim) {
		if (Direction[dim] > -LiteMath::EPSILON && Direction[dim] < LiteMath::EPSILON) {
			if (Position[dim] < Min[dim] || Position[dim] > Max[dim])
				return false;

			continue;
		}

		float t1 = (Min[dim] - Position[dim]) * InvDir[dim];
		float t2 = (Max[dim] - Position[dim]) * InvDir[dim];

		if (t1 > t2) {
			float temp_t = t1;
			t1 = t2;
			t2 = temp_t;
		}

		if (t1 > tMin)
			tMin = t1;

		if (t2 < tMax)
			tMax = t2;

		if (tMin > tMax || tMax < 0.0f)
			return false;
	}

	return tMin < tMax && tMax >= 0.0f;
}

CRT_Hit DSTree::traceTriangle(LiteMath::float3 Position, LiteMath::float3 Direction, unsigned trIndex) {
	CRT_Hit tempHitInfo;
	tempHitInfo.t = LiteMath::INF_POSITIVE;
	tempHitInfo.primId = trIndex;

	/*  ====  May be remade using getVertex()  ====  */
	/*
	const LiteMath::uint3 tr_index = LiteMath::uint3(indices[3 * instances[trIndex << 1]],
													 indices[3 * instances[trIndex << 1] + 1],
													 indices[3 * instances[trIndex << 1] + 2]);
	const LiteMath::float4 temp_a = getTransformMatrix(instances[(trIndex << 1) + 1]) *
									LiteMath::float4(vertices[tr_index.x].x, vertices[tr_index.x].y, vertices[tr_index.x].z, 0.f);
	const LiteMath::float4 temp_b = getTransformMatrix(instances[(trIndex << 1) + 1]) *
									LiteMath::float4(vertices[tr_index.y].x, vertices[tr_index.y].y, vertices[tr_index.y].z, 0.f) - temp_a;
	const LiteMath::float4 temp_c = getTransformMatrix(instances[(trIndex << 1) + 1]) *
									LiteMath::float4(vertices[tr_index.z].x, vertices[tr_index.z].y, vertices[tr_index.z].z, 0.f) - temp_a;

	const LiteMath::float3  a = LiteMath::float3(temp_a.x, temp_a.y, temp_a.z);
	const LiteMath::float3 E1 = LiteMath::float3(temp_b.x, temp_b.y, temp_b.z) - a;
	const LiteMath::float3 E2 = LiteMath::float3(temp_c.x, temp_c.y, temp_c.z) - a;
	*/
	const LiteMath::float3  a = getVertex(trIndex, instances.data(), 0u);
	const LiteMath::float3 E1 = getVertex(trIndex, instances.data(), 1u) - a;
	const LiteMath::float3 E2 = getVertex(trIndex, instances.data(), 2u) - a;
	/*  ===============  */
	const LiteMath::float3  P = cross(Direction, E2);
	const float d = dot(P, E1);

	if (d < -LiteMath::EPSILON || d > LiteMath::EPSILON) {
		LiteMath::float3 T = Position - a;
		const float f = 1.0f / d;
		const float temp_u = f * dot(P, T);
		T = cross(T, E1);
		const float temp_v = f * dot(Direction, T);
		const float temp_t = f * dot(E2, T);

		if (temp_u >= 0.0f && temp_u <= 1.0f && ((temp_u + temp_v) <= 1.0f) &&
			temp_v >= 0.0f && temp_t >= 0.0f && temp_t < tempHitInfo.t) {
			tempHitInfo.coords[0u] = temp_u;
			tempHitInfo.coords[1u] = temp_v;
			tempHitInfo.t = temp_t;
		}
	}
	return tempHitInfo;
}

CRT_Hit DSTree::RayQuery_NearestHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar) {
	return this->findHit(posAndNear, dirAndFar, false);
}

bool DSTree::RayQuery_AnyHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar) {
	CRT_Hit temp_hit = this->findHit(posAndNear, dirAndFar, true);

	return (temp_hit.t >= posAndNear.w) && (temp_hit.t <= dirAndFar.w);
}

CRT_Hit DSTree::findHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, bool findAny) {
	CRT_Hit temp_hit;
	temp_hit.t = LiteMath::INF_POSITIVE;
	LiteMath::float3 position  = LiteMath::float3(posAndNear.x, posAndNear.y, posAndNear.z);
	LiteMath::float3 direction = LiteMath::float3( dirAndFar.x,  dirAndFar.y,  dirAndFar.z);
	LiteMath::float3 invDir    = LiteMath::float3(1.f / dirAndFar.x, 1.f / dirAndFar.y, 1.f / dirAndFar.z);
	unsigned node = 0u;

	if (!traceAABB(position, direction, invDir)) return temp_hit;

	do {
		DSNode curr_node;
		unsigned header = 0u;
		bool3 parameter = bool3(false);

		do {
			curr_node = dst_nodes.nodes[node];
			node = curr_node.rightNode;
			header = curr_node.leftChild & 63u;
			curr_node.leftChild = curr_node.leftChild >> 6;
			parameter = bool3(bool(header & 1u),										//parameter[CHECK_LEAF_POLY]
							  bool(header & 2u),										//parameter[IS_CARVING_NODE]
							  bool(!((header & 48u) == 32u) && bool(header & 2u)));		//parameter[IS_DOUBLE_CARVE]

			if (!parameter[CHECK_LEAF_POLY] || parameter[IS_CARVING_NODE]) {
				LiteMath::uint2 planes = LiteMath::uint2((header & 12u) >> 2);
				bool2 is_normal_positive = bool2(!parameter[IS_CARVING_NODE], parameter[IS_CARVING_NODE]);

				if (parameter[IS_DOUBLE_CARVE]) {
					is_normal_positive = bool2(bool(planes.x & 2u), bool(planes.x & 1u));
					planes = LiteMath::uint2((header >> 5) & 1u, ((header >> 4) & 1u) + 1u);
				}

				bool2 planeTrav = bool2(is_normal_positive[0u] == (position[planes[0u]] < curr_node.planes[0u]),
										is_normal_positive[1u] == (position[planes[1u]] < curr_node.planes[1u]));

				if (!planeTrav.x && (direction[planes[0u]] < -LiteMath::EPSILON || direction[planes[0u]] > LiteMath::EPSILON)) {
					float temp_t = invDir[planes[0u]] * (curr_node.planes[0u] - position[planes[0u]]);
					planeTrav.x = temp_t >= 0.f && temp_t < temp_hit.t;
				}
				if (!planeTrav.y && (direction[planes[1u]] < -LiteMath::EPSILON || direction[planes[1u]] > LiteMath::EPSILON)) {
					float temp_t = invDir[planes[1u]] * (curr_node.planes[1u] - position[planes[1u]]);
					planeTrav.y = temp_t >= 0.f && temp_t < temp_hit.t;
				}
				if (planeTrav[0u] && (!parameter[IS_CARVING_NODE] || (planeTrav[1u] && !parameter[CHECK_LEAF_POLY])))
					node = curr_node.leftChild;
				if (!planeTrav[0u] && planeTrav[1u] && !parameter[IS_CARVING_NODE])
					node = dst_nodes.nodes[curr_node.leftChild].rightNode;

				parameter[CHECK_LEAF_POLY] = parameter[CHECK_LEAF_POLY] &&
					(!parameter[IS_CARVING_NODE] || (planeTrav[0u] && planeTrav[1u]));
			}
		} while (!parameter[CHECK_LEAF_POLY] && node != 0u);

		if (parameter[CHECK_LEAF_POLY]) {
			unsigned tr_num = DSTREE_MAX_POLY;
			if (!parameter[IS_CARVING_NODE])
				tr_num = unsigned(curr_node.planes[0]);

			tr_num += curr_node.leftChild;
			if (tr_num > indices.size()) tr_num = indices.size();
			for (unsigned j = curr_node.leftChild; j < tr_num; ++j) {
				CRT_Hit temp_t = traceTriangle(position, direction, j);

				if (temp_t.t < temp_hit.t && temp_t.t <= dirAndFar.w && temp_t.t >= posAndNear.w) {
					temp_hit = temp_t;
					if (findAny) return temp_hit;
				}
			}
		}
	} while (node != 0u);
	return temp_hit;
}


ISceneObject* CreateSceneRT(const char* a_impleName) { return new DSTree; }

void DeleteSceneRT(ISceneObject* a_pScene) { delete a_pScene; }