#include "dstree.h"
#include <cfloat>

float max(const LiteMath::float2& params) { return params[params.x < params.y]; }
float min(const LiteMath::float2& params) { return params[params.x > params.y]; }

float max(float par1, float par2) {
	float params[2] = { par1, par2 };
	return params[params[0] < params[1]];
}
float min(float par1, float par2) {
	float params[2] = { par1, par2 };
	return params[params[0] > params[1]];
}


unsigned DSTree::AddGeom_Triangles4f(const LiteMath::float4* a_vpos4f, size_t a_vertNumber, const unsigned* a_triIndices, size_t a_indNumber) {
	simpleMeshInfo tempInfo;
	tempInfo.firstVertID  = -1;
	tempInfo.firstIndID   = -1;
	tempInfo.lastVertID   = -1;
	tempInfo.lastIndID    = -1;
	tempInfo.DSTreeOffset = -1;

	instances.resize(a_indNumber / 3);
	for (size_t j = 0; j < a_indNumber / 3; ++j) {
		instances[j] = j;
	}

	if (a_vertNumber) {
		tempInfo.firstVertID = vertices.size();

		for (size_t j = 0; j < a_vertNumber; ++j) {
			vertices.push_back(a_vpos4f[j]);
		}
		tempInfo.lastVertID = vertices.size() - 1;
	}
	if (a_indNumber) {
		tempInfo.firstIndID = indices.size() / 3;

		for (size_t j = 0; j < instances.size(); ++j) {
			indices.push_back(a_triIndices[3 * j    ]);
			indices.push_back(a_triIndices[3 * j + 1]);
			indices.push_back(a_triIndices[3 * j + 2]);
		}

		tempInfo.lastIndID = indices.size() / 3 - 1;
	}
	meshes.push_back(tempInfo);
	if (a_vertNumber && a_indNumber) {
		dst_nodes.clear();

		dstBuilderRecur(0u, instances.size() - 1, SAH_MAX);

		meshes[meshes.size() - 1].DSTreeOffset = lower_tree.size();
		lower_tree.insert(lower_tree.end(), dst_nodes.begin(), dst_nodes.end());
		dst_nodes.clear();

		meshes[meshes.size() - 1].meshAABB = getAABB(0, a_vpos4f, a_triIndices);
		for (uint32_t j = 1; j < instances.size(); ++j) {
			meshes[meshes.size() - 1].meshAABB = mergeAABB(j, meshes[meshes.size() - 1].meshAABB, a_vpos4f, a_triIndices);
		}

		for (int i = 0; i < 3; ++i) {
			if (meshes[meshes.size() - 1].meshAABB.min[i] == meshes[meshes.size() - 1].meshAABB.max[i]) {
				meshes[meshes.size() - 1].meshAABB.min[i] -= AABBeps;
				meshes[meshes.size() - 1].meshAABB.max[i] += AABBeps;
			}
		}
	}
	if (a_indNumber) {
		for (size_t j = 0; j < instances.size(); ++j) {
			indices_sorted.push_back(instances[j]);
		}
	}
	return meshes.size() - 1;
}

//not supported
void DSTree::UpdateGeom_Triangles4f(unsigned a_geomId, const LiteMath::float4* a_vpos4f, size_t a_vertNumber, const unsigned* a_triIndices, size_t a_indNumber) {
	/*
	uint32_t firstVert = meshes[a_geomId].firstVertID;
	uint32_t firstInd = meshes[a_geomId].firstIndID;
	for (uint32_t i = 0; i < a_vertNumber; ++i) {
		vertices[i + firstVert] = a_vpos4f[i];
	}
	for (uint32_t i = 0; i < a_indNumber; ++i) {
		indices[3 * firstInd + i] = a_triIndices[i];
	}
	*/
}

uint32_t DSTree::AddInstance(uint32_t a_geomId, const LiteMath::float4x4& a_matrix) {
	simpleInstance tempInst;
	tempInst.transf = a_matrix;
	tempInst.geomID = a_geomId;
	tempInst.instAABB.min = a_matrix * meshes[a_geomId].meshAABB.min;
	tempInst.instAABB.max = a_matrix * meshes[a_geomId].meshAABB.max;

	float edges[6] = { meshes[a_geomId].meshAABB.min.x, meshes[a_geomId].meshAABB.min.y, meshes[a_geomId].meshAABB.min.z,
					   meshes[a_geomId].meshAABB.max.x, meshes[a_geomId].meshAABB.max.y, meshes[a_geomId].meshAABB.max.z };

	for (unsigned i = 0u; i < 8u; ++i) {
		LiteMath::float3 tempPos{ edges[3u * (i >> 2)], edges[3u * ((i >> 1) & 1u) + 1u], edges[3u * (i & 1u) + 2u] };
		tempPos = a_matrix * tempPos;

		tempInst.instAABB.min.x = min(tempInst.instAABB.min.x, tempPos.x);
		tempInst.instAABB.min.y = min(tempInst.instAABB.min.y, tempPos.y);
		tempInst.instAABB.min.z = min(tempInst.instAABB.min.z, tempPos.z);

		tempInst.instAABB.max.x = max(tempInst.instAABB.max.x, tempPos.x);
		tempInst.instAABB.max.y = max(tempInst.instAABB.max.y, tempPos.y);
		tempInst.instAABB.max.z = max(tempInst.instAABB.max.z, tempPos.z);
	}


	instances_info.push_back(tempInst);
	return instances_info.size() - 1;
}

void DSTree::UpdateInstance(uint32_t a_instanceId, const LiteMath::float4x4& a_matrix) {
	simpleInstance tempInst = instances_info[a_instanceId];
	simpleMeshInfo tempMesh = meshes[tempInst.geomID];

	tempInst.transf = a_matrix;
	tempInst.instAABB.min = a_matrix * tempMesh.meshAABB.min;
	tempInst.instAABB.max = a_matrix * tempMesh.meshAABB.max;

	float edges[6] = { tempMesh.meshAABB.min.x, tempMesh.meshAABB.min.y, tempMesh.meshAABB.min.z,
					   tempMesh.meshAABB.max.x, tempMesh.meshAABB.max.y, tempMesh.meshAABB.max.z };

	for (unsigned i = 0u; i < 8u; ++i) {
		LiteMath::float3 tempPos{ edges[3u * (i >> 2)], edges[3u * ((i >> 1) & 1u) + 1u], edges[3u * (i & 1u) + 2u] };
		tempPos = a_matrix * tempPos;

		tempInst.instAABB.min.x = min(tempInst.instAABB.min.x, tempPos.x);
		tempInst.instAABB.min.y = min(tempInst.instAABB.min.y, tempPos.y);
		tempInst.instAABB.min.z = min(tempInst.instAABB.min.z, tempPos.z);

		tempInst.instAABB.max.x = max(tempInst.instAABB.max.x, tempPos.x);
		tempInst.instAABB.max.y = max(tempInst.instAABB.max.y, tempPos.y);
		tempInst.instAABB.max.z = max(tempInst.instAABB.max.z, tempPos.z);
	}

	instances_info[a_instanceId] = tempInst;
}

void DSTree::ClearGeom() {
	vertices.clear();
	indices.clear();
	instances.clear();
	instances_info.clear();
}

void DSTree::ClearScene() {
	instances.clear();
	instances_info.clear();
}


/* ========  SAH calculation  ======== */

simpleAABB DSTree::getAABB(unsigned index, const LiteMath::float4* tempVertices, const unsigned* tempIndices) {
	simpleAABB tempAABB;
	LiteMath::float4 a = tempVertices[tempIndices[3 * index    ]];
	LiteMath::float4 b = tempVertices[tempIndices[3 * index + 1]];
	LiteMath::float4 c = tempVertices[tempIndices[3 * index + 2]];

	for (unsigned i = 0u; i < 3u; ++i) {
		tempAABB.min[i] = min(a[i], min(b[i], c[i]));
		tempAABB.max[i] = max(a[i], max(b[i], c[i]));
	}
	return tempAABB;
}

simpleAABB DSTree::mergeAABB(unsigned index, simpleAABB aabb, const LiteMath::float4* tempVertices, const unsigned* tempIndices) {
	LiteMath::float4 a = tempVertices[tempIndices[3 * index    ]];
	LiteMath::float4 b = tempVertices[tempIndices[3 * index + 1]];
	LiteMath::float4 c = tempVertices[tempIndices[3 * index + 2]];

	for (unsigned i = 0u; i < 3u; ++i) {
		aabb.min[i] = min(min(aabb.min[i], a[i]), min(b[i], c[i]));
		aabb.max[i] = max(max(aabb.max[i], a[i]), max(b[i], c[i]));
	}
	return aabb;
}


simpleAABB DSTree::getInstAABB(unsigned index) {
	return instances_info[index].instAABB;
}

simpleAABB DSTree::mergeInstAABB(unsigned index, simpleAABB aabb) {
	simpleAABB tempAABB = instances_info[index].instAABB;

	for (unsigned i = 0u; i < 3u; ++i) {
		aabb.min[i] = min(tempAABB.min[i], aabb.min[i]);
		aabb.max[i] = max(tempAABB.max[i], aabb.max[i]);
	}
	return aabb;
}

void DSTree::qsortUpper(unsigned* inst_arr, unsigned count, unsigned current_scene_axis) {
	if (!count) return;

	unsigned* stack = new unsigned[count];
	unsigned first = 0, last = count - 1;
	int top = -1;
	unsigned stackSize = 0;
	stack[++top] = 0;
	stack[++top] = last;

	while (top > 0) {
		last = stack[top--];
		first = stack[top--];

		last -= first;
		unsigned* temp_inst = &(inst_arr[first]);
		int i = -1;
		float max2 = instances_info[temp_inst[last]].instAABB.max[current_scene_axis];

		for (unsigned j = 0u; j < last; ++j) {
			float max1 = instances_info[temp_inst[j]].instAABB.max[current_scene_axis];

			int res = int(max1 > max2) - int(max1 < max2);
			if (res <= 0) {
				unsigned temp = temp_inst[++i];
				temp_inst[i] = temp_inst[j];
				temp_inst[j] = temp;
			}
		}

		unsigned temp = temp_inst[++i];
		temp_inst[i] = temp_inst[last];
		temp_inst[last] = temp;

		i += first;
		last += first;

		if (i - 1 > long(first)) {
			stack[++top] = first;
			stack[++top] = i - 1;
		}

		if (i + 1 < last) {
			stack[++top] = i + 1;
			stack[++top] = last;
		}
	}
	delete[] stack;
}

void DSTree::qsortLower(unsigned *inst_arr, unsigned count, unsigned current_scene_axis, LiteMath::float4* tempVertices, unsigned* tempIndices) {
	if (!count) return;

	unsigned* stack = new unsigned[count];
	unsigned first = 0, last = count - 1;
	int top = -1;
	unsigned stackSize = 0;
	stack[++top] = 0;
	stack[++top] = last;

	while (top > 0) {
		last = stack[top--];
		first = stack[top--];

		last -= first;
		unsigned* temp_inst = &(inst_arr[first]);
		int i = -1;

		unsigned tempInd = 3 * temp_inst[last];
		unsigned cmp_indices[3] = { tempIndices[tempInd], tempIndices[tempInd + 1u], tempIndices[tempInd + 2u] };
		float vert[3] = { tempVertices[cmp_indices[0]][current_scene_axis],
						  tempVertices[cmp_indices[1]][current_scene_axis],
						  tempVertices[cmp_indices[2]][current_scene_axis] };
		float max2 = max(vert[vert[0] < vert[1]], vert[2]);

		for (unsigned j = 0u; j < last; ++j) {
			tempInd = 3 * temp_inst[j];
			unsigned cmp_indices1[3] = { tempIndices[tempInd], tempIndices[tempInd + 1u], tempIndices[tempInd + 2u] };
			float vert1[3] = { tempVertices[cmp_indices1[0]][current_scene_axis],
							   tempVertices[cmp_indices1[1]][current_scene_axis],
							   tempVertices[cmp_indices1[2]][current_scene_axis] };
			vert1[1] = vert1[1u + unsigned(vert1[1] < vert1[2])];
			float max1 = vert1[vert1[0] < vert1[1]];

			int res = int(max1 > max2) - int(max1 < max2);
			if (res <= 0) {
				unsigned temp = temp_inst[++i];
				temp_inst[i] = temp_inst[j];
				temp_inst[j] = temp;
			}
		}

		unsigned temp = temp_inst[++i];
		temp_inst[i] = temp_inst[last];
		temp_inst[last] = temp;

		i += first;
		last += first;

		if (i - 1 > long(first)) {
			stack[++top] = first;
			stack[++top] = i - 1;
		}

		if (i + 1 < last) {
			stack[++top] = i + 1;
			stack[++top] = last;
		}
	}
	delete[] stack;
}

float DSTree::calculateSAH(unsigned& new_axis, unsigned& new_index, unsigned* instances_ptr,
	simpleAABB& parentAABB, simpleAABB& leftAABB, simpleAABB& rightAABB,
	unsigned elem_count, float parentSAH) {
	LiteMath::float4* tempVertices = &(vertices[meshes[meshes.size() - 1].firstVertID]);
	unsigned* tempIndices = &(indices[3u * meshes[meshes.size() - 1].firstIndID]);
	unsigned index = 0u, ind[3] = { 0 };

	unsigned* tempInstances = new unsigned[elem_count];
	new_axis = 3u;

	for (unsigned curr_axis = 0u; curr_axis < 3u; ++curr_axis) {
		std::vector <simpleAABB> right_surfs;

		qsortLower(instances_ptr, elem_count, curr_axis, tempVertices, tempIndices);
		right_surfs.push_back(getAABB(instances_ptr[elem_count - 1u], tempVertices, tempIndices));
		for (unsigned j = elem_count - 2u; j > 0u; --j) {
			right_surfs.push_back(mergeAABB(instances_ptr[j], right_surfs[right_surfs.size() - 1u], tempVertices, tempIndices));
		}
		parentAABB = mergeAABB(*instances_ptr, right_surfs[right_surfs.size() - 1u], tempVertices, tempIndices);

		simpleAABB temp_aabb(getAABB(*instances_ptr, tempVertices, tempIndices));

		for (unsigned k = 1u; k < elem_count - 1u; ++k) {
			temp_aabb = mergeAABB(instances_ptr[k - 1u], temp_aabb, tempVertices, tempIndices);

			float tempSAH = EMPTY_COST + k * calcSurf(temp_aabb) +
				(elem_count - k) * calcSurf(right_surfs[right_surfs.size() - k]);

			if (tempSAH < parentSAH) {
				parentSAH = tempSAH;
				rightAABB = right_surfs[right_surfs.size() - k];
				leftAABB = temp_aabb;
				new_index = k - 1u;
				ind[curr_axis] = k - 1u;
				new_axis = curr_axis;
			}
		}
		if (curr_axis < 2u && new_axis == curr_axis) {
			std::memcpy(tempInstances, instances_ptr, sizeof(unsigned) * elem_count);
		}
		right_surfs.clear();
	}

	if (new_axis < 2u) {
		std::memcpy(instances_ptr, tempInstances, sizeof(unsigned) * elem_count);
	}

	delete[] tempInstances;
	return parentSAH;
}

void DSTree::calculateUpperTreeSAH(unsigned& new_axis, unsigned& new_index, unsigned* instances_ptr,
	simpleAABB& parentAABB, simpleAABB& leftAABB, simpleAABB& rightAABB, unsigned elem_count) {
	unsigned index = 0u;
	new_axis = 3u;
	float parentSAH = FLT_MAX;
	unsigned* tempInstances = new unsigned[elem_count];

	for (unsigned curr_axis = 0u; curr_axis < 3u; ++curr_axis) {
		std::vector <simpleAABB> right_surfs;

		qsortUpper(instances_ptr, elem_count, curr_axis);

		right_surfs.push_back(getInstAABB(instances_ptr[elem_count - 1u]));
		for (unsigned j = elem_count - 2u; j > 0u; --j) {
			right_surfs.push_back(mergeInstAABB(instances_ptr[j], right_surfs[right_surfs.size() - 1u]));
		}
		parentAABB = mergeInstAABB(*instances_ptr, right_surfs[right_surfs.size() - 1u]);

		simpleAABB temp_aabb(getInstAABB(*instances_ptr));
		for (unsigned k = 1u; k < elem_count; ++k) {
			temp_aabb = mergeInstAABB(instances_ptr[k - 1u], temp_aabb);

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
		if (curr_axis < 2u && new_axis == curr_axis) {
			std::memcpy(tempInstances, instances_ptr, sizeof(unsigned) * elem_count);
		}

		right_surfs.clear();
	}

	if (new_axis < 2u) {
		std::memcpy(instances_ptr, tempInstances, sizeof(unsigned) * elem_count);
	}

	delete[] tempInstances;
	return;
}


/* ========  Tree builder  ======== */

unsigned DSTree::buildHeader(bool is_leaf, bool is_carve, unsigned info) {
	unsigned header = unsigned(is_leaf) + (unsigned(is_carve) << 1);
	return (info << 2u) + header;
}

unsigned DSTree::buildCarvingNodes(simpleAABB& parent_aabb, simpleAABB& child_aabb, unsigned axis,
									bool is_left_child, bool is_leaf, unsigned triangle_index) {
	unsigned size1 = dst_nodes.size();
	unsigned sortInd[6u] = { 0u, 1u, 2u, 3u, 4u, 5u };
	unsigned split_plane = axis + 3u * unsigned(is_left_child);
	float child[6u] = { child_aabb.min.x, child_aabb.min.y, child_aabb.min.z, child_aabb.max.x, child_aabb.max.y, child_aabb.max.z };
	float dif[6u] = { child_aabb.min.x - parent_aabb.min.x,
					  child_aabb.min.y - parent_aabb.min.y,
					  child_aabb.min.z - parent_aabb.min.z,
					  parent_aabb.max.x - child_aabb.max.x,
					  parent_aabb.max.y - child_aabb.max.y,
					  parent_aabb.max.z - child_aabb.max.z };

	{
		unsigned temp = sortInd[split_plane];
		sortInd[split_plane] = 5;
		sortInd[5] = temp;
	}
	for (unsigned i = 0u; i <= 5u; ++i) {
		for (unsigned j = i + 1u; j < 5u; ++j) {
			if (dif[sortInd[i]] < dif[sortInd[j]]) {
				unsigned temp = sortInd[i];
				sortInd[i] = sortInd[j];
				sortInd[j] = temp;
			}
		}
		unsigned tempInd = sortInd[i];
		if (dif[tempInd] < 0.01f) sortInd[i] += 6u;
		if (dif[tempInd] <= 0.f) sortInd[i] += 6u;
		if (tempInd == split_plane) sortInd[i] += 6u;
		if (i & 1u) {
			if (sortInd[i - 1] >= 6u || sortInd[i] >= 12u) {
				size1 = (i >> 1);
				break;
			}
			float planes[2];

			unsigned plane1 = sortInd[i - 1u];
			planes[0] = child[plane1];
			bool norm_pos1 = plane1 > 2u;
			if (norm_pos1) plane1 -= 3u;
			
			unsigned plane2 = sortInd[i];
			if (plane2 >= 6u) plane2 -= 6u;
			planes[1] = child[plane2];
			bool norm_pos2 = plane2 > 2u;
			if (norm_pos2) plane2 -= 3u;

			bool single_carve = plane1 == plane2;
			if (plane1 > plane2 || (single_carve && norm_pos1)) {
				unsigned temp1 = plane1;
				plane1 = plane2;
				plane2 = temp1;

				bool temp2 = norm_pos1;
				norm_pos1 = norm_pos2;
				norm_pos2 = temp2;

				float temp3 = planes[0];
				planes[0] = planes[1];
				planes[1] = temp3;
			}

			unsigned info = 8u + plane1;
			if (!single_carve) {
				info = (unsigned(plane1 == 1u) << 3) + (unsigned(plane2 == 2u) << 2) + (unsigned(norm_pos1) << 1) + unsigned(norm_pos2);
			}
			DSNode temp_node{ static_cast<unsigned>(((dst_nodes.size() + 1) << 6) + buildHeader(false, true, info)), 0u, {planes[0], planes[1]} };
			dst_nodes.push_back(temp_node);
		}
	}
	if (is_leaf && size1) {
		dst_nodes[dst_nodes.size() - 1u].leftChild = (dst_nodes[dst_nodes.size() - 1u].leftChild & 62u) + (triangle_index << 6) + 1u;
	}
	return size1;
}

void DSTree::dstBuilderRecur(unsigned first, unsigned last, float parentSAH) {
	unsigned curr_node_index = dst_nodes.size();
	unsigned new_index = last;

	dst_nodes.emplace_back();
	dst_nodes[curr_node_index].leftChild = (first << 6u) + buildHeader(true, false, 0u);
	dst_nodes[curr_node_index].planes[0] = static_cast<float>(last - first + 1u);

	if (last - first >= DSTREE_MAX_POLY) {
		simpleAABB left_aabb, right_aabb, curr_aabb;
		unsigned new_axis = 3u;
		float SAH = calculateSAH(new_axis, new_index, &(instances[first]),
								 curr_aabb, left_aabb, right_aabb, unsigned(last - first + 1u), parentSAH);

		if (new_axis < 3u && SAH <= parentSAH) {
			parentSAH = SAH;
			new_index += first;

			dst_nodes[curr_node_index].planes[0] =  left_aabb.max[new_axis];
			dst_nodes[curr_node_index].planes[1] = right_aabb.min[new_axis];
			dst_nodes[curr_node_index].leftChild = (dst_nodes.size() << 6u) + buildHeader(false, false, new_axis);

			unsigned temp1 = buildCarvingNodes(curr_aabb, left_aabb, new_axis, true, (new_index - first) < DSTREE_MAX_POLY, first);
			if (!((dst_nodes[dst_nodes.size() - 1].leftChild & 1u) && temp1))
				dstBuilderRecur(first, new_index, parentSAH);

			for (unsigned temp_index = curr_node_index + 1u;;) {
				dst_nodes[temp_index].rightNode = dst_nodes.size();
				if (dst_nodes[temp_index].leftChild & 1u)
					break;

				if (dst_nodes[temp_index].leftChild & 2u) {
					++temp_index;
				}
				else if (!(dst_nodes[temp_index].leftChild & 3u))
					temp_index = dst_nodes[dst_nodes[temp_index].leftChild >> 6].rightNode;
			}

			unsigned carve_nodes_num = buildCarvingNodes(curr_aabb, right_aabb, new_axis, false, (last - new_index) <= DSTREE_MAX_POLY, new_index + 1u);
			if (!((dst_nodes[dst_nodes.size() - 1].leftChild & 1u) && carve_nodes_num))
				dstBuilderRecur(new_index + 1u, last, parentSAH);
		}
	}
	return;
}

void DSTree::dstBuilderRecurUpper(unsigned first, unsigned last) {
	unsigned curr_node_index = dst_nodes.size();
	unsigned new_index = last;

	dst_nodes.emplace_back();
	dst_nodes[curr_node_index].leftChild = (instances[first] << 6u) + buildHeader(true, false, 0u);
	dst_nodes[curr_node_index].planes[0] = static_cast<float>(last - first + 1u);

	if (last - first) {
		simpleAABB left_aabb, right_aabb, curr_aabb;
		unsigned new_axis = 3u;
		calculateUpperTreeSAH(new_axis, new_index, &(instances[first]),
			curr_aabb, left_aabb, right_aabb, unsigned(last - first + 1u));

		if (new_axis < 3u) {
			new_index += first;
			dst_nodes[curr_node_index].planes[0] = left_aabb.max[new_axis];
			dst_nodes[curr_node_index].planes[1] = right_aabb.min[new_axis];
			dst_nodes[curr_node_index].leftChild = (dst_nodes.size() << 6u) + buildHeader(false, false, new_axis);

			buildCarvingNodes(curr_aabb, left_aabb, new_axis, true, new_index == first, instances[first]);
			if (!(new_index == first && (dst_nodes.size() - 1u != curr_node_index)))
				dstBuilderRecurUpper(first, new_index);

			for (unsigned temp_index = curr_node_index + 1u;;) {
				dst_nodes[temp_index].rightNode = dst_nodes.size();
				if (dst_nodes[temp_index].leftChild & 1u)
					break;

				if (dst_nodes[temp_index].leftChild & 2u) {
					++temp_index;
				}
				else if (!(dst_nodes[temp_index].leftChild & 3u))
					temp_index = dst_nodes[dst_nodes[temp_index].leftChild >> 6].rightNode;
			}

			unsigned carve_nodes_num = buildCarvingNodes(curr_aabb, right_aabb, new_axis, false, !bool(last - new_index - 1u), instances[new_index + 1u]);
			if (!(last == new_index + 1u && carve_nodes_num))
				dstBuilderRecurUpper(new_index + 1u, last);
		}
	}
	return;
}

void DSTree::CommitScene() {
	if (indices.empty())
		return;

	size_t instances_info_size = instances_info.size();
	//std::cout << "Triangles: " << indices.size() / 3 << std::endl;
	//std::cout << "Instances: " << instances_info_size << std::endl;
	/*std::fstream dstFile;
	dstFile.open("./dstree", std::fstream::in);
	
	if (dstFile.is_open() && !dstFile.eof()) {
		std::cout << "There is a file" << std::endl;
		dstFile >> scene_AABB.min.x;
		dstFile >> scene_AABB.min.y;
		dstFile >> scene_AABB.min.z;
		dstFile >> scene_AABB.max.x;
		dstFile >> scene_AABB.max.y;
		dstFile >> scene_AABB.max.z;

		uint32_t dst_size;
		dstFile >> dst_size;

		DSNode tempNode;
		for (uint32_t j = 0u; j < dst_size; ++j) {
			dstFile >> tempNode.leftChild;
			dstFile >> tempNode.rightNode;
			dstFile >> tempNode.planes[0];
			dstFile >> tempNode.planes[1];

			if (!dstFile.eof()) dst_nodes.push_back(tempNode);
		}

		for (uint32_t k = 0u; k < instances_size; ++k) {
			dstFile >> instances[k];
		}
	}
	if (!dst_nodes.size()) {
		std::cout << "The file was empty" << std::endl;

		dstFile.open("./dstree", std::fstream::out | std::fstream::trunc);*/

	/*  ========  Upper Tree  ========  */
	dst_nodes.clear();
	instances.clear();

	for (uint32_t i = 0u; i < instances_info_size; ++i)
		instances.push_back(instances.size());

	dstBuilderRecurUpper(0u, instances_info_size - 1u);

	upper_tree = dst_nodes;
	dst_nodes.clear();
	if (instances_info_size)
		scene_AABB = instances_info[0].instAABB;
	for (uint32_t i = 1u; i < instances_info_size; ++i) {
		scene_AABB = mergeInstAABB(i, scene_AABB);
	}
	/*
	size_t trNum = instances_info_size;
	for (uint32_t i = 0u; i < instances_info_size; ++i)
		trNum += (meshes[instances_info[i].geomID].lastIndID - meshes[instances_info[i].geomID].firstIndID);
	std::cout << "Triangles: " << trNum << std::endl;
	std::cout << "Upper tree: " << upper_tree.size() << "; Lower tree: " << lower_tree.size() << std::endl;
	std::cout << "AABB min: " << scene_AABB.min.x << ", " << scene_AABB.min.y << ", " << scene_AABB.min.z << std::endl;
	std::cout << "AABB max: " << scene_AABB.max.x << ", " << scene_AABB.max.y << ", " << scene_AABB.max.z << std::endl;*/

		/*
		dstFile << std::setprecision(16) << scene_AABB.min.x << " ";
		dstFile << std::setprecision(16) << scene_AABB.min.y << " ";
		dstFile << std::setprecision(16) << scene_AABB.min.z << " ";
		dstFile << std::setprecision(16) << scene_AABB.max.x << " ";
		dstFile << std::setprecision(16) << scene_AABB.max.y << " ";
		dstFile << std::setprecision(16) << scene_AABB.max.z << " \n\n";
		dstFile << dst_nodes.size() << " \n\n\n";

		for (uint32_t i = 0u; i < dst_nodes.size(); ++i) {
			dstFile << dst_nodes[i].leftChild << " ";
			dstFile << dst_nodes[i].rightNode << " ";
			dstFile << std::setprecision(16) << dst_nodes[i].planes[0] << " ";
			dstFile << std::setprecision(16) << dst_nodes[i].planes[1] << " \n\n";
		}

		for (uint32_t j = 0u; j < instances.size(); ++j) {
			dstFile << instances[j] << " \n";
		}
	}
	dstFile.close();
	std::cout << "Nodes: " << dst_nodes.size() << std::endl;*/
	return;
}


/* ========  Ray traverse  ======== */

bool DSTree::traceAABB(const simpleAABB &tempAABB, LiteMath::float3 Position, LiteMath::float3 Direction, LiteMath::float3 InvDir, LiteMath::float2 &tMinMax) {
	float tMin = LiteMath::INF_NEGATIVE, tMax = LiteMath::INF_POSITIVE;
	LiteMath::float3 Min = tempAABB.min, Max = tempAABB.max;

	for (unsigned dim = 0u; dim < 3u; ++dim) {
		if (Direction[dim] > -LiteMath::EPSILON && Direction[dim] < LiteMath::EPSILON) {
			if (Position[dim] < Min[dim] || Position[dim] > Max[dim])
				return false;

			continue;
		}
		float t1 = (Min[dim] /*- AABBeps*/ - Position[dim]) * InvDir[dim];
		float t2 = (Max[dim] /*+ AABBeps*/ - Position[dim]) * InvDir[dim];

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
	tMinMax = LiteMath::float2(max(tMin, 0.f), tMax);
	return tMin < tMax && tMax >= 0.0f;
}

CRT_Hit DSTree::traceTriangle(LiteMath::float3 Position, LiteMath::float3 Direction, LiteMath::float4* tempVertices, unsigned* tempIndices, uint32_t current_instance) {
	CRT_Hit tempHitInfo{ FLT_MAX, static_cast<unsigned>(-1), static_cast<unsigned>(-1), static_cast<unsigned>(-1), {0.f} };

	const LiteMath::float3  a = LiteMath::to_float3(tempVertices[tempIndices[0]]);
	const LiteMath::float3 E1 = LiteMath::to_float3(tempVertices[tempIndices[1]]) - a;
	const LiteMath::float3 E2 = LiteMath::to_float3(tempVertices[tempIndices[2]]) - a;

	const LiteMath::float3  P = cross(Direction, E2);
	const float d = dot(E1, P);

	if (d < -LiteMath::EPSILON || d > LiteMath::EPSILON) {
		const LiteMath::float3 T = Position - a;
		const float f = 1.0f / d;
		const float temp_u = f * dot(P, T);
		const LiteMath::float3 Q = cross(T, E1);
		const float temp_v = f * dot(Direction, Q);
		const float temp_t = f * dot(E2, Q);

		if (temp_u >= 0.0f && temp_u <= 1.0f && ((temp_u + temp_v) <= 1.0f + AABBeps) &&
			temp_v >= 0.0f && temp_t >= 0.0f) {

			tempHitInfo.coords[0] = temp_u;
			tempHitInfo.coords[1] = temp_v;
			tempHitInfo.t = temp_t;
		}
	}
	return tempHitInfo;
}

CRT_Hit DSTree::RayQuery_NearestHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar) {
	std::vector <uint32_t> insts;
	findInstHit(posAndNear, dirAndFar, insts);
	//for (uint32_t i = 0u; i < instances_info.size(); ++i)
	//	insts.push_back(i);

	CRT_Hit temp_hit{ FLT_MAX, static_cast<unsigned>(-1), static_cast<unsigned>(-1), static_cast<unsigned>(-1), {0.f} };

	for (uint32_t i = 0u; i < insts.size(); ++i) {
		//current_instance = insts[i];
		//if (dirAndFar.w > temp_hit.t) dirAndFar.w = temp_hit.t;
		CRT_Hit instHit = findHit(posAndNear, dirAndFar, false, insts[i], temp_hit.t);
		//temp_hit.coords[2] += instHit.coords[2];
		//temp_hit.coords[3] += instHit.coords[3];
		//instHit.coords[2] = temp_hit.coords[2];
		//instHit.coords[3] = temp_hit.coords[3];
		if (instHit.t < temp_hit.t)
			temp_hit = instHit;
	}

	return temp_hit;
}

bool DSTree::RayQuery_AnyHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar) {
	std::vector <uint32_t> insts;
	findInstHit(posAndNear, dirAndFar, insts);

	CRT_Hit temp_hit{ FLT_MAX, static_cast<unsigned>(-1), static_cast<unsigned>(-1), static_cast<unsigned>(-1), {0.f} };

	for (uint32_t i = 0u; i < insts.size(); ++i) {
		//current_instance = insts[i];

		//simpleMeshInfo tempMesh = meshes[instances_info[current_instance].geomID];
		//tempVertices = &(vertices[tempMesh.firstVertID]);
		//tempIndices = &(indices[3 * tempMesh.firstIndID]);
		//tempIndicesSorted = &(indices_sorted[tempMesh.firstIndID]);
		//dst_ptr = &(lower_tree[tempMesh.DSTreeOffset]);

		// NOT DONE YET
		//CRT_Hit instHit = findHit(posAndNear, dirAndFar, true, insts[i]);
		//if (instHit.t < temp_hit.t)
		//	temp_hit = instHit;
	}

	return (temp_hit.t >= posAndNear.w) && (temp_hit.t <= dirAndFar.w);
}

std::vector <unsigned> DSTree::TreePath(unsigned init, unsigned find) {
	std::vector <unsigned> tempRoute;
	for (unsigned tempInit = init; tempInit < find;) {
		tempRoute.push_back(tempInit);
		if (lower_tree[tempInit].rightNode > 0 && lower_tree[tempInit].rightNode < find)
			tempInit = lower_tree[tempInit].rightNode;
		else ++tempInit;
	}
	return tempRoute;
}

struct travStack {
	uint32_t node{ 0u };
	LiteMath::float2 t{ 0.f, 0.f };
};

void DSTree::findInstHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, std::vector <uint32_t>& insts) {
	LiteMath::float3 position(posAndNear.x, posAndNear.y, posAndNear.z);
	LiteMath::float3 direction(dirAndFar.x, dirAndFar.y, dirAndFar.z);
	LiteMath::float3 invDir(1.f / dirAndFar.x, 1.f / dirAndFar.y, 1.f / dirAndFar.z);
	bool isFin[3] = { std::isfinite(invDir.x), std::isfinite(invDir.y), std::isfinite(invDir.z) };

	travStack travInfo{ 0u, { 0.f, FLT_MAX } };

	if (!traceAABB(scene_AABB, position, direction, invDir, travInfo.t)) return;

	//std::vector <travStack> param_stack;
	travStack param_stack[64];
	int stack_ind = -1;
	do {
		DSNode curr_node;
		bool CHECK_LEAF_POLY(false), IS_CARVING_NODE(false), IS_DOUBLE_CARVE(false);

		do {
			curr_node = upper_tree[travInfo.node];

			unsigned header = curr_node.leftChild & 63u;
			curr_node.leftChild = curr_node.leftChild >> 6;
			CHECK_LEAF_POLY = header & 1u;
			IS_CARVING_NODE = header & 2u;
			IS_DOUBLE_CARVE = !((header & 48u) == 32u) && IS_CARVING_NODE;

			if (!CHECK_LEAF_POLY || IS_CARVING_NODE) {
				LiteMath::uint2 planes = LiteMath::uint2((header >> 2) & 3u);
				bool is_normal_positive[2] = { !IS_CARVING_NODE, IS_CARVING_NODE };

				if (IS_DOUBLE_CARVE) {
					is_normal_positive[0u] = bool(planes.x & 2u);
					is_normal_positive[1u] = bool(planes.x & 1u);
					planes = LiteMath::uint2((header >> 5) & 1u, ((header >> 4) & 1u) + 1u);
				}

				float t[2] = { -1.f, -1.f };

				if (isFin[planes.x]) t[0] = invDir[planes.x] * (curr_node.planes[0u] - position[planes.x]);
				if (isFin[planes.y]) t[1] = invDir[planes.y] * (curr_node.planes[1u] - position[planes.y]);

				if (!IS_CARVING_NODE) {
					bool is_min_left = direction[planes.x] <= 0.f;
					travStack t_info[2] = { travInfo, travInfo };

					t_info[0].node = curr_node.leftChild;
					t_info[1].node = upper_tree[curr_node.leftChild].rightNode;

					bool trav_correct[2] = { true, true };
					for (unsigned i = 0u; i < 2u; ++i) {
						if (t[i] > 0.f) {
							float temp_t[2] = { travInfo.t.x, travInfo.t.y };
							temp_t[is_min_left] = t[i];

							if (is_min_left) t_info[i].t.x = temp_t[temp_t[0] < temp_t[1]];
							if (!is_min_left) t_info[i].t.y = temp_t[temp_t[0] > temp_t[1]];
						}

						if (t_info[i].t.x > t_info[i].t.y || t_info[i].t.y < 0.f)
							trav_correct[i] = false;

						is_min_left = !is_min_left;
					}

					if (trav_correct[0] && trav_correct[1]) {
						travStack tempInfo = t_info[1];
						travInfo = t_info[0];

						if (tempInfo.t.x < travInfo.t.x) {
							tempInfo = t_info[0];
							travInfo = t_info[1];
						}

						//param_stack.push_back(tempInfo);
						param_stack[++stack_ind] = tempInfo;
					}
					if (trav_correct[0] != trav_correct[1]) {
						if (trav_correct[0]) travInfo = t_info[0];
						if (trav_correct[1]) travInfo = t_info[1];
					}
					if (!(trav_correct[0] || trav_correct[1])) {
						travInfo.node = 0u;

						if (stack_ind >= 0)
							travInfo = param_stack[stack_ind--];
						/*
						if (!param_stack.empty()) {
							travInfo = *param_stack.rbegin();
							param_stack.pop_back();
						}*/
					}
				}

				if (IS_CARVING_NODE) {
					bool is_min[2] = { is_normal_positive[0] != (direction[planes.x] > 0.f),
									   is_normal_positive[1] != (direction[planes.y] > 0.f) };

					if (is_min[0] == is_min[1]) {
						if (is_min[0]) {
							float temp_min = t[t[0] < t[1]];
							travInfo.t.x = max(travInfo.t.x, temp_min);
						}
						if (!is_min[0]) {
							float temp_max = t[t[0] > t[1]];
							travInfo.t.y = min(travInfo.t.y, temp_max);
						}
					}

					if (is_min[0] != is_min[1]) {
						if (!is_min[0]) {
							float temp = t[0];
							t[0] = t[1];
							t[1] = temp;
						}

						travInfo.t.x = max(travInfo.t.x, t[0]);
						travInfo.t.y = min(travInfo.t.y, t[1]);
					}

					bool skip = travInfo.t.x > travInfo.t.y || travInfo.t.y < 0.f;

					travInfo.node = curr_node.leftChild;
					if (skip || CHECK_LEAF_POLY) {
						travInfo.node = 0u;

						if (stack_ind >= 0)
							travInfo = param_stack[stack_ind--];
						/*
						if (!param_stack.empty()) {
							travInfo = *param_stack.rbegin();
							param_stack.pop_back();
						}*/
					}

					CHECK_LEAF_POLY = CHECK_LEAF_POLY && !skip;
				}
			}

		} while (!CHECK_LEAF_POLY && travInfo.node != 0u);

		if (CHECK_LEAF_POLY) {
			if (!IS_CARVING_NODE) {
				travInfo.node = 0u;

				if (stack_ind >= 0)
					travInfo = param_stack[stack_ind--];
				/*
				if (!param_stack.empty()) {
					travInfo = *param_stack.rbegin();
					param_stack.pop_back();
				}*/
			}

			insts.push_back(curr_node.leftChild);
		}

	} while (travInfo.node != 0u);

	return;
}

CRT_Hit DSTree::findHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, bool findAny, uint32_t current_instance, float triang_t) {
	simpleInstance tempInst = instances_info[current_instance];
	simpleMeshInfo tempMesh = meshes[tempInst.geomID];
	unsigned tree_tr_num = tempMesh.lastIndID - tempMesh.firstIndID + 1;	//used once

	LiteMath::float3 position = LiteMath::inverse4x4(tempInst.transf) * LiteMath::float3(posAndNear.x, posAndNear.y, posAndNear.z);
	LiteMath::float4 temp_direction = LiteMath::mul4x4x4(LiteMath::inverse4x4(tempInst.transf), LiteMath::float4(dirAndFar.x, dirAndFar.y, dirAndFar.z, 0.f));
	LiteMath::float3 direction = LiteMath::float3(temp_direction.x, temp_direction.y, temp_direction.z);
	LiteMath::float3 invDir = LiteMath::float3(1.f / direction.x, 1.f / direction.y, 1.f / direction.z);
	bool isFin[3] = { std::isfinite(invDir.x), std::isfinite(invDir.y), std::isfinite(invDir.z) };

	travStack travInfo{ 0u, { 0.f, FLT_MAX } };
	CRT_Hit temp_hit{ FLT_MAX, static_cast<unsigned>(-1), static_cast<unsigned>(-1), static_cast<unsigned>(-1), {0.f} };

	if (!traceAABB(tempMesh.meshAABB, position, direction, invDir, travInfo.t)) return temp_hit;
	if (triang_t < travInfo.t.x) return temp_hit;
	travInfo.t.y = min(travInfo.t.y, triang_t);

	LiteMath::float4* tempVertices = &(vertices[tempMesh.firstVertID]);
	unsigned* tempIndices = &(indices[3u * tempMesh.firstIndID]);
	unsigned* tempIndicesSorted = &(indices_sorted[tempMesh.firstIndID]);
	DSNode* dst_ptr = &(lower_tree[tempMesh.DSTreeOffset]);
	//std::vector <travStack> param_stack;
	travStack param_stack[64];
	int stack_ind = -1;

	do {
		DSNode curr_node;
		bool CHECK_LEAF_POLY(false), IS_CARVING_NODE(false), IS_DOUBLE_CARVE(false);

		do {
			curr_node = dst_ptr[travInfo.node];

			unsigned header = curr_node.leftChild & 63u;
			curr_node.leftChild = curr_node.leftChild >> 6;
			CHECK_LEAF_POLY = header & 1u;
			IS_CARVING_NODE = header & 2u;
			IS_DOUBLE_CARVE = !((header & 48u) == 32u) && IS_CARVING_NODE;

			if (!CHECK_LEAF_POLY || IS_CARVING_NODE) {
				LiteMath::uint2 planes = LiteMath::uint2((header >> 2) & 3u);
				bool is_normal_positive[2] = { !IS_CARVING_NODE, IS_CARVING_NODE };

				if (IS_DOUBLE_CARVE) {
					is_normal_positive[0u] = bool(planes.x & 2u);
					is_normal_positive[1u] = bool(planes.x & 1u);
					planes = LiteMath::uint2((header >> 5) & 1u, ((header >> 4) & 1u) + 1u);
				}

				float t[2] = { -1.f, -1.f };

				if (isFin[planes.x]) t[0] = invDir[planes.x] * (curr_node.planes[0u] - position[planes.x]);
				if (isFin[planes.y]) t[1] = invDir[planes.y] * (curr_node.planes[1u] - position[planes.y]);

				if (!IS_CARVING_NODE) {
					bool is_min_left = direction[planes.x] < 0.f;
					travStack t_info[2] = { travInfo, travInfo };

					t_info[0].node = curr_node.leftChild;
					t_info[1].node = dst_ptr[curr_node.leftChild].rightNode;

					bool trav_correct[2] = { true, true };
					for (unsigned i = 0u; i < 2u; ++i) {
						if (t[i] > 0.f) {
							float temp_t[2] = { travInfo.t.x, travInfo.t.y };
							temp_t[is_min_left] = t[i];

							if (is_min_left) t_info[i].t.x = temp_t[temp_t[0] < temp_t[1]];
							if (!is_min_left) t_info[i].t.y = temp_t[temp_t[0] > temp_t[1]];
						}

						if (t_info[i].t.x > t_info[i].t.y || t_info[i].t.y < 0.f)
							trav_correct[i] = false;

						is_min_left = !is_min_left;
					}

					if (trav_correct[0] && trav_correct[1]) {
						travStack tempInfo = t_info[1];
						travInfo = t_info[0];

						if (tempInfo.t.x < travInfo.t.x) {
							tempInfo = t_info[0];
							travInfo = t_info[1];
						}

						if (temp_hit.t >= tempInfo.t.x)
							//param_stack.push_back(tempInfo);
							param_stack[++stack_ind] = tempInfo;
					}
					if (trav_correct[0] != trav_correct[1]) {
						if (trav_correct[0]) travInfo = t_info[0];
						if (trav_correct[1]) travInfo = t_info[1];
					}
					if (!(trav_correct[0] || trav_correct[1])) {
						travInfo.node = 0u;

						if (stack_ind >= 0)
							travInfo = param_stack[stack_ind--];
						/*
						if (!param_stack.empty()) {
							travInfo = *param_stack.rbegin();
							param_stack.pop_back();
						}*/
					}
				}

				if (IS_CARVING_NODE) {
					bool is_min[2] = { is_normal_positive[0] != (direction[planes.x] > 0.f),
									   is_normal_positive[1] != (direction[planes.y] > 0.f) };

					travInfo.node = curr_node.leftChild;
					if (is_min[0] == is_min[1]) {
						if (is_min[0]) {
							float temp_min = t[t[0] < t[1]];
							travInfo.t.x = max(travInfo.t.x, temp_min);
						}
						if (!is_min[0]) {
							float temp_max = t[t[0] > t[1]];
							travInfo.t.y = min(travInfo.t.y, temp_max);
						}
					}

					if (is_min[0] != is_min[1]) {
						if (!is_min[0]) {
							float temp = t[0];
							t[0] = t[1];
							t[1] = temp;
						}

						travInfo.t.x = max(travInfo.t.x, t[0]);
						travInfo.t.y = min(travInfo.t.y, t[1]);
					}

					bool skip = travInfo.t.x > travInfo.t.y || travInfo.t.y < 0.f;

					if (skip || CHECK_LEAF_POLY) {
						travInfo.node = 0u;

						if (stack_ind >= 0)
							travInfo = param_stack[stack_ind--];
						/*
						if (!param_stack.empty()) {
							travInfo = *param_stack.rbegin();
							param_stack.pop_back();
						}*/
					}

					CHECK_LEAF_POLY = CHECK_LEAF_POLY && !skip;
				}
			}

		} while (!CHECK_LEAF_POLY && travInfo.node != 0u);

		if (CHECK_LEAF_POLY) {
			unsigned max_node_tr_index = DSTREE_MAX_POLY;

			if (!IS_CARVING_NODE) {
				max_node_tr_index = static_cast<unsigned>(curr_node.planes[0]);
				travInfo.node = 0u;

				if (stack_ind >= 0)
					travInfo = param_stack[stack_ind--];
				/*
				if (!param_stack.empty()) {
					travInfo = *param_stack.rbegin();
					param_stack.pop_back();
				}*/
			}

			max_node_tr_index += curr_node.leftChild;
			if (max_node_tr_index > tree_tr_num) max_node_tr_index = tree_tr_num;

			for (unsigned j = curr_node.leftChild; j < max_node_tr_index; ++j) {
				CRT_Hit temp_t = traceTriangle(position, direction, tempVertices, &(tempIndices[3 * tempIndicesSorted[j]]), current_instance);
				if (temp_t.t < temp_hit.t && temp_t.t >= 0.f) {
					temp_hit = temp_t;
					temp_hit.instId = current_instance;
					temp_hit.geomId = tempInst.geomID;
					temp_hit.primId = tempIndicesSorted[j];
					if (findAny) return temp_hit;
				}
			}
		}

	} while (travInfo.node != 0u);

	return temp_hit;
}


ISceneObject* CreateSceneRT(const char* a_impleName) { return new DSTree; }

void DeleteSceneRT(ISceneObject* a_pScene) { delete a_pScene; }