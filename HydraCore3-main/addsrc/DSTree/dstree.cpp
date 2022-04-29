#include "dstree.h"


DSTree* scene;

unsigned DSTree::AddGeom_Triangles4f(const LiteMath::float4* a_vpos4f, size_t a_vertNumber, const unsigned* a_triIndices, size_t a_indNumber) {
	simpleMeshInfo tempInfo;
	tempInfo.firstVertID  = -1;
	tempInfo.firstIndID   = -1;
	tempInfo.lastVertID   = -1;
	tempInfo.lastIndID    = -1;
	tempInfo.DSTreeOffset = -1;

	instances.resize(a_indNumber / 3);
	for (size_t j = 0; j < a_indNumber / 3; ++j) {
		instances[j] = j;/*tempInfo.firstVertID + a_triIndices[j];*/
		//std::cout << "Index " << j << std::endl;
	}

	if (a_vertNumber && a_indNumber) {
		tempVertices = a_vpos4f;
		tempIndices = a_triIndices;
		//getElemFunc = getMaxPos; //by default all function pointers are set to build lower_tree
		dst_nodes.clear();

		std::cout << "Mesh " << meshes.size() << "; Triangles: " << instances.size() << std::endl;
		std::cout << "Indices " << a_indNumber << "; Vertices: " << a_vertNumber << std::endl;
		scene = this;
		dstBuilderRecur(0u, instances.size() - 1, SAH_MAX);
		scene = nullptr;

		//for (size_t j = 0; j < a_indNumber; ++j)
		//	indices.push_back(instances[j]);

		tempInfo.DSTreeOffset = lower_tree.size();
		lower_tree.insert(lower_tree.end(), dst_nodes.begin(), dst_nodes.end());
		dst_nodes.clear();

		simpleAABB tempAABB = getAABB(0);
		for (uint32_t j = 1; j < instances.size(); ++j) {
			tempAABB = mergeAABB(j, tempAABB);
		}
		tempInfo.meshAABB = tempAABB;
		for (int i = 0; i < 3; ++i) {
			if (tempInfo.meshAABB.min[i] == tempInfo.meshAABB.max[i]) {
				tempInfo.meshAABB.min[i] -= AABBeps;
				tempInfo.meshAABB.max[i] += AABBeps;
			}
		}
	}
	if (a_vertNumber) {
		tempInfo.firstVertID = vertices.size();

		for (size_t j = 0; j < a_vertNumber; ++j) {
			vertices.push_back(a_vpos4f[j]);
			//if (!meshes.size()) std::cout << "Vert " << j << ": " << a_vpos4f[j].x << ", " << a_vpos4f[j].y << ", " << a_vpos4f[j].z << ", " << a_vpos4f[j].w << std::endl;
		}
		//if (!meshes.size()) {
		//	std::cout << "AABB min: " << tempInfo.meshAABB.min.x << ", " << tempInfo.meshAABB.min.y << ", " << tempInfo.meshAABB.min.z << std::endl;
		//	std::cout << "AABB max: " << tempInfo.meshAABB.max.x << ", " << tempInfo.meshAABB.max.y << ", " << tempInfo.meshAABB.max.z << std::endl;
		//}
		tempInfo.lastVertID = vertices.size() - 1;
	}
	if (a_indNumber) {
		tempInfo.firstIndID = indices.size() / 3;

		for (size_t j = 0; j < instances.size(); ++j) {
			indices.push_back(a_triIndices[3 * j    ]);
			indices.push_back(a_triIndices[3 * j + 1]);
			indices.push_back(a_triIndices[3 * j + 2]);

			indices_sorted.push_back(instances[j]);
		}

		tempInfo.lastIndID = indices.size() / 3 - 1;
	}
	meshes.push_back(tempInfo);
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

	for (uint32_t i = 0u; i < 3; ++i) {
		if (tempInst.instAABB.min[i] > tempInst.instAABB.max[i]) {
			float temp = tempInst.instAABB.min[i];
			tempInst.instAABB.min[i] = tempInst.instAABB.max[i];
			tempInst.instAABB.max[i] = temp;
		}
	}

	instances_info.push_back(tempInst);
	//std::cout << "Geom: " << a_geomId << std::endl;
	//std::cout << a_matrix.col(0)[0] << ", " << a_matrix.col(1)[0] << ", " << a_matrix.col(2)[0] << ", " << a_matrix.col(3)[0] << std::endl;
	//std::cout << a_matrix.col(0)[1] << ", " << a_matrix.col(1)[1] << ", " << a_matrix.col(2)[1] << ", " << a_matrix.col(3)[1] << std::endl;
	//std::cout << a_matrix.col(0)[2] << ", " << a_matrix.col(1)[2] << ", " << a_matrix.col(2)[2] << ", " << a_matrix.col(3)[2] << std::endl;
	//std::cout << a_matrix.col(0)[3] << ", " << a_matrix.col(1)[3] << ", " << a_matrix.col(2)[3] << ", " << a_matrix.col(3)[3] << std::endl;
	return instances_info.size() - 1;
}

void DSTree::UpdateInstance(uint32_t a_instanceId, const LiteMath::float4x4& a_matrix) {
	instances_info[a_instanceId].transf = a_matrix;
	instances_info[a_instanceId].instAABB.min = a_matrix * meshes[instances_info[a_instanceId].geomID].meshAABB.min;
	instances_info[a_instanceId].instAABB.max = a_matrix * meshes[instances_info[a_instanceId].geomID].meshAABB.max;
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


simpleAABB DSTree::getAABB(unsigned index) {
	simpleAABB tempAABB;
	LiteMath::float4 a = tempVertices[tempIndices[3 * index    ]];
	LiteMath::float4 b = tempVertices[tempIndices[3 * index + 1]];
	LiteMath::float4 c = tempVertices[tempIndices[3 * index + 2]];

	for (unsigned i = 0u; i < 3u; ++i) {
		tempAABB.min[i] = min(LiteMath::float2(a[i], min(LiteMath::float2(b[i], c[i]))));
		tempAABB.max[i] = max(LiteMath::float2(a[i], max(LiteMath::float2(b[i], c[i]))));
	}
	return tempAABB;
}

simpleAABB DSTree::mergeAABB(unsigned index, simpleAABB aabb) {
	LiteMath::float4 a = tempVertices[tempIndices[3 * index    ]];
	LiteMath::float4 b = tempVertices[tempIndices[3 * index + 1]];
	LiteMath::float4 c = tempVertices[tempIndices[3 * index + 2]];

	for (unsigned i = 0u; i < 3u; ++i) {
		aabb.min[i] = min(LiteMath::float2(min(LiteMath::float2(aabb.min[i], a[i])), min(LiteMath::float2(b[i], c[i]))));
		aabb.max[i] = max(LiteMath::float2(max(LiteMath::float2(aabb.max[i], a[i])), max(LiteMath::float2(b[i], c[i]))));
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

float DSTree::getMinPos(const void* index) {
	unsigned temp_index = *((unsigned*) index) * 3u;
	LiteMath::float3 coords = LiteMath::float3(tempVertices[tempIndices[temp_index     ]][current_scene_axis],
											   tempVertices[tempIndices[temp_index + 1u]][current_scene_axis],
											   tempVertices[tempIndices[temp_index + 2u]][current_scene_axis]);
	return min(LiteMath::float2(coords.x, min(LiteMath::float2(coords.y, coords.z))));
}

float DSTree::getMaxPos(const void* index) {
	unsigned temp_index = *((unsigned*) index) * 3u;
	DSTree* temp = this;
	LiteMath::float3 coords = LiteMath::float3(tempVertices[tempIndices[temp_index     ]][current_scene_axis],
											   tempVertices[tempIndices[temp_index + 1u]][current_scene_axis],
											   tempVertices[tempIndices[temp_index + 2u]][current_scene_axis]);

	return max(LiteMath::float2(coords.x, max(LiteMath::float2(coords.y, coords.z))));
}

float DSTree::getMinAABBpos(const void* index) {
	return instances_info[*((unsigned*)index)].instAABB.min[current_scene_axis];
}

float DSTree::getMaxAABBpos(const void* index) {
	return instances_info[*((unsigned*)index)].instAABB.max[current_scene_axis];
}

int cmpElemPos(const void* ind1, const void* ind2) {
	float max1 = scene->getMaxPos(ind1);
	float max2 = scene->getMaxPos(ind2);
	int res = int(max1 > max2) - int(max1 < max2);

	if (!res) res = int(ind1 > ind2) - int(ind1 < ind2);
	return res;
}

int cmpElemPosAABB(const void* ind1, const void* ind2) {
	float max1 = scene->getMaxAABBpos(ind1);
	float max2 = scene->getMaxAABBpos(ind2);
	int res = int(max1 > max2) - int(max1 < max2);

	if (!res) res = int(ind1 > ind2) - int(ind1 < ind2);
	return res;
}

float DSTree::calculateSAH(unsigned& new_axis, unsigned& new_index, unsigned *instances_ptr,
						   simpleAABB& parentAABB, simpleAABB& leftAABB, simpleAABB& rightAABB,
						   unsigned elem_count, float parentSAH) {
	unsigned index = 0u, ind[3] = { 0 };
	new_axis = 3u;
	//std::vector <unsigned> tempInds[4];
	std::vector <unsigned> tempInstances;
	tempInstances.resize(elem_count);
	//simpleAABB temps[4];
	for (unsigned curr_axis = 0u; curr_axis < 3u; ++curr_axis) {
		std::vector <simpleAABB> right_surfs;
		current_scene_axis = curr_axis;
		qsort(instances_ptr, elem_count, sizeof(unsigned), cmpElemPos);

		right_surfs.push_back(getAABB(instances_ptr[elem_count - 1u]));
		for (unsigned j = elem_count - 2u; j > 0u; --j) {
			right_surfs.push_back(mergeAABB(instances_ptr[j], right_surfs[right_surfs.size() - 1u]));
		}
		parentAABB = mergeAABB(*instances_ptr, right_surfs[right_surfs.size() - 1u]);

		simpleAABB temp_aabb(getAABB(*instances_ptr));
		//tempInds[curr_axis].push_back(*instances_ptr);
		for (unsigned k = 1u; k < elem_count - 1u; ++k) {
			temp_aabb = mergeAABB(instances_ptr[k - 1u], temp_aabb);

			//tempInds[curr_axis].push_back(instances_ptr[k-1]);
			float tempSAH = EMPTY_COST + k * calcSurf(temp_aabb) +
				(elem_count - k) * calcSurf(right_surfs[right_surfs.size() - k]);

			if (tempSAH < parentSAH) {
				parentSAH = tempSAH;
				rightAABB = right_surfs[right_surfs.size() - k];
				leftAABB = temp_aabb;
				//temps[curr_axis] = temp_aabb;
				new_index = k - 1u;
				ind[curr_axis] = k - 1u;
				new_axis = curr_axis;
			}
		}
		if (new_axis == curr_axis) {
			for (unsigned i = 0u; i < elem_count; ++i)
				tempInstances[i] = instances_ptr[i];
		}
		right_surfs.clear();
	}

	if (new_axis < 2u) {
		for (unsigned i = 0u; i < elem_count; ++i)
			instances_ptr[i] = tempInstances[i];
		//current_scene_axis = new_axis;
		//qsort(instances_ptr, elem_count, sizeof(unsigned), cmpElemPos);
	}

	//std::cout << "axis: " << new_axis << std::endl;
	/*tempInds[3].push_back(*instances_ptr);
	temps[3] = getAABB(*instances_ptr);
	for (unsigned k = 1u; k <= new_index; ++k) {
		temps[3] = mergeAABB(instances_ptr[k], temps[3]);
		tempInds[3].push_back(instances_ptr[k]);
	}
	
	if ((fabs(temps[3].min.x - temps[new_axis].min.x) > 0.001f) ||
		(fabs(temps[3].min.y - temps[new_axis].min.y) > 0.001f) ||
		(fabs(temps[3].min.z - temps[new_axis].min.z) > 0.001f) ||
		(fabs(temps[3].max.x - temps[new_axis].max.x) > 0.001f) ||
		(fabs(temps[3].max.y - temps[new_axis].max.y) > 0.001f) ||
		(fabs(temps[3].max.z - temps[new_axis].max.z) > 0.001f)) {
		std::cout << (fabs(temps[3].min.x - temps[new_axis].min.x) > 0.001f) << " " <<
					(fabs(temps[3].min.y - temps[new_axis].min.y) > 0.001f) << " " <<
					(fabs(temps[3].min.z - temps[new_axis].min.z) > 0.001f) << " " <<
					(fabs(temps[3].max.x - temps[new_axis].max.x) > 0.001f) << " " <<
					(fabs(temps[3].max.y - temps[new_axis].max.y) > 0.001f) << " " <<
					(fabs(temps[3].max.z - temps[new_axis].max.z) > 0.001f) << std::endl;
		int u = 0;
	}
	if (elem_count == 34)
		int u = 0;*/
	return parentSAH;
}

void DSTree::calculateUpperTreeSAH(unsigned& new_axis, unsigned& new_index, unsigned* instances_ptr,
	simpleAABB& parentAABB, simpleAABB& leftAABB, simpleAABB& rightAABB, unsigned elem_count) {
	unsigned index = 0u;
	new_axis = 3u;
	float parentSAH = FLT_MAX;
	std::vector <unsigned> tempInstances;
	tempInstances.resize(elem_count);

	for (unsigned curr_axis = 0u; curr_axis < 3u; ++curr_axis) {
		std::vector <simpleAABB> right_surfs;

		current_scene_axis = curr_axis;
		qsort(instances_ptr, elem_count, sizeof(unsigned), cmpElemPosAABB);

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
		if (new_axis == curr_axis) {
			for (unsigned i = 0u; i < elem_count; ++i)
				tempInstances[i] = instances_ptr[i];
		}

		right_surfs.clear();
	}

	if (new_axis < 2u) {
		for (unsigned i = 0u; i < elem_count; ++i)
			instances_ptr[i] = tempInstances[i];
		//current_scene_axis = new_axis;
		//qsort(instances_ptr, elem_count, sizeof(unsigned), cmpElemPos);
	}

	return;
}


/* ========  Tree builder  ======== */

unsigned DSTree::buildHeader(bool is_leaf, bool is_carve, unsigned info) {
	unsigned header = unsigned(is_leaf) + (unsigned(is_carve) << 1);
	//unsigned field1 = (info >> 2) & 3u, field2 = info & 3u;
	if (info >= 16u)
		int u = 0;
	return (info << 2u) + header;
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
	unsigned size1 = dst_nodes.size();
	if (is_leaf && !triangle_index)
		int u = 0;
	unsigned sortInd[6u] = { 0u, 1u, 2u, 3u, 4u, 5u };
	unsigned split_plane = axis + 3u * unsigned(is_left_child);
	float child[6u] = { child_aabb.min.x, child_aabb.min.y, child_aabb.min.z, child_aabb.max.x, child_aabb.max.y, child_aabb.max.z };
	float dif[6u] = { child_aabb.min.x - parent_aabb.min.x,
					  child_aabb.min.y - parent_aabb.min.y,
					  child_aabb.min.z - parent_aabb.min.z,
					  parent_aabb.max.x - child_aabb.max.x,
					  parent_aabb.max.y - child_aabb.max.y,
					  parent_aabb.max.z - child_aabb.max.z };
	if (dif[0] < 0.f || dif[1] < 0.f || dif[2] < 0.f || dif[3] < 0.f || dif[4] < 0.f || dif[5] < 0.f)
		std::cout << "OBJECTION!" << std::endl;
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
				//if (is_leaf && i > 1u) {
				//	dst_nodes[dst_nodes.size() - 1u].leftChild = (dst_nodes[dst_nodes.size() - 1u].leftChild & 62u) + (triangle_index << 6) + 1u;
				//}
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
			DSNode temp_node{ ((dst_nodes.size() + 1) << 6) + buildHeader(false, true, info), 0u, {planes[0], planes[1]} };
			dst_nodes.push_back(temp_node);
		}
	}
	if (is_leaf && size1) {
		dst_nodes[dst_nodes.size() - 1u].leftChild = (dst_nodes[dst_nodes.size() - 1u].leftChild & 62u) + (triangle_index << 6) + 1u;
	}
	return size1;
	/*
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
			if (!sortTempPlanes(new_planes, minDif, maxDif, curr_difference, i << 1, curr_difference >= 0.01f))
				sortTempPlanes(new_possible_planes, minDif, maxDif, curr_difference, i << 1,
								curr_difference > LiteMath::EPSILON && curr_difference < 0.01f);
		}

		if (!(is_left_child && i == axis)) {
			float curr_difference = parent_aabb.max[i] - child_aabb.max[i];
			if (!sortTempPlanes(new_planes, minDif, maxDif, curr_difference, (i << 1) + 1u, curr_difference >= 0.01f))
				sortTempPlanes(new_possible_planes, minDif, maxDif, curr_difference, (i << 1) + 1u,
								curr_difference > LiteMath::EPSILON && curr_difference < 0.01f);
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
		unsigned left_child = dst_nodes.size() + 1u;
		if (is_leaf && last_carve) {
			left_child = triangle_index;
		}

		if (double_carve) {
			field1 = (unsigned((plane1 >> 1) == 1u) << 1) + unsigned((plane2 >> 1) == 2u);
			field2 = ((plane1 & 1u) << 1) + (plane2 & 1u);
		}

		new_carving_node.leftChild = (left_child << 6) +
			buildHeader(is_leaf && last_carve, true, (field1 << 2) + field2);

		dst_nodes.push_back(new_carving_node);
	}
	return dst_nodes.size() - size1;*/
}

void DSTree::dstBuilderRecur(unsigned first, unsigned last, float parentSAH) {
	unsigned curr_node_index = dst_nodes.size();
	unsigned new_index = last;

	//std::cout << "Node " << curr_node_index << ", " << first << "-" << last << std::endl;
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

			/*
			simpleAABB tempAABB, tempAABB2;
			LiteMath::float4 a = tempVertices[tempIndices[3 * instances[last]]];
			LiteMath::float4 b = tempVertices[tempIndices[3 * instances[last] + 1]];
			LiteMath::float4 c = tempVertices[tempIndices[3 * instances[last] + 2]];
			float minX = a.x, minY = a.y, minZ = a.z, maxX = a.x, maxY = a.y, maxZ = a.z;
			if (b.x < minX) minX = b.x;
			if (c.x < minX) minX = c.x;
			if (b.y < minY) minY = b.y;
			if (c.y < minY) minY = c.y;
			if (b.z < minZ) minZ = b.z;
			if (c.z < minZ) minZ = c.z;

			if (b.x > maxX) maxX = b.x;
			if (c.x > maxX) maxX = c.x;
			if (b.y > maxY) maxY = b.y;
			if (c.y > maxY) maxY = c.y;
			if (b.z > maxZ) maxZ = b.z;
			if (c.z > maxZ) maxZ = c.z;

			tempAABB.min.x = minX;
			tempAABB.min.y = minY;
			tempAABB.min.z = minZ;
			tempAABB.max.x = maxX;
			tempAABB.max.y = maxY;
			tempAABB.max.z = maxZ;

			tempAABB2 = getAABB(instances[last]);
			std::vector <unsigned> inds = { instances[last] };
			for (unsigned i = last - 1; i >= new_index + 1; --i) {
				inds.push_back(instances[i]);
				a = tempVertices[tempIndices[3 * instances[i]]];
				b = tempVertices[tempIndices[3 * instances[i] + 1]];
				c = tempVertices[tempIndices[3 * instances[i] + 2]];

				if (a.x < tempAABB.min.x) tempAABB.min.x = a.x;
				if (b.x < tempAABB.min.x) tempAABB.min.x = b.x;
				if (c.x < tempAABB.min.x) tempAABB.min.x = c.x;
				if (a.y < tempAABB.min.y) tempAABB.min.y = a.y;
				if (b.y < tempAABB.min.y) tempAABB.min.y = b.y;
				if (c.y < tempAABB.min.y) tempAABB.min.y = c.y;
				if (a.z < tempAABB.min.z) tempAABB.min.z = a.z;
				if (b.z < tempAABB.min.z) tempAABB.min.z = b.z;
				if (c.z < tempAABB.min.z) tempAABB.min.z = c.z;

				if (a.x > tempAABB.max.x) tempAABB.max.x = a.x;
				if (b.x > tempAABB.max.x) tempAABB.max.x = b.x;
				if (c.x > tempAABB.max.x) tempAABB.max.x = c.x;
				if (a.y > tempAABB.max.y) tempAABB.max.y = a.y;
				if (b.y > tempAABB.max.y) tempAABB.max.y = b.y;
				if (c.y > tempAABB.max.y) tempAABB.max.y = c.y;
				if (a.z > tempAABB.max.z) tempAABB.max.z = a.z;
				if (b.z > tempAABB.max.z) tempAABB.max.z = b.z;
				if (c.z > tempAABB.max.z) tempAABB.max.z = c.z;

				tempAABB2 = mergeAABB(instances[i], tempAABB2);
			}
			if (tempAABB.min.x != right_aabb.min.x ||
				tempAABB.min.y != right_aabb.min.y ||
				tempAABB.min.z != right_aabb.min.z ||
				tempAABB.max.x != right_aabb.max.x ||
				tempAABB.max.y != right_aabb.max.y ||
				tempAABB.max.z != right_aabb.max.z) {
				std::cout << (fabs(tempAABB.min.x - left_aabb.min.x) > 0.001f) << " " <<
					(fabs(tempAABB.min.y - right_aabb.min.y) > 0.001f) << " " <<
					(fabs(tempAABB.min.z - right_aabb.min.z) > 0.001f) << " " <<
					(fabs(tempAABB.max.x - right_aabb.max.x) > 0.001f) << " " <<
					(fabs(tempAABB.max.y - right_aabb.max.y) > 0.001f) << " " <<
					(fabs(tempAABB.max.z - right_aabb.max.z) > 0.001f) << std::endl;
				int u = 0;
			}
			*/

			dst_nodes[curr_node_index].planes[0] =  left_aabb.max[new_axis];
			dst_nodes[curr_node_index].planes[1] = right_aabb.min[new_axis];
			dst_nodes[curr_node_index].leftChild = (dst_nodes.size() << 6u) + buildHeader(false, false, new_axis);

			unsigned temp1 = buildCarvingNodes(curr_aabb, left_aabb, new_axis, true, (new_index - first) < DSTREE_MAX_POLY, first);
			if (!((dst_nodes[dst_nodes.size() - 1].leftChild & 1u) && temp1))
				dstBuilderRecur(first, new_index, parentSAH);

			for (unsigned temp_index = curr_node_index + 1u;;) {
				if (dst_nodes[temp_index].rightNode)
					std::cout << "WROMH ALGORITHMMMM" << std::endl;
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

	if (curr_node_index == 3)
		int u = 0;	//used for breakpoints in the code
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
	std::cout << "Triangles: " << indices.size() / 3 << std::endl;
	std::cout << "Instances: " << instances_info_size << std::endl;
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

	scene = this;
	dstBuilderRecurUpper(0u, instances_info_size - 1u);
	scene = nullptr;
	upper_tree = dst_nodes;
	dst_nodes.clear();
	if (instances_info_size)
		scene_AABB = instances_info[0].instAABB;
	for (uint32_t i = 1u; i < instances_info_size; ++i) {
		scene_AABB = mergeInstAABB(i, scene_AABB);
	}
	std::cout << "Upper tree: " << upper_tree.size() << "; Lower tree: " << lower_tree.size() << std::endl;
	//std::cout << "AABB min: " << scene_AABB.min.x << ", " << scene_AABB.min.y << ", " << scene_AABB.min.z << std::endl;
	//std::cout << "AABB max: " << scene_AABB.max.x << ", " << scene_AABB.max.y << ", " << scene_AABB.max.z << std::endl;
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
		float t1 = (Min[dim] - AABBeps - Position[dim]) * InvDir[dim];
		float t2 = (Max[dim] + AABBeps - Position[dim]) * InvDir[dim];

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
	tMinMax = LiteMath::float2(max(LiteMath::float2(tMin, 0.f)), tMax);
	return tMin < tMax && tMax >= 0.0f;
}

CRT_Hit DSTree::traceTriangle(LiteMath::float3 Position, LiteMath::float3 Direction, LiteMath::float4* tempVertices, unsigned* tempIndices) {
	CRT_Hit tempHitInfo{ FLT_MAX, -1, -1, -1, {0.f} };
	//tempHitInfo.t = FLT_MAX;
	//tempHitInfo.primId = -1;
	//tempHitInfo.geomId = -1;
	//tempHitInfo.instId = -1;

	//LiteMath::float3x3 vertMat = getTrVertices(instances[trIndex]);
	//const LiteMath::float3  a = vertMat.row[0u];
	//const LiteMath::float3 E1 = vertMat.row[1u] - a;
	//const LiteMath::float3 E2 = vertMat.row[2u] - a;

	const LiteMath::float3  a = LiteMath::to_float3(tempVertices[tempIndices[0]]);
	const LiteMath::float3 E1 = LiteMath::to_float3(tempVertices[tempIndices[1]]) - a;
	const LiteMath::float3 E2 = LiteMath::to_float3(tempVertices[tempIndices[2]]) - a;

	/*  ===============  */
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
			//tempHitInfo.instId = getInstID(instances[trIndex]);
			//tempHitInfo.geomId = instances_info[tempHitInfo.instId].geomID;
			//tempHitInfo.primId = instances[trIndex] - instances_info[tempHitInfo.instId].firstInstID;
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

	CRT_Hit temp_hit{ FLT_MAX, -1, -1, -1, {0.f} };

	for (uint32_t i = 0u; i < insts.size(); ++i) {
		//current_instance = insts[i];
		//if (dirAndFar.w > temp_hit.t) dirAndFar.w = temp_hit.t;
		CRT_Hit instHit = findHit(posAndNear, dirAndFar, false, insts[i]);
		//temp_hit.coords[2] += instHit.coords[2];
		//temp_hit.coords[3] += instHit.coords[3];
		if (instHit.t < temp_hit.t)
			temp_hit = instHit;
	}

	return temp_hit;
}

bool DSTree::RayQuery_AnyHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar) {
	std::vector <uint32_t> insts;
	findInstHit(posAndNear, dirAndFar, insts);

	CRT_Hit temp_hit{ FLT_MAX, -1, -1, -1, {0.f} };

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

void DSTree::findInstHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, std::vector <uint32_t> &insts) {
	LiteMath::float3 position(posAndNear.x, posAndNear.y, posAndNear.z);
	LiteMath::float3 direction(dirAndFar.x, dirAndFar.y, dirAndFar.z);
	LiteMath::float3 invDir(1.f / dirAndFar.x, 1.f / dirAndFar.y, 1.f / dirAndFar.z);
	LiteMath::float2 tMinMax;

	CRT_Hit temp_hit{ FLT_MAX, -1, -1, -1, {0.f} };

	if (!traceAABB(scene_AABB, position, direction, invDir, tMinMax)) return;

	unsigned node = 0u;
	size_t instances_size = instances_info.size();
	bool isFin[3] = { std::isfinite(invDir.x), std::isfinite(invDir.y), std::isfinite(invDir.z) };

	do {
		DSNode curr_node;
		bool CHECK_LEAF_POLY, IS_CARVING_NODE, IS_DOUBLE_CARVE;
		do {
			curr_node = upper_tree[node];
			node = curr_node.rightNode;

			unsigned header = curr_node.leftChild & 63u;
			curr_node.leftChild = curr_node.leftChild >> 6;
			CHECK_LEAF_POLY = header & 1u;
			IS_CARVING_NODE = header & 2u;
			IS_DOUBLE_CARVE = !((header & 48u) == 32u) && IS_CARVING_NODE;

			if (!CHECK_LEAF_POLY || IS_CARVING_NODE) {
				LiteMath::uint2 planes = LiteMath::uint2((header & 12u) >> 2);
				bool is_normal_positive2 = IS_CARVING_NODE;
				bool is_normal_positive1 = !is_normal_positive2;

				if (IS_DOUBLE_CARVE) {
					is_normal_positive1 = bool(planes.x & 2u);
					is_normal_positive2 = bool(planes.x & 1u);
					planes = LiteMath::uint2((header >> 5) & 1u, ((header >> 4) & 1u) + 1u);
				}

				bool planeTrav[2] = { is_normal_positive1 == (position[planes.x] < curr_node.planes[0u]),
									  is_normal_positive2 == (position[planes.y] < curr_node.planes[1u]) };
				planeTrav[0] = planeTrav[0] || (planeTrav[1] && IS_CARVING_NODE);
				if (!planeTrav[0] && (IS_CARVING_NODE || !CHECK_LEAF_POLY)) {
					float t[2] = { -1.f }, tMin[2] = { 0.f }, tMax[2] = { tMinMax.y };

					if (isFin[planes.x]) t[0] = invDir[planes.x] * (curr_node.planes[0u] - position[planes.x]);
					if (isFin[planes.y]) t[1] = invDir[planes.y] * (curr_node.planes[1u] - position[planes.y]);
					/*
					if (rightNodes.size() < treeDepth) {
						if (!IS_CARVING_NODE) temp_hit.coords[3] += weight;
						//if (IS_CARVING_NODE && !IS_DOUBLE_CARVE) temp_hit.coords[2] += weight;
						if (IS_CARVING_NODE) temp_hit.coords[2] += weight;
					}*/
					if (isFin[planes.x] && isFin[planes.y]) {
						unsigned farChild = t[0] < t[1];
						if (!IS_CARVING_NODE) {
							planeTrav[farChild] = planeTrav[farChild] || t[farChild] <= tMinMax.y;
							planeTrav[1 - farChild] = planeTrav[1 - farChild] || t[1 - farChild] >= 0.f;
						}

						if (IS_CARVING_NODE) {

							if (!IS_DOUBLE_CARVE) {
								planeTrav[0] = t[farChild] >= 0.f;
							}
							if (IS_DOUBLE_CARVE) {
								bool is_norm_pos[2] = { is_normal_positive1, is_normal_positive2 };
								if (is_norm_pos[1 - farChild] == (direction[planes[1 - farChild]] < 0.f) && (t[farChild] <= tMinMax.y ||
									is_norm_pos[farChild] == (position[planes[farChild]] < curr_node.planes[farChild])))
									planeTrav[0] = true;
							}
							/*
							tMin[0] = t[1 - farChild];
							tMin[0] = tMin[tMin[0] < tMin[1]];
							tMax[0] = t[farChild];
							tMax[0] = tMax[tMax[0] > tMax[1]];

							if (IS_DOUBLE_CARVE) {
								tMin[0] = t[0];
								tMax[0] = tMinMax.y;
								if (direction[planes.x] < 0.f == is_normal_positive1) {
									tMax[0] = tMin[0];
									tMin[0] = 0.f;
								}

								tMin[1] = t[1];
								if (direction[planes.y] < 0.f == is_normal_positive2) {
									tMax[1] = tMin[1];
									tMin[1] = 0.f;
								}

								tMin[0] = tMin[tMin[0] < tMin[1]];
								tMin[1] = 0.f;
								tMin[0] = tMin[tMin[0] < tMin[1]];
								tMax[0] = tMax[tMax[0] > tMax[1]];
								tMax[1] = tMinMax.y;
								tMax[0] = tMax[tMax[0] > tMax[1]];
							}
							planeTrav[0] = tMin[0] <= tMax[0];*/
						}
					}
					//planeTrav[0] = true;
				}
				planeTrav[1] = !planeTrav[0] && planeTrav[1] && !IS_CARVING_NODE;
				/*
				if (planeTrav[1] && rightNodes.size() < treeDepth) {
					//if (!IS_CARVING_NODE) temp_hit.coords[3] += weight;
					//if (IS_CARVING_NODE && !IS_DOUBLE_CARVE) temp_hit.coords[2] += weight;
					//if (IS_CARVING_NODE) temp_hit.coords[2] += weight;
				}*/
				if (planeTrav[1]) {
					node = upper_tree[curr_node.leftChild].rightNode;
					//rightNodes.push_back(rightNodes[rightNodes.size() - 1]);
				}
				if (planeTrav[0] && !CHECK_LEAF_POLY) {
					node = curr_node.leftChild;
					//if (!IS_CARVING_NODE) rightNodes.push_back(dst_ptr[node].rightNode);
				}/*
				if ((!planeTrav[0] && !planeTrav[1]) || CHECK_LEAF_POLY) {
					while (rightNodes.size() && rightNodes[rightNodes.size() - 1] == curr_node.rightNode)
						rightNodes.pop_back();
					if (rightNodes.size()) rightNodes.push_back(rightNodes[rightNodes.size() - 1]);
				}*/

				CHECK_LEAF_POLY = CHECK_LEAF_POLY && planeTrav[0];

				/*bool planeTrav1 = is_normal_positive1 == (position[planes.x] < curr_node.planes[0u]);
				bool planeTrav2 = is_normal_positive2 == (position[planes.y] < curr_node.planes[1u]);
				float tMin0 = tMinMax.y, tMin1 = tMinMax.y, tMax0 = 0.f, tMax1 = 0.f;
				if (!planeTrav1 && (direction[planes.x] < -LiteMath::EPSILON || direction[planes.x] > LiteMath::EPSILON)) {
					float temp_t = invDir[planes.x] * (curr_node.planes[0u] - position[planes.x]);
					planeTrav1 = temp_t >= 0.f && temp_t <= tMinMax.y;
					if (IS_DOUBLE_CARVE) {
						//planeTrav1 = planeTrav1 || (is_normal_positive2 == (position[planes.y] + temp_t * direction[planes.y] < curr_node.planes[1u]));
						tMin0 = temp_t;
						tMax0 = tMinMax.y;

						if (direction[planes.x] < 0.f == is_normal_positive1) {
							tMax0 = tMin0;
							tMin0 = 0.f;
						}
					}
				}
				if (!planeTrav2 && (direction[planes.y] < -LiteMath::EPSILON || direction[planes.y] > LiteMath::EPSILON)) {
					float temp_t = invDir[planes.y] * (curr_node.planes[1u] - position[planes.y]);
					planeTrav2 = temp_t >= 0.f && temp_t <= tMinMax.y;
					if (IS_DOUBLE_CARVE) {
						//planeTrav2 = planeTrav2 || (is_normal_positive1 == (position[planes.x] + temp_t * direction[planes.x] < curr_node.planes[0u]));
						tMin1 = temp_t;
						tMax1 = tMinMax.y;

						if (direction[planes.y] < 0.f == is_normal_positive2) {
							tMax1 = tMin1;
							tMin1 = 0.f;
						}
					}
				}
				if (IS_DOUBLE_CARVE) {
					tMin0 = max(0.f, max(tMin0, tMin1));
					tMax0 = min(tMinMax.y, min(tMax0, tMax1));
					planeTrav1 = planeTrav1 || tMin0 <= tMax0;
					//planeTrav2 = planeTrav2 || tMin0 <= tMax0;
				}
				if (IS_DOUBLE_CARVE) {
					//planeTrav1 = true;
					//planeTrav2 = true;
				}
				if (!CHECK_LEAF_POLY && (planeTrav1 || (planeTrav2 && IS_CARVING_NODE)))
					node = curr_node.leftChild;
				if (!planeTrav1 && planeTrav2 && !IS_CARVING_NODE)
					node = upper_tree[curr_node.leftChild].rightNode;

				CHECK_LEAF_POLY = CHECK_LEAF_POLY && (planeTrav1 || planeTrav2);*/
			}
		} while (!CHECK_LEAF_POLY && node != 0u);

		if (CHECK_LEAF_POLY) {
			insts.push_back(curr_node.leftChild);
		}
	} while (node != 0u);
	return;
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

//std::vector <unsigned> visitedNodes[3];
CRT_Hit DSTree::findHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, bool findAny, uint32_t current_instance) {
	//unsigned treeDepth = 10u;
	//float weight = 0.5f / (treeDepth + instances_info.size());
	//uint32_t tempVisIndex = 0u;
	//if (visitedNodes[0].size())
	//	tempVisIndex = 1;
	//visitedNodes[2] = TreePath(0, 11113);
	//std::vector <DSNode> routeNode;
	//for (unsigned i = 0u; i < lower_tree.size(); ++i) {
	//	routeNode.push_back(lower_tree[i]);
	//	routeNode[i].leftChild = routeNode[i].leftChild >> 6;
	//}

	simpleInstance tempInst = instances_info[current_instance];
	LiteMath::float3 position  = LiteMath::inverse4x4(tempInst.transf) * LiteMath::float3(posAndNear.x, posAndNear.y, posAndNear.z);
	LiteMath::float4 temp_direction = LiteMath::mul4x4x4(LiteMath::inverse4x4(tempInst.transf), LiteMath::float4( dirAndFar.x,  dirAndFar.y,  dirAndFar.z, 0.f));
	LiteMath::float3 direction = LiteMath::float3(temp_direction.x, temp_direction.y, temp_direction.z);
	LiteMath::float3 invDir    = LiteMath::float3(1.f / direction.x, 1.f / direction.y, 1.f / direction.z);
	LiteMath::float2 tMinMax(0.f, FLT_MAX);
	//std::vector <unsigned> rightNodes = { 0 };

	CRT_Hit temp_hit{ FLT_MAX, -1, -1, -1, {0.f} };
	if (!traceAABB(meshes[tempInst.geomID].meshAABB, position, direction, invDir, tMinMax)) return temp_hit;

	simpleMeshInfo tempMesh = meshes[tempInst.geomID];
<<<<<<< Updated upstream
=======
	unsigned tempMeshSize = lower_tree.size();
	if (tempInst.geomID + 1 < meshes.size())
		tempMeshSize = meshes[tempInst.geomID + 1].DSTreeOffset;
	unsigned treeDepthMin = 1u, treeDepthMax = 10u; //root has depth == 1
	//float weight = 0.01f * (treeDepthMax - treeDepthMin + 2u) * lower_tree.size() / ((treeDepthMax + 1) * (1.f + ((tempMeshSize - tempMesh.DSTreeOffset))));
	//float weight = 0.1f * log10f(lower_tree.size()) * (treeDepthMax - treeDepthMin + 2u) / ((1.f + (tempMeshSize - tempMesh.DSTreeOffset)) * logf(treeDepthMax + 2.f));
	float weight = 0.01f;
>>>>>>> Stashed changes
	LiteMath::float4* tempVertices = &(vertices[tempMesh.firstVertID]);
	unsigned* tempIndices = &(indices[3 * tempMesh.firstIndID]);
	unsigned* tempIndicesSorted = &(indices_sorted[tempMesh.firstIndID]);
	DSNode* dst_ptr = &(lower_tree[tempMesh.DSTreeOffset]);
	unsigned node = 0u;
	bool isFin[3] = { std::isfinite(invDir.x), std::isfinite(invDir.y), std::isfinite(invDir.z) };
	//unsigned current_node = 0u;
	do {
		DSNode curr_node;
		bool CHECK_LEAF_POLY, IS_CARVING_NODE, IS_DOUBLE_CARVE;
		do {
			curr_node = dst_ptr[node];
			//if (!current_instance)
			//	visitedNodes[tempVisIndex].push_back(node);
			//if (current_node + 1 < node)
			//	int u = 0;
			//current_node = node;
			node = curr_node.rightNode;

			unsigned header = curr_node.leftChild & 63u;
			curr_node.leftChild = curr_node.leftChild >> 6;
			CHECK_LEAF_POLY = header & 1u;
			IS_CARVING_NODE = header & 2u;
			IS_DOUBLE_CARVE = !((header & 48u) == 32u) && IS_CARVING_NODE;

			if (!CHECK_LEAF_POLY || IS_CARVING_NODE) {
				LiteMath::uint2 planes = LiteMath::uint2((header >> 2) & 3u);
				bool is_normal_positive2 = IS_CARVING_NODE;
				bool is_normal_positive1 = !is_normal_positive2;

				if (IS_DOUBLE_CARVE) {
					is_normal_positive1 = bool(planes.x & 2u);
					is_normal_positive2 = bool(planes.x & 1u);
					planes = LiteMath::uint2((header >> 5) & 1u, ((header >> 4) & 1u) + 1u);
				}
				//curr_node.planes[0u] += 0.2f;
				//if (is_normal_positive1) curr_node.planes[0u] -= 0.4f;
				//curr_node.planes[1u] += 0.2f;
				//if (is_normal_positive2) curr_node.planes[1u] -= 0.4f;

				bool planeTrav[2] = { is_normal_positive1 == (position[planes.x] < curr_node.planes[0u]),
									  is_normal_positive2 == (position[planes.y] < curr_node.planes[1u]) };
<<<<<<< Updated upstream
				planeTrav[0] = planeTrav[0] || (planeTrav[1] && IS_CARVING_NODE);
=======
				//bool planeTrav[2] = { false };
				//planeTrav[0] = planeTrav[0] || (planeTrav[1] && IS_CARVING_NODE);
>>>>>>> Stashed changes
				if (!planeTrav[0] && (IS_CARVING_NODE || !CHECK_LEAF_POLY)) {
					float t[2] = { -1.f }, tMin[2] = { 0.f }, tMax[2] = { tMinMax.y };

					if (isFin[planes.x]) t[0] = invDir[planes.x] * (curr_node.planes[0u] - position[planes.x]);
					if (isFin[planes.y]) t[1] = invDir[planes.y] * (curr_node.planes[1u] - position[planes.y]);
<<<<<<< Updated upstream
					/*
					if (rightNodes.size() < treeDepth) {
						if (!IS_CARVING_NODE) temp_hit.coords[3] += weight;
						if (IS_CARVING_NODE) temp_hit.coords[2] += weight;
					}*/
=======

					if (rightNodes.size() >= treeDepthMin && rightNodes.size() <= treeDepthMax) {
						if (t[0] >= 0.f && t[0] <= tMinMax.y)
							temp_hit.coords[3] += weight;
						if (t[1] >= 0.f && t[0] <= tMinMax.y)
							temp_hit.coords[2] += weight;
					}
>>>>>>> Stashed changes
					if (isFin[planes.x] && isFin[planes.y]) {
						unsigned farChild = t[0] < t[1];
						if (!IS_CARVING_NODE) {
							planeTrav[farChild] = planeTrav[farChild] || t[farChild] <= tMinMax.y;
							planeTrav[1 - farChild] = planeTrav[1 - farChild] || t[1 - farChild] >= 0.f;
						}

						if (IS_CARVING_NODE) {
							
							if (!IS_DOUBLE_CARVE) {
								planeTrav[0] = t[farChild] >= 0.f;
							}
							if (IS_DOUBLE_CARVE) {
								bool is_norm_pos[2] = { is_normal_positive1, is_normal_positive2 };
								if (is_norm_pos[1 - farChild] == (direction[planes[1 - farChild]] < 0.f) && (t[farChild] <= tMinMax.y ||
									is_norm_pos[farChild] == (position[planes[farChild]] < curr_node.planes[farChild])))
									planeTrav[0] = true;
							}
							/*
							tMin[0] = t[1 - farChild];
							tMin[0] = tMin[tMin[0] < tMin[1]];
							tMax[0] = t[farChild];
							tMax[0] = tMax[tMax[0] > tMax[1]];

							if (IS_DOUBLE_CARVE) {
								tMin[0] = t[0];
								tMax[0] = tMinMax.y;
								if (direction[planes.x] < 0.f == is_normal_positive1) {
									tMax[0] = tMin[0];
									tMin[0] = 0.f;
								}

								tMin[1] = t[1];
								if (direction[planes.y] < 0.f == is_normal_positive2) {
									tMax[1] = tMin[1];
									tMin[1] = 0.f;
								}

								tMin[0] = tMin[tMin[0] < tMin[1]];
								tMin[1] = 0.f;
								tMin[0] = tMin[tMin[0] < tMin[1]];
								tMax[0] = tMax[tMax[0] > tMax[1]];
								tMax[1] = tMinMax.y;
								tMax[0] = tMax[tMax[0] > tMax[1]];
							}
							planeTrav[0] = tMin[0] <= tMax[0];*/
						}
					}
				}
				planeTrav[1] = !planeTrav[0] && planeTrav[1] && !IS_CARVING_NODE;
				//planeTrav[0] = true;
				/*
<<<<<<< Updated upstream
				if (planeTrav[1] && rightNodes.size() < treeDepth) {
					//if (!IS_CARVING_NODE) temp_hit.coords[3] += weight;
					//if (IS_CARVING_NODE) temp_hit.coords[2] += weight;
=======
				if ((planeTrav[0] || planeTrav[1]) && (!CHECK_LEAF_POLY || IS_CARVING_NODE) && rightNodes.size() >= treeDepthMin && rightNodes.size() <= treeDepthMax) {
					if (!IS_CARVING_NODE) temp_hit.coords[3] += weight;
					if (IS_CARVING_NODE) temp_hit.coords[2] += weight;
>>>>>>> Stashed changes
				}*/
				if (planeTrav[1]) {
					node = dst_ptr[curr_node.leftChild].rightNode;
					//rightNodes.push_back(rightNodes[rightNodes.size() - 1]);
				}
				if (planeTrav[0] && !CHECK_LEAF_POLY) {
					node = curr_node.leftChild;
					//if (!IS_CARVING_NODE) rightNodes.push_back(dst_ptr[node].rightNode);
				}/*
				if ((!planeTrav[0] && !planeTrav[1]) || CHECK_LEAF_POLY) {
					while (rightNodes.size() && rightNodes[rightNodes.size() - 1] == curr_node.rightNode)
						rightNodes.pop_back();
					if (rightNodes.size()) rightNodes.push_back(rightNodes[rightNodes.size() - 1]);
				}*/

				CHECK_LEAF_POLY = CHECK_LEAF_POLY && planeTrav[0];
				/*
				bool planeTrav1 = false;//is_normal_positive1 == (position[planes.x] < curr_node.planes[0u]);
				bool planeTrav2 = false;//is_normal_positive2 == (position[planes.y] < curr_node.planes[1u]);
				float tMin0 = tMinMax.y, tMin1 = tMinMax.y, tMax0 = 0.f, tMax1 = 0.f;
				if (!planeTrav1 && (direction[planes.x] < -LiteMath::EPSILON || direction[planes.x] > LiteMath::EPSILON)) {
					float temp_t = invDir[planes.x] * (curr_node.planes[0u] - position[planes.x]);
					planeTrav1 = temp_t >= 0.f && temp_t <= tMinMax.y;
					if (IS_DOUBLE_CARVE) {
						//planeTrav1 = planeTrav1 || (is_normal_positive2 == (position[planes.y] + temp_t * direction[planes.y] < curr_node.planes[1u]));
						tMin0 = temp_t;
						tMax0 = tMinMax.y;

						if (direction[planes.x] < 0.f == is_normal_positive1) {
							tMax0 = tMin0;
							tMin0 = 0.f;
						}
					}
				}
				if (!planeTrav2 && (direction[planes.y] < -LiteMath::EPSILON || direction[planes.y] > LiteMath::EPSILON)) {
					float temp_t = invDir[planes.y] * (curr_node.planes[1u] - position[planes.y]);
					planeTrav2 = temp_t >= 0.f && temp_t <= tMinMax.y;
					if (IS_DOUBLE_CARVE) {
						//planeTrav2 = planeTrav2 || (is_normal_positive1 == (position[planes.x] + temp_t * direction[planes.x] < curr_node.planes[0u]));
						tMin1 = temp_t;
						tMax1 = tMinMax.y;

						if (direction[planes.y] < 0.f == is_normal_positive2) {
							tMax1 = tMin1;
							tMin1 = 0.f;
						}
					}
				}
				if (IS_DOUBLE_CARVE) {
					tMin0 = max(0.f, max(tMin0, tMin1));
					tMax0 = min(tMinMax.y, min(tMax0, tMax1));
					planeTrav1 = planeTrav1 || tMin0 <= tMax0;
					//planeTrav2 = planeTrav2 || tMin0 <= tMax0;
				}
				if (planeTrav1 && rightNodes.size() < treeDepth) {
					if (!IS_CARVING_NODE) temp_hit.coords[1] += weight;
					//if (IS_CARVING_NODE && !IS_DOUBLE_CARVE) temp_hit.coords[2] += weight;
					if (IS_CARVING_NODE) temp_hit.coords[3] += weight;
				}
				if (planeTrav2 && rightNodes.size() < treeDepth) {
					if (!IS_CARVING_NODE) temp_hit.coords[1] += weight;
					//if (IS_CARVING_NODE && !IS_DOUBLE_CARVE) temp_hit.coords[2] += weight;
					if (IS_CARVING_NODE) temp_hit.coords[3] += weight;
				}
				
				if (IS_CARVING_NODE && !IS_DOUBLE_CARVE) {
					//planeTrav1 = true;
					//planeTrav2 = true;
				}
				if (!CHECK_LEAF_POLY && (planeTrav1 || (planeTrav2 && IS_CARVING_NODE))) {
					node = curr_node.leftChild;
					if (!IS_CARVING_NODE) rightNodes.push_back(dst_ptr[node].rightNode);
				}
				else if (!planeTrav1 && planeTrav2 && !IS_CARVING_NODE) {
					node = dst_ptr[curr_node.leftChild].rightNode;
					rightNodes.push_back(dst_ptr[node].rightNode);
				}
				else {
					while (rightNodes.size() && rightNodes[rightNodes.size() - 1] == curr_node.rightNode)
						rightNodes.pop_back();
					rightNodes.push_back(dst_ptr[node].rightNode);
				}

				CHECK_LEAF_POLY = CHECK_LEAF_POLY && (planeTrav1 || planeTrav2);*/
			}/*
			else {
				while (rightNodes.size() && rightNodes[rightNodes.size() - 1] == curr_node.rightNode)
					rightNodes.pop_back();
					if (rightNodes.size()) rightNodes.push_back(rightNodes[rightNodes.size() - 1]);
				}*/
		} while (!CHECK_LEAF_POLY && node != 0u);

		if (CHECK_LEAF_POLY) {
			/*unsigned tr_num = DSTREE_MAX_POLY + curr_node.leftChild;
			if (IS_CARVING_NODE) {
				unsigned i = current_node;
				while (dst_ptr[i].rightNode || !(dst_ptr[i].leftChild & 1u)) {
					++i;
					if (dst_ptr[i].leftChild & 1u) {
						tr_num = dst_ptr[i].leftChild >> 6;
						break;
					}
				}
			}
			else {
				tr_num = static_cast<unsigned>(curr_node.planes[0]) + curr_node.leftChild;
			}*/
			unsigned tr_num = DSTREE_MAX_POLY;
			if (!IS_CARVING_NODE)
				tr_num = static_cast<unsigned>(curr_node.planes[0]);
			tr_num += curr_node.leftChild;
			if (tr_num > tempMesh.lastIndID - tempMesh.firstIndID + 1) tr_num = tempMesh.lastIndID - tempMesh.firstIndID + 1;
			for (unsigned j = curr_node.leftChild; j < tr_num; ++j) {
				CRT_Hit temp_t = traceTriangle(position, direction, tempVertices, &(tempIndices[3 * tempIndicesSorted[j]]));
				temp_t.instId = current_instance;
				temp_t.geomId = tempInst.geomID;
				temp_t.primId = tempIndicesSorted[j];
				//temp_t.coords[2] = temp_hit.coords[2];
				//temp_t.coords[3] = temp_hit.coords[3];
				if (temp_t.t < temp_hit.t && temp_t.t <= tMinMax.y && temp_t.t >= tMinMax.x) {
				//if (temp_t.t < temp_hit.t && temp_t.t <= dirAndFar.w && temp_t.t >= posAndNear.w) {
					/*if (!(temp_t.t <= tMinMax.y && temp_t.t >= tMinMax.x)) {
						std::cout << "Bad triangle " << j << ", inst " << current_instance << "; AABB min t = " << tMinMax.x << "; AABB max t = " << tMinMax.y << "; trig = " << temp_t.t << std::endl;
						std::cout << "min index: " << tempMesh.firstIndID << ", max index: " << tempMesh.lastIndID << ", this index: " << tempIndicesSorted[j] << std::endl;
						std::cout << "AABB min: " << tempMesh.meshAABB.min.x << ", " << tempMesh.meshAABB.min.y << ", " << tempMesh.meshAABB.min.z << std::endl;
						std::cout << "AABB max: " << tempMesh.meshAABB.max.x << ", " << tempMesh.meshAABB.max.y << ", " << tempMesh.meshAABB.max.z << "\n" << std::endl;

						const LiteMath::float3 a = LiteMath::to_float3(tempVertices[tempIndices[3 * tempIndicesSorted[j]    ]]);
						const LiteMath::float3 b = LiteMath::to_float3(tempVertices[tempIndices[3 * tempIndicesSorted[j] + 1]]);
						const LiteMath::float3 c = LiteMath::to_float3(tempVertices[tempIndices[3 * tempIndicesSorted[j] + 2]]);
						std::cout << "Vert a: " << a.x << ", " << a.y << ", " << a.z << std::endl;
						std::cout << "Vert b: " << b.x << ", " << b.y << ", " << b.z << std::endl;
						std::cout << "Vert c: " << c.x << ", " << c.y << ", " << c.z << "\n" << std::endl;
					}*/
					///temp_hit.instId = current_instance;
					//temp_hit.geomId = instances_info[current_instance].geomID;
					//temp_hit.primId = tempIndicesSorted[j];
					//temp_hit.coords[0u] = temp_t.coords[0u];
					//temp_hit.coords[1u] = temp_t.coords[1u];
					temp_hit = temp_t;
					if (findAny) return temp_hit;
				}
			}
		}
	} while (node != 0u);
	/*if (rightNodes.size())
		std::cout << "SIZE " << rightNodes.size() << std::endl;
	if (temp_hit.geomId == 0)
		int u = 0;
	if (tempVisIndex && !current_instance) {
		for (unsigned i = 0u; i < visitedNodes[0].size(); ++i) {
			if (visitedNodes[0][i] != visitedNodes[1][i])
				std::cout << "wrong " << i << " : " << visitedNodes[0][i] << ", " << visitedNodes[1][i] << std::endl;
		}
	}*/
	return temp_hit;
}


ISceneObject* CreateSceneRT(const char* a_impleName) { return new DSTree; }

void DeleteSceneRT(ISceneObject* a_pScene) { delete a_pScene; }