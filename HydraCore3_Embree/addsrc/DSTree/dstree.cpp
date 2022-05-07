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
	}
	if (a_vertNumber) {
		tempInfo.firstVertID = vertices.size();

		for (size_t j = 0; j < a_vertNumber; ++j)
			vertices.push_back(a_vpos4f[j]);

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
	//uint32_t instID = instances_info.size();
	//tempInst.firstInstID = instances.size();
/*
	uint32_t firstIndex = 0u;
	
	for (uint32_t i = meshes[a_geomId].firstIndID; meshes[a_geomId].firstIndID <= meshes[a_geomId].lastIndID && i <= meshes[a_geomId].lastIndID; ++i) {
		instances.push_back(instances.size());
	}
	tempInst.lastInstID = instances.size() - 1;*/
	instances_info.push_back(tempInst);
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

			  /*	============ Don't need them, as vertices i/o is simplified ============	*/
/* 
uint32_t DSTree::getTrIndex(uint32_t instTrIndex) {
	uint32_t instID = 0;
	while (instTrIndex > instances_info[instID].lastInstID) ++instID;

	return instTrIndex - instances_info[instID].firstInstID + meshes[instances_info[instID].geomID].firstIndID;
}
//not needed, program always works with a certain instance "current_instance"
uint32_t DSTree::getInstID(uint32_t instTrIndex) {
	uint32_t instID = 0;
	while (instTrIndex > instances_info[instID].lastInstID) ++instID;

	return instID;
}

LiteMath::float3& DSTree::getVertex(unsigned instance_ind, unsigned* instances_ptr, unsigned geom_vert_ind) {
	uint32_t instTrIndex = instances_ptr[instance_ind];
	uint32_t instID = 0;
	while (instTrIndex > instances_info[instID].lastInstID) ++instID;

	simpleInstance temp_inst_info = instances_info[instID];
	uint32_t vecIndex = 3 * (instTrIndex - temp_inst_info.firstInstID + meshes[temp_inst_info.geomID].firstIndID) + geom_vert_ind;
	LiteMath::float4 coords = temp_inst_info.transf * vertices[indices[vecIndex]];

	return LiteMath::float3(coords.x, coords.y, coords.z);
}

LiteMath::float3x3& DSTree::getTrVertices(unsigned instTrIndex) {
	uint32_t instID = 0;
	while (instTrIndex > instances_info[instID].lastInstID) ++instID;

	simpleInstance temp_inst_info = instances_info[instID];
	uint32_t vecIndex = 3 * (instTrIndex - temp_inst_info.firstInstID + meshes[temp_inst_info.geomID].firstIndID);
	LiteMath::uint3 tempIndices = LiteMath::uint3(indices[vecIndex], indices[vecIndex + 1u], indices[vecIndex + 2u]);

	LiteMath::float4 coords1 = temp_inst_info.transf * vertices[tempIndices.x];
	LiteMath::float4 coords2 = temp_inst_info.transf * vertices[tempIndices.y];
	LiteMath::float4 coords3 = temp_inst_info.transf * vertices[tempIndices.z];

	LiteMath::float3x3 tempMat;
	tempMat.row[0] = LiteMath::float3(coords1.x, coords1.y, coords1.z);
	tempMat.row[1] = LiteMath::float3(coords2.x, coords2.y, coords2.z);
	tempMat.row[2] = LiteMath::float3(coords3.x, coords3.y, coords3.z);

	return tempMat;
}

LiteMath::float3& DSTree::getTrVerticesAxis(unsigned instTrIndex) {
	uint32_t instID = 0;
	while (instTrIndex > instances_info[instID].lastInstID) ++instID;

	simpleInstance temp_inst_info = instances_info[instID];
	uint32_t vecIndex = 3 * (instTrIndex - temp_inst_info.firstInstID + meshes[temp_inst_info.geomID].firstIndID);
	LiteMath::uint3 tempIndices = LiteMath::uint3(indices[vecIndex], indices[vecIndex + 1u], indices[vecIndex + 2u]);

	return LiteMath::float3((temp_inst_info.transf * vertices[tempIndices.x])[current_scene_axis],
							(temp_inst_info.transf * vertices[tempIndices.y])[current_scene_axis], 
							(temp_inst_info.transf * vertices[tempIndices.z])[current_scene_axis]);
}*/


float DSTree::getMinPos(const void* index) {
	unsigned temp_index = *((unsigned*) index) * 3u;

	//LiteMath::float3 coords = scene->getTrVerticesAxis(*temp_index);
	LiteMath::float3 coords = LiteMath::float3(tempVertices[tempIndices[temp_index     ]][current_scene_axis],
											   tempVertices[tempIndices[temp_index + 1u]][current_scene_axis],
											   tempVertices[tempIndices[temp_index + 2u]][current_scene_axis]);
	return min(LiteMath::float2(coords.x, min(LiteMath::float2(coords.y, coords.z))));
}

float DSTree::getMaxPos(const void* index) {
	//unsigned* temp_index = (unsigned*) index;
	//LiteMath::float3 coords = scene->getTrVerticesAxis(*temp_index);
	unsigned temp_index = *((unsigned*) index) * 3u;
	DSTree* temp = this;
	LiteMath::float3 coords = LiteMath::float3(tempVertices[tempIndices[temp_index     ]][current_scene_axis],
											   tempVertices[tempIndices[temp_index + 1u]][current_scene_axis],
											   tempVertices[tempIndices[temp_index + 2u]][current_scene_axis]);

	return max(LiteMath::float2(coords.x, max(LiteMath::float2(coords.y, coords.z))));
}

float DSTree::getMinAABBpos(const void* index) {
	//unsigned* temp_index = (unsigned*)index;
	return instances_info[*((unsigned*)index)].instAABB.min[current_scene_axis];
}

float DSTree::getMaxAABBpos(const void* index) {
	//unsigned* temp_index = (unsigned*)index;
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
	unsigned index = 0u;
	new_axis = 3u;

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
		for (unsigned k = 1u; k < elem_count - 1u; ++k) {
			temp_aabb = mergeAABB(instances_ptr[k - 1u], temp_aabb);

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
		qsort(instances_ptr, elem_count, sizeof(unsigned), cmpElemPos);
	}

	return parentSAH;
}

void DSTree::calculateUpperTreeSAH(unsigned& new_axis, unsigned& new_index, unsigned* instances_ptr,
	simpleAABB& parentAABB, simpleAABB& leftAABB, simpleAABB& rightAABB, unsigned elem_count) {
	unsigned index = 0u;
	new_axis = 3u;
	float parentSAH = FLT_MAX;

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

		right_surfs.clear();
	}

	if (new_axis < 2u) {
		current_scene_axis = new_axis;
		qsort(instances_ptr, elem_count, sizeof(unsigned), cmpElemPosAABB);
	}

	return;
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
	unsigned size1 = dst_nodes.size();
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
	return dst_nodes.size() - size1;
}

void DSTree::dstBuilderRecur(unsigned first, unsigned last, float parentSAH) {
	unsigned curr_node_index = dst_nodes.size();
	unsigned new_index = last;

	//std::cout << "Node " << curr_node_index << ", " << first << "-" << last << std::endl;
	dst_nodes.emplace_back();
	dst_nodes[curr_node_index].leftChild = (first << 6u) + 1u;
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

			buildCarvingNodes(curr_aabb, left_aabb, new_axis, true, (new_index - first) < DSTREE_MAX_POLY, first);
			if (!((new_index - first) < DSTREE_MAX_POLY && (dst_nodes.size() - 1u != curr_node_index)))
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

			unsigned carve_nodes_num = buildCarvingNodes(curr_aabb, right_aabb, new_axis, false, (last - new_index - 1u) < DSTREE_MAX_POLY, new_index + 1u);

			if (!((last - new_index - 1u) < DSTREE_MAX_POLY && carve_nodes_num))
				dstBuilderRecur(new_index + 1u, last, parentSAH);
		}
	}
	return;
}

void DSTree::dstBuilderRecurUpper(unsigned first, unsigned last) {
	unsigned curr_node_index = dst_nodes.size();
	unsigned new_index = last;

	if (curr_node_index == 3)	///////////
		int u = 0;
	dst_nodes.emplace_back();
	dst_nodes[curr_node_index].leftChild = (first << 6u) + 1u;
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

			buildCarvingNodes(curr_aabb, left_aabb, new_axis, true, new_index == first, first);
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

			unsigned carve_nodes_num = buildCarvingNodes(curr_aabb, right_aabb, new_axis, false, !bool(last - new_index - 1u), new_index + 1u);

			if (!(last == new_index + 1u && carve_nodes_num))
				dstBuilderRecurUpper(new_index + 1u, last);
		}
		/*else if (last - first == 1) {
			dstBuilderRecurUpper(first, first, parentSAH);
			dstBuilderRecurUpper(last, last, parentSAH);
		}*/
	}
	return;
}

void DSTree::CommitScene() {
	if (indices.empty())
		return;

	//size_t instances_size = instances.size();
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
	//getMeshAABB, mergeMeshAABB !!!
	if (instances_info_size)
		scene_AABB = instances_info[0].instAABB;
	for (uint32_t i = 1u; i < instances_info_size; ++i) {
		scene_AABB = mergeInstAABB(i, scene_AABB);
	}
	std::cout << "Upper: " << upper_tree.size() << "; Lower: " << lower_tree.size() << std::endl;
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
	dstFile.close();*/
	std::cout << "Nodes: " << dst_nodes.size() << std::endl;
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
	tMinMax = LiteMath::float2(max(LiteMath::float2(tMin, 0.f)), tMax);
	return tMin < tMax && tMax >= 0.0f;
}

CRT_Hit DSTree::traceTriangle(LiteMath::float3 Position, LiteMath::float3 Direction, LiteMath::float4* tempVertices, unsigned* tempIndices) {
	CRT_Hit tempHitInfo;
	tempHitInfo.t = FLT_MAX;
	tempHitInfo.primId = -1;
	tempHitInfo.geomId = -1;
	tempHitInfo.instId = -1;

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

		if (temp_u >= 0.0f && temp_u <= 1.0f && ((temp_u + temp_v) <= 1.0f) &&
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
	if (fabs(dirAndFar.x + 0.346856236) < 1e-4)
		int u = 0;
	findInstHit(posAndNear, dirAndFar, insts);
	//for (uint32_t i = 0u; i < instances_info.size(); ++i)
	//	insts.push_back(i);

	CRT_Hit temp_hit;
	temp_hit.t = FLT_MAX;
	temp_hit.primId = -1;
	temp_hit.geomId = -1;
	temp_hit.instId = -1;

	for (uint32_t i = 0u; i < insts.size(); ++i) {
		current_instance = insts[i];

		if (!current_instance)
			int u = 0;

		CRT_Hit instHit = findHit(posAndNear, dirAndFar, false, insts[i]);
		if (instHit.t < temp_hit.t)
			temp_hit = instHit;
	}

	return temp_hit;
}

bool DSTree::RayQuery_AnyHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar) {
	std::vector <uint32_t> insts;
	findInstHit(posAndNear, dirAndFar, insts);

	CRT_Hit temp_hit;
	temp_hit.t = FLT_MAX;
	temp_hit.primId = -1;
	temp_hit.geomId = -1;
	temp_hit.instId = -1;

	for (uint32_t i = 0u; i < insts.size(); ++i) {
		current_instance = insts[i];

		simpleMeshInfo tempMesh = meshes[instances_info[current_instance].geomID];
		tempVertices = &(vertices[tempMesh.firstVertID]);
		tempIndices = &(indices[3 * tempMesh.firstIndID]);
		tempIndicesSorted = &(indices_sorted[tempMesh.firstIndID]);
		dst_ptr = &(lower_tree[tempMesh.DSTreeOffset]);

		// NOT DONE YET
		CRT_Hit instHit = findHit(posAndNear, dirAndFar, true, insts[i]);
		if (instHit.t < temp_hit.t)
			temp_hit = instHit;
	}

	return (temp_hit.t >= posAndNear.w) && (temp_hit.t <= dirAndFar.w);
}

void DSTree::findInstHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, std::vector <uint32_t> &insts) {
	LiteMath::float3 position = LiteMath::float3(posAndNear.x, posAndNear.y, posAndNear.z);
	LiteMath::float3 direction = LiteMath::float3(dirAndFar.x, dirAndFar.y, dirAndFar.z);
	LiteMath::float3 invDir = LiteMath::float3(1.f / dirAndFar.x, 1.f / dirAndFar.y, 1.f / dirAndFar.z);
	LiteMath::float2 tMinMax;

	CRT_Hit temp_hit;
	temp_hit.t = FLT_MAX;
	temp_hit.primId = -1;
	temp_hit.geomId = -1;
	temp_hit.instId = -1;

	if (!traceAABB(scene_AABB, position, direction, invDir, tMinMax)) return;

	unsigned node = 0u;
	//unsigned nodesSize = upper_tree.size();
	size_t instances_size = instances_info.size();
	bool isFin[3] = { std::isfinite(invDir.x), std::isfinite(invDir.y), std::isfinite(invDir.z) };

	do {
		DSNode curr_node;
		bool CHECK_LEAF_POLY, IS_CARVING_NODE, IS_DOUBLE_CARVE;
		do {
			///TREE IS BUILT WRONG
			curr_node = upper_tree[node];
			node = curr_node.rightNode;
			//if (node >= nodesSize) node = 0u; //for the last leaf in the DSTree

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

				bool planeTrav1 = is_normal_positive1 == (position[planes.x] < curr_node.planes[0u]);
				bool planeTrav2 = is_normal_positive2 == (position[planes.y] < curr_node.planes[1u]);

				if (!planeTrav1 && (direction[planes.x] < -LiteMath::EPSILON || direction[planes.x] > LiteMath::EPSILON)) {
					float temp_t = invDir[planes.x] * (curr_node.planes[0u] - position[planes.x]);
					planeTrav1 = temp_t >= 0.f && temp_t < temp_hit.t;
				}
				if (!planeTrav2 && (direction[planes.y] < -LiteMath::EPSILON || direction[planes.y] > LiteMath::EPSILON)) {
					float temp_t = invDir[planes.y] * (curr_node.planes[1u] - position[planes.y]);
					planeTrav2 = temp_t >= 0.f && temp_t < temp_hit.t;
				}
				//planeTrav1 = true;
				//planeTrav2 = true;
				if (planeTrav1 && (!IS_CARVING_NODE || (planeTrav2 && !CHECK_LEAF_POLY)))
					node = curr_node.leftChild;
				if (!planeTrav1 && planeTrav2 && !IS_CARVING_NODE)
					node = upper_tree[curr_node.leftChild].rightNode;

				CHECK_LEAF_POLY = CHECK_LEAF_POLY && planeTrav1 && planeTrav2;
			}
		} while (!CHECK_LEAF_POLY && node != 0u);

		if (CHECK_LEAF_POLY) {
			insts.push_back(curr_node.leftChild);
		}
	} while (node != 0u);
	return;
}

CRT_Hit DSTree::findHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, bool findAny, uint32_t current_instance) {
	simpleMeshInfo tempMesh = meshes[instances_info[current_instance].geomID];
	//std::cout << "IND " << tempMesh.lastIndID << std::endl;
	LiteMath::float4* tempVertices = &(vertices[tempMesh.firstVertID]);
	unsigned* tempIndices = &(indices[3 * tempMesh.firstIndID]);
	unsigned* tempIndicesSorted = &(indices_sorted[tempMesh.firstIndID]);
	DSNode* dst_ptr = &(lower_tree[tempMesh.DSTreeOffset]);

	//change direction vector
	LiteMath::float3 position  = LiteMath::inverse4x4(instances_info[current_instance].transf) * LiteMath::float3(posAndNear.x, posAndNear.y, posAndNear.z);
	LiteMath::float4 temp_direction = LiteMath::mul4x4x4(LiteMath::inverse4x4(instances_info[current_instance].transf), LiteMath::float4( dirAndFar.x,  dirAndFar.y,  dirAndFar.z, 0.f));
	LiteMath::float3 direction = LiteMath::float3(temp_direction.x, temp_direction.y, temp_direction.z);
	LiteMath::float3 invDir    = LiteMath::float3(1.f / direction.x, 1.f / direction.y, 1.f / direction.z);
	LiteMath::float2 tMinMax;
	if (!current_instance)
		int u = 0;
	CRT_Hit temp_hit;
	temp_hit.t = FLT_MAX;
	temp_hit.primId = -1;
	temp_hit.geomId = -1;
	temp_hit.instId = -1;
	/*
	for (unsigned j = 0; j < 10; ++j) {
		CRT_Hit temp_t = traceTriangle(position, direction, j);

		if (temp_t.t < temp_hit.t && temp_t.t <= dirAndFar.w && temp_t.t >= posAndNear.w) {
			temp_hit = temp_t;
			if (findAny) return temp_hit;
		}
	}*/

	//dir if good: 0.346856147, -0.346856236, 0.871425033  //after matrix transform, what matters is they're different
	//dir if bad: direction {x=0.142260358 y=-0.323869467 z=0.935345173}
	//if (!traceAABB(meshes[instances_info[current_instance].geomID].meshAABB, position, direction, invDir, tMinMax)) return temp_hit;

	unsigned node = 0u;
	//unsigned nodesSize = dst_nodes.size();
	//size_t instances_size = instances.size();
	bool isFin[3] = { std::isfinite(invDir.x), std::isfinite(invDir.y), std::isfinite(invDir.z) };
	unsigned current_node = 0u;
	do {
		DSNode curr_node;
		bool CHECK_LEAF_POLY, IS_CARVING_NODE, IS_DOUBLE_CARVE;
		do {
			curr_node = dst_ptr[node];
			current_node = node;
			node = curr_node.rightNode;
			//++node;
			//if (node >= nodesSize) node = 0u; //for the last leaf in the DSTree

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
					//is_normal_positive = bool2(bool(planes.x & 2u), bool(planes.x & 1u));
					is_normal_positive1 = bool(planes.x & 2u);
					is_normal_positive2 = bool(planes.x & 1u);
					planes = LiteMath::uint2((header >> 5) & 1u, ((header >> 4) & 1u) + 1u);
				}

				bool planeTrav1 = is_normal_positive1 == (position[planes.x] < curr_node.planes[0u]);
				bool planeTrav2 = is_normal_positive2 == (position[planes.y] < curr_node.planes[1u]);
/*				bool nodeTrav = planeTrav1 && (!IS_CARVING_NODE || planeTrav2);
				bool splitTraceRight = planeTrav2 && !(IS_CARVING_NODE || CHECK_LEAF_POLY || planeTrav1);

				if (!nodeTrav && (isFin[planes.x] || isFin[planes.y])) {
					float temp_t[2] = { 0.f, 0.f };
					if (isFin[planes.x]) temp_t[0] = invDir[planes.x] * (curr_node.planes[0u] - position[planes.x]);
					if (isFin[planes.y]) temp_t[1] = invDir[planes.y] * (curr_node.planes[1u] - position[planes.y]);

					if (IS_CARVING_NODE && !IS_DOUBLE_CARVE) {
						nodeTrav = !(temp_t[temp_t[0] > temp_t[1]] > tMinMax[1] || temp_t[temp_t[0] < temp_t[1]] < tMinMax[0]);
					}

					if (!IS_CARVING_NODE) {
						nodeTrav = temp_t[0] < tMinMax[1];
						splitTraceRight = temp_t[1] < tMinMax[1];
					}

					if (IS_DOUBLE_CARVE) {
						nodeTrav = (!isFin[planes.x] && temp_t[1] < tMinMax[1]) || (!isFin[planes.y] && temp_t[0] < tMinMax[1]);

						if (!nodeTrav && isFin[planes.x] && isFin[planes.y])
							nodeTrav = is_normal_positive2 == ((position[planes.y] + temp_t[0] * direction[planes.y]) < curr_node.planes[1u]);
					}
				}

				if (!nodeTrav) node = curr_node.rightNode;
				if (splitTraceRight) node = dst_nodes.nodes[curr_node.leftChild].rightNode;

				CHECK_LEAF_POLY = CHECK_LEAF_POLY && nodeTrav;
				*/
				if (!planeTrav1 && (direction[planes.x] < -LiteMath::EPSILON || direction[planes.x] > LiteMath::EPSILON)) {
					float temp_t = invDir[planes.x] * (curr_node.planes[0u] - position[planes.x]);
					planeTrav1 = temp_t >= 0.f && temp_t < temp_hit.t;

					//if (planeTrav1 && IS_DOUBLE_CARVE) planeTrav1 = is_normal_positive2 == (position[planes.y] + temp_t * direction[planes.y] < curr_node.planes[1u]);
				}
				if (!planeTrav2 && (direction[planes.y] < -LiteMath::EPSILON || direction[planes.y] > LiteMath::EPSILON)) {
					float temp_t = invDir[planes.y] * (curr_node.planes[1u] - position[planes.y]);
					planeTrav2 = temp_t >= 0.f && temp_t < temp_hit.t;
				}
				if (planeTrav1 && (!IS_CARVING_NODE || (planeTrav2 && !CHECK_LEAF_POLY)))
					node = curr_node.leftChild;
				if (!planeTrav1 && planeTrav2 && !IS_CARVING_NODE)
					node = dst_ptr[curr_node.leftChild].rightNode;

				CHECK_LEAF_POLY = CHECK_LEAF_POLY && planeTrav1 && planeTrav2;
			}
		} while (!CHECK_LEAF_POLY && node != 0u);

		if (CHECK_LEAF_POLY) {
			unsigned tr_num = DSTREE_MAX_POLY;
			if (!IS_CARVING_NODE)
				tr_num = static_cast<unsigned>(curr_node.planes[0]);

			tr_num += curr_node.leftChild;
			if (tr_num > tempMesh.lastIndID + 1) tr_num = tempMesh.lastIndID + 1;
			for (unsigned j = curr_node.leftChild; j < tr_num; ++j) {
				CRT_Hit temp_t = traceTriangle(position, direction, tempVertices, &(tempIndices[3 * tempIndicesSorted[j]]));
				if (!current_instance)
					int u = 0;
				if (temp_t.t < temp_hit.t && temp_t.t <= dirAndFar.w && temp_t.t >= posAndNear.w) {

					temp_hit.instId = current_instance;
					temp_hit.geomId = instances_info[current_instance].geomID;
					temp_hit.primId = tempIndicesSorted[j];
					temp_hit.coords[0u] = temp_t.coords[0u];
					temp_hit.coords[1u] = temp_t.coords[1u];
					temp_hit.t = temp_t.t;
					if (findAny) return temp_hit;
				}
			}
		}
	} while (node != 0u);
	if (temp_hit.geomId == 0)
		int u = 0;
	return temp_hit;
}


ISceneObject* CreateSceneRT(const char* a_impleName) { return new DSTree; }

void DeleteSceneRT(ISceneObject* a_pScene) { delete a_pScene; }