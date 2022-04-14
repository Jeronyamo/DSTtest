#include "dstree.h"


unsigned DSTree::AddGeom_Triangles4f(const LiteMath::float4* a_vpos4f, size_t a_vertNumber, const unsigned* a_triIndices, size_t a_indNumber) {
	simpleMeshInfo tempInfo;
	tempInfo.firstVertID = vertices.size();
	tempInfo.firstIndID = indices.size() / 3;

	//std::cout << "Indices size 1: " << (indices.size()/3) << std::endl;
	for (size_t j = 0; j < a_vertNumber; ++j)
		vertices.push_back(a_vpos4f[j]);

	for (size_t j = 0; j < a_indNumber; ++j)
		indices.push_back(tempInfo.firstVertID + a_triIndices[j]);

	//unsigned index = indices.size() / 3;
	tempInfo.lastVertID = vertices.size() - 1;
	tempInfo.lastIndID = indices.size() / 3;
	if (tempInfo.lastIndID) --tempInfo.lastIndID;
	//std::cout << "Indices size 2: " << (indices.size()/3 - 1) << std::endl;
	//std::cout << "Mesh1:  " << tempInfo.firstIndID << " : " << tempInfo.lastIndID << std::endl;
	//std::cout << "Index last:  " << tempInfo.firstIndID << std::endl;
	meshes.push_back(tempInfo);
	//std::cout << "Mesh " << (meshes.size() - 1) << ": " << meshes[meshes.size() - 1].firstIndID << " - " << meshes[meshes.size() - 1].lastIndID << std::endl;
	//if (index) --index;
	return meshes.size() - 1;
}

void DSTree::UpdateGeom_Triangles4f(unsigned a_geomId, const LiteMath::float4* a_vpos4f, size_t a_vertNumber, const unsigned* a_triIndices, size_t a_indNumber) {
	/*for (size_t j = 0; j < (a_vertNumber > a_indNumber ? a_vertNumber : a_indNumber); ++j) {
		if (j < a_vertNumber)
			vertices[indices[a_vertNumber * a_geomId + j]] = a_vpos4f[j];
		if (j < a_indNumber)
			indices[a_indNumber * a_geomId + j] = a_triIndices[j];
	}*/
	uint32_t firstVert = meshes[a_geomId].firstVertID;
	uint32_t firstInd = meshes[a_geomId].firstIndID;
	for (uint32_t i = 0; i < a_vertNumber; ++i) {
		vertices[i + firstVert] = a_vpos4f[i];
	}
	for (uint32_t i = 0; i < a_indNumber; ++i) {
		indices[3 * firstInd + i] = a_triIndices[i];
	}
}

uint32_t DSTree::AddInstance(uint32_t a_geomId, const LiteMath::float4x4& a_matrix) {
	//if (!a_geomId) return -1;
	simpleInstance tempInst;
	//std::cout << "Adding instance... " << a_geomId << std::endl;
	tempInst.transf = a_matrix;
	tempInst.geomID = a_geomId;
	tempInst.firstInstID = instances.size();
	uint32_t instID = instances_info.size();

	uint32_t firstIndex = 0u;

	//std::cout << "Adding  " << instances.size() << std::endl;
	for (uint32_t i = meshes[a_geomId].firstIndID; meshes[a_geomId].firstIndID <= meshes[a_geomId].lastIndID && i <= meshes[a_geomId].lastIndID; ++i) {
		instances.push_back(instances.size());
		//std::cout << "Adding2 " << i - meshes[a_geomId].firstIndID << std::endl;
	}
	tempInst.lastInstID = instances.size() - 1;
	//instances.push_back(instID);
	instances_info.push_back(tempInst);

	return instID;
}

void DSTree::UpdateInstance(uint32_t a_instanceId, const LiteMath::float4x4& a_matrix) {
	instances_info[a_instanceId].transf = a_matrix;
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

DSTree *scene;

float max(const LiteMath::float2& params) { return params[params.x < params.y]; }
float min(const LiteMath::float2& params) { return params[params.x > params.y]; }


simpleAABB DSTree::getAABB(unsigned index) {
	simpleAABB tempAABB;
	LiteMath::float3x3 vertMat = getTrVertices(index);
	for (unsigned i = 0u; i < 3u; ++i) {
		tempAABB.min[i] = min(LiteMath::float2(vertMat.row[0u][i], min(LiteMath::float2(vertMat.row[1u][i], vertMat.row[2u][i]))));
		tempAABB.max[i] = max(LiteMath::float2(vertMat.row[0u][i], max(LiteMath::float2(vertMat.row[1u][i], vertMat.row[2u][i]))));
	}
	return tempAABB;
}

simpleAABB DSTree::mergeAABB(unsigned index, simpleAABB aabb) {
	LiteMath::float3x3 vertMat = getTrVertices(index);
	for (unsigned i = 0u; i < 3u; ++i) {
		aabb.min[i] = min(LiteMath::float2(min(LiteMath::float2(aabb.min[i], vertMat.row[0u][i])), min(LiteMath::float2(vertMat.row[1u][i], vertMat.row[2u][i]))));
		aabb.max[i] = max(LiteMath::float2(max(LiteMath::float2(aabb.max[i], vertMat.row[0u][i])), max(LiteMath::float2(vertMat.row[1u][i], vertMat.row[2u][i]))));
	}
	return aabb;
}

/*float* DSTree::getSceneAABB() {
	sceneAABB_arr[0u] = scene_AABB.min.x;
	sceneAABB_arr[1u] = scene_AABB.min.y;
	sceneAABB_arr[2u] = scene_AABB.min.z;
	sceneAABB_arr[3u] = scene_AABB.max.x;
	sceneAABB_arr[4u] = scene_AABB.max.y;
	sceneAABB_arr[5u] = scene_AABB.max.z;
	return sceneAABB_arr;
}*/

uint32_t DSTree::getTrIndex(uint32_t instTrIndex) {
	uint32_t instID = 0;
	while (instTrIndex > instances_info[instID].lastInstID) ++instID;

	return instTrIndex - instances_info[instID].firstInstID + meshes[instances_info[instID].geomID].firstIndID;
}

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
}


float getMinPos(const void* index) {
	unsigned* temp_index = (unsigned*) index;

	LiteMath::float3 coords = scene->getTrVerticesAxis(*temp_index);

	return min(LiteMath::float2(coords.x, min(LiteMath::float2(coords.y, coords.z))));
}

float getMaxPos(const void* index) {
	unsigned* temp_index = (unsigned*) index;

	LiteMath::float3 coords = scene->getTrVerticesAxis(*temp_index);

	return max(LiteMath::float2(coords.x, max(LiteMath::float2(coords.y, coords.z))));
}

int cmpMinPos(const void* ind1, const void* ind2) {
	float min1 = getMinPos(ind1);
	float min2 = getMinPos(ind2);
	int res = int(min1 > min2) - int(min1 < min2);

	if (!res) res = int(ind1 > ind2) - int(ind1 < ind2);
	return res;
}

int cmpMaxPos(const void* ind1, const void* ind2) {
	float max1 = getMaxPos(ind1);
	float max2 = getMaxPos(ind2);
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
		qsort(instances_ptr, elem_count, sizeof(unsigned), cmpMaxPos);

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
		qsort(instances_ptr, elem_count, sizeof(unsigned), cmpMaxPos);
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

	//std::cout << "Node " << curr_node_index << ", " << first << "-" << last << std::endl;
	dst_nodes.nodes.emplace_back();
	dst_nodes.nodes[curr_node_index].leftChild = (first << 6u) + 1u;
	dst_nodes.nodes[curr_node_index].planes[0] = static_cast<float>(last - first + 1u);

	if (last - first >= DSTREE_MAX_POLY) {
		simpleAABB left_aabb, right_aabb, curr_aabb;
		unsigned new_axis = 3u;
		float SAH = calculateSAH(new_axis, new_index, &(instances[first]),
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

	instances_size = instances.size();
	instances_info_size = instances_info.size();
	std::cout << "Triangles: " << indices.size() / 3 << std::endl;
	std::cout << "Instances: " << instances_size << std::endl;
	std::fstream dstFile;
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

			if (!dstFile.eof()) dst_nodes.nodes.push_back(tempNode);
		}

		for (uint32_t k = 0u; k < instances_size; ++k) {
			dstFile >> instances[k];
		}
	}
	if (!dst_nodes.nodes.size()) {
		std::cout << "The file was empty" << std::endl;

		dstFile.open("./dstree", std::fstream::out | std::fstream::trunc);
		scene = this;
		dstBuilderRecur(0u, instances_size - 1u, SAH_MAX);
		scene = nullptr;
		
		scene_AABB = getAABB(instances[0u]);
		for (uint32_t i = 1u; i < instances_size; ++i) {
			scene_AABB = mergeAABB(instances[i], scene_AABB);
		}

		dstFile << std::setprecision(16) << scene_AABB.min.x << " ";
		dstFile << std::setprecision(16) << scene_AABB.min.y << " ";
		dstFile << std::setprecision(16) << scene_AABB.min.z << " ";
		dstFile << std::setprecision(16) << scene_AABB.max.x << " ";
		dstFile << std::setprecision(16) << scene_AABB.max.y << " ";
		dstFile << std::setprecision(16) << scene_AABB.max.z << " \n\n";
		dstFile << dst_nodes.nodes.size() << " \n\n\n";

		for (uint32_t i = 0u; i < dst_nodes.nodes.size(); ++i) {
			dstFile << dst_nodes.nodes[i].leftChild << " ";
			dstFile << dst_nodes.nodes[i].rightNode << " ";
			dstFile << std::setprecision(16) << dst_nodes.nodes[i].planes[0] << " ";
			dstFile << std::setprecision(16) << dst_nodes.nodes[i].planes[1] << " \n\n";
		}

		for (uint32_t j = 0u; j < instances.size(); ++j) {
			dstFile << instances[j] << " \n";
		}
	}
	std::cout << "Nodes: " << dst_nodes.nodes.size() << std::endl;
	dstFile.close();
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
	tempHitInfo.t = FLT_MAX;
	tempHitInfo.primId = -1;
	tempHitInfo.geomId = -1;
	tempHitInfo.instId = -1;

	LiteMath::float3x3 vertMat = getTrVertices(instances[trIndex]);
	const LiteMath::float3  a = vertMat.row[0u];
	const LiteMath::float3 E1 = vertMat.row[1u] - a;
	const LiteMath::float3 E2 = vertMat.row[2u] - a;

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
			tempHitInfo.instId = getInstID(instances[trIndex]);
			tempHitInfo.geomId = instances_info[tempHitInfo.instId].geomID;
			tempHitInfo.primId = instances[trIndex] - instances_info[tempHitInfo.instId].firstInstID;
			tempHitInfo.coords[0u] = temp_u;
			tempHitInfo.coords[1u] = temp_v;
			tempHitInfo.t = temp_t;
		}
	}
	return tempHitInfo;
}

CRT_Hit DSTree::RayQuery_NearestHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar) {
	return findHit(posAndNear, dirAndFar, false);
}

bool DSTree::RayQuery_AnyHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar) {
	CRT_Hit temp_hit = findHit(posAndNear, dirAndFar, true);

	return (temp_hit.t >= posAndNear.w) && (temp_hit.t <= dirAndFar.w);
}

CRT_Hit DSTree::findHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, bool findAny) {
	LiteMath::float3 position  = LiteMath::float3(posAndNear.x, posAndNear.y, posAndNear.z);
	LiteMath::float3 direction = LiteMath::float3( dirAndFar.x,  dirAndFar.y,  dirAndFar.z);
	LiteMath::float3 invDir    = LiteMath::float3(1.f / dirAndFar.x, 1.f / dirAndFar.y, 1.f / dirAndFar.z);
	CRT_Hit temp_hit;
	temp_hit.t = FLT_MAX;
	temp_hit.primId = -1;
	temp_hit.geomId = -1;
	temp_hit.instId = -1;

	if (!traceAABB(position, direction, invDir)) return temp_hit;

	unsigned node = 0u;
	do {
		DSNode curr_node;
		bool CHECK_LEAF_POLY, IS_CARVING_NODE, IS_DOUBLE_CARVE;
		do {
			curr_node = dst_nodes.nodes[node];
			node = curr_node.rightNode;
			unsigned header = curr_node.leftChild & 63u;
			curr_node.leftChild = curr_node.leftChild >> 6;
			CHECK_LEAF_POLY = header & 1u;
			IS_CARVING_NODE = header & 2u;
			IS_DOUBLE_CARVE = !((header & 48u) == 32u) && IS_CARVING_NODE;

			if (!CHECK_LEAF_POLY || IS_CARVING_NODE) {
				LiteMath::uint2 planes = LiteMath::uint2((header & 12u) >> 2);
				//was bool is_normal_positive
				bool planeTrav2 = IS_CARVING_NODE;
				bool planeTrav1 = !planeTrav2;
				if (IS_DOUBLE_CARVE) {
					//is_normal_positive = bool2(bool(planes.x & 2u), bool(planes.x & 1u));
					planeTrav1 = bool(planes.x & 2u);
					planeTrav2 = bool(planes.x & 1u);
					planes = LiteMath::uint2((header >> 5) & 1u, ((header >> 4) & 1u) + 1u);
				}

				planeTrav1 = planeTrav1 == (position[planes.x] < curr_node.planes[0u]);
				planeTrav2 = planeTrav2 == (position[planes.y] < curr_node.planes[1u]);
				if (!planeTrav1 && (direction[planes.x] < -LiteMath::EPSILON || direction[planes.x] > LiteMath::EPSILON)) {
					float temp_t = invDir[planes.x] * (curr_node.planes[0u] - position[planes.x]);
					planeTrav1 = temp_t >= 0.f && temp_t < temp_hit.t;
				}
				if (!planeTrav2 && (direction[planes.y] < -LiteMath::EPSILON || direction[planes.y] > LiteMath::EPSILON)) {
					float temp_t = invDir[planes.y] * (curr_node.planes[1u] - position[planes.y]);
					planeTrav2 = temp_t >= 0.f && temp_t < temp_hit.t;
				}
				if (planeTrav1 && (!IS_CARVING_NODE || (planeTrav2 && !CHECK_LEAF_POLY)))
					node = curr_node.leftChild;
				if (!planeTrav1 && planeTrav2 && !IS_CARVING_NODE)
					node = dst_nodes.nodes[curr_node.leftChild].rightNode;

				CHECK_LEAF_POLY = CHECK_LEAF_POLY && planeTrav1 && planeTrav2;
			}
		} while (!CHECK_LEAF_POLY && node != 0u);

		if (CHECK_LEAF_POLY) {
			unsigned tr_num = DSTREE_MAX_POLY;
			if (!IS_CARVING_NODE)
				tr_num = static_cast<unsigned>(curr_node.planes[0]);

			tr_num += curr_node.leftChild;
			if (tr_num > instances_size) tr_num = instances_size;
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