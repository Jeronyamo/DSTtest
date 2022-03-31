#pragma once

#include <vector>


//16 bytes
struct DSNode {
	unsigned leftChild = 0u;
	unsigned rightNode = 0u;
	float    planes[2] = { 0.f, 0.f };
};

struct DSTNodeInfo {
	unsigned trOffset;
	unsigned length;
};

struct DST {
	std::vector <DSNode> nodes;
	std::vector <DSTNodeInfo> nodeTriangles;
};