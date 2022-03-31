#pragma once

#include "../../Include/glm/glm/vec2.hpp"
#include "../../Include/glm/glm/vec3.hpp"

/*
struct simpleVec2 {
	float x, y;

	simpleVec2() : x(0.f), y(0.f) {}
	simpleVec2(const float& new_x, const float& new_y) : x(new_x), y(new_y) {}
	simpleVec2(const simpleVec2& tex_vert) : x(tex_vert.x), y(tex_vert.y) {}
	float &operator[](const unsigned ind) { if (ind == 1u) return y; return x; }
	const float &operator[](const unsigned ind) const { if (ind == 1u) return y; return x; }
};

struct simpleVec3 {
	float x, y, z;

	simpleVec3() : x(0.f), y(0.f), z(0.f) {}
	simpleVec3(const float& new_x, const float& new_y, const float& new_z) : x(new_x), y(new_y), z(new_z) {}
	simpleVec3(const simpleVec3& tex_vert) : x(tex_vert.x), y(tex_vert.y), z(tex_vert.z) {}
	float &operator[](const unsigned ind) { if (ind == 1u) return y; if (ind == 2u) return z; return x; }
	const float &operator[](const unsigned ind) const { if (ind == 1u) return y; if (ind == 2u) return z; return x; }
};

struct simpleUVec3 {
	unsigned x, y, z;

	simpleUVec3() : x(0u), y(0u), z(0u) {}
	simpleUVec3(const unsigned& new_x, const unsigned& new_y, const unsigned& new_z) : x(new_x), y(new_y), z(new_z) {}
	simpleUVec3(const simpleUVec3& tex_vert) : x(tex_vert.x), y(tex_vert.y), z(tex_vert.z) {}
	unsigned &operator[](const unsigned ind) { if (ind == 1u) return y; if (ind == 2u) return z; return x; }
	const unsigned &operator[](const unsigned ind) const { if (ind == 1u) return y; if (ind == 2u) return z; return x; }
};

struct simpleAABB {
	glm::vec3 min;
	glm::vec3 max;

	simpleAABB() : min(), max() {}
	simpleAABB(const glm::vec3& new_min, const glm::vec3& new_max) : min(new_min), max(new_max) {}
};*/

struct simpleTriangle {
	glm::vec3 vertex1;
	glm::vec3 vertex2;
	glm::vec3 vertex3;

	simpleTriangle() : vertex1(), vertex2(), vertex3() {}
	simpleTriangle(const glm::vec3& new_vert1, const glm::vec3& new_vert2, const glm::vec3& new_vert3) : vertex1(new_vert1), vertex2(new_vert2), vertex3(new_vert3) {}
};