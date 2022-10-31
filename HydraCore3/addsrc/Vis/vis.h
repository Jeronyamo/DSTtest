#pragma once

#include <vector>

#include "../HydraCore3/ext/common.h"
#include "../../ext/imgui/imgui.h"
#include "../../ext/imgui/backends/imgui_impl_glfw.h"
#include "../../ext/imgui/backends/imgui_impl_opengl3.h"


class Visualizer {
	struct VisIvec2 { 
		int v[2];

		int& operator[](int ind) { return v[ind]; }
	};

	struct VisMesh {
		VisIvec2 triVAO;
		GLuint triVBO, triIBO;
		std::vector <VisIvec2> treeVAOs; //[0] - buffers for tree layers (squares/boxes); [1] - number of elements (indices) for glDrawElements()
		std::vector <GLuint> treeVBOs;
	};

	struct VisInstance {
		unsigned mesh;
		float M[16];
		bool visualize = false;
		unsigned min, max; //tree layers range to visualize
	};


	GLFWwindow *window = nullptr;
	GLuint triProg = 0u, boxProg = 0u;// , lineProg = 0u;
	std::vector <VisMesh> meshes;
	std::vector <VisInstance> instances;
	ImGuiIO imgui_io;


	GLuint create_buffer(GLenum buffer_type = GL_ARRAY_BUFFER, void *buffer_data = nullptr, size_t buf_size_bytes = 0, GLenum draw_type = GL_STATIC_DRAW, bool bind = false, GLuint bind_index = 0);
	GLuint create_array();
	GLuint create_shader_from_file(GLenum shader_type, const char* shader_path);
	GLuint create_shader(GLenum shader_type, char *shader, bool path_to_file);

	char *read_file(const char *path);
	int link_program(GLuint program);
	void ImGUI_initialization();

public:
	Visualizer();
	~Visualizer();

	unsigned addMesh(float* tri_buf, size_t tri_buf_size_bytes,
					 unsigned* tri_ind_buf, size_t tri_ind_buf_size_bytes,
					 std::vector <float> *tree_buf = nullptr);
	void addInstance(unsigned mesh_id, float* matrix);
	void start();
};