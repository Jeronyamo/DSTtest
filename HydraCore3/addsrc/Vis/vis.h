#pragma once

#include <vector>
#include "../HydraCore3/ext/common.h"
#include "../../ext/imgui/imgui.h"
#include "../../ext/imgui/backends/imgui_impl_glfw.h"
#include "../../ext/imgui/backends/imgui_impl_opengl3.h"


enum VIS_BUFFER_TYPE { LINE, BOX }; //squares are treated like boxes for now, will change if problems occur

class Visualizer {
	GLFWwindow *window;
	std::vector <GLuint> buffers;


	void create_buffer();
	void create_program();

	void ImGUI_initialization();

public:
	Visualizer();
	~Visualizer();

	void addMesh(float* tri_buf, size_t tri_buf_size_bytes, std::vector <std::vector <float>>& tree_buf, VIS_BUFFER_TYPE tree_buf_type);
	void start();
};