#pragma once

#include <vector>
#include "../HydraCore3/ext/common.h"


class Visualizer {
	GLFWwindow *window;
	std::vector <GLuint> buffers;


public:
	Visualizer();
	~Visualizer();

	void create_buffer();
	void create_program();
	void start();
};