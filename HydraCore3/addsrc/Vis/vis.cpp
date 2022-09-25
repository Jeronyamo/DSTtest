#include "vis.h"

#define MW_VER_MAJ 4
#define MW_VER_MIN 2
#define MW_RESIZABLE GL_TRUE
#define MW_WIDTH  1024
#define MW_HEIGHT 1024
#define MW_TITLE "DST Visualizer"


void fbSizeCallback(GLFWwindow* window, int width, int height) {
	glViewport(0, 0, width, height);
}

void mouseCallback(GLFWwindow* window, int button, int action, int mods) {
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		double x_pos, y_pos;

		glfwGetCursorPos(window, &x_pos, &y_pos);
	}
}

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, GL_TRUE);
	}
}


Visualizer::Visualizer() {
	if (!glfwInit()) {
		std::cout << "ERROR:: Failed to init GLFW" << std::endl;
		exit(-1);
	}

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, MW_VER_MAJ);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, MW_VER_MIN);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, MW_RESIZABLE);

	window = glfwCreateWindow(MW_WIDTH, MW_HEIGHT, MW_TITLE, NULL, NULL);
	if (window == nullptr) {
		std::cout << "ERROR:: Failed to create GLFW window" << std::endl;
		exit(-1);
	}
	glfwMakeContextCurrent(window);

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
		std::cout << "ERROR:: Failed to initialize OpenGL context" << std::endl;
		exit(-1);
	}

	int fb_width, fb_height;
	glfwSetMouseButtonCallback(window, mouseCallback);
	glfwSetKeyCallback(window, keyCallback);
	glfwGetFramebufferSize(window, &fb_width, &fb_height);
	glViewport(0, 0, fb_width, fb_height);
}

Visualizer::~Visualizer() {
	glfwDestroyWindow(window);
}

void Visualizer::create_buffer() {

}

void Visualizer::create_program() {

}

void Visualizer::start() {
	while (!glfwWindowShouldClose(window)) {
		glfwSwapBuffers(window);
		glfwPollEvents();
	}
}