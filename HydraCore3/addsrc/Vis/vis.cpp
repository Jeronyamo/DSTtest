#include "vis.h"

#define MW_WIDTH  1024
#define MW_HEIGHT 1024
#define MW_RESIZABLE GL_TRUE

#define MW_VER_MAJ 4
#define MW_VER_MIN 2
#define MW_GLSL_VERSION "#version 420"  //should be changed together with MW_VER_MAJ, MW_VER_MIN

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
	glfwTerminate();
}


void Visualizer::addMesh(float* tri_buf, size_t tri_buf_size_bytes, std::vector <std::vector <float>> &tree_buf, VIS_BUFFER_TYPE tree_buf_type) {

}


void Visualizer::create_buffer() {

}


void Visualizer::create_program() {

}


void Visualizer::ImGUI_initialization() {
	const char* glsl_version = MW_GLSL_VERSION;

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;

	ImGui::StyleColorsDark();

	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init(glsl_version);
}


void Visualizer::start() {
	ImGUI_initialization();

	ImVec4 clear_color = ImVec4(0.1f, 0.1f, 0.1f, 1.0f);
	unsigned num_meshes = 11u;
	bool *active_meshes = new bool[num_meshes];
	for (unsigned i = 0u; i < num_meshes; ++i)
		active_meshes[i] = false;

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();

		//Draw meshes here

		{
			ImGui_ImplOpenGL3_NewFrame();
			ImGui_ImplGlfw_NewFrame();
			ImGui::NewFrame();

			ImGui::Begin("DST visualizer");
			ImGui::Text("Active meshes:   ");

			for (unsigned i = 0u; i < num_meshes; ++i) {
				if (ImGui::Selectable((std::string("   ") + std::to_string(i) + std::string("   ")).c_str(), active_meshes[i], 0, ImVec2(55, 30)))
					active_meshes[i] = !active_meshes[i];
			}

			ImGui::SameLine();
			ImGui::End();
		}

		ImGui::Render();
		int display_w, display_h;
		glfwGetFramebufferSize(window, &display_w, &display_h);
		glViewport(0, 0, display_w, display_h);
		glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
		glClear(GL_COLOR_BUFFER_BIT);
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
	}

	delete[] active_meshes;

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
}