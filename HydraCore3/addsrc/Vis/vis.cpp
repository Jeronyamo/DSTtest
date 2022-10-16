#include "vis.h"

#define VIS_WIDTH  1024
#define VIS_HEIGHT 1024
#define VIS_NEAR 0.1f
#define VIS_FAR 300.f
#define VIS_RESIZABLE GL_TRUE

#define VIS_VER_MAJ 4
#define VIS_VER_MIN 2
#define VIS_IMGUI_GLSL_VERSION "#version 420"  //should be changed together with MW_VER_MAJ, MW_VER_MIN

#define VIS_TITLE "DST Visualizer"


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

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, VIS_VER_MAJ);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, VIS_VER_MIN);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, VIS_RESIZABLE);

	window = glfwCreateWindow(VIS_WIDTH, VIS_HEIGHT, VIS_TITLE, NULL, NULL);
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
	glEnable(GL_CULL_FACE); GL_CHECK_ERRORS;

	//load shaders here
	triProg = glCreateProgram();

	GLuint triVert = create_shader(GL_VERTEX_SHADER, "./shaders/triVert.glsl", true); GL_CHECK_ERRORS;
	GLuint triFrag = create_shader(GL_FRAGMENT_SHADER, "./shaders/triFrag.glsl", true); GL_CHECK_ERRORS;

	glAttachShader(triProg, triVert); GL_CHECK_ERRORS;
	glAttachShader(triProg, triFrag); GL_CHECK_ERRORS;
	glLinkProgram(triProg); GL_CHECK_ERRORS;

	glDeleteShader(triVert); GL_CHECK_ERRORS;
	glDeleteShader(triFrag); GL_CHECK_ERRORS;

	//projection matrix
	float tanHalfFovy = tan(20.f);
	P[0] = 1.f / (VIS_WIDTH / VIS_HEIGHT * tanHalfFovy);
	P[5] = 1.f / tanHalfFovy;
	P[11] = -1.f;
	P[10] = -(VIS_FAR + VIS_NEAR) / (VIS_FAR - VIS_NEAR);
	P[14] = -(2.f * VIS_FAR * VIS_NEAR) / (VIS_FAR - VIS_NEAR);
	P[1] = P[2] = P[3] = P[4] = P[6] = P[7] = P[8] = P[9] = P[12] = P[13] = P[15] = 0.f;
}


Visualizer::~Visualizer() {
	glfwDestroyWindow(window);
	glfwTerminate();
}


unsigned Visualizer::addMesh(float* tri_buf, size_t tri_buf_size_bytes, unsigned *tri_ind_buf, size_t tri_ind_buf_size_bytes, std::vector <float> *tree_buf) {
	VisMesh new_mesh;

	new_mesh.triVAO[0] = create_array();
	new_mesh.triVAO[1] = tri_ind_buf_size_bytes / (3 * sizeof(unsigned));
	new_mesh.triVBO = create_buffer(GL_ARRAY_BUFFER, tri_buf, tri_buf_size_bytes, GL_STATIC_DRAW);
	new_mesh.triIBO = create_buffer(GL_ELEMENT_ARRAY_BUFFER, tri_ind_buf, tri_ind_buf_size_bytes, GL_STATIC_DRAW);


	glBindVertexArray(new_mesh.triVAO[0]);

	glBindBuffer(GL_ARRAY_BUFFER, new_mesh.triVBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, new_mesh.triIBO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), nullptr);
	glEnableVertexAttribArray(0);

	glBindVertexArray(0); GL_CHECK_ERRORS;
	/*
	for (unsigned i = 0u; i < tree_buf.size(); ++i) {
		GLuint VAO = create_array();
		GLuint VBO = create_buffer(GL_ARRAY_BUFFER, tree_buf[i].data(), tree_buf[i].size() * sizeof(float), GL_STATIC_DRAW);
		//CONTINUE HERE - LOAD DUAL-SPLIT TREE
	}
	*/

	meshes.push_back(new_mesh);
	return meshes.size() - 1;
}

void Visualizer::addInstance(unsigned mesh_id, float *matrix) {
	VisInstance new_instance;
	new_instance.mesh = mesh_id;
	new_instance.visualize = false;
	new_instance.min = 0;
	new_instance.max = 0;

	for (int i = 0; i < 16; ++i) {
		new_instance.M[i] = matrix[i];
	}
	instances.push_back(new_instance);
	return;
}


GLuint Visualizer::create_buffer(GLenum buffer_type, void* buffer_data, size_t buf_size_bytes, GLenum draw_type, bool bind, GLuint bind_index) {
	GLuint BO = 0;
	glCreateBuffers(1, &BO);

	if (buffer_data != nullptr) {
		glBindBuffer(buffer_type, BO);
		glBufferData(buffer_type, buf_size_bytes, buffer_data, draw_type);
		if (bind) glBindBufferBase(buffer_type, bind_index, BO);
		glBindBuffer(buffer_type, 0);
	}

	return BO;
}

GLuint Visualizer::create_array() {
	GLuint VAO = 0;
	glGenVertexArrays(1, &VAO);
	return VAO;
}


char* Visualizer::read_file(char* path) {
	std::ifstream shader_file(path);

	if (!shader_file.is_open()) {
		std::cout << "ERROR:: Could not open a file: " << path << std::endl;
		return nullptr;
	}

	shader_file.seekg(0, shader_file.end);

	unsigned length = shader_file.tellg();
	char* file = new char[length + 1u];

	for (unsigned i = 0; i <= length; ++i)
		file[i] = '\0';

	shader_file.seekg(0, shader_file.beg);
	shader_file.clear();
	shader_file.read(file, length);
	shader_file.close();
	file[length] = '\0';
	return file;
}

GLuint Visualizer::create_shader(GLenum shader_type, char* shader_inp, bool path_to_file) {
	char *shader = shader_inp;
	if (path_to_file) shader = read_file(shader_inp);


	GLuint shader_obj = glCreateShader(shader_type);
	glShaderSource(shader_obj, 1, &shader, NULL);
	glCompileShader(shader_obj); GL_CHECK_ERRORS;

	int success = 0;
	char infoLog[1024];
	glGetShaderiv(shader_obj, GL_COMPILE_STATUS, &success);
	if (!success) {
		glGetShaderInfoLog(shader_obj, 1024, NULL, infoLog);
		std::cout << "ERROR:: Shader compilation error:\n" << infoLog << std::endl;
		return -1;
	}

	return shader_obj;
}

int Visualizer::link_program(GLuint program) {
	glLinkProgram(program);

	int success = 0;
	char infoLog[512];
	glGetProgramiv(program, GL_LINK_STATUS, &success);
	if (!success) {
		glGetProgramInfoLog(program, 512, NULL, infoLog);
		std::cout << "ERROR:: Program linking error:\n" << infoLog << std::endl;
		return -1;
	}
	return 0;
}


void Visualizer::ImGUI_initialization() {
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;

	ImGui::StyleColorsDark();

	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init(VIS_IMGUI_GLSL_VERSION);
}


void Visualizer::start() {
	ImGUI_initialization();

	ImVec4 clear_color = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);
	unsigned num_meshes = instances.size();
	bool *active_meshes = new bool[num_meshes];
	for (unsigned i = 0u; i < num_meshes; ++i)
		active_meshes[i] = false;

	glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();

		//Draw meshes here
		glUseProgram(triProg);
		glUniformMatrix4fv(glGetUniformLocation(triProg, "P"), 1, true, P);
		glUniformMatrix4fv(glGetUniformLocation(triProg, "V"), 1, false, camera.V);

		for (unsigned i = 0u; i < num_meshes; ++i) {
			if (active_meshes[i]) {
				VisIvec2 tmpVAO{ meshes[instances[i].mesh].triVAO[0], meshes[instances[i].mesh].triVAO[1] };
				glBindVertexArray(meshes[i].triVAO[0]);
				glUniformMatrix4fv(glGetUniformLocation(triProg, "M"), 1, false, instances[i].M);
				glDrawElements(GL_TRIANGLES, meshes[i].triVAO[1], GL_UNSIGNED_INT, nullptr);
				glBindVertexArray(0);
			}
		}
		glUseProgram(0);

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
		glViewport(0, 0, display_w, display_h); GL_CHECK_ERRORS;
		glClear(GL_COLOR_BUFFER_BIT);
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
	}

	delete[] active_meshes;

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
}

void Visualizer::VisCam::LookAtUpd() {
	ImVec4 x, y, z;
	z.x = CamPos.x - LookAt.x;
	z.y = CamPos.y - LookAt.y;
	z.z = CamPos.z - LookAt.z;

	float norm = sqrt(z.x * z.x + z.y * z.y + z.z * z.z);
	if (norm > 10e-6 || norm < -10e-6)
		z = ImVec4(z.x / norm, z.y / norm, z.z / norm, 0.f);
	y = UpVect;

	x = ImVec4(y.y * z.z - y.z * z.y, y.z * z.x - y.x * z.z, y.x * z.y - y.y * z.x, 0.f);
	y = ImVec4(z.y * x.z - z.z * x.y, z.z * x.x - z.x * x.z, z.x * x.y - z.y * x.x, 0.f);

	norm = sqrt(x.x * x.x + x.y * x.y + x.z * x.z);
	if (norm > 10e-6 || norm < -10e-6)
		x = ImVec4(x.x / norm, x.y / norm, x.z / norm, 0.f);
	
	norm = sqrt(y.x * y.x + y.y * y.y + y.z * y.z);
	if (norm > 10e-6 || norm < -10e-6)
		y = ImVec4(y.x / norm, y.y / norm, y.z / norm, 0.f);

	float M[16] = { x.x, y.x, z.x, 0.0f,
					x.y, y.y, z.y, 0.0f,
					x.z, y.z, z.z, 0.0f,
					-x.x * CamPos.x - x.y * CamPos.y - x.z * CamPos.z,
					-y.x * CamPos.x - y.y * CamPos.y - y.z * CamPos.z,
					-z.x * CamPos.x - z.y * CamPos.y - z.z * CamPos.z,
					1.0f };
	
	for (int i = 0; i < 16; ++i)
		V[i] = M[i];
	return;
}
