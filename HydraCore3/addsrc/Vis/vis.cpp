#include "vis.h"
#include <math.h>

#define VIS_VER_MAJ 4
#define VIS_VER_MIN 2
#define VIS_IMGUI_GLSL_VERSION "#version 420"  //should be changed together with MW_VER_MAJ, MW_VER_MIN

#define VIS_TITLE "DST Visualizer"



struct VisCamera {
	float V[16] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1. };
	ImVec4 CamPos{ 0.f, 0.f, 1.f, 0.f };
	ImVec4 LookAt{ 0.f, 0.f, 0.f, 0.f };
	ImVec4 UpVect{ 0.f, 1.f, 0.f, 0.f };
	float dist = 16.f, yaw = 0.f, pitch = 0.f;

	bool update = true;
	void CamPosUpd();
	void LookAtUpd();
} camera;

struct VisProj {
	float P[16] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
	unsigned width = 1024u, height = 1024u;
	float near = 0.5f, far = 300.f, fov_deg = 50.f;

	bool update = true;
	void ProjUpd();
} proj;

struct VisMouse {
	bool LMB = false, RMB = false, MMB = false;
	double xpos = 0., ypos = 0.;
	double timedel = 0.;
} mouse;


void fbSizeCallback(GLFWwindow* window, int width, int height) {
	//glViewport(0, 0, width, height);
	//proj.width = width;
	//proj.height = height;
	//proj.update = true;
}

void mouseCallback(GLFWwindow* window, int button, int action, int mods) {
	if ((button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_RIGHT || button == GLFW_MOUSE_BUTTON_MIDDLE)
		&& action == GLFW_PRESS && !(mouse.LMB || mouse.RMB || mouse.MMB) && !ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow)) {

		if (button == GLFW_MOUSE_BUTTON_LEFT)
			mouse.LMB = true;
		if (button == GLFW_MOUSE_BUTTON_RIGHT)
			mouse.RMB = true;
		if (button == GLFW_MOUSE_BUTTON_MIDDLE)
			mouse.MMB = true;

		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
		glfwGetCursorPos(window, &mouse.xpos, &mouse.ypos);
	}
	if ((button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_RIGHT || button == GLFW_MOUSE_BUTTON_MIDDLE)
		&& action == GLFW_RELEASE) {

		if (button == GLFW_MOUSE_BUTTON_LEFT)
			mouse.LMB = false;
		if (button == GLFW_MOUSE_BUTTON_RIGHT)
			mouse.RMB = false;
		if (button == GLFW_MOUSE_BUTTON_MIDDLE)
			mouse.MMB = false;

		if (!(mouse.LMB || mouse.RMB || mouse.MMB))
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
	}
}

void mouseMoveCallback(GLFWwindow* window, double xpos, double ypos) {
	double xdel = mouse.xpos - xpos, ydel = mouse.ypos - ypos;

	if (mouse.RMB) {
		camera.LookAt.y += ydel * mouse.timedel;
		camera.dist -= xdel * mouse.timedel;

		if (camera.dist < 2.f)
			camera.dist = 2.f + 10e-6;

		glfwSetCursorPos(window, mouse.xpos, mouse.ypos);
		camera.update = true;
	}

	if (mouse.LMB){
		camera.yaw -= xdel * mouse.timedel * 0.1f;
		camera.pitch -= ydel * mouse.timedel * 0.1f;

		if (camera.pitch >= 3.1415f * 0.49f)
			camera.pitch = 3.1415f * 0.49f;
		if (camera.pitch <= 3.1415f * -0.49f)
			camera.pitch = 3.1415f * -0.49f;
		if (camera.yaw <= -3.1415f * 2.f || camera.yaw >= 3.1415f * 2.f)
			camera.yaw = 0.f;
		//if (camera.yaw >= 360.0f)
		//	camera.yaw -= 360.f;
		//if (camera.yaw < -360.0f)
		//	camera.yaw += 360.0f;

		glfwSetCursorPos(window, mouse.xpos, mouse.ypos);
		camera.update = true;
	}

	if (mouse.MMB) {
		ImVec2 forwVec{ camera.LookAt.x - camera.CamPos.x, camera.LookAt.z - camera.CamPos.z};
		float norm = sqrt(forwVec.x * forwVec.x + forwVec.y * forwVec.y);
		forwVec = ImVec2(forwVec.x / norm * mouse.timedel, forwVec.y / norm * mouse.timedel);

		camera.LookAt.x += ydel * forwVec.x + xdel * forwVec.y;
		camera.LookAt.z += ydel * forwVec.y - xdel * forwVec.x;

		glfwSetCursorPos(window, mouse.xpos, mouse.ypos);
		camera.update = true;
	}
}

void scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
	if (!ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow)) {
		proj.fov_deg -= (float)yoffset;
		if (proj.fov_deg < 10.0f)
			proj.fov_deg = 10.0f;
		if (proj.fov_deg > 65.0f)
			proj.fov_deg = 65.0f;
		proj.update = true;
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
	glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

	window = glfwCreateWindow(proj.width, proj.height, VIS_TITLE, NULL, NULL);
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
	glfwSetCursorPosCallback(window, mouseMoveCallback);
	glfwSetScrollCallback(window, scrollCallback);
	glfwSetKeyCallback(window, keyCallback);
	//glfwSetFramebufferSizeCallback(window, fbSizeCallback);
	glfwGetFramebufferSize(window, &fb_width, &fb_height);
	glViewport(0, 0, fb_width, fb_height);
	//glEnable(GL_CULL_FACE); GL_CHECK_ERRORS;

	//load shaders here
	triProg = glCreateProgram();
	GLuint triVert = create_shader_from_file(  GL_VERTEX_SHADER, "./shaders/triVert.glsl"); GL_CHECK_ERRORS;
	GLuint triFrag = create_shader_from_file(GL_FRAGMENT_SHADER, "./shaders/triFrag.glsl"); GL_CHECK_ERRORS;

	glAttachShader(triProg, triVert); GL_CHECK_ERRORS;
	glAttachShader(triProg, triFrag); GL_CHECK_ERRORS;
	glLinkProgram (triProg); GL_CHECK_ERRORS;

	glDeleteShader(triVert); GL_CHECK_ERRORS;
	glDeleteShader(triFrag); GL_CHECK_ERRORS;
}


Visualizer::~Visualizer() {
	glfwDestroyWindow(window);
	glfwTerminate();
}


void Visualizer::addTree(std::vector <std::vector <float>>* tree_buf, std::vector <std::vector <unsigned>>* layers_inds) {
	tree_buffer = tree_buf;
	layers_indices = layers_inds;

	boxProg = glCreateProgram();
	GLuint treeVert = create_shader_from_file(GL_VERTEX_SHADER, "./shaders/treeVert.glsl"); GL_CHECK_ERRORS;
	GLuint treeGeom = create_shader_from_file(GL_FRAGMENT_SHADER, "./shaders/treeGeom.glsl"); GL_CHECK_ERRORS;
	GLuint treeFrag = create_shader_from_file(GL_FRAGMENT_SHADER, "./shaders/treeFrag.glsl"); GL_CHECK_ERRORS;

	glAttachShader(boxProg, treeVert); GL_CHECK_ERRORS;
	glAttachShader(boxProg, treeGeom); GL_CHECK_ERRORS;
	glAttachShader(boxProg, treeFrag); GL_CHECK_ERRORS;
	glLinkProgram (boxProg); GL_CHECK_ERRORS;

	glDeleteShader(treeVert); GL_CHECK_ERRORS;
	glDeleteShader(treeGeom); GL_CHECK_ERRORS;
	glDeleteShader(treeFrag); GL_CHECK_ERRORS;

		//CONTINUE HERE - LOAD DUAL-SPLIT TREE
		//CREATE ONE BUFFER AND MANY VAO's FOR DIFFERENT PARTS OF IT!!!
	GLuint VAO = create_array();
	for (unsigned i = 0u; i < tree_buffer->size(); ++i) {
		GLuint VBO = create_buffer(GL_ARRAY_BUFFER, tree_buf[i].data(), tree_buf[i].size() * sizeof(float), GL_STATIC_DRAW);
		
	}
}

unsigned Visualizer::addMesh(float* tri_buf, size_t tri_buf_size_bytes, unsigned *tri_ind_buf, size_t tri_ind_buf_size_bytes) {
	VisMesh new_mesh;

	new_mesh.triVAO[0] = create_array();
	new_mesh.triVAO[1] = tri_ind_buf_size_bytes / (sizeof(unsigned));
	new_mesh.triVBO = create_buffer(GL_ARRAY_BUFFER, tri_buf, tri_buf_size_bytes, GL_STATIC_DRAW);
	new_mesh.triIBO = create_buffer(GL_ELEMENT_ARRAY_BUFFER, tri_ind_buf, tri_ind_buf_size_bytes, GL_STATIC_DRAW);


	glBindVertexArray(new_mesh.triVAO[0]);

	glBindBuffer(GL_ARRAY_BUFFER, new_mesh.triVBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, new_mesh.triIBO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), nullptr);
	glEnableVertexAttribArray(0);

	glBindVertexArray(0); GL_CHECK_ERRORS;

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
	glGenBuffers(1, &BO);

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


char* Visualizer::read_file(const char* path) {
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

GLuint Visualizer::create_shader_from_file(GLenum shader_type, const char* shader_path) {
	return create_shader(shader_type, read_file(shader_path), false);
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
	imgui_io = ImGui::GetIO();

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
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();
		if (camera.update) camera.LookAtUpd();
		if (proj.update) proj.ProjUpd();

		{
			glUseProgram(triProg); GL_CHECK_ERRORS;

			glUniformMatrix4fv(glGetUniformLocation(triProg, "P"), 1, false, proj.P); GL_CHECK_ERRORS;
			glUniformMatrix4fv(glGetUniformLocation(triProg, "V"), 1, false, camera.V); GL_CHECK_ERRORS;
			for (unsigned i = 0u; i < num_meshes; ++i) {
				if (active_meshes[i]) {
					VisIvec2 tmpVAO{ meshes[instances[i].mesh].triVAO[0], meshes[instances[i].mesh].triVAO[1] };
					glBindVertexArray(tmpVAO[0]);
					glUniformMatrix4fv(glGetUniformLocation(triProg, "M"), 1, false, instances[i].M);
					glDrawElements(GL_TRIANGLES, tmpVAO[1], GL_UNSIGNED_INT, nullptr);
					glBindVertexArray(0);
				}
			}

			glUseProgram(0);
		}
		{
			ImGui_ImplOpenGL3_NewFrame();
			ImGui_ImplGlfw_NewFrame();
			ImGui::NewFrame();
			mouse.timedel = imgui_io.DeltaTime;

			ImGui::Begin("DST visualizer");
			ImGui::Text("Active meshes:   ");

			for (unsigned i = 0u; i < num_meshes; ++i) {
				if (ImGui::Selectable((std::string("   ") + std::to_string(i) + std::string("   ")).c_str(), active_meshes[i], 0, ImVec2(55, 30)))
					active_meshes[i] = !active_meshes[i];
			}

			ImGui::SameLine();
			ImGui::End();
			ImGui::Render();
			ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		}

		glfwSwapBuffers(window);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	}

	delete[] active_meshes;

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
}


void VisCamera::CamPosUpd() {
	CamPos.x = cos(yaw) * cos(pitch);
	CamPos.y = sin(pitch);
	CamPos.z = sin(yaw) * cos(pitch);

	float norm = sqrt(CamPos.x * CamPos.x + CamPos.y * CamPos.y + CamPos.z * CamPos.z);
	CamPos.x *= dist / norm;
	CamPos.y *= dist / norm;
	CamPos.z *= dist / norm;

	CamPos.x += LookAt.x;
	CamPos.y += LookAt.y;
	CamPos.z += LookAt.z;
}

void VisCamera::LookAtUpd() {
	ImVec4 f, s, u;

	update = false;
	CamPosUpd();

	f = ImVec4(LookAt.x - CamPos.x, LookAt.y - CamPos.y, LookAt.z - CamPos.z, 0.f);

	float norm = sqrt(f.x * f.x + f.y * f.y + f.z * f.z);
	if (norm > 10e-6 || norm < -10e-6)
		f = ImVec4(f.x / norm, f.y / norm, f.z / norm, 0.f);

	s = ImVec4(f.y * UpVect.z - f.z * UpVect.y, f.z * UpVect.x - f.x * UpVect.z, f.x * UpVect.y - f.y * UpVect.x, 0.f);

	norm = sqrt(s.x * s.x + s.y * s.y + s.z * s.z);
	if (norm > 10e-6 || norm < -10e-6)
		s = ImVec4(s.x / norm, s.y / norm, s.z / norm, 0.f);

	u = ImVec4(s.y * f.z - s.z * f.y, s.z * f.x - s.x * f.z, s.x * f.y - s.y * f.x, 0.f);

	float dots[3] = { s.x * CamPos.x + s.y * CamPos.y + s.z * CamPos.z,
					  u.x * CamPos.x + u.y * CamPos.y + u.z * CamPos.z,
					  f.x * CamPos.x + f.y * CamPos.y + f.z * CamPos.z };

	V[0] =       s.x;
	V[4] =       s.y;
	V[8] =       s.z;
	V[12] = -dots[0];

	V[1] =       u.x;
	V[5] =       u.y;
	V[9] =       u.z;
	V[13] = -dots[1];

	V[2] =      -f.x;
	V[6] =      -f.y;
	V[10] =     -f.z;
	V[14] =  dots[2];
	return;
}

void VisProj::ProjUpd() {
	update = false;

	float tanHalfFovy = tan(3.1415f / 180.f * fov_deg);
	P[0] = 1.f / (width / height * tanHalfFovy);
	P[5] = 1.f / tanHalfFovy;
	P[11] = -1.f;
	P[10] = -(far + near) / (far - near);
	P[14] = -(2.f * far * near) / (far - near);
}