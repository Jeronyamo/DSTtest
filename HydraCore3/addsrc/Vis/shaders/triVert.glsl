#version 420 core

layout(location = 0) in vec3 position;

uniform mat4 M, V, P;


void main() {
	gl_Position = P * V * M * vec4(position, 1.0);
}