#version 420 core

layout(location = 0) in vec4 position;

uniform mat4 M, V, P;


void main() {
	gl_Position = P * V * M * vec4(position.x, position.y, position.z, 1.0);
	gl_Position = gl_Position / gl_Position.w;
	gl_Position.w = position.w;
}