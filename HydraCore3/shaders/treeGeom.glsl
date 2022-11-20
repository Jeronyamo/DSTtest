#version 420 core

layout(lines) in;
layout(triangle_strip, max_vertices = 8) out;

//out vec4 in_color;

uniform mat4 M, V, P;



void build_rectangle(float plane, vec2 axis1_minmax, vec2 axis2_minmax, uvec3 axes) {
	for (uint i = 0u; i < 4; ++i) {
		vec4 temp = vec4(0.f, 0.f, 0.f, 1.f);

		temp[axes.x] = plane;
		temp[axes.y] = axis1_minmax[i >> 1];
		temp[axes.z] = axis2_minmax[i & 1];

		gl_Position = P * V * M * temp;
		EmitVertex();
	}

	EndPrimitive();
}


void main() {
	int flag1 = int(gl_in[0].gl_Position.w), flag2 = int(gl_in[1].gl_Position.w);
	float plane = gl_in[0].gl_Position[0];
	uvec3 axes = uvec3(0, 1, 2);
	vec2 axis1_minmax = vec2(gl_in[0].gl_Position[axes.y], gl_in[1].gl_Position[axes.y]);
	vec2 axis2_minmax = vec2(gl_in[0].gl_Position[axes.z], gl_in[1].gl_Position[axes.z]);

	if (flag1 == -1) {
		build_rectangle(plane, axis1_minmax, axis2_minmax, axes);

		//in_color = vec4(0.f, 0.f, 0.f, 0.f);
	}
	if (flag1 == 0 || flag1 == 1) {
		plane = gl_in[1 - flag1].gl_Position[flag2];
		axes = uvec3(flag2, (flag2 + 1) % 3, (flag2 + 2) % 3);
		axis1_minmax = vec2(gl_in[0].gl_Position[axes.y], gl_in[1].gl_Position[axes.y]);
		axis2_minmax = vec2(gl_in[0].gl_Position[axes.z], gl_in[1].gl_Position[axes.z]);

		build_rectangle(plane, axis1_minmax, axis2_minmax, axes);

		vec4 colors[2] = { vec4(0.f, 0.8f, 0.8f, 1.f), vec4(0.f, 0.f, 0.8f, 1.f) };
		//in_color = colors[flag1];
	}
	if (flag1 == 2) {
		axes = uvec3(flag2, (flag2 + 1) % 3, (flag2 + 2) % 3);
		axis1_minmax = vec2(gl_in[0].gl_Position[axes.y], gl_in[1].gl_Position[axes.y]);
		axis2_minmax = vec2(gl_in[0].gl_Position[axes.z], gl_in[1].gl_Position[axes.z]);

		build_rectangle(gl_in[0].gl_Position[flag2], axis1_minmax, axis2_minmax, axes);
		build_rectangle(gl_in[1].gl_Position[flag2], axis1_minmax, axis2_minmax, axes);

		//in_color = vec4(0.8f, 0.f, 0.f, 1.f);
	}
	if (flag1 == 3) {
		uvec2 plane_ind = uvec2(flag2 >> 4, (flag2 >> 2) & 3);
		uvec2 is_max = uvec2((flag2 >> 1) & 1, flag2 & 1);

		for (uint i = 0u; i < 2u; ++i) {
			plane = gl_in[is_max[i]].gl_Position[plane_ind[i]];
			axes = uvec3(plane_ind[i], (plane_ind[i] + 1) % 3, (plane_ind[i] + 2) % 3);
			axis1_minmax = vec2(gl_in[0].gl_Position[axes.y], gl_in[1].gl_Position[axes.y]);
			axis2_minmax = vec2(gl_in[0].gl_Position[axes.z], gl_in[1].gl_Position[axes.z]);

			build_rectangle(plane, axis1_minmax, axis2_minmax, axes);
		}

		//in_color = vec4(0.8f, 0.3f, 0.f, 1.f);
	}
}

/*
	* flag1 =-1 -- one-noded tree, don't know what to do
	* flag1 = 0 -- split node first
			flag2 = plane (max)
			in_color = vec4(0.f, 0.8f, 0.8f, 1.f);
		-> 1 plane

	* flag1 = 1 -- split node second
			flag2 = plane (min)
			in_color = vec4(0.f, 0.f, 0.8f, 1.f);
		-> 1 plane

	* flag1 = 2 -- single-carve node
			flag2 = plane (min/max)
			in_color = vec4(0.8f, 0.f, 0.f, 1.f);
		-> 2 planes

	* flag1 = 3 -- double-carve node
			(flag2 >> 4) = plane 1
			(flag2 >> 2) & 3 = plane 2
			bool(flag2 & 2) = plane 1 min/max
			bool(flag2 & 1) = plane 2 min/max
			in_color = vec4(0.8f, 0.3f, 0.f, 1.f);
		-> 2 planes
*/