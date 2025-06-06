#pragma once

//------------------------------------------------------------------------------
// This file contains simple classes for storing geomtery on the CPU and the GPU
// Later assignments will require you to expand these classes or create your own
// similar classes with the needed functionality
//------------------------------------------------------------------------------

#include "VertexArray.h"
#include "VertexBuffer.h"

#include <glad/glad.h>
#include <glm/glm.hpp>

#include <vector>


// List of vertices and texture coordinates using std::vector and glm::vec3
struct CPU_Geometry {
	// for index buffer
	std::vector<unsigned int> indices;
	std::vector<glm::vec3> verts;
	std::vector<glm::vec3> cols;
};


// VAO and two VBOs for storing vertices and texture coordinates, respectively
class GPU_Geometry {
public:
	GPU_Geometry();
	// Public interface
	void bind() {
		vao.bind();
	}
	void setVerts(const std::vector<glm::vec3>& verts);
	void setCols(const std::vector<glm::vec3>& cols);
protected:
	// note: due to how OpenGL works, vao needs to be
// defined and initialized before the vertex buffers
	VertexArray vao;

	VertexBuffer vertBuffer;
	VertexBuffer colorsBuffer;
private:

};
