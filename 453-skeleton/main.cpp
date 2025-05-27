#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Geometry.h"
#include "SceneNode.h"
#include "ShaderLoader.h"

// DEBUGGING PURPOSES
void printVectorOfPairs(const std::vector<std::pair<glm::vec3, glm::vec3>>& vec) {
	for (const auto& pair : vec) {
		std::cout << "Pair: " << std::endl;
		std::cout << "  First (vec3): (" << pair.first.x << ", " << pair.first.y << ", " << pair.first.z << ")" << std::endl;
		std::cout << "  Second (vec3): (" << pair.second.x << ", " << pair.second.y << ", " << pair.second.z << ")" << std::endl;
	}
}
void fillMappingGeometry(
	const std::vector<std::pair<glm::vec3, glm::vec3>>& mappings,
	CPU_Geometry& mappingLines
) {
	for (const auto& [contourPoint, closestPoint] : mappings) {
		// add both vertices (pair)
		mappingLines.verts.push_back(contourPoint);
		mappingLines.verts.push_back(closestPoint);

		// indices for the line
		int startIndex = mappingLines.verts.size() - 2;
		mappingLines.indices.push_back(startIndex);
		mappingLines.indices.push_back(startIndex + 1);
	}
}

GLuint vao, vboPos, vboColor, ebo;

void setupBuffers() {
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	glGenBuffers(1, &vboPos);
	glBindBuffer(GL_ARRAY_BUFFER, vboPos);
	glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);  
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	glEnableVertexAttribArray(0);

	glGenBuffers(1, &vboColor);
	glBindBuffer(GL_ARRAY_BUFFER, vboColor);
	glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);  
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	glEnableVertexAttribArray(1);

	glGenBuffers(1, &ebo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);

	glBindVertexArray(0);
}

void updateBuffers(const std::vector<glm::vec3>& verts,
	const std::vector<glm::vec3>& colors,
	const std::vector<unsigned int>& indices) {
	glBindVertexArray(vao);

	glBindBuffer(GL_ARRAY_BUFFER, vboPos);
	glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(glm::vec3), verts.data(), GL_DYNAMIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, vboColor);
	glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(glm::vec3), colors.data(), GL_DYNAMIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_DYNAMIC_DRAW);
}

void draw(GLenum primitive, GLsizei vertexCount, GLsizei indexCount) {
	glBindVertexArray(vao);
	glDrawArrays(GL_POINTS, 0, vertexCount);
	glDrawElements(primitive, indexCount, GL_UNSIGNED_INT, 0);
}

int main() {
	glfwInit();
	GLFWwindow* window = glfwCreateWindow(800, 800, "Branching Structure", NULL, NULL);
	glfwMakeContextCurrent(window);
	gladLoadGL();

	glEnable(GL_DEPTH_TEST);
	GLuint shader = ShaderLoader("D:/Program/C++/NewPhytologist2017/articulated-structure/articulated-structure/assets/shaders/test.vert", "D:/Program/C++/NewPhytologist2017/articulated-structure/articulated-structure/assets/shaders/test.frag").ID;

	// create and bind VAO and VBO
	setupBuffers();


	// branch initialization
	CPU_Geometry branchGeometry;
	CPU_Geometry branchUpdates;
	SceneNode* root = SceneNode::createBranch(0, 1, 45.0f, 1.0f, false);
	root->updateBranch(glm::mat4(1.0f), glm::mat4(1.0f), branchGeometry);
	// contour initialization
	CPU_Geometry contourGeometry;
	root->generateInitialContourControlPoints(root);
	//std::vector<glm::vec3> contour = SceneNode::bSplineCurve(0, root);
	std::vector<glm::vec3> contour = SceneNode::contourCatmullRom(root, 8);
	// branch-contour mapping
	std::vector<std::pair<SceneNode*, SceneNode*>> pairs;
	std::vector<std::pair<SceneNode*, SceneNode*>> branchPairs = root->getBranches(root, pairs);
	std::vector<ContourBinding> bindings = root->bindContourToBranches(contour, root, branchPairs);
	CPU_Geometry mappingLines;

	float lastTime = glfwGetTime();

	while (!glfwWindowShouldClose(window)) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// animation
		float currentTime = glfwGetTime();
		float deltaTime = (currentTime - lastTime)/10;
		lastTime = currentTime;
		root->animate(deltaTime);

		// need to clear geometry before calling update to draw the new positions
		branchGeometry.verts.clear();
		branchGeometry.cols.clear();
		branchGeometry.indices.clear();
		contourGeometry.verts.clear();
		contourGeometry.cols.clear();
		branchUpdates.verts.clear();
		branchUpdates.cols.clear();
		branchUpdates.indices.clear();

		// camera setup
		glm::mat4 view = glm::lookAt(glm::vec3(0, 0, 6), glm::vec3(0, 1, 0), glm::vec3(0, 1, 0));
		glm::mat4 proj = glm::perspective(glm::radians(45.0f), 800.f / 800.f, 0.1f, 100.f);
		glm::mat4 viewProj = proj * view;
		glUseProgram(shader);
		glUniformMatrix4fv(glGetUniformLocation(shader, "viewProj"), 1, GL_FALSE, glm::value_ptr(viewProj));

		// update branch position
		root->updateBranch(glm::mat4(1.0f), glm::mat4(1.0f), branchGeometry);
		root->interpolateBranch(bindings, branchUpdates);
		for (int i = 0; i < branchUpdates.verts.size(); ++i) branchUpdates.cols.push_back(glm::vec3(1.0f));

		// contour
		contour = root->animateContour(bindings);
		for (int i = 0; i < contour.size(); ++i) {
			contourGeometry.verts.push_back(contour[i]);
			contourGeometry.cols.push_back(glm::vec3(1.0f, 0.f, 0.f));
		}

		mappingLines.verts.clear();
		mappingLines.indices.clear();

		// UPDATE ONCE BRANCHES INTERPOLATE
		int i = 0;
		for (const auto& binding : bindings) {
			int startIdx = mappingLines.verts.size();
			mappingLines.verts.push_back(contour[i]);
			mappingLines.verts.push_back(binding.closestPoint);
			mappingLines.cols.push_back(glm::vec3(0.f, 0.f, 1.0f));
			mappingLines.cols.push_back(glm::vec3(0.f, 0.f, 1.0f));
			mappingLines.indices.push_back(startIdx);     // from contour
			mappingLines.indices.push_back(startIdx + 1); // to closest branch point
			++i;
		}

		// Branch
		updateBuffers(branchGeometry.verts, branchGeometry.cols, branchGeometry.indices);
		draw(GL_LINES, branchGeometry.verts.size(), branchGeometry.indices.size());

		// Interpolated branch
		updateBuffers(branchUpdates.verts, branchUpdates.cols, branchUpdates.indices);
		draw(GL_LINES, branchUpdates.verts.size(), branchUpdates.indices.size());

		// Contour
		updateBuffers(contourGeometry.verts, contourGeometry.cols, {});
		glPointSize(5);
		glDrawArrays(GL_POINTS, 0, contourGeometry.verts.size());
		glDrawArrays(GL_LINE_STRIP, 0, contourGeometry.verts.size());

		// Mapping (DEBUGGING PURPOSES)
		//updateBuffers(mappingLines.verts, mappingLines.cols, mappingLines.indices);
		//draw(GL_LINES, mappingLines.verts.size(), mappingLines.indices.size());

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwTerminate();
	return 0;
}

