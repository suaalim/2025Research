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

int main() {
	glfwInit();
	GLFWwindow* window = glfwCreateWindow(800, 800, "Branching Structure", NULL, NULL);
	glfwMakeContextCurrent(window);
	gladLoadGL();

	glEnable(GL_DEPTH_TEST);
	GLuint shader = ShaderLoader("D:/Program/C++/NewPhytologist2017/articulated-structure/articulated-structure/assets/shaders/test.vert", "D:/Program/C++/NewPhytologist2017/articulated-structure/articulated-structure/assets/shaders/test.frag").ID;

	// create and bind VAO and VBO
	GLuint vaoBranch, vaoContour, vboBranch, vboContour, vboBranchColor, vboContourColor, vboMappingColor;
	glGenVertexArrays(1, &vaoBranch);
	glGenVertexArrays(1, &vaoContour);
	glGenBuffers(1, &vboBranch);
	glGenBuffers(1, &vboContour);
	glGenBuffers(1, &vboBranchColor);
	glGenBuffers(1, &vboContourColor);
	glGenBuffers(1, &vboMappingColor);

	// branch initialization
	CPU_Geometry branchGeometry;
	SceneNode* root = SceneNode::createBranch(0, 1, 45.0f, 1.0f, false);
	root->updateBranches(glm::mat4(1.0f), glm::mat4(1.0f), branchGeometry);
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

		// camera setup
		glm::mat4 view = glm::lookAt(glm::vec3(0, 0, 6), glm::vec3(0, 1, 0), glm::vec3(0, 1, 0));
		glm::mat4 proj = glm::perspective(glm::radians(45.0f), 800.f / 800.f, 0.1f, 100.f);
		glm::mat4 viewProj = proj * view;
		glUseProgram(shader);
		glUniformMatrix4fv(glGetUniformLocation(shader, "viewProj"), 1, GL_FALSE, glm::value_ptr(viewProj));

		// update branch position
		root->updateBranches(glm::mat4(1.0f), glm::mat4(1.0f), branchGeometry);

		// contour
		contour = root->animateContourPoints(bindings);
		for (int i = 0; i < contour.size(); ++i) {
			contourGeometry.verts.push_back(contour[i]);
			contourGeometry.cols.push_back(glm::vec3(1.0f, 0.f, 0.f));
		}

		mappingLines.verts.clear();
		mappingLines.indices.clear();

		for (const auto& binding : bindings) {
			int startIdx = mappingLines.verts.size();
			mappingLines.verts.push_back(binding.contourPoint);
			mappingLines.verts.push_back(binding.closestPoint);
			mappingLines.cols.push_back(glm::vec3(0.f, 0.f, 1.0f));
			mappingLines.cols.push_back(glm::vec3(0.f, 0.f, 1.0f));
			mappingLines.indices.push_back(startIdx);     // from contour
			mappingLines.indices.push_back(startIdx + 1); // to closest branch point
		}

		// draw branch
		glBindVertexArray(vaoBranch);
		glBindBuffer(GL_ARRAY_BUFFER, vboBranch);
		glBufferData(GL_ARRAY_BUFFER, branchGeometry.verts.size() * sizeof(glm::vec3), branchGeometry.verts.data(), GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vboBranchColor);
		glBufferData(GL_ARRAY_BUFFER, branchGeometry.cols.size() * sizeof(glm::vec3), branchGeometry.cols.data(), GL_DYNAMIC_DRAW);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		glEnableVertexAttribArray(1);
		GLuint elementbuffer;
		glGenBuffers(1, &elementbuffer);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, branchGeometry.indices.size() * sizeof(unsigned int), branchGeometry.indices.data(), GL_DYNAMIC_DRAW);
		glPointSize(5);
		glDrawArrays(GL_POINTS, 0, branchGeometry.verts.size());
		glDrawElements(GL_LINES, branchGeometry.indices.size(), GL_UNSIGNED_INT, 0);
		// draw contour
		glBindVertexArray(vaoContour);
		glBindBuffer(GL_ARRAY_BUFFER, vboContour);
		glBufferData(GL_ARRAY_BUFFER, contourGeometry.verts.size() * sizeof(glm::vec3), contourGeometry.verts.data(), GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vboContourColor);
		glBufferData(GL_ARRAY_BUFFER, contourGeometry.cols.size() * sizeof(glm::vec3), contourGeometry.cols.data(), GL_DYNAMIC_DRAW);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		glEnableVertexAttribArray(1);
		glDrawArrays(GL_POINTS, 0, contourGeometry.verts.size());
		glDrawArrays(GL_LINE_STRIP, 0, contourGeometry.verts.size());
		// draw mapping (DEBUGGING PURPOSES)
		GLuint vaoMap, vboMap, eboMap;
		glGenVertexArrays(1, &vaoMap);
		glGenBuffers(1, &vboMap);
		glGenBuffers(1, &eboMap);
		glBindVertexArray(vaoMap);
		glBindBuffer(GL_ARRAY_BUFFER, vboMap);
		glBufferData(GL_ARRAY_BUFFER, mappingLines.verts.size() * sizeof(glm::vec3), mappingLines.verts.data(), GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eboMap);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, mappingLines.indices.size() * sizeof(unsigned int), mappingLines.indices.data(), GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, vboMappingColor);
		glBufferData(GL_ARRAY_BUFFER, mappingLines.cols.size() * sizeof(glm::vec3), mappingLines.cols.data(), GL_DYNAMIC_DRAW);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		glEnableVertexAttribArray(1);
		//glDrawElements(GL_LINES, mappingLines.indices.size(), GL_UNSIGNED_INT, 0);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwTerminate();
	return 0;
}

