#pragma once

#include <memory>

#include "Geometry.h"
#include "Panel.h"
#include "ShaderProgram.h"
#include "Window.h"
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

// define the class before the struct since it uses the class
class SceneNode;

struct ContourBinding {
	SceneNode* parentNode;   
	SceneNode* childNode;
	glm::vec3 contourPoint;
	float t;                 
	glm::vec3 closestPoint;
};

// SceneNode for Scene Graph
class SceneNode {
public:
	SceneNode();
	void addChild(SceneNode* child);
	static SceneNode* createBranch(int depth, int maxDepth, float angle, float length, bool alternate);
	void updateBranches(const glm::mat4& parentTransform, const glm::mat4& parentRest, CPU_Geometry& outGeometry);
	void animate(float deltaTime);
	void deleteSceneGraph(SceneNode* node);
	static glm::vec3 intersectionPoint(glm::vec3 P, glm::vec3 Q, glm::vec3 R);
	std::vector<std::pair<SceneNode*, SceneNode*>> getBranches(SceneNode* node, std::vector<std::pair<SceneNode*, SceneNode*>>& segments);
	std::vector<ContourBinding> bindContourToBranches(
		const std::vector<glm::vec3>& contourPoints,
		SceneNode* root,
		std::vector<std::pair<SceneNode*, SceneNode*>>&
	);
	std::vector<glm::vec3> animateContourPoints(const std::vector<ContourBinding>& bindings);
	static void getLeafNodes(SceneNode* node, std::vector<SceneNode*>& leaves);
	static std::vector<glm::vec3> generateInitialContourControlPoints(SceneNode* root);
	static std::vector<glm::vec3> bSplineCurve(int iterations, SceneNode* root);
	static std::vector<glm::vec3> contourCatmullRom(SceneNode* root, int points);
	// global transformation A = T*V
	glm::mat4 globalTransformation = glm::mat4(1.0f);
	// global to local transformation for rest pose 
	glm::mat4 restPose;

private:
	SceneNode* parent;
	std::vector<SceneNode*> children;
	// T trasformation (rest pose)
	glm::mat4 localTranslation;
	glm::quat localRotation;
	glm::mat4 localScaling;
	// V transformation (animation)
	glm::mat4 animateTranslation = glm::mat4(1.0f);
	glm::quat animateRotation = glm::quat(1.0f, 0.0f, 0.0f, 0.0f);
	glm::mat4 animateScaling = glm::mat4(1.0f);
	// interpolated matrix for the contour
	glm::mat4 interpolatedTranslation = glm::mat4(1.0f);
	glm::quat interpolatedRotation = glm::quat(1.0f, 0.0f, 0.0f, 0.0f);
	glm::mat4 interpolatedScaling = glm::mat4(1.0f);
	// animation variables
	float deltatime = 0.0f;
	float animationDirection = 1.0f; // control how left and right branches move differently (+angle, -angle)
	float animationAngle = 0.0f;
	float animationScaling = 1.0f;
	float animationTime = 0.0f;
	float animationDuration = 0.3f; // how long the animation lasts

};
