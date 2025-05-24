#include "SceneNode.h"

#include "Geometry.h"
#include "Log.h"
#include "Panel.h"
#include "ShaderProgram.h"
#include "Window.h"

#include <imgui.h>
#include <memory>
#include <glm/glm.hpp>
#include <iostream>
#include <glm/gtc/type_ptr.hpp>
#include <cstdlib>

#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/spline.hpp>

// just a helper function to print the matrices for debugging purposes
void printMat4(const glm::mat4& mat) {
	for (int i = 0; i < 4; ++i) {
		std::cout << "| ";
		for (int j = 0; j < 4; ++j) {
			std::cout << mat[j][i] << "\t";
		}
		std::cout << "|\n";
	}
	std::cout << std::endl;
}

// Scene Graph Structure
SceneNode::SceneNode() : localTranslation(1.0f), localRotation(1.0f, 0.0f, 0.0f, 0.0f), localScaling(1.0f), parent(nullptr) {}

// function to add children to current node
void SceneNode::addChild(SceneNode* child) {
	child->parent = this;
	children.push_back(child);
}

// generating branches recursively using SceneNodes
// create scene graph (parent-child relationship)
// creating all the local transformations
SceneNode* SceneNode::createBranch(int depth, int maxDepth, float angle, float length, bool alternate) {
	// base case
	if (depth > maxDepth) return nullptr;

	// recursive case
	SceneNode* branch = new SceneNode();
	branch->localTranslation = glm::mat4(1.0f);
	branch->localRotation = glm::quat(1.0f, 0.f, 0.f, 0.f);
	branch->localScaling = glm::mat4(1.0f);

	float childLength = length * 0.5f;
	std::vector<float> angles = { angle };
	//std::vector<float> angles = { angle, 0.0f, -angle };

	// selecting angles based on alternating structure or symmetric structure
	std::vector<float> selectedAngles;
	if (!alternate) {
		selectedAngles = angles;
	}
	else {
		selectedAngles = (depth % 2 == 0)
			? std::vector<float>{angles[0], angles[1]}
		: std::vector<float>{ angles[1], angles[2] };
	}

	// children
	for (float a : selectedAngles) {
		SceneNode* child = createBranch(depth + 1, maxDepth, angle, childLength, alternate);
		if (!child) continue;

		glm::mat4 scale = glm::scale(glm::mat4(1.0f), glm::vec3(childLength));
		glm::quat rotQuat = glm::toQuat(glm::rotate(glm::mat4(1.0f), glm::radians(a), glm::vec3(0, 0, 1)));
		glm::mat4 translation = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 1.0f, 0.0f));

		child->localTranslation = translation;
		child->localRotation = rotQuat;
		child->localScaling = scale;
		//child->animationDirection = (a > 0) ? 1.0f : -1.0f; // left and right move opposite
		branch->addChild(child);
	}

	return branch;
}

// animation
void SceneNode::animate(float deltaTime) {
	// stop animation after certain time
	if (animationTime >= animationDuration) {
		return;
	}
	animationTime += deltaTime;

	animationAngle += deltaTime * 30.0f * animationDirection; 
	animateRotation = glm::toQuat(glm::rotate(glm::mat4(1.0f), glm::radians(animationAngle), glm::vec3(0, 0, 1)));

	animationScaling += deltaTime * animationScaling;
	animateScaling = glm::scale(glm::mat4(1.0f), glm::vec3(animationScaling));

	for (SceneNode* child : children) {
		child->animate(deltaTime);
	}
}

// compute global transformation for all nodes recursively
// to get the final position and draw
void SceneNode::updateBranches(const glm::mat4& parentAnimate, const glm::mat4& parentRest, CPU_Geometry& outGeometry) {
	// convert rotation quaternion back to matrix form
	glm::mat4 animateRotationMatrix = glm::toMat4(animateRotation);
	glm::mat4 localRotationMatrix = glm::toMat4(localRotation);
	// local to global animated matrix: A = T*V
	glm::mat4 global = parentAnimate * animateRotationMatrix * localRotationMatrix * localTranslation * localScaling *  animateScaling;
	globalTransformation = global;
	// global to local rest post matrix
	glm::mat4 local = glm::inverse(parentRest * localRotationMatrix * localTranslation * localScaling);
	restPose = local;
	// global position of node
	// for drawing purposes
	glm::vec3 rootPos = glm::vec3(global[3]);

	// parent index
	unsigned int currentIndex = outGeometry.verts.size();
	// parent geometry
	outGeometry.verts.push_back(rootPos);
	outGeometry.cols.push_back(glm::vec3(0.f, 1.f, 0.f));

	for (SceneNode* child : children) {
		// child index to draw line segment from parent to child pair
		unsigned int childIndex = outGeometry.verts.size();

		// parent and child pair (to draw line segment for each pair)
		outGeometry.indices.push_back(currentIndex);
		outGeometry.indices.push_back(childIndex);

		// recurse
		child->updateBranches(global, local, outGeometry);
	}
}

// destructor
void SceneNode::deleteSceneGraph(SceneNode* node) {
	if (!node) return;
	for (SceneNode* child : node->children) {
		deleteSceneGraph(child);
	}
	delete node;
}

// get all leaf nodes
void SceneNode::getLeafNodes(SceneNode* node, std::vector<SceneNode*>& leaves) {
	if (!node) return;
	if (node->children.empty()) {
		leaves.push_back(node);
	}
	else {
		for (SceneNode* child : node->children) {
			getLeafNodes(child, leaves);
		}
	}
}

// initial control points for contour using root and leaf nodes
std::vector<glm::vec3> SceneNode::generateInitialContourControlPoints(SceneNode* root) {
	std::vector<glm::vec3> controlPoints;

	// root
	glm::vec3 rootPos = glm::vec3(root->globalTransformation[3]);

	glm::vec3 leftOffset = rootPos - glm::vec3(0.5f, 0.0f, 0.0f);
	glm::vec3 rightOffset = rootPos + glm::vec3(0.5f, 0.0f, 0.0f);

	controlPoints.push_back(leftOffset); 

	// leaf
	std::vector<SceneNode*> leaves;
	getLeafNodes(root, leaves);

	for (SceneNode* leaf : leaves) {
		glm::vec3 leafPos = leaf->globalTransformation[3];
		glm::vec3 leafParentPos = leaf->parent->globalTransformation[3];
		glm::vec3 dir = glm::normalize(leafPos - leafParentPos);
		glm::vec3 offsetPos = leafPos + dir * 0.5f;

		controlPoints.push_back(offsetPos);
	}

	controlPoints.push_back(rightOffset);  

	return controlPoints;
}

// contour using b spline curve using chaikin's algorithm 
std::vector<glm::vec3> SceneNode::bSplineCurve(int subdivision, SceneNode* root) {
	std::vector<glm::vec3> curve = generateInitialContourControlPoints(root);
	// special case for 3 initial points
	if (curve.size() == 3) {
		// take the average of x values
		float newX = abs(curve.at(1).x - curve.at(0).x) / 2;
		std::vector <glm::vec3> newInitial;
		newInitial.push_back(curve.at(0));
		newInitial.push_back(glm::vec3(curve.at(1).x - newX, curve.at(1).y, curve.at(1).z));
		newInitial.push_back(glm::vec3(curve.at(1).x + newX, curve.at(1).y, curve.at(1).z));
		newInitial.push_back(curve.at(2));
		curve = newInitial;
		
	}
	// based on the number of iterations to repeat the chaikin's algorithm:
	for (int iter = 0; iter < subdivision; ++iter) {
		std::vector<glm::vec3> newCurve;

		// special mask the first control point
		if (!curve.empty()) {
			newCurve.push_back(curve.front());
		}
		newCurve.push_back(0.5f * curve[0] + 0.5f * curve[1]);

		// chaikin's algorithm  
		for (int i = 1; i < curve.size() - 2; ++i) {
			glm::vec3 p1 = 0.75f * curve[i] + 0.25f * curve[i + 1];
			glm::vec3 p2 = 0.25f * curve[i] + 0.75f * curve[i + 1];

			newCurve.push_back(p1);
			newCurve.push_back(p2);
		}

		// special mask for the last point
		newCurve.push_back(0.5f * curve[curve.size() - 2] + 0.5f * curve[curve.size() - 1]);
		if (!curve.empty()) {
			newCurve.push_back(curve.back());
		}

		// update the curve for the next iteration
		curve = newCurve;
	}
	return curve;
}

// contour using catmull rom spline 
std::vector<glm::vec3> SceneNode::contourCatmullRom(SceneNode* root, int points) {
	std::vector<glm::vec3> controlPoints = generateInitialContourControlPoints(root);
	// add first and last points
	std::vector<glm::vec3> paddedPoints;
	glm::vec3 first = controlPoints[0] + (controlPoints[0] - controlPoints[1]);
	glm::vec3 last = controlPoints.back() + (controlPoints.back() - controlPoints[controlPoints.size() - 1]);

	paddedPoints.push_back(first);
	paddedPoints.insert(paddedPoints.end(), controlPoints.begin(), controlPoints.end());
	paddedPoints.push_back(last);

	// generate contour points
	std::vector<glm::vec3> curvePoints;
	// minus 3 because we have i + 3
	for (size_t i = 0; i < paddedPoints.size() - 3; ++i) {
		for (int j = 0; j <= points; ++j) {
			float t = float(j) / points;
			glm::vec3 pt = glm::catmullRom(
				paddedPoints[i],
				paddedPoints[i + 1],
				paddedPoints[i + 2],
				paddedPoints[i + 3],
				t
			);
			curvePoints.push_back(pt);
		}
	}
	return curvePoints;
}

// finding the closest point from R (on contour) to the branch segment PQ 
glm::vec3 SceneNode::intersectionPoint(glm::vec3 P, glm::vec3 Q, glm::vec3 R) {
	float t = dot(Q - P, R - P) / dot(Q - P, Q - P);
	if (t <= 0) return P;
	else if (t >= 1) return Q;
	else {
		return P + t * (Q - P);
	}	
}

// extract parent child pair (endpoints of the branches)
std::vector<std::pair<SceneNode*, SceneNode*>> SceneNode::getBranches(SceneNode* node, std::vector<std::pair<SceneNode*, SceneNode*>>& segments) {
	for (SceneNode* child : node->children) {
		segments.push_back({ node, child });
		getBranches(child, segments);
	}
	return segments;
}

// binding for each contour point
std::vector<ContourBinding> SceneNode::bindContourToBranches(const std::vector<glm::vec3>& contourPoints, SceneNode* root, std::vector<std::pair<SceneNode*, SceneNode*>>& segments) {
	std::vector<ContourBinding> bindings;

	// for each contour point
	for (const glm::vec3& contourPoint : contourPoints) {
		float minDist = FLT_MAX;
		ContourBinding bestBinding;

		for (auto& [parent, child] : segments) {
			// position of endpoints of the branches
			glm::vec3 P = parent->globalTransformation[3]; 
			glm::vec3 Q = child->globalTransformation[3];    
			glm::vec3 closest = SceneNode::intersectionPoint(P, Q, contourPoint);
			float t = glm::dot(Q - P, contourPoint - P) / glm::dot(Q - P, Q - P);
			t = glm::clamp(t, 0.0f, 1.0f);

			float dist = glm::length(closest - contourPoint);
			if (dist < minDist) {
				minDist = dist;
				bestBinding = { parent, child, contourPoint, t, closest };
			}
		}
		bindings.push_back(bestBinding);
	}

	return bindings;
}

// interpolate the branch transformations and apply to contour points
std::vector<glm::vec3> SceneNode::animateContour(const std::vector<ContourBinding>& bindings) {
	std::vector<glm::vec3> animatedPoints;

	for (const auto& binding : bindings) {
		// NEED TO SEPARATE 
		glm::mat4 parentT = binding.parentNode->globalTransformation;
		glm::mat4 childT = binding.childNode->globalTransformation;

		// linearly interpolate matrix, then apply to the contour point
		// NOTE: blending matrices don't work well, so need to multiply by the inverse first then blend
		glm::mat4 animatedPosMat = binding.t * childT * binding.childNode->restPose + (1 - binding.t) * (parentT* binding.parentNode->restPose);

		// P' = A'A-1P
		animatedPoints.push_back(animatedPosMat * glm::vec4(binding.contourPoint, 1.0f));
		//glm::vec3 finalPos = glm::vec3(interpolatedMat * restPoseMat * glm::vec4(binding.contourPoint, 1.0f));
		//float elapsedTime = glm::clamp(animationTime / animationDuration, 0.0f, 1.0f);
		//// interpolate original and final pos using time elapsed for animation
		//glm::vec3 currentPos = glm::mix(binding.contourPoint, finalPos, elapsedTime);
		//animatedPoints.push_back(currentPos);
	}

	return animatedPoints;
}
