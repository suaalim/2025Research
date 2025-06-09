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
#include <map>
#include <numeric>
#include <cmath>
#include <tuple>
#include <vector>

// just a helper function to print the matrices for debugging purposes
void printMat4(const glm::mat4& mat) {
	for (int i = 0; i < 4; i++) {
		std::cout << "| ";
		for (int j = 0; j < 4; j++) {
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
//SceneNode* SceneNode::createBranch(int depth, int maxDepth, float angle, float length, bool alternate, std::vector<float> selectedAngles) {
//	// base case
//	if (depth > maxDepth) return nullptr;
//
//	// recursive case
//	SceneNode* branch = new SceneNode();
//	branch->localTranslation = glm::mat4(1.0f);
//	branch->localRotation = glm::quat(1.0f, 0.f, 0.f, 0.f);
//	branch->localScaling = glm::mat4(1.0f);
//
//	float childLength = length * 0.5f;
//
//	//// selecting angles based on alternating structure or symmetric structure
//	//std::vector<float> selectedAngles;
//	//if (!alternate) {
//	//	selectedAngles = angles;
//	//}
//	//else {
//	//	selectedAngles = (depth % 2 == 0)
//	//		? std::vector<float>{angles[0], angles[1]}
//	//	: std::vector<float>{ angles[1], angles[2] };
//	//}
//
//	// children
//	std::vector<float> angles = { 45.0f, 0.0f, -45.0f };
//	for (float a : selectedAngles) {
//		SceneNode* child = createBranch(depth + 1, maxDepth, angle, childLength, alternate, angles);
//		if (!child) continue;
//
//		// uniform scaling
//		glm::mat4 scale = glm::scale(glm::mat4(1.0f), glm::vec3(childLength));
//		// non-uniform scaling
//		//glm::mat4 scale = glm::scale(glm::mat4(1.0f), glm::vec3(1.0f, childLength, 1.0f));
//		glm::quat rotQuat = glm::toQuat(glm::rotate(glm::mat4(1.0f), glm::radians(a), glm::vec3(0, 0, 1)));
//		glm::mat4 translation = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 1.0f, 0.0f));
//
//		child->localTranslation = translation;
//		child->localRotation = rotQuat;
//		child->localScaling = scale;
//		//child->animationDirection = (a > 0) ? 1.0f : -1.0f; // left and right move opposite
//		branch->addChild(child);
//	}
//
//	return branch;
//}

SceneNode* SceneNode::createBranch(int depth, int maxDepth, float angle, float length, bool alternate, std::vector<float> selectedAngles) {
	// base case
	if (depth > maxDepth) return nullptr;

	// create this branch node
	SceneNode* branch = new SceneNode();
	branch->localTranslation = glm::mat4(1.0f);
	branch->localRotation = glm::quat(1.0f, 0.f, 0.f, 0.f);
	branch->localScaling = glm::mat4(1.0f);

	float childLength = length * 0.5f;

	// if we're at the root (depth 0), create an identity child
	if (depth == 0) {
		SceneNode* child = createBranch(depth + 1, maxDepth, angle, childLength, alternate, selectedAngles);
		if (child) {
			child->localTranslation = glm::mat4(1.0f);
			child->localRotation = glm::mat4(1.0f);
			child->localScaling = glm::mat4(1.0f);
			branch->addChild(child);
		}
	}
	else {
		// from depth 1 onward, apply the selected angles
		for (float a : selectedAngles) {
			SceneNode* child = createBranch(depth + 1, maxDepth, angle, childLength, alternate, selectedAngles);
			if (!child) continue;

			glm::mat4 scale = glm::scale(glm::mat4(1.0f), glm::vec3(childLength));
			glm::quat rotQuat = glm::toQuat(glm::rotate(glm::mat4(1.0f), glm::radians(a), glm::vec3(0, 0, 1)));
			glm::mat4 translation = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 1.0f, 0.0f));

			child->localTranslation = translation;
			child->localRotation = rotQuat;
			child->localScaling = scale;
			branch->addChild(child);
		}
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
	// uniform scaling
	animateScaling = glm::scale(glm::mat4(1.0f), glm::vec3(animationScaling));
	// non-uniform scaling
	//animateScaling = glm::scale(glm::mat4(1.0f), glm::vec3(1.0f, animationScaling, 1.0f));

	for (SceneNode* child : children) {
		child->animate(deltaTime);
	}
}

// compute global transformation for all nodes recursively
// to get the final position and draw
void SceneNode::updateBranch(const glm::mat4& parentTransform, const glm::mat4& parentRestInverse, const glm::mat4& parentRest, CPU_Geometry& outGeometry) {
	// convert rotation quaternion back to matrix form
	glm::mat4 animateRotationMatrix = glm::toMat4(animateRotation);
	glm::mat4 localRotationMatrix = glm::toMat4(localRotation);
	// local to global animated matrix: A = T*V
	globalTransformation = parentTransform * animateScaling * localRotationMatrix * localTranslation * localScaling;  // scale in the local coordinate
	// global to local rest post matrix
	// need to apply parentRest outside because if not it will be double inversed (inverse every call)
	restPoseInverse = glm::inverse(localRotationMatrix * localTranslation * localScaling) * parentRestInverse;
	// rest pose matrix
	restPose = parentRest * localRotationMatrix * localTranslation * localScaling;
	// global position of node
	// for drawing purposes
	glm::vec3 rootPos = glm::vec3(globalTransformation[3]);

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
		child->updateBranch(globalTransformation, restPoseInverse, restPose, outGeometry);
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

	glm::vec3 leftOffset = rootPos - glm::vec3(0.5f, 0.25f, 0.0f);
	glm::vec3 rightOffset = rootPos + glm::vec3(0.5f, -0.25f, 0.0f);

	controlPoints.push_back(leftOffset);

	// leaf
	std::vector<SceneNode*> leaves;
	getLeafNodes(root, leaves);

	for (SceneNode* leaf : leaves) {
		glm::vec3 leafPos = leaf->globalTransformation[3];
		glm::vec3 leafParentPos = leaf->parent->globalTransformation[3];
		glm::vec3 dir = glm::normalize(leafPos - leafParentPos);
		glm::vec3 offsetPos = leafPos + dir * 0.15f;
		controlPoints.push_back(offsetPos);
	}

	controlPoints.push_back(rightOffset);

	return controlPoints;
}

std::vector<glm::vec3> SceneNode::midPoints(std::vector<glm::vec3>& contourPoints) {
	std::vector<glm::vec3> copyOfContour;
	copyOfContour.push_back(contourPoints[0]);
	for (int i = 1; i < contourPoints.size() - 2; i++) {
		copyOfContour.push_back(contourPoints[i]);
		copyOfContour.push_back((contourPoints[i] + contourPoints[i + 1]) / 2.f);
	}
	copyOfContour.push_back(contourPoints[contourPoints.size() - 2]);
	copyOfContour.push_back(contourPoints[contourPoints.size() - 1]);
	return copyOfContour;
}

std::vector<std::vector<glm::vec3>> SceneNode::contourCatmullRomGrouped(std::vector<glm::vec3> controlPoints, int pointsPerSegment) {
	std::vector<std::vector<glm::vec3>> groupedContourPoints;

	std::vector<glm::vec3> paddedPoints;
	glm::vec3 first = controlPoints[0] + (controlPoints[0] - controlPoints[1]);
	glm::vec3 last = controlPoints.back() + (controlPoints.back() - controlPoints[controlPoints.size() - 2]);

	paddedPoints.push_back(first);
	paddedPoints.insert(paddedPoints.end(), controlPoints.begin(), controlPoints.end());
	paddedPoints.push_back(last);

	for (size_t i = 0; i < paddedPoints.size() - 3; i++) {
		std::vector<glm::vec3> segmentPoints;

		for (int j = 0; j <= pointsPerSegment; j++) {
			float t = float(j) / pointsPerSegment;
			glm::vec3 pt = glm::catmullRom(
				paddedPoints[i],
				paddedPoints[i + 1],
				paddedPoints[i + 2],
				paddedPoints[i + 3],
				t
			);
			segmentPoints.push_back(pt);
		}

		groupedContourPoints.push_back(segmentPoints);
	}

	return groupedContourPoints;
}

// contour points are grouped by segments
// now we need to bind
// get each group, map the first and last points to the tip, the middle point to the tip, and interpolate the rest
std::vector<ContourBinding> SceneNode::bindInterpolatedContourToBranches(const std::vector<std::vector<glm::vec3>>& contourPoints, SceneNode* root, std::vector<std::pair<SceneNode*, SceneNode*>>& segments) {
	std::vector<ContourBinding> bindings;
	glm::vec3 rootPos = root->globalTransformation[3];

	// segment
	for (int i = 0; i < contourPoints.size(); i++) {
		// individual points in the segment
		for (int j = 0; j < contourPoints[i].size(); j++) {
			float minDist = FLT_MAX;
			ContourBinding bestBinding;

			// different scenarios depending on the index of i
			for (auto& [parent, child] : segments) {
				glm::vec3 P = parent->globalTransformation[3];
				glm::vec3 Q = child->globalTransformation[3];

				glm::vec3 closest = SceneNode::intersectionPoint(P, Q, contourPoints[i][j]);
				float t = glm::clamp(glm::dot(Q - P, contourPoints[i][j] - P) / glm::dot(Q - P, Q - P), 0.0f, 1.0f);
				float dist = glm::length(closest - contourPoints[i][j]);

				if (dist < minDist) {
					minDist = dist;
					// instead of taking t directly, interpolate along the branch
					if (i % 2 == 1) t = 1 - (j / ((float)contourPoints[i].size() - 1));
					else t = (j / ((float)contourPoints[i].size() - 1));
					bestBinding = { parent, child, contourPoints[i][j], t, t * Q + (1 - t) * P, glm::inverse(t * child->globalTransformation + (1 - t) * parent->globalTransformation)};

				}
			}
			bindings.push_back(bestBinding);
		}
	}

	// first and last contour points are mapped to the root
	bindings[0].parentNode = std::get<0>(segments[1]);
	bindings[0].childNode = std::get<1>(segments[1]);
	bindings[0].t = 0.f;
	bindings[0].closestPoint = (bindings[0].t * bindings[0].childNode->globalTransformation[3] + (1 - bindings[0].t) * bindings[0].parentNode->globalTransformation[3]);
	bindings[0].previousAnimateInverse = glm::inverse(bindings[0].t * bindings[0].childNode->globalTransformation + (1 - bindings[0].t) * bindings[0].parentNode->globalTransformation);

	bindings[bindings.size() - 1].parentNode = std::get<0>(segments[1]);
	bindings[bindings.size() - 1].childNode = std::get<1>(segments[1]);
	bindings[bindings.size() - 1].t = 0.f;
	bindings[bindings.size() - 1].closestPoint = (bindings[bindings.size() - 1].t * bindings[bindings.size() - 1].childNode->globalTransformation[3] + (1 - bindings[bindings.size() - 1].t) * bindings[0].parentNode->globalTransformation[3]);
	bindings[bindings.size() - 1].previousAnimateInverse = glm::inverse(bindings[bindings.size() - 1].t * bindings[bindings.size() - 1].childNode->globalTransformation + (1 - bindings[bindings.size() - 1].t) * bindings[bindings.size() - 1].parentNode->globalTransformation);

	return bindings;
}

// contour using catmull rom spline 
std::vector<glm::vec3> SceneNode::contourCatmullRom(std::vector<glm::vec3> controlPoints, int points) {
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
	for (size_t i = 0; i < paddedPoints.size() - 3; i++) {
		for (int j = 0; j <= points; j++) {
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
void SceneNode::getBranches(SceneNode* node, std::vector<std::pair<SceneNode*, SceneNode*>>& segments) {
	for (SceneNode* child : node->children) {
		segments.push_back({ node, child });
		getBranches(child, segments);
	}
}

std::vector<ContourBinding> SceneNode::bindContourToBranches(const std::vector<glm::vec3>& contourPoints, SceneNode* root, std::vector<std::pair<SceneNode*, SceneNode*>>& segments) {
	std::vector<ContourBinding> bindings;
	glm::vec3 rootPos = root->globalTransformation[3];

	for (const glm::vec3& contourPoint : contourPoints) {
		float minDist = FLT_MAX;
		ContourBinding bestBinding;

		for (auto& [parent, child] : segments) {
			glm::vec3 P = parent->globalTransformation[3];
			glm::vec3 Q = child->globalTransformation[3];

			glm::vec3 closest = SceneNode::intersectionPoint(P, Q, contourPoint);
			float t = glm::dot(Q - P, contourPoint - P) / glm::dot(Q - P, Q - P);
			t = glm::clamp(t, 0.0f, 1.0f);
			float dist = glm::length(closest - contourPoint);

			if (dist < minDist) {
				minDist = dist;
				bestBinding = { parent, child, contourPoint, t, closest, glm::inverse(t * child->globalTransformation + (1 - t) * parent->globalTransformation) };

			}
		}

		bindings.push_back(bestBinding);
	}

	// first and last contour points are mapped to the root
	bindings[0].parentNode = std::get<0>(segments[1]);
	bindings[0].childNode = std::get<1>(segments[1]);
	bindings[0].t = 0.f;
	bindings[0].closestPoint = (bindings[0].t * bindings[0].childNode->globalTransformation[3] + (1 - bindings[0].t) * bindings[0].parentNode->globalTransformation[3]);
	bindings[0].previousAnimateInverse = glm::inverse(bindings[0].t * bindings[0].childNode->globalTransformation + (1 - bindings[0].t) * bindings[0].parentNode->globalTransformation);

	bindings[bindings.size() - 1].parentNode = std::get<0>(segments[1]);
	bindings[bindings.size() - 1].childNode = std::get<1>(segments[1]);
	bindings[bindings.size() - 1].t = 0.f;
	bindings[bindings.size() - 1].closestPoint = (bindings[bindings.size() - 1].t * bindings[bindings.size() - 1].childNode->globalTransformation[3] + (1 - bindings[bindings.size() - 1].t) * bindings[0].parentNode->globalTransformation[3]);
	bindings[bindings.size() - 1].previousAnimateInverse = glm::inverse(bindings[bindings.size() - 1].t * bindings[bindings.size() - 1].childNode->globalTransformation + (1 - bindings[bindings.size() - 1].t) * bindings[bindings.size() - 1].parentNode->globalTransformation);

	return bindings;
}

// get midpoint of the contourPoints and map to branching points
std::vector<ContourBinding> SceneNode::branchingPointMap(std::vector<ContourBinding>& contourPoints) {
	std::vector<ContourBinding> copyOfContour;
	copyOfContour.push_back(contourPoints[0]);
	for (int i = 1; i < contourPoints.size() - 2; i++) {
		// midpoint is binded to the parent branch of the original points
		copyOfContour.push_back(contourPoints[i]);
		ContourBinding newContour;
		newContour.parentNode = contourPoints[i].parentNode->parent;
		newContour.childNode = contourPoints[i].parentNode;
		newContour.contourPoint = (contourPoints[i].contourPoint + contourPoints[i + 1].contourPoint) / 2.f;
		newContour.t = 1.f;
		newContour.closestPoint = (newContour.t * newContour.childNode->globalTransformation[3] + (1 - newContour.t) * newContour.parentNode->globalTransformation[3]);
		newContour.previousAnimateInverse = glm::inverse((newContour.t * newContour.childNode->globalTransformation + (1 - newContour.t) * newContour.parentNode->globalTransformation));
		copyOfContour.push_back(newContour);
	}
	copyOfContour.push_back(contourPoints[contourPoints.size() - 2]);
	copyOfContour.push_back(contourPoints[contourPoints.size() - 1]);
	return copyOfContour;
}

// bind the interpolated points to the branch using blending
std::vector<ContourBinding> SceneNode::interpolateBetweenContour(std::vector<ContourBinding>& contourPoints) {
	std::vector<ContourBinding> copyOfContour;
	copyOfContour.push_back(contourPoints[0]);
	for (int i = 0; i < contourPoints.size() - 1; i++) {
		copyOfContour.push_back(contourPoints[i]);
		std::vector<glm::vec3> points = contourCatmullRom({ contourPoints[i].contourPoint, contourPoints[i + 1].contourPoint }, 6);
		for (int j = 1; j < points.size() - 1; j++) {
			ContourBinding newContour;
			if (i % 2 == 1) {
				newContour.parentNode = contourPoints[i].parentNode;
				newContour.childNode = contourPoints[i].childNode;
				newContour.t = 1 - (j / (float)points.size());
			}
			else {
				newContour.parentNode = contourPoints[i + 1].parentNode;
				newContour.childNode = contourPoints[i + 1].childNode;
				newContour.t = j / (float)points.size();
			}
			newContour.contourPoint = points[j];
			newContour.closestPoint = (newContour.t * newContour.childNode->globalTransformation[3] + (1 - newContour.t) * newContour.parentNode->globalTransformation[3]);
			newContour.previousAnimateInverse = glm::inverse((newContour.t * newContour.childNode->globalTransformation + (1 - newContour.t) * newContour.parentNode->globalTransformation));
			copyOfContour.push_back(newContour);
		}
	}
	copyOfContour.push_back(contourPoints[contourPoints.size() - 1]);
	return copyOfContour;
}

// interpolate branches 
void SceneNode::interpolateBranchTransforms(std::vector<std::pair<SceneNode*, SceneNode*>>& pair, std::vector<CPU_Geometry>& outGeometry) {
	for (auto& [parent, child] : pair) {
		glm::mat4 T1 = parent->globalTransformation;
		glm::mat4 T2 = child->globalTransformation;
		CPU_Geometry geom;

		for (int i = 0; i <= 7; i++) {
			float t = static_cast<float>(i) / 7.0f;

			glm::mat4 animatedMat =
				t * child->globalTransformation * child->restPoseInverse +
				(1 - t) * parent->globalTransformation * parent->restPoseInverse;

			glm::vec3 pos = t * child->restPose[3] + (1 - t) * parent->restPose[3];

			geom.verts.push_back(glm::vec3(animatedMat * glm::vec4(pos, 1.0f)));
			geom.cols.push_back(glm::vec3(1.0f));
		}

		outGeometry.push_back(geom);
	}
}

// check distance bewteen each contour points
// if distance is longer than a threshold, add a new point in bewteen those points
std::vector<glm::vec3> SceneNode::distanceBetweenContourPoints(std::vector<glm::vec3> contourPoints) {
	std::vector<glm::vec3> newContourPoints;
	// arbitrary threshold
	float threshold = 0.3f;
	for (int i = 0; i < contourPoints.size() - 1; i++) {
		float distance = glm::length(contourPoints[i + 1] - contourPoints[i]);
		if (distance >= threshold) {
			// add a point in between the two original points
			newContourPoints.push_back(contourPoints[i]);
			newContourPoints.push_back(glm::mix(contourPoints[i], contourPoints[i + 1], 0.5f));
			contourChanged = true;
		}
		else {
			newContourPoints.push_back(contourPoints[i]);
		}
	}
	newContourPoints.push_back(contourPoints[contourPoints.size() - 1]);
	return newContourPoints;
}

// add a new point to the curve
std::vector<ContourBinding> SceneNode::addContourPoints(std::vector<ContourBinding>& bindings) {
	std::vector<ContourBinding> newBindingSet;
	// arbitrary threshold
	float threshold = 0.2f;
	for (int i = 0; i < bindings.size() - 1; i++) {
		float distance = glm::length(bindings[i + 1].contourPoint - bindings[i].contourPoint);
		if (distance >= threshold) {
			// add a point in between the two original points
			glm::vec3 newPoint = glm::mix(bindings[i].contourPoint, bindings[i + 1].contourPoint, 0.5f);

			// which branch to bind to?
			glm::vec3 firstNeighborParent = bindings[i].parentNode->globalTransformation[3];
			glm::vec3 firstNeighborChild = bindings[i].childNode->globalTransformation[3];
			glm::vec3 proj1 = intersectionPoint(firstNeighborParent, firstNeighborChild, newPoint);
			float dist1 = glm::length(newPoint - proj1);

			glm::vec3 secondNeighborParent = bindings[i + 1].parentNode->globalTransformation[3];
			glm::vec3 secondNeighborChild = bindings[i + 1].childNode->globalTransformation[3];
			glm::vec3 proj2 = intersectionPoint(secondNeighborParent, secondNeighborChild, newPoint);
			float dist2 = glm::length(newPoint - proj2);

			if (dist1 < dist2) {
				float t1 = (bindings[i].t + bindings[i + 1].t) / 2.f;
				glm::mat4 previousiInverseAnimationMat = glm::inverse(t1 * bindings[i].childNode->globalTransformation + (1 - t1) * bindings[i].parentNode->globalTransformation);
				newBindingSet.push_back(bindings[i]);
				newBindingSet.push_back({ bindings[i].parentNode, bindings[i].childNode, newPoint, t1, glm::mix(bindings[i].closestPoint, bindings[i + 1].closestPoint, 0.5f), previousiInverseAnimationMat });
			}
			else {
				float t2 = (bindings[i].t + bindings[i + 1].t) / 2.f;
				glm::mat4 previousiInverseAnimationMat = glm::inverse(t2 * bindings[i + 1].childNode->globalTransformation + (1 - t2) * bindings[i + 1].parentNode->globalTransformation);
				newBindingSet.push_back(bindings[i]);
				newBindingSet.push_back({ bindings[i + 1].parentNode, bindings[i + 1].childNode, newPoint, t2, glm::mix(bindings[i].closestPoint, bindings[i + 1].closestPoint, 0.5f), previousiInverseAnimationMat });
			}
		}
		else {
			newBindingSet.push_back(bindings[i]);
		}
	}
	newBindingSet.push_back(bindings[bindings.size() - 1]);

	return newBindingSet;
}

// relative transformation between frame (global frame)
void SceneNode::animationPerFrame(std::vector<ContourBinding>& bindings) {
	// need the previous frame's animation matrix and current frame's animation matrix
	// update the contour point every frame (in the ContourBinding) so that you just apply the matrix to the contour point
	// this is in the "global" frame
	for (auto& binding : bindings) {
		glm::mat4 animatedPosMat = binding.t * binding.childNode->globalTransformation + (1 - binding.t) * (binding.parentNode->globalTransformation);
		binding.contourPoint = animatedPosMat * binding.previousAnimateInverse * glm::vec4(binding.contourPoint, 1.0f);
		binding.previousAnimateInverse = glm::inverse(binding.t * binding.childNode->globalTransformation + (1 - binding.t) * binding.parentNode->globalTransformation);
	}
}

// label the branches by indices
void SceneNode::labelBranches(SceneNode* node, std::vector<std::tuple<SceneNode*, SceneNode*, int>>& segments, int& i) {
	for (SceneNode* child : node->children) {
		segments.push_back({ node, child, i });
		i++;
		labelBranches(child, segments, i);
	}
}

void SceneNode::handleMouseClick(double xpos, double ypos, int screenWidth, int screenHeight, glm::mat4 view, glm::mat4 projection, std::vector<glm::vec3> contourPoints, CPU_Geometry geom) {
	glm::vec2 clickPos;
	clickPos = glm::vec2(xpos, ypos);
	clickPos += glm::vec2(0.5f, 0.5f);
	clickPos /= (glm::vec2(screenWidth, screenHeight));
	clickPos = glm::vec2(clickPos.x, 1.0f - clickPos.y);
	clickPos *= 2.0f;
	clickPos -= glm::vec2(1.0f, 1.0f);

	contourPoints.push_back(glm::vec3(clickPos, 0.0f));
	std::cout << glm::to_string(clickPos) << std::endl;

}
