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
		glm::vec3 offsetPos = leafPos + dir * 0.15f;
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
		for (int i = 1; i < curve.size() - 2; i++) {
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

		glm::vec3 rootToContour = glm::normalize(contourPoint - rootPos);

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

	return bindings;
}

// interpolate the branch transformations and apply to contour points
std::vector<glm::vec3> SceneNode::animateContour(std::vector<ContourBinding>& bindings) {
	std::vector<glm::vec3> animatedPoints;

	for (auto& binding : bindings) {
		// linearly interpolate matrix, then apply to the contour point
		// NOTE: blending matrices don't work well, so need to multiply by the inverse first then blend
		// P' = A'A-1P but before we interpolate

		// my understanding: restPoseInverse moves the contour point wrt to the local coordinate frame that the root node is in
		// then globalTransformation takes the contour point, all the way to where the node that it is binded to is, then apply the same transformation as the binded node
		glm::mat4 animatedPosMat = binding.t * binding.childNode->globalTransformation * binding.childNode->restPoseInverse + (1 - binding.t) * (binding.parentNode->globalTransformation * binding.parentNode->restPoseInverse);
		animatedPoints.push_back(animatedPosMat * glm::vec4(binding.contourPoint, 1.0f));
	}

	return animatedPoints;
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

// inverse transform the deformed contour 
void SceneNode::inverseTransform(std::vector<ContourBinding>& bindings) {
	for (auto& binding : bindings) {
		glm::mat4 transformInverseMat = binding.t * glm::inverse(binding.childNode->globalTransformation * binding.childNode->restPoseInverse) + (1 - binding.t) * glm::inverse(binding.parentNode->globalTransformation * binding.parentNode->restPoseInverse);
		binding.contourPoint = glm::vec3(transformInverseMat * glm::vec4(binding.contourPoint, 1.0f));
	}
}

// add a new point to the curve
// when adding, get its neighbors and the branch that it belongs to
// calculate the closest point using those branches
// then you can get the animation matrix from the deformed branch
// will need to update the position of the contour point each frame

// change so that instead of copying the vector everytime, just insert to the right index
// need a way to initialize the weight vector and previous animation matrix list
// for the weights, take the average with its neighbors
// for animation inverse matrices, isert the inverse for this new added point to the previous list -> keep track of the index of the added point
std::vector<ContourBinding> SceneNode::addContourPoints(std::vector<ContourBinding>& bindings) {
	std::vector<ContourBinding> newBindingSet;
	// arbitrary threshold
	float threshold = 0.2f;
	//for (int i = 1; i < bindings.size() - 2; i++) {
	for (int i = 0; i < bindings.size() - 1; i++) {
		float distance = glm::length(bindings[i + 1].contourPoint - bindings[i].contourPoint);
		if (distance >= threshold) {
			// add a point in between the two original points
			glm::vec3 newPoint = glm::mix(bindings[i].contourPoint, bindings[i + 1].contourPoint, 0.5f);

			//// catmullrom interpolation
			//glm::vec3 p0 = bindings[i - 1].contourPoint;
			//glm::vec3 p1 = bindings[i].contourPoint;
			//glm::vec3 p2 = bindings[i + 1].contourPoint;
			//glm::vec3 p3 = bindings[i + 2].contourPoint;
			//float t = 0.5f; 
			//glm::vec3 newPoint = glm::catmullRom(p0, p1, p2, p3, t);

			glm::vec3 firstNeighborParent = bindings[i].parentNode->globalTransformation[3];
			glm::vec3 firstNeighborChild = bindings[i].childNode->globalTransformation[3];
			glm::vec3 proj1 = intersectionPoint(firstNeighborParent, firstNeighborChild, newPoint);
			float dist1 = glm::length(newPoint - proj1);

			glm::vec3 secondNeighborParent = bindings[i + 1].parentNode->globalTransformation[3];
			glm::vec3 secondNeighborChild = bindings[i + 1].childNode->globalTransformation[3];
			glm::vec3 proj2 = intersectionPoint(secondNeighborParent, secondNeighborChild, newPoint);
			float dist2 = glm::length(newPoint - proj2);


			if (dist1 < dist2) {
				float t1 = glm::dot(firstNeighborChild - firstNeighborParent, newPoint - firstNeighborParent) / glm::dot(firstNeighborChild - firstNeighborParent, firstNeighborChild - firstNeighborParent);
				//t1 = glm::clamp(t1, 0.0f, 1.0f);
				t1 = 1.f;
				glm::mat4 previousiInverseAnimationMat = glm::inverse(t1 * bindings[i].childNode->globalTransformation + (1 - t1) * bindings[i].parentNode->globalTransformation);
				newBindingSet.push_back(bindings[i]);
				newBindingSet.push_back({ bindings[i].parentNode, bindings[i].childNode, newPoint, t1, proj1, previousiInverseAnimationMat, bindings[i].weights, bindings[i].previousAnimateInverseMat });
			}
			else {
				float t2 = glm::dot(secondNeighborChild - secondNeighborParent, newPoint - secondNeighborParent) / glm::dot(secondNeighborChild - secondNeighborParent, secondNeighborChild - secondNeighborParent);
				//t2 = glm::clamp(t2, 0.0f, 1.0f);
				t2 = 1.f;
				glm::mat4 previousiInverseAnimationMat = glm::inverse(t2 * bindings[i + 1].childNode->globalTransformation + (1 - t2) * bindings[i + 1].parentNode->globalTransformation);
				newBindingSet.push_back(bindings[i]);
				newBindingSet.push_back({ bindings[i + 1].parentNode, bindings[i + 1].childNode, newPoint, t2, proj2, previousiInverseAnimationMat, bindings[i + 1].weights, bindings[i + 1].previousAnimateInverseMat });
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
		binding.previousAnimateInverse = glm::inverse(binding.t * binding.childNode->globalTransformation + (1 - binding.t) *binding.parentNode->globalTransformation);
	}
}

// MULTIPLE BRANCHES

// label the branches by indices
void SceneNode::labelBranches(SceneNode* node, std::vector<std::tuple<SceneNode*, SceneNode*, int>>& segments, int& i) {
	for (SceneNode* child : node->children) {
		segments.push_back({ node, child, i });
		i++;
		labelBranches(child, segments, i);
	}
}

// initialize a vector of weights for blending with multiple branches
std::vector<ContourBinding> SceneNode::bindContourToMultipleBranches(const std::vector<glm::vec3>& contourPoints, SceneNode* root, std::vector<std::tuple<SceneNode*, SceneNode*, int>>& segments) {
	std::vector<ContourBinding> bindings;
	std::vector<float> weights(segments.size());
	std::vector<glm::mat4> previousAnimateInverseVec(segments.size(), glm::mat4(1.0f));
	glm::vec3 rootPos = root->globalTransformation[3];

	for (const glm::vec3& contourPoint : contourPoints) {
		float minDist = FLT_MAX;
		ContourBinding bestBinding;

		glm::vec3 rootToContour = glm::normalize(contourPoint - rootPos);

		for (auto& [parent, child, index] : segments) {
			glm::vec3 P = parent->globalTransformation[3];
			glm::vec3 Q = child->globalTransformation[3];

			glm::vec3 closest = SceneNode::intersectionPoint(P, Q, contourPoint);
			float t = glm::dot(Q - P, contourPoint - P) / glm::dot(Q - P, Q - P);
			//t = glm::clamp(t, 0.0f, 1.0f);
			t = 1.0f;
			float dist = glm::length(closest - contourPoint);

			if (dist < minDist) {
				minDist = dist;
				std::fill(weights.begin(), weights.end(), 0);
				// initialize weight to 1 for the branch is it binded to
				weights[index] = 1.0f;
				bestBinding = { parent, child, contourPoint, t, closest, glm::inverse(t * child->globalTransformation + (1 - t) * parent->globalTransformation), weights };
			}
		}

		bindings.push_back(bestBinding);
	}

	// previous animatrion matrix for all branches
	for (int i = 0; i < segments.size(); i++) {
		previousAnimateInverseVec[i] = glm::inverse(bindings[i].t * std::get<1>(segments[i])->globalTransformation + (1 - bindings[i].t) * std::get<0>(segments[i])->globalTransformation);
	}

	// for each contour, store the previous animation inverse matrix list (for all branches)
	for (int i = 0; i < bindings.size(); i++) {
		bindings[i].previousAnimateInverseMat = previousAnimateInverseVec;
	}

	return bindings;
}

// go through contour points in a loop, from index 1 to n - 1 (not from 0 to n)
// for each contour point, get its neighbors
// take the average weight for each weight in the weights vector and update the current contour point weight
void SceneNode::multipleWeights(std::vector<ContourBinding>& bindings) {
	std::vector<std::vector<float>> weightCopies;
	for (int i = 0; i < bindings.size(); i++) {
		weightCopies.push_back(bindings[i].weights);
	}
	// for each contour point
	for (int i = 1; i < bindings.size() - 1; i++) {
		// for each weight of the contour point
		for (int j = 0; j < bindings[i].weights.size(); j++) {
			// take the average (direct averaging or using rate of blending)
			//float average = (weightCopies[i - 1][j] + weightCopies[i][j] + weightCopies[i + 1][j]) / 3.0f;
			float average = 0.2*(((weightCopies[i - 1][j] + weightCopies[i][j] + weightCopies[i + 1][j]) / 3.0f) - weightCopies[i][j]) + weightCopies[i][j];
			bindings[i].weights[j] = average;
		}
	}
}

// map middle point to branching point
void SceneNode::bindToBranchingPoint(std::vector<ContourBinding>& bindings, std::vector<std::tuple<SceneNode*, SceneNode*, int>>& segments) {
	// get all the "leaf" branches
	std::vector<std::tuple<SceneNode*, SceneNode*, int>> endSegments;
	for (int i = 0; i < segments.size(); i++) {
		if (std::get<1>(segments[i])->children.empty()) {
			endSegments.push_back(segments[i]);
		}
	}

	// for all leaf branches
	for (int i = 1; i < endSegments.size(); i++) {
		std::vector<int> contourPointsToConsider;
		// if they are two consecutive branches (same parent -> index differ by 1)
		if (abs(std::get<2>(endSegments[i - 1]) - std::get<2>(endSegments[i])) == 1 && std::get<0>(endSegments[i - 1]) == std::get<0>(endSegments[i])) {
			// out of all the contour points that are binded to either of these two segments
			for (int j = 0; j < bindings.size() - 1; j++) {
				// if the previous binding point of the first contour point is the tip of the first branch and the previous binding point of the second contour point is the tip of the second point
				if (glm::vec4(bindings[j].closestPoint, 1.f) == std::get<1>(endSegments[i - 1])->globalTransformation[3] && glm::vec4(bindings[j + 1].closestPoint, 1.f) == std::get<1>(endSegments[i])->globalTransformation[3]) {
					// bind the middle point to the branching point
					for (int k = 0; k < segments.size(); k++) {
						// if the current branch's child is the same as the parent of the currently connected branch (so we have found the ancestor)
						if (std::get<1>(segments[k]) == std::get<0>(endSegments[i - 1])) {
							bindings[j].childNode = std::get<0>(endSegments[i - 1]);
							bindings[j].parentNode = std::get<0>(segments[k]);
							bindings[j].t = 1.f;
							// closest point is now the branching point
							bindings[j].closestPoint = (bindings[j].t * bindings[j].childNode->globalTransformation[3] + (1 - bindings[j].t) * bindings[j].parentNode->globalTransformation[3]);
							// update weight to reflect the change in binding
							std::fill(bindings[j].weights.begin(), bindings[j].weights.end(), 0);
							bindings[j].weights[std::get<2>(segments[k])] = 1.f;
							bindings[j].previousAnimateInverse = glm::inverse(std::get<0>(endSegments[i - 1])->globalTransformation);
							bindings[j].previousAnimateInverseMat[std::get<2>(segments[k])] = glm::inverse(std::get<0>(endSegments[i - 1])->globalTransformation);
							break;
						}
					}
				}
			}
		}
	}

	// first and last contour points are mapped to the root
	bindings[0].parentNode = std::get<0>(segments[0]);
	bindings[0].childNode = std::get<1>(segments[0]);
	bindings[0].t = 1.f;
	bindings[0].closestPoint = (bindings[0].t * bindings[0].childNode->globalTransformation[3] + (1 - bindings[0].t) * bindings[0].parentNode->globalTransformation[3]);
	std::fill(bindings[0].weights.begin(), bindings[0].weights.end(), 0);
	// assign weight to the first branch (root stem branch)
	bindings[0].weights[0] = 1.f;
	bindings[0].previousAnimateInverse = glm::inverse(bindings[0].t * bindings[0].childNode->globalTransformation + (1 - bindings[0].t) * bindings[0].parentNode->globalTransformation);
	bindings[0].previousAnimateInverseMat[0] = glm::inverse(bindings[0].t * bindings[0].childNode->globalTransformation + (1 - bindings[0].t) * bindings[0].parentNode->globalTransformation);

	bindings[bindings.size() - 1].parentNode = std::get<0>(segments[0]);
	bindings[bindings.size() - 1].childNode = std::get<1>(segments[0]);
	bindings[bindings.size() - 1].t = 1.f;
	bindings[bindings.size() - 1].closestPoint = (bindings[bindings.size() - 1].t * bindings[bindings.size() - 1].childNode->globalTransformation[3] + (1 - bindings[bindings.size() - 1].t) * bindings[0].parentNode->globalTransformation[3]);
	std::fill(bindings[bindings.size() - 1].weights.begin(), bindings[bindings.size() - 1].weights.end(), 0);
	// assign weight to the first branch (root stem branch)
	bindings[bindings.size() - 1].weights[0] = 1.f;
	bindings[bindings.size() - 1].previousAnimateInverse = glm::inverse(bindings[bindings.size() - 1].t * bindings[0].childNode->globalTransformation + (1 - bindings[bindings.size() - 1].t) * bindings[0].parentNode->globalTransformation);
	bindings[bindings.size() - 1].previousAnimateInverseMat[0] = glm::inverse(bindings[bindings.size() - 1].t * bindings[bindings.size() - 1].childNode->globalTransformation + (1 - bindings[0].t) * bindings[bindings.size() - 1].parentNode->globalTransformation);
}

// relative transformation between frame (global frame) but taking into account blending of weights for multiple branches
void SceneNode::animationPerFrameUsingMultipleWeights(std::vector<ContourBinding>& bindings, std::vector<std::tuple<SceneNode*, SceneNode*, int>>& segments) {
	glm::vec3 animatedPoint;

	// for each contour point, blend the transformation using the weight vector
	for (auto& binding : bindings) {
		animatedPoint = glm::vec3(0.f);
		// for each branch, get the branch trasnformation and blend using weight vector
		for (int i = 0; i < segments.size(); i++) {
			// weight * A' * inverseA * contour point
			animatedPoint += glm::vec3(binding.weights[i] * (binding.t * std::get<1>(segments[i])->globalTransformation + (1 - binding.t) * std::get<0>(segments[i])->globalTransformation) * binding.previousAnimateInverseMat[i] * glm::vec4(binding.contourPoint, 1.0f));
			binding.previousAnimateInverseMat[i] = glm::inverse((binding.t * std::get<1>(segments[i])->globalTransformation) + (1 - binding.t) * (std::get<0>(segments[i])->globalTransformation));
		}
		binding.contourPoint = animatedPoint;
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
