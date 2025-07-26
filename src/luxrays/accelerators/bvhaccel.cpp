/***************************************************************************
 * Copyright 1998-2020 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxCoreRender.                                   *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

// Boundary Volume Hierarchy accelerator
// Based of "Efficiency Issues for Ray Tracing" by Brian Smits
// Available at http://www.cs.utah.edu/~bes/papers/fastRT/paper.html

#include <iostream>
#include <functional>
#include <algorithm>
#include <limits>

#include "luxrays/accelerators/bvhaccel.h"
#include "luxrays/utils/utils.h"
#include "luxrays/core/context.h"
#include "luxrays/core/geometry/vector.h"
#include "luxrays/core/geometry/transform.h"

// 🔥GRIN
#include "luxrays/core/geometry/xprimeray.h"
// 🔥GRIN
#include <cmath>  // for sinf, cosf, etc.


using namespace std;

namespace luxrays {

// 🔥GRIN Trial Function
bool GRINRK4_Intersect(
    const Ray &ray,
    const Point &p0,
    const Point &p1,
    const Point &p2,
    const float stepSize,
    const int maxSteps,
    float *hitT,
    float *b1,
    float *b2
) {
    const Vector dir0 = Normalize(ray.d);
    Point pos = ray.o;
	const float bendScale = 0.01f; // or 0.05f for less curvature

	for (int i = 0; i < maxSteps; ++i) {
		const Vector posVec = Vector(pos);
		const float r = luxrays::Max(0.001f, posVec.Length());

		const Vector bend = Normalize(posVec) * (bendScale / r);

		const Vector k1 = dir0;
		const Vector k2 = Normalize(dir0 + bend * stepSize * 0.5f);

		const Point midPos = pos + k2 * (stepSize * 0.5f);
		const float r2 = luxrays::Max(0.001f, Vector(midPos).Length());

		const Vector bend2 = Normalize(Vector(midPos)) * (bendScale / r2);
		const Vector k3 = Normalize(dir0 + bend2 * stepSize * 0.5f);

		const Vector k4 = Normalize(dir0 + bend * stepSize);

		const Vector avgDir = Normalize((k1 + 2.f * k2 + 2.f * k3 + k4) / 6.f);
		const Point next = pos + avgDir * stepSize;

		const Ray segment(pos, next - pos, 0.f, 1.f);

		float localT, localB1, localB2;
		if (Triangle::Intersect(segment, p0, p1, p2, &localT, &localB1, &localB2)) {
			if (hitT) *hitT = localT;
			if (b1) *b1 = localB1;
			if (b2) *b2 = localB2;
			return true;
		}

		pos = next;
	}
    return false;
}
// 🔥GRIN Trial Function


// BVHAccel Method Definitions

BVHAccel::BVHAccel(const Context *context) : ctx(context) {
	params = ToBVHParams(ctx->GetConfig());

	initialized = false;
}

BVHAccel::~BVHAccel() {
	if (initialized)
		delete[] bvhTree;
}

BVHParams BVHAccel::ToBVHParams(const Properties &props) {
	// Tree type to generate (2 = binary, 4 = quad, 8 = octree)
	const int treeType = props.Get(Property("accelerator.bvh.treetype")(4)).Get<int>();
	// Samples to get for cost minimization
	const int costSamples = props.Get(Property("accelerator.bvh.costsamples")(0)).Get<int>();
	const int isectCost = props.Get(Property("accelerator.bvh.isectcost")(80)).Get<int>();
	const int travCost = props.Get(Property("accelerator.bvh.travcost")(10)).Get<int>();
	const float emptyBonus = props.Get(Property("accelerator.bvh.emptybonus")(.5)).Get<float>();
	
	BVHParams params;
	// Make sure treeType is 2, 4 or 8
	if (treeType <= 2) params.treeType = 2;
	else if (treeType <= 4) params.treeType = 4;
	else params.treeType = 8;

	params.costSamples = costSamples;
	params.isectCost = isectCost;
	params.traversalCost = travCost;
	params.emptyBonus = emptyBonus;

	return params;
}

void BVHAccel::Init(const deque<const Mesh *> &ms, const u_longlong totVert,
		const u_longlong totTri) {
	assert (!initialized);

	meshes = ms;
	totalVertexCount = totVert;
	totalTriangleCount = totTri;

	// Handle the empty DataSet case
	if (totalTriangleCount == 0) {
		LR_LOG(ctx, "Empty BVH");
		nNodes = 0;
		bvhTree = NULL;
		initialized = true;

		return;
	}

	const double t0 = WallClockTime();

	//--------------------------------------------------------------------------
	// Build the list of triangles
	//--------------------------------------------------------------------------
	
	vector<BVHTreeNode> bvNodes(totalTriangleCount);
	vector<BVHTreeNode *> bvList(totalTriangleCount, NULL);
	u_int meshIndex = 0;
	u_int bvListIndex = 0;
	BOOST_FOREACH(const Mesh *mesh, meshes) {
		const Triangle *p = mesh->GetTriangles();
		const u_int triangleCount = mesh->GetTotalTriangleCount();

		#pragma omp parallel for
		for (
				// Visual C++ 2013 supports only OpenMP 2.5
#if _OPENMP >= 200805
				unsigned
#endif
				int i = 0; i < triangleCount; ++i) {
			const int index = bvListIndex + i;
			BVHTreeNode *node = &bvNodes[index];

			// This is a fast path because I know Mesh can be only TYPE_TRIANGLE/TYPE_EXT_TRIANGLE
			// in BVH
			node->bbox = Union(
					BBox(
						mesh->GetVertex(Transform::TRANS_IDENTITY, p[i].v[0]),
						mesh->GetVertex(Transform::TRANS_IDENTITY, p[i].v[1])),
						mesh->GetVertex(Transform::TRANS_IDENTITY, p[i].v[2]));
			// NOTE - Ratow - Expand bbox a little to make sure rays collide
			node->bbox.Expand(MachineEpsilon::E(node->bbox));
			node->triangleLeaf.meshIndex = meshIndex;
			node->triangleLeaf.triangleIndex = i;

			node->leftChild = NULL;
			node->rightSibling = NULL;

			bvList[index] = node;
		}
		
		bvListIndex += triangleCount;
		++meshIndex;
	}

	LR_LOG(ctx, "BVH Dataset preprocessing time: " << int((WallClockTime() - t0) * 1000) << "ms");

	//--------------------------------------------------------------------------
	// Build the BVH hierarchy
	//--------------------------------------------------------------------------

	const double t1 = WallClockTime();

	const string builderType = ctx->GetConfig().Get(Property("accelerator.bvh.builder.type")(
		"EMBREE_BINNED_SAH"
		)).Get<string>();

	LR_LOG(ctx, "BVH builder: " << builderType);
	if (builderType == "CLASSIC")
		bvhTree = BuildBVH(params, &nNodes, &meshes, bvList);
	else if (builderType == "EMBREE_BINNED_SAH")
		bvhTree = BuildEmbreeBVHBinnedSAH(params, &nNodes, &meshes, bvList);
	else if (builderType == "EMBREE_MORTON")
		bvhTree = BuildEmbreeBVHMorton(params, &nNodes, &meshes, bvList);
	else
		throw runtime_error("Unknown BVH builder type in BVHAccel::Init(): " + builderType);

	LR_LOG(ctx, "BVH build hierarchy time: " << int((WallClockTime() - t1) * 1000) << "ms");

	//--------------------------------------------------------------------------
	// Done
	//--------------------------------------------------------------------------

	LR_LOG(ctx, "BVH total build time: " << int((WallClockTime() - t0) * 1000) << "ms");
	LR_LOG(ctx, "Total BVH memory usage: " << nNodes * sizeof(luxrays::ocl::BVHArrayNode) / 1024 << "Kbytes");

	initialized = true;
}


// 🔥GRIN Straight-Line Ray Decisions
//bool BVHAccel::Intersect(const Ray *initialRay, RayHit *rayHit, const float beta, const luxrays::Vector &gamma) const {
bool BVHAccel::Intersect(const Ray *initialRay, RayHit *rayHit) const {
	assert (initialized);

	rayHit->t = initialRay->maxt;
	rayHit->SetMiss();
	if (!nNodes)
		return false;

	Ray ray(*initialRay);
	// -- Convert Ray to xPRIMEray --
	xPRIMEray xPRIMEray(
			initialRay->o,                 // origin
			initialRay->d,                 // direction
			luxrays::Point(0.f, 0.f, 0.f), // center of curvature (can customize)
			1.0f,                          // beta (GRIN intensity scalar)
			luxrays::Vector(1.f,1.f,1.f),  // gamma (exponents per axis)
			xPRIMErayType::POWER,          // curvature model
			initialRay->mint,
			initialRay->maxt
	);
	
	// 🔥GRIN
	//std::cout <<  "🔥[GRIN] Ray origin: " << ray.o << ", dir: " << ray.d << ", maxt: " << ray.maxt << std::endl;
	//LR_LOG(ctx, "🔥[GRIN] Ray origin: " << ray.o << ", dir: " << ray.d << ", maxt: " << ray.maxt);

	u_int currentNode = 0; // Root Node
	const u_int stopNode = BVHNodeData_GetSkipIndex(bvhTree[0].nodeData); // Non-existent

	float t, b1, b2;
	while (currentNode < stopNode) {
		const luxrays::ocl::BVHArrayNode &node = bvhTree[currentNode];

		const u_int nodeData = node.nodeData;
		if (BVHNodeData_IsLeaf(nodeData)) {
			// It is a leaf, check the triangle
			const Mesh *mesh = meshes[node.triangleLeaf.meshIndex];

			// This is a fast path because I know Mesh can be only TYPE_TRIANGLE/TYPE_EXT_TRIANGLE
			// in BVH
			const Point p0 = mesh->GetVertex(Transform::TRANS_IDENTITY, node.triangleLeaf.v[0]);
			const Point p1 = mesh->GetVertex(Transform::TRANS_IDENTITY, node.triangleLeaf.v[1]);
			const Point p2 = mesh->GetVertex(Transform::TRANS_IDENTITY, node.triangleLeaf.v[2]);


			// 🔥GRIN Straight-Line Hit Operation
			//if (Triangle::xPRIMEIntersect(xPRIMEray, p0, p1, p2, &t, &b1, &b2)) {
			//if (Triangle::xPRIMEIntersect(xPRIMEray, p0, p1, p2, &t, &b1, &b2)) {
			if (Triangle::Intersect(ray, p0, p1, p2, &t, &b1, &b2)) {
				if (t < rayHit->t) {
					ray.maxt = t;
					rayHit->t = t;
					rayHit->b1 = b1;
					rayHit->b2 = b2;
					rayHit->meshIndex = node.triangleLeaf.meshIndex;
					rayHit->triangleIndex = node.triangleLeaf.triangleIndex;
					// Continue testing for closer intersections
				}
			}
			++currentNode;
		} else {
			// It is a node, check the bounding box
			if (BBox::IntersectP(ray,
					*reinterpret_cast<const Point *>(&node.bvhNode.bboxMin[0]),
					*reinterpret_cast<const Point *>(&node.bvhNode.bboxMax[0])))
				++currentNode;
			else {
				// I don't need to use BVHNodeData_GetSkipIndex() here because
				// I already know the leaf flag is 0
				currentNode = nodeData;
			}
		}
	}

	return !rayHit->Miss();
}

bool BVHAccel::xPRIMEIntersect(const Ray *initialRay, RayHit *rayHit,
               const float beta, const luxrays::Vector &gamma,
               const luxrays::Point &grinCenter,
			   const float stepSize, const int numSteps) const {
	assert (initialized);

	rayHit->t = initialRay->maxt;
	rayHit->SetMiss();
	if (!nNodes)
		return false;

	Ray ray(*initialRay);
	// -- Convert Ray to xPRIMEray --
	xPRIMEray xPRIMEray(
		initialRay->o,              // origin
		initialRay->d,				// direction
		grinCenter,					// center of curvature (can customize)
		beta,                       // beta (GRIN intensity scalar)
		gamma,						// gamma (exponents per axis)
		xPRIMErayType::POWER,       // curvature model
		initialRay->mint,
		initialRay->maxt,
		stepSize,
		numSteps
	);

	// Use a scaled step count when computing the curved-ray envelope to
	// capture sharper bends for long rays.
	const int envelopeSteps = xPRIMEray.numSteps * 2;
	const BBox sweptBBox = ComputeSweptBBox(xPRIMEray, envelopeSteps);

	// 🔥GRIN
	//std::cout <<  "🔥[GRIN] Ray origin: " << ray.o << ", dir: " << ray.d << ", maxt: " << ray.maxt << std::endl;
	//LR_LOG(ctx, "🔥[GRIN] Ray origin: " << ray.o << ", dir: " << ray.d << ", maxt: " << ray.maxt);

	u_int currentNode = 0; // Root Node
	const u_int stopNode = BVHNodeData_GetSkipIndex(bvhTree[0].nodeData); // Non-existent

	float t, b1, b2;
	while (currentNode < stopNode) {
		const luxrays::ocl::BVHArrayNode &node = bvhTree[currentNode];

		const u_int nodeData = node.nodeData;
		if (BVHNodeData_IsLeaf(nodeData)) {
			// It is a leaf, check the triangle
			const Mesh *mesh = meshes[node.triangleLeaf.meshIndex];

			// This is a fast path because I know Mesh can be only TYPE_TRIANGLE/TYPE_EXT_TRIANGLE
			// in BVH
			const Point p0 = mesh->GetVertex(Transform::TRANS_IDENTITY, node.triangleLeaf.v[0]);
			const Point p1 = mesh->GetVertex(Transform::TRANS_IDENTITY, node.triangleLeaf.v[1]);
			const Point p2 = mesh->GetVertex(Transform::TRANS_IDENTITY, node.triangleLeaf.v[2]);

			// 🔥GRIN Straight-Line Hit Operation
			if (Triangle::xPRIMEIntersect(xPRIMEray, p0, p1, p2, grinCenter, &t, &b1, &b2)) {
			//if (Triangle::Intersect(ray, p0, p1, p2, &t, &b1, &b2)) {
				if (t < rayHit->t) {
					ray.maxt = t;
					rayHit->t = t;
					rayHit->b1 = b1;
					rayHit->b2 = b2;
					rayHit->meshIndex = node.triangleLeaf.meshIndex;
					rayHit->triangleIndex = node.triangleLeaf.triangleIndex;
					// Continue testing for closer intersections
				}
			}
			++currentNode;
		} else {
			// It is a node, check against the curved-ray envelope
			const BBox nodeBBox(
						*reinterpret_cast<const Point *>(&node.bvhNode.bboxMin[0]),
						*reinterpret_cast<const Point *>(&node.bvhNode.bboxMax[0]));
			if (sweptBBox.Overlaps(nodeBBox))
				++currentNode;
			else {
				// I don't need to use BVHNodeData_GetSkipIndex() here because
				// I already know the leaf flag is 0
				currentNode = nodeData;
			}
		}
	}

	return !rayHit->Miss();
}

}
