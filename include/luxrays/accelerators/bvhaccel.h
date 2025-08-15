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

#ifndef _LUXRAYS_BVHACCEL_H
#define	_LUXRAYS_BVHACCEL_H

#include <vector>
#include <boost/foreach.hpp>

#include "luxrays/luxrays.h"
#include "luxrays/core/accelerator.h"
#include "slg/scene/stitchhint.h"
#include "luxrays/core/bvh/bvhbuild.h"
#include "luxrays/core/geometry/triangle.h"

namespace luxrays {

class HardwareIntersectionDevice;

// BVHAccel Declarations
class BVHAccel : public Accelerator {
public:
	// BVHAccel Public Methods
	BVHAccel(const Context *context);
	virtual ~BVHAccel();

	virtual AcceleratorType GetType() const { return ACCEL_BVH; }

	virtual bool HasNativeSupport(const IntersectionDevice &device) const;
	virtual bool HasHWSupport(const IntersectionDevice &device) const;

	virtual HardwareIntersectionKernel *NewHardwareIntersectionKernel(HardwareIntersectionDevice &device) const;

	virtual void Init(const std::deque<const Mesh *> &meshes,
		const u_longlong totalVertexCount,
		const u_longlong totalTriangleCount);

	virtual bool Intersect(const Ray *ray, RayHit *hit) const;

virtual bool xPRIMEIntersect(const Ray *ray, RayHit *hit,
							const float beta, const luxrays::Vector &gamma,
							const luxrays::Point &grinCenter,
							const float rInner, const float rOuter,
							const float stepSize, const int numSteps,
							const bool invert = false,
							const float insightCurvatureThreshold = 1e-6f,
							const float barycentricEpsilon = 0.03f,
							const float rk4PlaneThreshold = 1e-4f,
							const float uvSeamTolerance = 1e-6f,
							const UVCrossPolicy uvPolicy = UV_REJECT,
							slg::StitchHint *stitchHint = nullptr,
							float stitchPlaneFactor = 2.f,
							float stitchBaryMargin  = 0.02f,
							// Adaptive controls
							bool  adaptiveEnable = false,
							float adaptivePlaneTriggerFactor = 1.0f,
							float adaptiveCurvatureTrigger = 0.2f,
							int   adaptiveMaxSubdiv = 2,
							int   adaptiveBisectIters = 5,
							float adaptiveMinStep = 1e-5f,
							float adaptiveInsightAcceptMargin = 0.0f,
							float adaptiveRate = 0.25f,
							float adaptiveMaxScale = 4.0f) const;

	static BVHParams ToBVHParams(const Properties &props);

	friend class BVHKernel;
	friend class MBVHKernel;
	friend class MBVHAccel;

private:
	BVHParams params;

	u_int nNodes;
	luxrays::ocl::BVHArrayNode *bvhTree;

	const Context *ctx;
	std::deque<const Mesh *> meshes;
	u_longlong totalVertexCount, totalTriangleCount;

	bool initialized;
};

}

#endif	/* _LUXRAYS_BVHACCEL_H */
