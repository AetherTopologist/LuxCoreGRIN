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

#ifndef _LUXRAYS_TRIANGLE_H
#define	_LUXRAYS_TRIANGLE_H

#include <algorithm>
#include <cmath>


#include "luxrays/luxrays.h"
#include "luxrays/core/geometry/point.h"
#include "luxrays/core/geometry/vector.h"
#include "luxrays/core/geometry/normal.h"
#include "luxrays/core/geometry/ray.h"
#include "luxrays/core/geometry/bbox.h"
#include "luxrays/utils/mc.h"
#include "luxrays/utils/serializationutils.h"

#include "luxrays/core/geometry/xprimeray.h" // GRIN PRIME Ray Class


namespace luxrays {

// UV seam handling policy
enum UVCrossPolicy { UV_REJECT = 0, UV_EDGE_PROJECT = 1 };

struct BaryResult {
    bool inside;          // strictly inside triangle (with epsilon)
    bool nearEdge;        // within tolerance of any edge
    float b1, b2;         // barycentrics (b0 = 1 - b1 - b2)
    u_int edgeId;         // 0:(p0,p1) 1:(p1,p2) 2:(p2,p0), or 3 if none
};

// Project a 3D point (curved hit) to triangle plane and compute robust barycentrics.
// Uses an orthonormal basis, returns nearEdge if within edgeTol of an edge.
static inline BaryResult ProjectToTriangleBary(
    const Point &p0, const Point &p1, const Point &p2,
    const Point &hit, const float insideTol, const float edgeTol) {
    BaryResult br;
    br.inside = false;
    br.nearEdge = false;
    br.b1 = br.b2 = 0.f;
    br.edgeId = 3;

    // Build orthonormal basis on the triangle plane
    const Vector e1 = p1 - p0;
    const Vector e2 = p2 - p0;
    const Vector Nraw = Cross(e1, e2);
    const float nlen2 = Nraw.LengthSquared();
    if (nlen2 < 1e-20f)
        return br; // Degenerate triangle, leave br.inside=false
    const Vector N = Normalize(Nraw);
    const Vector u = Normalize(e1);
    const Vector v = Normalize(Cross(N, u));

    // Project hit to plane
    const Vector hpv = hit - p0;
    const float dist = Dot(hpv, N);
    const Point hp = hit - dist * N;

    // 2D coordinates in basis
    const double x1 = Dot(e1, u);
    const double y1 = Dot(e1, v);
    const double x2 = Dot(e2, u);
    const double y2 = Dot(e2, v);
    const double xh = Dot(hp - p0, u);
    const double yh = Dot(hp - p0, v);

    const double denom = x1 * y2 - x2 * y1;
    if (denom != 0.0) {
        const double b1d = (xh * y2 - x2 * yh) / denom;
        const double b2d = (x1 * yh - xh * y1) / denom;
        br.b1 = static_cast<float>(b1d);
        br.b2 = static_cast<float>(b2d);
        const double b0d = 1.0 - b1d - b2d;

        // inside check with tolerance
        if (b1d >= -insideTol && b2d >= -insideTol && b0d >= -insideTol &&
            b1d <= 1.0 + insideTol && b2d <= 1.0 + insideTol && b0d <= 1.0 + insideTol)
            br.inside = true;
    }

    // Edge distance check for seam handling
    const Point pts[3] = { p0, p1, p2 };
    for (u_int i = 0; i < 3; ++i) {
        const Point &a = pts[i];
        const Point &b = pts[(i + 1) % 3];
        const Vector ab = b - a;
        const double abLen2 = ab.LengthSquared();
        if (abLen2 == 0.0)
            continue;
        double t = Dot(hp - a, ab) / abLen2;
        t = std::clamp(t, 0.0, 1.0);
        const Point proj = a + ab * t;
        const double d = Distance(hp, proj);
        if (d <= edgeTol) {
            br.nearEdge = true;
            br.edgeId = i;
            break;
        }
    }

    return br;
}

// If near an edge, project 'hit' onto closest edge segment in-plane,
// recompute barycentrics restricted to that edge (one coord = 0).
static inline bool EdgeProjectBary(
    const Point &p0, const Point &p1, const Point &p2,
    const Point &hit, u_int edgeId, float *b1, float *b2) {
    // Compute plane normal and project hit onto the plane
    const Vector e10 = p1 - p0;
    const Vector e20 = p2 - p0;
    const Vector Nraw = Cross(e10, e20);
    const float nlen2 = Nraw.LengthSquared();
    if (nlen2 < 1e-20f)
        return false;
    const Vector N = Normalize(Nraw);
    const Point hp = hit - Dot(hit - p0, N) * N;

    const Point *a = nullptr;
    const Point *b = nullptr;
    switch (edgeId) {
        case 0: a = &p0; b = &p1; break;
        case 1: a = &p1; b = &p2; break;
        case 2: a = &p2; b = &p0; break;
        default: return false;
    }

    const Vector ab = *b - *a;
    const double abLen2 = ab.LengthSquared();
    if (abLen2 == 0.0)
        return false;
    double t = Dot(hp - *a, ab) / abLen2;
    t = std::clamp(t, 0.0, 1.0);
    const double wA = 1.0 - t;
    const double wB = t;

    switch (edgeId) {
        case 0: // p0-p1
            *b1 = static_cast<float>(wB);
            *b2 = 0.f;
            break;
        case 1: // p1-p2
            *b1 = static_cast<float>(wA);
            *b2 = static_cast<float>(wB);
            break;
        case 2: // p2-p0
            *b1 = 0.f;
            *b2 = static_cast<float>(wA);
            break;
    }
    return true;
}

// OpenCL data types
namespace ocl {
#include "luxrays/core/geometry/triangle_types.cl"
}

class Triangle {
public:
	Triangle() { }
	Triangle(const unsigned int v0, const unsigned int v1, const unsigned int v2) {
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
	}

	BBox WorldBound(const Point *verts) const {
		const Point &p0 = verts[v[0]];
		const Point &p1 = verts[v[1]];
		const Point &p2 = verts[v[2]];

		return Union(BBox(p0, p1), p2);
	}

	// Analytic intersection of POWER curved rays with a triangle's plane using symbolic per-axis gamma
	static bool IntersectPlaneSymbolic(const xPRIMEray &ray, const Point &p0, const Vector &normal,
		float *tHit, Point *hitPoint, const float insightCurvatureThreshold = 1e-6f) {
		
		if (ray.type != xPRIMErayType::POWER)
			return false;

		const Vector rayToPlane = ray.origin - p0;
		const float A = Dot(rayToPlane, normal);
		const float B = Dot(ray.direction, normal);
		const float beta = ray.beta;

		// Early out: Avoid division by zero or negligible curvature
		const float denom = beta * B;
		if (std::fabs(denom) < insightCurvatureThreshold)
			return false;
		
		const float tBase = -A / denom;
		if (tBase <= 0.f) {
			if (ray.beta < 0.f && B < 0.f) {
				// Allow negative beta with ray pointing toward the plane
			} else
				return false;
		}

		// Use axis of maximum curvature to determine t exponent
		const float gammaMax = std::max(1e-6f,
				std::max(ray.gamma.x, std::max(ray.gamma.y, ray.gamma.z)));
		const float t = std::pow(tBase, 1.f / gammaMax);

		if (!std::isfinite(t) || t < ray.mint || t > ray.maxt)
			return false;

		// Compute curved offset from origin using per-axis gamma
		const Vector curveOffset(
			ray.direction.x * std::pow(t, ray.gamma.x),
			ray.direction.y * std::pow(t, ray.gamma.y),
			ray.direction.z * std::pow(t, ray.gamma.z));

		*hitPoint = ray.origin + beta * curveOffset;
		*tHit = t;

		return true;
	}

	static Vector ComputeGRINField(
			const Point &pos,
			const float beta,
			const Vector &gamma,
			const Point &GRINCenter,
			const float rInner,
			const float rOuter,
			const bool invert = false) {

		const Vector offset = pos - GRINCenter;
		const float r = offset.Length();

		// Straight-line region inside rInner
		if (r < rInner)
			return Vector(0.f, 0.f, 0.f);

		// Normalize r into t ∈ [0, 1]
		const float t = std::clamp((r - rInner) / (rOuter - rInner), 0.f, 1.f);
		const float t_x = std::pow(t, gamma.x);
		const float t_y = std::pow(t, gamma.y);
		const float t_z = std::pow(t, gamma.z);

		const float effectiveBeta = invert ? -beta : beta;

		return Vector(
				effectiveBeta * t_x * offset.x / r,
				effectiveBeta * t_y * offset.y / r,
				effectiveBeta * t_z * offset.z / r);
	}

	static bool IntersectINSIGHT(
			const xPRIMEray &ray,
			const Point &p0,
			const Vector &normal,
			float *tINSIGHT,
			Point *approxHit,
			const float insightCurvatureThreshold = 1e-6f) {
		
		if (ray.type != xPRIMErayType::POWER)
			return false;

		const Vector rayToPlane = ray.origin - p0;
		const float A = Dot(rayToPlane, normal);
		const float B = Dot(ray.direction, normal);
		const float beta = ray.beta;

		// Early out for negligible curvature
		const float denom = beta * B;
		if (std::fabs(denom) < insightCurvatureThreshold)
			return false;
		
		const float tBase = -A / denom;
		if (tBase <= 0.f) {
			if (ray.beta < 0.f && B < 0.f) {
				// Allow negative beta with ray pointing toward the plane
			} else
				return false;
		}

		// Use strongest curvature axis for exponent
		const float gammaMax = std::max(1e-6f,
				std::max(ray.gamma.x, std::max(ray.gamma.y, ray.gamma.z)));
		const float t = std::pow(tBase, 1.f / gammaMax);

		if (!std::isfinite(t) || t < ray.mint || t > ray.maxt)
			return false;

		// Compute approximate curved hit point
		const Vector curveOffset(
			ray.direction.x * std::pow(t, ray.gamma.x),
			ray.direction.y * std::pow(t, ray.gamma.y),
			ray.direction.z * std::pow(t, ray.gamma.z));

		*approxHit = ray.origin + beta * curveOffset;
		*tINSIGHT = t;

		return true;
	}

	// Project point P to the closest point on triangle (p0,p1,p2). Returns clamped barycentrics.
	static inline void ClosestPointBarycentric(const Point &p,
					const Point &p0, const Point &p1, const Point &p2,
					float *b1, float *b2) {
		// Ericson-style closest-point on triangle
		const Vector v0 = p1 - p0;
		const Vector v1 = p2 - p0;
		const Vector v2 = p - p0;

		const float d00 = Dot(v0, v0);
		const float d01 = Dot(v0, v1);
		const float d11 = Dot(v1, v1);
		const float d20 = Dot(v2, v0);
		const float d21 = Dot(v2, v1);
		const float denom = d00 * d11 - d01 * d01;

		float v = (d11 * d20 - d01 * d21) / denom;
		float w = (d00 * d21 - d01 * d20) / denom;
		float u = 1.f - v - w;

		// If inside, done
		if (u >= 0.f && v >= 0.f && w >= 0.f) {
			*b1 = v; *b2 = w;
			return;
		}

		// Otherwise clamp to edges
		auto clamp01 = [](float x) { return (x < 0.f) ? 0.f : (x > 1.f ? 1.f : x); };

		// Edge p0-p1
		{
			const Vector e = v0;
			float t = clamp01(Dot(v2, e) / Dot(e, e));
			Point q = p0 + t * e;
			Vector rq = p - q;
			float d = Dot(rq, rq);
			// store best in locals
		}

		// We’ll do simple bary clamp: move to nearest of 3 edges by analytic projection
		// Edge p0-p1
		{
			const Vector e = v0;
			float t = clamp01(Dot(v2, e) / Dot(e, e));
			float vv = t, ww = 0.f;
			float uu = 1.f - vv - ww;
			if (uu >= 0.f && vv >= 0.f && ww >= 0.f) { *b1 = vv; *b2 = ww; return; }
		}
		// Edge p0-p2
		{
			const Vector e = v1;
			float t = clamp01(Dot(v2, e) / Dot(e, e));
			float vv = 0.f, ww = t;
			float uu = 1.f - vv - ww;
			if (uu >= 0.f && vv >= 0.f && ww >= 0.f) { *b1 = vv; *b2 = ww; return; }
		}
		// Edge p1-p2
		{
			const Vector e = v1 - v0;
			const Vector w2 = p - p1;
			float t = clamp01(Dot(w2, e) / Dot(e, e));
			// On edge p1->p2, bary: u=0, v=1-t, w=t
			float vv = 1.f - t, ww = t;
			if (vv >= 0.f && ww >= 0.f) { *b1 = vv; *b2 = ww; return; }
		}

		// Fallback clamp to [0,1] and renormalize
		*b1 = clamp01(v);
		*b2 = clamp01(w);
		float u2 = 1.f - *b1 - *b2;
		if (u2 < 0.f) { // push back into simplex
			if (*b1 > *b2) *b1 = clamp01(1.f - *b2); else *b2 = clamp01(1.f - *b1);
		}
	}

	static inline void ClampBarycentricSoft(const Point &p,
					const Point &p0, const Point &p1, const Point &p2,
					float tol, float *b1, float *b2) {
		const float u = 1.f - *b1 - *b2;
		if (u >= -tol && *b1 >= -tol && *b2 >= -tol)
			ClosestPointBarycentric(p, p0, p1, p2, b1, b2);
	}

	static bool RK4_GRINIntersect(
							const xPRIMEray &ray,
							const Point &p0,
							const Point &p1,
							const Point &p2,
							const Point &grinCenter,
							const float rInner,
							const float rOuter,
							float *tHit,
							Point *rk4Hit,
							float *b1,
							float *b2,
							const bool invert = false,
							const float barycentricEpsilon = 0.03f,
							const float rk4PlaneThreshold = 1e-4f,
							float *finalPlaneDist = nullptr,
							bool *nearBary = nullptr,
							const float uvSeamTolerance = 1e-6f,
							const UVCrossPolicy uvPolicy = UV_REJECT,
							const float side0 = 0.f) {
		
		// Triangle plane setup
		const Vector edge1 = p1 - p0;
		const Vector edge2 = p2 - p0;
		const Vector N = Normalize(Cross(edge1, edge2));

		// Optional INSIGHT straight-line shortcut
		const float r0 = (ray.origin - grinCenter).Length();
		if (r0 < rInner) {
			Ray linearRay(ray.origin, ray.direction, ray.mint, ray.maxt);
			return Intersect(linearRay, p0, p1, p2, tHit, b1, b2);
		}

		// RK4 setup
		const float stepSize = ray.stepSize;
		const int maxSteps = ray.numSteps;

		Point pos = ray.origin;
		Vector dir = ray.direction;

		float tAccum = 0.f;

		// For sign-change detection
		Point prevPos = ray.origin;
		float prevDist = Dot(prevPos - p0, N);
		float currDist = prevDist;

		bool nearEdge = false;

		for (int i = 0; i < maxSteps; ++i) {
			// Compute GRIN curvature at current position
			Vector k1 = ComputeGRINField(pos, ray.beta, ray.gamma, grinCenter, rInner, rOuter, invert);
			Vector k2 = ComputeGRINField(pos + 0.5f * stepSize * k1, ray.beta, ray.gamma, grinCenter, rInner, rOuter, invert);
			Vector k3 = ComputeGRINField(pos + 0.5f * stepSize * k2, ray.beta, ray.gamma, grinCenter, rInner, rOuter, invert);
			Vector k4 = ComputeGRINField(pos + stepSize * k3, ray.beta, ray.gamma, grinCenter, rInner, rOuter, invert);

			// RK4 update
			dir += (stepSize / 6.f) * (k1 + 2.f * k2 + 2.f * k3 + k4);
			pos += stepSize * dir;
			tAccum += stepSize;

			// Current signed distance to plane
			currDist = Dot(pos - p0, N);

			// Crossing detected or near the plane: refine within the step
			if ((prevDist * currDist <= 0.f) || (std::fabs(currDist) < rk4PlaneThreshold)) {
				Point A = prevPos, B = pos;
				float da = prevDist, db = currDist;
				for (int it = 0; it < 5; ++it) {
					const Point M = (A + B) * 0.5f;
					const float dm = Dot(M - p0, N);
					if (dm == 0.f) { A = B = M; break; }
					if (da * dm <= 0.f) { B = M; db = dm; }
					else { A = M; da = dm; }
				}
				const Point planePoint = (A + B) * 0.5f;

				// Precise tHit inside the step
				const float segLen = (pos - prevPos).Length();
				const float hitLen = (planePoint - prevPos).Length();
				const float frac = (segLen > 0.f) ? (hitLen / segLen) : 0.f;
				const float tStepStart = tAccum - stepSize;
				const float tPlane = tStepStart + stepSize * frac;

				// Backface guard inside the RK4 step (reject if we didn't actually cross)
				{
						const float sidePlane = Dot(planePoint - p0, N);
						if (!invert && (side0 * sidePlane > std::max(1e-8f, 0.5f * std::max(barycentricEpsilon, uvSeamTolerance) + 0.25f * rk4PlaneThreshold)))
								return false;
				}

				const BaryResult br = ProjectToTriangleBary(p0, p1, p2, planePoint,
												barycentricEpsilon, uvSeamTolerance);
				if (br.inside) {
					*tHit = tPlane;
					*rk4Hit = planePoint;
					*b1 = br.b1;
					*b2 = br.b2;
					ClampBarycentricSoft(planePoint, p0, p1, p2, 1e-6f, b1, b2);
					return true;
				}
				if (br.nearEdge && uvPolicy == UV_EDGE_PROJECT &&
					EdgeProjectBary(p0, p1, p2, planePoint, br.edgeId, b1, b2)) {
					*tHit = tPlane;
					*rk4Hit = planePoint;
					ClampBarycentricSoft(planePoint, p0, p1, p2, 1e-6f, b1, b2);
					return true;
				}
				if (br.nearEdge)
					nearEdge = true;
			}

			prevPos = pos;
			prevDist = currDist;
		}

		if (finalPlaneDist)
			*finalPlaneDist = std::fabs(currDist);
		if (nearBary)
			*nearBary = nearEdge;

		return false; // No intersection found
	}

	// 🔥GRIN Curved Path Additions 
        static bool xPRIMEIntersect(
						const xPRIMEray &ray,
						const Point &p0,
						const Point &p1,
						const Point &p2,
						const Point &grinCenter,
						const float rInner,
						const float rOuter,
						float *tHit,
						float *b1,
						float *b2,
						const bool invert = false,
						const float insightCurvatureThreshold = 1e-6f,
						const float barycentricEpsilon = 0.03f,
						const float rk4PlaneThreshold = 1e-4f,
						const float uvSeamTolerance = 1e-6f,
						const UVCrossPolicy uvPolicy = UV_REJECT,
						float *finalPlaneDist = nullptr,
						bool *nearBary = nullptr,
						const float stitchBaryMargin = 0.02f,
						Point *outHitPos = nullptr) {

		// Compute plane normal
		const Vector edge1 = p1 - p0;
		const Vector edge2 = p2 - p0;
		const Vector N = Normalize(Cross(edge1, edge2));
		const float side0 = Dot(ray.origin - p0, N);

		// Adaptive-ish guard tolerance (keeps backfaces out but allows tiny slack)
		const float effBaryEps = std::max(barycentricEpsilon, uvSeamTolerance);
		const float sideTol = std::max(1e-8f, 0.5f * effBaryEps + 0.25f * rk4PlaneThreshold);

		Point approxHit;
		float tINSIGHT;

		// STEP 1: Do INSIGHT symbolic intersection with the triangle's plane
		if (!IntersectINSIGHT(ray, p0, N, &tINSIGHT, &approxHit, insightCurvatureThreshold))
			return false;

		// Early backface guard at the symbolic plane hit (unless invert)
		const float sideApprox = Dot(approxHit - p0, N);
		if (!invert && (side0 * sideApprox > sideTol))
			return false;

		// STEP 2: Quick barycentric test
		const BaryResult br0 = ProjectToTriangleBary(p0, p1, p2, approxHit,
						barycentricEpsilon, uvSeamTolerance);
		if (!br0.inside) {
			if (br0.nearEdge && uvPolicy == UV_EDGE_PROJECT) {
				EdgeProjectBary(p0, p1, p2, approxHit, br0.edgeId, b1, b2);
			} else {
				if (nearBary)
						*nearBary = br0.nearEdge ||
								ProjectToTriangleBary(p0, p1, p2, approxHit,
										barycentricEpsilon + stitchBaryMargin,
										uvSeamTolerance).nearEdge;
				return false;
			}
		} else if (nearBary)
			*nearBary = false;
		
		// STEP 3: RK4 refinement through GRIN field (like your Blender)
		Point rk4Hit;
		float tRK4;

		if (!RK4_GRINIntersect(ray, p0, p1, p2, grinCenter, rInner, rOuter,
										&tRK4, &rk4Hit, b1, b2, invert,
										barycentricEpsilon, rk4PlaneThreshold,
										finalPlaneDist, nearBary,
										uvSeamTolerance, uvPolicy,
										side0))
			return false;

		// Clamp barycentrics softly to the triangle interior before returning
		ClampBarycentricSoft(rk4Hit, p0, p1, p2, /*tol=*/0.f, b1, b2);
		// Same-side guard: reject hits that didn't cross the plane
		const float sideH = Dot(rk4Hit - p0, N);
		if (!invert && (side0 * sideH > sideTol)) {
			if (nearBary) *nearBary = false;
			return false;
		}

		if (outHitPos)
			*outHitPos = rk4Hit;

		*tHit = tRK4;

		// Final guard: tiny numerical excursions can still happen → clamp softly
		ClampBarycentricSoft(rk4Hit, p0, p1, p2, 1e-6f, b1, b2);
		if (nearBary)
			*nearBary = false;
		return true;
	}

	static bool Intersect(const Ray &ray, const Point &p0, const Point &p1, const Point &p2,
		float *t, float *b1, float *b2) {
		const Vector e1 = p1 - p0;
		const Vector e2 = p2 - p0;
		const Vector s1 = Cross(ray.d, e2);

		const float divisor = Dot(s1, e1);
		if (divisor == 0.f)
			return false;

		const float invDivisor = 1.f / divisor;

		// Compute first barycentric coordinate
		const Vector d = ray.o - p0;
		*b1 = Dot(d, s1) * invDivisor;
		if (*b1 < 0.f)
			return false;

		// Compute second barycentric coordinate
		const Vector s2 = Cross(d, e1);
		*b2 = Dot(ray.d, s2) * invDivisor;
		if (*b2 < 0.f)
			return false;

		const float b0 = 1.f - *b1 - *b2;
		if (b0 < 0.f)
			return false;

		// Compute _t_ to intersection point
		*t = Dot(e2, s2) * invDivisor;
		if (*t < ray.mint || *t > ray.maxt)
			return false;

		return true;
	}

	bool Intersect(const Ray &ray, const Point *verts, float *t, float *b1, float *b2) const {
		const Point &p0 = verts[v[0]];
		const Point &p1 = verts[v[1]];
		const Point &p2 = verts[v[2]];

		return Intersect(ray, p0, p1, p2, t, b1, b2);
	}

	float Area(const Point *verts) const {
		const Point &p0 = verts[v[0]];
		const Point &p1 = verts[v[1]];
		const Point &p2 = verts[v[2]];

		return Area(p0, p1, p2);
	}

	Normal GetGeometryNormal(const Point *verts) const {
		const Point &p0 = verts[v[0]];
		const Point &p1 = verts[v[1]];
		const Point &p2 = verts[v[2]];

		return Normal(Normalize(Cross(p1 - p0, p2 - p0)));
	}

	void Sample(const Point *verts, const float u0,
		const float u1, Point *p, float *b0, float *b1, float *b2) const {
		// Old triangle uniform sampling
		// UniformSampleTriangle(u0, u1, b0, b1);

		// This new implementation samples from a one dimensional sample
		LowDiscrepancySampleTriangle(u0, b0, b1);

		// Get triangle vertices in _p1_, _p2_, and _p3_
		const Point &p0 = verts[v[0]];
		const Point &p1 = verts[v[1]];
		const Point &p2 = verts[v[2]];
		*b2 = 1.f - (*b0) - (*b1);
		*p = (*b0) * p0 + (*b1) * p1 + (*b2) * p2;
	}

	static float Area(const Point &p0, const Point &p1, const Point &p2) {
		return .5f * Cross(p1 - p0, p2 - p0).Length();
	}

	bool GetBaryCoords(const Point *verts, const Point &hitPoint, float *b1, float *b2) const {
		const Point &p0 = verts[v[0]];
		const Point &p1 = verts[v[1]];
		const Point &p2 = verts[v[2]];

		return GetBaryCoords(p0, p1, p2, hitPoint, b1, b2);
	}

	static bool GetBaryCoords(const Point &p0, const Point &p1, const Point &p2,
			const Point &hitPoint, float *b1, float *b2) {
		const Vector u = p1 - p0;
		const Vector v = p2 - p0;
		const Vector w = hitPoint - p0;

		const Vector vCrossW = Cross(v, w);
		const Vector vCrossU = Cross(v, u);

		if (Dot(vCrossW, vCrossU) < 0.f)
			return false;

		const Vector uCrossW = Cross(u, w);
		const Vector uCrossV = Cross(u, v);

		if (Dot(uCrossW, uCrossV) < 0.f)
			return false;

		const float denom = uCrossV.Length();
		const float r = vCrossW.Length() / denom;
		const float t = uCrossW.Length() / denom;
		
		*b1 = r;
		*b2 = t;

		return ((r <= 1.f) && (t <= 1.f) && (r + t <= 1.f));
	}

	bool GetBaryCoordsSoft(const Point *verts, const Point &hitPoint,
			float *b1, float *b2, const float epsilon = 0.01f) const {
		const Point &p0 = verts[v[0]];
		const Point &p1 = verts[v[1]];
		const Point &p2 = verts[v[2]];

		return GetBaryCoordsSoft(p0, p1, p2, hitPoint, b1, b2, epsilon);
	}

	static bool GetBaryCoordsSoft(const Point &p0, const Point &p1, const Point &p2,
			const Point &hitPoint, float *b1, float *b2, const float epsilon = 0.1f) {
		const Vector u = p1 - p0;
		const Vector v = p2 - p0;
		const Vector w = hitPoint - p0;

		const Vector vCrossW = Cross(v, w);
		const Vector vCrossU = Cross(v, u);

		if (Dot(vCrossW, vCrossU) < 0.f)
			return false;

		const Vector uCrossW = Cross(u, w);
		const Vector uCrossV = Cross(u, v);

		if (Dot(uCrossW, uCrossV) < 0.f)
			return false;

		const float denom = uCrossV.Length();
		const float r = vCrossW.Length() / denom;
		const float t = uCrossW.Length() / denom;

		*b1 = r;
		*b2 = t;

		// Epsilon-softened bounds
		return ((r >= -epsilon) && (t >= -epsilon) && ((r + t) <= (1.f + epsilon)));
	}

	
	static float GetHeight(const float a, const float b, const float c) {
		// Heron's formula for triangle area
		const float s = (a + b + c) * .5f;
		const float area = sqrtf(s * (s - a) * (s - b) * (s - c));

		// h = (A / a) * 2
		return (area / a) * 2.f;
	}


	unsigned int v[3];

	friend class boost::serialization::access;

private:
	template<class Archive>	void serialize(Archive & ar, const unsigned int version) {
		ar & v[0];
		ar & v[1];
		ar & v[2];
	}
};

inline std::ostream & operator<<(std::ostream &os, const Triangle &tri) {
	os << "Triangle[" << tri.v[0] << ", " << tri.v[1] << ", " << tri.v[2] << "]";
	return os;
}

}

// Eliminate serialization overhead at the cost of
// never being able to increase the version.
BOOST_CLASS_IMPLEMENTATION(luxrays::Triangle, boost::serialization::object_serializable)
BOOST_CLASS_EXPORT_KEY(luxrays::Triangle)

#endif	/* _LUXRAYS_TRIANGLE_H */
