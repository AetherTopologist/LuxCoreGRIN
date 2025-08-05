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
		const float gammaMax = std::max(ray.gamma.x, std::max(ray.gamma.y, ray.gamma.z));
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
		const float gammaMax = std::max(ray.gamma.x, std::max(ray.gamma.y, ray.gamma.z));
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
				const float rk4PlaneThreshold = 1e-4f) {		
		
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
		float prevDist = Dot(pos - p0, N);

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
			float currDist = Dot(pos - p0, N);

			// 🔥 Two tests: crossing the plane OR direct near-plane hit
			if ((prevDist * currDist < 0.f) || (std::fabs(currDist) < rk4PlaneThreshold)) {
				if (GetBaryCoordsSoft(p0, p1, p2, pos, b1, b2, barycentricEpsilon)) {
					*tHit = tAccum;
					*rk4Hit = pos;
					return true;
				}
			}

			prevDist = currDist;
		}

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
					const float rk4PlaneThreshold = 1e-4f) {

		// Compute plane normal
		const Vector edge1 = p1 - p0;
		const Vector edge2 = p2 - p0;
		const Vector N = Normalize(Cross(edge1, edge2));

		Point approxHit;
		float tINSIGHT;

		// STEP 1: Do INSIGHT symbolic intersection with the triangle's plane
		if (!IntersectINSIGHT(ray, p0, N, &tINSIGHT, &approxHit, insightCurvatureThreshold))
			return false;

		// STEP 2: Quick barycentric test
		if (!GetBaryCoordsSoft(p0, p1, p2, approxHit, b1, b2, barycentricEpsilon))
			return false;

		// STEP 3: RK4 refinement through GRIN field (like your Blender)
		Point rk4Hit;
		float tRK4;

		if (!RK4_GRINIntersect(ray, p0, p1, p2, grinCenter, rInner, rOuter,
								&tRK4, &rk4Hit, b1, b2, invert,
								barycentricEpsilon, rk4PlaneThreshold))
			return false;

		*tHit = tRK4;
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
