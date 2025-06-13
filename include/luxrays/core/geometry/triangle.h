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

	
	// 🔥[GRIN] Curved Path Additions
	static bool xPRIMEIntersect(
		const xPRIMEray &ray,
		const Point &p0,
		const Point &p1,
		const Point &p2,
		float *tHit,
		float *b1,
		float *b2) {

		// Step 1: Triangle plane definition
		const Vector edge1 = p1 - p0;
		const Vector edge2 = p2 - p0;
		//const Normal N = Normalize(Cross(edge1, edge2));
		const Vector N = Normalize(Cross(edge1, edge2));

		// Step 2: Solve intersection with curved ray path R(t)
		// R(t) = O + beta * D * t^gamma  ⇒ we solve Dot(R(t) - p0, N) = 0
		// Let R(t) = origin + beta * vec(t^gamma.x, t^gamma.y, t^gamma.z)

		const Point &O = ray.origin;
		const Vector &D = ray.direction;
		const float beta = ray.beta;
		const Vector &gamma = ray.gamma;

		// Build vector-valued power-law ray term
		auto CurveAt = [&](float t) -> Point {
			return O + beta * Vector(
				D.x * std::pow(t, gamma.x),
				D.y * std::pow(t, gamma.y),
				D.z * std::pow(t, gamma.z)
			);
		};

		// Step 3: Numerical root-finding: Solve f(t) = Dot(R(t) - p0, N) = 0
		float t = 1e-4f; // initial guess
		const float tMax = 100.f;
		const int maxIter = 50;
		const float epsilon = 1e-5f;

		bool solved = false;
		for (int i = 0; i < maxIter; ++i) {
			const Point Rt = CurveAt(t);
			const Vector diff = Rt - p0;
			const float f = Dot(diff, N);

			// Derivative estimate via finite difference
			const float dt = t * 0.001f + 1e-6f;
			const Point Rt2 = CurveAt(t + dt);
			const float f2 = Dot(Rt2 - p0, N);
			const float dfdt = (f2 - f) / dt;

			if (std::fabs(f) < epsilon) {
				solved = true;
				break;
			}

			if (dfdt == 0.f)
				break;

			t = t - f / dfdt;
			if (t < 0.f || t > tMax)
				break;
		}

		if (!solved)
			return false;

		// Step 4: Compute hit point and barycentric coordinates
		const Point hitPoint = CurveAt(t);
		const Vector u = p1 - p0;
		const Vector v = p2 - p0;
		const Vector w = hitPoint - p0;

		const float denom = Dot(u, u) * Dot(v, v) - Dot(u, v) * Dot(u, v);
		if (denom == 0.f)
			return false;

		const float s = (Dot(u, u) * Dot(w, v) - Dot(u, v) * Dot(w, u)) / denom;
		const float t_param = (Dot(v, v) * Dot(w, u) - Dot(u, v) * Dot(w, v)) / denom;

		if (s < 0.f || t_param < 0.f || (s + t_param) > 1.f)
			return false;

		*b1 = s;
		*b2 = t_param;
		*tHit = t;
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

	bool GetBaryCoords(const Point *verts, const Point &hitPoint, float *b1, float *b2) const {
		const Point &p0 = verts[v[0]];
		const Point &p1 = verts[v[1]];
		const Point &p2 = verts[v[2]];

		return GetBaryCoords(p0, p1, p2, hitPoint, b1, b2);
	}

	static float Area(const Point &p0, const Point &p1, const Point &p2) {
		return .5f * Cross(p1 - p0, p2 - p0).Length();
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
