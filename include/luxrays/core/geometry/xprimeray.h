// 🔥GRIN XPRIMERAY Header
#ifndef _LUXRAYS_XPRIMERAY_H
#define _LUXRAYS_XPRIMERAY_H

#include "luxrays/core/geometry/point.h"
#include "luxrays/core/geometry/vector.h"
#include "luxrays/core/geometry/bbox.h"
#include "luxrays/core/geometry/triangle.h"

namespace luxrays {

// Enumeration for different symbolic curvature types (power, exp, etc.)
enum class xPRIMErayType {
	POWER,
	LOG,
	EXP,
	SIN,
	LINEAR,     // optional fallback
	CUSTOM      // future symbolic override
};


// GRIN-style curvature-aware ray definition
class xPRIMEray {
public:
	Point origin;        // Starting point of ray
	Vector direction;    // Initial direction (should be normalized)
	Point center;        // Curvature center (origin of field)
	float beta;          // Scalar modulation (GRIN field intensity)
	Vector gamma;        // Power or exponential scaling per axis
	xPRIMErayType type;  // Curve type

	// Optional t-range for validity (for compatibility with Ray)
	float mint, maxt;
	// Integration controls for RK4 path solving
	float stepSize;
	int numSteps;


	// Constructor
	xPRIMEray(const Point &o, const Vector &d,
	          const Point &c, float b,
	          const Vector &g,
	          xPRIMErayType t = xPRIMErayType::POWER,
	          float minT = 0.0001f, float maxT = 1e30f,
			  float step = 0.01f, int steps = 64)
              : origin(o), direction(Normalize(d)), center(c), beta(b), gamma(g),
                type(t), mint(minT), maxt(maxT), stepSize(step), numSteps(steps) { }

	// Default constructor
	xPRIMEray()
		: origin(Point()), direction(Vector(0, 0, 1)),
			center(Point()), beta(2.f),
			gamma(Vector(1.f, 1.f, 1.6f)),
			type(xPRIMErayType::POWER), mint(0.0001f), maxt(1e30f),
			stepSize(0.01f), numSteps(64) { }

};

//------------------------------------------------------------------------------
// Utility helpers
//------------------------------------------------------------------------------

inline Point EvaluateCurvePoint(const xPRIMEray &ray, const float t) {
    return ray.origin + ray.beta * Vector(
            ray.direction.x * std::pow(t, ray.gamma.x),
            ray.direction.y * std::pow(t, ray.gamma.y),
            ray.direction.z * std::pow(t, ray.gamma.z));
}

inline BBox ComputeSweptBBox(const xPRIMEray &ray, const int steps = 16) {
    BBox bbox(EvaluateCurvePoint(ray, ray.mint));
    const float step = (ray.maxt - ray.mint) / static_cast<float>(steps);
    for (int i = 1; i <= steps; ++i) {
        const float t = ray.mint + step * i;
        bbox = Union(bbox, EvaluateCurvePoint(ray, t));
    }
    return bbox;
}

} // namespace luxrays

#endif	/* _LUXRAYS_XPRIMERAY_H */
