#ifndef _LUXRAYS_XPRIMERAY_H
#define _LUXRAYS_XPRIMERAY_H

#include "luxrays/core/geometry/point.h"
#include "luxrays/core/geometry/vector.h"

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

	// Constructor
	xPRIMEray(const Point &o, const Vector &d,
	          const Point &c, float b,
	          const Vector &g,
	          xPRIMErayType t = xPRIMErayType::POWER,
	          float minT = 0.0001f, float maxT = 1e30f)
		: origin(o), direction(Normalize(d)), center(c), beta(b), gamma(g),
		  type(t), mint(minT), maxt(maxT) { }

	// Default constructor
	xPRIMEray()
		: origin(Point()), direction(Vector(0, 0, 1)),
		  center(Point()), beta(1.f),
		  gamma(Vector(1.f, 1.f, 1.f)),
		  type(xPRIMErayType::POWER), mint(0.0001f), maxt(1e30f) { }

};

} // namespace luxrays

#endif	/* _LUXRAYS_XPRIMERAY_H */
