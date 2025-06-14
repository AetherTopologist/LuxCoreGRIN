#ifndef _LUXRAYS_XPRIMERAYCONTEXT_H
#define _LUXRAYS_XPRIMERAYCONTEXT_H

#include "luxrays/core/geometry/vector.h"
#include "luxrays/core/geometry/point.h"
#include "luxrays/core/geometry/xprimeray.h"

namespace luxrays {

// Context object for passing per-scene GRIN field parameters to Intersect()
class xPRIMErayContext {
public:
	bool enabled;
	Point fieldCenter;
	float beta;
	Vector gamma;
	xPRIMErayType type;

	xPRIMErayContext()
		: enabled(false), fieldCenter(Point(0, 0, 0)),
		  beta(2.f), gamma(Vector(1.f, 1.f, 1.6f)),
		  type(xPRIMErayType::POWER) { }

	xPRIMErayContext(const Point &c, float b, const Vector &g, xPRIMErayType t)
		: enabled(true), fieldCenter(c), beta(b), gamma(g), type(t) { }

	// Optional: convert to actual xPRIMEray
	xPRIMEray GenerateRay(const Point &origin, const Vector &direction) const {
		return xPRIMEray(origin, direction, fieldCenter, beta, gamma, type);
	}
};

}

#endif // _LUXRAYS_XPRIMERAYCONTEXT_H
