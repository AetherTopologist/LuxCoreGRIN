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

#ifndef _SLG_GRIN_UV_H
#define _SLG_GRIN_UV_H

#include "luxrays/luxrays.h"

#ifndef LUXCORE_GRIN_USE_GRAM_PROJECTION
#define LUXCORE_GRIN_USE_GRAM_PROJECTION 1
#endif

#ifndef SLG_FORCE_INLINE
#define SLG_FORCE_INLINE inline
#endif

namespace slg {

SLG_FORCE_INLINE static void ProjectTangentToUV(
        const Vector &tangent, const Vector &dpdu, const Vector &dpdv,
        float &du, float &dv) {
#if LUXCORE_GRIN_USE_GRAM_PROJECTION
        const float a = Dot(dpdu, dpdu);
        const float b = Dot(dpdu, dpdv);
        const float d = Dot(dpdv, dpdv);
        const float e = Dot(tangent, dpdu);
        const float f = Dot(tangent, dpdv);
        const float det = a * d - b * b;

        if (det > 1e-20f) {
                const float invDet = 1.f / det;
                du = (e * d - f * b) * invDet;
                dv = (f * a - e * b) * invDet;
        } else {
                du = (a > 0.f) ? (e / a) : 0.f;
                dv = (d > 0.f) ? (f / d) : 0.f;
        }
#else
        const float a = Dot(dpdu, dpdu);
        const float d = Dot(dpdv, dpdv);
        const float e = Dot(tangent, dpdu);
        const float f = Dot(tangent, dpdv);
        du = (a > 0.f) ? (e / a) : 0.f;
        dv = (d > 0.f) ? (f / d) : 0.f;
#endif
}

}

#endif /* _SLG_GRIN_UV_H */