#ifndef _SLG_STITCHHINT_H
#define _SLG_STITCHHINT_H

#include <limits>
#include "luxrays/luxrays.h"

namespace slg {

// Hint describing the best near-miss candidate during curved-ray traversal
struct StitchHint {
    u_int meshIndex = 0xffffffffu;
    u_int triIndex  = 0xffffffffu;
    float finalPlaneDist = std::numeric_limits<float>::infinity();
    bool nearBary = false;

    bool valid() const {
        return meshIndex != 0xffffffffu && triIndex != 0xffffffffu;
    }

    void consider(u_int mesh, u_int tri, float planeDist, bool nearB) {
        if (planeDist < finalPlaneDist || (nearB && !nearBary)) {
            meshIndex = mesh;
            triIndex = tri;
            finalPlaneDist = planeDist;
            nearBary = nearB;
        }
    }
};

} // namespace slg

#endif // _SLG_STITCHHINT_H