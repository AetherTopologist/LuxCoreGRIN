#line 2 "grin_uv_funcs.cl"

#ifndef LUXCORE_GRIN_USE_GRAM_PROJECTION
#define LUXCORE_GRIN_USE_GRAM_PROJECTION 1
#endif

OPENCL_FORCE_INLINE float3 Triangle_ComputeGRINField(
    const float3 pos,
    const float beta,
    const float3 gamma,
    const float3 center,
    const float rInner,
    const float rOuter,
    const bool invert) {
    const float3 offset = pos - center;
    const float r = length(offset);

    if (r < rInner)
        return MAKE_FLOAT3(0.f, 0.f, 0.f);

    const float t = clamp((r - rInner) / (rOuter - rInner), 0.f, 1.f);
    const float t_x = pow(t, gamma.x);
    const float t_y = pow(t, gamma.y);
    const float t_z = pow(t, gamma.z);

    const float effectiveBeta = invert ? -beta : beta;

    return MAKE_FLOAT3(
        effectiveBeta * t_x * offset.x / r,
        effectiveBeta * t_y * offset.y / r,
        effectiveBeta * t_z * offset.z / r);
}

OPENCL_FORCE_INLINE void ProjectTangentToUV_cl(
    const float3 tangent, const float3 dpdu, const float3 dpdv,
    float *du, float *dv) {
#if LUXCORE_GRIN_USE_GRAM_PROJECTION
    const float a = dot(dpdu, dpdu);
    const float b = dot(dpdu, dpdv);
    const float d = dot(dpdv, dpdv);
    const float e = dot(tangent, dpdu);
    const float f = dot(tangent, dpdv);
    const float det = a * d - b * b;

    if (det > 1e-20f) {
        const float invDet = 1.f / det;
        *du = (e * d - f * b) * invDet;
        *dv = (f * a - e * b) * invDet;
    } else {
        *du = (a > 0.f) ? (e / a) : 0.f;
        *dv = (d > 0.f) ? (f / d) : 0.f;
    }
#else
    const float a = dot(dpdu, dpdu);
    const float d = dot(dpdv, dpdv);
    const float e = dot(tangent, dpdu);
    const float f = dot(tangent, dpdv);
    *du = (a > 0.f) ? (e / a) : 0.f;
    *dv = (d > 0.f) ? (f / d) : 0.f;
#endif
    const float mag = hypot(*du, *dv);
    if (mag > 1e3f) {
        const float s = 1e3f / mag;
        *du *= s;
        *dv *= s;
    }
}