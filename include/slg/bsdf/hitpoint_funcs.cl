#line 2 "hitpoint_funcs.cl"

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

#ifndef GRIN_UV_USE_RAW_BARY
#define GRIN_UV_USE_RAW_BARY 1
#endif

// Used when hitting a surface
OPENCL_FORCE_INLINE void HitPoint_Init(__global HitPoint *hitPoint, const bool throughShadowTransp,
                const uint meshIndex, const uint triIndex,
                const float3 pnt, const float3 fixedDir,
                const float b1, const float b2,
				const float passThroughEvnt
				MATERIALS_PARAM_DECL) {
	hitPoint->throughShadowTransparency = throughShadowTransp;
	hitPoint->passThroughEvent = passThroughEvnt;

	VSTORE3F(pnt, &hitPoint->p.x);
	VSTORE3F(fixedDir, &hitPoint->fixedDir.x);

	hitPoint->objectID = sceneObjs[meshIndex].objectID;

	// Interpolate face normal
	const float3 geometryN = ExtMesh_GetGeometryNormal(&hitPoint->localToWorld, meshIndex, triIndex EXTMESH_PARAM);
	VSTORE3F(geometryN,  &hitPoint->geometryN.x);
	const float3 interpolatedN = ExtMesh_GetInterpolateNormal(&hitPoint->localToWorld, meshIndex, triIndex, b1, b2 EXTMESH_PARAM);
	VSTORE3F(interpolatedN,  &hitPoint->interpolatedN.x);
	const float3 shadeN = interpolatedN;
	VSTORE3F(shadeN,  &hitPoint->shadeN.x);

	hitPoint->intoObject = (dot(-fixedDir, geometryN) < 0.f);

	// Interpolate UV coordinates
#if GRIN_UV_USE_RAW_BARY
	float3 localHit = Transform_InvApplyPoint(&hitPoint->localToWorld, pnt);
	float rb0, rb1, rb2;
	{
			__global const ExtMesh* restrict meshDesc = &meshDescs[meshIndex];
			__global const Triangle* restrict tri = &triangles[meshDesc->trisOffset + triIndex];
			__global const Point* restrict verts = &vertices[meshDesc->vertsOffset];
			const float3 p0 = VLOAD3F(&verts[tri->v[0]].x);
			const float3 p1 = VLOAD3F(&verts[tri->v[1]].x);
			const float3 p2 = VLOAD3F(&verts[tri->v[2]].x);
			const float3 v0 = p1 - p0;
			const float3 v1 = p2 - p0;
			const float3 v2 = localHit - p0;
			const float d00 = dot(v0, v0);
			const float d01 = dot(v0, v1);
			const float d11 = dot(v1, v1);
			const float d20 = dot(v2, v0);
			const float d21 = dot(v2, v1);
			const float denom = d00 * d11 - d01 * d01;
			rb1 = (d11 * d20 - d01 * d21) / denom;
			rb2 = (d00 * d21 - d01 * d20) / denom;
			rb0 = 1.f - rb1 - rb2;
	}
	rb0 = clamp(rb0, -1e-6f, 1.f + 1e-6f);
	rb1 = clamp(rb1, -1e-6f, 1.f + 1e-6f);
	rb2 = 1.f - rb1 - rb0;
	const float2 defaultUV = ExtMesh_GetInterpolateUV(meshIndex, triIndex, rb1, rb2, 0 EXTMESH_PARAM);
#else
	const float2 defaultUV = ExtMesh_GetInterpolateUV(meshIndex, triIndex, b1, b2, 0 EXTMESH_PARAM);
#endif
	VSTORE2F(defaultUV, &hitPoint->defaultUV.u);
	// GRIN-UV DISABLED (minimal fast revert): UV mapping mirrors non-GRIN behavior.
	// Reason: adaptive barycentric epsilon + smart stepping fixed geometry; UV distortion is unnecessary and caused apparent "zoom".
	// TODO: re-enable via a runtime toggle if we want to experiment later.
	VSTORE2F(MAKE_FLOAT2(0.f, 0.f), &hitPoint->grinUvDelta.u);

	hitPoint->meshIndex = meshIndex;
	hitPoint->triangleIndex = triIndex;
	hitPoint->triangleBariCoord1 = b1;
	hitPoint->triangleBariCoord2 = b2;

	// Compute geometry differentials
	float3 dndu, dndv, dpdu, dpdv;
	ExtMesh_GetDifferentials(
			&hitPoint->localToWorld,
			meshIndex,
			triIndex,
			shadeN, 0,
			&dpdu, &dpdv,
			&dndu, &dndv
			EXTMESH_PARAM);
	VSTORE3F(dpdu, &hitPoint->dpdu.x);
	VSTORE3F(dpdv, &hitPoint->dpdv.x);
	VSTORE3F(dndu, &hitPoint->dndu.x);
	VSTORE3F(dndv, &hitPoint->dndv.x);
}

// Initialize all fields
OPENCL_FORCE_INLINE void HitPoint_InitDefault(__global HitPoint *hitPoint) {
	hitPoint->meshIndex = NULL_INDEX;

	hitPoint->throughShadowTransparency = false;
	hitPoint->passThroughEvent = 0.f;

	VSTORE3F(MAKE_FLOAT3(0.f, 0.f, 0.f), &hitPoint->p.x);
	VSTORE3F(MAKE_FLOAT3(0.f, 0.f, 0.f), &hitPoint->fixedDir.x);

	hitPoint->objectID = 0;
	
	VSTORE3F(MAKE_FLOAT3(0.f, 0.f, 0.f),  &hitPoint->geometryN.x);
	VSTORE3F(MAKE_FLOAT3(0.f, 0.f, 0.f),  &hitPoint->interpolatedN.x);
	VSTORE3F(MAKE_FLOAT3(0.f, 0.f, 0.f),  &hitPoint->shadeN.x);

	hitPoint->intoObject = true;

	VSTORE2F(MAKE_FLOAT2(0.f, 0.f), &hitPoint->defaultUV.u);
	// GRIN-UV DISABLED (minimal fast revert): UV mapping mirrors non-GRIN behavior.
	// Reason: adaptive barycentric epsilon + smart stepping fixed geometry; UV distortion is unnecessary and caused apparent "zoom".
	// TODO: re-enable via a runtime toggle if we want to experiment later.
	VSTORE2F(MAKE_FLOAT2(0.f, 0.f), &hitPoint->grinUvDelta.u);

	VSTORE3F(MAKE_FLOAT3(0.f, 0.f, 0.f), &hitPoint->dpdu.x);
	VSTORE3F(MAKE_FLOAT3(0.f, 0.f, 0.f), &hitPoint->dpdv.x);
	VSTORE3F(MAKE_FLOAT3(0.f, 0.f, 0.f), &hitPoint->dndu.x);
	VSTORE3F(MAKE_FLOAT3(0.f, 0.f, 0.f), &hitPoint->dndv.x);

	Transform_Init(&hitPoint->localToWorld);

	hitPoint->interiorVolumeIndex = NULL_INDEX;
	hitPoint->exteriorVolumeIndex = NULL_INDEX;
	hitPoint->interiorIorTexIndex = NULL_INDEX;
	hitPoint->exteriorIorTexIndex = NULL_INDEX;
}

OPENCL_FORCE_INLINE void HitPoint_GetFrame(__global const HitPoint *hitPoint, Frame *frame) {
	Frame_Set_Private(frame, VLOAD3F(&hitPoint->dpdu.x), VLOAD3F(&hitPoint->dpdv.x), VLOAD3F(&hitPoint->shadeN.x));
}

OPENCL_FORCE_INLINE float3 HitPoint_GetGeometryN(__global const HitPoint *hitPoint) {
	return (hitPoint->intoObject ? 1.f : -1.f) * VLOAD3F(&hitPoint->geometryN.x);
}

OPENCL_FORCE_INLINE float3 HitPoint_GetInterpolatedN(__global const HitPoint *hitPoint) {
	return (hitPoint->intoObject ? 1.f : -1.f) * VLOAD3F(&hitPoint->interpolatedN.x);
}

OPENCL_FORCE_INLINE float3 HitPoint_GetShadeN(__global const HitPoint *hitPoint) {
	return (hitPoint->intoObject ? 1.f : -1.f) * VLOAD3F(&hitPoint->shadeN.x);
}

OPENCL_FORCE_INLINE float2 HitPoint_GetUV(__global const HitPoint *hitPoint, const uint dataIndex EXTMESH_PARAM_DECL) {
	const uint meshIndex = hitPoint->meshIndex;

	if (meshIndex != NULL_INDEX) {
		return (dataIndex == 0) ?
				// GRIN-UV DISABLED (minimal fast revert): UV mapping mirrors non-GRIN behavior.
				// Reason: adaptive barycentric epsilon + smart stepping fixed geometry; UV distortion is unnecessary and caused apparent "zoom".
				// TODO: re-enable via a runtime toggle if we want to experiment later.
				VLOAD2F(&hitPoint->defaultUV.u) :
				ExtMesh_GetInterpolateUV(meshIndex, hitPoint->triangleIndex, hitPoint->triangleBariCoord1, hitPoint->triangleBariCoord2, dataIndex EXTMESH_PARAM);
	} else
		return MAKE_FLOAT2(0.f, 0.f);
}

OPENCL_FORCE_INLINE float3 HitPoint_GetColor(__global const HitPoint *hitPoint, const uint dataIndex EXTMESH_PARAM_DECL) {
	const uint meshIndex = hitPoint->meshIndex;

	if (meshIndex != NULL_INDEX)
		return ExtMesh_GetInterpolateColor(meshIndex, hitPoint->triangleIndex, hitPoint->triangleBariCoord1, hitPoint->triangleBariCoord2, dataIndex EXTMESH_PARAM);
	else
		return WHITE;
}	

OPENCL_FORCE_INLINE float HitPoint_GetAlpha(__global const HitPoint *hitPoint, const uint dataIndex EXTMESH_PARAM_DECL) {
	const uint meshIndex = hitPoint->meshIndex;

	if (meshIndex != NULL_INDEX)
		return ExtMesh_GetInterpolateAlpha(meshIndex, hitPoint->triangleIndex, hitPoint->triangleBariCoord1, hitPoint->triangleBariCoord2, dataIndex EXTMESH_PARAM);
	else
		return 1.f;
}

OPENCL_FORCE_INLINE float HitPoint_GetVertexAOV(__global const HitPoint *hitPoint, const uint dataIndex EXTMESH_PARAM_DECL) {
	const uint meshIndex = hitPoint->meshIndex;

	if (meshIndex != NULL_INDEX)
		return ExtMesh_GetInterpolateVertexAOV(meshIndex, hitPoint->triangleIndex, hitPoint->triangleBariCoord1, hitPoint->triangleBariCoord2, dataIndex EXTMESH_PARAM);
	else
		return 0.f;
}

OPENCL_FORCE_INLINE float HitPoint_GetTriAOV(__global const HitPoint *hitPoint, const uint dataIndex EXTMESH_PARAM_DECL) {
	const uint meshIndex = hitPoint->meshIndex;

	if (meshIndex != NULL_INDEX)
		return ExtMesh_GetTriAOV(meshIndex, hitPoint->triangleIndex, dataIndex EXTMESH_PARAM);
	else
		return 0.f;
}
