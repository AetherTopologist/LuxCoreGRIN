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

#include "slg/bsdf/bsdf.h"
#include "slg/scene/scene.h"
#include "slg/materials/glass.h"
#include "luxrays/core/geometry/triangle.h"

using namespace luxrays;
using namespace slg;
using namespace std;

// Used when hitting a surface
void BSDF::Init(const bool fixedFromLight, const bool throughShadowTransparency,
		const Scene &scene, const Ray &ray, const RayHit &rayHit,
		const float passThroughEvent, const PathVolumeInfo *volInfo) {
	// Get the scene object
	sceneObject = scene.objDefs.GetSceneObject(rayHit.meshIndex);

	//SLG_LOG("🔥GRIN [BSDF::Init()] Surface");

	// Get the mesh
	const ExtMesh *mesh = sceneObject->GetExtMesh();
	mesh->GetLocal2World(ray.time, hitPoint.localToWorld);

	Point hp;
	// ------------------------------------------------------------
	// Path Logic Branch: Curved vs. Straight Ray Behavior
	//
	// If ray.rayType == RAYTYPE_CURVED:
	//     - Apply RK4 or symbolic curved ray stepping
	//     - Interpolate hitpoint via barycentric mesh vertices
	//     - Compute GRIN field distortion, non-linear IOR paths
	//
	// Else (RAYTYPE_DEFAULT):
	//     - Standard linear ray marching and hit resolution
	//     - Use ray(rayHit.t) for hitpoint position
	//     - Mesh UVs, normals, and shading all follow original logic
	// ------------------------------------------------------------
	if (ray.IsCurved()) {
		const Triangle *tris = mesh->GetTriangles();
		const Point *verts = mesh->GetVertices();
		const Triangle &tri = tris[rayHit.triangleIndex];
		const float b0 = 1.f - rayHit.b1 - rayHit.b2;
		const Point lp = b0 * verts[tri.v[0]] + rayHit.b1 * verts[tri.v[1]] + rayHit.b2 * verts[tri.v[2]];
		hp = hitPoint.localToWorld * lp;
	} else
		hp = ray(rayHit.t);

	hitPoint.Init(fixedFromLight, throughShadowTransparency,
					scene, rayHit.meshIndex, rayHit.triangleIndex,
					hp, -ray.d,
					rayHit.b1, rayHit.b2,
					passThroughEvent);

	//hitPoint.Init(fixedFromLight, throughShadowTransparency, scene, rayHit.meshIndex, rayHit.triangleIndex, ray(rayHit.t), -ray.d, rayHit.b1, rayHit.b2, passThroughEvent);
	
        // Apply GRIN-based UV distortion if enabled
        // ------------------------------------------------------------
        // Path Logic Branch: Curved vs. Straight Ray Behavior
        //
        // If ray.rayType == RAYTYPE_CURVED:
        //     - Apply RK4 or symbolic curved ray stepping
        //     - Interpolate hitpoint via barycentric mesh vertices
        //     - Compute GRIN field distortion, non-linear IOR paths
        //
        // Else (RAYTYPE_DEFAULT):
        //     - Standard linear ray marching and hit resolution
        //     - Use ray(rayHit.t) for hitpoint position
        //     - Mesh UVs, normals, and shading all follow original logic
        // ------------------------------------------------------------
        if (ray.IsCurved()) {
		const Vector field = Triangle::ComputeGRINField(
								hp,
								scene.worldGrinInfo.beta,
								scene.worldGrinInfo.gamma,
								scene.worldGrinInfo.center,
								scene.worldGrinInfo.rInner,
								scene.worldGrinInfo.rOuter,
								scene.worldGrinInfo.invert);

		// Project the distortion to be tangent to the surface
		const Vector tangent = field - Dot(field, hitPoint.shadeN) * Vector(hitPoint.shadeN);

		const float du = Dot(tangent, hitPoint.dpdu) / hitPoint.dpdu.LengthSquared();
		const float dv = Dot(tangent, hitPoint.dpdv) / hitPoint.dpdv.LengthSquared();

		const float strength = scene.grinUVDistortionStrength;
		hitPoint.grinUvDelta.u = strength * du;
		hitPoint.grinUvDelta.v = strength * dv;
	}

	// Get the material
	material = sceneObject->GetMaterial();

	// Set interior and exterior volumes
	volInfo->SetHitPointVolumes(hitPoint,
			material->GetInteriorVolume(hitPoint, hitPoint.passThroughEvent),
			material->GetExteriorVolume(hitPoint, hitPoint.passThroughEvent),
			scene.defaultWorldVolume);

	// Check if it is a light source
	if (material->IsLightSource())
		triangleLightSource = scene.lightDefs.GetLightSourceByMeshAndTriIndex(rayHit.meshIndex, rayHit.triangleIndex);
	else
		triangleLightSource = NULL;

	// Apply bump or normal mapping
	material->Bump(&hitPoint);

	// Build the local reference system
	frame = hitPoint.GetFrame();
}

// Used when have a point of a surface
void BSDF::Init(const Scene &scene,
		const u_int meshIndex, const u_int triangleIndex,
		const Point &surfacePoint,
		const float surfacePointBary1, const float surfacePointBary2, 
		const float time,
		const float passThroughEvent, const PathVolumeInfo *volInfo) {
	// Get the scene object
	sceneObject = scene.objDefs.GetSceneObject(meshIndex);

	//SLG_LOG("🔥GRIN [BSDF::Init()] Point of Surface");

	// Get the mesh
	const ExtMesh *mesh = sceneObject->GetExtMesh();
	mesh->GetLocal2World(time, hitPoint.localToWorld);

	const Vector fixedDir = Vector(mesh->GetGeometryNormal(hitPoint.localToWorld, triangleIndex));
	hitPoint.Init(false, false,
				scene, meshIndex, triangleIndex,
				surfacePoint, fixedDir,
				surfacePointBary1, surfacePointBary2, passThroughEvent);

	// Apply GRIN-based UV distortion if enabled
	if (scene.worldGrinInfo.enabled) {
		const Vector field = Triangle::ComputeGRINField(
									surfacePoint,
									scene.worldGrinInfo.beta,
									scene.worldGrinInfo.gamma,
									scene.worldGrinInfo.center,
									scene.worldGrinInfo.rInner,
									scene.worldGrinInfo.rOuter,
									scene.worldGrinInfo.invert);
		const Vector tangent = field -
						Dot(field, hitPoint.shadeN) * Vector(hitPoint.shadeN);

		const float du = Dot(tangent, hitPoint.dpdu) / hitPoint.dpdu.LengthSquared();
		const float dv = Dot(tangent, hitPoint.dpdv) / hitPoint.dpdv.LengthSquared();

		const float strength = scene.grinUVDistortionStrength;
		hitPoint.grinUvDelta.u = strength * du;
		hitPoint.grinUvDelta.v = strength * dv;
	}
	
	// Get the material
	material = sceneObject->GetMaterial();

	// Set interior and exterior volumes
	volInfo->SetHitPointVolumes(hitPoint,
			material->GetInteriorVolume(hitPoint, hitPoint.passThroughEvent),
			material->GetExteriorVolume(hitPoint, hitPoint.passThroughEvent),
			scene.defaultWorldVolume);

	// Check if it is a light source
	if (material->IsLightSource())
		triangleLightSource = scene.lightDefs.GetLightSourceByMeshAndTriIndex(meshIndex, triangleIndex);
	else
		triangleLightSource = NULL;

	// Apply bump or normal mapping
	material->Bump(&hitPoint);

	// Build the local reference system
	frame = hitPoint.GetFrame();
}

// Used when hitting a volume scatter point
void BSDF::Init(const bool fixedFromLight, const bool throughShadowTransparency,
		const Scene &scene, const luxrays::Ray &ray,
		const Volume &volume, const float t, const float passThroughEvent) {
	hitPoint.fromLight = fixedFromLight;
	hitPoint.throughShadowTransparency = throughShadowTransparency;
	hitPoint.passThroughEvent = passThroughEvent;

	//SLG_LOG("🔥GRIN [BSDF::Init()] Volume Scatter Point");

	hitPoint.p = ray(t);
	hitPoint.fixedDir = -ray.d;

	sceneObject = NULL;
	material = &volume;

	hitPoint.geometryN = Normal(-ray.d);
	hitPoint.interpolatedN = hitPoint.geometryN;
	hitPoint.shadeN = hitPoint.geometryN;

	hitPoint.intoObject = true;
	hitPoint.interiorVolume = &volume;
	hitPoint.exteriorVolume = &volume;

	triangleLightSource = NULL;

	hitPoint.defaultUV = UV(0.f, 0.f);
	hitPoint.grinUvDelta.u = 0.f;
	hitPoint.grinUvDelta.v = 0.f;

	CoordinateSystem(Vector(hitPoint.shadeN), &hitPoint.dpdu, &hitPoint.dpdv);
	hitPoint.dndu = Normal();
	hitPoint.dndv = Normal();

	hitPoint.mesh = nullptr;
	hitPoint.triangleIndex = NULL_INDEX;
	hitPoint.triangleBariCoord1 = 0.f;
	hitPoint.triangleBariCoord2 = 0.f;

	hitPoint.objectID = NULL_INDEX;

	// Build the local reference system
	frame.SetFromZ(hitPoint.shadeN);
}

void BSDF::MoveHitPoint(const Point &p, const Normal &n) {
	hitPoint.p = p;
	hitPoint.geometryN = n;
	hitPoint.interpolatedN = n;
	hitPoint.shadeN = n;

	Vector x, y;
	CoordinateSystem(Vector(n), &x, &y);
	frame = Frame(x, y, n);
}

bool BSDF::IsAlbedoEndPoint(const AlbedoSpecularSetting albedoSpecularSetting,
		const float albedoSpecularGlossinessThreshold) const {
	const BSDFEvent eventTypes = GetEventTypes();
	if (!IsDelta() && !((eventTypes & GLOSSY) && (GetGlossiness() < albedoSpecularGlossinessThreshold)))
		return true;

	switch (albedoSpecularSetting) {
		case NO_REFLECT_TRANSMIT:
			return true;
		case ONLY_REFLECT:
			return !((eventTypes & REFLECT) && !(eventTypes & TRANSMIT));
		case ONLY_TRANSMIT:
			return !(!(eventTypes & REFLECT) && (eventTypes & TRANSMIT));
		case REFLECT_TRANSMIT:
			return !((eventTypes & REFLECT) || (eventTypes & TRANSMIT));
		default:
			throw runtime_error("Unknown AlbedoSpecularSetting in BSDF::IsAlbedoEndPoint(): " + ToString(albedoSpecularSetting));
	}
}

bool BSDF::IsCameraInvisible() const {
	return (sceneObject) ? sceneObject->IsCameraInvisible() : false;
}

u_int BSDF::GetObjectID() const {
	//SLG_LOG("🔥GRIN [BSDF::GetObjectID())] ");
	return (sceneObject) ? sceneObject->GetID() : std::numeric_limits<u_int>::max();
}

static string MaterialNULLptrName = "NULL pointer";

const string &BSDF::GetMaterialName() const {
	//SLG_LOG("🔥GRIN [BSDF::GetMaterialName())] ");
	if (material)
		return material->GetName();
	else
		return MaterialNULLptrName;
}

Spectrum BSDF::Albedo() const {
	return material->Albedo(hitPoint);
}

Spectrum BSDF::EvaluateTotal() const {
	return material->EvaluateTotal(hitPoint);
}

bool BSDF::HasBakeMap(const BakeMapType type) const {
	return sceneObject && sceneObject->HasBakeMap(type);
}

Spectrum BSDF::GetBakeMapValue() const {
	return sceneObject->GetBakeMapValue(hitPoint.GetUV(sceneObject->GetBakeMapUVIndex()));
}

//------------------------------------------------------------------------------
// "A Microfacet-Based Shadowing Function to Solve the Bump Terminator Problem"
// by Alejandro Conty Estevez, Pascal Lecocq, and Clifford Stein
// http://www.aconty.com/pdf/bump-terminator-nvidia2019.pdf
//------------------------------------------------------------------------------

// Return alpha ^2 parameter from normal divergence
//static float ShadowTerminatorAvoidanceAlpha2(const Normal &Ni, const Normal &Ns) {
//	const float cos_d = Min(fabsf(Dot(Ni , Ns)), 1.f);
//	const float tan2_d = (1.f - cos_d * cos_d) / (cos_d * cos_d);
//
//	return Clamp(.125f * tan2_d, 0.f, 1.f);
//}

// Shadowing factor
//static float ShadowTerminatorAvoidanceFactor(const Normal &Ni, const Normal &Ns,
//		const Vector &lightDir) {
//	const float alpha2 = ShadowTerminatorAvoidanceAlpha2(Ni, Ns);
//	if (alpha2 > 0.f) {
//		const float cos_i = Max(fabsf(Dot(Ni , lightDir)), 1e-6f);
//		const float tan2_i = (1.f - cos_i * cos_i) / (cos_i * cos_i);
//
//		return 2.f / (1.f + sqrtf(1.f + alpha2 * tan2_i));
//	} else
//		return 1.f;
//}

//------------------------------------------------------------------------------
// "Taming the Shadow Terminator"
// by Matt Jen-Yuan Chiang, Yining Karl Li and Brent Burley
// https://www.yiningkarlli.com/projects/shadowterminator.html
//------------------------------------------------------------------------------

static float ShadowTerminatorAvoidanceFactor(const Normal &Ni, const Normal &Ns,
		const Vector &lightDir) {
	const float dotNsLightDir = Dot(Ns, lightDir);
	if (dotNsLightDir <= 0.f)
		return 0.f;

	const float dotNiNs = Dot(Ni, Ns);
	if (dotNiNs <= 0.f)
		return 0.f;

	const float G = Min(1.f, Dot(Ni, lightDir) / (dotNsLightDir * dotNiNs));
	if (G <= 0.f)
		return 0.f;
	
	const float G2 = G * G;
	const float G3 = G2 * G;

	return -G3 + G2 + G;
}

Spectrum BSDF::Evaluate(const Vector &generatedDir,
		BSDFEvent *event, float *directPdfW, float *reversePdfW) const {
	const Vector &eyeDir = hitPoint.fromLight ? generatedDir : hitPoint.fixedDir;
	const Vector &lightDir = hitPoint.fromLight ? hitPoint.fixedDir : generatedDir;

	const float dotLightDirNG = Dot(lightDir, hitPoint.geometryN);
	const float absDotLightDirNG = fabsf(dotLightDirNG);
	const float dotEyeDirNG = Dot(eyeDir, hitPoint.geometryN);
	const float absDotEyeDirNG = fabsf(dotEyeDirNG);
	
	if (!IsVolume()) {
		// These kind of tests make sense only for materials

		// Avoid glancing angles
		if ((absDotLightDirNG < DEFAULT_COS_EPSILON_STATIC) ||
				(absDotEyeDirNG < DEFAULT_COS_EPSILON_STATIC))
			return Spectrum();

		// Check geometry normal and light direction side
		const float sideTestNG = dotEyeDirNG * dotLightDirNG;
		const BSDFEvent matEvents = material->GetEventTypes();
		if (((sideTestNG > 0.f) && !(matEvents & REFLECT)) ||
				((sideTestNG < 0.f) && !(matEvents & TRANSMIT)))
			return Spectrum();

		// Check shading normal and light direction side
		const float sideTestIS = Dot(eyeDir, hitPoint.interpolatedN) * Dot(lightDir, hitPoint.interpolatedN);
		if (((sideTestIS > 0.f) && !(matEvents & REFLECT)) ||
				((sideTestIS < 0.f) && !(matEvents & TRANSMIT)))
			return Spectrum();
	}

	//SLG_LOG("🔥GRIN [BSDF::Evaluate()]");

	const Vector localLightDir = frame.ToLocal(lightDir);
	const Vector localEyeDir = frame.ToLocal(eyeDir);
	Spectrum result = material->Evaluate(hitPoint, localLightDir, localEyeDir,
			event, directPdfW, reversePdfW);
	assert (!result.IsNaN() && !result.IsInf());
	if (result.Black())
		return result;

	if (!IsVolume()) {
		// Shadow terminator artefact avoidance
		if ((*event & REFLECT) &&
				(*event & (DIFFUSE | GLOSSY)) &&
				(hitPoint.shadeN != hitPoint.interpolatedN))
			result *= ShadowTerminatorAvoidanceFactor(hitPoint.GetLandingInterpolatedN(),
					hitPoint.GetLandingShadeN(), lightDir);

		// Adjoint BSDF (not for volumes)
		if (hitPoint.fromLight)
			result *= (absDotEyeDirNG / absDotLightDirNG);
	}

	return result;
}

Spectrum BSDF::ShadowCatcherSample(Vector *sampledDir,
		float *pdfW, float *absCosSampledDir, BSDFEvent *event) const {
	// Just continue to trace the ray
	*sampledDir = -hitPoint.fixedDir;
	*absCosSampledDir = AbsDot(*sampledDir, hitPoint.geometryN);

	*pdfW = 1.f;
	*event = SPECULAR | TRANSMIT;
	const Spectrum result(1.f);

	// Adjoint BSDF
	if (hitPoint.fromLight) {
		const float absDotFixedDirNG = AbsDot(hitPoint.fixedDir, hitPoint.geometryN);
		const float absDotSampledDirNG = AbsDot(*sampledDir, hitPoint.geometryN);
		return result * (absDotSampledDirNG / absDotFixedDirNG);
	} else
		return result;
}

Spectrum BSDF::Sample(Vector *sampledDir,
		const float u0, const float u1,
		float *pdfW, float *absCosSampledDir,
		BSDFEvent *event) const {
	Vector localFixedDir = frame.ToLocal(hitPoint.fixedDir);
	Vector localSampledDir;
	
	Spectrum result = material->Sample(hitPoint,
			localFixedDir, &localSampledDir, u0, u1, hitPoint.passThroughEvent,
			pdfW, event);
	if (result.Black())
		return result;

	//SLG_LOG("🔥GRIN [BSDF::Sample()]");

	*absCosSampledDir = fabsf(CosTheta(localSampledDir));
	*sampledDir = frame.ToWorld(localSampledDir);

	// Shadow terminator artefact avoidance
	if ((*event & REFLECT) &&
			(*event & (DIFFUSE | GLOSSY))
			&& (hitPoint.shadeN != hitPoint.interpolatedN)) {
		const Vector &lightDir = hitPoint.fromLight ? hitPoint.fixedDir : (*sampledDir);

		result *= ShadowTerminatorAvoidanceFactor(hitPoint.GetLandingInterpolatedN(),
				hitPoint.GetLandingShadeN(), lightDir);
	}

	// Adjoint BSDF
	if (hitPoint.fromLight) {
		const float absDotFixedDirNG = AbsDot(hitPoint.fixedDir, hitPoint.geometryN);
		const float absDotSampledDirNG = AbsDot(*sampledDir, hitPoint.geometryN);
		result *= (absDotSampledDirNG / absDotFixedDirNG);
	}

	return result;
}

void BSDF::Pdf(const Vector &sampledDir, float *directPdfW, float *reversePdfW) const {
	const Vector &eyeDir = hitPoint.fromLight ? sampledDir : hitPoint.fixedDir;
	const Vector &lightDir = hitPoint.fromLight ? hitPoint.fixedDir : sampledDir;
	Vector localLightDir = frame.ToLocal(lightDir);
	Vector localEyeDir = frame.ToLocal(eyeDir);

	material->Pdf(hitPoint, localLightDir, localEyeDir, directPdfW, reversePdfW);
	//SLG_LOG("🔥GRIN [BSDF::Pdf()]");
}

Spectrum BSDF::GetPassThroughTransparency(const bool backTracing) const {
	const Vector localFixedDir = frame.ToLocal(hitPoint.fixedDir);
	
	//SLG_LOG("🔥GRIN [BSDF::GetPassThroughTransparency()]");

	return material->GetPassThroughTransparency(hitPoint, localFixedDir,
			hitPoint.passThroughEvent, backTracing);
}

Spectrum BSDF::GetEmittedRadiance(float *directPdfA, float *emissionPdfW) const {
	return triangleLightSource ?
		triangleLightSource->GetRadiance(hitPoint, directPdfA, emissionPdfW) :
		Spectrum();
}

AlbedoSpecularSetting slg::String2AlbedoSpecularSetting(const string &type) {
	if (type == "NO_REFLECT_TRANSMIT")
		return NO_REFLECT_TRANSMIT;
	else if (type == "ONLY_REFLECT")
		return ONLY_REFLECT;
	else if (type == "ONLY_TRANSMIT")
		return ONLY_TRANSMIT;
	else if (type == "REFLECT_TRANSMIT")
		return REFLECT_TRANSMIT;
	else
		throw runtime_error("Unknown albedo specular setting in String2AlbedoSpecularSetting(): " + type);
}

const string slg::AlbedoSpecularSetting2String(const AlbedoSpecularSetting type) {
	switch (type) {
		case NO_REFLECT_TRANSMIT:
			return "NO_REFLECT_TRANSMIT";
		case ONLY_REFLECT:
			return "ONLY_REFLECT";
		case ONLY_TRANSMIT:
			return "ONLY_TRANSMIT";
		case REFLECT_TRANSMIT:
			return "REFLECT_TRANSMIT";
		default:
			throw runtime_error("Unknown albedo specular setting in AlbedoSpecularSetting2String(): " + ToString(type));
	}
}
