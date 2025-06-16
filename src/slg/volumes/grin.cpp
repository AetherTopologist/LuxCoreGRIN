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

#include <cstddef>

#include "slg/volumes/grin.h"
#include "slg/bsdf/bsdf.h"

using namespace std;
using namespace luxrays;
using namespace slg;

//------------------------------------------------------------------------------
// GRINVolume
//------------------------------------------------------------------------------
GRINVolume::GRINVolume(const Texture *iorTex, const Texture *emiTex, const Texture *a,
                       const luxrays::Spectrum &minIor, const luxrays::Spectrum &maxIor,
                       const luxrays::Vector &stretchVec, const std::string &profileType,
					   const float beta, const Vector &gamma)
    : Volume(iorTex, emiTex), iorMin(minIor), iorMax(maxIor),
      stretch(stretchVec), profile(profileType), beta(beta), gamma(gamma) {
    sigmaA = a;
}

Spectrum GRINVolume::SigmaA(const HitPoint &hitPoint) const {
	return sigmaA->GetSpectrumValue(hitPoint).Clamp();
}

Spectrum GRINVolume::SigmaS(const HitPoint &hitPoint) const {
	return Spectrum();
}

float GRINVolume::Scatter(const Ray &ray, const float u,
		const bool scatteredStart, Spectrum *connectionThroughput,
		Spectrum *connectionEmission) const {

	// TODO BB GRIN
	// Convert ray.o to object-local space
    // Evaluate radial distance
    // Get IOR from gradient logic (Lerp from min→max)
    // Apply basic exponential decay for transmittance
    // Apply curvature step for the ray (for future RK4)
	// TODO BB GRIN

	SLG_LOG("🔥GRIN [GRINVolume::Scatter()]");

	// Point where to evaluate the volume
	HitPoint hitPoint;
	hitPoint.Init();
	hitPoint.fixedDir = ray.d;
	hitPoint.p = ray.o;
	hitPoint.geometryN = hitPoint.interpolatedN = hitPoint.shadeN = Normal(-ray.d);
	hitPoint.passThroughEvent = u;
	
	const float distance = ray.maxt - ray.mint;
	Spectrum transmittance(1.f);

	const Spectrum sigma = SigmaT(hitPoint);
	if (!sigma.Black()) {
		const Spectrum tau = (distance * sigma).Clamp();
		transmittance = Exp(-tau);
	}

	// Apply volume transmittance
	*connectionThroughput *= transmittance;

	// Apply volume emission
	if (volumeEmissionTex)
		*connectionEmission += *connectionThroughput * distance * volumeEmissionTex->GetSpectrumValue(hitPoint).Clamp();
	
	return -1.f;
}

Spectrum GRINVolume::Albedo(const HitPoint &hitPoint) const {
	throw runtime_error("Internal error: called GRINVolume::Albedo()");
}

Spectrum GRINVolume::Evaluate(const HitPoint &hitPoint,
		const Vector &localLightDir, const Vector &localEyeDir, BSDFEvent *event,
		float *directPdfW, float *reversePdfW) const {
	throw runtime_error("Internal error: called GRINVolume::Evaluate()");
}

Spectrum GRINVolume::Sample(const HitPoint &hitPoint,
		const Vector &localFixedDir, Vector *localSampledDir,
		const float u0, const float u1, const float passThroughEvent,
		float *pdfW, BSDFEvent *event) const {
	throw runtime_error("Internal error: called GRINVolume::Sample()");
	SLG_LOG("🔥GRIN [GRINVolume::Sample()]");
}

void GRINVolume::Pdf(const HitPoint &hitPoint,
		const Vector &localLightDir, const Vector &localEyeDir,
		float *directPdfW, float *reversePdfW) const {
	throw runtime_error("Internal error: called GRINVolume::Pdf()");
}

void GRINVolume::AddReferencedTextures(boost::unordered_set<const Texture *> &referencedTexs) const {
	Volume::AddReferencedTextures(referencedTexs);

	sigmaA->AddReferencedTextures(referencedTexs);
}

void GRINVolume::UpdateTextureReferences(const Texture *oldTex, const Texture *newTex) {
	Volume::UpdateTextureReferences(oldTex, newTex);

	if (sigmaA == oldTex)
		sigmaA = newTex;
}

Properties GRINVolume::ToProperties() const {
	Properties props;

	const string name = GetName();
	props.Set(Property("scene.volumes." + name + ".type")("grin"));
	props.Set(Property("scene.volumes." + name + ".absorption")(sigmaA->GetSDLValue()));
	props.Set(Volume::ToProperties());
	props.Set(Property("scene.volumes." + name + ".grin.beta")(beta));
	props.Set(Property("scene.volumes." + name + ".grin.gamma")(gamma));

	return props;
}
