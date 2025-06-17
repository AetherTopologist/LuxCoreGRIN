// 🔥GRIN MATERIAL VOLUME NODE HEADER
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

#ifndef _SLG_GRINVOL_H
#define	_SLG_GRINVOL_H

#include "slg/volumes/volume.h"

namespace slg {

//------------------------------------------------------------------------------
// GRINVolume
//------------------------------------------------------------------------------
class GRINVolume : public Volume {
public:
	GRINVolume(const Texture *iorTex, const Texture *emiTex, const Texture *a,
			const luxrays::Spectrum &minIor, const luxrays::Spectrum &maxIor,
			const luxrays::Vector &stretchVec, const std::string &profileType,
			const float beta, const Vector &gamma);

	virtual float Scatter(const luxrays::Ray &ray, const float u, const bool scatteredStart,
		luxrays::Spectrum *connectionThroughput, luxrays::Spectrum *connectionEmission) const;

	// Material interface

	virtual MaterialType GetType() const { return GRIN_VOL; }
	virtual BSDFEvent GetEventTypes() const { return DIFFUSE | REFLECT; };

	virtual luxrays::Spectrum Albedo(const HitPoint &hitPoint) const;

	virtual luxrays::Spectrum Evaluate(const HitPoint &hitPoint,
		const luxrays::Vector &localLightDir, const luxrays::Vector &localEyeDir, BSDFEvent *event,
		float *directPdfW = NULL, float *reversePdfW = NULL) const;
	virtual luxrays::Spectrum Sample(const HitPoint &hitPoint,
		const luxrays::Vector &localFixedDir, luxrays::Vector *localSampledDir,
		const float u0, const float u1, const float passThroughEvent,
		float *pdfW, BSDFEvent *event) const;
	virtual void Pdf(const HitPoint &hitPoint,
		const luxrays::Vector &localLightDir, const luxrays::Vector &localEyeDir,
		float *directPdfW, float *reversePdfW) const;

	virtual void AddReferencedTextures(boost::unordered_set<const Texture *> &referencedTexs) const;
	virtual void UpdateTextureReferences(const Texture *oldTex, const Texture *newTex);

	virtual luxrays::Properties ToProperties() const;

	const Texture *GetSigmaA() const { return sigmaA; }

	luxrays::Spectrum iorMin;
	luxrays::Spectrum iorMax;
	luxrays::Vector stretch; // e.g., Vector(1.0f, 2.0f, 1.0f)
	std::string profile; // "radial", "axial", "shell", etc.
    float beta;
    luxrays::Vector gamma;

	const luxrays::Spectrum &GetIORMin() const { return iorMin; }
	const luxrays::Spectrum &GetIORMax() const { return iorMax; }
	const luxrays::Vector &GetStretch() const { return stretch; }
	const std::string &GetProfile() const { return profile; }
	//const std::string &GetBeta() const { return ; }
	const float &GetBeta() const { return beta; }
	//const std::string &GetGamma() const { return gamma; }
	const luxrays::Vector &GetGamma() const { return gamma; }


protected:
	virtual luxrays::Spectrum SigmaA(const HitPoint &hitPoint) const;
	virtual luxrays::Spectrum SigmaS(const HitPoint &hitPoint) const;

private:
	const Texture *sigmaA;
};

}

#endif	/* _SLG_GRINVOL_H */
