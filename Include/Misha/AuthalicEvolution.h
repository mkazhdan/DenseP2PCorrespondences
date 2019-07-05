/*
Copyright (c) 2018, Michael Kazhdan and Sing Chun Lee
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef AUTHALIC_EVOLUTION_INCLUDED
#define AUTHALIC_EVOLUTION_INCLUDED

#include "Misha/SphericalGeometry.h"
#include "Misha/MoebiusCentering.h"

namespace SphericalGeometry
{
	template< typename Real >
	struct AuthalicEvolution
	{
		struct Stats{ Real min , max , avg , dev; };
		struct Parameters
		{
			unsigned int resolution;
			unsigned int subSteps;
			unsigned int centeringIterations;
			Real stepSize;
			Real smoothness;
			Real epsilon;
			bool useGoldenSectionSearch;
			bool useSemiImplicit;
			Parameters( void )
				: resolution(256) , subSteps(1) , centeringIterations(3) , stepSize( (Real)0.005 ) , smoothness( (Real)2.5e-3 ) , epsilon( (Real)1e-16 ) , useGoldenSectionSearch(false) , useSemiImplicit(true)
			{}
		};
		AuthalicEvolution( SphericalGeometry::Mesh< Real > &sMesh , const std::vector< Point3D< Real > > &meshVertices , Parameters params );
		void advance( void );

		Point3D< Real > flowField( Point3D< Real > p ) const;

		Real flippedTriangleFraction( void ) const;
		Stats getStats( std::vector< Real > &logScaleFactors ) const;
		Stats getStats( void ) const;

	protected:
		SphericalGeometry::Mesh< Real > &_sMesh;
		const std::vector< Point3D< Real > > &_meshVertices;
		Parameters _params;
		std::vector< TriangleIndex > _triangles;
		SphericalGrid< Real > _sGrid;
		typename SphericalGeometry::RegularGrid< Real >::template PrimalSignal< Real > _logScaleFactors;
		typename SphericalGeometry::RegularGrid< Real >::template PrimalSignal< Point3D< Real > > _flowField;
		typename SphericalGeometry::RegularGrid< Real >::template DualSignal< Point3D< Real > > _gradientField;
		HarmonicTransform< Real > _hForm; 
		FourierKeyS2< Real > _key;
		Real _areaScale;

		void _advect( void );
		void _setFlowField( void );
		void _setLogScaleFactors( std::vector< Real > &logScaleFactors ) const;
	};
}
#include "Misha/AuthalicEvolution.inl"
#endif // AUTHALIC_EVOLUTION_INCLUDED