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

#ifndef OPTICAL_FLOW_INCLUDED
#define OPTICAL_FLOW_INCLUDED

#include "Misha/SphericalGeometry.h"
#include "Misha/MoebiusCentering.h"

namespace SphericalGeometry
{
	template< typename Real , bool DivergenceFree >
	struct OpticalFlow
	{
		struct Stats { Real l1Norm , l2Norm , l1Difference , l2Difference; };

		struct Parameters
		{
			unsigned int resolution;
			unsigned int subSteps;
			unsigned int integrationSamples;
			Real stepSize;
			Real smoothingWeight;
			Real epsilon;
			Real bandLimitDampening;
			bool useSemiImplicit;
			bool symmetric;
			bool rescaleFlow;
			Parameters( void )
				: resolution(256) , subSteps(1) , integrationSamples(1024) , stepSize( (Real)0.005 ) , smoothingWeight( (Real)0.05 ) , epsilon( (Real)1e-16 ) , useSemiImplicit(true) , symmetric(true) , rescaleFlow(false) , bandLimitDampening( (Real)1 )
			{}
			Real scalarSigma( void ) const { return - (Real)log(bandLimitDampening) / ( resolution/2 * ( resolution/2 + 1 ) ); }
#if 1
#pragma message( "[WARNING] Forcing smoother flow field" )
			Real flowFieldSigma( void ) const { return - (Real)log(bandLimitDampening/resolution) / ( resolution/2 * ( resolution/2 + 1 ) ); }
#else
			Real flowFieldSigma( void ) const { return - (Real)log(bandLimitDampening) / ( resolution/2 * ( resolution/2 + 1 ) ); }
#endif
		};
		OpticalFlow( SphericalGeometry::Mesh< Real > &source , SphericalGeometry::Mesh< Real > &target , const std::vector< Real > &sourceSignal , const std::vector< Real > &targetSignal , Parameters params );
		~OpticalFlow( void );
		void advance( void );
		Stats getStats( void ) const;
		Point3D< Real > flowField( Point3D< Real > p ) const;
	protected:
		typedef EigenSolverCholeskyLDLt< Real > Solver;
		typedef typename SphericalGeometry::RegularGrid< Real >::template PrimalSignal< Point3D< Real > > SphericalFlowField;
		typedef typename SphericalGeometry::RegularGrid< Real >::template PrimalSignal< Real > SphericalSignal;
		typedef typename SphericalGeometry::RegularGrid< Real >::template PrimalVectorField< true , false > SphericalGradientField;
		typedef typename SphericalGeometry::RegularGrid< Real >::template PrimalVectorField< !DivergenceFree , true > FlowField;


		SphericalGeometry::Mesh< Real > &_source , &_target;
		const std::vector< Real > &_sourceSignal , &_targetSignal;
		Parameters _params;

		SphericalGrid< Real > _sGrid , _sGrids[2];
		SphericalSignal _signals[2] , _dSignal;
		SphericalGradientField _gradients[2] , _sGradients;
		SphericalFlowField _flowField;
		HarmonicTransform< Real > _hForm; 
		FourierKeyS2< Real > _key;
		SparseMatrix< Real , int > _M , _FlowFieldMass , _MassAndStiffness;
		std::vector< Real > _b;
		FlowField _x;
		SphericalSignal __x;

		Solver *_solver;

		void _rasterizeSignals( void );
		void _setFlowField( void );
		void _advect( void );

		template< typename Data , bool _DivergenceFree=DivergenceFree >
		typename std::enable_if< _DivergenceFree >::type _setWeights( typename FlowField::template MassAndStiffnessOperator< Data > &massAndStiffness , Real mGradientWeight , Real mJGradientWeight , Real sWeight ) const
		{
			massAndStiffness.mWeightJGradient() = mJGradientWeight;
			massAndStiffness.sWeightJGradient() = sWeight;
		}
		template<  typename Data ,bool _DivergenceFree=DivergenceFree >
		typename std::enable_if< !_DivergenceFree >::type _setWeights( typename FlowField::template MassAndStiffnessOperator< Data > &massAndStiffness , Real mGradientWeight , Real mJGradientWeight , Real sWeight ) const
		{
			massAndStiffness.template mWeightGradient< !DivergenceFree , true >() = mGradientWeight;
			massAndStiffness.mWeightJGradient() = mJGradientWeight;
			massAndStiffness.template sWeightGradient< !DivergenceFree , true >() = sWeight;
			massAndStiffness.sWeightJGradient() = sWeight;
		}

	};
}
#include "Misha/OpticalFlow.inl"
#endif // OPTICAL_FLOW_INCLUDED