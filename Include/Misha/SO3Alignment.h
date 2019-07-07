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

#ifndef SO3_ALIGNMENT_INCLUDED
#define SO3_ALIGNMENT_INCLUDED

#include "Misha/SphericalGeometry.h"

namespace SphericalGeometry
{
	template< typename Real >
	struct SO3Alignment
	{
		struct Parameters
		{
			unsigned int resolution;
			unsigned int maxRotations;
			Real bandLimitDampening;
			Real rotationSeparation;
			Real correlationFraction;
			Real epsilon;
			bool removeDC;
			Parameters( void )
				: resolution(128) , bandLimitDampening( (Real)1 ) , removeDC(true) , rotationSeparation( (Real)7.9 ) , correlationFraction( (Real)7.5 ) , maxRotations(1) , epsilon( (Real)1e-16 )
			{}
			Real sigma( void ) const { return - (Real)log(bandLimitDampening) / ( resolution/2 * ( resolution/2 + 1 ) ); }
		};

		static void SetCorrelation( unsigned int gridNum , SphericalGrid< Real > source[] , SphericalGrid< Real > target[] , Parameters parameters , RotationGrid< Real > &correlation );

		static unsigned int TopNRotations( const RotationGrid< Real > &correlation , Parameters parameters , std::vector< Real > &correlationValues , std::vector< SquareMatrix< Real , 3 > > &rotations );
		static unsigned int TopNRotations( unsigned int gridNum , SphericalGrid< Real > source[] , SphericalGrid< Real > target[] , Parameters parameters , std::vector< Real > &correlationValues , std::vector< SquareMatrix< Real , 3 > > &rotations );

		static void AlignVertexValues ( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > sourceSignals[] , const std::vector< Real > targetSignals[] , unsigned int signalNum , Parameters parameters , RotationGrid< Real > &correlation );
		static void AlignFaceIntegrals( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > sourceSignals[] , const std::vector< Real > targetSignals[] , unsigned int signalNum , Parameters parameters , RotationGrid< Real > &correlation );
		static void AlignVertexValues ( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > sourceSignals[] , const std::vector< Real > targetSignals[] , unsigned int signalNum , Parameters parameters , std::vector< SquareMatrix< Real , 3 > > &rotations , std::vector< Real > &correlationValues );
		static void AlignFaceIntegrals( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > sourceSignals[] , const std::vector< Real > targetSignals[] , unsigned int signalNum , Parameters parameters , std::vector< SquareMatrix< Real , 3 > > &rotations , std::vector< Real > &correlationValues );

		static SquareMatrix< Real , 3 > AlignVertexValues( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > sourceSignals[] , const std::vector< Real > targetSignals , unsigned int signalNum , Parameters parameters )
		{
			parameters.maxRotation = 1;
			std::vector< SquareMatrix< Real , 3 > > rotations;
			std::vector< Real > correlationValues;
			AlignVertexValues( source , target , sourceSignals , targetSignals , signalNum , parameters , rotations , correlationValues );
			return rotations[0];
		}
		static SquareMatrix< Real , 3 > AlignFaceIntegrals( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > sourceSignals[] , const std::vector< Real > targetSignals , unsigned int signalNum , Parameters parameters )
		{
			parameters.maxRotation = 1;
			std::vector< SquareMatrix< Real , 3 > > rotations;
			std::vector< Real > correlationValues;
			AlignFaceIntegrals( source , target , sourceSignals , targetSignals , signalNum , parameters , rotations , correlationValues );
			return rotations[0];
		}

		static void AlignVertexValues ( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > &sourceSignal , const std::vector< Real > &targetSignal , Parameters parameters , std::vector< SquareMatrix< Real , 3 > > &rotations , std::vector< Real > &correlationValues ){ return AlignVertexValues ( source , target , &sourceSignal , &targetSignal , 1 , parameters , rotations , correlationValues ); }
		static void AlignFaceIntegrals( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > &sourceSignal , const std::vector< Real > &targetSignal , Parameters parameters , std::vector< SquareMatrix< Real , 3 > > &rotations , std::vector< Real > &correlationValues ){ return AlignFaceIntegrals( source , target , &sourceSignal , &targetSignal , 1 , parameters , rotations , correlationValues ); }

		static SquareMatrix< Real , 3 > AlignVertexValues ( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > &sourceSignal , const std::vector< Real > &targetSignal , Parameters parameters ){ return AlignVertexValues ( source , target , &sourceSignal , &targetSignal , 1 , parameters ); }
		static SquareMatrix< Real , 3 > AlignFaceIntegrals( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > &sourceSignal , const std::vector< Real > &targetSignal , Parameters parameters ){ return AlignFaceIntegrals( source , target , &sourceSignal , &targetSignal , 1 , parameters ); }
	};
}
#include "Misha/SO3Alignment.inl"
#endif // SO3_ALIGNMENT_INCLUDED