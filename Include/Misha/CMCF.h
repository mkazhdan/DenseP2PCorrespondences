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

#ifndef CMCF_INCLUDED
#define CMCF_INCLUDED

#include "Misha/SphericalGeometry.h"
#include "Misha/SparseMatrix.h"
#include "Misha/Solver.h"

namespace SphericalGeometry
{
	template< typename Real >
	struct CMCF
	{
		struct Stats{ Real deformationScale , quasiConformalRatio , radialDeviation; };

		struct Parameters
		{
			Real stepSize;
			bool lump;
			Parameters( void )
				: stepSize( (Real)0.1 ) , lump(false)
			{}
		};
		CMCF( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , SphericalGeometry::Mesh< Real > &sMesh , Parameters params );
		~CMCF( void );
		void advance( void );

		Stats getStats( void ) const;

	protected:
		typedef EigenSolverCholeskyLLt< Real > Solver;

		static const Real _QuasiConformalCutOff;
		SphericalGeometry::Mesh< Real > &_sMesh;
		Parameters _params;
		SparseMatrix< Real , int > _M , _S , _A;
		std::vector< Point3D< Real > > _b , _x;
		std::vector< Real > __b , __x;
		Solver *_solver;
		std::vector< SquareMatrix< Real , 2 > > _massInvs;
		std::vector< Real > _areas;

		void _normalize( void );
		Real _deformationScale( void ) const;
		Real _radialDeviation( void ) const;
		Real _quasiConformalRatio( void ) const;

		static SquareMatrix< Real , 2 > __TriangleMassMatrix( const Point3D< Real > v[] );
		static SquareMatrix< Real , 3 > _TriangleMassMatrix( const Point3D< Real > vertices[] );
		static SquareMatrix< Real , 3 > _TriangleStiffnessMatrix( const Point3D< Real > vertices[] );
		static void _ReorderMatrixEntries( SparseMatrix< Real , int > &M );
		template< typename TriangleMatrixFunctor > SparseMatrix< Real , int > _systemMatrix( TriangleMatrixFunctor F ) const;
		SparseMatrix< Real , int > _massMatrix( void ) const { return _systemMatrix( _TriangleMassMatrix ); }
		SparseMatrix< Real , int > _stiffnessMatrix( void ) const { return _systemMatrix( _TriangleStiffnessMatrix ); }
	};
}
#include "Misha/CMCF.inl"
#endif // CMCF_INCLUDED