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

#ifndef MOEBIUS_CENTERING_INCLUDED
#define MOEBIUS_CENTERING_INCLUDED

#include "Misha/SphericalGeometry.h"

namespace SphericalGeometry
{
	template< typename Real >
	struct MoebiusCentering
	{
		struct Stats{ Point3D< Real > center; };

		struct Parameters
		{
			enum
			{
				PARAMETERS_CENTER_TO_INVERSION ,
				PARAMETERS_GOLDEN_SECTION_SEARCH ,
				PARAMETERS_POINCARE ,
				PARAMETERS_COUNT
			};
			union
			{
				Real gssTolerance;
				Real poincareMaxNorm;
			};
			bool gaussNewton;
			int type;
			Parameters( int pType=PARAMETERS_POINCARE )
				: type(pType) , gaussNewton(true)
			{
				if( type==PARAMETERS_GOLDEN_SECTION_SEARCH ) gssTolerance = (Real)1e-6;
				else if( type==PARAMETERS_POINCARE ) poincareMaxNorm = (Real)2.;
			}
		};

		MoebiusCentering( SphericalGeometry::Mesh< Real > &sMesh , Parameters params );
		void advance( void );

		Stats getStats( void ) const;
	protected:
		SphericalGeometry::Mesh< Real > &_sMesh;
		Parameters _params;
	};
}
#include "Misha/MoebiusCentering.inl"
#endif // MOEBIUS_CENTERING_INCLUDED