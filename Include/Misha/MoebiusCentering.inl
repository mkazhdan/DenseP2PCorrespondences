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

namespace SphericalGeometry
{
	template< typename Real >
	MoebiusCentering< Real >::MoebiusCentering( SphericalGeometry::Mesh< Real > &sMesh , Parameters params )
		: _sMesh(sMesh) , _params(params)
	{
	}

	template< typename Real >
	void MoebiusCentering< Real >::advance( void )
	{
		Point3D< Real > c = _sMesh.center( );
		SquareMatrix< Real , 3 > D = _sMesh.dCenter( );
		if( _params.gaussNewton )
		{
			SquareMatrix< Real , 3 > D2 = ( D.transpose() * D ).inverse();
			c = - D2 * ( D * c );
		}
		else c = - D * c * 2;
		SphericalInversion< Real > inv;
		switch( _params.type )
		{
		case Parameters::PARAMETERS_CENTER_TO_INVERSION:   inv = typename SphericalGeometry::Mesh< Real >::CenterToInversion()                                 ( _sMesh , c ) ; break;
		case Parameters::PARAMETERS_GOLDEN_SECTION_SEARCH: inv = typename SphericalGeometry::Mesh< Real >::GSSCenterToInversion( _params.gssTolerance )        ( _sMesh , c ) ; break;
		case Parameters::PARAMETERS_POINCARE:              inv = typename SphericalGeometry::Mesh< Real >::PoincareCenterToInversion( _params.poincareMaxNorm )( _sMesh , c ) ; break;
		default: fprintf( stderr , "[ERROR] Unrecognized parameter type: %d\n" , _params.type ) , exit( 0 );
		}
#pragma omp parallel for
		for( int i=0 ; i<_sMesh.vertices.size() ; i++ ) _sMesh.vertices[i] = inv( _sMesh.vertices[i] );
	}

	template< typename Real >
	typename MoebiusCentering< Real >::Stats MoebiusCentering< Real >::getStats( void ) const
	{
		Stats stats;
		stats.center = _sMesh.center();
		return stats;
	}

}
