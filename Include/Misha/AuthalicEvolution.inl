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
	AuthalicEvolution< Real >::AuthalicEvolution( SphericalGeometry::Mesh< Real > &sMesh , const std::vector< Point3D< Real > > &meshVertices , Parameters params )
		: _sMesh(sMesh) , _meshVertices( meshVertices ) , _params(params)
	{
		auto Triangulate = []( const std::vector< Point3D< Real > > &vertices , const std::vector< std::vector< int > > &polygons )
		{
			std::vector< TriangleIndex > triangles;

			int count = 0;
			for( int i=0 ; i<polygons[i].size() ; i++ ) count += (int)polygons[i].size()-2;
			triangles.reserve( count );
#pragma omp parallel for
			for( int i=0 ; i<polygons.size() ; i++ )
			{
				MinimalAreaTriangulation< Real > MAT;
				const std::vector< int >& poly = polygons[i];
				std::vector< Point3D< Real > > _vertices( poly.size() );
				std::vector< TriangleIndex > _triangles;
				for( int j=0 ; j<poly.size() ; j++ ) _vertices[j] = vertices[ poly[j] ];
				MAT.GetTriangulation( _vertices , _triangles );
				for( int i=0 ; i<_triangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) _triangles[i][j] = poly[ _triangles[i][j] ];
#pragma omp critical
				triangles.insert( triangles.end() , _triangles.begin() , _triangles.end() );
			}
			return triangles;
		};

		// Triangulate the mesh
		_triangles = Triangulate( _meshVertices , _sMesh.polygons );

		// Compute the scale factor makeing the original mesh have area 4\pi
		{
			Real area = 0;
			for( int i=0 ; i<_triangles.size() ; i++ )
			{
				Point3D< Real > v[] = { _meshVertices[ _triangles[i][0] ] , _meshVertices[ _triangles[i][1] ] , _meshVertices[ _triangles[i][2] ] };
				area += Point3D< Real >::Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2;
			}
			_areaScale = (Real)sqrt( 4. * M_PI / area );
		}

		_sGrid.resize( _params.resolution );

		_flowField.setResolution( _sGrid.resolution() );
		_logScaleFactors.setResolution( _sGrid.resolution() );
		_hForm.resize( _sGrid.resolution() );
		_key.resize( _sGrid.resolution() );

		_setFlowField();
	}

	template< typename Real >
	void AuthalicEvolution< Real >::advance( void )
	{
		// Advect the vertex positions
		_advect();

		// Moebius center the mesh
		typename MoebiusCentering< Real >::Parameters parameters( _params.useGoldenSectionSearch ? MoebiusCentering< Real >::Parameters::PARAMETERS_GOLDEN_SECTION_SEARCH : MoebiusCentering< Real >::Parameters::PARAMETERS_POINCARE );
		MoebiusCentering< Real > mc( _sMesh , parameters );
		for( unsigned int i=0 ; i<_params.centeringIterations ; i++ ) mc.advance();

		// Compute the new scale factors and flow field
		_setFlowField();
	}

	template< typename Real >
	Real AuthalicEvolution< Real >::flippedTriangleFraction( void ) const
	{
		Real flippedArea = 0 , totalArea = 0;
#pragma omp parallel for reduction( + : flippedArea , totalArea )
		for( int t=0 ; t<_triangles.size() ; t++ )
		{
			Point3D< Real > v[] = { _sMesh.vertices[ _triangles[t][0] ] , _sMesh.vertices[ _triangles[t][1] ] , _sMesh.vertices[ _triangles[t][2] ] };
			Point3D< Real > c = ( v[0] + v[1] + v[2] ) / 3;
			Point3D< Real > n = Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
			Real a = (Real)sqrt( Point3D< Real >::SquareNorm( n ) );
			totalArea += a;
			if( Point3D< Real >::Dot( c , n )<0 ) flippedArea += a;
		}
		return flippedArea / totalArea;
	}

	template< typename Real >
	typename AuthalicEvolution< Real >::Stats AuthalicEvolution< Real >::getStats( void ) const
	{
		static std::vector< Real > logScaleFactors( _sMesh.vertices.size() );
		return getStats( logScaleFactors );
	}

	template< typename Real >
	typename AuthalicEvolution< Real >::Stats AuthalicEvolution< Real >::getStats( std::vector< Real > &logScaleFactors ) const
	{
		auto SetVertexWeights = []( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , std::vector< Real > &vertexWeights )
		{
#pragma omp parallel for
			for( int i=0 ; i<vertices.size() ; i++ ) vertexWeights[i] = 0;

			Real area = 0;
#pragma omp parallel for reduction( + : area )
			for( int i=0 ; i<triangles.size() ; i++ )
			{
				Point3D< Real > v[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
				Real a = Point3D< Real >::Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2;
				area += a;
				for( int j=0 ; j<3 ; j++ )
#pragma omp atomic
					vertexWeights[ triangles[i][j] ] += a/3;
			}
#pragma omp parallel for
			for( int i=0 ; i<vertices.size() ; i++ ) vertexWeights[i] /= area;
		};

		auto GetStats = []( const std::vector< Real > &signal , const std::vector< Real > &vWeights )
		{
			Stats stats;

			stats.min = stats.max = signal[0];
			stats.avg = stats.dev = 0;
			for( int i=0 ; i<signal.size() ; i++ )
			{
				stats.min = std::min< Real >( stats.min , signal[i] ) , stats.max = std::max< Real >( stats.max , signal[i] );
				stats.avg += signal[i] * vWeights[i];
				stats.dev += signal[i] * signal[i] * vWeights[i];
			}
			stats.dev = (Real)sqrt( stats.dev );
			return stats;
		};

		static std::vector< Real > vWeights( _sMesh.vertices.size() );

		_setLogScaleFactors( logScaleFactors );
		SetVertexWeights( _sMesh.vertices , _triangles , vWeights );
		return GetStats( logScaleFactors , vWeights );
	}

	template< typename Real >
	Point3D< Real > AuthalicEvolution< Real >::flowField( Point3D< Real > p ) const { return _flowField(p); }

	template< typename Real >
	void AuthalicEvolution< Real >::_setLogScaleFactors( std::vector< Real > &logScaleFactors ) const
	{
		auto SetScaleFactors = []( const std::vector< Point3D< Real > > &paramVertices , const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , std::vector< Real > &scaleFactors , Real areaScale )
		{
			scaleFactors.resize( vertices.size() );
			static std::vector< Real > paramAreas( vertices.size() );
#pragma omp parallel for
			for( int i=0 ; i<scaleFactors.size() ; i++ ) scaleFactors[i] = paramAreas[i] = 0;
			std::vector< Real > &areas = scaleFactors;
#pragma omp parallel for
			for( int i=0 ; i<triangles.size() ; i++ )
			{
				Point3D< Real > pv[] = { paramVertices[ triangles[i][0] ] , paramVertices[ triangles[i][1] ] , paramVertices[ triangles[i][2] ] };
				Point3D< Real >  v[] = {      vertices[ triangles[i][0] ] ,      vertices[ triangles[i][1] ] ,      vertices[ triangles[i][2] ] };
				Real pArea = Point3D< Real >::Length( Point3D< Real >::CrossProduct( pv[1]-pv[0] , pv[2]-pv[0] ) ) / 6;
				Real  area = Point3D< Real >::Length( Point3D< Real >::CrossProduct(  v[1]- v[0] ,  v[2]- v[0] ) ) / 6;
				for( int j=0 ; j<3 ; j++ )
				{
#pragma omp atomic
					paramAreas[ triangles[i][j] ] += pArea;
#pragma omp atomic
					areas[ triangles[i][j] ] += area;
				}
			}
#pragma omp parallel for
			for( int i=0 ; i<vertices.size() ; i++ ) scaleFactors[i] = areas[i] * areaScale / paramAreas[i];
		};

		// Compute the scale factors
		SetScaleFactors( _sMesh.vertices , _meshVertices , _triangles , logScaleFactors , _areaScale );

		// Clamp the scale factors
		if( _params.epsilon>0 )
#pragma omp parallel for
			for( int i=0 ; i<logScaleFactors.size() ; i++ )
			{
				if( logScaleFactors[i]!=logScaleFactors[i] ) fprintf( stderr , "[ERROR] Undefined scale factor (%d) %g\n" , i , logScaleFactors[i] ) , exit( 0 );
				else if( logScaleFactors[i]<   _params.epsilon ) logScaleFactors[i] = _params.epsilon;
				else if( logScaleFactors[i]>1./_params.epsilon ) logScaleFactors[i] = (Real)( 1./_params.epsilon );
			}

		// Compute the logarithm
#pragma omp parallel for
		for( int i=0 ; i<logScaleFactors.size() ; i++ ) logScaleFactors[i] = (Real)log( logScaleFactors[i] );
	}

	template< typename Real >
	void AuthalicEvolution< Real >::_advect( void )
	{
		auto Advect = []( Point3D< Real > p , Point3D< Real > v )
		{
			// Project into the tangent plane
			v -= p * Point3D< Real >::Dot( v , p );
			Real l = Point3D< Real >::Length( v );
			if( l ) p = p * cos( l ) + v / l * sin( l );
			return p;
		};
		struct Frame
		{
			Point3D< Real > p , t , n;
			Frame( void ){}
			Frame( Point3D< Real > point , Point3D< Real > tangent ) : p( point ) , t( tangent ) { n = Point3D< Real >::CrossProduct( p , t ); }
			static Point3D< Real > Transport( const Frame &start , const Frame &end , Point3D< Real > d )
			{
				return end.t * Point3D< Real >::Dot( start.t , d ) + end.n * Point3D< Real >::Dot( start.n , d );
			};
		};
		auto AdvectFrame = []( Point3D< Real > p , Point3D< Real > d , Frame &startFrame , Frame &endFrame )
		{
			d -= p * Point3D< Real >::Dot( d , p );
			Real l = Point3D< Real >::Length( d );
			if( l )
			{
				d /= l;
				Real c = (Real)cos(l) , s = (Real)sin(l);

				startFrame = Frame( p , d );
				endFrame = Frame( p * c + d * s , d * c - p * s );

				return true;
			}
			else return false;
		};
		auto SampleFlow = [&]( Point3D< Real > p )
		{
			Real x , y;
			SphericalGrid< Real >::SetCoordinates( &(p[0]) , x , y );
			x *= (Real)( _sGrid.resolution() / ( 2.*M_PI ) );
			y *= (Real)( _sGrid.resolution() / M_PI );
			int x0 = (int)floor(x) , y0 = (int)floor(y);
			int x1 = x0+1 , y1 = y0+1;
			y0 = std::max< int >( 0 , y0 ) , y1 = std::min< int >( y1 , _sGrid.resolution() );
			Real dx = x-x0 , dy = y-y0;
			return
				_flowField[ SphericalGeometry::Tessellation< Real >::GridVertexIndex( x0 , y0 , _sGrid.resolution() ) ] * (Real)( ( 1.-dx ) * ( 1.-dy ) ) +
				_flowField[ SphericalGeometry::Tessellation< Real >::GridVertexIndex( x0 , y1 , _sGrid.resolution() ) ] * (Real)( ( 1.-dx ) * (    dy ) ) +
				_flowField[ SphericalGeometry::Tessellation< Real >::GridVertexIndex( x1 , y0 , _sGrid.resolution() ) ] * (Real)( (    dx ) * ( 1.-dy ) ) +
				_flowField[ SphericalGeometry::Tessellation< Real >::GridVertexIndex( x1 , y1 , _sGrid.resolution() ) ] * (Real)( (    dx ) * (    dy ) ) ;
		};

#pragma omp parallel for
		for( int i=0 ; i<_sMesh.vertices.size() ; i++ )
		{
			Point3D< Real > p = _sMesh.vertices[i];
			if( _params.useSemiImplicit )
			{
				for( unsigned int s=0 ; s<_params.subSteps ; s++ )
				{
					Frame start , end;
					Point3D< Real > v1 , v2;
					// Get the vector at the start position
					v1 = SampleFlow( p ) * _params.stepSize / _params.subSteps;
					if( AdvectFrame( p , v1 , start , end ) )
					{
						// Get the vector at the end and transport it back
						v2 = SampleFlow( end.p ) * _params.stepSize / _params.subSteps;
						v2 = Frame::Transport( end , start , v2 );
						v1 = (v1 + v2) / 2;
					}
					p = Advect( p , v1 );
				}
			}
			else for( unsigned int i=0 ; i<_params.subSteps ; i++ ) p = Advect( p , SampleFlow( p ) * _params.stepSize / _params.subSteps );
			_sMesh.vertices[i] = p;
		}
	}

	template< typename Real >
	void AuthalicEvolution< Real >::_setFlowField( void )
	{
		// Compute the scale factors for every spherical cell
		SphericalGeometry::Tessellation< Real >::SampleFaceIntegrals( _sMesh , _sMesh.masses , _sGrid , _params.epsilon );
		// Take the logarithm
#pragma omp parallel for
		for( int i=0 ; i<_sGrid.resolution() ; i++ ) for( int j=0 ; j<_sGrid.resolution() ; j++ ) _sGrid(i,j) = (Real)log( _sGrid(i,j)*_areaScale );
		_logScaleFactors.set( _sGrid );

		// Smooth
		{
			_hForm.ForwardFourier( _sGrid , _key );
			for( int i=0 ; i<_key.bandWidth() ; i++ ) for( int j=0 ; j<=i ; j++ ) _key(i,j) *= (Real)exp( -_params.smoothness * i * ( (Real)(i+1) ) * 1. );
			_hForm.InverseFourier( _key , _sGrid );

		}
		// Sample the smoothed (log) scale factors to the sphere vertices
		_logScaleFactors.set( _sGrid );

		{
			// Start by computing the per-face flow field
			_logScaleFactors.setGradients( _gradientField );

			// Negate the flow field
#pragma omp parallel for
			for( int i=0 ; i<_gradientField.size() ; i++ ) _gradientField[i] = -_gradientField[i];

			// Sample to the vertices
			_flowField.set( _gradientField );
		}
	}
}
