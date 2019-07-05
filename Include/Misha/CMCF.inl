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
	template< typename Real > const Real CMCF< Real >::_QuasiConformalCutOff = (Real)1e-15;

	template< typename Real >
	CMCF< Real >::CMCF( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , SphericalGeometry::Mesh< Real > &sMesh , Parameters params )
		: _sMesh(sMesh) , _params(params)
	{
		_sMesh.vertices = vertices;
		_sMesh.polygons.resize( triangles.size() );
		for( int i=0 ; i<triangles.size() ; i++ )
		{
			_sMesh.polygons[i].resize(3);
			for( int j=0 ; j<3 ; j++ ) _sMesh.polygons[i][j] = triangles[i][j];
		}

		_normalize();

		_A = _S = _stiffnessMatrix();
		_solver = new Solver( _A , true );

		_b.resize( _sMesh.vertices.size() );
		_x.resize( _sMesh.vertices.size() );
		__b.resize( _sMesh.vertices.size() );
		__x.resize( _sMesh.vertices.size() );

		_massInvs.resize( _sMesh.polygons.size() );
		_areas.resize( _sMesh.polygons.size() , 0 );

		Real area = 0;
#pragma omp parallel for reduction( + : area )
		for( int i=0 ; i<_sMesh.polygons.size() ; i++ )
		{
			Point3D< Real > v[] = { _sMesh.vertices[ _sMesh.polygons[i][0] ] , _sMesh.vertices[ _sMesh.polygons[i][1] ] , _sMesh.vertices[ _sMesh.polygons[i][2] ] };
			_areas[i] = (Real)Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2;
			area += _areas[i];
			if( _areas[i]<=_QuasiConformalCutOff ) continue;
			_massInvs[i] = __TriangleMassMatrix( v ).inverse();
		}
#pragma omp parallel for
		for( int i=0 ; i<_areas.size() ; i++ ) _areas[i] /= area;
	}

	template< typename Real >
	CMCF< Real >::~CMCF( void ){ delete _solver; }

	template< typename Real >
	void CMCF< Real >::advance( void )
	{
		_M = _massMatrix();
		_M.Multiply( &_sMesh.vertices[0] , &_b[0] );
#pragma omp parallel for
		for( int j=0 ; j<_A.rows ; j++ ) for( int k=0 ; k<_A.rowSizes[j] ; k++ ) _A[j][k].Value = _S[j][k].Value * _params.stepSize + _M[j][k].Value;
		_solver->update( _A );

		for( int d=0 ; d<3 ; d++ )
		{
#pragma omp parallel for
			for( int j=0 ; j<_b.size() ; j++ ) __b[j] = _b[j][d] , __x[j] = _x[j][d] = _sMesh.vertices[j][d];
			_solver->solve( &__b[0] , &__x[0] );
#pragma omp parallel for
			for( int j=0 ; j<_sMesh.vertices.size() ; j++ ) _sMesh.vertices[j][d] = __x[j];
		}
		_normalize();
	}

	template< typename Real >
	typename CMCF< Real >::Stats CMCF< Real >::getStats( void ) const
	{
		Stats stats;
		stats.deformationScale = _deformationScale();
		stats.quasiConformalRatio = _quasiConformalRatio();
		stats.radialDeviation = _radialDeviation();
		return stats;
	}

	template< typename Real >
	Real CMCF< Real >::_deformationScale( void )  const
	{
		Real differenceNorm=0 , oldNorm=0;
#pragma omp parallel for reduction( + : differenceNorm , oldNorm )
		for( int i=0 ; i<_M.rows ; i++ ) for( int j=0 ; j<_M.rowSizes[i] ; j++ ) 
		{
			differenceNorm += Point3D< Real >::Dot( _x[i]-_sMesh.vertices[_M[i][j].N] , _x[i]-_sMesh.vertices[_M[i][j].N] ) * _M[i][j].Value;
			oldNorm += Point3D< Real >::Dot( _x[i] , _x[_M[i][j].N] ) * _M[i][j].Value;
		}
		return (Real)sqrt( differenceNorm / oldNorm );
	}
	template< typename Real >
	Real CMCF< Real >::_radialDeviation( void ) const
	{
		Real area = 0;
		std::vector< Real > vertexAreas( _sMesh.vertices.size() , 0 );
#pragma omp parallel for
		for( int i=0 ; i<_sMesh.polygons.size() ; i++ )
		{
			Point3D< Real > v[] = { _sMesh.vertices[ _sMesh.polygons[i][0] ] , _sMesh.vertices[ _sMesh.polygons[i][1] ] , _sMesh.vertices[ _sMesh.polygons[i][2] ] };
			Real a = (Real)Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2;
			area += a;
			a /= 3;
			for( int j=0 ; j<3 ; j++ )
#pragma omp atomic
				vertexAreas[ _sMesh.polygons[i][j] ] += a;
		}
		Point3D< Real > center;
		for( int i=0 ; i<_sMesh.vertices.size() ; i++ ) center += _sMesh.vertices[i] * vertexAreas[i];
		center /= area;

		// set v = \sum_i[ w_i * v_i ]
		// var = \sum_i[ w_i * ( v_i -v )^2 ]
		//     = \sum_i[ w_i * ( v_i^2 + v^2 - 2 * v_i * v ) ]
		//     = \sum_i[ w_i * v_i^2 ] + \sum_i[ w_i * v^2 ] - \sum_i[ 2 * w_i * v_i * v ) ]
		//     = \sum_i[ w_i * v_i^2 ] + v^2 - 2 * v^2
		//     = \sum_i[ w_i * v_i^2 ] - v^2
		Real average = 0 , var = 0;
#pragma omp parallel for reduction( + : average , var )
		for( int i=0 ; i<_sMesh.vertices.size() ; i++ )
		{
			Real r = (Real)sqrt( Point3D< Real >::SquareNorm( Point3D< Real >( _sMesh.vertices[i]-center ) ) );
			average += r * vertexAreas[i];
			var += r * r * vertexAreas[i];
		}
		average /= area , var /= area;
		var -= average * average;
		var = (Real)sqrt( var );
		return var / average;
	}

	template< typename Real >
	void CMCF< Real >::_normalize( void )
	{
		Real area = 0 , centerX = 0 , centerY = 0 , centerZ = 0;

#pragma omp parallel for reduction( + : area , centerX , centerY , centerZ )
		for( int i=0 ; i<_sMesh.polygons.size() ; i++ )
		{
			Point3D< Real > v[] = { _sMesh.vertices[ _sMesh.polygons[i][0] ] , _sMesh.vertices[ _sMesh.polygons[i][1] ] , _sMesh.vertices[ _sMesh.polygons[i][2] ] };
			Point3D< Real > c = ( v[0] + v[1] + v[2] ) / 3;
			Real a = (Real)Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) );
			centerX += c[0] * a;
			centerY += c[1] * a;
			centerZ += c[2] * a;
			area += a;
		}
		Point3D< Real > translate( -centerX/area , -centerY/area , -centerZ/area );
		Real scale = (Real)1./sqrt(area);
#pragma omp parallel for
		for( int i=0 ; i<_sMesh.vertices.size() ; i++ ) _sMesh.vertices[i] = ( _sMesh.vertices[i] + translate ) * scale;
	}

	template< typename Real >
	template< typename TriangleMatrixFunctor >
	SparseMatrix< Real , int > CMCF< Real >::_systemMatrix( TriangleMatrixFunctor F ) const
	{
		struct Entry
		{
			Entry( void ) : row(-1) , col(-1) , value(0){}
			Entry( int r , int c , Real v ) : row(r) , col(c) , value(v){}
			int row , col;
			Real value;
		};
		SparseMatrix< Real , int > M;
		M.resize( (int)_sMesh.vertices.size() );
		std::vector< std::vector< Entry > > entries( omp_get_max_threads() );
#pragma omp parallel for
		for( int t=0 ; t<_sMesh.polygons.size() ; t++ )
		{
			int thread = omp_get_thread_num();
			Point3D< Real > v[] = { _sMesh.vertices[ _sMesh.polygons[t][0] ] , _sMesh.vertices[ _sMesh.polygons[t][1] ] , _sMesh.vertices[ _sMesh.polygons[t][2] ] };
			SquareMatrix< Real , 3 > m = F( v );
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) entries[thread].push_back( Entry( _sMesh.polygons[t][i] , _sMesh.polygons[t][j] , m(i,j) ) );
		}
		for( int i=0 ; i<entries.size() ; i++ ) for( int j=0 ; j<entries[i].size() ; j++ )	M.rowSizes[ entries[i][j].row ]++;
#pragma omp parallel for
		for( int i=0 ; i<M.rows ; i++ )
		{
			int rowSize = (int)M.rowSizes[i];
			M.rowSizes[i] = 0;
			M.SetRowSize( i , rowSize );
			M.rowSizes[i] = 0;
		}
		for( int i=0 ; i<entries.size() ; i++ ) for( int j=0 ; j<entries[i].size() ; j++ ) M[ entries[i][j].row ][ M.rowSizes[entries[i][j].row]++ ] = MatrixEntry< Real , int >( entries[i][j].col , entries[i][j].value );
#pragma omp parallel for
		for( int i=0 ; i<M.rows ; i++ )
		{
			std::unordered_map< int , Real > row;
			for( int j=0 ; j<M.rowSizes[i] ; j++ ) row[ M[i][j].N ] += M[i][j].Value;
			M.SetRowSize( i , (int)row.size() );
			int j=0;
			for( typename std::unordered_map< int , Real >::const_iterator iter=row.begin() ; iter!=row.end() ; iter++ ) M[i][j++] = MatrixEntry< Real , int >( iter->first , iter->second );
		}
		return M;
	}

	template< typename Real >
	SquareMatrix< Real , 2 > CMCF< Real >::__TriangleMassMatrix( const Point3D< Real > v[] )
	{
		Point3D< Real > t[] = { v[1]-v[0] , v[2]-v[0] };
		SquareMatrix< Real , 2 > mass;
		for( int  j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ ) mass(j,k) = Point3D< Real >::Dot( t[j] , t[k] );
		return mass;
	}

	template< typename Real >
	SquareMatrix< Real , 3 > CMCF< Real >::_TriangleMassMatrix( const Point3D< Real > vertices[] )
	{
		SquareMatrix< Real , 3 > mass;
		Real area = (Real)Length( Point3D< Real >::CrossProduct( vertices[1]-vertices[0] , vertices[2]-vertices[0] ) ) / 2;
		mass(1,0) = mass(0,1) = mass(2,0) = mass(0,1) = mass(0,2) = mass(1,2) = mass(2,1) = area / 12;
		mass(0,0) = mass(1,1) = mass(2,2) = area / 6;
		return mass;
	}

	template< typename Real >
	SquareMatrix< Real , 3 > CMCF< Real >::_TriangleStiffnessMatrix( const Point3D< Real > vertices[] )
	{
		SquareMatrix< Real , 3 > stiffness;
		Real area = (Real)Length( Point3D< Real >::CrossProduct( vertices[1]-vertices[0] , vertices[2]-vertices[0] ) ) / 2;
		for( int i=0 ; i<3 ; i++ )
		{
			int i0 = (i+1)%3 , i1 = (i+2)%3;
			Real dot = Point3D< Real >::Dot( vertices[i0]-vertices[i] , vertices[i1]-vertices[i] );
			stiffness(i0,i1) = stiffness(i1,i0) = - dot / area;
		}
		stiffness(0,0) = - ( stiffness(0,1) + stiffness(0,2) );
		stiffness(1,1) = - ( stiffness(1,2) + stiffness(1,0) );
		stiffness(2,2) = - ( stiffness(2,0) + stiffness(2,1) );
		return stiffness;
	}

	template< typename Real >
	Real CMCF< Real >::_quasiConformalRatio( void ) const
	{
		Real error = 0 , area = 0;
#pragma omp parallel for reduction( + : error , area )
		for( int i=0 ; i<_sMesh.polygons.size() ; i++ )
		{
			Point3D< Real > v[] = { _sMesh.vertices[ _sMesh.polygons[i][0] ] , _sMesh.vertices[ _sMesh.polygons[i][1] ] , _sMesh.vertices[ _sMesh.polygons[i][2] ] };
			SquareMatrix< Real , 2 > m = _massInvs[i] * __TriangleMassMatrix( v );
			Real det = m.determinant() , trace = m.trace();
			trace /= 2;
			if( det/4<=_QuasiConformalCutOff*_QuasiConformalCutOff ) continue;
			// The eigenvectors of this matrix are the roots of
			// P(x) = [x-d(0,0)]*[x-d(1,1)] - d(1,0)*d(0,1)
			//      = x^2 - x * Tr(d) + Det(d)
			// x = Tr(d)/2 +/- sqrt( Tr^2(d)/4 - Det(d) )
			Real disc  = trace*trace-det;
			if( disc<=Real(0) ) disc = 0;
			else                disc = sqrt( disc );
			Real x1 = trace - disc;
			Real x2 = trace + disc;
			if( x1<=0 ) continue;
			Real tError = sqrt( x2/x1 ) * _areas[i];
			error += tError;
			area += _areas[i];
		}
		return error / area;
	}

}
