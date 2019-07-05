/*
Copyright (c) 2018, Michael Kazhdan, Alex Baden, and Keenan Crane
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

#ifndef SIGNALS_INCLUDED
#define SIGNALS_INCLUDED


#include "Misha/SourceToTargetCorrespondence.h"
#include "Misha/FEM.h"
#include "Misha/Solver.h"

template< typename Real >
struct Signal
{
	struct Smoother
	{
		typedef EigenSolverCholeskyLDLt< Real > Solver;
		Smoother( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles );
		Smoother( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , Real sWeight );
		~Smoother( void );
		void setStiffnessWeight( Real sWeight );
		void smooth( ConstPointer( Real ) in , Pointer( Real ) out );
	protected:
		void _init( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles );
		SparseMatrix< Real , int  > _M , _S , _L;
		std::vector< Real > _b;
		Solver *_solver;
	};

	std::vector< std::vector< Real > > values;
	void read( const char *fileName );
	void write( const char *fileName ) const;

	static void SetHKS( const Spectrum< Real > &spectrum , Real time , std::vector< Real > &hks );
	static void Normalize( std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles );
	static void SetNormals( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , std::vector< Point3D< Real > > &normals );
	static void SmoothNormals( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , std::vector< Point3D< Real > > &normals , Real sWeight , int iters );
	static void SetGaussAndMeanCurvatures( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , const std::vector< Point3D< Real > > &normals , std::vector< Point2D< Real > > &curvatures );
protected:
	static const unsigned long long _MAGIC_NUMBER;
};

////////////
// Signal //
////////////
template< typename Real > const unsigned long long Signal< Real >::_MAGIC_NUMBER = 0x0516ull;

template< class Real >
void Signal< Real >::SetHKS( const Spectrum< Real > &spectrum , Real time , std::vector< Real > &hks )
{
	hks.resize( spectrum.eVector(0).size() , 0 );
#pragma omp parallel for
	for( int i=0 ; i<hks.size() ; i++ ) for( int j=0 ; j<spectrum.size() ; j++ ) hks[i] += (Real)exp( - spectrum.eValue(j) * time ) * spectrum.eVector(j)[i] * spectrum.eVector(j)[i];
}

template< class Real >
void Signal< Real >::Normalize( std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles )
{
	auto Area = [&]( int idx )
	{
		return (Real)Length( Point3D< Real >::CrossProduct( vertices[ triangles[idx][1] ] - vertices[ triangles[idx][0] ] , vertices[ triangles[idx][2] ] - vertices[ triangles[idx][0] ] ) ) / 2;
	};
	Real area = 0;
	Point3D< Real > center;
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Real a = Area( i );
		Point3D< Real > c = ( vertices[ triangles[i][0] ] + vertices[ triangles[i][1] ] + vertices[ triangles[i][2] ]) / 3;
		center += c * a;
		area += a;
	}
	center /= area;
	Real scl = (Real)sqrt( 1./area );
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = ( vertices[i] - center ) * scl;
}

template< class Real >
void Signal< Real >::SetNormals( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , std::vector< Point3D< Real > > &normals )
{
	normals.resize( vertices.size() );
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< Real > v[] = { Point3D< Real >( vertices[ triangles[i][0] ] ) , Point3D< Real >( vertices[ triangles[i][1] ] ) , Point3D< Real >( vertices[ triangles[i][2] ] ) };
		Point3D< Real > n = Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
		for( int j=0 ; j<3 ; j++ ) normals[ triangles[i][j] ] += n;
	}
	for( int i=0 ; i<normals.size() ; i++ ) normals[i] /= (Real)Length( normals[i] );
}

template< class Real >
void Signal< Real >::SmoothNormals( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , std::vector< Point3D< Real > > &normals , Real sWeight , int iters )
{
	typedef EigenSolverCholeskyLDLt< Real > Solver;

	std::vector< TriangleIndex > _triangles = triangles;
	FEM::RiemannianMesh< Real > mesh( GetPointer( _triangles ) , triangles.size() );
	mesh.setMetricFromEmbedding( GetPointer( vertices ) );
	mesh.makeUnitArea();
	SparseMatrix< Real , int > M , _M = mesh.template massMatrix< FEM::BASIS_0_WHITNEY >() , _S = mesh.template stiffnessMatrix< FEM::BASIS_0_WHITNEY >();
	M.resize( 2*vertices.size() );
#pragma omp parallel for
	for( int i=0 ; i<vertices.size() ; i++ ) for( int ii=0 ; ii<2 ; ii++ )
	{
		M.SetRowSize( 2*i+ii , 2*_M.rowSizes[i] );
		for( int j=0 ; j<_M.rowSizes[i] ; j++ ) for( int jj=0 ; jj<2 ; jj++ ) M[2*i+ii][2*j+jj].N = _M[i][j].N*2+jj;
	}
	std::vector< Point3D< Real > > tangents( vertices.size()*2 );
	std::vector< Real > b( vertices.size()*2 ) , o( vertices.size()*2 );

	Solver solver( M , true );

	Real mWeight = (Real)1.;
	if( sWeight<=0 )
	{
		fprintf( stderr , "[WARNING] Ascribing no weight to the mass term\n" );
		mWeight = 0 , sWeight = (Real)1.;
	}

	for( int iter=0 ; iter<iters ; iter++ )
	{
		Timer t;

		// Set the tangent directions
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ )
		{
			Point3D< Real > v( 1 , 0 , 0 );
			if( fabs( Point3D< Real >::Dot( v , normals[i] ) )>0.99 ) v = Point3D< Real >( 0 , 1 , 0 );
			tangents[2*i+0] = Point3D< Real >::CrossProduct( normals[i] , v               ) ; tangents[2*i+0] /= Length( tangents[2*i+0] );
			tangents[2*i+1] = Point3D< Real >::CrossProduct( normals[i] , tangents[2*i+0] ) ; tangents[2*i+1] /= Length( tangents[2*i+1] );
		}

		// Solve for the tangent offsets minimizing the dirichlet energy:
		// E( o1 , o2 ) = || \sum o[i] * T[i] ||^2 + e * || \nabla( \sum n[i] + o[i] * T[i] ) ||^2
		//              = o^t * T^t * M * T * o + e * [ o^t * T^t * S * T * o + 2 * o^t * T^t * S * n + n^t * S * n ]
		// \nabla E = 0:
		// 0 = T^t * ( M + e * S ) * T * o + e * T^t * S * n
		{
			Timer t;
#pragma omp parallel for 
			for( int i=0 ; i<vertices.size() ; i++ ) for( int ii=0 ; ii<2 ; ii++ ) 
			{
				b[2*i+ii] = 0;
				for( int j=0 ; j<_M.rowSizes[i] ; j++ )
				{
					for( int jj=0 ; jj<2 ; jj++ ) M[2*i+ii][2*j+jj].Value = ( _M[i][j].Value*mWeight + _S[i][j].Value * sWeight ) * Point3D< Real >::Dot( tangents[2*i+ii] , tangents[ 2*_M[i][j].N+jj ] );
					b[2*i+ii] -= _S[i][j].Value * Point3D< Real >::Dot( normals[ _S[i][j].N ] , tangents[2*i+ii] ) * sWeight;
				}
			}
		}
		{
			Timer t;
			solver.update( M );
			solver.solve( GetPointer( b ) , GetPointer( o ) );
#pragma omp parallel for
			for( int i=0 ; i<vertices.size() ; i++ ) normals[i] += tangents[2*i+0] * o[2*i+0] + tangents[2*i+1] * o[2*i+1] , normals[i] /= Length( normals[i] );
		}
	}
}

template< typename Real >
void Signal< Real >::SetGaussAndMeanCurvatures( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , const std::vector< Point3D< Real > > &normals , std::vector< Point2D< Real > > &curvatures )
{
	auto SecondFundamentalForm = []( const Point3D< Real > v[3] , const Point3D< Real > n[3] )
	{
		SquareMatrix< Real , 2 > II;
		Point3D< Real > dv[] = { v[1]-v[0] , v[2]-v[0] } , dn[] = { n[1]-n[0] , n[2]-n[0] };
		for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) II( i , j ) = Point3D< Real >::Dot( dv[i] , dn[j] );
		II(1,0) = II (0,1) = ( II(0,1) + II(1,0) ) / (Real)2.;
		return II;
	};

	std::vector< TriangleIndex > _triangles = triangles;
	FEM::RiemannianMesh< Real > mesh( GetPointer( _triangles ) , triangles.size() );
	mesh.setMetricFromEmbedding( GetPointer( vertices ) );
	Real area = mesh.area();
	mesh.makeUnitArea();

	curvatures.resize( vertices.size() );
	std::vector< Point2D< Real > > _curvatures( triangles.size() );
	{
		Real s = (Real)( 1. / sqrt(area) );
#pragma omp parallel for
		for( int i=0 ; i<triangles.size() ; i++ )
		{
			Point3D< Real > v[] = { vertices[ triangles[i][0] ] * s , vertices[ triangles[i][1] ] * s , vertices[ triangles[i][2] ] * s };
			Point3D< Real > n[] = { normals [ triangles[i][0] ]     , normals [ triangles[i][1] ]     , normals [ triangles[i][2] ]     };
			SquareMatrix< Real , 2 > II = SecondFundamentalForm( v , n );

			// The columns of A_inverse given an orthonormal frame with respect to g.
			SquareMatrix< Real , 2 > A = FEM::TensorRoot( mesh.g(i) ) , A_inverse = A.inverse() , D , R;
			// The matrix A_inverse^t * II * A_inverse gives the second fundamental form with respect to a basis that is orthonormal w.r.t. g
			SquareMatrix< Real , 2 > _II = A_inverse.transpose() * II * A_inverse;
			_curvatures[i] = Point2D< Real >( _II.determinant() , _II.trace() );
		}

		{
			std::vector< Real > vAreas( vertices.size() , 0 );
			for( int i=0 ; i<triangles.size() ; i++ )
			{
				Real area = mesh.area(i);
				for( int j=0 ; j<3 ; j++ ) vAreas[ triangles[i][j] ] += area , curvatures[ triangles[i][j] ] += _curvatures[i] * area;
			}
			for( int i=0 ; i<vertices.size() ; i++ ) curvatures[i] /= vAreas[i];
		}
	}
}

template< typename Real >
void Signal< Real >::read( const char *fileName )
{
	FILE* fp = fopen( fileName , "rb" );
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open file for reading: %s\n" , fileName ) , exit( 0 );

	/// Read simple header
	// Read Magic Number
	unsigned long long magicNum;
	if( fread( &magicNum , sizeof(unsigned long long) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] Failed to read magic number\n" ) , exit( 0 );
	if( magicNum!=_MAGIC_NUMBER ) fprintf( stderr , "[ERROR] Bad magin number: %s\n" , fileName ) , exit( 0 );

	// type size
	unsigned int typeSize;
	if( fread( &typeSize , sizeof(unsigned int) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] Failed to read type size\n" ) , exit( 0 );
	bool needCasting = typeSize!=sizeof(Real);
	if( needCasting ) fprintf( stderr , "[WARNING] File type size does not match: %d != %d.  Will cast.\n" , (int)typeSize , (int)sizeof(Real) );

	auto Read = [&]( void )
	{
		if( needCasting )
		{
			switch( typeSize )
			{
			case 4:
			{
				float _value;
				if( fread( &_value , sizeof(float) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] Failed to read value\n" );
				return (Real)_value;
			}
			case 8:
			{
				double _value;
				if( fread( &_value , sizeof(double) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] Failed to read value\n" );
				return (Real)_value;
			}
			case 10:
			{
				long double _value;
				if( fread( &_value , sizeof(long double) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] Failed to read value\n" );
				return (Real)_value;
			}
			default:
				fprintf( stderr , "[ERROR] Unrecognized type\n" ) , exit( 0 );
				return (Real)-1;
			}
		}
		else
		{
			Real value;
			if( fread( &value , sizeof(Real) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] Failed to read value\n" );
			return value;
		}
	};

	// Num of Channels
	int numChannels , dimension;
	if( fread( &numChannels , sizeof(int) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] Failed to read number of channels\n" );
	if( fread( &dimension , sizeof(int) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] Failed to read dimension\n" );

	values.resize( numChannels );
	for( int c=0 ; c<numChannels ; c++ ) values[c].resize( dimension );

	/// Content
	if ( needCasting )
	{
		// Signals
		unsigned char* tempMemory = new unsigned char[ typeSize * dimension ];
		for( int c=0 ; c<numChannels ; c++ )
		{
			if( fread( tempMemory , typeSize , dimension , fp )!=dimension ) fprintf( stderr , "[ERROR] Failed to read temporary memory\n" ) , exit( 0 );
#pragma omp parallel for
			for( int j=0 ; j<dimension ; j++ )
			{
				switch( typeSize )
				{
				case  4: values[c][j] = (Real)reinterpret_cast< float *       >( tempMemory )[j] ; break;
				case  8: values[c][j] = (Real)reinterpret_cast< double *      >( tempMemory )[j] ; break;
				case 10: values[c][j] = (Real)reinterpret_cast< long double * >( tempMemory )[j] ; break;
				}
			}
		}
		delete[] tempMemory;
	}
	else for( int c=0 ; c<numChannels ; c++ ) if( fread( values[c].data() , sizeof(Real) , dimension , fp )!=dimension ) fprintf( stderr , "[ERROR] Failed to read values\n" ) , exit( 0 );

	fclose( fp );
}
template< class Real >
void Signal< Real >::write( const char *fileName ) const
{
	FILE* fp = fopen( fileName , "wb" );
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open file for writing: %s\n" , fileName ) , exit( 0 );

	// Read simple header
	// Read Magic Number
	fwrite( &_MAGIC_NUMBER , sizeof(unsigned long long) , 1 , fp );

	// type size
	unsigned int typeSize = sizeof(Real);
	fwrite( &typeSize , sizeof(unsigned int) , 1 , fp );

	// Num of Channels
	int numChannels = (int)values.size() , dimension = (int)values[0].size();
	fwrite( &numChannels , sizeof(int) , 1 , fp );
	fwrite( &dimension , sizeof(int) , 1 , fp );
	for( int c=0 ; c<numChannels ; c++ ) fwrite( values[c].data() , sizeof(Real) , dimension , fp );

	fclose( fp );
}

//////////////////////
// Signal::Smoother //
//////////////////////
template< typename Real >
Signal< Real >::Smoother::~Smoother( void ){ delete _solver; }
template< typename Real >
Signal< Real >::Smoother::Smoother( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles )
{
	_init( vertices , triangles );
}
template< typename Real >
Signal< Real >::Smoother::Smoother( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , Real sWeight )
{
	_init( vertices , triangles );
	setStiffnessWeight( sWeight );
}
template< typename Real >
void Signal< Real >::Smoother::_init( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles )
{
	std::vector< TriangleIndex > _triangles = triangles;
	FEM::RiemannianMesh< Real > mesh( GetPointer( _triangles ) , triangles.size() );
	mesh.setMetricFromEmbedding( GetPointer( vertices ) );
	mesh.makeUnitArea();
	_M = mesh.template massMatrix< FEM::BASIS_0_WHITNEY >();
	_S = mesh.template stiffnessMatrix< FEM::BASIS_0_WHITNEY >();
	_L.resize( _M.rows );
#pragma omp parallel for
	for( int i=0 ; i<_S.rows ; i++ )
	{
		_L.SetRowSize( i , _M.rowSizes[i] );
		for( int j=0 ; j<_M.rowSizes[i] ; j++ ) _L[i][j].N = _M[i][j].N;
	}
	_b.resize( _M.rows );
	_solver = new Solver( _L , true );
}
template< typename Real >
void Signal< Real >::Smoother::setStiffnessWeight( Real sWeight )
{
#pragma omp parallel for
	for( int i=0 ; i<_S.rows ; i++ ) for( int j=0 ; j<_S.rowSizes[i] ; j++ ) _L[i][j].Value = _S[i][j].Value * sWeight + _M[i][j].Value;
	_solver->update( _L );
}
template< typename Real >
void Signal< Real >::Smoother::smooth( ConstPointer( Real ) in , Pointer( Real ) out )
{
	_M.Multiply( in , GetPointer(_b) );
	_solver->solve( GetPointer(_b) , out );
}

#endif // SIGNALS_INCLUDED