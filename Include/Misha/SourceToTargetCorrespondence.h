/*
Copyright (c) 2019
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

#ifndef SOURCE_TO_TARGET_CORRESPONDENCE_INCLUDED
#define SOURCE_TO_TARGET_CORRESPONDENCE_INCLUDED

#include "Misha/Ply.h"
#include "Misha/SphericalGeometry.h"

#include <Spectra/SymGEigsSolver.h>
#include "Misha/FEM.h"
#include "Misha/Solver.h"


template< typename Real >
struct SourceTargetCorrespondences
{
	struct BarycentricCoordinates{ size_t tIndex ; Point3D< Real > weights; };

	static std::vector< BarycentricCoordinates > ReadCorrespondences( const char *fileName );
	static void WriteCorrespondences( const char *fileName , const std::vector< BarycentricCoordinates > &correspondences );

	std::vector< Point3D< Real > > vertices[2];
	std::vector< TriangleIndex > triangles[2];
	std::vector< BarycentricCoordinates > correspondences[2];
	Point3D< Real > sourceToTarget( size_t idx ) const
	{
		Point3D< Real > weights = correspondences[0][idx].weights;
		size_t v[] = { triangles[1][ correspondences[0][idx].tIndex ][0] , triangles[1][ correspondences[0][idx].tIndex ][1] , triangles[1][ correspondences[0][idx].tIndex ][2] };
		return vertices[1][ v[0] ] * weights[0] + vertices[1][ v[1] ] * weights[1] + vertices[1][ v[2] ] * weights[2];
	}
	Point3D< Real > targetToSource( size_t idx ) const
	{
		Point3D< Real > weights = correspondences[1][idx].weights;
		size_t v[] = { triangles[0][ correspondences[1][idx].tIndex ][0] , triangles[0][ correspondences[1][idx].tIndex ][1] , triangles[0][ correspondences[1][idx].tIndex ][2] };
		return vertices[0][ v[0] ] * weights[0] + vertices[0][ v[1] ] * weights[1] + vertices[0][ v[2] ] * weights[2];
	}

	SourceTargetCorrespondences( const SphericalGeometry::Mesh< Real > &source , const std::vector< Point3D< Real > > &sourceMeshVertices , const SphericalGeometry::Mesh< Real > &target , const std::vector< Point3D< Real > > &targetMeshVertices , unsigned int res );
	SourceTargetCorrespondences( const char *sourceFile , const char *targetFile , unsigned int res );
	SourceTargetCorrespondences( const char *sourceFile , const char *targetFile , const char *sourceToTargetCorrespondenceFile , const char *targetToSourceCorrespondenceFile );
	SourceTargetCorrespondences( const SourceTargetCorrespondences &correspondences )
	{
		for( int m=0 ; m<2 ; m++ )
		{
			vertices[m] = correspondences.vertices[m];
			triangles[m] = correspondences.triangles[m];
			this->correspondences[m] = correspondences.correspondences[m];
		}
	}
protected:
	unsigned int _res;
	void _init( const SphericalGeometry::Mesh< Real > &source , const std::vector< Point3D< Real > > &sourceMeshVertices , const SphericalGeometry::Mesh< Real > &target , const std::vector< Point3D< Real > > &targetMeshVertices , unsigned int res );
	static Point3D< Real > _NearestPointOnTriangle( Point3D< Real > point , const Point3D< Real > triangle[3] , Real* b );
	static BarycentricCoordinates _NearestPoint( Point3D< Real > p , const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , const std::vector< std::vector< SphericalGeometry::Polygon< Real > > > &polygonCells , unsigned int res );
	static Point3D< Real > _NearestPointOnEdge( Point3D< Real > point , const Point3D< Real > edge[2] , Real& b0 , Real& b1 );
	static std::vector< TriangleIndex > _Triangulate( const std::vector< Point3D< Real > > &vertices , const std::vector< std::vector< int > > &polygons );
};

template< typename Real >
struct Spectrum
{
	void set( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , unsigned int dimension , Real offset , bool lump );
	void read( const std::string &fileName );
	void write( const std::string &fileName ) const;
	Real biharmonicDistance( unsigned int i , unsigned int j ) const;
	Real biharmonicDistance( unsigned int i , TriangleIndex tri , Point3D< Real> weights ) const;
	size_t size( void ) const { return _eigenvalues.size(); }
	const Real &eValue( unsigned int idx ) const { return _eigenvalues[idx]; }
	const std::vector< Real > &eVector( unsigned int idx ) const { return _eigenvectors[idx]; }
protected:
	static const unsigned long long _MAGIC_NUMBER;
	std::vector< Real > _eigenvalues;
	std::vector< std::vector< Real > > _eigenvectors;
};

////////////////////////////////
// SourceTargetCorrespondence //
////////////////////////////////
template< typename Real >
std::vector< typename SourceTargetCorrespondences< Real >::BarycentricCoordinates > SourceTargetCorrespondences< Real >::ReadCorrespondences( const char *fileName )
{
	std::vector< BarycentricCoordinates > correspondences;
	FILE *fp = fopen( fileName , "r" );
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open correspondence file for reading: %s\n" , fileName ) , exit( 0 );
	int idx;
	float w[3];
	BarycentricCoordinates bc;
	while( fscanf( fp , " %d %f %f %f" , &idx , w+0 , w+1 , w+2 )!=4 )
	{
		bc.tIndex = (size_t)idx;
		bc.weights[0] = (Real)w[0];
		bc.weights[1] = (Real)w[1];
		bc.weights[2] = (Real)w[2];
		correspondences.push_back( bc );
	}
	fclose( fp );
	return bc;
}

template< typename Real >
void SourceTargetCorrespondences< Real >::WriteCorrespondences( const char *fileName , const std::vector< BarycentricCoordinates > &correspondences )
{
	FILE *fp = fopen( fileName , "w" );
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open correspondence file for writing: %s\n" , fileName ) , exit( 0 );
	int idx;
	float w[3];
	for( int i=0 ; i<correspondences.size() ; i++ ) fprintf( fp , " %d %f %f %f\n" , (int)correspondences[i].tIndex , (float)correspondences[i].weights[0] , (float)correspondences[i].weights[1] , (float)correspondences[i].weights[2] );
	fclose( fp );
}

template< typename Real >
SourceTargetCorrespondences< Real >::SourceTargetCorrespondences( const SphericalGeometry::Mesh< Real > &source , const std::vector< Point3D< Real > > &sourceMeshVertices , const SphericalGeometry::Mesh< Real > &target , const std::vector< Point3D< Real > > &targetMeshVertices , unsigned int res )
{
	_init( source , sourceMeshVertices , target , targetMeshVertices , res );
}
template< typename Real >
void SourceTargetCorrespondences< Real >::_init( const SphericalGeometry::Mesh< Real > &source , const std::vector< Point3D< Real > > &sourceMeshVertices , const SphericalGeometry::Mesh< Real > &target , const std::vector< Point3D< Real > > &targetMeshVertices , unsigned int res )
{
	_res = res;
	vertices[0] = sourceMeshVertices , vertices[1] = targetMeshVertices;
	SphericalGeometry::Mesh< Real > sMeshes[] = { source , target };
	for( int m=0 ; m<2 ; m++ )
	{
		triangles[m] = _Triangulate( vertices[m] , sMeshes[m].polygons );

		sMeshes[m].polygons.resize( triangles[m].size() );
		for( int i=0 ; i<triangles[m].size() ; i++ )
		{
			sMeshes[m].polygons[i].resize(3);
			for( int j=0 ; j<3 ; j++ ) sMeshes[m].polygons[i][j] = triangles[m][i][j];
		}
	}
	for( int m=0 ; m<2 ; m++ )
	{
		std::vector< std::vector< SphericalGeometry::Polygon< Real > > > polygonCells( _res * _res );

		SphericalGeometry::Tessellation< Real > tessellator( sMeshes[(m+1)%2] , _res );
		const std::vector< SphericalGeometry::Polygon< Real > > p = tessellator.polygons();
		for( int j=0 ; j<tessellator.polygons().size() ; j++ ) polygonCells[ p[j].theta*_res + p[j].phi ].push_back( p[j] );

		correspondences[m].resize( vertices[m].size() );
#pragma omp parallel for
		for( int i=0 ; i<vertices[m].size() ; i++ ) correspondences[m][i] = _NearestPoint( sMeshes[m].vertices[i] , sMeshes[(m+1)%2].vertices , triangles[(m+1)%2] , polygonCells , res );
	}
}
template< typename Real >
SourceTargetCorrespondences< Real >::SourceTargetCorrespondences( const char *sourceFile , const char *targetFile , unsigned int res )
{
	SphericalGeometry::Mesh< Real > sMeshes[2];
	std::vector< Point3D< Real > > meshVertices[2];
	sMeshes[0].read( sourceFile , meshVertices[0] );
	sMeshes[1].read( targetFile , meshVertices[1] );

	_init( sMeshes[0] , meshVertices[0] , sMeshes[1] , meshVertices[1] , res );
}

template< typename Real >
SourceTargetCorrespondences< Real >::SourceTargetCorrespondences( const char *sourceFile , const char *targetFile , const char *sourceToTargetCorrespondenceFile , const char *targetToSourceCorrespondenceFile )
{
	auto ReadMesh = []( const char *fileName , std::vector< Point3D< Real > > &vertices , std::vector< TriangleIndex > &triangles )
	{
		std::vector< PlyVertex< float > > _vertices;
		int fileType;
		PlyReadTriangles( fileName , _vertices , triangles , PlyVertex< float >::ReadProperties , NULL , PlyVertex< float >::ReadComponents , fileType );
		vertices.resize( _vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = Point3D< Real >( _vertices[i] );
	};
	ReadMesh( sourceFile , vertices[0] , triangles[0] );
	ReadMesh( targetFile , vertices[1] , triangles[1] );

	for( int m=0 ; m<2 ; m++ )
	{
		const char *correspondenceFile = m==0 ? sourceToTargetCorrespondenceFile : targetToSourceCorrespondenceFile;
		FILE *fp = fopen( correspondenceFile , "r" );
		if( !fp ) fprintf( stderr , "[ERROR] Failed to open correspondence file for reading: %s\n" , correspondenceFile ) , exit( 0 );
		correspondences[m].resize( vertices[m].size() );
		for( int i=0 ; i<vertices[m].size() ; i++ )
		{
			int idx;
			float w[3];
			if( fscanf( fp , " %d %f %f %f" , &idx , w+0 , w+1 , w+2 )!=4 ) fprintf( stderr , "[ERROR] Failed to read correspondence from file: %d\n" , i ) , exit( 0 );
			else
			{
				correspondences[m][i].tIndex = idx;
				correspondences[m][i].weights = Point3D< Real >( w[0] , w[1] , w[2] );
			}
		}
		fclose( fp );
	}
}

template< typename Real >
typename SourceTargetCorrespondences< Real >::BarycentricCoordinates SourceTargetCorrespondences< Real >::_NearestPoint( Point3D< Real > p , const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , const std::vector< std::vector< SphericalGeometry::Polygon< Real > > > &polygonCells , unsigned int res )
{
	BarycentricCoordinates bc;

	Real theta , phi; 

	SphericalGrid< Real >::SetCoordinates( p.coords , theta , phi );
	if( theta<0 ) theta += (Real)( 2. * M_PI );

	Real x = (Real)( (theta*res) / (2.*M_PI) );
	Real y = (Real)( (phi*res) / M_PI );
	int _x = (int)floor( x ) , _y = (int)floor( y );

	int idx = _x *  res + _y;

	TriangleIndex tri;
	Real dist = -1;
	Point3D< Real > b;
	for( int j=0 ; j<polygonCells[idx].size() ; j++ )
	{
		TriangleIndex tri = triangles[ polygonCells[idx][j].sourceID ];
		Point3D< Real > triangle[] = { vertices[ tri[0] ] , vertices[ tri[1] ] , vertices[ tri[2] ] };

		Point3D< Real > weights;
		Point3D< Real > q = _NearestPointOnTriangle( p , triangle , &weights[0] );

		Real d2 = Point3D< Real >::SquareNorm( q - p );
		if( dist<0 || d2<dist )
		{
			dist = d2;
			bc.tIndex = polygonCells[idx][j].sourceID;
			bc.weights = weights;
		}
	}
	if( dist<0 ) fprintf( stderr , "[ERROR] Could not find a triangle in cell: %d %d\n" , _x , _y ) , exit( 0 );

	return bc;
}

template< class Real >
Point3D< Real > SourceTargetCorrespondences< Real >::_NearestPointOnTriangle( Point3D< Real > point , const Point3D< Real > triangle[3] , Real* b )
{
	b[0] = b[1] = b[2] = 0;
	Point3D< Real > d[] = { triangle[1]-triangle[0] , triangle[2]-triangle[0] } , p = point - triangle[0] , n = CrossProduct( d[0] , d[1] );

	if( !Length(n) ) return ( triangle[0] + triangle[1] + triangle[2] ) / (Real)3.;


	if     ( Point3D< Real >::Dot( point-triangle[0] , CrossProduct( n , triangle[1]-triangle[0] ) )<0 ){ Point3D< Real > edge[] = { triangle[0] , triangle[1] } ; return _NearestPointOnEdge( point , edge , b[0] , b[1] ); }
	else if( Point3D< Real >::Dot( point-triangle[1] , CrossProduct( n , triangle[2]-triangle[1] ) )<0 ){ Point3D< Real > edge[] = { triangle[1] , triangle[2] } ; return _NearestPointOnEdge( point , edge , b[1] , b[2] ); }
	else if( Point3D< Real >::Dot( point-triangle[2] , CrossProduct( n , triangle[0]-triangle[2] ) )<0 ){ Point3D< Real > edge[] = { triangle[2] , triangle[0] } ; return _NearestPointOnEdge( point , edge , b[2] , b[0] ); }
	else
	{
		// Solve for the minimizer of:
		//                            E(s,t) = || p - s*d[0]-t*d[1] ||^2
		//                                   = || p ||^2 - 2 s < p , d[0] > - 2 t < p , d[1] > + 2 s t < d[0] , d[1] > + s^2 || d[0] ||^2 + t^2 || d[1] ||^2
		//   =>  (0,0) = ( -< p , d[0] > + t < d[0] , d[1] > + s || d[0] ||^2 , -< p , d[1] > + s < d[0] , d[1] > + t || d[1] ||^2
		//            <=> | < p , d[0] > | = | < d[0] , d[0] >   < d[0] , d[1] > | | s |
		//                | < p , d[1] > |   | < d[0] , d[1] >   < d[1] , d[1] > | | t |
		SquareMatrix< Real , 2 > M , M_inverse;
		M(0,0) = Point3D< Real >::SquareNorm( d[0] ) , M(1,0) = M(0,1) = Point3D< Real >::Dot( d[0] , d[1] ) , M(1,1) = Point3D< Real >::SquareNorm( d[1] );
		Real det = M(0,0)*M(1,1) - M(0,1)*M(1,0);
		M_inverse(0,0) = M(1,1) , M_inverse(0,1) = -M(0,1) , M_inverse(1,0) = -M(1,0) , M_inverse(1,1) = M(0,0);
		M_inverse /= det;
		Point2D< Real > st = M_inverse * Point2D< Real >( Point3D< Real >::Dot( p , d[0] ) , Point3D< Real >::Dot( p , d[1] ) );
		// Sanity check
		b[0] = 1. - st[0] - st[1] , b[1] = st[0] , b[2] = st[1];
		Point3D< Real > ret = triangle[0] * ( Real )( 1. - st[0] - st[1] ) + d[0] * st[0] + d[1] * st[1];

		return triangle[0] * ( Real )( 1. - st[0] - st[1] ) + d[0] * st[0] + d[1] * st[1];
	}
}

template< class Real >
Point3D< Real > SourceTargetCorrespondences< Real >::_NearestPointOnEdge( Point3D< Real > point , const Point3D< Real > edge[2] , Real& b0 , Real& b1 )
{
	Point3D< Real > d = edge[1] - edge[0] , p = point - edge[0];
	Real dot = Point3D< Real >::Dot( p , d );
	if( dot<0 ) 
	{
		b0 = 1.;
		return edge[0];
	}
	else if( dot>Point3D< Real >::SquareNorm( d ) ) { 
		b1 = 1.; 
		return edge[1];
	}
	else
	{
		// Solve for the minimizer of:
		//                            E(t) = || p - t*d ||^2
		//                                 = || p ||^2 - 2 t < p , d > + t^2 || d ||^2
		//            =>             0 = -< p , d > + t || d ||^2
		//            <=>    t = < p , d > / || d ||^2
		Real t = dot / Point3D< Real >::SquareNorm( d );
		b0 = 1.-t , b1 = t;
		return edge[0] + d * t;
	}
}

template< typename Real >
std::vector< TriangleIndex > SourceTargetCorrespondences< Real >::_Triangulate( const std::vector< Point3D< Real > > &vertices , const std::vector< std::vector< int > > &polygons )
{
	std::vector< TriangleIndex > triangles;

	int count = 0;
	for( int i=0 ; i<polygons[i].size() ; i++ ) count += (int)polygons[i].size()-2;
	triangles.reserve( count );
//#pragma omp parallel for
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
}
//////////////
// Spectrum //
//////////////
template< typename Real > const unsigned long long Spectrum< Real >::_MAGIC_NUMBER = 0x2019ull;

template< typename Real >
void Spectrum< Real >::set( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , unsigned int dimension , Real offset , bool lump )
{
	// [Definition]
	//	We define the generalized eigensystem (A,B) to be the system A v = \lambda B v
	// [Fact 1]
	// If (v,\lambda) is a solution to the ges (A,B) then (v,\lambda+\epsilon) is a solution to the ges (A+\epsilon*B,B):
	//		(A + \epsilon B) v = (\lambda + \epsilon) B v
	//	<=>	A v + \epsilon B v = \lambda B v + \epsilon B V
	//	<=>                A v = \lambda B v
	// [Fact 2]
	// If (w,\delta) is a solution to the ges (A^{-1},B^{-1}) then (A^{-1}w,1/\delta) is a solution to the ges (A,B):
	//		A^{-1} w = \delta B^{-1} w
	//	<=> v = \delta B^{-1} A v
	//	<=> 1/\delta B v = A v
	// [Corollary]
	// If (w,\delta) is a solution to the ges ( (A+\epsilon*B)^{-1} , B^{-1} ) then:
	//	=> ( (A+\epsilon*B)^{-1} w , 1/\delta ) is a solution to the ges (A+\epsilon*B,B)
	//	=> ( (A+\epsilon*B)^{-1} w , 1\delta-\epsilon ) is a solution to the ges (A,B)

	typedef EigenSolverCholeskyLDLt< Real > Solver;
	struct InverseOperator
	{
		Solver solver;
		InverseOperator( const SparseMatrix< Real , int > &M ) : _M(M) , solver( M ){}
		int rows( void ) const { return (int)_M.rows; }
		int cols( void ) const { return (int)_M.rows; }
		void perform_op( const Real *in , Real *out ) const { const_cast< Solver & >(solver).solve( in , out ); };
	protected:
		const SparseMatrix< Real , int > &_M;
	};

	struct InverseBOperator
	{
		Solver solver;
		InverseBOperator( const SparseMatrix< Real , int > &M ) : _M(M) , solver( M ){}
		int rows( void ) const { return (int)_M.rows; }
		int cols( void ) const { return (int)_M.rows; }
		void solve( const Real *in , Real *out ) const { _M.Multiply( in , out ); }
		void mat_prod( const Real *in , Real *out ) const { const_cast< Solver & >(solver).solve( in , out ); };
	protected:
		const SparseMatrix< Real , int > &_M;
	};

	std::vector< TriangleIndex > _triangles = triangles;
	FEM::RiemannianMesh< Real > mesh( GetPointer( _triangles ) , _triangles.size() );
	mesh.setMetricFromEmbedding( GetPointer( vertices ) );
	mesh.makeUnitArea();
	SparseMatrix< Real , int > M = mesh.template massMatrix< FEM::BASIS_0_WHITNEY >( lump ) , S = mesh.template stiffnessMatrix< FEM::BASIS_0_WHITNEY >();

	// Offset the stiffness matrix so that it becomes positive definite
#pragma omp parallel for
	for( int i=0 ; i<M.rows ; i++ ) for( int j=0 ; j<M.rowSizes[i] ; j++ ) S[i][j].Value += M[i][j].Value * offset;

	InverseOperator op( S );
	InverseBOperator Bop( M );

	Spectra::SymGEigsSolver< Real , Spectra::LARGEST_ALGE , InverseOperator , InverseBOperator , Spectra::GEIGS_REGULAR_INVERSE > geigs( &op , &Bop , dimension , 2*dimension );
	geigs.init();
	int nconv = geigs.compute();
	if( nconv!=dimension ) fprintf( stderr , "[WARNING] Number of converged is not equal to dimension: %d != %d\n" , nconv , dimension );
	Eigen::VectorXd evalues;
	Eigen::MatrixXd evecs;
	if( geigs.info()==Spectra::SUCCESSFUL )
	{
		evalues = geigs.eigenvalues();
		evecs = geigs.eigenvectors();
	}
	else if( geigs.info()==Spectra::NOT_COMPUTED )    fprintf( stderr , "[ERROR] Not computed\n" ) , exit(0);
	else if( geigs.info()==Spectra::NOT_CONVERGING 	) fprintf( stderr , "[ERROR] Not converging\n" ) , exit(0);
	else if( geigs.info()==Spectra::NUMERICAL_ISSUE ) fprintf( stderr , "[ERROR] Numerical issue\n" ) , exit(0);
	else                                              fprintf( stderr , "[ERROR] Failed\n" ) , exit(0);

	_eigenvalues.resize( evecs.cols() );
	_eigenvectors.resize( evecs.cols() );

	for( int i=0 ; i<evecs.cols() ; i++ )
	{
		_eigenvectors[i].resize( vertices.size() );
		_eigenvalues[i] = (Real)1./evalues[i] - offset;
		std::vector< Real > w( vertices.size() );
		for( int j=0 ; j<evecs.rows() ; j++ ) w[j] = evecs(j,i);
		op.perform_op( &w[0] , &_eigenvectors[i][0] );
		Real l2Norm = 0;
#pragma omp parallel for reduction( + : l2Norm )
		for( int j=0 ; j<M.rows ; j++ ) for( int k=0 ; k<M.rowSizes[j] ; k++ ) l2Norm += M[j][k].Value * _eigenvectors[i][j] * _eigenvectors[i][ M[j][k].N ];
		l2Norm = (Real)sqrt( l2Norm );
#pragma omp parallel for
		for( int j=0 ; j<_eigenvectors[i].size() ; j++ ) _eigenvectors[i][j] /= l2Norm;
	}
}
template< typename Real >
void Spectrum< Real >::read( const std::string &fileName )
{
	FILE *fp = fopen( fileName.c_str() , "rb" );
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open file for reading: %s\n" , fileName.c_str() ) , exit( 0 );
	/// Read simple header
	// Read Magic Number
	unsigned long long magicNum;
	if( fread( &magicNum , sizeof(unsigned long long) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] Failed to read magic number\n" ) , exit( 0 );
	if( magicNum!=_MAGIC_NUMBER ) fprintf( stderr , "[ERROR] Bad magic number: %s\n" , fileName.c_str() ) , exit( 0 );

	// type size
	unsigned int typeSize;
	if( fread( &typeSize , sizeof(unsigned int) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] Failed to read type size\n" ) , exit( 0 );
	bool needCasting = typeSize != sizeof(Real);
	if( needCasting ) fprintf( stderr , "[WARNING] Types sizes don't match: %d != %d\n" , (int)typeSize , (int)sizeof(Real) );
	// Num of eigenvalues
	unsigned int numOfValues;
	if( fread( &numOfValues , sizeof(unsigned int) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] Failed to read number of values\n" ) , exit( 0 );
	_eigenvalues.resize( numOfValues );
	_eigenvectors.resize( numOfValues );
	// Dimensions
	unsigned int dimension;
	if( fread( &dimension , sizeof(unsigned int) , 1,  fp )!=1 ) fprintf( stderr , "[ERROR] Failed to read dimension\n" ) , exit( 0 );
	for( int i=0 ; i<_eigenvectors.size() ; i++ ) _eigenvectors[i].resize( dimension );
	/// Content
	if( needCasting )
	{
		// Eigenvalues
		unsigned char *tempMemory = new unsigned char[ typeSize*numOfValues ];
		if( fread( tempMemory , typeSize , numOfValues , fp )!=numOfValues ) fprintf( stderr , "[ERROR] Failed to read values\n" ) , exit( 0 );
#pragma omp parallel for
		for( int i=0 ; i<_eigenvalues.size() ; i++ )
		{
			switch( typeSize )
			{
				case sizeof(float):  _eigenvalues[i] = (Real)reinterpret_cast< float       * >(tempMemory)[i] ; break;
					case sizeof(double): _eigenvalues[i] = (Real)reinterpret_cast< double      * >(tempMemory)[i] ; break;
					default: fprintf( stderr , "[ERROR] Could not determine type from size: %d\n" , (int) typeSize ) , exit( 0 );
			}
		}
		delete[] tempMemory;
		// Eigenvectors
		tempMemory = new unsigned char[ typeSize*dimension ];
		for( int i=0 ; i<_eigenvectors.size() ; i++ )
		{
			if( fread( tempMemory , typeSize , dimension , fp )!=dimension ) fprintf( stderr , "[ERROR] Failed to read eigenvectors\n" ) , exit( 0 );
#pragma omp parallel for
			for( int j=0 ; j<_eigenvectors[i].size() ; j++ )
			{
				switch(typeSize )
				{
					case sizeof(float):  _eigenvectors[i][j] = (Real)reinterpret_cast< float       * >(tempMemory)[j] ; break;
						case sizeof(double): _eigenvectors[i][j] = (Real)reinterpret_cast< double      * >(tempMemory)[j] ; break;
						default: fprintf( stderr , "[ERROR] Could not determine type from size: %d\n" , (int) typeSize ) , exit( 0 );
				}
			}
		}
		delete[] tempMemory;
	}
	else
	{
		// Eigenvalues
		if( fread( _eigenvalues.data() , sizeof(Real) , numOfValues , fp )!=numOfValues ) fprintf( stderr , "[ERROR] Failed to read eigenvalues\n" ) , exit( 0 );
		// Eigenvectors
		for( int i=0 ; i<_eigenvectors.size() ; i++ ) if( fread( _eigenvectors[i].data() , sizeof(Real) , dimension , fp )!=dimension ) fprintf( stderr , "[ERROR] Failed to read eigenvectors\n" ) , exit( 0 );
	}
	fclose( fp );
}
template< typename Real >
void Spectrum< Real >::write( const std::string &fileName ) const
{
	FILE *fp = fopen( fileName.c_str() , "wb" );
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open file for writing: %s\n" , fileName.c_str() ) , exit( 0 );
	// Write simple header
	// Write Magic Number
	fwrite( &_MAGIC_NUMBER , sizeof(unsigned long long) , 1 , fp );

	// type size
	unsigned int typeSize = sizeof(Real);
	fwrite( &typeSize , sizeof(unsigned int) , 1 , fp );
	// Num of eigenvalues
	unsigned int numOfValues = (unsigned int)_eigenvectors.size();
	fwrite( &numOfValues , sizeof(unsigned int) , 1 , fp );
	// Dimensions
	unsigned int dimension = (unsigned int)_eigenvectors[0].size();
	fwrite( &dimension , sizeof(unsigned int) , 1,  fp );
	// Eigenvalues
	fwrite( _eigenvalues.data() , sizeof(Real) , numOfValues , fp );
	// Eigenvectors
	for( int i=0 ; i<_eigenvectors.size() ; i++ ) fwrite( _eigenvectors[i].data() , sizeof(Real) , dimension , fp );
	fclose( fp );
}

template< typename Real >
Real Spectrum< Real >::biharmonicDistance( unsigned int i , unsigned int j ) const
{
	Real distance = (Real)0;
	for( unsigned int k=1 ; k<_eigenvectors.size() ; k++ )
	{
		auto& v = _eigenvectors[k];
		Real temp = (Real) pow( v[i]-v[j] , 2) / (Real)pow( _eigenvalues[k] , 2 );
		distance += temp;
	}
	return (Real)sqrt( distance );
}
template< typename Real >
Real Spectrum< Real >::biharmonicDistance( unsigned int i , TriangleIndex tri , Point3D< Real> weights ) const
{
	Real distance = (Real)0;
	for( unsigned int k=1 ; k<_eigenvectors.size() ; k++ )
	{
		auto& v = _eigenvectors[k];
		Real v1 = v[i] , v2 = v[ tri[0] ] * weights[0] + v[ tri[1] ] * weights[1] + v[ tri[2] ] * weights[2];
		Real temp = (Real) pow( v1-v2 , 2 ) / (Real)pow( _eigenvalues[k] , 2 );
		distance += temp;
	}
	return (Real)sqrt( distance );
}

#endif // SOURCE_TO_TARGET_CORRESPONDENCE_INCLUDED