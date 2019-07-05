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

#include <typeinfo>
#include <algorithm>
#include <unordered_set>
#include <random>
#include "Misha/Timer.h"
#include "Misha/Fourier.h"
#include "Misha/SphericalSignals.h"

#define ALIGN_INVERSION_TO_CENTER

namespace SphericalGeometry
{
	inline unsigned long long Key( int v1 , int v2 )
	{
		if( v1<v2 ) return ( (unsigned long long)v1 )<<32 | ( (unsigned long long )v2);
		else        return ( (unsigned long long)v2 )<<32 | ( (unsigned long long )v1);
	}

	template<>
	inline void AddAtomic( float &a , const float &b )
	{
#pragma omp atomic
		a += b;
	}
	template<>
	inline void AddAtomic( double &a , const double &b )
	{
#pragma omp atomic
		a += b;
	}
	template< typename Data >
	inline void AddAtomic( Point2D< Data > &a , const Point2D< Data > &b )
	{
		for( int d=0 ; d<3; d++ ) AddAtomic( a[d] , b[d] );
	}
	template< typename Data >
	inline void AddAtomic( Point3D< Data > &a , const Point3D< Data > &b )
	{
		for( int d=0 ; d<3; d++ ) AddAtomic( a[d] , b[d] );
	}
	template< typename Data , unsigned int Dim >
	inline void AddAtomic( Point< Data , Dim > &a , const Point< Data , Dim > &b )
	{
		for( unsigned int d=0 ; d<Dim ; d++ ) AddAtomic( a[d] , b[d] );
	}
	template< typename Data , unsigned int Dim >
	inline void AddAtomic( SquareMatrix< Data , Dim > &a , const SquareMatrix< Data , Dim > &b )
	{
		for( unsigned int d1=0 ; d1<Dim ; d1++ ) for( int d2=0 ; d2<Dim ; d2++ ) AddAtomic( a(d1,d2) , b(d1,d2) );
	}

	template< typename Real , typename Data >
	inline void AddAtomic( VectorSpaceDirectSumElement< Real , Data , Real > &a , const VectorSpaceDirectSumElement< Real , Data , Real > &b )
	{
		AddAtomic( std::get< 0 >(a) , std::get< 0 >(b) );
		AddAtomic( std::get< 1 >(a) , std::get< 1 >(b) );
	}

	template< typename Data >
	inline void AddAtomic( Data &a , const Data &b )
	{
#pragma message( "[WARNING] This is very inefficient" )
#pragma omp critical ( add_atomic_critical )
		{
			a += b;
		}
	}
}

template< class Real >
SquareMatrix< Real , 3 > SphericalGeometry::Correlate( SphericalGrid< Real >& source , SphericalGrid< Real >& target , Real& error , bool gradientDomain , bool invert )
{
	return Correlate( 1 , &source , &target , error , gradientDomain , invert );
}

template< class Real >
SquareMatrix< Real , 3 > SphericalGeometry::Correlate( unsigned int gridNum , SphericalGrid< Real > source[] , SphericalGrid< Real > target[] , Real& error , bool gradientDomain , bool invert )
{
	FourierKeySO3< Real> keySO3; 
	HarmonicTransform< Real > hForm; 
	WignerTransform< Real > wForm; 
	FourierKeyS2<Real> sourceKey, targetKey; 
	RotationGrid<Real> euler; 
	Real matrix[3][3];
	for( unsigned int g=0 ; g<gridNum ; g++ ) if( source[g].resolution()!=target[g].resolution() || source[g].resolution()!=source[0].resolution() ) fprintf( stderr , "[ERROR] Resolutions differ\n" ) , exit( 0 );

	int resolution = source[0].resolution();

	// Allocate memory 
	if( !keySO3.resize( resolution ) ) fprintf( stderr , "[ERROR] Correlate: Could not allocate key: %d\n" , resolution ) , exit( 0 ); 
	for( int i=0 ; i<keySO3.bandWidth(); i++ ) for( int j=0 ; j<=i ; j++ ) for( int k=0 ; k<=i ; k++ ) keySO3(i,j,k) = 0;

	for( unsigned int g=0 ; g<gridNum ; g++ )
	{
		hForm.ForwardFourier( source[g] , sourceKey );
		hForm.ForwardFourier( target[g] , targetKey );

		if( gradientDomain )
#pragma omp parallel for
			for( int i=0 ; i<sourceKey.bandWidth() ; i++ )
			{
				Real scl = (Real)sqrt( i*(i+1) );
				for( int j=0 ; j<=i ; j++ ) sourceKey(i,j) *= scl , targetKey(i,j) *= scl;
			}

		// Compute the aligning rotation
#pragma omp parallel for
		for( int i=0 ; i<keySO3.bandWidth(); i++ ) for( int j=0 ; j<=i ; j++ ) for( int k=0 ; k<=i ; k++ ) 
		{
			// set each Wigner-D coefficient to be the cross multiplication 
			// within each frequency 
			int sign = 1;
			if( invert && (i%2) )sign = -1;
			keySO3(i,j,k) += sourceKey(i,k) * targetKey(i,j).conjugate() * sign;
			if( k ) keySO3(i,j,-k) += sourceKey(i,k).conjugate() * targetKey(i,j).conjugate() * sign;
		}
	}

	// Take the inverse Wigner D transform 
	wForm.InverseFourier( keySO3 , euler );

	Real max = 0;
	int maxcoords[] = { 0 , 0 , 0 };
	for( int x=0 ; x<euler.resolution() ; x++ ) for( int y=0 ; y<euler.resolution() ; y++ ) for( int z=0 ; z<euler.resolution() ; z++ ) if( euler(x,y,z)>max )
	{
		max = euler(x,y,z);
		maxcoords[0] = x;
		maxcoords[1] = y;
		maxcoords[2] = z; 
	}

	euler.setCoordinates( maxcoords[0] , maxcoords[1] , maxcoords[2] , matrix );
	SquareMatrix< Real , 3 > m;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) m(i,j) = matrix[i][j];
	error = ( sourceKey.squareNorm() + targetKey.squareNorm() - 2.*max ) / ( sourceKey.squareNorm() + targetKey.squareNorm() );
	if( invert ) return -m;
	else         return  m;
}

template< typename Real >
void SphericalGeometry::SetCorrelation( unsigned int gridNum , SphericalGrid< Real > source[] , SphericalGrid< Real > target[] , bool gradientDomain , bool invert , RotationGrid< Real > &correlation )
{
	FourierKeySO3< Real> keySO3; 
	HarmonicTransform< Real > hForm; 
	WignerTransform< Real > wForm; 
	FourierKeyS2<Real> sourceKey, targetKey; 
	for( unsigned int g=0 ; g<gridNum ; g++ ) if( source[g].resolution()!=target[g].resolution() || source[g].resolution()!=source[0].resolution() ) fprintf( stderr , "[ERROR] Resolutions differ\n" ) , exit( 0 );

	int resolution = source[0].resolution();

	// Allocate memory 
	if( !keySO3.resize( resolution ) ) fprintf( stderr , "[ERROR] Correlate: Could not allocate key: %d\n" , resolution ) , exit( 0 ); 
	for( int i=0 ; i<keySO3.bandWidth(); i++ ) for( int j=0 ; j<=i ; j++ ) for( int k=0 ; k<=i ; k++ ) keySO3(i,j,k) = 0;

	for( unsigned int g=0 ; g<gridNum ; g++ )
	{
		hForm.ForwardFourier( source[g] , sourceKey );
		hForm.ForwardFourier( target[g] , targetKey );
		if( gradientDomain )
#pragma omp parallel for
			for( int i=0 ; i<sourceKey.bandWidth() ; i++ )
			{
				Real scl = (Real)sqrt( i*(i+1) );
				for( int j=0 ; j<=i ; j++ ) sourceKey(i,j) *= scl , targetKey(i,j) *= scl;
			}

		// Compute the aligning rotation
#pragma omp parallel for
		for( int i=0 ; i<keySO3.bandWidth(); i++ ) for( int j=0 ; j<=i ; j++ ) for( int k=0 ; k<=i ; k++ ) 
		{
			// set each Wigner-D coefficient to be the cross multiplication 
			// within each frequency 
			int sign = 1;
			if( invert && (i%2) ) sign = -1;
			keySO3(i,j,k) += sourceKey(i,k) * targetKey(i,j).conjugate() * sign;
			if( k ) keySO3(i,j,-k) += sourceKey(i,k).conjugate() * targetKey(i,j).conjugate() * sign;
		}
	}

	// Take the inverse Wigner D transform 
	wForm.InverseFourier( keySO3 , correlation );
}

template< class Real >
unsigned int SphericalGeometry::TopNRotations( unsigned int N , unsigned int gridNum , SphericalGrid< Real > source[] , SphericalGrid< Real > target[] , bool gradientDomain , bool invert , Real separation , Real correlationFraction, std::vector< Real > &correlationValues , std::vector< SquareMatrix< Real , 3 > > &rotations )
{
	correlationValues.resize( 0 );
	rotations.resize( 0 );

	RotationGrid< Real > correlation;
	SetCorrelation( gridNum , source , target , gradientDomain , invert , correlation );

	// Find the top correlation values
	auto RotationMatrix = [&]( int x , int y , int z )
	{
		Real _R[3][3];
		correlation.setCoordinates( x , y , z , _R );
		SquareMatrix< Real , 3 > R;
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) R(i,j) = _R[i][j];
		return R;
	};
	auto ValidRotation = [&]( SquareMatrix< Real , 3 > R , unsigned int n , Real correlation )
	{
		for( unsigned int i=0 ; i<n ; i++ ) if( SquareMatrix< Real , 3 >::SquareNorm( rotations[i].inverse() * R - SquareMatrix< Real , 3 >::Identity() )<separation ) return false;
		if( n && correlation<correlationValues[0]*correlationFraction ) return false;
		return true;
	};

	for( unsigned int n=0 ; n<N ; n++ )
	{
		bool found = false;
		Real maxCorrelation = 0;
		SquareMatrix< Real , 3 > R;
		for( int x=0 ; x<correlation.resolution() ; x++ ) for( int y=0 ; y<correlation.resolution() ; y++ ) for( int z=0 ; z<correlation.resolution() ; z++ ) if( correlation(x,y,z)>maxCorrelation )
		{
			SquareMatrix< Real , 3 > _R = RotationMatrix( x , y , z );
			if( correlationFraction<0 || ValidRotation( _R , n , correlation(x,y,z) ) )
			{
				maxCorrelation = correlation( x , y , z );
				R = _R;
				found = true;
			}
		}
		if( !found ) return n;
		correlationValues.push_back( maxCorrelation );
		rotations.push_back( invert ? -R : R );
	}
	return N;
}

///////////////////////////////////////////
// SphericalGeometry::SphericalInversion //
///////////////////////////////////////////
template< class Real >
Point3D< Real > SphericalGeometry::SphericalInversion< Real >::operator() ( Point3D< Real > p ) const
{
#ifdef ALIGN_INVERSION_TO_CENTER
	p += center;
#else // !ALIGN_INVERSION_TO_CENTER
	p -= center;
#endif // ALIGN_INVERSION_TO_CENTER
	p /= Point3D< Real >::SquareNorm( p );
	p *= 1 - Point3D< Real >::SquareNorm( center );
#ifdef ALIGN_INVERSION_TO_CENTER
	p += center;
#else // !ALIGN_INVERSION_TO_CENTER
	p -= center;
#endif // ALIGN_INVERSION_TO_CENTER
	return p;
}

///////////////////////////////////////////////////////
// SphericalGeometry::FractionalLinearTransformation //
///////////////////////////////////////////////////////
template< class Real >
SphericalGeometry::FractionalLinearTransformation< Real >::FractionalLinearTransformation( SphericalInversion< Real > si )
{
	matrix = _ITransformation( StereographicProjection( si(_ZERO) ) , StereographicProjection( si(_INFINITY) ) , StereographicProjection( si(_ONE) ) );
}
template< class Real >
SphericalGeometry::FractionalLinearTransformation< Real >::FractionalLinearTransformation( SquareMatrix< Real , 3 > m )
{
	matrix = _ITransformation( StereographicProjection( m * _ZERO ) , StereographicProjection( m * _INFINITY ) , StereographicProjection( m * _ONE ) );
}
template< class Real >
SphericalGeometry::FractionalLinearTransformation< Real >::FractionalLinearTransformation( const std::pair< Point2D< Real > , Point2D< Real > > pairs[3] )
{
	matrix = _ITransformation( pairs[0].second , pairs[1].second , pairs[2].second ) * _Transformation( pairs[0].first , pairs[1].first , pairs[2].first );
}

template< class Real >
Point2D< Real > SphericalGeometry::FractionalLinearTransformation< Real >::operator()( Point2D< Real > p ) const
{
	std::complex< Real > c( p[0] , p[1] );
	c = ( c * matrix(0,0) + matrix(1,0) ) / ( c * matrix(0,1) + matrix(1,1) );
	return Point2D< Real >( c.real() , c.imag() );
}
template< class Real >
Point3D< Real > SphericalGeometry::FractionalLinearTransformation< Real >::operator()( Point3D< Real > p ) const
{
	return IStereographicProjection( (*this)( StereographicProjection( p ) ) );
}
template< class Real >
SphericalGeometry::FractionalLinearTransformation< Real > SphericalGeometry::FractionalLinearTransformation< Real >::operator * ( SphericalGeometry::FractionalLinearTransformation< Real > flt ) const
{
	FractionalLinearTransformation _flt = *this;
	return _flt *= flt;
}
template< class Real >
SphericalGeometry::FractionalLinearTransformation< Real >& SphericalGeometry::FractionalLinearTransformation< Real >::operator *= ( SphericalGeometry::FractionalLinearTransformation< Real > flt )
{
	matrix *= flt.matrix;
	return *this;
};
template< class Real >
SphericalGeometry::FractionalLinearTransformation< Real > SphericalGeometry::FractionalLinearTransformation< Real >::inverse( void ) const
{
	FractionalLinearTransformation flt;
	std::complex< Real > det = matrix(0,0) * matrix(1,1) - matrix(1,0) * matrix(0,1);
	flt.matrix(0,0) =  matrix(1,1) / det , flt.matrix(1,1) =  matrix(0,0) / det;
	flt.matrix(1,0) = -matrix(1,0) / det , flt.matrix(0,1) = -matrix(0,1) / det;
	return flt;
}

template< class Real >
SquareMatrix< std::complex< Real > , 2 > SphericalGeometry::FractionalLinearTransformation< Real >::_Transformation( Point2D< Real > p1 , Point2D< Real > p2 , Point2D< Real > p3 )
{
	SquareMatrix< std::complex< Real > , 2 > matrix;
	std::complex< Real > _p1( p1[0] , p1[1] ) , _p2( p2[0] , p2[1] ) , _p3( p3[0] , p3[1] );
	if( std::isfinite( p2[0] ) && std::isfinite( p2[1] ) )
	{
		matrix(0,0) =         ( _p3 - _p2 );
		matrix(1,0) = - _p1 * ( _p3 - _p2 );
		matrix(0,1) =         ( _p3 - _p1 );
		matrix(1,1) = - _p2 * ( _p3 - _p1 );
	}
	else
	{
		matrix(0,0) =     1;
		matrix(1,0) = - _p1;
		matrix(0,1) =     0;
		matrix(1,1) = _p3 - _p1;
	}

	return matrix;
}
template< class Real >
SquareMatrix< std::complex< Real > , 2 > SphericalGeometry::FractionalLinearTransformation< Real >::_ITransformation( Point2D< Real > p1 , Point2D< Real > p2 , Point2D< Real > p3 )
{
	SquareMatrix< std::complex< Real > , 2 > matrix , _matrix = _Transformation( p1 , p2 , p3 );
	std::complex< Real > det = _matrix(0,0) * _matrix(1,1) - _matrix(1,0) * _matrix(0,1);
	matrix(0,0) = _matrix(1,1) , matrix(1,1) = _matrix(0,0);
	matrix(0,1) = -_matrix(0,1) , matrix(1,0) = -_matrix(1,0);
	return matrix / det;
}
template< class Real > const Point3D< Real > SphericalGeometry::FractionalLinearTransformation< Real >::_ZERO    (  0 ,  0 , -1 );
template< class Real > const Point3D< Real > SphericalGeometry::FractionalLinearTransformation< Real >::_INFINITY(  0 ,  0 ,  1 );
template< class Real > const Point3D< Real > SphericalGeometry::FractionalLinearTransformation< Real >::_ONE     (  1 ,  0 ,  0 );

template< class Real >
Real SphericalGeometry::Area( const std::vector< Point3D< Real > > &vertices , const std::vector< int > &polygon )
{
	// For a triangle:
	//		N = ( v1-v0 ) x ( v2-v0 )
	//        = v0 x v1 + v1 x v2 + v2 x v0 
	//        = [ ( v1-v0 ) x ( v1+v0 ) +  ( v2-v1 ) x ( v2+v1 ) +  ( v0-v2 ) x ( v0+v2 ) ] / 2
	Point3D< Real > n;
	for( int i=0 ; i<polygon.size() ; i++ )
	{
		Point3D< Real > b = vertices[ polygon[(i+1)%polygon.size() ] ] + vertices[ polygon[i] ];
		Point3D< Real > e = vertices[ polygon[(i+1)%polygon.size() ] ] - vertices[ polygon[i] ];
		n += Point3D< Real >::CrossProduct( b , e ) / 2;
	}
	return (Real)Length( n ) / 2;
}
template< class Real >
Real SphericalGeometry::Area( const std::vector< Point3D< Real > > &vertices ){ return Area( &vertices[0] , vertices.size() ); }
template< class Real >
Real SphericalGeometry::Area( const Point3D< Real > *vertices , size_t vNum )
{
	// For a triangle:
	//		N = ( v1-v0 ) x ( v2-v0 )
	//        = v0 x v1 + v1 x v2 + v2 x v0 
	//        = [ ( v1-v0 ) x ( v1+v0 ) +  ( v2-v1 ) x ( v2+v1 ) +  ( v0-v2 ) x ( v0+v2 ) ] / 2
	Point3D< Real > n;
	for( int i=0 ; i<vNum ; i++ )
	{
		Point3D< Real > b = vertices[ (i+1)%vNum ] + vertices[i];
		Point3D< Real > e = vertices[ (i+1)%vNum ] - vertices[i];
		n += Point3D< Real >::CrossProduct( b , e ) / 2;
	}
	return (Real)Length( n ) / 2;
}
/////////////////////////////
// SphericalGeometry::Mesh //
/////////////////////////////
template< class Real >
SphericalGeometry::Mesh< Real >::Mesh( const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& sVertices , const std::vector< std::vector< int > > &polygons )
{
	this->vertices = sVertices;
	this->polygons = polygons;
	masses.resize( polygons.size() );
#pragma omp parallel for
	for( int i=0 ; i<polygons.size() ; i++ ) masses[i] = Area( vertices , polygons[i] );
}

template< class Real >
std::vector< TriangleIndex > SphericalGeometry::Mesh< Real >::triangulate( const std::vector< Point3D< Real > > &vertices ) const
{
	std::vector< TriangleIndex > triangles;

	int count = 0;
	for( int i=0 ; i<polygons.size() ; i++ ) count += (int)polygons[i].size()-2;
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
}
template< class Real >
void SphericalGeometry::Mesh< Real >::write( const char* fileName , const std::vector< Point3D< Real > >& vertices , bool binary ) const
{
	std::vector< PlyParametrizedVertex< float , Real > > _vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].point = Point3D< float >( vertices[i] ) , _vertices[i].param = Point3D< float >( this->vertices[i] );
	PlyWritePolygons( fileName , _vertices , polygons , PlyParametrizedVertex< float , Real >::WriteProperties , PlyParametrizedVertex< float , Real >::WriteComponents , binary ? PLY_BINARY_NATIVE : PLY_ASCII , NULL , 0 );
}
template< class Real >
void SphericalGeometry::Mesh< Real >::write( const char* fileName , const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& colors , bool binary ) const
{
	std::vector< PlyParametrizedColorVertex< float , Real > > _vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].point = Point3D< float >( vertices[i] ) , _vertices[i].param = Point3D< float >( this->vertices[i] ) , _vertices[i].color = Point3D< float >( colors[i] );
	PlyWritePolygons( fileName , _vertices , polygons , PlyParametrizedColorVertex< float , Real >::WriteProperties , PlyParametrizedColorVertex< float , Real >::WriteComponents , binary ? PLY_BINARY_NATIVE : PLY_ASCII , NULL , 0 );
}
template< class Real >
void SphericalGeometry::Mesh< Real >::setMasses( const std::vector< Point3D< Real > > &meshVertices )
{
	masses.resize( polygons.size() );
#pragma omp parallel for
	for( int i=0 ; i<polygons.size() ; i++ ) masses[i] = Area( meshVertices , polygons[i] );
}
template< class Real >
void SphericalGeometry::Mesh< Real >::read( const char* fileName , std::vector< Point3D< Real > >& vertices , bool verbose )
{
	int plyType;
	std::vector< PlyParametrizedVertex< float , Real > > _vertices;
	PlyReadPolygons( fileName , _vertices , polygons , PlyParametrizedVertex< float , Real >::ReadProperties , NULL , PlyParametrizedVertex< float , Real >::ReadComponents , plyType , NULL , 0 );
	this->vertices.resize( _vertices.size() ) , vertices.resize( _vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) this->vertices[i] = Point3D< Real >( _vertices[i].param ) , vertices[i] = Point3D< Real >( _vertices[i].point );
	masses.resize( polygons.size() );
#pragma omp parallel for
	for( int i=0 ; i<polygons.size() ; i++ ) masses[i] = Area( vertices , polygons[i] );
}
template< class Real >
bool SphericalGeometry::Mesh< Real >::read( const char* fileName , std::vector< Point3D< Real > >& vertices , std::vector< Point3D< Real > >& colors , bool verbose  )
{
	int plyType;
	bool propertyFlags[ PlyParametrizedColorVertex< float , Real >::ReadComponents ];
	std::vector< PlyParametrizedColorVertex< float , Real > > _vertices;
	PlyReadPolygons( fileName , _vertices , polygons , PlyParametrizedColorVertex< float , Real >::ReadProperties , propertyFlags , PlyParametrizedColorVertex< float , Real >::ReadComponents , plyType , NULL , 0 );
	bool hasColor = (propertyFlags[6]||propertyFlags[9]) && (propertyFlags[7]||propertyFlags[10]) && (propertyFlags[8]||propertyFlags[11]);

	this->vertices.resize( _vertices.size() ) , vertices.resize( _vertices.size() ) , colors.resize( _vertices.size() );
	for( int i=0 ; i<_vertices.size() ; i++ ) this->vertices[i] = Point3D< Real >( _vertices[i].param ) , vertices[i] = Point3D< Real >( _vertices[i].point ) , colors[i] = Point3D< Real >( _vertices[i].color );
	masses.resize( polygons.size() );
#pragma omp parallel for
	for( int i=0 ; i<polygons.size() ; i++ ) masses[i] = Area( vertices , polygons[i] );
	return hasColor;
}

template< class Real >
double SphericalGeometry::Mesh< Real >::area( int p ) const { return SphericalGeometry::Area( vertices , polygons[p] ); }
template< class Real > template< typename F >
double SphericalGeometry::Mesh< Real >::area( int p , F f ) const
{
	std::vector< Point3D< Real > > v( polygons[p].size() );
	for( int i=0 ; i<polygons[p].size() ; i++ ) v[i] = f( vertices[ polygons[p][i] ] );
	return SphericalGeometry::Area( v );
}
template< class Real >
Point3D< Real > SphericalGeometry::Mesh< Real >::center( int p ) const
{
	Point3D< Real > c;
	for( int i=0 ; i<polygons[p].size() ; i++ ) c += vertices[ polygons[p][i] ];
	return c / Length( c );
}
template< class Real > template< typename F >
Point3D< Real > SphericalGeometry::Mesh< Real >::center( int p , F f ) const
{
	Point3D< Real > c;
	for( int i=0 ; i<polygons[p].size() ; i++ ) c += f( vertices[ polygons[p][i] ] );
	return c / Length( c );
}
template< class Real >
Real SphericalGeometry::Mesh< Real >::constant( void ) const
{
	Real c = 0;
	for( int i=0 ; i<masses.size() ; i++ ) c += masses[i];
	return c;
}
template< class Real >
Point3D< Real > SphericalGeometry::Mesh< Real >::center( void ) const
{
	// We are looking for the vector v that minimizes
	// E(v) = \sum_T ( < v , p_T > - f(p_T) )^2 * |T|
	//      = v^t ( \sum_T (p_T) * (p_T)^T * |T| ) v - v^t ( 2 * \sum_T (p_T) * f(p_T) * |T| ) + ...
	// => ( \sum_T (p_T) * (p_T)^T * |T| ) v = \sum_T (p_T) * f(p_T) * |T|
	// Noting that the integral of p * p_T over the sphere is 4*Pi/3 * Id., this gives:
	//    v = 3 / (4*PI) * \sum_T (p_T) * f(p_T) * |T|
	Point3D< Real > c;
	for( int p=0 ; p<polygons.size() ; p++ ) c += center(p) * masses[p]; // = p * ( mA / sA ) * sA;
	return c;
}
template< class Real > template< typename F >
Point3D< Real > SphericalGeometry::Mesh< Real >::center( F f ) const
{
	Point3D< Real > c;
	for( int p=0 ; p<polygons.size() ; p++ ) c += center( p , f ) * masses[p];
	return c;
}
template< class Real > template< unsigned int SHDegree >
Point< Real , SphericalHarmonics::Dimension< SHDegree >() > SphericalGeometry::Mesh< Real >::centerSH( void ) const
{
	static_assert( SHDegree<=SphericalHarmonics::MaxDegree , "[ERROR] Degree exceeds maximum spherical harmonic degree" );
	std::vector< Point< double , SphericalHarmonics::Dimension< SHDegree >() > > projections( omp_get_max_threads() );
#pragma omp parallel for
	for( int t=0 ; t<masses.size() ; t++ )
	{
		int thread = omp_get_thread_num();
		double _p[ SphericalHarmonics::Dimension< SHDegree >() ];
		SphericalHarmonics::HarmonicValues< SHDegree >( Point3D< double >( center(t) ) , _p );
		for( unsigned int i=0 ; i<SphericalHarmonics::Dimension< SHDegree >() ; i++ ) projections[thread][i] += _p[i] * masses[t];
	}
	Point< Real , SphericalHarmonics::Dimension< SHDegree >() > p;
	for( int i=0 ; i<projections.size() ; i++ ) for( unsigned int j=0 ; j<SphericalHarmonics::Dimension< SHDegree >() ; j++ ) p[j] += (Real)projections[i][j];
	return p;
}
template< class Real > template< unsigned int SHDegree >
SquareMatrix< Real , SphericalHarmonics::Dimension< SHDegree >() > SphericalGeometry::Mesh< Real >::dCenterSH( void ) const
{
	static_assert( SHDegree<=SphericalHarmonics::MaxDegree , "[ERROR] Degree exceeds maximum spherical harmonic degree" );
	std::vector< SquareMatrix< double , SphericalHarmonics::Dimension< SHDegree >() > > dProjections( omp_get_max_threads() );
#pragma omp parallel for
	for( int t=0 ; t<masses.size() ; t++ )
	{
		int thread = omp_get_thread_num();
		Point3D< double > _p[ SphericalHarmonics::Dimension< SHDegree >() ];
		SphericalHarmonics::HarmonicGradients< SHDegree >( Point3D< double >( center(t) ) , _p );
		for( unsigned int i=0 ; i<SphericalHarmonics::Dimension< SHDegree >() ; i++ ) for( unsigned int j=0 ; j<=i ; j++ ) dProjections[thread](i,j) += Point3D< Real >::Dot( _p[i] , _p[j] ) * masses[t];
	}
	SquareMatrix< Real , SphericalHarmonics::Dimension< SHDegree >() > d;
	for( int i=0 ; i<dProjections.size() ; i++ ) for( unsigned int j=0 ; j<SphericalHarmonics::Dimension< SHDegree >() ; j++ ) for( unsigned int k=0 ; k<=j ; k++ ) d(j,k) += (Real)dProjections[i](j,k);
	for( unsigned int i=0 ; i<SphericalHarmonics::Dimension< SHDegree >() ; i++ ) for( unsigned int j=0 ; j<=i ; j++ ) d(j,i) = d(i,j);
	return d;
}
template< class Real > template< typename VF >
void SphericalGeometry::Mesh< Real >::advect( VF vf , int steps )
{
#pragma omp parallel for
	for( int i=0 ; i<vertices.size() ; i++ )
	{
		Point3D< Real >& v = vertices[i];
		for( int j=0 ; j<steps ; j++ )
		{
			Point3D< Real > _vf = Point3D< Real >( vf( v ) );
			Real l = Length( _vf );
			if( l ) v = v * cos( l ) + _vf / l * sin( l );
		}
	}
}

template< class Real >
SquareMatrix< Real , 3 > SphericalGeometry::Mesh< Real >::dCenter( void ) const
{
	SquareMatrix< Real , 3 > D;
	for( int p=0 ; p<polygons.size() ; p++ )
	{
		SquareMatrix< Real , 3 > _D = SquareMatrix< Real , 3 >::Identity();
		Point3D< Real > c = center( p );
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) _D(i,j) -= c[i] * c[j];
		D += _D * masses[p];
	}
#ifdef ALIGN_INVERSION_TO_CENTER
	return D * (Real)2;
#else // !ALIGN_INVERSION_TO_CENTER
	return -D * (Real)2;
#endif // ALIGN_INVERSION_TO_CENTER
};
template< class Real > template< typename F >
SquareMatrix< Real , 3 > SphericalGeometry::Mesh< Real >::dCenter( F f ) const
{
	SquareMatrix< Real , 3 > D;
	for( int p=0 ; p<polygons.size() ; p++ )
	{
		SquareMatrix< Real , 3 > _D = SquareMatrix< Real , 3 >::Identity();
		Point3D< Real > c = center( p , f );
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) _D(i,j) -= c[i] * c[j];
		D += _D * masses[p];
	}
#ifdef ALIGN_INVERSION_TO_CENTER
	return D * (Real)2;
#else // !ALIGN_INVERSION_TO_CENTER
	return -D * (Real)2;
#endif // ALIGN_INVERSION_TO_CENTER
};

template< class Real > template< typename F >
SphericalGeometry::Mesh< Real >& SphericalGeometry::Mesh< Real >::operator *= ( F f )
{
#pragma omp parallel for
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = f( vertices[i] );
	return *this;
}

template< class Real >
void SphericalGeometry::Mesh< Real >::makeUnitMass( void )
{
	Real m = constant();
	for( int i=0 ; i<masses.size() ; i++ ) masses[i] /= m;
}

template< class Real >
SphericalGeometry::FractionalLinearTransformation< Real > SphericalGeometry::Mesh< Real >::normalizer( int iters , double cutOff , bool gaussNewton , bool verbose ) const
{
	FractionalLinearTransformation< Real > flt;
	for( int i=0 ; i<iters ; i++ )
	{
		Point3D< Real > c = center( flt );
		if( verbose ) printf( "Moebius register[%d]: %8.2e\t(%8.1e %8.1e %8.1e)\n" , i , Length( c ) , c[0] , c[1] , c[2] );
		if( Point3D< Real >::SquareNorm( c )<=cutOff*cutOff ) return flt;
		SquareMatrix< Real , 3 > D = dCenter( flt );
		Point3D< Real > center;
		if( gaussNewton )
		{
			SquareMatrix< Real , 3 > D2 = ( D.transpose() * D ).inverse();
			center = - D2 * ( D * c );
		}
		else center = - D * c * 2;
		// [WARNING] Accumulating double inversions through fractional linear transformations can result in a loss of precision
		flt = FractionalLinearTransformation< Real >( SphericalInversion< Real >( center ) ) * flt;
	}
	if( verbose )
	{
		Point3D< Real > c = center( flt );
		printf( "\tMoebius register[%d]: %8.2e\t(%8.1e %8.1e %8.1e)\n" , iters , Length( c ) , c[0] , c[1] , c[2] );
	}
	return flt;
}

template< class Real >
int SphericalGeometry::Mesh< Real >::normalize( int iters , double cutOff , bool gaussNewton , const CenterToInversion& c2i , int verbose )
{
	Timer timer;
	for( int i=0 ; i<iters ; i++ )
	{
		Point3D< Real > c = center( );
		if( verbose>1 ) printf( "M-center[%d]: %8.2e\t(%8.1e %8.1e %8.1e)      \r" , i , Length( c ) , c[0] , c[1] , c[2] );
		if( Point3D< Real >::SquareNorm( c )<=cutOff*cutOff )
		{
			if( verbose==1 )
			{
				Point3D< Real > c = center( );
				printf( "\nM-center[%d] %.2f(s): %8.2e\t(%8.1e %8.1e %8.1e)\n" , i , timer.elapsed() , Length( c ) , c[0] , c[1] , c[2] );
			}
			return i;
		}
		SquareMatrix< Real , 3 > D = dCenter( );
		if( gaussNewton )
		{
			SquareMatrix< Real , 3 > D2 = ( D.transpose() * D ).inverse();
			c = - D2 * ( D * c );
		}
		else c = - D * c * 2;
		SphericalInversion< Real > inv = c2i( *this , c );
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = inv( vertices[i] );
	}
	if( verbose )
	{
		if( iters ) printf( "\n" );
		Point3D< Real > c = center( );
		printf( "M-center[%d] %.2f(s): %8.2e\t(%8.1e %8.1e %8.1e)\n" , iters , timer.elapsed() , Length( c ) , c[0] , c[1] , c[2] );
	}
	return iters;
}

template< class Real >
template< unsigned int SHDegree >
int SphericalGeometry::Mesh< Real >::normalizeSH( int iters , int advectionSteps , Real advectionStepSize , double cutOff , bool gaussNewton , int verbose )
{
	static const unsigned int Dim = SphericalHarmonics::Dimension< SHDegree >() - 1;
	auto SubPoint = [&]( const Point< Real , SphericalHarmonics::Dimension< SHDegree >() >& p )
	{
		Point< Real , Dim > _center , _p;
		for( int i=0 ; i<Dim ; i++ ) _p[i] = p[i+1];
		return _p;
	};
	auto SubMatrix = [&]( const SquareMatrix< Real , SphericalHarmonics::Dimension< SHDegree >() >& M )
	{
		SquareMatrix< Real , Dim > _M;
		for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) _M(i,j) = M(i+1,j+1);
		return _M;
	};
	
	Timer timer;
	for( int i=0 ; i<iters ; i++ )
	{
		Point< Real , SphericalHarmonics::Dimension< SHDegree >() > center;
		Point< Real , Dim > c = SubPoint ( centerSH< SHDegree >() );

		if( verbose>1 )
		{
			printf( "M-center[%d]: %8.2e\t" , i , sqrt( Point< Real , Dim >::SquareNorm( c ) ) );
			printf( "( " );
			for( int i=0 ; i<Dim ; i++ ) printf( " %8.1e" , c[i] );
			printf( " )\n" );
		}
		if( Point< Real , Dim >::SquareNorm( c )<=cutOff*cutOff )
		{
			if( verbose==1 )
			{
				Point< Real , Dim > c = SubPoint ( centerSH< SHDegree >() );
				printf( "M-center[%d] %.2f(s): %8.2e\t" , i , timer.elapsed() , sqrt( Point< Real , Dim >::SquareNorm( c ) ) );
				printf( "( " );
				for( int i=0 ; i<Dim ; i++ ) printf( " %8.1e" , c[i] );
				printf( " )\n" );
			}
			return i;
		}

		SquareMatrix< Real , Dim > D = SubMatrix( dCenterSH< SHDegree >() );
		if( gaussNewton )
		{
			SquareMatrix< Real , Dim > D2 = ( D.transpose() * D ).inverse();
			c = - D2 * ( D * c );
		}
		else c = - D * c * 2;

		for( int i=0 ; i<Dim ; i++ ) center[i+1] = c[i];

		advect( [&]( Point3D< Real > p ){ return SphericalHarmonics::Gradient< SHDegree >( center*advectionStepSize , p ); } , advectionSteps );
	}
	if( verbose )
	{
		Point< Real , Dim > c = SubPoint ( centerSH< SHDegree >() );
		printf( "M-center[%d] %.2f(s): %8.2e\t" , iters , timer.elapsed() , sqrt( Point< Real , Dim >::SquareNorm( c ) ) );
		printf( "( " );
		for( int i=0 ; i<Dim ; i++ ) printf( " %8.1e" , c[i] );
		printf( " )\n" );
	}
	return iters;
}

////////////////////////////////
// SphericalGeometry::Polygon //
////////////////////////////////
template< class Real > template< class VertexData >
void SphericalGeometry::Polygon< Real >::_split( Point3D< Real > pNormal , Real pOffset , int pIdx , Polygon &back , Polygon &front , std::vector< Point3D< Real > >& vertices , std::vector< PlanarSources > &pSources , std::vector< VertexData >* vData , std::unordered_map< unsigned long long , int >& vMap ) const
{
	std::vector< Real > values( _vertices.size() );
	bool frontSet = false , backSet = false;
	front.sourceID = back.sourceID = sourceID;

	// Evaluate the plane's function at the vertices and mark if front/back vertices have been found
	for( int i=0 ; i<size() ; i++ )
	{
		values[i] = Point3D< Real >::Dot( vertices[ _vertices[i] ] , pNormal ) - pOffset;

		// Don't split on edges near the pole
		if( pIdx>=0 && pSources[ _vertices[i] ].full() ) values[i] = 0;

		backSet |= ( values[i]<0 ) , frontSet |= ( values[i]>0 );
	}

	// If no vertices are in front or no vertices are behind, life is easy
	if( !frontSet ){  back._vertices = _vertices ;  back.mass = mass ; return; }
	if( !backSet  ){ front._vertices = _vertices ; front.mass = mass ; return; }

	enum
	{
		PLANE_STATE_BACK ,
		PLANE_STATE_FRONT ,
		PLANE_STATE_BOUNDARY
	};
	auto State = []( Real value )
	{
		if     ( value<0 ) return PLANE_STATE_BACK;
		else if( value>0 ) return PLANE_STATE_FRONT;
		else               return PLANE_STATE_BOUNDARY;
	};
	struct _VertexIndex
	{
		int idx , state;
		_VertexIndex( int i , int s ) : idx(i) , state(s) {}
	};

	// Compute the indices of the refined polygons and their associated states
	std::vector< _VertexIndex > vertexIndices;
	for( int i0=0 ; i0<size() ; i0++ )
	{
		int i1 = (i0+1)%size();
		vertexIndices.push_back( _VertexIndex( _vertices[i0] , State( values[i0] ) ) );
		// If the edge crosses the plane, we need to create a new vertex
		if( values[i0]*values[i1]<0 )
		{
			int vIdx;
			unsigned long long key = Key( _vertices[i0] , _vertices[i1] );
			if( vMap.find( key )==vMap.end() )
			{
				Real t = values[i0] / ( values[i0] - values[i1] );
				Point3D< Real > v = vertices[ _vertices[i1] ]*t + vertices[ _vertices[i0] ]*(Real)(1.-t);
				v /= (Real)Length( v );
				vIdx = (int)vertices.size();
				vertices.push_back( v );
				if( vData ) vData->push_back( (*vData)[ _vertices[i1] ]*t + (*vData)[ _vertices[i0] ]*(Real)(1.-t) );
				{
					PlanarSources ps;
					ps.addIndex( PlanarSources::SharedIndex( pSources[ _vertices[i0] ] , pSources[ _vertices[i1] ] ) );
					ps.addIndex( pIdx );
					pSources.push_back( ps );
				}
				vMap[key] = vIdx;
			}
			else vIdx = vMap[key];
			vertexIndices.push_back( _VertexIndex( vIdx , PLANE_STATE_BOUNDARY ) );
		}
	}

	for( int i=0 ; i<vertexIndices.size() ; i++ )
	{
		if( vertexIndices[i].state!=PLANE_STATE_BACK ) front._vertices.push_back( vertexIndices[i].idx );
		if( vertexIndices[i].state!=PLANE_STATE_FRONT ) back._vertices.push_back( vertexIndices[i].idx );
	}

	Real a0 = back.area( vertices ) , a1 = front.area( vertices );
	Real a = a0 + a1;
	if( a ) a0 /= a , a1 /= a;
	else a0 = a1 = (Real)1./2;
	back.mass = mass * a0 , front.mass = mass * a1;
}
template< class Real >
double SphericalGeometry::Polygon< Real >::area( const std::vector< Point3D< Real > >& vertices ) const { return SphericalGeometry::Area( vertices , _vertices ); }

////////////////////////////////////////
// SphericalGeometry::ClippedTriangle //
////////////////////////////////////////
template< class Real >
Real SphericalGeometry::ClippedTriangle< Real >::area( void ) const
{
	Point3D< Real > n;
	for( int i=0 ; i<_vNum ; i++ )
	{
		Point3D< Real > b = _vertices[ (i+1)%_vNum ] + _vertices[i];
		Point3D< Real > e = _vertices[ (i+1)%_vNum ] - _vertices[i];
		n += Point3D< Real >::CrossProduct( b , e ) / 2;
	}
	return (Real)Length( n ) / 2;
}

//////////////////////////////////////////////
// SphericalGeometry::ClippedPrimalTriangle //
//////////////////////////////////////////////
template< class Real , typename Data >
void SphericalGeometry::ClippedPrimalTriangle< Real , Data >::split( Point3D< Real > pNormal , Real pOffset , ClippedPrimalTriangle &back , ClippedPrimalTriangle &front , double eps ) const
{
	enum
	{
		PLANE_STATE_BACK ,
		PLANE_STATE_FRONT ,
		PLANE_STATE_BOUNDARY
	};
	Real values[ _MAX_VERTICES ];
	bool frontSet = false , backSet = false;

	// Evaluate the plane's function at the vertices and mark if front/back vertices have been found
	for( int i=0 ; i<size() ; i++ )
	{
		values[i] = Point3D< Real >::Dot( _vertices[i] , pNormal ) - pOffset;
		if( values[i]>=-eps && values[i]<=eps ) values[i] = 0;
		backSet |= ( values[i]<0 ) , frontSet |= ( values[i]>0 );
	}

	// If no vertices are in front or no vertices are behind, life is easy
	if( !frontSet ){  back = *this ; return; }
	if( !backSet  ){ front = *this ; return; }

	auto State = []( Real value , double eps )
	{
		if     ( value<0 ) return PLANE_STATE_BACK;
		else if( value>0 ) return PLANE_STATE_FRONT;
		else               return PLANE_STATE_BOUNDARY;
	};
	struct _VertexIndex
	{
		Point3D< Real > vertex;
		Data data;
		int state;
		_VertexIndex( void ){}
		_VertexIndex( Point3D< Real > v , Data d , int s ) : vertex(v) , data(d) , state(s) {}
	};

	// Compute the indices of the refined polygons and their associated states
	_VertexIndex vertexIndices[ 2*_MAX_VERTICES ];
	unsigned int vertexIndexCount = 0;
	for( int i0=0 ; i0<size() ; i0++ )
	{
		int i1 = (i0+1)%size();
		vertexIndices[ vertexIndexCount++ ] = _VertexIndex( _vertices[i0] , _data[i0] , State( values[i0] , eps ) );

		// If the edge crosses the plane, we need to create a new vertex
		if( values[i0]*values[i1]<0 )
		{
			Real t = values[i0] / ( values[i0] - values[i1] );
			Point3D< Real > v = _vertices[i1]*t + _vertices[i0]*(Real)(1.-t);
			Data d = _data[i1]*t + _data[i0]*(Real)(1.-t);
			v /= (Real)Length( v );
			vertexIndices[ vertexIndexCount++ ] = _VertexIndex( v , d , PLANE_STATE_BOUNDARY );
		}
	}
	for( unsigned int i=0 ; i<vertexIndexCount ; i++ )
	{
		if( vertexIndices[i].state!=PLANE_STATE_BACK ) front.push_back( vertexIndices[i].vertex , vertexIndices[i].data );
		if( vertexIndices[i].state!=PLANE_STATE_FRONT ) back.push_back( vertexIndices[i].vertex , vertexIndices[i].data );
	}
}

template< class Real , typename Data >
Data SphericalGeometry::ClippedPrimalTriangle< Real , Data >::integral( void ) const
{
	Data d{};
	for( int i=1 ; i<_vNum-1 ; i++ )
	{
		Data _d = ( _data[0] + _data[i] + _data[i+1] ) / (Real)3;
		Real _a = Point3D< Real >::Length( Point3D< Real >::CrossProduct( _vertices[i] - _vertices[0] , _vertices[i+1] - _vertices[0] ) ) / (Real)2;
		d += _d * _a;
	}
	return d;
}

////////////////////////////////////////////
// SphericalGeometry::ClippedDualTriangle //
////////////////////////////////////////////
template< class Real , typename Data >
void SphericalGeometry::ClippedDualTriangle< Real , Data >::split( Point3D< Real > pNormal , Real pOffset , ClippedDualTriangle &back , ClippedDualTriangle &front , double eps ) const
{
	enum
	{
		PLANE_STATE_BACK ,
		PLANE_STATE_FRONT ,
		PLANE_STATE_BOUNDARY
	};
	Real values[ _MAX_VERTICES ];
	bool frontSet = false , backSet = false;

	// Evaluate the plane's function at the vertices and mark if front/back vertices have been found
	for( int i=0 ; i<size() ; i++ )
	{
		values[i] = Point3D< Real >::Dot( _vertices[i] , pNormal ) - pOffset;
		if( values[i]>=-eps && values[i]<=eps ) values[i] = 0;
		backSet |= ( values[i]<0 ) , frontSet |= ( values[i]>0 );
	}

	// If no vertices are in front or no vertices are behind, life is easy
	if( !frontSet ){  back = *this ; return; }
	if( !backSet  ){ front = *this ; return; }

	auto State = []( Real value , double eps )
	{
		if     ( value<0 ) return PLANE_STATE_BACK;
		else if( value>0 ) return PLANE_STATE_FRONT;
		else               return PLANE_STATE_BOUNDARY;
	};
	struct _VertexIndex
	{
		Point3D< Real > vertex;
		int state;
		_VertexIndex( void ){}
		_VertexIndex( Point3D< Real > v , int s ) : vertex(v) , state(s) {}
	};

	// Compute the indices of the refined polygons and their associated states
	_VertexIndex vertexIndices[ 2*_MAX_VERTICES ];
	unsigned int vertexIndexCount = 0;
	for( int i0=0 ; i0<size() ; i0++ )
	{
		int i1 = (i0+1)%size();
		vertexIndices[ vertexIndexCount++ ] = _VertexIndex( _vertices[i0] , State( values[i0] , eps ) );

		// If the edge crosses the plane, we need to create a new vertex
		if( values[i0]*values[i1]<0 )
		{
			Real t = values[i0] / ( values[i0] - values[i1] );
			Point3D< Real > v = _vertices[i1]*t + _vertices[i0]*(Real)(1.-t);
			v /= (Real)Length( v );
			vertexIndices[ vertexIndexCount++ ] = _VertexIndex( v , PLANE_STATE_BOUNDARY );
		}
	}
	for( unsigned int i=0 ; i<vertexIndexCount ; i++ )
	{
		if( vertexIndices[i].state!=PLANE_STATE_BACK ) front.push_back( vertexIndices[i].vertex );
		if( vertexIndices[i].state!=PLANE_STATE_FRONT ) back.push_back( vertexIndices[i].vertex );
	}

	// Set the integrated value
	back.data = front.data = data;
	if( front.size() && back.size() )
	{
		Real fArea = front.area() , bArea = back.area();
		if( fArea<=0 && bArea<=0 ) front.data /= 2 , back.data /= 2;
		else
		{
			front.data *= fArea / ( fArea + bArea );
			back.data *= bArea / ( fArea + bArea );
		}
	}
}

/////////////////////////////////////
// SphericalGeometry::Tessellation //
/////////////////////////////////////
template< class Real > template< class VertexData >
void SphericalGeometry::Tessellation< Real >::_splitPolygonPhi( Polygon< Real > p , int theta , std::unordered_map< unsigned long long , int >& vMap , std::vector< typename SphericalGeometry::Polygon< Real >::PlanarSources > &pSources , std::vector< VertexData >* vData )
{
	auto Phi = [&] ( Point3D< Real > p )
	{
		if     ( p[1]<=-1 ) return M_PI;
		else if( p[1]>= 1 ) return 0.;
		else                return acos( p[1] );
	};
	Real minPhi = Phi( _vertices[ p[0] ] ) , maxPhi = Phi( _vertices[ p[0] ] );
	for( int j=1 ; j<p.size() ; j++ ) minPhi = std::min< Real >( minPhi , Phi( _vertices[ p[j] ] ) ) , maxPhi = std::max< Real >( maxPhi , Phi( _vertices[ p[j] ] ) );
	minPhi *= (Real)_resolution / M_PI , maxPhi *= (Real)_resolution / M_PI;
	int start = (int)floor( minPhi ) , end = (int)ceil( maxPhi );

	for( int i=start+1 ; i<=end-1 ; i++ )
	{
		Real phi = ( (Real)i ) * M_PI / _resolution;
		Point3D< Real > normal( 0 , -1 , 0 );
		Real offset = -cos( phi );
		Polygon< Real > back , front;

		if( vData ) p.split( normal , offset , -1 , back , front , _vertices , pSources , *vData , vMap );
		else        p.split( normal , offset , -1 , back , front , _vertices , pSources ,          vMap );
		back.theta = theta<0 ? theta+_resolution : theta;
		back.phi = i-1;
		_polygons.push_back( back );
		p = front;
	}
	p.theta = theta<0 ? theta+_resolution : theta;
	p.phi = end-1;
	_polygons.push_back( p );
}

template< class Real > template< class VertexData >
void SphericalGeometry::Tessellation< Real >::_splitPolygon( SphericalGeometry::Polygon< Real > p , std::unordered_map< unsigned long long , int >& vMap , std::vector< typename SphericalGeometry::Polygon< Real >::PlanarSources >& pSources , std::vector< VertexData >* vData )
{
	Polygon< Real > _p[2];
	if( vData ) p.split( Point3D< Real >( 0 , 0 , 1 ) , 0. , 0 , _p[0] , _p[1] , _vertices , pSources , *vData , vMap );
	else        p.split( Point3D< Real >( 0 , 0 , 1 ) , 0. , 0 , _p[0] , _p[1] , _vertices , pSources ,          vMap );

	for( int w=0 ; w<2 ; w++ ) if( _p[w].size() )
	{
		// Recall that the spherical parameterization is given by:
		//		Phi( theta , phi ) = ( cos(theta) * sin(phi) , cos(phi) , sin(theta) * sin(phi) )
		// with theta in [0,2*pi] and phi in [0,pi]
		auto Theta = [&] ( int vIdx )
		{
			Point3D< Real > v = _vertices[ vIdx ];
			if( ( w==0 && v[2]>0 ) || ( w==1 && v[2]<0 ) ) v[2] = 0;
			Real theta = atan2( v[2] , v[0] );
			if( w==0 ) // [-Pi,0]
			{
				if( theta>0 ) theta -= 2 * M_PI;
			}
			else // [0,Pi]
			{
				if( theta<0 ) theta += 2 * M_PI;
			}
			return theta;
		};
		Polygon< Real > p = _p[w];

		Real minTheta = Theta( p[0] ) , maxTheta = Theta( p[0] );
		int vSkip = -1;
		if( pSources[ p[0] ].full() ) minTheta = Theta( p[1] ) , maxTheta = Theta( p[1] ) , vSkip = 0;
		else for( int j=1 ; j<p.size() ; j++ ) if( pSources[ p[j] ].full() ) vSkip = j;

		for( int j=1 ; j<p.size() ; j++ ) 
		{
			if( j==vSkip ) continue; 
			minTheta = std::min< Real >( minTheta , Theta( p[j] ) );
			maxTheta = std::max< Real >( maxTheta , Theta( p[j] ) );	
		}

		minTheta *= (Real)_resolution / ( 2.*M_PI ) , maxTheta *= (Real)_resolution / ( 2.*M_PI );
		int start = (int)floor( minTheta ) , end = (int)ceil( maxTheta );
		for( int i=start+1 ; i<=end-1 ; i++ )
		{
			Real theta = ( (Real)i ) * 2. * M_PI / _resolution;
			Polygon< Real > back , front;
			if( vData ) p.split( Point3D< Real >( -sin(theta) , 0 , cos(theta) ) , 0. , ( i + _resolution ) % _resolution , back , front, _vertices , pSources , *vData , vMap );
			else        p.split( Point3D< Real >( -sin(theta) , 0 , cos(theta) ) , 0. , ( i + _resolution ) % _resolution , back , front, _vertices , pSources ,          vMap );
			_splitPolygonPhi( back , i-1 , vMap , pSources , vData );
			p = front;
		}
		_splitPolygonPhi( p , end-1 , vMap , pSources , vData );
	}
}

template< class Real > template< class VertexData >
void SphericalGeometry::Tessellation< Real >::_init( const SphericalGeometry::Mesh< Real >& mesh , std::vector< VertexData >* vData , int resolution )
{
	_resolution = resolution;
	_vertices = mesh.vertices;
	_polygons.resize( mesh.polygons.size() );
#pragma omp parallel for
	for( int i=0 ; i<mesh.polygons.size() ; i++ ) _polygons[i] = Polygon< Real >( i , mesh.masses[i] , mesh.polygons[i] );
#pragma omp parallel for
	for( int i=0 ; i<_vertices.size() ; i++ ) _vertices[i] /= (Real)Length( _vertices[i] );

	std::vector< typename SphericalGeometry::Polygon< Real >::PlanarSources > pSources( _vertices.size() );
	std::unordered_map< unsigned long long , int > vMap;
	for( int i=(int)_polygons.size()-1 ; i>=0 ; i-- )
	{
		Polygon< Real > p = _polygons[i];
		_polygons[i] = _polygons.back() , _polygons.pop_back();
		_splitPolygon( p , vMap , pSources , vData );
	}
}

template< class Real >
void SphericalGeometry::Tessellation< Real >::createSGrid( SphericalGrid< Real >& sGrid , Real smooth ) const
{
	if( sGrid.resolution()!=_resolution ) 
	{
		fprintf( stderr , "[WARNING] Input sGrid resolution (%d) did not match resolution of tessellation (%d). Resizing.\n" , sGrid.resolution() , _resolution );
		sGrid.resize( _resolution ); 
	}

	auto Theta = [&]( Real i ){ return Real( 2.0*M_PI*i/sGrid.resolution() ); };
	auto Phi   = [&]( Real i ){ return Real( M_PI*(2.0*i+1)/(2.0*sGrid.resolution()) ); };

	for( int x=0 ; x<sGrid.resolution() ; x++ ) for( int y=0 ; y<sGrid.resolution() ; y++ ) sGrid(x,y) = 0;
	for( int i=0 ; i<_polygons.size() ; i++ ) sGrid( _polygons[i].theta , _polygons[i].phi ) += _polygons[i].mass;
	for( int x=0 ; x<sGrid.resolution() ; x++ )  for( int y=0 ; y<sGrid.resolution() ; y++ ) 
	{
		// Normalize the bucket values by area 
		Real theta[] = { Theta(x-0.5) , Theta(x+0.5) };
		Real phi[] = { Phi(y-0.5) , Phi(y+0.5) };
		Real area = ( theta[1] - theta[0] ) * ( -cos(phi[1]) + cos(phi[0]) );
		sGrid(x,y) = sqrt( sGrid(x,y) / area );
	}

	if( smooth>0 )
	{
		HarmonicTransform< Real > hForm; 
		FourierKeyS2< Real > key;

		hForm.ForwardFourier( sGrid , key );
		for( int i=0 ; i<key.bandWidth() ; i++ ) for( int j=0 ; j<=i ; j++ ) key(i,j) *= exp( -smooth * i * ( (Real)(i+1) ) * 1. );
		hForm.InverseFourier( key , sGrid );
	}
}

template< class Real > template< class VertexData >
void SphericalGeometry::Tessellation< Real >::_collapseCells( std::vector< VertexData >* vData )
{
	std::vector< std::vector< Polygon< Real > > > polygons( _resolution*_resolution );
	for( int i=(int)_polygons.size()-1 ; i>=0 ; i-- ) 
	{
		if( _polygons[i].size() ) polygons[ _polygons[i].theta*_resolution + _polygons[i].phi ].push_back( _polygons[i] );
		_polygons.pop_back();
	}
	auto EdgeKey = []( int v1 , int v2 )
	{
		unsigned long long _v1 = v1 , _v2 = v2;
		return ( _v1<<32 ) | _v2;
	};
	auto FactorEdgeKey = []( unsigned long long key , int& v1 , int& v2 )
	{
		v1 = (int)((key>>32)    );
		v2 = (int)((key<<32)>>32);
	};
#pragma omp parallel for
	for( int i=0 ; i<polygons.size() ; i++ )
	{
		Real mass = 0;
		std::unordered_set< unsigned long long > eSet;
		// Compute the total mass
		for( int j=0 ; j<polygons[i].size() ; j++ ) mass += polygons[i][j].mass;
		// Get the set of half-edges for all polygons within the face
		for( int j=0 ; j<polygons[i].size() ; j++ ) for( int k=0 ; k<polygons[i][j].size() ; k++ )
		{
			unsigned long long key = EdgeKey( polygons[i][j][k] , polygons[i][j][k+1] );
			if( eSet.find( key )!=eSet.end() ) fprintf( stderr , "[WARNING] Edge is already occupied\n" );
			eSet.insert(key);
		}

		// Extract the subset of boundary edges
		std::vector< std::pair< int , int > > bEdges;
		for( std::unordered_set< unsigned long long >::iterator i=eSet.begin() ; i!=eSet.end() ; i++ )
		{
			int v1 , v2;
			FactorEdgeKey( *i , v1 , v2 );
			if( eSet.find( EdgeKey( v2 , v1 ) )==eSet.end() ) bEdges.push_back( std::pair< int , int >( v1 , v2 ) );
		}
		int pCount = 0;

		// Stitch together the boundary edges (assuming each vertex has valence of exactly two)
		while( bEdges.size() )
		{
			std::pair< int , int > e = bEdges.back() ; bEdges.pop_back();
			int first = e.first , last = e.second;
			Polygon< Real > p( -1 , mass );
			p.theta = polygons[i][0].theta;
			p.phi = polygons[i][0].phi;
			p.push_back( e.first );
			while( first!=last )
			{
				bool found = false;
				for( int b=0 ; b<bEdges.size() ; b++ ) if( bEdges[b].first==last )
				{
					p.push_back( last );
					last = bEdges[b].second;
					bEdges[b] = bEdges.back() ; bEdges.pop_back();
					found = true;
					break;
				}
				if( !found ) fprintf( stderr , "[ERROR] Couldn't find next edge\n" ) , exit( 0 );
			}
#pragma omp critical
			{
				_polygons.push_back( p );
			}
			pCount++;
		}
	}

	std::vector< int > vMap( _vertices.size() , -1 );
	int idx = 0;
	for( int i=0 ; i<_polygons.size() ; i++ ) for( int j=0 ; j<_polygons[i].size() ; j++ ) if( vMap[ _polygons[i][j] ]==-1 ) vMap[ _polygons[i][j] ] = idx++;

	// Because the mesh could have holes, we can only remove valence two vertices if they are not on the boundary
	std::unordered_set< unsigned long long > eSet;
	for( int i=0 ; i<_polygons.size() ; i++ ) for( int j=0 ; j<_polygons[i].size() ; j++ )
	{
		unsigned long long key = EdgeKey( _polygons[i][j] , _polygons[i][j+1] );
		if( eSet.find( key )!=eSet.end() ) fprintf( stderr , "[WARNING] Edge is already occupied\n" );
		eSet.insert(key);
	}
	std::vector< bool > boundaryVertex( _vertices.size() , false );
	for( std::unordered_set< unsigned long long >::iterator i=eSet.begin() ; i!=eSet.end() ; i++ )
	{
		int v1 , v2;
		FactorEdgeKey( *i , v1 , v2 );
		if( eSet.find( EdgeKey( v2 , v1 ) )==eSet.end() ) boundaryVertex[v1] = boundaryVertex[v2] = true;
	}

	// Remove valence-2 vertices
	std::vector< int > vCount( _vertices.size() , 0 );
	for( int i=0 ; i<_polygons.size() ; i++ ) for( int j=0 ; j<_polygons[i].size() ; j++ ) vCount[ _polygons[i][j] ]++;
	for( int i=0 ; i<_polygons.size() ; i++ )
	{
		Polygon< Real > p( -1 , _polygons[i].mass );
		p.theta = _polygons[i].theta , p.phi = _polygons[i].phi;
		for( int j=0 ; j<_polygons[i].size() ; j++ ) if( vCount[ _polygons[i][j] ]>2 || boundaryVertex[ _polygons[i][j] ] ) p.push_back( _polygons[i][j] );
		_polygons[i] = p;
	}

	// Remove unused vertices
	std::vector< Point3D< Real > > vertices( idx );
	for( int i=0 ; i<vMap.size() ; i++ ) if( vMap[i]>=0 ) vertices[ vMap[i] ] = _vertices[i];
	for( int i=0 ; i<_polygons.size() ; i++ ) for( int j=0 ; j<_polygons[i].size() ; j++ ) _polygons[i][j] = vMap[ _polygons[i][j] ];
	_vertices = vertices;
	if( vData )
	{
		std::vector< VertexData > _vData( idx );
		for( int i=0 ; i<vMap.size() ; i++ ) if( vMap[i]>=0 ) _vData[ vMap[i] ] = (*vData)[i];
		*vData = _vData;
	}
}
template< class Real >
void SphericalGeometry::Tessellation< Real >::write( const char* fileName , bool binary ) const
{
	std::vector< PlyVertex< float > > v( _vertices.size() );
	std::vector< std::vector< int > > p( _polygons.size() );
	for( int i=0 ; i<_vertices.size() ; i++ ) v[i].point = Point3D< float >( _vertices[i] );
	for( int i=0 ; i<_polygons.size() ; i++ )
	{
		p[i].resize( _polygons[i].size() );
		for( int j=0 ; j<_polygons[i].size() ; j++ ) p[i][j] = _polygons[i][j];
	}
	std::vector< float > m;
	m.resize( _polygons.size() );
	for( int i=0 ; i<_polygons.size() ; i++ ) m[i] = (float)_polygons[i].mass;
	PlyProperty massProperties[] = { { "mass" , PLY_FLOAT , PLY_FLOAT , 0 , 0 , 0 , 0 , 0 } };
	PlyWritePolygons( fileName , v , p , m , PlyVertex< float >::WriteProperties , PlyVertex< float >::WriteComponents , massProperties , 1 , PLY_BINARY_NATIVE , NULL , 0 );
}
template< class Real >
void SphericalGeometry::Tessellation< Real >::write( const char* fileName , const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& colors , bool binary ) const
{
	std::vector< PlyParametrizedColorVertex< float , Real > > v( _vertices.size() );
	std::vector< std::vector< int > > p( _polygons.size() );
	for( int i=0 ; i<_vertices.size() ; i++ )
	{
		v[i].param = Point3D< float >( _vertices[i] );
		v[i].point = Point3D< float >(  vertices[i] );
		for( int j=0 ; j<3 ; j++ ) v[i].color[j] = (float)colors[i][j];
	}
	for( int i=0 ; i<_polygons.size() ; i++ )
	{
		p[i].resize( _polygons[i].size() );
		for( int j=0 ; j<_polygons[i].size() ; j++ ) p[i][j] = _polygons[i][j];
	}
	PlyWritePolygons( fileName , v , p , PlyParametrizedColorVertex< float , Real >::WriteProperties , PlyParametrizedColorVertex< float , Real >::WriteComponents , PLY_BINARY_NATIVE , NULL , 0 );
}

template< class Real > template< typename PointScaleFunction >
void SphericalGeometry::Tessellation< Real >::writeGrid( const char* fileName , Real smooth , PointScaleFunction psFunction , bool binary ) const
{
	SphericalGrid< Real > sGrid( _resolution );
	createSGrid( sGrid , smooth );

	Real max = 0;
	for( int i=0 ; i<_resolution ; i++ ) for( int j=0 ; j<_resolution ; j++ ) max = std::max< Real >( max , (Real)fabs( sGrid(i,j) ) );

	WriteGrid( fileName , sGrid , psFunction , binary );
}

template< class Real >
size_t SphericalGeometry::Tessellation< Real >::GridVertexNum( unsigned int resolution ){ return resolution*(resolution-1)+2; }
template< class Real >
size_t SphericalGeometry::Tessellation< Real >::GridTriangleNum( unsigned int resolution ){ return (resolution-1) * resolution * 2; }
template< class Real >
size_t SphericalGeometry::Tessellation< Real >::GridFaceNum( unsigned int resolution ){ return resolution * resolution; }
template< class Real >
Point3D< Real > SphericalGeometry::Tessellation< Real >::GridVertex( unsigned int i , unsigned int resolution )
{
	int ii , jj;
	if( i==0 ) ii = jj = 0;
	else if( i==resolution*(resolution-1)+1 ) ii=0 , jj=resolution;
	else ii = (i-1) % resolution , jj = (i-1)/resolution+1;

	Real theta = (2.*M_PI*ii)/resolution;
	Real phi = (PI*jj)/resolution;
	Point3D< Real > p;
	SphericalGrid< Real >::SetCoordinates( theta , phi , &p[0] );
	return p;
}
template< class Real >
unsigned int SphericalGeometry::Tessellation< Real >::GridVertexIndex( unsigned int i , unsigned int j , unsigned int resolution )
{
	i %= resolution;
	if( j==0 ) return 0u;
	else if( j==resolution ) return resolution*(resolution-1)+1;
	else return 1 + (j-1)*resolution+i;
}
template< class Real >
std::vector< unsigned int > SphericalGeometry::Tessellation< Real >::GridFace( unsigned int i , unsigned int j , unsigned int resolution )
{
	if( j==0 )
	{
		std::vector< unsigned int > p(3);
		p[0] = GridVertexIndex( 0   , j  , resolution );
		p[1] = GridVertexIndex( i+1 , j+1 , resolution );
		p[2] = GridVertexIndex( i   , j+1 , resolution );
		return p;
	}
	else if( j==resolution-1 )
	{
		std::vector< unsigned int > p(3);
		p[0] = GridVertexIndex( 0   , j+1 , resolution );
		p[1] = GridVertexIndex( i   , j   , resolution );
		p[2] = GridVertexIndex( i+1 , j   , resolution );
		return p;
	}
	else
	{
		std::vector< unsigned int > p(4);
		p[0] = GridVertexIndex( i   , j   , resolution );
		p[1] = GridVertexIndex( i+1 , j   , resolution );
		p[2] = GridVertexIndex( i+1 , j+1 , resolution );
		p[3] = GridVertexIndex( i   , j+1 , resolution );
		return p;
	}
}
template< class Real >
TriangleIndex SphericalGeometry::Tessellation< Real >::GridTriangle( unsigned int i , unsigned int resolution )
{
	TriangleIndex tri;
	if( i<resolution )
	{
		tri[0] = GridVertexIndex( 0   , 0 , resolution );
		tri[1] = GridVertexIndex( i+1 , 1 , resolution );
		tri[2] = GridVertexIndex( i   , 1 , resolution );
	}
	else if( i>=(resolution-1) * resolution * 2 - resolution )
	{
		i -= (resolution-1) * resolution * 2 - resolution;
		tri[0] = GridVertexIndex( 0   , resolution   , resolution );
		tri[1] = GridVertexIndex( i   , resolution-1 , resolution );
		tri[2] = GridVertexIndex( i+1 , resolution-1 , resolution );
	}
	else
	{
		i -= resolution;
		int j = i/(2*resolution)+1;
		i %= (2*resolution);
		if( i&1 )
		{
			tri[0] = GridVertexIndex( i/2   , j   , resolution );
			tri[1] = GridVertexIndex( i/2+1 , j   , resolution );
			tri[2] = GridVertexIndex( i/2+1 , j+1 , resolution );
		}
		else
		{
			tri[0] = GridVertexIndex( i/2   , j   , resolution );
			tri[1] = GridVertexIndex( i/2+1 , j+1 , resolution );
			tri[2] = GridVertexIndex( i/2   , j+1 , resolution );
		}
	}
	return tri;
}
template< class Real >
Real SphericalGeometry::Tessellation< Real >::GridFaceArea( unsigned int i , unsigned int j , unsigned int resolution )
{
	std::vector< unsigned int > face = GridFace( i , j , resolution );
	std::vector< Point3D< Real > > v( face.size() );
	for( int i=0 ; i<face.size() ; i++ ) v[i] = GridVertex( face[i] , resolution );
	return SphericalGeometry::Area( v );
}

template< class Real >
Real SphericalGeometry::Tessellation< Real >::GridTriangleArea( unsigned int i , unsigned int resolution )
{
	TriangleIndex tri = GridTriangle( i , resolution );
	std::vector< Point3D< Real > > v( 3 );
	for( int i=0 ; i<3 ; i++ ) v[i] = GridVertex( tri[i] , resolution );
	return SphericalGeometry::Area( v );
}



template< class Real > template< typename PointScaleFunction >
void SphericalGeometry::Tessellation< Real >::WriteGrid( const char* fileName , const SphericalGrid< Real >& sGrid , PointScaleFunction psFunction , bool binary )
{
	int res = sGrid.resolution();
	std::vector< PlyColorVertex< float > > v( res*(res-1)+2 );
	std::vector< std::vector< int > > polygons;

	auto Index = [&]( int i , int j )
	{
		i %= res;
		if( j==0 ) return 0;
		else if( j==res ) return res*(res-1)+1;
		else return 1 + (j-1)*res+i;
	};
	auto Position = [&]( int i , int j )
	{
		Real theta = (2.*M_PI*i)/res;
		Real phi = (PI*j)/res;
		Point3D< Real > p;
		SphericalGrid< Real >::SetCoordinates( theta , phi , &p[0] );
		return p;
	};
	auto Color = []( Real value )
	{
		Point3D< Real > blue( 0 , 0 , 255 ) , red( 255 , 0 , 0 ) , gray( 128 , 128 , 128 );
		value = std::max< Real >( -1 , std::min< Real >( value , 1 ) );
		if( value<0 ) return - blue * value + gray * ( 1 + value );
		else          return    red * value + gray * ( 1 - value );
	};

	for( int i=0 ; i<res ; i++ ) for( int j=0 ; j<=res ; j++ ) v[ Index(i,j) ].point = Position( i , j );
	// North triangles
	{
		int j=0;
		for( int i=0 ; i<res ; i++ )
		{
			std::vector< int > triangle(3);
			triangle[0] = Index(0,j);
			triangle[1] = Index(i+1,j+1);
			triangle[2] = Index(i,j+1);
			polygons.push_back( triangle );
		}
	}
	// South triangles
	{
		int j=res;
		for( int i=0 ; i<res ; i++ )
		{
			std::vector< int > triangle(3);
			triangle[0] = Index(0,j);
			triangle[1] = Index(i,j-1);
			triangle[2] = Index(i+1,j-1);
			polygons.push_back( triangle );
		}
	}
	// Interior quads
	for( int j=1 ; j<res-1 ; j++ ) for( int i=0 ; i<res ; i++ )
	{
		std::vector< int > quad(4);
		quad[0] = Index(i,j);
		quad[1] = Index(i+1,j);
		quad[2] = Index(i+1,j+1);
		quad[3] = Index(i,j+1);
		polygons.push_back( quad );
	}

	std::vector< Real > values( v.size() );
	for( int i=0 ; i<values.size() ; i++ )
	{
		Point3D< Real > p = v[i].point;
		Real x , y;
		sGrid.setCoordinates( &p[0] , x , y );
		values[i] = sGrid(x,y);
	}

	for( int i=0 ; i<v.size() ; i++ )
	{
		v[i].point *= (float)psFunction( values[i] );
		v[i].color = Point3D< float >( Color( values[i] ) );
	}

	PlyWritePolygons( fileName , v , polygons , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , PLY_ASCII , NULL , 0 );
}

template< class Real >
void SphericalGeometry::Tessellation< Real >::SampleScaleFactors( const std::vector< Point3D< Real > > &sphereVertices , const std::vector< Point3D< Real > > &meshVertices , const std::vector< TriangleIndex > &triangles , SphericalGrid< Real > &sGrid , double eps )
{
	std::vector< Real > meshTriangleAreas( triangles.size() );
#pragma omp parallel for
	for( int t=0 ; t<triangles.size() ; t++ )
	{
		Point3D< Real > v[] = { meshVertices[ triangles[t][0] ] , meshVertices[ triangles[t][1] ] , meshVertices[ triangles[t][2] ] };
		meshTriangleAreas[t] = (Real)Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2;
	}
	_SampleFaceIntegrals( sphereVertices , triangles , []( TriangleIndex ){ return (unsigned int)3; } , meshTriangleAreas , sGrid , eps );
}

template< class Real >
void SphericalGeometry::Tessellation< Real >::SampleMetric( const std::vector< Point3D< Real > > &sphereVertices , const std::vector< Point3D< Real > > &meshVertices , const std::vector< TriangleIndex > &triangles , SphericalGrid< SquareMatrix< Real , 2 > > &sGrid , double eps )
{
	auto TriangleMetric = []( const Point3D< Real > sVertices[] , const Point3D< Real > mVertices[] )
	{
		Point3D< Real > c = sVertices[0] + sVertices[1] + sVertices[2];
		c /= (Real)Length(c);
		std::pair< Point3D< Real > , Point3D< Real > > tangents = SphericalGeometry::RegularGrid< Real >::SphereTangents( c );
		Matrix< Real , 2 , 3 > T;
		for( int i=0 ; i<3 ; i++ ) T(0,i) = tangents.first[i] , T(1,i) = tangents.second[i];

		// Solve for the linear transformation taking sVertices[i]-sVertices[0] -> mVertices[i]-mVertices[0]
		Point3D< Real > sv[] = { sVertices[1]-sVertices[0] , sVertices[2]-sVertices[0] , Point3D< Real >::CrossProduct( sVertices[1]-sVertices[0] , sVertices[2]-sVertices[0] ) };
		Point3D< Real > mv[] = { mVertices[1]-mVertices[0] , mVertices[2]-mVertices[0] , Point3D< Real >::CrossProduct( mVertices[1]-mVertices[0] , mVertices[2]-mVertices[0] ) };
		SquareMatrix< Real , 3 > sM , mM;
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) sM(i,j) = sv[i][j] , mM(i,j) = mv[i][j];
		SquareMatrix< Real , 3 > D = mM * sM.inverse();

		return SquareMatrix< Real , 2 >( T.transpose() * D.transpose() * D * T );
	};
	auto MeshMetric = [&]( const ClippedPrimalTriangle< Real , Point3D< Real > > &p )
	{
		SquareMatrix< Real , 2 > M;
		Point3D< Real > sCenter , mCenter , n;
		for( int i=0 ; i<p.size() ; i++ )
		{
			sCenter += p[i].first , mCenter += p[i].second;
			Point3D< Real > b = p[i+1].second + p[i].second;
			Point3D< Real > e = p[i+1].second - p[i].second;
			n += Point3D< Real >::CrossProduct( b , e ) / 2;
		}
		sCenter /= (Real)Length(sCenter) , mCenter /= (Real)p.size();
		Point3D< Real > sVertices[3] , mVertices[3];
		sVertices[0] = sCenter , mVertices[0] = mCenter;
		for( int i=0 ; i<p.size() ; i++ )
		{
			sVertices[1] = p[i  ].first , mVertices[1] = p[i  ].second;
			sVertices[2] = p[i+1].first , mVertices[2] = p[i+1].second;
			M += TriangleMetric( sVertices , mVertices );
		}
		return M * (Real)Length( n ) / 2;
	};

	unsigned int threads = omp_get_max_threads();
	SphericalGrid< SquareMatrix< Real , 2 > > *sGrids = new SphericalGrid< SquareMatrix< Real , 2 > >[ threads ];
#pragma omp parallel for
	for( int t=0 ; t<(int)threads ; t++ )
	{
		sGrids[t].resize( sGrid.resolution() );
		for( int i=0 ; i<sGrid.resolution() ; i++ ) for( int j=0 ; j<sGrid.resolution() ; j++ ) sGrids[t](i,j) = SquareMatrix< Real , 2 >();
	}
#pragma omp parallel for
	for( int i=0 ; i<sGrid.resolution() ; i++ ) for( int j=0 ; j<sGrid.resolution() ; j++ ) sGrid(i,j) = SquareMatrix< Real , 2 >();
#pragma omp parallel for
	for( int t=0 ; t<triangles.size() ; t++ )
	{
		ClippedPrimalTriangle< Real , Point3D< Real > > p;
		for( int i=0 ; i<3 ; i++ ) p.push_back( std::pair< Point3D< Real > , Point3D< Real > >( sphereVertices[ triangles[t][i] ] , meshVertices[ triangles[t][i] ] ) );
		_Accumulate( p , sGrids[ omp_get_thread_num() ] , MeshMetric , eps );
	}
#pragma omp parallel for
	for( int i=0 ; i<sGrid.resolution() ; i++ ) for( int j=0 ; j<sGrid.resolution() ; j++ ) for( unsigned int t=0 ; t<threads ; t++ ) sGrid(i,j) += sGrids[t](i,j);
	delete[] sGrids;

#pragma omp parallel for
	for( int j=0 ; j<sGrid.resolution() ; j++ )
	{
		Real scaleFactor = (Real)( 1. / GridFaceArea( 0 , j , sGrid.resolution() ) );
		for( int i=0 ; i<sGrid.resolution() ; i++ ) sGrid(i,j) *= scaleFactor;
	}
}

template< class Real >
template< typename Data , typename Face , typename FaceSizeFunctor >
void SphericalGeometry::Tessellation< Real >::_SampleVertexValues( const std::vector< Point3D< Real > > &sphereVertices , const std::vector< Face > &faces , const FaceSizeFunctor &FSF , const std::vector< Data > &vertexValues , SphericalGrid< Data > &sGrid , double eps )
{
	typedef VectorSpaceDirectSumElement< Real , Data , Real > _Data;
	std::vector< _Data > _vertexValues( vertexValues.size() );
#pragma omp parallel for
	for( int i=0 ; i<vertexValues.size() ; i++ )
	{
		std::get< 0 >( _vertexValues[i] ) = vertexValues[i];
		std::get< 1 >( _vertexValues[i] ) = 1;
	}

	auto AccumulationFunction = []( const ClippedPrimalTriangle< Real , _Data > &p ){ return p.integral(); };

	SphericalGrid< _Data > _sGrid;
	_sGrid.resize( sGrid.resolution() );

	unsigned int threads = omp_get_max_threads();
	SphericalGrid< _Data > *_sGrids = new SphericalGrid< _Data >[ threads ];
#pragma omp parallel for
	for( int t=0 ; t<(int)threads ; t++ )
	{
		_sGrids[t].resize( _sGrid.resolution() );
		for( int i=0 ; i<_sGrid.resolution() ; i++ ) for( int j=0 ; j<_sGrid.resolution() ; j++ ) std::get< 0 >( _sGrids[t](i,j) ) = Data() , std::get< 1 >( _sGrids[t](i,j) ) = 0;
	}

#pragma omp parallel for
	for( int i=0 ; i<_sGrid.resolution() ; i++ ) for( int j=0 ; j<_sGrid.resolution() ; j++ ) std::get< 0 >( _sGrid(i,j) ) = Data() , std::get< 1 >( _sGrid(i,j) ) = 0;

#pragma omp parallel for
	for( int f=0 ; f<faces.size() ; f++ )
	{
		if( FSF( faces[f] )==3 )
		{
			ClippedPrimalTriangle< Real , _Data > p;
			for( int i=0 ; i<3 ; i++ ) p.push_back( sphereVertices[ faces[f][i] ] , _vertexValues[ faces[f][i] ] );
			_Accumulate( p , _sGrids[ omp_get_thread_num() ] , AccumulationFunction , eps );
		}
		else
		{
			MinimalAreaTriangulation< Real > MAT;
			std::vector< Point3D< Real > > _vertices( FSF( faces[f] ) );
			std::vector< TriangleIndex > _triangles;
			for( unsigned int j=0 ; j<FSF( faces[f] ) ; j++ ) _vertices[j] = sphereVertices[ faces[f][j] ];
			MAT.GetTriangulation( _vertices , _triangles );
			for( int i=0 ; i<_triangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) _triangles[i][j] = faces[f][ _triangles[i][j] ];
			for( int t=0 ; t<_triangles.size() ; t++ )
			{
				ClippedPrimalTriangle< Real , _Data > p;
				for( int i=0 ; i<3 ; i++ ) p.push_back( sphereVertices[ _triangles[t][i] ] , _vertexValues[ _triangles[t][i] ] );
				_Accumulate( p , _sGrids[ omp_get_thread_num() ] , AccumulationFunction , eps );
			}
		}
	}

#pragma omp parallel for
	for( int i=0 ; i<_sGrid.resolution() ; i++ ) for( int j=0 ; j<_sGrid.resolution() ; j++ ) for( unsigned int t=0 ; t<threads ; t++ ) _sGrid(i,j) += _sGrids[t](i,j);
	delete[] _sGrids;

#pragma omp parallel for
	for( int i=0 ; i<sGrid.resolution() ; i++ ) for( int j=0 ; j<sGrid.resolution() ; j++ ) sGrid(i,j) = std::get< 0 >( _sGrid(i,j) ) / std::get< 1 >( _sGrid(i,j) );
}

template< class Real >
template< typename Data , typename Face , typename FaceSizeFunctor >
void SphericalGeometry::Tessellation< Real >::_SampleFaceIntegrals( const std::vector< Point3D< Real > > &sphereVertices , const std::vector< Face > &faces , const FaceSizeFunctor &FSF , const std::vector< Data > &faceIntegrals , SphericalGrid< Data > &sGrid , double eps )
{
	auto AccumulationFunction = []( const ClippedDualTriangle< Real , Data > &p ){ return p.integral(); };

	unsigned int threads = omp_get_max_threads();
	SphericalGrid< Data > *_sGrids = new SphericalGrid< Data >[ threads ];
#pragma omp parallel for
	for( int t=0 ; t<(int)threads ; t++ )
	{
		_sGrids[t].resize( sGrid.resolution() );
		for( int i=0 ; i<sGrid.resolution() ; i++ ) for( int j=0 ; j<sGrid.resolution() ; j++ ) _sGrids[t](i,j) = Data();
	}

#pragma omp parallel for
	for( int i=0 ; i<sGrid.resolution() ; i++ ) for( int j=0 ; j<sGrid.resolution() ; j++ ) sGrid(i,j) = Data();

#pragma omp parallel for
	for( int f=0 ; f<faces.size() ; f++ )
	{
		if( FSF( faces[f] )==3 )
		{
			ClippedDualTriangle< Real , Data > p;
			for( int i=0 ; i<3 ; i++ ) p.push_back( sphereVertices[ faces[f][i] ] );
			p.data = faceIntegrals[f];
			_Accumulate( p , _sGrids[ omp_get_thread_num() ] , AccumulationFunction , eps );
		}
		else
		{
			MinimalAreaTriangulation< Real > MAT;
			std::vector< Point3D< Real > > _vertices( FSF( faces[f] ) );
			std::vector< TriangleIndex > _triangles;
			for( unsigned int j=0 ; j<FSF( faces[f] ) ; j++ ) _vertices[j] = sphereVertices[ faces[f][j] ];
			MAT.GetTriangulation( _vertices , _triangles );
			std::vector< ClippedDualTriangle< Real , Data > > p( _triangles.size() );
			Real area = 0;
			for( int t=0 ; t<_triangles.size() ; t++ )
			{
				for( int i=0 ; i<3 ; i++ ) p[t].push_back( sphereVertices[ faces[f][ _triangles[t][i] ] ] );
				Real a = p[t].area();
				p[t].data = faceIntegrals[f] * a;
				area += a;
			}
			for( int t=0 ; t<_triangles.size() ; t++ )
			{
				p[t].data /= area;
				_Accumulate( p[t] , _sGrids[ omp_get_thread_num() ] , AccumulationFunction , eps );
			}
		}
	}
#pragma omp parallel for
	for( int i=0 ; i<sGrid.resolution() ; i++ ) for( int j=0 ; j<sGrid.resolution() ; j++ ) for( unsigned int t=0 ; t<threads ; t++ ) sGrid(i,j) += _sGrids[t](i,j);
	delete[] _sGrids;

#pragma omp parallel for
	for( int j=0 ; j<sGrid.resolution() ; j++ )
	{
		Real scaleFactor = (Real)( 1. / GridFaceArea( 0 , j , sGrid.resolution() ) );
		for( int i=0 ; i<sGrid.resolution() ; i++ ) sGrid(i,j) *= scaleFactor;
	}
}

template< class Real >
template< typename Polygon , typename Data , typename PolygonDataFunctor >
void SphericalGeometry::Tessellation< Real >::_Accumulate( Polygon p , SphericalGrid< Data > &sGrid , const PolygonDataFunctor &polygonDataFunctor , double eps )
{
	int resolution = sGrid.resolution();
	Real a = p.area();

	Polygon _p[2];
	p.split( Point3D< Real >( 0 , 0 , 1 ) , 0. ,  _p[0] , _p[1] , eps );

	for( int w=0 ; w<2 ; w++ ) if( _p[w].size() )
	{
		// Recall that the spherical parameterization is given by:
		//		Phi( theta , phi ) = ( cos(theta) * sin(phi) , cos(phi) , sin(theta) * sin(phi) )
		// with theta in [0,2*pi] and phi in [0,pi]
		auto Theta = [&] ( Point3D< Real > v )
		{
			if( ( w==0 && v[2]>0 ) || ( w==1 && v[2]<0 ) ) v[2] = 0;
			Real theta = atan2( v[2] , v[0] );
			if( w==0 ) // [-Pi,0]
			{
				if( theta>0 ) theta -= 2 * M_PI;
			}
			else // [0,Pi]
			{
				if( theta<0 ) theta += 2 * M_PI;
			}
			return theta;
		};
		Polygon p = _p[w];

		Real minTheta = Theta( p[0] ) , maxTheta = Theta( p[0] );
		for( int j=1 ; j<p.size() ; j++ ) 
		{
			minTheta = std::min< Real >( minTheta , Theta( p[j] ) );
			maxTheta = std::max< Real >( maxTheta , Theta( p[j] ) );	
		}

		minTheta *= (Real)resolution / ( 2.*M_PI ) , maxTheta *= (Real)resolution / ( 2.*M_PI );
		int start = (int)floor( minTheta ) , end = (int)ceil( maxTheta );
		for( int i=start+1 ; i<=end-1 ; i++ )
		{
			Real theta = ( (Real)i ) * 2. * M_PI / resolution;
			Polygon back , front;
			p.split( Point3D< Real >( -sin(theta) , 0 , cos(theta) ) , 0. , back , front , eps );
			_AccumulatePhi( i-1 , back , sGrid , polygonDataFunctor , eps );
			p = front;
		}
		_AccumulatePhi( end-1 , p , sGrid , polygonDataFunctor , eps );
	}
}

template< class Real >
template< typename Polygon , typename Data , typename PolygonDataFunctor >
void SphericalGeometry::Tessellation< Real >::_AccumulatePhi( int theta , Polygon p , SphericalGrid< Data > &sGrid , const PolygonDataFunctor &polygonDataFunctor , double eps )
{
	int resolution = sGrid.resolution();
	theta = theta<0 ? theta+resolution : theta;
	auto Phi = [&] ( Point3D< Real > p )
	{
		if     ( p[1]<=-1 ) return (double)M_PI;
		else if( p[1]>= 1 ) return (double)0.;
		else                return (double)acos( p[1] );
	};
	Real minPhi = Phi( p[0] ) , maxPhi = Phi( p[0] );
	for( int j=1 ; j<p.size() ; j++ ) minPhi = std::min< Real >( minPhi , Phi( p[j] ) ) , maxPhi = std::max< Real >( maxPhi , Phi( p[j] ) );
	minPhi *= (Real)resolution / M_PI , maxPhi *= (Real)resolution / M_PI;
	int start = (int)floor( minPhi ) , end = (int)ceil( maxPhi );

	for( int i=start+1 ; i<=end-1 ; i++ )
	{
		Real phi = ( (Real)i ) * M_PI / resolution;
		Point3D< Real > normal( 0 , -1 , 0 );
		Real offset = -cos( phi );
		Polygon back , front;

		p.split( normal , offset , back , front , eps );
		AddAtomic( sGrid( theta , i-1 ) , polygonDataFunctor( back ) );
		p = front;
	}
	AddAtomic( sGrid( theta , end-1 ) , polygonDataFunctor( p ) );
}


////////////////////////////////////////////
// SphericalGeometry::CircumscribedSphere //
////////////////////////////////////////////
template< class Real >
void SphericalGeometry::CircumscribedSphere< Real >::_RandomCenter( const std::vector< Point3D< Real > >& vertices , Point3D< Real >& c , Real& a )
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, (int)vertices.size() - 1);

	Point3D< Real > v[4]; 
	for( int i=0 ; i<4 ; i++ ) v[i] = vertices[dis(gen)];

	// Compute determinant matrices 
	SquareMatrix< Real , 4 > A , Dx , Dy , Dz; 
	Real dx , dy , dz;

	for( int i=0 ; i<4 ; i++ ) 
	{
		for( int j=0 ; j<3 ; j++ ) A(i, j) = v[i][j]; 
		A( i , 3 ) = 1; 
	}
	a = A.determinant(); 

	for( int i=0 ; i<4 ; i++ )
	{
		Dx( i , 0 ) = v[i].squareNorm(); 
		Dx( i , 1 ) = v[i][1];
		Dx( i , 2 ) = v[i][2];
		Dx( i , 3 ) = 1; 
	}
	dx = Dx.determinant(); 

	for( int i=0 ; i<4 ; i++ )
	{
		Dy( i , 0 ) = v[i].squareNorm(); 
		Dy( i , 1 ) = v[i][0];
		Dy( i , 2 ) = v[i][2];
		Dy( i , 3 ) = 1; 
	}
	dy = -Dy.determinant(); 

	for( int i=0 ; i<4 ; i++ )
	{
		Dz( i , 0 ) = v[i].squareNorm(); 
		Dz( i , 1 ) = v[i][0];
		Dz( i , 2 ) = v[i][1];
		Dz( i , 3 ) = 1; 
	}
	dz = Dz.determinant(); 

	c[0] = dx/(2.*a); 
	c[1] = dy/(2.*a);
	c[2] = dz/(2.*a);
}

template< class Real >
Point3D< Real > SphericalGeometry::CircumscribedSphere< Real >::Center( const std::vector< Point3D< Real > >& vertices , int iters )
{
	Real w=0 , x=0 , y=0 , z=0;
#pragma omp parallel for reduction( + : w , x , y , z )
	for( int i=0 ; i<iters ; i++ )
	{
		Point3D< Real > c;
		Real _w;
		_RandomCenter( vertices , c , _w );
		if( _w>0 ) x += c[0] *_w , y += c[1] * _w , z += c[2] * _w , w +=_w;
	}
	return Point3D< Real >( x/w , y/w , z/w );
}
