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

///////////////////////////
// MonomialIntegralTable //
///////////////////////////
template< typename Real , unsigned int D >
void SphericalGeometry::MonomialIntegralTable< Real , D >::SetIntegralTable( Real start , Real end , Real table[D+1] )
{
	Real _start=start , _end=end;
	for( unsigned int d=0 ; d<=D ; d++ , _start *= start , _end *= end ) table[d] = ( _end - _start ) / ( d+1 );
}
template< typename Real , unsigned int D >
void SphericalGeometry::MonomialIntegralTable< Real , D >::init( Real phiStart , Real phiEnd , unsigned int mcSamples )
{
	Real dPhi = ( phiEnd - phiStart ) / mcSamples;
	for( unsigned int d=0 ; d<=D ; d++ ) _sinTable[d] = _cscTable[d] = ____table[d] = 0;
	for( unsigned int j=0 ; j<mcSamples ; j++ )
	{
		Real phi = phiStart + dPhi * (Real)(j+0.5);
		Real a1 = (Real)( dPhi * sin(phi) );
		Real a2 = (Real)( dPhi / sin(phi) );
		Real _phi = 1;
		for( int d=0 ; d<=D ; d++ , _phi *= phi  ) _sinTable[d] += _phi * a1 , _cscTable[d] += _phi * a2;
	}
	SetIntegralTable( phiStart , phiEnd , ____table );
}

////////////////
// Polynomial //
////////////////
template< typename Real , unsigned int D1 , unsigned int D2 >
Real SphericalGeometry::Polynomial< Real , D1 , D2 >::operator()( Real x , Real y ) const
{
	Real v = 0 , _x = 1;
	for( int d1=0 ; d1<=D1 ; d1++ , _x *= x )
	{
		Real _v = 0 , _y = 1;
		for( int d2=0 ; d2<=D2 ; d2++ , _y *= y ) _v += _y * _coefficients[d1][d2];
		v += _v * _x;
	}
	return v;
}
template< typename Real , unsigned int D1 , unsigned int D2 >
SphericalGeometry::Polynomial< Real , D1-1 , D2 > SphericalGeometry::Polynomial< Real , D1 , D2 >::dx( void ) const
{
	SphericalGeometry::Polynomial< Real , D1-1 , D2 > p;
	for( int d1=1 ; d1<=D1 ; d1++ ) for( int d2=0 ; d2<=D2 ; d2++ ) p[d1-1][d2] = _coefficients[d1][d2] * d1;
	return p;
}
template< typename Real , unsigned int D1 , unsigned int D2 >
SphericalGeometry::Polynomial< Real , D1 , D2-1 > SphericalGeometry::Polynomial< Real , D1 , D2 >::dy( void ) const
{
	SphericalGeometry::Polynomial< Real , D1 , D2-1 > p;
	for( int d1=0 ; d1<=D1 ; d1++ ) for( int d2=1 ; d2<=D2 ; d2++ ) p[d1][d2-1] = _coefficients[d1][d2] * d2;
	return p;
}

template< typename Real , unsigned int D1 , unsigned int D2 >
template< unsigned int _D1 , unsigned int _D2 >
SphericalGeometry::Polynomial< Real , D1+_D1 , D2+_D2 > SphericalGeometry::Polynomial< Real , D1 , D2 >::operator * ( const Polynomial< Real , _D1 , _D2 > &p ) const
{
	Polynomial< Real , D1+_D1 , D2+_D2 > q;
	for( unsigned int d1=0 ; d1<=D1 ; d1++ ) for( unsigned int d2=0 ; d2<=D2 ; d2++ ) for( unsigned int _d1=0 ; _d1<=_D1 ; _d1++ ) for( unsigned int _d2=0 ; _d2<=_D2 ; _d2++ )
		q[d1+_d1][d2+_d2] += _coefficients[d1][d2] * p[_d1][_d2];
	return q;
}

template< typename Real , unsigned int D1 , unsigned int D2 >
SphericalGeometry::Polynomial< Real , D1 , D2 > SphericalGeometry::Polynomial< Real , D1 , D2 >::operator - ( void ) const
{
	Polynomial< Real , D1 , D2 > q;
	for( unsigned int d1=0 ; d1<=D1 ; d1++ ) for( unsigned int d2=0 ; d2<=D2 ; d2++ ) q._coefficients[d1][d2] -= _coefficients[d1][d2];
	return q;
}

template< typename Real , unsigned int D1 , unsigned int D2 >
SphericalGeometry::Polynomial< Real , D1 , D2 > SphericalGeometry::Polynomial< Real , D1 , D2 >::operator + ( const Polynomial< Real , D1 , D2 > &p ) const
{
	Polynomial< Real , D1 , D2 > q;
	for( unsigned int d1=0 ; d1<=D1 ; d1++ ) for( unsigned int d2=0 ; d2<=D2 ; d2++ ) q._coefficients[d1][d2] = _coefficients[d1][d2] + p._coefficients[d1][d2];
	return q;
}

template< typename Real , unsigned int D1 , unsigned int D2 >
SphericalGeometry::Polynomial< Real , D1 , D2 >& SphericalGeometry::Polynomial< Real , D1 , D2 >::operator += ( const Polynomial< Real , D1 , D2 > &p )
{
	for( unsigned int d1=0 ; d1<=D1 ; d1++ ) for( unsigned int d2=0 ; d2<=D2 ; d2++ ) _coefficients[d1][d2] += p._coefficients[d1][d2];
	return *this;
}

template< typename Real , unsigned int D1 , unsigned int D2 >
SphericalGeometry::Polynomial< Real , D1 , D2 > SphericalGeometry::Polynomial< Real , D1 , D2 >::operator * ( Real s ) const
{
	Polynomial< Real , D1 , D2 > q;
	for( unsigned int d1=0 ; d1<=D1 ; d1++ ) for( unsigned int d2=0 ; d2<=D2 ; d2++ ) q._coefficients[d1][d2] = _coefficients[d1][d2] * s;
	return q;
}

template< typename Real , unsigned int D1 , unsigned int D2 >
SphericalGeometry::Polynomial< Real , D1 , D2 >& SphericalGeometry::Polynomial< Real , D1 , D2 >::operator *= ( Real s )
{
	for( unsigned int d1=0 ; d1<=D1 ; d1++ ) for( unsigned int d2=0 ; d2<=D2 ; d2++ ) _coefficients[d1][d2] *= s;
	return *this;
}

template< typename Real , unsigned int D1 , unsigned int D2 >
Real SphericalGeometry::Polynomial< Real , D1 , D2 >::integrate( Real thetaStart , Real thetaEnd , Real phiStart , Real phiEnd ) const
{
	Real thetaIntegrals[D1+1] , phiIntegrals[D2+1];
	MonomialIntegralTable< Real , D1 >::SetIntegralTable( thetaStart , thetaEnd , thetaIntegrals );
	MonomialIntegralTable< Real , D2 >::SetIntegralTable( phiStart , phiEnd , phiIntegrals );
	Real value = 0;
	for( unsigned int d1=0 ; d1<=D1 ; d1++ ) for( unsigned int d2=0 ; d2<=D2 ; d2++ ) value += _coefficients[d1][d2] * thetaIntegrals[d1] * phiIntegrals[d2];
	return value;
}

template< typename Real , unsigned int D1 , unsigned int D2 >
Real SphericalGeometry::Polynomial< Real , D1 , D2 >::integrateSine( Real thetaStart , Real thetaEnd , Real phiStart , Real phiEnd , unsigned int mcSamples ) const
{
	Real value = 0;
	Real dTheta = ( thetaEnd - thetaStart ) / mcSamples , dPhi = ( phiEnd - phiStart ) / mcSamples;
	for( unsigned int j=0 ; j<mcSamples ; j++ )
	{
		Real phi = phiStart + dPhi * (Real)(j+0.5);
		Real area = (Real)( dPhi * dTheta * sin(phi) );
		for( unsigned int i=0 ; i<mcSamples ; i++ )
		{
			Real theta = thetaStart + dTheta * (Real)(i+0.5);
			value += this->operator()( theta , phi ) * area;
		}
	}
	return value;
}

template< typename Real , unsigned int D1 , unsigned int D2 >
Real SphericalGeometry::Polynomial< Real , D1 , D2 >::integrateCosecant( Real thetaStart , Real thetaEnd , Real phiStart , Real phiEnd , unsigned int mcSamples ) const
{
	Real value = 0;
	Real dTheta = ( thetaEnd - thetaStart ) / mcSamples , dPhi = ( phiEnd - phiStart ) / mcSamples;
	for( unsigned int j=0 ; j<mcSamples ; j++ )
	{
		Real phi = phiStart + dPhi * (Real)(j+0.5);
		Real area = (Real)( dPhi * dTheta / sin(phi) );
		for( unsigned int i=0 ; i<mcSamples ; i++ )
		{
			Real theta = thetaStart + dTheta * (Real)(i+0.5);
			value += this->operator()( theta , phi ) * area;
		}
	}
	return value;
}

template< typename Real , unsigned int D1 , unsigned int D2 >
template< unsigned int D >
Real SphericalGeometry::Polynomial< Real , D1 , D2 >::integrateSine( Real thetaStart , Real thetaEnd , const MonomialIntegralTable< Real , D > &mit ) const
{
	static_assert( D>=D2 , "[ERROR] Monomial integration table degree must be at least as large as the phi degree" );
	Real value = 0;
	Real thetaIntegrals[D1+1];
	MonomialIntegralTable< Real , D1 >::SetIntegralTable( thetaStart , thetaEnd , thetaIntegrals );
	for( unsigned int d1=0 ; d1<=D1 ; d1++ ) for( unsigned int d2=0 ; d2<=D2 ; d2++ ) value += _coefficients[d1][d2] * thetaIntegrals[d1] * mit.table< 1 >(d2);
	return value;
}

template< typename Real , unsigned int D1 , unsigned int D2 >
template< unsigned int D >
Real SphericalGeometry::Polynomial< Real , D1 , D2 >::integrateCosecant( Real thetaStart , Real thetaEnd , const MonomialIntegralTable< Real , D > &mit ) const
{
	static_assert( D>=D2 , "[ERROR] Monomial integration table degree must be at least as large as the phi degree" );
	Real value = 0;
	Real thetaIntegrals[D1+1];
	MonomialIntegralTable< Real , D1 >::SetIntegralTable( thetaStart , thetaEnd , thetaIntegrals );
	for( unsigned int d1=0 ; d1<=D1 ; d1++ ) for( unsigned int d2=0 ; d2<=D2 ; d2++ ) value += _coefficients[d1][d2] * thetaIntegrals[d1] * mit.table< -1 >(d2);
	return value;
}

template< typename Real , unsigned int D1 , unsigned int D2 >
std::ostream &SphericalGeometry::operator << ( std::ostream &os , const SphericalGeometry::Polynomial< Real , D1 , D2 > &p )
{
	bool first = true;
	for( int d1=0 ; d1<=D1 ; d1++ ) for( int d2=0 ; d2<=D2 ; d2++ ) if( p[d1][d2] )
	{
		if( p[d1][d2]>0 )
		{
			if( first ) os << p[d1][d2];
			else        os << " + " <<  p[d1][d2];
		}
		else
		{
			if( first ) os << "- " << -p[d1][d2];
			else        os << " - " << -p[d1][d2];
		}
		if( d1==0 ) ;
		else if( d1==1 ) os << " x";
		else             os << " x^" << d1;
		if( d2==0 ) ;
		else if( d2==1 ) os << " y";
		else             os << " y^" << d2;
		first = false;
	}
	return os;
}

template< typename Real , unsigned int D1 , unsigned int D2 , int SinYExponent >
Real SphericalGeometry::SinYPolynomial< Real , D1 , D2 , SinYExponent >::integrate( Real xStart , Real xEnd , Real yStart , Real yEnd , unsigned int mcSamples ) const
{
	Real value = 0;
	Real dx = ( xEnd - xStart ) / mcSamples , dy = ( yEnd - yStart ) / mcSamples;
	for( unsigned int j=0 ; j<mcSamples ; j++ )
	{
		Real y = yStart + dy * (Real)(j+0.5);
		for( unsigned int i=0 ; i<mcSamples ; i++ )
		{
			Real x = xStart + dx * (Real)(i+0.5);
			value += this->operator()( x , y );
		}
	}
	return value * dx * dy;
}

template< typename Real , unsigned int D1 , unsigned int D2 , int SinYExponent >
template< unsigned int D >
Real SphericalGeometry::SinYPolynomial< Real , D1 , D2 , SinYExponent >::integrate( Real xStart , Real xEnd , const MonomialIntegralTable< Real , D > &mit ) const
{
	static_assert( D>=D2 , "[ERROR] Monomial integration table degree must be at least as large as the y degree" );
	Real value = 0;
	Real xIntegrals[D1+1];
	MonomialIntegralTable< Real , D1 >::SetIntegralTable( xStart , xEnd , xIntegrals );
	for( unsigned int d1=0 ; d1<=D1 ; d1++ ) for( unsigned int d2=0 ; d2<=D2 ; d2++ ) value += _poly[d1][d2] * xIntegrals[d1] * mit.template table< SinYExponent >(d2);
	return value;
}

template< typename Real , unsigned int D1 , unsigned int D2 , int SinYExponent >
std::ostream &SphericalGeometry::operator << ( std::ostream &os , const SphericalGeometry::SinYPolynomial< Real , D1 , D2 , SinYExponent > &p )
{
	if     ( SinYExponent== 0 ) return os << p();
	else if( SinYExponent== 1 ) return os << "( " << p() << " ) * sin( y )" ;
	else if( SinYExponent==-1 ) return os << "( " << p() << " ) / sin( y )" ;
	else if( SinYExponent>  0 ) return os << "( " << p() << " ) * sin^" <<  SinYExponent << "( y )" ;
	else if( SinYExponent<  0 ) return os << "( " << p() << " ) / sin^" << -SinYExponent << "( y )" ;
}

/////////////////
// RegularGrid //
/////////////////
template< typename Real >
template< unsigned int Cols , unsigned int Rows >
SparseMatrix< Real , int > SphericalGeometry::RegularGrid< Real >::ExpandMatrix( const SparseMatrix< Matrix< Real , Cols , Rows > , int > &M )
{
	SparseMatrix< Real , int > _M;
	_M.resize( M.rows * Rows );
#pragma omp parallel for
	for( int i=0 ; i<M.rows ; i++ )
	{
		for( int r=0 ; r<Rows ; r++ ) _M.SetRowSize( i*Rows+r , M.rowSizes[i]*Cols );
		for( int j=0 ; j<M.rowSizes[i] ; j++ ) for( int r=0 ; r<Rows ; r++ ) for( int c=0 ; c<Cols ; c++ ) _M[i*Rows+r][j*Cols+c] = MatrixEntry< Real , int >( M[i][j].N*Cols+c , M[i][j].Value(c,r) );
	}
	return _M;
}
template< typename Real >
template< unsigned int Dim >
SparseMatrix< Real , int > SphericalGeometry::RegularGrid< Real >::ExpandMatrix( const SparseMatrix< SquareMatrix< Real , Dim > , int > &M )
{
	SparseMatrix< Real , int > _M;
	_M.resize( M.rows * Dim );
#pragma omp parallel for
	for( int i=0 ; i<M.rows ; i++ )
	{
		for( int r=0 ; r<Dim ; r++ ) _M.SetRowSize( i*Dim+r , M.rowSizes[i]*Dim );
		for( int j=0 ; j<M.rowSizes[i] ; j++ ) for( int r=0 ; r<Dim ; r++ ) for( int c=0 ; c<Dim ; c++ ) _M[i*Dim+r][j*Dim+c] = MatrixEntry< Real , int >( M[i][j].N*Dim+c , M[i][j].Value(c,r) );
	}
	return _M;
}

template< typename Real >
Point3D< Real > SphericalGeometry::RegularGrid< Real >::SpherePoint( Point2D< Real > thetaPhi )
{
	return Point3D< Real >( (Real)( sin(thetaPhi[1]) * cos(thetaPhi[0]) ) , (Real)cos(thetaPhi[1]) , (Real)( sin(thetaPhi[1]) * sin(thetaPhi[0]) ) );
}

template< typename Real >
Point2D< Real > SphericalGeometry::RegularGrid< Real >::SphereParameters( Point3D< Real > p )
{
	Point2D< Real > thetaPhi;
	thetaPhi[1] = (Real)acos( std::max< Real >( std::min< Real >( p[1] , 1 ) , -1 ) );
	if( p[1]>-1 && p[1]<1 ) thetaPhi[0] = (Real)atan2( p[2] , p[0] );
	else thetaPhi[0] = 0;
	return thetaPhi;
}
template< typename Real >
std::pair< Point3D< Real > , Point3D< Real > > SphericalGeometry::RegularGrid< Real >::SphereTangents( Real theta , Real phi )
{
	Point3D< Real > dTheta( -(Real)( sin(theta) ) , 0 , (Real)( cos(theta) ) );
	Point3D< Real > dPhi( (Real)( cos(phi) * cos(theta) ) , -(Real)sin(phi) , (Real)( cos(phi) * sin(theta) ) );
	return std::pair< Point3D< Real > , Point3D< Real > >( dTheta , dPhi );
}
template< typename Real >
std::pair< Point3D< Real > , Point3D< Real > > SphericalGeometry::RegularGrid< Real >::SphereTangents( Point3D< Real > p )
{
	Point2D< Real > thetaPhi = SphereParameters( p );
	return SphereTangents( thetaPhi[0] , thetaPhi[1] );
}

template< typename Real >
SquareMatrix< Real , 2 > SphericalGeometry::RegularGrid< Real >::J( void )
{
	SquareMatrix< Real , 2 > j;
	j(0,1) = 1 , j(1,0) = -1;
	return j;
}

template< typename Real >
SquareMatrix< Real , 2 > SphericalGeometry::RegularGrid< Real >::J( const SquareMatrix< Real , 2 > &g )
{
	// The rotation by 90 degrees operator, w.r.t. the metric g, is given by:
	//		J_g = J * g / sqrt( det(g) )
	return J() * g / (Real)sqrt( g.determinant() );
}

template< typename Real > Point2D< Real > SphericalGeometry::RegularGrid< Real >::J(                                     Point2D< Real > p ){ return J( )*p; }
template< typename Real > Point2D< Real > SphericalGeometry::RegularGrid< Real >::J( const SquareMatrix< Real , 2 > &g , Point2D< Real > p ){ return J(g)*p; }

template< typename Real >
size_t SphericalGeometry::RegularGrid< Real >::VertexNum( size_t resolution ){ return resolution*(resolution-1)+2; }
template< typename Real >
size_t SphericalGeometry::RegularGrid< Real >::FaceNum( size_t resolution ){ return resolution * resolution; }
template< typename Real >
size_t SphericalGeometry::RegularGrid< Real >::VertexIndex( int i , int j , size_t resolution )
{
	i = (i+(int)resolution) % (int)resolution;
	if( j==0 ) return 0u;
	else if( j==resolution ) return resolution*(resolution-1)+1;
	else return 1 + (j-1)*resolution+i;
}
template< typename Real >
size_t SphericalGeometry::RegularGrid< Real >::FaceIndex( int i , int j , size_t resolution )
{
	i = (i+(int)resolution) % (int)resolution;
	return i + j*resolution;
}
template< typename Real >
void SphericalGeometry::RegularGrid< Real >::FactorVertexIndex( size_t idx , int &i , int &j , size_t resolution )
{
	if( idx==0 ) i = j = 0;
	else if( idx==VertexNum(resolution)-1 ) i = 0 , j = (int)resolution;
	else
	{
		j = 1 + ( (int)idx-1 ) / (int)resolution;
		i = ((int)idx-1) % (int)resolution;
	}
}
template< typename Real >
void SphericalGeometry::RegularGrid< Real >::FactorFaceIndex( size_t idx , int &i , int &j , size_t resolution )
{
	i = (int)( idx % resolution );
	j = (int)( idx / resolution );
}
template< typename Real >
Point3D< Real > SphericalGeometry::RegularGrid< Real >::Vertex( size_t idx , size_t resolution )
{
	int i , j;
	FactorVertexIndex( idx , i , j , resolution );

	Real theta = (Real)( (2.*M_PI*i)/resolution );
	Real phi = (Real)( (PI*j)/resolution );
	return SpherePoint( Point2D< Real >( theta , phi ) );
}
template< typename Real >
Point3D< Real > SphericalGeometry::RegularGrid< Real >::FaceCenter( size_t idx , size_t resolution )
{
	int i , j;
	FactorFaceIndex( idx , i , j , resolution );

	Real theta = (Real)( (2.*M_PI*(i+0.5))/resolution );
	Real phi = (Real)( (PI*(j+0.5))/resolution );
	return SpherePoint( Point2D< Real >( theta , phi ) );
}
template< typename Real >
std::pair< Point3D< Real > , Point3D< Real > > SphericalGeometry::RegularGrid< Real >::VertexTangents( size_t idx , size_t resolution )
{
	int i , j;
	FactorVertexIndex( idx , i , j , resolution );

	if( j==0 || j==resolution )
	{
		Point3D< Real > v = Vertex( 0 , j , resolution );
		Point3D< Real > t1 = Point3D< Real >( 1 , 0 , 0 );
		Point3D< Real > t2 = Point3D< Real >::CrossProduct( v , t1 );
		return std::pair< Point3D< Real > , Point3D< Real > >( t1 , t2 );
	}
	else
	{
		Real theta = (Real)( (2.*M_PI*i)/resolution );
		Real phi   = (Real)( (   M_PI*j)/resolution );
		return SphereTangents( theta , phi );
	}
}

template< typename Real >
std::pair< Point3D< Real > , Point3D< Real > > SphericalGeometry::RegularGrid< Real >::FaceTangents( size_t idx , size_t resolution )
{
	int i , j;
	FactorFaceIndex( idx , i , j , resolution );

	Real theta = (Real)( (2.*M_PI*(i+0.5))/resolution );
	Real phi   = (Real)( (   M_PI*(j+0.5))/resolution );
	return SphereTangents( theta , phi );
}

template< typename Real >
typename SphericalGeometry::RegularGrid< Real >::FaceIndices SphericalGeometry::RegularGrid< Real >::Face( size_t idx , size_t resolution )
{
	FaceIndices face;
	int i , j;
	FactorFaceIndex( idx , i , j , resolution );

	if( j==0 )
	{
		face.size = 3;
		face[0] = VertexIndex( 0   , j  , resolution );
		face[1] = VertexIndex( i+1 , j+1 , resolution );
		face[2] = VertexIndex( i   , j+1 , resolution );
	}
	else if( j==resolution-1 )
	{
		face.size = 3;
		face[0] = VertexIndex( i   , j   , resolution );
		face[1] = VertexIndex( i+1 , j   , resolution );
		face[2] = VertexIndex( 0   , j+1 , resolution );
	}
	else
	{
		face.size = 4;
		face[0] = VertexIndex( i   , j   , resolution );
		face[1] = VertexIndex( i+1 , j   , resolution );
		face[2] = VertexIndex( i+1 , j+1 , resolution );
		face[3] = VertexIndex( i   , j+1 , resolution );
	}
	return face;
}
template< typename Real >
Real SphericalGeometry::RegularGrid< Real >::FaceMass( size_t idx , size_t resolution )
{
	int i , j;
	FactorFaceIndex( idx , i , j , resolution );
	Real dTheta = (Real)( (2.*M_PI)/resolution ) , dPhi = (Real)( PI/resolution );
	// dA = sin(\phi)

	// I = \int_0^dTheta \int_dPhi*j^dPhi*(j+1) sin(\phi) d phi d theta
	//   = dTheta [ cos( dPhi*j ) - cos( dPhi*(j+1) ) ]
	return (Real)( cos( dPhi*j ) - cos( dPhi*(j+1) ) ) * dTheta;
}
template< typename Real >
Real SphericalGeometry::RegularGrid< Real >::VertexMass( size_t idx , size_t resolution )
{
	int i , j;
	FactorVertexIndex( idx , i , j , resolution );
	if( j==0 || j==resolution ) return ( FaceMass( FaceIndex( 0 , 0 , resolution ) , resolution ) * resolution ) / 3;
	else return ( FaceMass( FaceIndex( 0 , j , resolution ) , resolution ) + FaceMass( FaceIndex( 0 , j-1 , resolution ) , resolution ) ) / 2;
}

template< typename Real >
template< unsigned int D >
SphericalGeometry::MonomialIntegralTable< Real , D > SphericalGeometry::RegularGrid< Real >::_GetMonomialIntegralTable( int j , size_t resolution , unsigned int mcSamples )
{
	MonomialIntegralTable< Real , D > mit;
	Real dPhi = (Real)( M_PI / resolution );
	mit.init( dPhi*j , dPhi*(j+1) , mcSamples );
	return mit;
}

template< typename Real >
void SphericalGeometry::RegularGrid< Real >::_BilinearFaceBasis( int j , size_t resolution , SinYPolynomial< Real , 1 , 1 , 0 > *polynomials )
{
	auto XFunction = []( Real x0 , Real x1 , Real v0 , Real v1 )
	{
		Polynomial< Real , 1 , 0 > p;
		Real a = ( v1 - v0 ) / ( x1 - x0 );
		Real b = v0 - a * x0;
		p[0][0] = b;
		p[1][0] = a;
		return p;
	};
	auto YFunction = []( Real y0 , Real y1 , Real v0 , Real v1 )
	{
		Polynomial< Real , 0 , 1 > p;
		Real a = ( v1 - v0 ) / ( y1 - y0 );
		Real b = v0 - a * y0;
		p[0][0] = b;
		p[0][1] = a;
		return p;
	};
	Real dTheta = (Real)( 2.*M_PI / resolution ) , dPhi = (Real)( M_PI / resolution );
	Polynomial< Real , 1 , 0 > x0 = XFunction( 0 , dTheta , 1 , 0 ) , x1 = XFunction( 0 , dTheta , 0 , 1 );
	Polynomial< Real , 0 , 1 > y0 = YFunction( dPhi*j , dPhi*(j+1) , 1 , 0 ) , y1 = YFunction( dPhi*j , dPhi*(j+1) , 0 , 1 );
	if( j==0 || j==resolution-1 )
	{
		if     ( j==0            ) polynomials[0] = y0*Polynomial< Real , 1 , 0 >(1) , polynomials[1] = y1*x1 , polynomials[2] = y1*x0;
		else if( j==resolution-1 ) polynomials[0] = y0*x0 , polynomials[1] = y0*x1 , polynomials[2] = y1*Polynomial< Real , 1 , 0 >(1);
	}
	else polynomials[0] = y0*x0 , polynomials[1] = y0*x1 , polynomials[2] = y1*x1 , polynomials[3] = y1*x0;
}
template< typename Real >
void SphericalGeometry::RegularGrid< Real >::_SetFaceMass( int j , size_t resolution , const MonomialIntegralTable< Real , 2 > &mit , Real mass[4][4] )
{
	BilinearFunctionType< Real > p[4];
	SinYPolynomial< Real , 0 , 0 , 1 > sine(1);
	_BilinearFaceBasis( j , resolution , p );
	unsigned int dim = ( j==0 || j==resolution-1 ) ? 3 : 4;
	for( unsigned int ii=0 ; ii<dim ; ii++ ) for( unsigned int jj=0 ; jj<dim ; jj++ ) mass[ii][jj] = (  p[ii] * p[jj] * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit );
}

template< typename Real >
void SphericalGeometry::RegularGrid< Real >::_SetFaceStiffness( int j , size_t resolution , const MonomialIntegralTable< Real , 2 > &mit , Real stiffness[4][4] )
{
	BilinearFunctionType< Real > p[4];
	BilinearGradientType< Real > g[4];
	SinYPolynomial< Real , 0 , 0 , 1 > sine(1);
	_BilinearFaceBasis( j , resolution , p );
	unsigned int dim = ( j==0 || j==resolution-1 ) ? 3 : 4;
	for( unsigned int d=0 ; d<dim ; d++ ) g[d] = p[d].gradient();
	for( unsigned int ii=0 ; ii<dim ; ii++ ) for( unsigned int jj=0 ; jj<dim ; jj++ )
		stiffness[ii][jj] = ( g[ii].x * g[jj].x * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit ) + ( g[ii].y * g[jj].y * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit );
}

template< typename Real >
void SphericalGeometry::RegularGrid< Real >::_SetFaceStiffnesses( int j , size_t resolution , const MonomialIntegralTable< Real , 2 > &mit , SquareMatrix< Real , 2 > stiffnesses[4][4] )
{
	BilinearFunctionType< Real > p[4];
	BilinearGradientType< Real > g[4];
	SinYPolynomial< Real , 0 , 0 , 1 > sine(1);
	_BilinearFaceBasis( j , resolution , p );
	unsigned int dim = ( j==0 || j==resolution-1 ) ? 3 : 4;
	for( unsigned int d=0 ; d<dim ; d++ ) g[d] = p[d].gradient();
	for( unsigned int ii=0 ; ii<dim ; ii++ ) for( unsigned int jj=0 ; jj<dim ; jj++ )
	{
		stiffnesses[ii][jj](0,0) = ( g[ii].x * g[jj].x * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit );
		stiffnesses[ii][jj](1,0) = ( g[ii].y * g[jj].x * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit );
		stiffnesses[ii][jj](0,1) = ( g[ii].x * g[jj].y * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit );
		stiffnesses[ii][jj](1,1) = ( g[ii].y * g[jj].y * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit );
	}
}

template< typename Real >
Real SphericalGeometry::RegularGrid< Real >::_Dot( const float &d1 , const float &d2 ){ return (Real)( d1*d2 ); }

template< typename Real >
Real SphericalGeometry::RegularGrid< Real >::_Dot( const double &d1 , const double &d2 ){ return ( Real )( d1*d2 ); }

template< typename Real >
template< typename Data >
Real SphericalGeometry::RegularGrid< Real >::_Dot( const Data &d1 , const Data &d2 ){ return Data::Dot( d1 , d2 ); }

/////////////////////////////
// RegularGrid::DualSignal //
/////////////////////////////
template< typename Real >
template< typename Data >
void SphericalGeometry::RegularGrid< Real >::DualSignal< Data >::_setFromPrimal( const Data *values , size_t resolution )
{
	setResolution( resolution );
#pragma omp parallel for
	for( int i=0 ; i<RegularGrid< Real >::FaceNum( resolution ) ; i++ )
	{
		RegularGrid< Real >::FaceIndices face = RegularGrid< Real >::Face( i , resolution );
		Data value{};
		for( int i=0 ; i<face.size ; i++ ) value += values[ face[i] ];
		_signal[i] = value / (Real)face.size;
	}
}

template< typename Real >
template< typename Data >
SparseMatrix< Real , int > SphericalGeometry::RegularGrid< Real >::DualSignal< Data >::MassMatrix( size_t resolution )
{
	SparseMatrix< Real , int > M;
	M.resize( FaceNum( resolution ) );
#pragma omp parallel for
	for( int j=0 ; j<resolution ; j++ )
	{
		Real mass = RegularGrid< Real >::FaceMass( 0 , j , resolution );
		for( int i=0 ; i<resolution ; i++ )
		{
			size_t idx = RegularGrid< Real >::FaceIndex( i , j , resolution );
			M.SetRowSize( idx , 1 );
			M[idx][0] = MatrixEntry< Real , int >( (int)idx , mass );
		}
	}
	return M;
}

template< typename Real >
template< typename Data >
Real SphericalGeometry::RegularGrid< Real >::DualSignal< Data >::squareNorm( void ) const
{
	Real l2 = 0;
#pragma omp parallel for reduction( + : l2 )
	for( int j=0 ; j<_resolution ; j++ )
	{
		Real mass = RegularGrid< Real >::FaceMass( 0 , j , _resolution );
		for( int i=0 ; i<_resolution ; i++ ) l2 += _Dot( _signal[i] , _signal[i] ) * mass;
	}
	return l2;
}

template< typename Real >
template< typename Data >
Data SphericalGeometry::RegularGrid< Real >::DualSignal< Data >::operator() ( Point3D< Real > p ) const
{
	Point2D< Real > thetaPhi = SphereParameters( p );
	thetaPhi[0] *= (Real)( _resolution / ( 2.*M_PI ) );
	thetaPhi[1] *= (Real)( _resolution / M_PI );

	int theta = (int)floor(thetaPhi[0]) , phi = (int)floor(thetaPhi[1]);
	return _signal[ FaceIndex( theta , phi , _resolution ) ];
}

///////////////////////////////
// RegularGrid::PrimalSignal //
///////////////////////////////
template< typename Real >
template< typename Data >
Data SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::operator() ( Point3D< Real > p ) const
{
	Point2D< Real > thetaPhi = SphereParameters( p );
	thetaPhi[0] *= (Real)( _resolution / ( 2.*M_PI ) );
	thetaPhi[1] *= (Real)( _resolution / M_PI );

	int theta = (int)floor(thetaPhi[0]) , phi = (int)floor(thetaPhi[1]);
	if( phi==_resolution ) phi--;
	Real dTheta = thetaPhi[0]-theta , dPhi = thetaPhi[1]-phi;
	return
		_signal[ VertexIndex( theta+0 , phi+0 , _resolution ) ] * (Real)( ( 1.-dTheta ) * ( 1.-dPhi ) ) +
		_signal[ VertexIndex( theta+0 , phi+1 , _resolution ) ] * (Real)( ( 1.-dTheta ) * (    dPhi ) ) +
		_signal[ VertexIndex( theta+1 , phi+0 , _resolution ) ] * (Real)( (    dTheta ) * ( 1.-dPhi ) ) +
		_signal[ VertexIndex( theta+1 , phi+1 , _resolution ) ] * (Real)( (    dTheta ) * (    dPhi ) ) ;
}

template< typename Real >
template< typename Data >
typename SphericalGeometry::RegularGrid< Real >::template PrimalSignal< Data > SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::upSample( void ) const
{
	PrimalSignal up( 2 *_resolution );
	up._signal[0] = _signal[0];
	up._signal.back() = _signal.back();
#pragma omp parallel for
	for( int i=1 ; i<up.size()-1 ; i++ )
	{
		int ii=(i-1)%(int)up._resolution , jj = (i-1)/(int)up._resolution + 1;
		if( ii&1 && jj&1 ) up[i] = ( _signal[ RegularGrid< Real >::VertexIndex(ii/2,jj/2,_resolution) ] + _signal[ RegularGrid< Real >::VertexIndex(ii/2+1,jj/2,_resolution) ] + _signal[ RegularGrid< Real >::VertexIndex(ii/2+1,jj/2+1,_resolution) ] + _signal[ VertexIndex(ii/2,jj/2+1,_resolution) ] ) / 4;
		else if( ii&1 ) up[i] = ( _signal[ RegularGrid< Real >::VertexIndex(ii/2,jj/2,_resolution) ] + _signal[ RegularGrid< Real >::VertexIndex(ii/2+1,jj/2,_resolution) ] ) / 2;
		else if( jj&1 ) up[i] = ( _signal[ RegularGrid< Real >::VertexIndex(ii/2,jj/2,_resolution) ] + _signal[ RegularGrid< Real >::VertexIndex(ii/2,jj/2+1,_resolution) ] ) / 2;
		else up[i] = _signal[ RegularGrid< Real >::VertexIndex(ii/2,jj/2,_resolution) ];
	}

	return up;
}

template< typename Real >
template< typename Data >
typename SphericalGeometry::RegularGrid< Real >::template PrimalSignal< Data > SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::downSample( void ) const
{
	if( _resolution&1 ) fprintf( stderr , "[ERROR] Cannot down-sample odd resolution: %d\n" , (int)_resolution ) , exit( 0 );
	PrimalSignal down( _resolution/2 );
	{
		down[0] = 0;
		for( int i=0 ; i<_resolution ; i++ ) down[0] += _signal[1+i];
		down[0] /= 2*_resolution;
		down[0] += _signal[0] / 2;
	}
	{
		down._signal.back() = 0;
		for( int i=0 ; i<_resolution ; i++ ) down._signal.back() += _signal[ RegularGrid< Real >::VertexIndex(i,(int)_resolution-1,_resolution) ];
		down._signal.back() /= 2*_resolution;
		down._signal.back() += _signal.back() / 2;
	}
#pragma omp parallel for
	for( int i=1 ; i<down.size()-1 ; i++ )
	{
		int ii=(i-1)%(int)down._resolution , jj = (i-1)/(int)down._resolution + 1;
		down[i] = ( _signal[ RegularGrid< Real >::VertexIndex(2*ii,2*jj,_resolution) ]*4 + _signal[ RegularGrid< Real >::VertexIndex(2*ii+1,2*jj,_resolution) ] + _signal[ RegularGrid< Real >::VertexIndex(2*ii-1+(int)_resolution,2*jj,_resolution) ] + _signal[ RegularGrid< Real >::VertexIndex(2*ii,2*jj-1,_resolution) ] + _signal[ RegularGrid< Real >::VertexIndex(2*ii,2*jj+1,_resolution) ] ) / 8;
	}

	return down;
}

template< typename Real >
template< typename Data >
template< typename ColorFunctor >
void SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::write( const char *fileName , const ColorFunctor &colorFunctor ) const
{
	std::vector< std::vector< int > > polygons( RegularGrid< Real >::FaceNum( _resolution ) );
	std::vector< PlyColorVertex< float > > vertices( RegularGrid< Real >::VertexNum( _resolution ) );
	for( int i=0 ; i<RegularGrid< Real >::VertexNum( _resolution ) ; i++ )
	{
		vertices[i].point = Point3D< float >( Vertex( i , _resolution ) );
		vertices[i].color = Point3D< float >( colorFunctor( _signal[i] ) );
	}
	for( int i=0 ; i<RegularGrid< Real >::FaceNum( _resolution ) ; i++ )
	{
		FaceIndices face = RegularGrid< Real >::Face( i , _resolution );
		polygons[i].resize( face.size );
		for( int j=0 ; j<face.size ; j++ ) polygons[i][j] = (int)face[j];
	}
	PlyWritePolygons( fileName , vertices , polygons , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , PLY_BINARY_NATIVE );
}

template< typename Real >
template< typename Data >
template< typename ColorFunctor >
void SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::write( const char *fileName , const ColorFunctor &colorFunctor , const DualSignal< Point3D< Real > > &vectorField ) const
{
	std::vector< PlyVFFace< float > > polygons( RegularGrid< Real >::FaceNum( _resolution ) );
	std::vector< PlyColorVertex< float > > vertices( RegularGrid< Real >::VertexNum( _resolution ) );
	for( int i=0 ; i<RegularGrid< Real >::VertexNum( _resolution ) ; i++ )
	{
		vertices[i].point = Point3D< float >( Vertex( i , _resolution ) );
		vertices[i].color = Point3D< float >( colorFunctor( _signal[i] ) );
	}
	for( int i=0 ; i<RegularGrid< Real >::FaceNum( _resolution ) ; i++ )
	{
		FaceIndices face = RegularGrid< Real >::Face( i , _resolution );
		polygons[i].resize( (int)face.size );
		for( int j=0 ; j<face.size ; j++ ) polygons[i][j] = (int)face[j];
		polygons[i].v = Point3D< float >( vectorField[i] );
	}
	PlyWritePolygons( fileName , vertices , polygons , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , PlyVFFace< float >::WriteProperties , PlyVFFace< float >::WriteComponents , PLY_BINARY_NATIVE );
}

template< typename Real >
template< typename Data >
template< typename ColorFunctor , bool Gradient , bool JGradient >
void SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::write( const char *fileName , const ColorFunctor &colorFunctor , const PrimalVectorField< Gradient , JGradient > &vectorField ) const
{
	std::vector< PlyVFFace< float > > polygons( RegularGrid< Real >::FaceNum( _resolution ) );
	std::vector< PlyColorVertex< float > > vertices( RegularGrid< Real >::VertexNum( _resolution ) );
	for( int i=0 ; i<RegularGrid< Real >::VertexNum( _resolution ) ; i++ )
	{
		vertices[i].point = Point3D< float >( Vertex( i , _resolution ) );
		vertices[i].color = Point3D< float >( colorFunctor( _signal[i] ) );
	}
	for( int i=0 ; i<RegularGrid< Real >::FaceNum( _resolution ) ; i++ )
	{
		FaceIndices face = RegularGrid< Real >::Face( i , _resolution );
		polygons[i].resize( (int)face.size );
		for( int j=0 ; j<face.size ; j++ ) polygons[i][j] = (int)face[j];
		polygons[i].v = Point3D< float >( vectorField.faceValue(i) );
	}
	PlyWritePolygons( fileName , vertices , polygons , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , PlyVFFace< float >::WriteProperties , PlyVFFace< float >::WriteComponents , PLY_BINARY_NATIVE );
}

template< typename Real >
template< typename Data >
void SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::_setFromDual( const Data *values , size_t resolution )
{
	setResolution( resolution );
#pragma omp parallel for
	for( int i=0 ; i<_signal.size() ; i++ ) _signal[i] *= 0;
#pragma omp parallel for
	for( int j=0 ; j<resolution ; j++ )
	{
		Real faceArea = RegularGrid< Real >::FaceMass( RegularGrid< Real >::FaceIndex( 0 , j , resolution ) , resolution );
		for( int i=0 ; i<resolution ; i++ )
		{
			size_t f = RegularGrid< Real >::FaceIndex( i , j , resolution );
			FaceIndices face = RegularGrid< Real >::Face( f , resolution );
			Data d = values[f] * faceArea;
			for( int k=0 ; k<face.size ; k++ ) AddAtomic( _signal[ face[k] ] , d );
		}
	}
	_signal[0] /= RegularGrid< Real >::VertexMass( 0 , resolution ) * 3;
	_signal.back() /= RegularGrid< Real >::VertexMass( VertexIndex( 0 , (int)resolution , resolution ) , resolution ) * 3;
#pragma omp parallel for
	for( int j=1 ; j<resolution ; j++ )
	{
		Real vertexArea = RegularGrid< Real >::VertexMass( RegularGrid< Real >::VertexIndex( 0 , j , resolution ) , resolution ) * 4;
		for( int i=0 ; i<resolution ; i++ ) _signal[ RegularGrid< Real >::VertexIndex( i , j , resolution ) ] /= vertexArea;
	}
}

template< typename Real >
template< typename Data >
Matrix< Real , 4 , 2 > SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::_GradientOperator( int j , size_t resolution )
{
	SquareMatrix< Real , 2 > M;
	Matrix< Real , 4 , 2 > B;
	// Given vertices {v[]}, values {w[]}, and tangents {t1,t2}, we want to solve for a and b minimizing:
	// E(a,b) = \sum_i ( w[i+1]-w[i] - < v[i+1]-v[i] , a*t1 + b*t2 > )^2
	//        = \sum_i ( w[i+1]-w[i] )^2 + < v[i+1]-v[i] , a*t1 + b*t2 >^2 - 2 * (w[i+1]-w[i])*< v[i+1]-v[i] , a*t1 + b*t2 >
	// Taking the partial derivatives:
	// dE/da  = \sum_i 2 < v[i+1]-v[i] , a*t1 + b*t2 > * < v[i+1]-v[i] , t1 > - 2 * (w[i+1]-w[i])*< v[i+1]-v[i] , t1 >
	// dE/db  = \sum_i 2 < v[i+1]-v[i] , a*t1 + b*t2 > * < v[i+1]-v[i] , t2 > - 2 * (w[i+1]-w[i])*< v[i+1]-v[i] , t2 >
	// Setting to zero gives:
	//       \sum_i (w[i+1]-w[i])*< v[i+1]-v[i] , t1 > = \sum_i < v[i+1]-v[i] , a*t1 + b*t2 > * < v[i+1]-v[i] , t1 >
	//       \sum_i (w[i+1]-w[i])*< v[i+1]-v[i] , t2 > = \sum_i < v[i+1]-v[i] , a*t1 + b*t2 > * < v[i+1]-v[i] , t2 >

	Point3D< Real > v[4];
	FaceIndices face = RegularGrid< Real >::Face( RegularGrid< Real >::FaceIndex( 0 , j , resolution ) , resolution );
	for( int j=0 ; j<face.size ; j++ ) v[j] = RegularGrid< Real >::Vertex( face[j] , resolution );


	std::pair< Point3D< Real > , Point3D< Real > > tangents = RegularGrid< Real >::FaceTangents( 0 , j , resolution );
	Point3D< Real > v1 = tangents.first , v2 = tangents.second;
	for( int i=0 ; i<face.size ; i++ )
	{
		Point3D< Real > dv =  v[(i+1)%face.size]-v[i];
		M(0,0) += Point3D< Real >::Dot( dv , v1 ) * Point3D< Real >::Dot( dv , v1 );
		M(0,1) += Point3D< Real >::Dot( dv , v2 ) * Point3D< Real >::Dot( dv , v1 );
		M(1,1) += Point3D< Real >::Dot( dv , v2 ) * Point3D< Real >::Dot( dv , v2 );
		B(i,0) -= Point3D< Real >::Dot( dv , v1 );
		B(i,1) -= Point3D< Real >::Dot( dv , v2 );
		B((i+1)%face.size,0) += Point3D< Real >::Dot( dv , v1 );
		B((i+1)%face.size,1) += Point3D< Real >::Dot( dv , v2 );
	}
	M(1,0) = M(0,1);
	return M.inverse() * B;
}

template< typename Real >
template< typename Data >
void SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::_setGradients( Point3D< Data > *gradientValues ) const
{
#pragma omp parallel for
	for( int j=0 ; j<_resolution ; j++ )
	{
		Point2D< Data > ab;
		Matrix< Real , 4 ,  2 > G = _GradientOperator( j , _resolution );
		for( int i=0 ; i<_resolution ; i++ )
		{
			size_t idx = RegularGrid< Real >::FaceIndex( i , j , _resolution );
			FaceIndices face = RegularGrid< Real >::Face( idx , _resolution );
			std::pair< Point3D< Real > , Point3D< Real > > tangents = RegularGrid< Real >::FaceTangents( i , j , _resolution );
			Point< Data , 4 > p;
			for( int f=0 ; f<face.size ; f++ ) p[f] = _signal[ face[f] ];
			ab = G * p;
			Point3D< Real > vertices[3];
			for( int j=0 ; j<3 ; j++ ) vertices[j] = RegularGrid< Real >::Vertex( face[j] , _resolution );
			for( int j=0 ; j<3 ; j++ ) gradientValues[idx][j] = ab[0] * tangents.first[j] + ab[1] * tangents.second[j];
		}
	}
}


template< typename Real >
template< typename Data >
SparseMatrix< Real , int > SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::RestrictionMatrix( size_t coarseResolution )
{
	size_t fineResolution = 2*coarseResolution;
	size_t lowSZ = RegularGrid< Real >::VertexNum( coarseResolution );
	size_t highSZ = RegularGrid< Real >::VertexNum( fineResolution );

	SparseMatrix< Real , int > R;
	R.resize( lowSZ );
	// Do the poles
	{
		R.SetRowSize( 0 , 1+fineResolution );
		R[0][0] = MatrixEntry< Real , int >( 0 , 1 );
		for( int i=0 ; i<fineResolution ; i++ ) R[0][i+1] = MatrixEntry< Real , int >( i+1 , (Real)0.5 );

		R.SetRowSize( lowSZ-1 , 1+fineResolution );
		R[lowSZ-1][0] = MatrixEntry< Real , int >( (int)highSZ-1 , 1 );
		for( int i=0 ; i<fineResolution ; i++ ) R[lowSZ-1][i+1] = MatrixEntry< Real , int >( (int)highSZ-2-i, (Real)0.5 );
	}
#pragma omp parallel for
	for( int j=1 ; j<coarseResolution ; j++ ) for( int i=0 ; i<coarseResolution ; i++ )
	{
		size_t idx = VertexIndex( i , j , coarseResolution );
		R.SetRowSize( idx , 9 );
		R[idx][0] = MatrixEntry< Real , int >( (int)VertexIndex( 2*i+0 , 2*j+0 , fineResolution ) , (Real)1.00 );

		R[idx][1] = MatrixEntry< Real , int >( (int)VertexIndex( 2*i-1 , 2*j+0 , fineResolution ) , (Real)0.50 );
		R[idx][2] = MatrixEntry< Real , int >( (int)VertexIndex( 2*i+1 , 2*j+0 , fineResolution ) , (Real)0.50 );
		R[idx][3] = MatrixEntry< Real , int >( (int)VertexIndex( 2*i+0 , 2*j-1 , fineResolution ) , (Real)0.50 );
		R[idx][4] = MatrixEntry< Real , int >( (int)VertexIndex( 2*i+0 , 2*j+1 , fineResolution ) , (Real)0.50 );

		R[idx][5] = MatrixEntry< Real , int >( (int)VertexIndex( 2*i-1 , 2*j-1 , fineResolution ) , (Real)0.25 );
		R[idx][6] = MatrixEntry< Real , int >( (int)VertexIndex( 2*i+1 , 2*j-1 , fineResolution ) , (Real)0.25 );
		R[idx][7] = MatrixEntry< Real , int >( (int)VertexIndex( 2*i-1 , 2*j+1 , fineResolution ) , (Real)0.25 );
		R[idx][8] = MatrixEntry< Real , int >( (int)VertexIndex( 2*i+1 , 2*j+1 , fineResolution ) , (Real)0.25 );
	}

	return R;
}

template< typename Real >
template< typename Data >
template< typename Entry >
void SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::_InitSystemMatrix( size_t resolution , SparseMatrix< Entry , int > &M )
{
	size_t sz = RegularGrid< Real >::VertexNum( resolution );

	// Allocte the memory for the sparse matrices
	M.resize( sz );

#pragma omp parallel for
	for( int j=0 ; j<=resolution ; j++ )
	{
		if( j==0 )
		{
			size_t idx = VertexIndex( 0 , j , resolution ) , _idx = 0;
			M.SetRowSize( (int)idx , resolution+1 );
			M[idx][_idx++] = MatrixEntry< Entry , int >( (int)idx , Entry() );
			for( int i=0 ; i<resolution ; i++ ) M[idx][_idx++] = MatrixEntry< Entry , int >( (int)VertexIndex( i , j+1 , resolution ) , Entry() );
		}
		else if( j==resolution )
		{
			size_t idx = VertexIndex( 0 , j , resolution ) , _idx = 0;
			M.SetRowSize( (int)idx , resolution+1 );
			for( int i=0 ; i<resolution ; i++ ) M[idx][_idx++] = MatrixEntry< Entry , int >( (int)VertexIndex( i , j-1 , resolution ) , Entry() );
			M[idx][_idx++] = MatrixEntry< Entry , int >( (int)idx , Entry() );
		}
		else if( j==1 )
		{
			for( int i=0 ; i<resolution ; i++ )
			{
				size_t idx = VertexIndex( i , j , resolution ) , _idx = 0;
				M.SetRowSize( (int)idx , 7 );
				M[idx][_idx++] = MatrixEntry< Entry , int >( (int)VertexIndex( i , j-1 , resolution ) , Entry() );
				for( int jj=0 ; jj<=1 ; jj++ ) for( int ii=-1 ; ii<=1 ; ii++ ) M[idx][_idx++] = MatrixEntry< Entry , int >( (int)VertexIndex( i+ii , j+jj , resolution ) , Entry() );
			}
		}
		else if( j==resolution-1 )
		{
			for( int i=0 ; i<resolution ; i++ )
			{
				size_t idx = VertexIndex( i , j , resolution ) , _idx = 0;
				M.SetRowSize( (int)idx , 7 );
				for( int jj=-1 ; jj<=0 ; jj++ ) for( int ii=-1 ; ii<=1 ; ii++ ) M[idx][_idx++] = MatrixEntry< Entry , int >( (int)VertexIndex( i+ii , j+jj , resolution ) , Entry() );
				M[idx][_idx++] = MatrixEntry< Entry , int >( (int)VertexIndex( i , j+1 , resolution ) , Entry() );
			}
		}
		else
		{
			for( int i=0 ; i<resolution ; i++ )
			{
				size_t idx = VertexIndex( i , j , resolution ) , _idx = 0;
				M.SetRowSize( (int)idx , 9 );
				for( int jj=-1 ; jj<=1 ; jj++ ) for( int ii=-1 ; ii<=1 ; ii++ ) M[idx][_idx++] = MatrixEntry< Entry , int >( (int)VertexIndex( i+ii , j+jj , resolution ) , Entry() );
			}
		}
	}
}

template< typename Real >
template< typename Data >
size_t SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::_RowOffset( size_t idx1 , size_t idx2 , size_t resolution )
{
	auto ThetaNeighbors = [] ( int i , int j , size_t resolution )
	{
		if( ( j>=i-1 && j<=i+1 ) || ( ( j>=i+resolution-1 && j<=i+resolution+1 ) ) || ( j+resolution>=i-1 && j+resolution<=i+1 ) ) return true;
		return false;
	};

	int i1 , i2 , j1 , j2;
	FactorVertexIndex( idx1 , i1 , j1 , resolution );
	FactorVertexIndex( idx2 , i2 , j2 , resolution );
	if( j1==0 )
	{
		if     ( j2==0 ) return 0;
		else if( j2==1 ) return i2+1;
	}
	else if( j1==resolution )
	{
		if     ( j2==resolution-1 ) return i2;
		else if( j2==resolution   ) return (int)resolution;
	}
	else if( j1==1 )
	{
		if     ( j2==0 ) return 0;
		else if( j2==1 ){ if( ThetaNeighbors( i1 , i2 , resolution ) ) return (int)( ( (i2-i1+1+resolution )%resolution ) + 1 ); }
		else if( j2==2 ){ if( ThetaNeighbors( i1 , i2 , resolution ) ) return (int)( ( (i2-i1+1+resolution )%resolution ) + 4 ); }
	}
	else if( j1==resolution-1 )
	{
		if     ( j2==resolution-2 ){ if( ThetaNeighbors( i1 , i2 , resolution ) ) return (int)( ( (i2-i1+1+resolution )%resolution ) + 0 ); }
		else if( j2==resolution-1 ){ if( ThetaNeighbors( i1 , i2 , resolution ) ) return (int)( ( (i2-i1+1+resolution )%resolution ) + 3 ); }
		else if( j2==resolution   ) return 6;
	}
	else{ if( j2>=j1-1 && j2<=j1+1 && ThetaNeighbors( i1 , i2 , resolution ) ) return (j2-j1+1)*3 + (int)( ( (i2-i1+1)+resolution )%resolution ); }
	fprintf( stderr , "[ERROR] Vertices are not neighbors: (%d %d) (%d %d)\n" , i1 , j1 , i2 , j2 ) , exit( 0 );
	return -1;
}

template< typename Real >
template< typename Data >
SparseMatrix< Real , int > SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::MassMatrix( size_t resolution , unsigned int mcSamples )
{
	SparseMatrix< Real , int > M;
	size_t sz = VertexNum( resolution );

	_InitSystemMatrix( resolution , M );

	// Set the matrix coefficients by iterating over the faces
#pragma omp parallel for
	for( int j=0 ; j<resolution ; j++ )
	{
		MonomialIntegralTable< Real , 2 > mit = _GetMonomialIntegralTable< 2 >( j , resolution , mcSamples );
		Real m[4][4];
		_SetFaceMass( j , resolution , mit , m );
		for( int i=0 ; i<resolution ; i++ )
		{
			FaceIndices face = Face( i , j , resolution );
			for( int ii=0 ; ii<face.size ; ii++ ) for( int jj=0 ; jj<face.size ; jj++ )
			{
				size_t offset = _RowOffset( face[ii] , face[jj] , resolution );
#pragma omp atomic
				M[ face[ii] ][ offset ].Value += m[ii][jj];
			}
		}
	}
	return M;
}
template< typename Real >
template< typename Data >
SparseMatrix< Real , int > SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::StiffnessMatrix( size_t resolution , unsigned int mcSamples )
{
	SparseMatrix< Real , int > S;
	size_t sz = VertexNum( resolution );

	_InitSystemMatrix( resolution , S );

	// Set the matrix coefficients by iterating over the faces
#pragma omp parallel for
	for( int j=0 ; j<resolution ; j++ )
	{
		MonomialIntegralTable< Real , 2 > mit = _GetMonomialIntegralTable< 2 >( j , resolution , mcSamples );
		Real s[4][4];
		_SetFaceStiffness( j , resolution , mit , s );
		for( int i=0 ; i<resolution ; i++ )
		{
			FaceIndices face = Face( i , j , resolution );
			for( int ii=0 ; ii<face.size ; ii++ ) for( int jj=0 ; jj<face.size ; jj++ )
			{
				size_t offset = _RowOffset( face[ii] , face[jj] , resolution );
#pragma omp atomic
				S[ face[ii] ][ offset ].Value += s[ii][jj];
			}
		}
	}
	return S;
}

template< typename Real >
template< typename Data >
template< typename MetricTensor2x2 >
SparseMatrix< Real , int > SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::MassMatrix( size_t resolution , unsigned int mcSamples , const MetricTensor2x2 &metricTensor )
{
	SparseMatrix< Real , int > M;
	size_t sz = VertexNum( resolution );

	_InitSystemMatrix( resolution , M );

	// Set the matrix coefficients by iterating over the faces
#pragma omp parallel for
	for( int j=0 ; j<resolution ; j++ )
	{
		MonomialIntegralTable< Real , 2 > mit = _GetMonomialIntegralTable< 2 >( j , resolution , mcSamples );
		Real m[4][4];
		_SetFaceMass( j , resolution , mit , m );
		for( int i=0 ; i<resolution ; i++ )
		{
			Real scl = (Real)sqrt( metricTensor(i,j).determinant() );
			FaceIndices face = Face( i , j , resolution );
			for( int ii=0 ; ii<face.size ; ii++ ) for( int jj=0 ; jj<face.size ; jj++ )
			{
				size_t offset = _RowOffset( face[ii] , face[jj] , resolution );
#pragma omp atomic
				M[ face[ii] ][ offset ].Value += m[ii][jj]*scl;
			}
		}
	}
	return M;
}

template< typename Real >
template< typename Data >
template< typename MetricTensor2x2 >
SparseMatrix< Real , int > SphericalGeometry::RegularGrid< Real >::PrimalSignal< Data >::StiffnessMatrix( size_t resolution , unsigned int mcSamples , const MetricTensor2x2 &metricTensor )
{
	SparseMatrix< Real , int > S;
	size_t sz = VertexNum( resolution );

	_InitSystemMatrix( resolution , S );

	// Set the matrix coefficients by iterating over the faces
#pragma omp parallel for
	for( int j=0 ; j<resolution ; j++ )
	{
		MonomialIntegralTable< Real , 2 > mit = _GetMonomialIntegralTable< 2 >( j , resolution , mcSamples );
		SquareMatrix< Real , 2 > s[4][4];
		_SetFaceStiffnesses( j , resolution , mit , s );
		for( int i=0 ; i<resolution ; i++ )
		{
			SquareMatrix< Real , 2 > g = metricTensor(i,j);
			FaceIndices face = Face( i , j , resolution );
			for( int ii=0 ; ii<face.size ; ii++ ) for( int jj=0 ; jj<face.size ; jj++ )
			{
				size_t offset = _RowOffset( face[ii] , face[jj] , resolution );
#pragma omp atomic
				S[ face[ii] ][ offset ].Value += SquareMatrix< Real , 2  >::Dot( s[ii][jj] , g );
			}
		}
	}
	return S;
}


//////////////////////////////////////////////////////////////
// RegularGrid::PrimalVectorField::MassAndStiffnessOperator //
//////////////////////////////////////////////////////////////
template< typename Real >
template< bool Gradient , bool JGradient >
template< typename Data >
SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::MassAndStiffnessOperator< Data >::MassAndStiffnessOperator( size_t resolution , unsigned int mcSamples )
{
	for( int b=0 ; b<Blocks ; b++ ) _mWeights[b] = _sWeights[b] = 0;
	SparseMatrix< Real , int > M = PrimalSignal< Real >::MassMatrix( resolution , mcSamples );
	_S = PrimalSignal< Real >::StiffnessMatrix( resolution , mcSamples );
	_invM.resize( M.rows );
	_scratch.resize( M.rows );
#pragma omp parallel for
	for( int i=0 ; i<M.rows ; i++ )
	{
		_invM.SetRowSize( i , 1 );
		Real rowSum = 0;
		for( int j=0 ; j<M.rowSizes[i] ; j++ ) rowSum += M[i][j].Value;
		_invM[i][0] = MatrixEntry< Real , int >( i , (Real)1./rowSum );
	}
}

template< typename Real >
template< bool Gradient , bool JGradient >
template< typename Data >
template< typename MetricTensor2x2 > 
SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::MassAndStiffnessOperator< Data >::MassAndStiffnessOperator( size_t resolution , unsigned int mcSamples , const MetricTensor2x2 &metricTensor )
{
	for( int b=0 ; b<Blocks ; b++ ) _mWeights[b] = _sWeights[b] = 0;
	SparseMatrix< Real , int > M = PrimalSignal< Real >::MassMatrix( resolution , mcSamples , metricTensor );
	_S = PrimalSignal< Real >::StiffnessMatrix( resolution , mcSamples , metricTensor );
	_invM.resize( M.rows );
	_scratch.resize( M.rows );
#pragma omp parallel for
	for( int i=0 ; i<M.rows ; i++ )
	{
		_invM.SetRowSize( i , 1 );
		Real rowSum = 0;
		for( int j=0 ; j<M.rowSizes[i] ; j++ ) rowSum += M[i][j].Value;
		_invM[i][0] = MatrixEntry< Real , int >( i , (Real)1./rowSum );
	}
}
template< typename Real >
template< bool Gradient , bool JGradient >
template< typename Data >
SparseMatrix< Real , int > SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::MassAndStiffnessOperator< Data >::toMatrix( void ) const
{
	if( Gradient!=JGradient ) return _S * _mWeights[0] + _S * _invM * _S * _sWeights[0];
	else if( Gradient && JGradient )
	{
		SparseMatrix< Real , int > M , divM , curlM;
		divM = _S * _mWeights[0] + _S * _invM * _S * _sWeights[0];
		curlM = _S * _mWeights[1] + _S * _invM * _S * _sWeights[1];
		M.resize( divM.rows + curlM.rows );
		for( int i=0 ; i<divM.rows ; i++ )
		{
			M.SetRowSize( i , divM.rowSizes[i] );
			for( int j=0 ; j<divM.rowSizes[i] ; j++ ) M[i][j] = divM[i][j];
		}
		for( int i=0 ; i<curlM.rows ; i++ )
		{
			M.SetRowSize( i+divM.rows , curlM.rowSizes[i] );
			for( int j=0 ; j<curlM.rowSizes[i] ; j++ ) M[i+divM.rows][j] = MatrixEntry< Real , int >( curlM[i][j].N + (int)divM.rows , curlM[i][j].Value );
		}
		return M;
	}
}

template< typename Real >
template< bool Gradient , bool JGradient >
template< typename Data >
void SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::MassAndStiffnessOperator< Data >::Multiply( ConstPointer( Data ) x , Pointer( Data ) Mx ) const
{
	auto SetSmoothness = []( const SparseMatrix< Real , int > &S , const SparseMatrix< Real , int > &invM , ConstPointer( Data ) x , Pointer( Data ) Mx , Pointer( Data ) scratch , Real weight )
	{
#pragma omp parallel for
		for( int i=0 ; i<S.rows ; i++ )
		{
			Data d{};
			for( int j=0 ; j<S.rowSizes[i] ; j++ ) d += x[ S[i][j].N ] * S[i][j].Value;
			scratch[i] = d * invM[i][0].Value;
		}
#pragma omp parallel for
		for( int i=0 ; i<S.rows ; i++ )
		{
			Data d{};
			for( int j=0 ; j<S.rowSizes[i] ; j++ ) d += scratch[ S[i][j].N ] * S[i][j].Value;
			Mx[i] = d * weight;
		}
	};
	auto AddMass = []( const SparseMatrix< Real , int > &S , ConstPointer( Data ) x , Pointer( Data ) Mx , Real weight )
	{
#pragma omp parallel for
		for( int i=0 ; i<S.rows ; i++ )
		{
			Data d{};
			for( int j=0 ; j<S.rowSizes[i] ; j++ ) d += x[ S[i][j].N ] * S[i][j].Value;
			Mx[i] += d * weight;
		}
	};

	size_t sz = _S.rows;

	if( Gradient )
	{
		unsigned int idx = 0;
		size_t off = idx*sz;
		if( _sWeights[idx] ) SetSmoothness( _S , _invM , x+off , Mx+off , GetPointer( const_cast< std::vector< Data > & >( _scratch ) ) , _sWeights[idx] );
		else
#pragma omp parallel for
			for( int i=0 ; i<sz ; i++ ) Mx[i+off] = Data{};
		if( _mWeights[idx] ) AddMass( _S , x+off , Mx+off , _mWeights[idx]  );
	}
	if( JGradient )
	{
		unsigned int idx = Blocks-1;
		size_t off = idx*sz;
		if( _sWeights[idx] ) SetSmoothness( _S , _invM , x+off , Mx+off , GetPointer( const_cast< std::vector< Data > & >( _scratch ) ) , _sWeights[idx] );
		else
#pragma omp parallel for
			for( int i=0 ; i<sz ; i++ ) Mx[i+off] = Data{};
		if( _mWeights[idx] ) AddMass( _S , x+off , Mx+off , _mWeights[idx]  );
	}
}

////////////////////////////////////
// RegularGrid::PrimalVectorField //
////////////////////////////////////
template< typename Real >
template< bool Gradient , bool JGradient >
void SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::setVertexValues( PrimalSignal< Point3D< Real > > &vf ) const
{
	vf.setResolution( _resolution );
#pragma omp parallel for
	for( int i=0 ; i<vf.size() ; i++ ) vf[i] = vertexValue(i);
}

template< typename Real >
template< bool Gradient , bool JGradient >
Point3D< Real > SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::_sample( int i , int j , Real di , Real dj ) const
{
	Real theta = (Real)( (i+di) * 2 * M_PI / _resolution );
	Real phi = (Real)( (j+dj) * M_PI / _resolution );
	std::pair< Point3D< Real > , Point3D< Real > > tangents = SphereTangents( theta , phi );

	Real sinPhi = (Real)sin( phi );
	Real thetaScale = (Real)( _resolution / (2*M_PI) ) / sinPhi , phiScale = (Real)( _resolution / M_PI );

	// Functions -> Gradients:
	// (1-theta*resolution/(2PI)) * (1-phi*resolution/PI) -> -(resolution / (2PI) ) * (1-phi) * T / sin(phi) - (resolution / PI ) * (1-theta) * P
	// (  theta*resolution/(2PI)) * (1-phi*resolution/PI) ->  (resolution / (2PI) ) * (1-phi) * T / sin(phi) - (resolution / PI ) * (  theta) * P
	// (1-theta*resolution/(2PI)) * (  phi*resolution/PI) -> -(resolution / (2PI) ) * (  phi) * T / sin(phi) + (resolution / PI ) * (1-theta) * P
	// (  theta*resolution/(2PI)) * (  phi*resolution/PI) ->  (resolution / (2PI) ) * (  phi) * T / sin(phi) + (resolution / PI ) * (  theta) * P

	auto GradientValue = [&]( const Real *signal )
	{
		return
			Point2D< Real >( - (Real)(1.-dj) * thetaScale , - (Real)(1.-di) * phiScale ) * signal[ VertexIndex( i+0 , j+0 , _resolution ) ] +
			Point2D< Real >(   (Real)(1.-dj) * thetaScale , - (Real)(   di) * phiScale ) * signal[ VertexIndex( i+1 , j+0 , _resolution ) ] +
			Point2D< Real >( - (Real)(   dj) * thetaScale ,   (Real)(1.-di) * phiScale ) * signal[ VertexIndex( i+0 , j+1 , _resolution ) ] +
			Point2D< Real >(   (Real)(   dj) * thetaScale ,   (Real)(   di) * phiScale ) * signal[ VertexIndex( i+1 , j+1 , _resolution ) ] ;
	};
	Point2D< Real > g;
	if(  Gradient ) g +=    GradientValue( &_signal[0]                                     );
	if( JGradient ) g += J( GradientValue( &_signal[ VertexNum(_resolution)*(Blocks-1) ] ) );
	return tangents.first * g[0] + tangents.second * g[1];
}

template< typename Real >
template< bool Gradient , bool JGradient >
Point3D< Real > SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::_sample( const SquareMatrix< Real , 2 > &g , int i , int j , Real di , Real dj ) const
{
	Real theta = (Real)( (i+di) * 2 * M_PI / _resolution );
	Real phi = (Real)( (j+dj) * M_PI / _resolution );
	std::pair< Point3D< Real > , Point3D< Real > > tangents = SphereTangents( theta , phi );

	Real sinPhi = (Real)sin( phi );
	Real thetaScale = (Real)( _resolution / (2*M_PI) ) / sinPhi , phiScale = (Real)( _resolution / M_PI );

	auto GradientValue = [&]( const Real *signal )
	{
		return
			Point2D< Real >( - (Real)(1.-dj) * thetaScale , - (Real)(1.-di) * phiScale ) * signal[ VertexIndex( i+0 , j+0 , _resolution ) ] +
			Point2D< Real >(   (Real)(1.-dj) * thetaScale , - (Real)(   di) * phiScale ) * signal[ VertexIndex( i+1 , j+0 , _resolution ) ] +
			Point2D< Real >( - (Real)(   dj) * thetaScale ,   (Real)(1.-di) * phiScale ) * signal[ VertexIndex( i+0 , j+1 , _resolution ) ] +
			Point2D< Real >(   (Real)(   dj) * thetaScale ,   (Real)(   di) * phiScale ) * signal[ VertexIndex( i+1 , j+1 , _resolution ) ] ;
	};
	Point2D< Real > v;
	if(  Gradient ) v +=        GradientValue( &_signal[0]                                     );
	if( JGradient ) v += J( g , GradientValue( &_signal[ VertexNum(_resolution)*(Blocks-1) ] ) );
	v = g.inverse() * v;
	return tangents.first * v[0] + tangents.second * v[1];
}

template< typename Real >
template< bool Gradient , bool JGradient >
Point3D< Real > SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::operator() ( Point3D< Real > p ) const
{
	Point2D< Real > thetaPhi = SphereParameters( p );
	thetaPhi[0] *= (Real)( _resolution / ( 2.*M_PI ) );
	thetaPhi[1] *= (Real)( _resolution / M_PI );
	int theta = (int)floor(thetaPhi[0]) , phi = (int)floor(thetaPhi[1]);
	if( phi==_resolution ) phi--;
	Real dTheta = thetaPhi[0]-theta , dPhi = thetaPhi[1]-phi;
	return _sample( theta , phi , dTheta , dPhi );
}

template< typename Real >
template< bool Gradient , bool JGradient >
Point3D< Real > SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::faceValue( size_t idx ) const
{
	int i , j;
	FactorFaceIndex( idx , i , j , _resolution );
	return _sample( i , j , (Real)0.5 , (Real)0.5 );
}

template< typename Real >
template< bool Gradient , bool JGradient >
Point3D< Real > SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::vertexValue( size_t idx ) const
{
	int i , j;
	FactorVertexIndex( idx , i , j , _resolution );
	Point3D< Real > v;
	const Real *divSignal = NULL , *curlSignal = NULL;
	if( Gradient )   divSignal = &_signal[0];
	if( JGradient ) curlSignal = &_signal[ VertexNum(_resolution)*(Blocks-1) ];
	if( j==0 || j==_resolution )
	{
		Point3D< Real > p = SpherePoint( Point2D< Real >( 0 , j==0 ? 0 : M_PI ) );
		Real thetaScale = (Real)( _resolution / (2*M_PI) ) , phiScale = (Real)( _resolution / M_PI );
#if 0
		Real theta0 = (Real)( (i+0) * 2 * M_PI / _resolution );
		Real theta1 = (Real)( (i+1) * 2 * M_PI / _resolution );
		Real phi = (Real)( j==0 ? 0 : M_PI );
		std::pair< Point3D< Real > , Point3D< Real > > tangents0 = SphereTangents( theta0 , phi );
		std::pair< Point3D< Real > , Point3D< Real > > tangents1 = SphereTangents( theta1 , phi );
		auto GradientValue = [&]( const Real *signal , int i , int j , Real di )
		{
			if( j==0 )
			{
				return
					Point2D< Real >( - thetaScale * phiScale , (Real)(1.-di) * phiScale ) * signal[ VertexIndex( i+0 , ( j==0 ? 1 : (int)_resolution-1 ) , _resolution ) ] +
					Point2D< Real >(   thetaScale * phiScale , (Real)(   di) * phiScale ) * signal[ VertexIndex( i+1 , ( j==0 ? 1 : (int)_resolution-1 ) , _resolution ) ] ;
			}
			else
			{
				return
					Point2D< Real >( - thetaScale * phiScale , - (Real)(1.-di) * phiScale ) * signal[ VertexIndex( i+0 , j+0 , _resolution ) ] +
					Point2D< Real >(   thetaScale * phiScale , - (Real)(   di) * phiScale ) * signal[ VertexIndex( i+1 , j+0 , _resolution ) ] ;
			}
		};
		Point2D< Real > g0 , g1;
		if( j==0 )
		{
			for( int i=0 ; i<_resolution ; i++ )
			{
				if(  Gradient ) g0 +=    GradientValue(  divSignal , i , j , 0 )   , g1 +=    GradientValue(  divSignal , i , j , 1 )  ;
				if( JGradient ) g0 -= J( GradientValue( curlSignal , i , j , 0 ) ) , g1 -= J( GradientValue( curlSignal , i , j , 1 ) );
			}
		}
		else
		{
			for( int i=0 ; i<_resolution ; i++ )
			{
				if(  Gradient ) g0 +=    GradientValue(  divSignal , i , j-1 , 0 )   , g1 +=    GradientValue(  divSignal , i , j-1 , 1 )  ;
				if( JGradient ) g0 -= J( GradientValue( curlSignal , i , j-1 , 0 ) ) , g1 -= J( GradientValue( curlSignal , i , j-1 , 1 ) );
			}
		}
		v += tangents0.first * g0[0] + tangents0.second * g0[1];
		v += tangents1.first * g1[0] + tangents1.second * g1[1];
#else
		auto GradientValue = [&]( const Real *signal , int i , int j , Real di )
		{
			Real theta = (Real)( (i+di) * 2 * M_PI / _resolution );
			Real phi = (Real)( j==0 ? 0 : M_PI );
			std::pair< Point3D< Real > , Point3D< Real > > tangents = SphereTangents( theta , phi );
			if( j==0 )
			{
				return
					( - tangents.first * thetaScale * phiScale + tangents.second * (Real)(1.-di) * phiScale ) * signal[ VertexIndex( i+0 , ( j==0 ? 1 : (int)_resolution-1 ) , _resolution ) ] +
					(   tangents.first * thetaScale * phiScale + tangents.second * (Real)(   di) * phiScale ) * signal[ VertexIndex( i+1 , ( j==0 ? 1 : (int)_resolution-1 ) , _resolution ) ] ;
			}
			else
			{
				return
					( - tangents.first * thetaScale * phiScale - tangents.second * (Real)(1.-di) * phiScale ) * signal[ VertexIndex( i+0 , j+0 , _resolution ) ] +
					(   tangents.first * thetaScale * phiScale - tangents.second * (Real)(   di) * phiScale ) * signal[ VertexIndex( i+1 , j+0 , _resolution ) ] ;
			}
		};
		if( j==0 )
		{
			for( int i=0 ; i<_resolution ; i++ )
			{
				if(  Gradient ) v += GradientValue( divSignal , i , j , 0 ) +  GradientValue( divSignal , i , j , 1 );
				//				if( JGradient ) v -= Point3D< Real >::CrossProduct( p , GradientValue( curlSignal , i , j , 0 ) + GradientValue( curlSignal , i , j , 1 ) );
				if( JGradient ) v += Point3D< Real >::CrossProduct( p , GradientValue( curlSignal , i , j , 0 ) + GradientValue( curlSignal , i , j , 1 ) );
			}
		}
		else
		{
			for( int i=0 ; i<_resolution ; i++ )
			{
				if(  Gradient ) v += GradientValue( divSignal , i , j-1 , 0 ) + GradientValue( divSignal , i , j-1 , 1 );
				//				if( JGradient ) v -= Point3D< Real >::CrossProduct( p , GradientValue( curlSignal , i , j-1 , 0 ) + GradientValue( curlSignal , i , j-1 , 1 ) );
				if( JGradient ) v += Point3D< Real >::CrossProduct( p , GradientValue( curlSignal , i , j-1 , 0 ) + GradientValue( curlSignal , i , j-1 , 1 ) );
			}
		}
#endif

		v /= 2 * (Real)_resolution;
	}
	else v = ( _sample( i , j , 0 , 0 ) + _sample( i-1 , j , 1 , 0 ) + _sample( i , j-1 , 0 , 1 ) + _sample( i-1 , j-1 , 1 , 1 ) ) / 4;
	return v;
}

template< typename Real >
template< bool Gradient , bool JGradient >
template< typename MetricTensor3x3 >
SparseMatrix< Real , int > SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::MassMatrix( size_t resolution , unsigned int mcSamples , const MetricTensor3x3 &metricTensor )
{
	auto ReducedTensor = []( const std::pair< Point3D< Real > , Point3D< Real > > &tangents , const SquareMatrix< Real , 3 > &tensor )
	{
		SquareMatrix< Real , 2 > _tensor;
		_tensor(0,0) = Point3D< Real >::Dot( tangents.first  , tensor * tangents.first  );
		_tensor(1,0) = Point3D< Real >::Dot( tangents.first  , tensor * tangents.second );
		_tensor(0,1) = Point3D< Real >::Dot( tangents.second , tensor * tangents.first  );
		_tensor(1,1) = Point3D< Real >::Dot( tangents.second , tensor * tangents.second );
		return _tensor;
	};
	return _MassMatrix( resolution , mcSamples , [&]( int i , int j ){ return ReducedTensor( FaceTangents( i , j , resolution ) , metricTensor(i,j) ); } );
}

template< typename Real >
template< bool Gradient , bool JGradient >
void SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::_InitSystemMatrix( size_t resolution , SparseMatrix< Real , int > &M )
{
	size_t sz = RegularGrid< Real >::VertexNum( resolution );

	// Allocte the memory for the sparse matrices
	M.resize( sz*Blocks );

	auto InitEntry = [&]( size_t idx , size_t _idx ){ M[idx][ _RowOffset( idx , _idx , resolution ) ] = MatrixEntry< Real , int >( (int)_idx , 0 ); };
#pragma omp parallel for
	for( int j=0 ; j<=resolution ; j++ )
	{
		if( j==0 )
		{
			for( int k=0 ; k<Blocks ; k++ )
			{
				size_t idx = VertexIndex( 0 , j , resolution ) + sz*k , _idx;
				M.SetRowSize( (int)idx , Blocks*(resolution+1) );
				_idx = VertexIndex( 0 , j , resolution );
				for( int kk=0 ; kk<Blocks ; kk++ ) InitEntry( idx , _idx+kk*sz );
				for( int i=0 ; i<resolution ; i++ )
				{
					_idx = VertexIndex( i , j+1 , resolution );
					for( int kk=0 ; kk<Blocks ; kk++ ) InitEntry( idx , _idx+kk*sz );
				}
			}
		}
		else if( j==resolution )
		{
			for( int k=0 ; k<Blocks ; k++ )
			{
				size_t idx = VertexIndex( 0 , j , resolution ) + sz*k , _idx;
				M.SetRowSize( (int)idx , Blocks*(resolution+1) );
				for( int i=0 ; i<resolution ; i++ )
				{
					_idx = VertexIndex( i , j-1 , resolution );
					for( int kk=0 ; kk<Blocks ; kk++ ) InitEntry( idx , _idx+kk*sz );
				}
				_idx = VertexIndex( 0 , j , resolution );
				for( int kk=0 ; kk<Blocks ; kk++ ) InitEntry( idx , _idx+kk*sz );
			}
		}
		else if( j==1 )
		{
			for( int i=0 ; i<resolution ; i++ ) for( int k=0 ; k<Blocks ; k++ )
			{
				size_t idx = VertexIndex( i , j , resolution ) + sz*k , _idx;
				M.SetRowSize( (int)idx , Blocks*7 );
				_idx = VertexIndex( 0 , j-1 , resolution );
				for( int kk=0 ; kk<Blocks ; kk++ ) InitEntry( idx , _idx+kk*sz );
				for( int jj=0 ; jj<=1 ; jj++ ) for( int ii=-1 ; ii<=1 ; ii++ )
				{
					_idx = VertexIndex( i+ii , j+jj , resolution );
					for( int kk=0 ; kk<Blocks ; kk++ ) InitEntry( idx , _idx+kk*sz );
				}
			}
		}
		else if( j==resolution-1 )
		{
			for( int i=0 ; i<resolution ; i++ ) for( int k=0 ; k<Blocks ; k++ )
			{
				size_t idx = VertexIndex( i , j , resolution ) + sz*k , _idx;
				M.SetRowSize( (int)idx , Blocks*7 );
				for( int jj=-1 ; jj<=0 ; jj++ ) for( int ii=-1 ; ii<=1 ; ii++ )
				{
					_idx = VertexIndex( i+ii , j+jj , resolution );
					for( int kk=0 ; kk<Blocks ; kk++ ) InitEntry( idx , _idx+kk*sz );
				}
				_idx = VertexIndex( 0 , j+1 , resolution );
				for( int kk=0 ; kk<Blocks ; kk++ ) InitEntry( idx , _idx+kk*sz );
			}
		}
		else
		{
			for( int i=0 ; i<resolution ; i++ ) for( int k=0 ; k<Blocks ; k++ )
			{
				size_t idx = VertexIndex( i , j , resolution ) + sz*k , _idx;
				M.SetRowSize( (int)idx , Blocks*9 );
				for( int jj=-1 ; jj<=1 ; jj++ ) for( int ii=-1 ; ii<=1 ; ii++ )
				{
					_idx = VertexIndex( i+ii , j+jj , resolution );
					for( int kk=0 ; kk<Blocks ; kk++ ) InitEntry( idx , _idx+kk*sz );
				}
			}
		}
	}
}

template< typename Real >
template< bool Gradient , bool JGradient >
size_t SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::_RowOffset( size_t idx1 , size_t idx2 , size_t resolution )
{
	size_t sz = RegularGrid< Real >::VertexNum( resolution );

	auto ThetaNeighbors = [] ( int i , int j , size_t resolution )
	{
		if( ( j>=i-1 && j<=i+1 ) || ( ( j>=i+resolution-1 && j<=i+resolution+1 ) ) || ( j+resolution>=i-1 && j+resolution<=i+1 ) ) return true;
		return false;
	};

	int i1 , i2 , j1 , j2 , off=-1;
	FactorVertexIndex( idx1%sz , i1 , j1 , resolution );
	FactorVertexIndex( idx2%sz , i2 , j2 , resolution );
	if( j1==0 )
	{
		if     ( j2==0 ) off = 0;
		else if( j2==1 ) off = i2+1;
	}
	else if( j1==resolution )
	{
		if     ( j2==resolution-1 ) off = i2;
		else if( j2==resolution   ) off = (int)resolution;
	}
	else if( j1==1 )
	{
		if     ( j2==0 ) off = 0;
		else if( j2==1 ){ if( ThetaNeighbors( i1 , i2 , resolution ) ) off = (int)( ( (i2-i1+1+resolution )%resolution ) + 1 ); }
		else if( j2==2 ){ if( ThetaNeighbors( i1 , i2 , resolution ) ) off = (int)( ( (i2-i1+1+resolution )%resolution ) + 4 ); }
	}
	else if( j1==resolution-1 )
	{
		if     ( j2==resolution-2 ){ if( ThetaNeighbors( i1 , i2 , resolution ) ) off = (int)( ( (i2-i1+1+resolution )%resolution ) + 0 ); }
		else if( j2==resolution-1 ){ if( ThetaNeighbors( i1 , i2 , resolution ) ) off = (int)( ( (i2-i1+1+resolution )%resolution ) + 3 ); }
		else if( j2==resolution   ) off = 6;
	}
	else{ if( j2>=j1-1 && j2<=j1+1 && ThetaNeighbors( i1 , i2 , resolution ) ) off = (j2-j1+1)*3 + (int)( ( (i2-i1+1)+resolution )%resolution ); }
	if( off>=0 )
	{
		if( Blocks==2 ) return idx2<sz ? 2*off : 2*off+1;
		else            return off;
	}
	fprintf( stderr , "[ERROR] Vertices are not neighbors: [%d](%d %d) [%d](%d %d) @ %d\n" , (int)idx1 , i1 , j1 , (int)idx2 , i2 , j2 , (int)resolution ) , exit( 0 );
	return -1;
}

template< typename Real >
template< bool Gradient , bool JGradient >
template< typename MetricTensor2x2 >
SparseMatrix< Real , int > SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::_MassMatrix( size_t resolution , unsigned int mcSamples , const MetricTensor2x2 &metricTensor )
{
	size_t sz = RegularGrid< Real >::VertexNum( resolution );

	SparseMatrix< Real , int > M;
	_InitSystemMatrix( resolution , M );

	// Set the matrix coefficients by iterating over the faces
#pragma omp parallel for
	for( int j=0 ; j<resolution ; j++ )
	{
		// Compute the integrals of monomials in \phi
		MonomialIntegralTable< Real , 2 > mit = _GetMonomialIntegralTable< 2 >( j , resolution , mcSamples );

		BilinearFunctionType< Real > F[4];
		BilinearGradientType< Real > G[4];
		RotatedBilinearGradientType< Real > JG[4];

		_BilinearFaceBasis( j , resolution , F );
		for( int ii=0 ; ii<4 ; ii++ )
		{
			G[ii] = F[ii].gradient();
			JG[ii] = G[ii].rotate90();
		}
		SinYPolynomial< Real , 0 , 0 , 1 > sine(1);

		auto PartialDot = [&]( unsigned int ii , unsigned int jj )
		{
			Real dTheta = (Real)( 2. * M_PI / resolution );
			SquareMatrix< Real , 2 > dot;
			if( ii<4 )
			{
				if( jj<4 )
				{
					dot(0,0) = ( G[ii].x * G[jj].x * sine ).integrate( 0 , dTheta , mit );
					dot(1,0) = ( G[ii].y * G[jj].x * sine ).integrate( 0 , dTheta , mit );
					dot(0,1) = ( G[ii].x * G[jj].y * sine ).integrate( 0 , dTheta , mit );
					dot(1,1) = ( G[ii].y * G[jj].y * sine ).integrate( 0 , dTheta , mit );
				}
				else
				{
					dot(0,0) = ( G[ii].x * JG[jj-4].x * sine ).integrate( 0 , dTheta , mit );
					dot(1,0) = ( G[ii].y * JG[jj-4].x * sine ).integrate( 0 , dTheta , mit );
					dot(0,1) = ( G[ii].x * JG[jj-4].y * sine ).integrate( 0 , dTheta , mit );
					dot(1,1) = ( G[ii].y * JG[jj-4].y * sine ).integrate( 0 , dTheta , mit );
				}
			}
			else
			{
				if( jj<4 )
				{
					dot(0,0) = ( JG[ii-4].x * G[jj].x * sine ).integrate( 0 , dTheta , mit );
					dot(1,0) = ( JG[ii-4].y * G[jj].x * sine ).integrate( 0 , dTheta , mit );
					dot(0,1) = ( JG[ii-4].x * G[jj].y * sine ).integrate( 0 , dTheta , mit );
					dot(1,1) = ( JG[ii-4].y * G[jj].y * sine ).integrate( 0 , dTheta , mit );
				}
				else
				{
					dot(0,0) = ( JG[ii-4].x * JG[jj-4].x * sine ).integrate( 0 , dTheta , mit );
					dot(1,0) = ( JG[ii-4].y * JG[jj-4].x * sine ).integrate( 0 , dTheta , mit );
					dot(0,1) = ( JG[ii-4].x * JG[jj-4].y * sine ).integrate( 0 , dTheta , mit );
					dot(1,1) = ( JG[ii-4].y * JG[jj-4].y * sine ).integrate( 0 , dTheta , mit );
				}
			}

			return dot;
		};

		SquareMatrix< Real , 2 > s[8][8];
#if 0
		SquareMatrix< Real , 2 > _s[4][4];
		_SetFaceStiffnesses( j , resolution , mit , _s );
		for( int ii=0 ; ii<4; ii++ ) for( int jj=0 ; jj<4 ; jj++ )
		{
			s[ii+0][jj+0] =                   _s[ii][jj]      ;
			s[ii+4][jj+0] =                   _s[ii][jj] * J();
			s[ii+0][jj+4] = J().transpose() * _s[ii][jj]      ;
			s[ii+4][jj+4] = J().transpose() * _s[ii][jj] * J();
		}
#else
		for( int ii=0 ; ii<8 ; ii++ ) for( int jj=0 ; jj<8 ; jj++ ) s[ii][jj] = PartialDot(ii,jj);
#endif
#if 0
SquareMatrix< Real , 2 > _s[4][4];
_SetFaceStiffnesses( j , resolution , mit , _s );
#pragma omp critical
{
	printf( "%d] " , j );
	{
		double foo = 0;
		for( int i=0 ; i<4 ; i++ ) for( int j=0 ; j<4 ; j++ )
		{
			for( int ii=0 ; ii<2 ; ii++ ) for( int jj=0 ; jj<2 ; jj++ ) foo += ( s[i][j](ii,jj) - _s[i][j](ii,jj) ) * ( s[i][j](ii,jj) - _s[i][j](ii,jj) );
		}
		printf( " %g" , foo );
	}
	{
		double foo = 0;
		for( int i=0 ; i<4 ; i++ ) for( int j=0 ; j<4 ; j++ )
		{
			SquareMatrix< Real , 2 > __s = _s[i][j] * J().transpose();
			for( int ii=0 ; ii<2 ; ii++ ) for( int jj=0 ; jj<2 ; jj++ ) foo += ( s[i+4][j](ii,jj) - __s(ii,jj) ) * ( s[i+4][j](ii,jj) - __s(ii,jj) );
		}
		printf( " %g" , foo );
	}
	{
		double foo = 0;
		for( int i=0 ; i<4 ; i++ ) for( int j=0 ; j<4 ; j++ )
		{
			SquareMatrix< Real , 2 > __s = J() * _s[i][j];
			for( int ii=0 ; ii<2 ; ii++ ) for( int jj=0 ; jj<2 ; jj++ ) foo += ( s[i][j+4](ii,jj) - __s(ii,jj) ) * ( s[i][j+4](ii,jj) - __s(ii,jj) );
		}
		printf( " %g" , foo );
	}
	{
		double foo = 0;
		for( int i=0 ; i<4 ; i++ ) for( int j=0 ; j<4 ; j++ )
		{
			SquareMatrix< Real , 2 > __s = J() * _s[i][j] * J().transpose();
			for( int ii=0 ; ii<2 ; ii++ ) for( int jj=0 ; jj<2 ; jj++ ) foo += ( s[i+4][j+4](ii,jj) - __s(ii,jj) ) * ( s[i+4][j+4](ii,jj) - __s(ii,jj) );
		}
		printf( " %g" , foo );
	}
	printf( "\n" );
}
#endif


		for( int i=0 ; i<resolution ; i++ )
		{
			// Get the vertices on the face
			FaceIndices face = Face( i , j , resolution );

			SquareMatrix< Real , 2 > tensor = metricTensor(i,j);

			// Accumulate into the matrix
			for( int ii=0 ; ii<face.size ; ii++ ) for( int jj=0 ; jj<face.size ; jj++ ) for( int iii=0 ; iii<Blocks ; iii++ ) for( int jjj=0 ; jjj<Blocks ; jjj++ )
			{
				size_t idx1 = face[ii] + sz*iii , idx2 = face[jj] + sz*jjj;
				size_t offset = _RowOffset( idx1 , idx2 , resolution );
				Real value = SquareMatrix< Real , 2 >::Dot( tensor , Gradient ? s[ii+4*iii][jj+4*jjj] : s[ii+4][jj+4] );
				AddAtomic( M[ idx1 ][ offset ].Value , value );
			}
		}
	}
	return M;
}

template< typename Real >
template< bool Gradient , bool JGradient >
template< bool _Gradient , bool _JGradient >
std::vector< Real > SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::WeightedMass( unsigned int Dim , const PrimalSignal< Real > *weights , const PrimalVectorField< _Gradient , _JGradient > *vf , unsigned int mcSamples )
{
	size_t resolution = weights[0].resolution();
	for( unsigned int d=0 ; d<Dim ; d++ ) if( weights[d].resolution()!=resolution || vf[d].resolution()!=resolution ) fprintf( stderr , "[ERROR] Resolutions don't match\n" ) , exit( 0 );
	size_t sz = RegularGrid< Real >::VertexNum( resolution );
	std::vector< Real > mass( Blocks*sz , 0 );

	// Set the matrix coefficients by iterating over the faces
#pragma omp parallel for
	for( int j=0 ; j<resolution ; j++ )
	{
		// Compute the integrals of monomials in \phi
		MonomialIntegralTable< Real , 3 > mit = _GetMonomialIntegralTable< 3 >( j , resolution , mcSamples );

		BilinearFunctionType< Real > F[4];
		BilinearGradientType< Real > G[4];
		RotatedBilinearGradientType< Real > JG[4];

		_BilinearFaceBasis( j , resolution , F );
		for( int ii=0 ; ii<4 ; ii++ )
		{
			G[ii] = F[ii].gradient();
			JG[ii] = G[ii].rotate90();
		}
		SinYPolynomial< Real , 0 , 0 , 1 > sine(1);

		auto Dot = [&]( unsigned int ii , unsigned int jj , unsigned int kk )
		{
			if( jj<4 )
			{
				if( kk<4 )
				{
					return
						( F[ii] * G[jj].x * G[kk].x * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit ) +
						( F[ii] * G[jj].y * G[kk].y * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit ) ;
				}
				else
				{
					return
						( F[ii] * G[jj].x * JG[kk-4].x * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit ) +
						( F[ii] * G[jj].y * JG[kk-4].y * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit ) ;
				}
			}
			else
			{
				if( kk<4 )
				{
					return
						( F[ii] * JG[jj-4].x * G[kk].x * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit ) +
						( F[ii] * JG[jj-4].y * G[kk].y * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit ) ;
				}
				else
				{
					return
						( F[ii] * JG[jj-4].x * JG[kk-4].x * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit ) +
						( F[ii] * JG[jj-4].y * JG[kk-4].y * sine ).integrate( 0 , (Real)( 2.*M_PI / resolution ) , mit ) ;
				}
			}
		};

		Real integrals[4][8][8];
		for( int ii=0 ; ii<4 ; ii++ ) for( int jj=0 ; jj<8 ; jj++ ) for( int kk=0 ; kk<8 ; kk++ ) integrals[ii][jj][kk] = Dot( ii , jj , kk );

		for( int i=0 ; i<resolution ; i++ )
		{
			// Get the vertices on the face
			FaceIndices face = Face( i , j , resolution );

			// Accumulate into the constraint vector
			static const unsigned int _Blocks = PrimalVectorField< _Gradient , _JGradient >::Blocks;
			for( int ii=0 ; ii<face.size ; ii++ ) for( int jj=0 ; jj<face.size ; jj++ ) for( int kk=0 ; kk<face.size ; kk++ ) for( int iii=0 ; iii<Blocks ; iii++ ) for( int jjj=0 ; jjj<_Blocks ; jjj++ )
			{
				size_t idx1 = face[ii] + sz*iii , idx2 = face[jj] + sz*jjj , idx3 = face[kk];
				Real value=0; 
				for( unsigned int d=0 ; d<Dim ; d++ ) value += vf[d][idx2] * weights[d][idx3];
				if     ( Gradient && _Gradient ) value *= integrals[kk][jj+4*jjj][ii+4*iii];
				else if( Gradient              ) value *= integrals[kk][jj+4][ii+4*iii];
				else if(             _Gradient ) value *= integrals[kk][jj+4*jjj][ii+4];
				else                             value *= integrals[kk][jj+4][ii+4];
				AddAtomic( mass[ idx1 ] , value );
			}
		}
	}
	return mass;
}

template< typename Real >
template< bool Gradient , bool JGradient >
std::vector< Real > SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::Mass( const DualSignal< Point3D< Real > > &V , unsigned int mcSamples )
{
	size_t resolution = V.resolution();
	size_t sz = RegularGrid< Real >::VertexNum( resolution );
	std::vector< Real > mass( Blocks*sz , 0 );

	// Set the vector coefficients by iterating over the faces
#pragma omp parallel for
	for( int j=0 ; j<resolution ; j++ )
	{
		// Compute the integrals of monomials in \phi
		MonomialIntegralTable< Real , 3 > mit = _GetMonomialIntegralTable< 3 >( j , resolution , mcSamples );

		BilinearFunctionType< Real > F[4];
		BilinearGradientType< Real > G[4];
		RotatedBilinearGradientType< Real > JG[4];

		_BilinearFaceBasis( j , resolution , F );
		for( int ii=0 ; ii<4 ; ii++ )
		{
			G[ii] = F[ii].gradient();
			JG[ii] = G[ii].rotate90();
		}
		SinYPolynomial< Real , 0 , 0 , 1 > sine(1);

		auto Dot = [&]( unsigned int ii )
		{
			Real dTheta = (Real)( 2. * M_PI / resolution );
			if( ii<4 ) return Point2D< Real >( (  G[ii  ].x * sine ).integrate( 0 , dTheta , mit ) , (  G[ii  ].y * sine ).integrate( 0 , dTheta , mit ) );
			else       return Point2D< Real >( ( JG[ii-4].x * sine ).integrate( 0 , dTheta , mit ) , ( JG[ii-4].y * sine ).integrate( 0 , dTheta , mit ) );
		};

		Point2D< Real > integrals[8];
		for( int ii=0 ; ii<8 ; ii++ ) integrals[ii] = Dot( ii );

		for( int i=0 ; i<resolution ; i++ )
		{
			// Get the vertices on the face
			FaceIndices face = Face( i , j , resolution );
			std::pair< Point3D< Real > , Point3D< Real > > tangents = FaceTangents( i , j , resolution );

			Point3D< Real > v = V(i,j);
			Point2D< Real > _v( Point3D< Real >::Dot( v , tangents.first ) , Point3D< Real >::Dot( v , tangents.second ) );

			// Accumulate into the constraint vector
			for( int ii=0 ; ii<face.size ; ii++ ) for( int iii=0 ; iii<Blocks ; iii++ )
			{
				size_t idx = face[ii] + sz*iii;
				AddAtomic( mass[ idx ] , Point2D< Real >::Dot( _v , Gradient ? integrals[ii+4*iii] : integrals[ii+4] ) );
			}
		}
	}
	return mass;
}

template< typename Real >
template< bool Gradient , bool JGradient >
template< typename MetricTensor3x3 > 
std::vector< Real > SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::mass( unsigned int mcSamples , const MetricTensor3x3 &metricTensor ) const
{
	auto ReducedTensor = []( const std::pair< Point3D< Real > , Point3D< Real > > &tangents , const SquareMatrix< Real , 3 > &tensor )
	{
		SquareMatrix< Real , 2 > _tensor;
		_tensor(0,0) = Point3D< Real >::Dot( tangents.first  , tensor * tangents.first  );
		_tensor(1,0) = Point3D< Real >::Dot( tangents.first  , tensor * tangents.second );
		_tensor(0,1) = Point3D< Real >::Dot( tangents.second , tensor * tangents.first  );
		_tensor(1,1) = Point3D< Real >::Dot( tangents.second , tensor * tangents.second );
		return _tensor;
	};
	return _mass( mcSamples , [&]( int i , int j ){ return ReducedTensor( FaceTangents( i , j , _resolution ) , metricTensor(i,j) ); } );
}

template< typename Real >
template< bool Gradient , bool JGradient >
template< typename MetricTensor2x2 > 
std::vector< Real > SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::_mass( unsigned int mcSamples , const MetricTensor2x2 &metricTensor ) const
{
	size_t sz = RegularGrid< Real >::VertexNum( _resolution );
	std::vector< Real > mass( Blocks*sz , 0 );

	// Set the coefficients coefficients by iterating over the faces
#pragma omp parallel for
	for( int j=0 ; j<_resolution ; j++ )
	{
		// Compute the integrals of monomials in \phi
		MonomialIntegralTable< Real , 2 > mit = _GetMonomialIntegralTable< 2 >( j , _resolution , mcSamples );

		BilinearFunctionType< Real > F[4];
		BilinearGradientType< Real > G[4];
		RotatedBilinearGradientType< Real > JG[4];

		_BilinearFaceBasis( j , _resolution , F );
		for( int ii=0 ; ii<4 ; ii++ )
		{
			G[ii] = F[ii].gradient();
			JG[ii] = G[ii].rotate90();
		}
		SinYPolynomial< Real , 0 , 0 , 1 > sine(1);

		auto Dot = [&]( unsigned int ii , unsigned int jj )
		{
			SquareMatrix< Real , 2 > dot;
			if( ii<4 )
			{
				if( jj<4 )
				{
					dot(0,0) = ( G[ii].x * G[jj].x * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
					dot(1,0) = ( G[ii].y * G[jj].x * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
					dot(0,1) = ( G[ii].x * G[jj].y * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
					dot(1,1) = ( G[ii].y * G[jj].y * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
				}
				else
				{
					dot(0,0) = ( G[ii].x * JG[jj-4].x * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
					dot(1,0) = ( G[ii].y * JG[jj-4].x * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
					dot(0,1) = ( G[ii].x * JG[jj-4].y * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
					dot(1,1) = ( G[ii].y * JG[jj-4].y * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
				}
			}
			else
			{
				if( jj<4 )
				{
					dot(0,0) = ( JG[ii-4].x * G[jj].x * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
					dot(1,0) = ( JG[ii-4].y * G[jj].x * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
					dot(0,1) = ( JG[ii-4].x * G[jj].y * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
					dot(1,1) = ( JG[ii-4].y * G[jj].y * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
				}
				else
				{
					dot(0,0) = ( JG[ii-4].x * JG[jj-4].x * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
					dot(1,0) = ( JG[ii-4].y * JG[jj-4].x * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
					dot(0,1) = ( JG[ii-4].x * JG[jj-4].y * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
					dot(1,1) = ( JG[ii-4].y * JG[jj-4].y * sine ).integrate( 0 , (Real)( 2.*M_PI / _resolution ) , mit );
				}
			}

			return dot;
		};

		SquareMatrix< Real , 2 > s[8][8];
		for( int ii=0 ; ii<8 ; ii++ ) for( int jj=0 ; jj<8 ; jj++ ) s[ii][jj] = Dot( ii , jj );

		for( int i=0 ; i<_resolution ; i++ )
		{
			// Get the vertices on the face
			FaceIndices face = Face( i , j , _resolution );

			SquareMatrix< Real , 2 > tensor = metricTensor( i , j );

			// Accumulate into the matrix
			for( int ii=0 ; ii<face.size ; ii++ ) for( int jj=0 ; jj<face.size ; jj++ ) for( int iii=0 ; iii<Blocks ; iii++ ) for( int jjj=0 ; jjj<Blocks ; jjj++ )
			{
				size_t idx1 = face[ii] + sz*iii , idx2 = face[jj] + sz*jjj;
				Real value = SquareMatrix< Real , 2 >::Dot( tensor , Gradient ? s[ii+4*iii][jj+4*jjj] : s[ii+4][jj+4] ) * _signal[idx2];
				AddAtomic( mass[ idx1 ] , value );
			}
		}
	}
	return mass;
}

template< typename Real >
template< bool Gradient , bool JGradient >
SparseMatrix< Real , int > SphericalGeometry::RegularGrid< Real >::PrimalVectorField< Gradient , JGradient >::RestrictionMatrix( size_t coarseResolution )
{
	SparseMatrix< Real , int > R;
	if( Gradient!=JGradient ) R = PrimalSignal< Real >::RestrictionMatrix( coarseResolution );
	else if( Gradient && !Gradient )
	{
		size_t lowOff = VertexNum( coarseResolution ) , highOff = VertexNum( 2*coarseResolution );
		SparseMatrix< Real , int > _R = PrimalSignal< Real >::RestrictionMatrix( coarseResolution );
		R.resize( 2*_R.rows );
#pragma omp parallel for
		for( int i=0 ; i<_R.rows ; i++ )
		{
			R.SetRowSize( i         , _R.rowSizes[i] );
			R.SetRowSize( i+_R.rows , _R.rowSizes[i] );
			for( int j=0 ; j<_R.rowSizes[i] ; j++ )
			{
				R[i       ][j] = MatrixEntry< Real , int >( _R[i][j].N              , _R[i][j].Value );
				R[i+lowOff][j] = MatrixEntry< Real , int >( _R[i][j].N+(int)highOff , _R[i][j].Value );
			}
		}
	}
	return R;
}
