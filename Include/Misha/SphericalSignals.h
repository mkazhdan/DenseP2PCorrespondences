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
// To do:
// 1. Functors for tensors
// 2. Dual vector fields
// 3. Point-wise evaluation of primal signals
#ifndef SPHERICAL_SIGNALS_INCLUDED
#define SPHERICAL_SIGNALS_INCLUDED

#include "Misha/SphericalGeometry.h"

namespace SphericalGeometry
{
	// This structure stores integrals of the form x^d * sin(x) and x^d / sin(x) over a given domain
	// Things we want to integrate:
	//	1. Functions:
	//		I(f,g) = \int \int f(\theta,\phi) * g(\theta,\phi) * sin(\phi) d\theta d\phi
	//	2. Gradients:
	//		I(f,g) = \int \int < \nabla f(\theta,\phi) , \nabla g(\theta,\phi) > * sin(\phi) d\theta d\phi
	//		       = \int \int < ( df/d\theta * 1/sin(\phi) , df/d\phi ) , ( dg/d\theta * 1/sin(\phi) , dg/d\phi ) > * sin(\phi) d\theta d\phi
	//		       = \int \int ( df/d\theta * dg/d\theta * 1/sin^2(phi) + df/d\phi * dg/d\phi ) * sin(phi) d\theta d\phi
	//		       = \int \int   df/d\theta * dg/d\theta * 1/sin(phi) + df/d\phi * dg/d\phi * sin(phi) d\theta d\phi
	//	3. Gradients relative to a metric H:
	//		I_H(f,g) = \int \int < \nabla f(\theta,\phi) , H * \nabla g(\theta,\phi) > * sin(\phi) d\theta d\phi
	//		         = \int \int < ( df/d\theta * 1/sin(\phi) , df/d\phi ) , H * ( dg/d\theta * 1/sin(\phi) , dg/d\phi ) > * sin(\phi) d\theta d\phi
	//		         = \int \int	( 
	//								H_00 * df/d\theta * dg/d\theta * 1/sin(\phi) * 1/sin(\phi) +
	//								H_10 * df/d\theta * dg/d\phi * 1/sin(\phi) +
	//								H_01 * df/d\phi * dg/d\theta * 1/sin(\phi) +
	//								H_11 * df/d\phi * dg/d\phi
	//								) * sin(\phi) d\theta d\phi
	//		         = \int \int	( 
	//								H_00 * df/d\theta * dg/d\theta * 1/sin(\phi) +
	//								H_10 * df/d\theta * dg/d\phi +
	//								H_01 * df/d\phi * dg/d\theta +
	//								H_11 * df/d\phi * dg/d\phi * sin(\phi)
	//								) d\theta d\phi

	template< typename Real , unsigned int D >
	struct MonomialIntegralTable
	{
		static void SetIntegralTable( Real start , Real end , Real table[D+1] );
		void init( Real startPhi , Real endPhi , unsigned int mcSamples );
		template< int Exp > typename std::enable_if< Exp==-1 , Real >::type table( unsigned int idx ) const{ return _cscTable[idx]; }
		template< int Exp > typename std::enable_if< Exp== 0 , Real >::type table( unsigned int idx ) const{ return ____table[idx]; }
		template< int Exp > typename std::enable_if< Exp== 1 , Real >::type table( unsigned int idx ) const{ return _sinTable[idx]; }
	protected:
		Real ____table[D+1] , _sinTable[D+1] , _cscTable[D+1];
	};
	template< typename Real , unsigned int D1 , unsigned int D2 , int SinYExponent > struct SinYPolynomial;
	template< typename Real , unsigned int D1 , unsigned int D2 , int SinYExponent , unsigned int _D1 , unsigned int _D2 , int _SinYExponent > struct SinYPolynomialVector;

	template< typename Real , unsigned int D1 , unsigned int D2 >
	struct Polynomial
	{
		Polynomial( void ){ memset( _coefficients , 0 , sizeof(_coefficients) ); }
		Polynomial( Real c ){ memset( _coefficients , 0 , sizeof(_coefficients) ) ; _coefficients[0][0] = c; }
		Real *operator[]( unsigned int d ){ return _coefficients[d]; }
		const Real *operator[]( unsigned int d ) const { return _coefficients[d]; }

		Polynomial operator - ( void ) const;

		Polynomial  operator +  ( const Polynomial< Real , D1 , D2 > &p ) const;
		Polynomial &operator += ( const Polynomial< Real , D1 , D2 > &p );
		Polynomial  operator *  ( Real s ) const;
		Polynomial &operator *= ( Real s );

		Polynomial  operator -  ( const Polynomial< Real , D1 , D2 > &p ) const { return *this + (-p); }
		Polynomial &operator -= ( const Polynomial< Real , D1 , D2 > &p ){ return *this += -p; }
		Polynomial &operator /= ( Real s ){ return *this *= ( (Real)1./s ); }
		Polynomial  operator /  ( Real s ) const { return *this * ( (Real)1./s ); }

		template< unsigned int _D1 , unsigned int _D2 >
		Polynomial< Real , D1+_D1 , D2+_D2 > operator * ( const Polynomial< Real , _D1 , _D2 > &p ) const;

		Real operator()( Real x , Real y ) const;
		Polynomial< Real , D1-1 , D2 > dx( void ) const;
		Polynomial< Real , D1 , D2-1 > dy( void ) const;

		Real integrate( Real thetaStart , Real thetaEnd , Real phiStart , Real phiEnd ) const;
		Real integrateSine( Real thetaStart , Real thetaEnd , Real phiStart , Real phiEnd , unsigned int mcSamples ) const;
		Real integrateCosecant( Real thetaStart , Real thetaEnd , Real phiStart , Real phiEnd , unsigned int mcSamples ) const;
		template< unsigned int D > Real integrateSine( Real thetaStart , Real thetaEnd , const MonomialIntegralTable< Real , D > &mit ) const;
		template< unsigned int D > Real integrateCosecant( Real thetaStart , Real thetaEnd , const MonomialIntegralTable< Real , D > &mit ) const;
	protected:
		Real _coefficients[D1+1][D2+1];
	};
	template< typename Real , unsigned int D1 , unsigned int D2 >
	Polynomial< Real , D1 , D2 > operator * ( Real s , const Polynomial< Real , D1 , D2 > &p ){ return p*s; }
	template< typename Real , unsigned int D1 , unsigned int D2 >
	std::ostream & operator << ( std::ostream &os , const Polynomial< Real , D1 , D2 > &p );

	template< typename Real > using BilinearFunctionType = SinYPolynomial< Real , 1 , 1 , 0 >;
	template< typename Real > using BilinearGradientType = SinYPolynomialVector< Real , 0 , 1 , -1 , 1 , 0 , 0 >;
	template< typename Real > using RotatedBilinearGradientType = typename BilinearGradientType< Real >::RotatedType;

	template< typename Real , unsigned int D1 , unsigned int D2 , int SinYExponent >
	struct SinYPolynomial
	{
		template< typename _Real , unsigned int _D1 , unsigned int _D2 , int _SinYExponent > friend struct SinYPolynomial;

		SinYPolynomial( void ){}
		SinYPolynomial( Real c ) : _poly( c ){}
		SinYPolynomial( const Polynomial< Real , D1 , D2 > &poly ) : _poly( poly ){}
		template< int _SinYExponent > SinYPolynomial( const SinYPolynomial< Real , D1 , D2 , _SinYExponent > &p ) : SinYPolynomial( p._poly ){}
		Real *operator[]( unsigned int d ){ return _poly[d]; }
		const Real *operator[]( unsigned int d ) const { return _poly[d]; }

		SinYPolynomial operator - ( void ) const { return SinYPolynomial( -_poly ); }

		SinYPolynomial  operator +  ( const Polynomial< Real , D1 , D2 > &p ) const { return SinYPolynomial( _poly+p ); }
		SinYPolynomial &operator += ( const Polynomial< Real , D1 , D2 > &p ){ _poly += p ; return *this; }
		SinYPolynomial  operator *  ( Real s ) const { return SinYPolynomial( _poly*s ); }
		SinYPolynomial &operator *= ( Real s ){ _poly *= s ; return *this; }

		SinYPolynomial  operator -  ( const Polynomial< Real , D1 , D2 > &p ) const { return *this + (-p); }
		SinYPolynomial &operator -= ( const Polynomial< Real , D1 , D2 > &p ){ return *this += -p; }
		SinYPolynomial &operator /= ( Real s ){ return *this *= ( (Real)1./s ); }
		SinYPolynomial  operator /  ( Real s ) const { return *this * ( (Real)1./s ); }

		template< unsigned int _D1 , unsigned int _D2 , int _SinYExponent >
		SinYPolynomial< Real , D1+_D1 , D2+_D2 , SinYExponent+_SinYExponent > operator * ( const SinYPolynomial< Real , _D1 , _D2 , _SinYExponent > &p ) const { return SinYPolynomial< Real , D1+_D1 , D2+_D2 , SinYExponent+_SinYExponent >( _poly * p._poly ); }

		Real operator()( Real x , Real y ) const { return (Real)( _poly(x,y) * pow( sin(y) , SinYExponent ) ); }

		Real integrate( Real xStart , Real xEnd , Real yStart , Real yEnd , unsigned int mcSamples ) const;
		template< unsigned int D > Real integrate( Real xStart , Real xEnd , const MonomialIntegralTable< Real , D > &mit ) const;

		typedef SinYPolynomialVector< Real , D1-1 , D2 , -1 , D1 , D2-1 , 0 > GradientType;
		template< int _SinYExponent=SinYExponent >
		typename std::enable_if< _SinYExponent==0 , GradientType >::type gradient( void ){ return GradientType( _poly.dx() , _poly.dy() ); }

		Polynomial< Real , D1 , D2 > &operator()( void ){ return _poly; }
		const Polynomial< Real , D1 , D2 > &operator()( void ) const { return _poly; }
	protected:
		Polynomial< Real , D1 , D2 > _poly;
	};
	template< typename Real , unsigned int D1 , unsigned int D2 , int SinYExponent >
	std::ostream & operator << ( std::ostream &os , const SinYPolynomial< Real , D1 , D2 , SinYExponent > &p );

	template< typename Real , unsigned int D1 , unsigned int D2 , int SinYExponent , unsigned int _D1 , unsigned int _D2 , int _SinYExponent >
	struct SinYPolynomialVector
	{
		typedef SinYPolynomialVector< Real , _D1 , _D2 , _SinYExponent , D1 , D2 , SinYExponent > RotatedType;
		typedef SinYPolynomial< Real ,  D1 ,  D2 ,  SinYExponent > XType;
		typedef SinYPolynomial< Real , _D1 , _D2 , _SinYExponent > YType;
		XType x;
		YType y;

		SinYPolynomialVector( void ){}
		SinYPolynomialVector( const XType &x , const YType &y ){ this->x = x , this->y = y; }
		RotatedType rotate90( void ) const { return RotatedType( -y , x ); }
	};

	template< typename Real >
	struct RegularGrid
	{
		struct FaceIndices
		{
			size_t idx[4];
			size_t size;
			size_t &operator[] ( size_t i ){ return idx[i]; }
			const size_t &operator[] ( size_t i )const { return idx[i]; }
		};

		// \Phi( theta , phi ) = ( sin(phi) * cos(theta)) ) , cos(phi) , sin(phi) * sin(theta) ) | \theta \in [ 0 , 2\pi ) , \phi \in [0,\pi]
		// d\Phi / d\theta = ( -sin(phi) * sin(theta) , 0         , sin(phi) * cos(theta) ) = sin(phi) * ( -sin(theta) , 0 , cos(theta) ) = sin(phi) * v_\theta
		// d\Phi / d\phi   = (  cos(phi) * cos(theta) , -sin(phi) , cos(phi) * sin(theta) ) = v_\phi
		// I = | sin^2(\phi) 0 |
		//     | 0           1 |
		// \nabla f = df/d\theta * 1/sin(\phi) * v_\theta + df/d\phi * v\_phi


		// Given \delta, define f(\theta,\phi) to be the function that is 1 at (\delta,\delta) and zero at (0,0) and (0,\delta)
		//		f(\theta,\phi) = \phi / delta * \theta / \delta
		// => \nabla f = (1/\delta^2)[ \phi / sin(\phi) * v_\theta + \theta * v_\phi ]
		// Taking the limit as \phi -> 0 gives:
		// lim_{\phi \rightarrow 0}\nabla f(0,phi) = (1/\delta^2) * v_\theta

		//		f(\theta,\phi) = a \theta + b \phi + c
		// Using the fact that it's zero at (0,0)
		//		f(\theta,\phi) = a \theta + b \phi
		// Using the fact that it's zero (\delta,\delta)
		//		f(\theta,\phi) = a ( \theta - \phi )
		// Using the fact that it's one at (0,\delta)
		//		f(\theta,\phi) = -1/\delta ( \theta - \phi )
		// Then:
		//		\nabla f = -1/\delta * 1 / \sin(\phi) * v_\theta + 1/\delta * v_\phi

		// If f is the linear function:
		//		f(\theta,\phi) = \theta * \alpha + \phi * \beta
		// \nabla f = \beta * v_\theta / \sin(\phi) 

		static Point3D< Real > SpherePoint( Point2D< Real > thetaPhi );
		static Point2D< Real > SphereParameters( Point3D< Real > p );
		static std::pair< Point3D< Real > , Point3D< Real > > SphereTangents( Real theta , Real phi );
		static std::pair< Point3D< Real > , Point3D< Real > > SphereTangents( Point3D< Real > p );
		static Point2D< Real > J( Point2D< Real > t );
		static Point2D< Real > J( const SquareMatrix< Real , 2 > &g , Point2D< Real > t );
		static SquareMatrix< Real , 2 > J( void );
		static SquareMatrix< Real , 2 > J( const SquareMatrix< Real , 2 > &g );

		static size_t VertexNum( size_t resolution );
		static size_t FaceNum( size_t resolution );
		static size_t VertexIndex( int i , int j , size_t resolution );
		static size_t FaceIndex( int i , int j , size_t resolution );
		static void FactorVertexIndex( size_t idx , int &i , int &j , size_t resolution );
		static void FactorFaceIndex( size_t idx , int &i , int &j , size_t resolution );

		static Real VertexMass( size_t idx , size_t resolution );
		static Real VertexMass( int i , int j , size_t resolution ){ return VertexMass( VertexIndex( i , j , resolution ) , resolution ); }
		static Real FaceMass( size_t idx , size_t resolution );
		static Real FaceMass( int i , int j , size_t resolution ){ return FaceMass( FaceIndex( i , j , resolution ) , resolution ); }
		static Point3D< Real > Vertex( size_t idx , size_t resolution );
		static Point3D< Real > Vertex( int i , int j , size_t resolution ){ return Vertex( VertexIndex( i , j , resolution ) , resolution ); }
		static Point3D< Real > FaceCenter( size_t idx , size_t resolution );
		static Point3D< Real > FaceCenter( int i , int j , size_t resolution ){ return FaceCenter( FaceIndex( i , j , resolution ) , resolution ); }
		static FaceIndices Face( size_t idx , size_t resolution );
		static FaceIndices Face( int i , int j , size_t resolution ){ return Face( FaceIndex( i , j , resolution ) , resolution ); }
		static std::pair< Point3D< Real > , Point3D< Real > > FaceTangents( size_t idx , size_t resolution );
		static std::pair< Point3D< Real > , Point3D< Real > > FaceTangents( int i , int j , size_t resolution ){ return FaceTangents( FaceIndex( i , j , resolution ) , resolution ); }
		static std::pair< Point3D< Real > , Point3D< Real > > VertexTangents( size_t idx , size_t resolution );
		static std::pair< Point3D< Real > , Point3D< Real > > VertexTangents( int i , int j , size_t resolution ){ return VertexTangents( VertexIndex( i , j , resolution ) , resolution ); }

		template< unsigned int Cols , unsigned int Rows >
		static SparseMatrix< Real , int > ExpandMatrix( const SparseMatrix< Matrix< Real , Cols , Rows > , int > &M );
		template< unsigned int Dim >
		static SparseMatrix< Real , int > ExpandMatrix( const SparseMatrix< SquareMatrix< Real , Dim > , int > &M );

		template< typename Data > struct PrimalSignal;
		template< typename Data > struct DualSignal;
		template< bool Gradient , bool JGradient > struct PrimalVectorField;

		template< typename Data >
		struct DualSignal
		{
			DualSignal( void ) : _resolution(0){};
			DualSignal( size_t resolution ) : _resolution(resolution) , _signal( FaceNum( resolution ) ){}
			DualSignal( size_t resolution , Data value ) : _resolution(resolution) , _signal( FaceNum( resolution ) , value ){}
			DualSignal( const PrimalSignal< Data > &pSignal ){ _setFromPrimal( &pSignal[0] , pSignal.resolution() ); }

			void set( const PrimalSignal< Data > &pSignal ){ _setFromPrimal( &pSignal[0] , pSignal.resolution() ); }

			size_t resolution( void ) const { return _resolution; }
			size_t size( void ) const { return _signal.size(); }
			void setResolution( size_t resolution ){ _resolution = resolution ; _signal.resize( FaceNum( resolution ) ); }
			Data &operator[]( size_t i ){ return _signal[i]; }
			const Data &operator[] ( size_t i ) const { return _signal[i]; }
			Data &operator()( int i , int j ){ return _signal[ FaceIndex(i,j,_resolution) ]; }
			const Data &operator() ( int i , int j ) const { return _signal[ FaceIndex(i,j,_resolution) ]; }
			std::vector< Data > &operator()( void ){ return _signal; }
			const std::vector< Data > &operator()( void ) const { return _signal; }

			Data operator()( Point3D< Real > p ) const;

			Real squareNorm( void ) const;

			size_t index( int i , int j ) const { return FaceIndex( i , j , _resolution ); }
			void factorIndex( size_t idx , int &i , int &j ) const { FactorFaceIndex( idx , i , j , _resolution ); }
			Point3D< Real > center( size_t idx ) const { return FaceCenter( idx , _resolution ); }
			Point3D< Real > center( int i , int j ) const { return center( FaceIndex(i,j,_resolution) ); }

			static SparseMatrix< Real , int > MassMatrix( size_t resolution );
		protected:
			void _setFromPrimal( const Data *values , size_t resolution );
			std::vector< Data > _signal;
			size_t _resolution;
		};

		template< typename Data >
		struct PrimalSignal
		{
			PrimalSignal( void ) : _resolution(0){};
			PrimalSignal( const SphericalGrid< Data > &sGrid ){ _setFromDual( &sGrid(0,0) , sGrid.resolution() ); }
			PrimalSignal( const DualSignal< Data > &dSignal ){ _setFromDual( &dSignal[0] , dSignal.resolution() ); }
			PrimalSignal( size_t resolution ) : _resolution(resolution) , _signal( VertexNum( resolution ) ){}
			PrimalSignal( size_t resolution , Data value ) : _resolution(resolution) , _signal( VertexNum( resolution ) , value ){}

			void set( const SphericalGrid< Data > &sGrid ){ _setFromDual( &sGrid(0,0) , sGrid.resolution() ); }
			void set( const DualSignal< Data > &dSignal ){ _setFromDual( &dSignal[0] , dSignal.resolution() ); }

			size_t resolution( void ) const { return _resolution; }
			size_t size( void ) const { return _signal.size(); }
			void setResolution( size_t resolution ){ _resolution = resolution ; _signal.resize( VertexNum( resolution ) ); }
			Data &operator[]( size_t i ){ return _signal[i]; }
			const Data &operator[] ( size_t i ) const { return _signal[i]; }
			Data &operator()( int i , int j ){ return _signal[ VertexIndex(i,j,_resolution) ]; }
			const Data &operator() ( int i , int j ) const { return _signal[ VertexIndex(i,j,_resolution) ]; }
			std::vector< Data > &operator()( void ){ return _signal; }
			const std::vector< Data > &operator()( void ) const { return _signal; }
			Point3D< Real > vertex( size_t idx ) const { return Vertex( idx , _resolution ); }
			Point3D< Real > vertex( int i , int j ) const { return vertex( VertexIndex(i,j,_resolution) ); }
			size_t index( int i , int j ) const { return VertexIndex( i , j , _resolution ); }
			void factorIndex( size_t idx , int &i , int &j ) const { FactorVertexIndex( idx , i , j , _resolution ); }

			Data operator()( Point3D< Real > p ) const;

			PrimalSignal upSample( void ) const;
			PrimalSignal downSample( void ) const;

			void setGradients( SphericalGrid< Point3D< Data > > &gradients ) const{ gradients.resize( _resolution ) ; _setGradients( &gradients[0][0] ); }
			void setGradients( DualSignal< Point3D< Data > > &gradients ) const{ gradients.setResolution( _resolution ) ; _setGradients( &gradients[0] ); }

			template< typename ColorFunctor >
			void write( const char *fileName , const ColorFunctor &colorFunctor ) const;
			template< typename ColorFunctor >
			void write( const char *fileName , const ColorFunctor &colorFunctor , const DualSignal< Point3D< Real > > &vectorField ) const;
			template< typename ColorFunctor , bool Divergence , bool Curl >
			void write( const char *fileName , const ColorFunctor &colorFunctor , const PrimalVectorField< Divergence , Curl > &vectorField ) const;

			static SparseMatrix< Real , int > MassMatrix( size_t resolution , unsigned int mcSamples );
			static SparseMatrix< Real , int > StiffnessMatrix( size_t resolution , unsigned int mcSamples );

			// MetricTensor2x2 is a functor taking the indices of a face and returning a 2x2 matrix
			template< typename MetricTensor2x2 > static SparseMatrix< Real , int > MassMatrix( size_t resolution , unsigned int mcSamples , const MetricTensor2x2 &metricTensor );
			template< typename MetricTensor2x2 > static SparseMatrix< Real , int > StiffnessMatrix( size_t resolution , unsigned int mcSamples , const MetricTensor2x2 &metricTensor );

			static SparseMatrix< Real , int > RestrictionMatrix( size_t coarseResolution );
		protected:
			template< typename Entry > static void _InitSystemMatrix( size_t resolution , SparseMatrix< Entry , int > &M );
			static size_t _RowOffset( size_t idx1 , size_t idx2 , size_t resolution );
			static Matrix< Real , 4 , 2 > _GradientOperator( int j , size_t resolution );
			void _setFromDual( const Data *values , size_t resolution );
			void _setGradients( Point3D< Data > *gradientValues ) const;
			std::vector< Data > _signal;
			size_t _resolution;
		};

		template< bool Gradient , bool JGradient >
		struct PrimalVectorField
		{
			static const unsigned int Blocks = ( Gradient ? 1 : 0 ) + ( JGradient ? 1 : 0 );
			template< typename Data >
			struct MassAndStiffnessOperator
			{
				template< bool _Gradient=Gradient , bool _JGradient=JGradient > typename std::enable_if< _Gradient ,       Real & >::type mWeightGradient( void )       { return _mWeights[0]; }
				template< bool _Gradient=Gradient , bool _JGradient=JGradient > typename std::enable_if< _Gradient , const Real & >::type mWeightGradient( void ) const { return _mWeights[0]; }
				template< bool _Gradient=Gradient , bool _JGradient=JGradient > typename std::enable_if< _Gradient ,       Real & >::type sWeightGradient( void )       { return _sWeights[0]; }
				template< bool _Gradient=Gradient , bool _JGradient=JGradient > typename std::enable_if< _Gradient , const Real & >::type sWeightGradient( void ) const { return _sWeights[0]; }

				template< bool _Gradient=Gradient , bool _JGradient=JGradient >	typename std::enable_if< _JGradient ,       Real & >::type mWeightJGradient( void )       { return _mWeights[ _Gradient ? 1 : 0 ]; }
				template< bool _Gradient=Gradient , bool _JGradient=JGradient > typename std::enable_if< _JGradient , const Real & >::type mWeightJGradient( void ) const { return _mWeights[ _Gradient ? 1 : 0 ]; }
				template< bool _Gradient=Gradient , bool _JGradient=JGradient > typename std::enable_if< _JGradient ,       Real & >::type sWeightJGradient( void )       { return _sWeights[ _Gradient ? 1 : 0 ]; }
				template< bool _Gradient=Gradient , bool _JGradient=JGradient > typename std::enable_if< _JGradient , const Real & >::type sWeightJGradient( void ) const { return _sWeights[ _Gradient ? 1 : 0 ]; }

				MassAndStiffnessOperator( size_t resolution , unsigned int mcSamples );
				// MetricTensor2x2 is a functor taking the indices of a face and returning a 2x2 matrix
				template< typename MetricTensor2x2 > MassAndStiffnessOperator( size_t resolution , unsigned int mcSamples , const MetricTensor2x2 &metricTensor );


				void Multiply( ConstPointer( Data ) x , Pointer( Data ) Mx ) const;
				SparseMatrix< Real , int > toMatrix( void ) const;
			protected:
				Real _mWeights[Blocks] , _sWeights[Blocks];
				SparseMatrix< Real , int > _S , _invM;
				std::vector< Data > _scratch;
			};

			PrimalVectorField( void ) : _resolution(0){};
			PrimalVectorField( size_t resolution ) : _resolution(resolution) , _signal( Blocks*VertexNum( resolution ) ){}
			PrimalVectorField( size_t resolution , Real value ) : _resolution(resolution) , _signal( Blocks*VertexNum( resolution ) , value ){}

			void setVertexValues( PrimalSignal< Point3D< Real > > &vf ) const;

			size_t resolution( void ) const { return _resolution; }
			size_t size( void ) const { return _signal.size(); }
			void setResolution( size_t resolution ){ _resolution = resolution ; _signal.resize( Blocks*VertexNum( resolution ) ); }
			Real &operator[]( size_t i ){ return _signal[i]; }
			const Real &operator[] ( size_t i ) const { return _signal[i]; }
			std::vector< Real > &operator()( void ){ return _signal; }
			const std::vector< Real > &operator()( void ) const { return _signal; }

			Point3D< Real > operator()( Point3D< Real > p ) const;
			Point3D< Real > vertexValue( size_t idx ) const;
			Point3D< Real > vertexValue( int i , int j ) const { return vertexValue( VertexIndex( i , j , _resolution ) ); }
			Point3D< Real > faceValue( size_t idx ) const;
			Point3D< Real > faceValue( int i , int j ) const { return faceValue( FaceIndex( i , j , _resolution ) ); }

			static SparseMatrix< Real , int > MassMatrix( size_t resolution , unsigned int mcSamples ){ return MassMatrix( resolution , mcSamples , []( int , int ){ return SquareMatrix< Real , 3 >::Identity(); } ); }
			// MetricTensor3x3 is a function taking the indices of a face and returning a 3x3 matrix
			template< typename MetricTensor3x3 > static SparseMatrix< Real , int > MassMatrix( size_t resolution , unsigned int mcSamples , const MetricTensor3x3 &metricTensor );

			std::vector< Real > mass( unsigned int mcSamples ) const { return mass( mcSamples , []( int , int ){ return SquareMatrix< Real , 2 >::IdentityMatrix(); } ); }
			// MetricTensor3x3 is a function taking the indices of a face and returning a 3x3 matrix
			template< typename MetricTensor3x3 > std::vector< Real > mass( unsigned int mcSamples , const MetricTensor3x3 &metricTensor ) const;

			template< bool _Gradient , bool _JGradient >
			static std::vector< Real > WeightedMass( const PrimalSignal< Real > &weights , const PrimalVectorField< _Gradient , _JGradient > &vf , unsigned int mcSamples ){ return WeightedMass( 1 , &weights , &vf , mcSamples ); }
			template< bool _Gradient , bool _JGradient >
			static std::vector< Real > WeightedMass( unsigned int Dim , const PrimalSignal< Real > *weights , const PrimalVectorField< _Gradient , _JGradient > *vf , unsigned int mcSamples );

			static std::vector< Real > Mass( const DualSignal< Point3D< Real > > &V , unsigned int mcSamples );

			static SparseMatrix< Real , int > RestrictionMatrix( size_t coarseResolution );
		protected:
			// MetricTensor2x2 is a functor taking the indices of a face and returning a 2x2 matrix
			template< typename MetricTensor2x2 > static SparseMatrix< Real , int > _MassMatrix( size_t resolution , unsigned int mcSamples , const MetricTensor2x2 &metricTensor );
			// MetricTensor2x2 is a functor taking the indices of a face and returning a 2x2 matrix
			template< typename MetricTensor2x2 > std::vector< Real > _mass( unsigned int mcSamples , const MetricTensor2x2 &metricTensor ) const;

			Point3D< Real > _sample( int i , int j , Real di , Real dj ) const;
			Point3D< Real > _sample( const SquareMatrix< Real , 2 > &g , int i , int j , Real di , Real dj ) const;
			static void _InitSystemMatrix( size_t resolution , SparseMatrix< Real , int > &M );
			static size_t _RowOffset( size_t idx1 , size_t idx2 , size_t resolution );
			std::vector< Real > _signal;
			size_t _resolution;
		};
	protected:
		template< unsigned int D >
		static MonomialIntegralTable< Real , D > _GetMonomialIntegralTable( int j , size_t resolution , unsigned int mcSamples );
		static void _BilinearFaceBasis( int j , size_t resolution , BilinearFunctionType< Real > polynomials[4] );
		static void _SetFaceMass( int j , size_t resolution , const MonomialIntegralTable< Real , 2 > &mit , Real mass[4][4] );
		static void _SetFaceStiffness( int j , size_t resolution , const MonomialIntegralTable< Real , 2 > &mit , Real stiffness[4][4] );
		static void _SetFaceStiffnesses( int j , size_t resolution , const MonomialIntegralTable< Real , 2 > &mit , SquareMatrix< Real , 2 > stiffnesses[4][4] );
		static Real _Dot( const float &d1 , const float &d2 );
		static Real _Dot( const double &d1 , const double &d2 );
		template< typename Data > static Real _Dot( const Data &d1 , const Data &d2 );
	};
}
#include "Misha/SphericalSignals.inl"
#endif // SPHERICAL_SIGNALS_INCLUDED