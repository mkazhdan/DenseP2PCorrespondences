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
	void SO3Alignment< Real >::SetCorrelation( unsigned int gridNum , SphericalGrid< Real > source[] , SphericalGrid< Real > target[] , Parameters parameters , RotationGrid< Real > &correlation )
	{
		FourierKeySO3< Real> keySO3; 
		HarmonicTransform< Real > hForm; 
		WignerTransform< Real > wForm; 
		FourierKeyS2<Real> sourceKey , targetKey; 
		for( unsigned int g=0 ; g<gridNum ; g++ ) if( source[g].resolution()!=target[g].resolution() || source[g].resolution()!=source[0].resolution() ) fprintf( stderr , "[ERROR] Resolutions differ\n" ) , exit( 0 );

		int resolution = source[0].resolution();

		// Allocate memory 
		if( !keySO3.resize( resolution ) ) fprintf( stderr , "[ERROR] SO3Alignment< Real >::SetCorrelation: Could not allocate key: %d\n" , resolution ) , exit( 0 ); 
		for( int i=0 ; i<keySO3.bandWidth(); i++ ) for( int j=0 ; j<=i ; j++ ) for( int k=0 ; k<=i ; k++ ) keySO3(i,j,k) = 0;

		Real sigma = parameters.sigma();
		for( unsigned int g=0 ; g<gridNum ; g++ )
		{
			hForm.ForwardFourier( source[g] , sourceKey );
			hForm.ForwardFourier( target[g] , targetKey );
			if( sigma>0 )
			{
#pragma omp parallel for
				for( int i=0 ; i<sourceKey.bandWidth() ; i++ )
				{
					Real scl = (Real)exp( -sigma * i * ( (Real)(i+1) ) );
					for( int j=0 ; j<=i ; j++ ) sourceKey(i,j) *= scl , targetKey(i,j) *= scl;
				}
			}
			if( parameters.removeDC ) sourceKey(0,0) *= (Real)0 , targetKey(0,0) *= (Real)0;

			// Compute the Wigner-D coefficients
#pragma omp parallel for
			for( int i=0 ; i<keySO3.bandWidth() ; i++ ) for( int j=0 ; j<=i ; j++ ) for( int k=0 ; k<=i ; k++ ) 
			{
				// set each Wigner-D coefficient to be the cross multiplication within each frequency 
				keySO3(i,j,k) += sourceKey(i,k) * targetKey(i,j).conjugate();
				if( k ) keySO3(i,j,-k) += sourceKey(i,k).conjugate() * targetKey(i,j).conjugate();
			}
		}

		// Take the inverse Wigner D transform 
		wForm.InverseFourier( keySO3 , correlation );
	}

	template< typename Real >
	unsigned int SO3Alignment< Real >::TopNRotations( unsigned int gridNum , SphericalGrid< Real > source[] , SphericalGrid< Real > target[] , Parameters parameters , std::vector< Real > &correlationValues , std::vector< SquareMatrix< Real , 3 > > &rotations )
	{
		RotationGrid< Real > correlation;
		SetCorrelation( gridNum , source , target , parameters , correlation );
		return TopNRotations( correlation , parameters , correlationValues , rotations );
	}

	template< typename Real >
	unsigned int SO3Alignment< Real >::TopNRotations( const RotationGrid< Real > &correlation , Parameters parameters , std::vector< Real > &correlationValues , std::vector< SquareMatrix< Real , 3 > > &rotations )
	{
		correlationValues.resize( 0 );
		rotations.resize( 0 );

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
			for( unsigned int i=0 ; i<n ; i++ ) if( SquareMatrix< Real , 3 >::SquareNorm( rotations[i].inverse() * R - SquareMatrix< Real , 3 >::Identity() )<parameters.rotationSeparation ) return false;
			if( n && correlation<correlationValues[0]*parameters.correlationFraction ) return false;
			return true;
		};

		for( unsigned int n=0 ; n<parameters.maxRotations ; n++ )
		{
			bool found = false;
			Real maxCorrelation = 0;
			SquareMatrix< Real , 3 > R;
			for( int x=0 ; x<correlation.resolution() ; x++ ) for( int y=0 ; y<correlation.resolution() ; y++ ) for( int z=0 ; z<correlation.resolution() ; z++ ) if( correlation(x,y,z)>maxCorrelation )
			{
				SquareMatrix< Real , 3 > _R = RotationMatrix( x , y , z );
				if( parameters.correlationFraction<0 || ValidRotation( _R , n , correlation(x,y,z) ) )
				{
					maxCorrelation = correlation( x , y , z );
					R = _R;
					found = true;
				}
			}
			if( !found ) return n;
			correlationValues.push_back( maxCorrelation );
			rotations.push_back( R );
		}
		return (unsigned int)correlationValues.size();
	}


	template< typename Real >
	void SO3Alignment< Real >::AlignVertexValues( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > sourceSignals[] , const std::vector< Real > targetSignals[] , unsigned int signalNum , Parameters parameters , RotationGrid< Real > &correlation )
	{
		std::vector< SphericalGrid< Real > > sourceSGrids , targetSGrids;
		for( int m=0 ; m<2 ; m++ )
		{
			const SphericalGeometry::Mesh< Real > &sMesh = m==0 ? source : target; 
			std::vector< SphericalGrid< Real > > &sGrids = m==0 ? sourceSGrids : targetSGrids;
			const std::vector< Real > *signals = m==0 ? sourceSignals : targetSignals;
			sGrids.resize( signalNum );
			for( unsigned int i=0 ; i<signalNum ; i++ )
			{
				sGrids[i].resize( (int)parameters.resolution );
				SphericalGeometry::Tessellation< Real >::SampleVertexValues( sMesh.vertices , sMesh.polygons , signals[i] , sGrids[i] , parameters.epsilon );
			}
		}
		SetCorrelation( signalNum , &sourceSGrids[0] , &targetSGrids[0] , parameters , correlation );
	}

	template< typename Real >
	void SO3Alignment< Real >::AlignFaceIntegrals( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > sourceSignals[] , const std::vector< Real > targetSignals[] , unsigned int signalNum , Parameters parameters , RotationGrid< Real > &correlation )
	{
		std::vector< SphericalGrid< Real > > sourceSGrids , targetSGrids;
		for( int m=0 ; m<2 ; m++ )
		{
			const SphericalGeometry::Mesh< Real > &sMesh = m==0 ? source : target; 
			std::vector< SphericalGrid< Real > > &sGrids = m==0 ? sourceSGrids : targetSGrids;
			const std::vector< Real > *signals = m==0 ? sourceSignals : targetSignals;
			sGrids.resize( signalNum );
			for( unsigned int i=0 ; i<signalNum ; i++ )
			{
				sGrids[i].resize( (int)parameters.resolution );
				SphericalGeometry::Tessellation< Real >::SampleFaceIntegrals( sMesh.vertices , sMesh.polygons , signals[i] , sGrids[i] , parameters.epsilon );
			}
		}
		SetCorrelation( signalNum , &sourceSGrids[0] , &targetSGrids[0] , parameters , correlation );
	}

	template< typename Real >
	void SO3Alignment< Real >::AlignVertexValues( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > sourceSignals[] , const std::vector< Real > targetSignals[] , unsigned int signalNum , Parameters parameters , std::vector< SquareMatrix< Real , 3 > > &rotations , std::vector< Real > &correlationValues )
	{
#if 1
		RotationGrid< Real > correlation;
		AlignVertexValues( source , target , sourceSignals , targetSignals , signalNum , parameters , correlation );
		TopNRotations( correlation , parameters , correlationValues , rotations );
#else
		std::vector< SphericalGrid< Real > > sourceSGrids , targetSGrids;
		for( int m=0 ; m<2 ; m++ )
		{
			const SphericalGeometry::Mesh< Real > &sMesh = m==0 ? source : target; 
			std::vector< SphericalGrid< Real > > &sGrids = m==0 ? sourceSGrids : targetSGrids;
			const std::vector< Real > *signals = m==0 ? sourceSignals : targetSignals;
			sGrids.resize( signalNum );
			for( unsigned int i=0 ; i<signalNum ; i++ )
			{
				sGrids[i].resize( (int)parameters.resolution );
				SphericalGeometry::Tessellation< Real >::SampleVertexValues( sMesh.vertices , sMesh.polygons , signals[i] , sGrids[i] , parameters.epsilon );
			}
		}

		TopNRotations( parameters.maxRotations , signalNum , &sourceSGrids[0] , &targetSGrids[0] , parameters.sigma() , parameters.removeDC , parameters.rotationSeparation , parameters.correlationFraction , correlationValues , rotations );
#endif
	}

	template< typename Real >
	void SO3Alignment< Real >::AlignFaceIntegrals( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > sourceSignals[] , const std::vector< Real > targetSignals[] , unsigned int signalNum , Parameters parameters , std::vector< SquareMatrix< Real , 3 > > &rotations , std::vector< Real > &correlationValues )
	{
#if 1
		RotationGrid< Real > correlation;
		AlignFaceIntegrals( source , target , sourceSignals , targetSignals , signalNum , parameters , correlation );
		TopNRotations( correlation , parameters , correlationValues , rotations );
#else
		std::vector< SphericalGrid< Real > > sourceSGrids , targetSGrids;
		for( int m=0 ; m<2 ; m++ )
		{
			const SphericalGeometry::Mesh< Real > &sMesh = m==0 ? source : target; 
			std::vector< SphericalGrid< Real > > &sGrids = m==0 ? sourceSGrids : targetSGrids;
			const std::vector< Real > *signals = m==0 ? sourceSignals : targetSignals;
			sGrids.resize( signalNum );
			for( unsigned int i=0 ; i<signalNum ; i++ )
			{
				sGrids[i].resize( (int)parameters.resolution );
				SphericalGeometry::Tessellation< Real >::SampleFaceIntegrals( sMesh.vertices , sMesh.polygons , signals[i] , sGrids[i] , parameters.epsilon );
			}
		}

		TopNRotations( parameters.maxRotations , signalNum , &sourceSGrids[0] , &targetSGrids[0] , parameters.sigma , parameters.removeDC , parameters.rotationSeparation , parameters.correlationFraction , correlationValues , rotations );
#endif
	}

}
