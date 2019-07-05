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
	template< typename Real , bool DivergenceFree >
	OpticalFlow< Real , DivergenceFree >::OpticalFlow( SphericalGeometry::Mesh< Real > &source , SphericalGeometry::Mesh< Real > &target , const std::vector< Real > &sourceSignal , const std::vector< Real > &targetSignal , Parameters params )
		: _source(source) , _target(target) , _sourceSignal(sourceSignal) , _targetSignal(targetSignal) , _params(params)
	{
		_sGrid.resize( _params.resolution );
		_sGrids[0].resize( _params.resolution );
		_sGrids[1].resize( _params.resolution );
		_signals[0].setResolution( _params.resolution );
		_signals[1].setResolution( _params.resolution );
		_gradients[0].setResolution( _params.resolution );
		_gradients[1].setResolution( _params.resolution );
		_dSignal.setResolution( _params.resolution );
		_sGradients.setResolution( _params.resolution );
		_hForm.resize( (int)_params.resolution );
		_key.resize( (int)_params.resolution );
		_x.setResolution( (int)_params.resolution );
		__x.setResolution( (int)_params.resolution );

		typename FlowField::template MassAndStiffnessOperator< Real > massAndStiffness( _params.resolution , _params.integrationSamples );
		_setWeights( massAndStiffness , 0 , 0 , _params.smoothingWeight );
		_MassAndStiffness = massAndStiffness.toMatrix();

		_solver = NULL;

		_rasterizeSignals();
	}

	template< typename Real , bool DivergenceFree >
	OpticalFlow< Real , DivergenceFree >::~OpticalFlow( void )
	{
		if( _solver ) delete _solver;
		_solver = NULL;
	}

	template< typename Real , bool DivergenceFree >
	void OpticalFlow< Real , DivergenceFree >::advance( void )
	{
		_setFlowField();
		_advect();
		_rasterizeSignals();
	}

	template< typename Real , bool DivergenceFree >
	typename OpticalFlow< Real , DivergenceFree >::Stats OpticalFlow< Real , DivergenceFree >::getStats( void ) const
	{
		Stats stats;

		stats.l1Norm = _sGrids[0].l1Norm() + _sGrids[1].l1Norm();
		stats.l2Norm = _sGrids[0].squareL2Norm() + _sGrids[1].squareL2Norm();
		stats.l1Difference = SphericalGrid< Real >::L1Difference( _sGrids[0] , _sGrids[1] );
		stats.l2Difference = SphericalGrid< Real >::SquareL2Difference( _sGrids[0] , _sGrids[1] );
		stats.l2Norm = (Real)sqrt( stats.l2Norm );
		stats.l2Difference = (Real)sqrt( stats.l2Difference );
		return stats;
	}

	template< typename Real , bool DivergenceFree >
	Point3D< Real > OpticalFlow< Real , DivergenceFree >::flowField( Point3D< Real > p ) const { return _flowField(p); }

	template< typename Real , bool DivergenceFree >
	void OpticalFlow< Real , DivergenceFree >::_rasterizeSignals( void )
	{
		auto SetSphericalGrid = []( const SphericalGeometry::Mesh< Real > &sMesh , const std::vector< Real > &values , SphericalGrid< Real > &sGrid , Real epsilon , Real sigma , HarmonicTransform< Real > &hForm , FourierKeyS2< Real > &key )
		{
			SphericalGeometry::Tessellation< Real >::SampleVertexValues( sMesh , values , sGrid , epsilon );
			if( sigma>0 )
			{
				hForm.ForwardFourier( sGrid , key );
				for( int i=0 ; i<key.bandWidth() ; i++ ) for( int j=0 ; j<=i ; j++ ) key(i,j) *= (Real)exp( -sigma * i * ( (Real)(i+1) ) );
				hForm.InverseFourier( key , sGrid );
			}
		};

		for( int i=0 ; i<2 ; i++ )
		{
			SetSphericalGrid( i==0 ? _source : _target , i==0 ? _sourceSignal : _targetSignal , _sGrids[i] , _params.epsilon , _params.scalarSigma() , _hForm , _key );
			_signals[i].set( _sGrids[i] );
#pragma omp parallel for
			for( int j=0 ; j<_signals[i].size() ; j++ ) _gradients[i][j] = _signals[i][j];
		}
#pragma omp parallel for
		for( int i=0 ; i<_dSignal.size() ; i++ ) _dSignal[i] = _signals[1][i] - _signals[0][i] , _sGradients[i] = _gradients[1][i] + _gradients[0][i];
	}

	template< typename Real , bool DivergenceFree >
	void OpticalFlow< Real , DivergenceFree >::_setFlowField( void )
	{
		auto MetricTensor = [&]( int i , int j )
		{
			SquareMatrix< Real , 3 > tensor;
			tensor += OuterProduct( _gradients[0].faceValue(i,j) ) + OuterProduct( _gradients[1].faceValue(i,j) );
			return tensor;
		};

		// Noting that infinitesimal advection of a scalar field f by a vector field W acts as -< \nabla f , V > we are looking for the vector field V minimizing:
		//		E(V) = || - < V , \nabla s[0] > - (s[1]-s[0]) ||^2 + || < V , \nabla s[1] > - (s[0]-s[1]) ||^2 + \epsilon || V ||^2 + \delta || \nabla V ||^2
		//		     = || < V , \nabla s[0] > + (s[1]-s[0]) ||^2 + || < V , \nabla s[1] > + (s[1]-s[0]) ||^2 + \epsilon || V ||^2 + \delta || \nabla V ||^2
		//		     = V^t ( \nabla s[0] \nabla s[0]^2 ) V + V^t ( \nabla s[1] \nabla s[1]^t ) V + 2 V^t \nabla (s[0]+s[1])*(s[1]-s[0]) + \epsilon || V ||^2 + \delta || \nabla V ||^2 + ...
		//		     = V^t ( \nabla s[0] \nabla s[0]^2 + \nabla s[1] \nabla s[1]^2 + \epsilon + \delta \Delta ) V + 2 V^t \nabla (s[0]+s[1])*(s[1]-s[0]) ...
		// Taking the derivative w.r.t. V gives:
		//		   0 = 2 ( \nabla s[0] \nabla s[0]^t + \nabla s[1] \nabla s[1]^t + \epsilon + \delta \Delta ) V + 2 \nabla(s[0]+s[1])*(s[1]-s[0])
		// <=>     -\nabla(s[0]+s[1])*(s[1]-s[0]) = ( \nabla s[0] \nabla s[0]^t + \nabla s[1] \nabla s[1]^t + \epsilon + \delta \Delta ) V

		_b = FlowField::WeightedMass( 1 , &_dSignal , &_sGradients , _params.integrationSamples );
#pragma omp parallel for
		for( int i=0 ; i<_b.size() ; i++ ) _b[i] = -_b[i];

		_FlowFieldMass = FlowField::MassMatrix( _params.resolution , _params.integrationSamples , MetricTensor );
		SparseMatrix< Real , int > M = _MassAndStiffness + _FlowFieldMass;

#if 1
		// [WARNING] Massaging matrix to make it definite
		for( int i=0 ; i<M.rows ; i++ ) for( int j=0 ; j<M.rowSizes[i] ; j++ ) if( M[i][j].N==i ) M[i][j].Value += 1e-8;
#endif

		if( !_solver ) _solver = new Solver( M , true );
		_solver->update( M );
		_solver->solve( GetPointer( _b ) , GetPointer( _x() ) );

		if( _params.bandLimitDampening>0 )
		{
			Real sigma = _params.flowFieldSigma();
			{
#pragma omp parallel for
				for( int i=0 ; i<__x.size() ; i++ ) __x[i] = _x[i];
#pragma omp parallel for
				for( int j=0 ; j<(int)_params.resolution ; j++ ) for( int i=0 ; i<(int)_params.resolution ; i++ ) _sGrid(i,j) = __x( SphericalGeometry::RegularGrid< Real >::FaceCenter(i,j,_params.resolution) );
				_hForm.ForwardFourier( _sGrid , _key );
				for( int i=0 ; i<_key.bandWidth() ; i++ ) for( int j=0 ; j<=i ; j++ ) _key(i,j) *= (Real)exp( -sigma * i * ( (Real)(i+1) ) );
				_hForm.InverseFourier( _key , _sGrid );
				__x.set( _sGrid );
#pragma omp parallel for
				for( int i=0 ; i<__x.size() ; i++ ) _x[i] = __x[i];
			}
			if( !DivergenceFree )
			{
#pragma omp parallel for
				for( int i=0 ; i<__x.size() ; i++ ) __x[i] = _x[i+__x.size()];
#pragma omp parallel for
				for( int j=0 ; j<(int)_params.resolution ; j++ ) for( int i=0 ; i<(int)_params.resolution ; i++ ) _sGrid(i,j) = __x( SphericalGeometry::RegularGrid< Real >::FaceCenter(i,j,_params.resolution) );
				_hForm.ForwardFourier( _sGrid , _key );
				for( int i=0 ; i<_key.bandWidth() ; i++ ) for( int j=0 ; j<=i ; j++ ) _key(i,j) *= (Real)exp( -sigma * i * ( (Real)(i+1) ) );
				_hForm.InverseFourier( _key , _sGrid );
				__x.set( _sGrid );
#pragma omp parallel for
				for( int i=0 ; i<__x.size() ; i++ ) _x[i+__x.size()] = __x[i];
			}
		}

		// Solve for the optimal scaling term assuming no smoothing
		if( _params.rescaleFlow )
		{
			Real v1=0 , v2=0;
			for( int i=0 ; i<_FlowFieldMass.rows ; i++ ) for( int j=0 ; j<_FlowFieldMass.rowSizes[i] ; j++ ) v1 += _FlowFieldMass[i][j].Value * _x[i] * _x[ _FlowFieldMass[i][j].N ];
			for( int i=0 ; i<_FlowFieldMass.rows ; i++ ) v2 += _b[i] * _x[i];
			if( v1 )
			{
				Real scale = v2/v1;
#pragma omp parallel for
				for( int i=0 ; i<_x.size() ; i++ ) _x[i] *= scale;
			}
		}
		_x.setVertexValues( _flowField );
#pragma omp parallel for
		for( int i=0 ; i<_flowField.size() ; i++ ) _flowField[i] *= _params.stepSize;
	}

	template< typename Real , bool DivergenceFree >
	void OpticalFlow< Real , DivergenceFree >::_advect( void )
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

		Real scale = _params.symmetric ? (Real)1./2 : (Real)1.;
		scale /= _params.subSteps;
#pragma omp parallel for
		for( int i=0 ; i<_source.vertices.size() ; i++ )
		{
			Point3D< Real > p = _source.vertices[i];
			if( _params.useSemiImplicit )
			{
				for( unsigned int s=0 ; s<_params.subSteps ; s++ )
				{
					Frame start , end;
					Point3D< Real > v1 , v2;
					// Get the vector at the start position
					v1 = _flowField( p ) * scale;
					if( AdvectFrame( p , v1 , start , end ) )
					{
						// Get the vector at the end and transport it back
						v2 = _flowField( end.p ) * scale;
						v2 = Frame::Transport( end , start , v2 );
						v1 = (v1 + v2) / 2;
					}
					p = Advect( p , v1 );
				}
			}
			else for( unsigned int s=0 ; s<_params.subSteps ; s++ ) p = Advect( p , _flowField( p ) * scale );
			_source.vertices[i] = p;
		}
		if( _params.symmetric )
			for( int i=0 ; i<_target.vertices.size() ; i++ )
			{
				Point3D< Real > p = _target.vertices[i];
				if( _params.useSemiImplicit )
				{
					for( unsigned int s=0 ; s<_params.subSteps ; s++ )
					{
						Frame start , end;
						Point3D< Real > v1 , v2;
						// Get the vector at the start position
						v1 = - _flowField( p ) * scale;
						if( AdvectFrame( p , v1 , start , end ) )
						{
							// Get the vector at the end and transport it back
							v2 = - _flowField( end.p ) * scale;
							v2 = Frame::Transport( end , start , v2 );
							v1 = (v1 + v2) / 2;
						}
						p = Advect( p , v1 );
					}
				}
				else for( unsigned int s=0 ; s<_params.subSteps ; s++ ) p = Advect( p , - _flowField( p ) * scale );
				_target.vertices[i] = p;
			}
	}
}
