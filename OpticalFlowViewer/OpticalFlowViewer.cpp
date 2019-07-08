/*
Copyright (c) 2019, Michael Kazhdan and Sing Chun Lee
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

#include <cstdlib>
#include <vector>
#include <omp.h>
#include "Misha/CmdLineParser.h"
#include "Misha/SourceToTargetCorrespondence.h"
#include "Misha/Signals.h"
#include "Misha/SO3Alignment.h"
#ifdef NO_OPEN_GL
#include "Misha/OpticalFlow.h"
#else // !NO_OPEN_GL
#include "OpticalFlowViewer.inl"
#endif // NO_OPEN_GL


cmdLineParameterArray< char * , 2 >
	In( "in" ) ,										// Input spherical meshes with signals encoded in the normals
	InSignal( "signal" ) ,								// Input signals (assumed to be 80 channels)
	Out( "out" );										// Aligned spherical mesh headers
cmdLineParameter< char * >
	InCamera( "camera" );								// The initial camera
cmdLineParameter< float >
	SmoothingWeight( "sWeight" , 0.05f ) ,				// Flow field smoothness weight
	RotationSeparation( "rSeparation" , 3.f ) ,			// Minimum frobenius norm for rotation separation
	CorrelationFraction( "cFraction" , 0.9f ) ,			// Correlation proximity to peak
	BandLimitDampening( "filterScale" , 1e-2f ) ,		// Scale for band-pass filtering
	StepSize( "stepSize" , 1.f );						// Fraction of the flow step to take
cmdLineParameter< int >
	BaseResolution( "res" , 16 ) ,						// Resolution of the coarsest sphere
	RotationalAlignmentResolution( "rotRes" ) ,			// Resolution for performing rotational alignment
	SolvesPerLevel( "solvesPerLevel" , 6 ) ,			// Number of optical flow iterations per level of the hierarchy
	AdvectionSubSteps( "subSteps" , 1 ) ,				// Number of sub-steps to use for performing advection
	IntegrationSamples( "iSamples" , 1024 ) ,			// Number of samples to use in performing Monte-Carlo integration
	Rotate( "rotate" , -1 ) ,							// Use the n-th best rotation
	MaximumCandidateRotations( "maxRot" , 4 ) ,			// Maximum number of rotations
	Threads( "threads" , omp_get_num_procs() );			// Number of threads to use for parallelization
cmdLineReadable
	UseL2Difference( "l2Difference" ) ,					// Select optimal alignment by L2 difference
	NoDivergenceFree( "noDivFree" ) ,					// Use both curl-free and divergence-free components for the flow field
	RescaleFlow( "rescaleFlow" ) ,						// Adjust the advection step-size by estimating the step-size without smoothness constraints
	Animate( "animate" ) ,								// Start the visualization running the optical flow
	ShowEdges( "edges" ) ,								// Start the visualization running in edge mode
	ShowColors( "color" ) ,								// Start the visualization running in color mode
	Verbose( "verbose" ) ,								// Print out general statistics to the command line
	OutputAllRotations( "allRotations" ) ,				// Output the result from all rotation candidates
	FullVerbose( "fullVerbose" );						// Print out per-solve statistics to the command line

void Usage( const char* ex ) 
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input source/target spherical meshes>\n" , In.name );
	printf( "\t --%s <input source/target signals>\n" , InSignal.name );
	printf( "\t[--%s <output source/target mesh headers>]\n" , Out.name );
	printf( "\t[--%s <initial camera>]\n" , InCamera.name );
	printf( "\t[--%s <starting spherical grid resolution>=%d]\n" , BaseResolution.name , BaseResolution.value );
	printf( "\t[--%s <rotational alignment resolution>]\n" , RotationalAlignmentResolution.name );
	printf( "\t[--%s <smoothing weight>=%g]\n" , SmoothingWeight.name , SmoothingWeight.value );
	printf( "\t[--%s <integration samples>=%d]\n" , IntegrationSamples.name , IntegrationSamples.value );
	printf( "\t[--%s <number of threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s <band limit dampening factor>=%g]\n" , BandLimitDampening.name , BandLimitDampening.value );
	printf( "\t[--%s <step size>=%g]\n" , StepSize.name , StepSize.value );
	printf( "\t[--%s <solves at each levels>=%d]\n" , SolvesPerLevel.name , SolvesPerLevel.value );
	printf( "\t[--%s <number of advection sub-steps>=%d]\n" , AdvectionSubSteps.name , AdvectionSubSteps.value );
	printf( "\t[--%s <minimal rotation separation>=%f]\n" , RotationSeparation.name , RotationSeparation.value );
	printf( "\t[--%s <correlation fraction>=%f]\n" , CorrelationFraction.name , CorrelationFraction.value );
	printf( "\t[--%s <best rotation>]\n" , Rotate.name );
	printf( "\t[--%s <maximum candidate rotations>=%d]\n" , MaximumCandidateRotations.name , MaximumCandidateRotations.value );
	printf( "\t[--%s]\n" , UseL2Difference.name );
	printf( "\t[--%s]\n" , NoDivergenceFree.name );
	printf( "\t[--%s]\n" , RescaleFlow.name );
	printf( "\t[--%s]\n" , Animate.name );
	printf( "\t[--%s]\n" , ShowEdges.name );
	printf( "\t[--%s]\n" , ShowColors.name );
	printf( "\t[--%s]\n" , OutputAllRotations.name );
	printf( "\t[--%s]\n" , Verbose.name );
	printf( "\t[--%s]\n" , FullVerbose.name );
}
cmdLineReadable* params[] = { &In , &InSignal , &Out , &InCamera , &BaseResolution , &RotationalAlignmentResolution , &SmoothingWeight , &IntegrationSamples , &Threads , &StepSize , &Rotate , &NoDivergenceFree , &AdvectionSubSteps , &RescaleFlow , &SolvesPerLevel , &Verbose , &FullVerbose , &Animate , &ShowEdges , &ShowColors , &RotationSeparation , &CorrelationFraction , &MaximumCandidateRotations , &UseL2Difference , &BandLimitDampening , &OutputAllRotations , NULL };

#ifdef NO_OPEN_GL
#else // !NO_OPEN_GL
template< bool DivFree >
struct OpticalFlowVisualizationWrapper
{
	static OpticalFlowVisualization< double , DivFree > *ofv;
	static void Idle             ( void )                                  { return ofv->Idle(); }
	static void Display          ( void )                                  { return ofv->Display(); }
	static void KeyboardFunc     ( unsigned char key , int x , int y )     { return ofv->KeyboardFunc( key , x , y ); }
	static void SpecialFunc      ( int key , int x , int y )               { return ofv->SpecialFunc( key , x , y ); }
	static void Reshape          ( int w , int h )                         { return ofv->Reshape( w , h ); }
	static void MouseFunc        ( int button , int state , int x , int y ){ return ofv->MouseFunc( button , state , x , y ); }
	static void MotionFunc       ( int x , int y )                         { return ofv->MotionFunc( x , y ); }
	static void PassiveMotionFunc( int x , int y )                         { return ofv->PassiveMotionFunc( x , y ); }
};

template< bool DivFree >
OpticalFlowVisualization< double , DivFree > *OpticalFlowVisualizationWrapper< DivFree >::ofv;
#endif // NO_OPEN_GL

template< class Real >
Point3D< Real > NormalColor( Point3D< Real > n )
{
//	if( InvertNormal.set )  n = - n;
	Point3D< Real > c = ( -n + Point3D< Real >( 1 , 1 , 1 ) ) * Real(128);
	for( int d=0 ; d<3 ; d++ )
	{
		if( c[d]>Real(255) ) c[d] = Real(255);
		if( c[d]<Real(  0) ) c[d] = Real(  0);
	}
	return c;
}

template< typename Real >
std::vector< SquareMatrix< Real , 3 > > GetAligningRotations( const SphericalGeometry::Mesh< Real > &source , const SphericalGeometry::Mesh< Real > &target , const std::vector< Real > &sourceValues , const std::vector< Real > &targetValues , typename SphericalGeometry::SO3Alignment< Real >::Parameters parameters , bool verbose )
{
	std::vector< SquareMatrix< Real , 3 > > candidateRotations;
	std::vector< Real > l2CorrelationValues;
	SphericalGeometry::SO3Alignment< Real >::AlignVertexValues( source , target , sourceValues , targetValues , parameters , candidateRotations , l2CorrelationValues );

#if 1
	if( verbose ) for( int i=0 ; i<l2CorrelationValues.size() ; i++ )
	{
		printf( "Correlation values[%d] %g (%f)\t" , i+1 , l2CorrelationValues[i] , l2CorrelationValues[i]/l2CorrelationValues[0] );
		for( int j=0 ; j<i ; j++ ) printf( " %g" , SquareMatrix< Real , 3 >::SquareNorm( candidateRotations[i].inverse() * candidateRotations[j] - SquareMatrix< Real , 3 >::Identity() ) );
		printf( "\n" );
	}
#else
	if( verbose ) for( int i=0 ; i<l2CorrelationValues.size() ; i++ ) printf( "Correlation values[%d] %g (%f)\n" , i+1 , l2CorrelationValues[i] , l2CorrelationValues[i]/l2CorrelationValues[0] );
#endif
	return candidateRotations;
}

template< typename Real >
void ApplyRotation( SphericalGeometry::Mesh< Real > &source , SphericalGeometry::Mesh< Real > &target , SquareMatrix< Real , 3 > R )
{
	R = RotationSquareRoot< Real >( R );
#pragma omp parallel for
	for( int i=0 ; i<source.vertices.size() ; i++ ) source.vertices[i] = R * source.vertices[i];
	R = R.transpose();
#pragma omp parallel for
	for( int i=0 ; i<target.vertices.size() ; i++ ) target.vertices[i] = R * target.vertices[i];
}

template< typename Real , bool DivFree >
int _main( int argc , char *argv[] )
{
	// Read the mesh in and set the signal
	SphericalGeometry::Mesh< Real > sMeshes[2];
	std::vector< Point3D< Real > > meshVertices[2];
	Signal< Real > signals[2];
	std::vector< Point3D< Real > > colors[2];

	for( int i=0 ; i<2 ; i++ )
	{
		if( !sMeshes[i].read( In.values[i] , meshVertices[i] , colors[i] ) )
		{
			colors[i].resize( sMeshes[i].vertices.size() );

			for( int p=0 ; p<sMeshes[i].polygons.size() ; p++ )
			{
				Point3D< Real > n;
				const std::vector< int > &poly = sMeshes[i].polygons[p];
				for( int j=0 ; j<poly.size() ; j++ )
				{
					Point3D< Real > b = meshVertices[i][ poly[ (j+1)%poly.size() ] ] + meshVertices[i][ poly[j] ];
					Point3D< Real > e = meshVertices[i][ poly[ (j+1)%poly.size() ] ] - meshVertices[i][ poly[j] ];
					n += Point3D< Real >::CrossProduct( b , e ) / 2;
				}
				for( int j=0 ; j<poly.size() ; j++ ) colors[i][ poly[j] ] += n;
			}
#pragma omp parallel for
			for( int v=0 ; v<colors[i].size() ; v++ ) colors[i][v] = NormalColor( -colors[i][v]/(Real)Length( colors[i][v] ) );
		}
		signals[i].read( InSignal.values[i] );
	}
	if( signals[0].values.size()!=signals[1].values.size() ) fprintf( stderr , "[ERROR] Source and target signals have different dimensions: %d != %d\n" , (int)signals[0].values.size() , (int)signals[1].values.size() ) , exit( 0 );
	if( signals[0].values.size()<3 ) fprintf( stderr , "[ERROR] Expected at least three channels in the signal: %d >= 3\n" , (int)signals[0].values.size() ) , exit( 0 );
	if( !RotationalAlignmentResolution.set ) RotationalAlignmentResolution.value = BaseResolution.value << ( signals[0].values.size()-3 );

	// Get aliging rotations
	std::vector< SquareMatrix< Real , 3 > > candidateRotations;
	{
		Timer timer;
		typename SphericalGeometry::SO3Alignment< Real >::Parameters parameters;
		parameters.resolution = RotationalAlignmentResolution.value;
		parameters.bandLimitDampening = BandLimitDampening.value;
		parameters.maxRotations = MaximumCandidateRotations.value;
		parameters.rotationSeparation = RotationSeparation.value;
		parameters.correlationFraction = CorrelationFraction.value;
		parameters.removeDC = true;

		candidateRotations = GetAligningRotations< Real >( sMeshes[0] , sMeshes[1] , signals[0].values[0] , signals[1].values[0] , parameters , Verbose.set );
		if( Verbose.set ) printf( "Rotational alignment time: %.2f(s)\n" , timer.elapsed() );
	}

#ifdef NO_OPEN_GL
#else // !NO_OPEN_GL
	if( Out.set )
#endif // NO_OPEN_GL
	{
		auto OutputAdvectedMesh = []( const char *inFileName , const char *outFileName , const std::vector< Point3D< Real > > &advectedVertices )
		{
			SphericalGeometry::Mesh< Real > sMesh;
			std::vector< Point3D< Real > > meshVertices , colors;
			sMesh.read( inFileName , meshVertices , colors );
			sMesh.vertices = advectedVertices;
			sMesh.write( outFileName , meshVertices , colors );
		};
		auto OutputAdvectedMeshes = [&]( int rot , const std::vector< Point3D< Real > > &advectedSourceVertices , const std::vector< Point3D< Real > > &advectedTargetVertices )
		{
			char outFileNames[2][1024];
			if( rot>=0 ) sprintf( outFileNames[0] , "%s_%d.ply" , Out.values[0] , rot ) , sprintf( outFileNames[1] , "%s_%d.ply" , Out.values[1] , rot );
			else         sprintf( outFileNames[0] , "%s.ply"    , Out.values[0]       ) , sprintf( outFileNames[1] , "%s.ply"    , Out.values[1]       );
			OutputAdvectedMesh( In.values[0] , outFileNames[0] , advectedSourceVertices );
			OutputAdvectedMesh( In.values[1] , outFileNames[1] , advectedTargetVertices );
		};

		typename SphericalGeometry::OpticalFlow< double , DivFree >::Parameters params;
		params.subSteps = AdvectionSubSteps.value;
		params.integrationSamples = IntegrationSamples.value;
		params.stepSize = StepSize.value;
		params.smoothingWeight = SmoothingWeight.value;
		params.bandLimitDampening = BandLimitDampening.value;
		params.rescaleFlow = RescaleFlow.set;

		Timer timer;
		typename SphericalGeometry::OpticalFlow< double , DivFree >::Stats advanceState;
		if( !Rotate.set ) // Try all candidate rotations
		{
			std::vector< Point3D< Real > > advectedSourceVertices , advectedTargetVertices;
			int bestRotationIndex = -1;
			Timer timer;
			for( int r=0 ; r<candidateRotations.size() ; r++ )
			{
				typename SphericalGeometry::OpticalFlow< double , DivFree >::Stats _advanceState;
				SphericalGeometry::Mesh< Real > _sMeshes[] = { sMeshes[0] , sMeshes[1] };
				ApplyRotation( _sMeshes[0] , _sMeshes[1] , candidateRotations[r] );

				unsigned int levels = (unsigned int)signals[0].values.size()-2;
				for( unsigned int l=0 ; l<levels ; l++ )
				{
					params.resolution = BaseResolution.value<<l;
					typename SphericalGeometry::OpticalFlow< double , DivFree > opticalFlow( _sMeshes[0] , _sMeshes[1] , signals[0].values[l+1] , signals[1].values[l+1] , params );
					for( int i=0 ; i<SolvesPerLevel.value ; i++ ) opticalFlow.advance();
				}
				{
					typename SphericalGeometry::OpticalFlow< double , DivFree > opticalFlow( _sMeshes[0] , _sMeshes[1] , signals[0].values[levels+1] , signals[1].values[levels+1] , params );
					_advanceState = opticalFlow.getStats();
				}
				bool foundBetter = UseL2Difference.set ? ( _advanceState.l2Difference<advanceState.l2Difference ) : ( _advanceState.l1Difference<advanceState.l1Difference );
				if( r==0 || foundBetter )
				{
					bestRotationIndex = r;
					advanceState = _advanceState;
					advectedSourceVertices = _sMeshes[0].vertices;
					advectedTargetVertices = _sMeshes[1].vertices;
				}
				if( Verbose.set ) printf( "Alignment[%d]: L1-Difference: %g (%g); L2-Difference: %g (%g)\n" , r+1 , _advanceState.l1Difference , _advanceState.l1Difference / _advanceState.l1Norm , _advanceState.l2Difference , _advanceState.l2Difference / _advanceState.l2Norm );
				if( Out.set && OutputAllRotations.set ) OutputAdvectedMeshes( r , _sMeshes[0].vertices , _sMeshes[1].vertices );
			}
			if( Verbose.set )
			{
				printf( "Best rotation: %d / %d\n" , bestRotationIndex+1 , (int)candidateRotations.size() );
				printf( "Optical flow time: %.2f(s)\n" , timer.elapsed() );
			}
			if( Out.set )
				if( OutputAllRotations.set ) OutputAdvectedMeshes( -1 , advectedSourceVertices , advectedTargetVertices );
				else
				{
					OutputAdvectedMesh( In.values[0] , Out.values[0] , advectedSourceVertices );
					OutputAdvectedMesh( In.values[1] , Out.values[1] , advectedTargetVertices );
				}
		}
		else // Try a single rotation
		{
			Timer timer;
			if( Rotate.value>=0 )
			{
				if( Rotate.value<candidateRotations.size() ) ApplyRotation( sMeshes[0] , sMeshes[1] , candidateRotations[ Rotate.value ] );
				else fprintf( stderr , "[ERROR] Rotation index out of bounds: %d <= %d < %d\n" , 0 , Rotate.value , (int)candidateRotations.size() ) , exit( 0 );
			}
			unsigned int levels = (unsigned int)signals[0].values.size()-2;
			for( unsigned int l=0 ; l<levels ; l++ )
			{
				params.resolution = BaseResolution.value<<l;
				typename SphericalGeometry::OpticalFlow< double , DivFree > opticalFlow( sMeshes[0] , sMeshes[1] , signals[0].values[l+1] , signals[1].values[l+1] , params );
				for( int i=0 ; i<SolvesPerLevel.value ; i++ ) opticalFlow.advance();
			}
			{
				typename SphericalGeometry::OpticalFlow< double , DivFree > opticalFlow( sMeshes[0] , sMeshes[1] , signals[0].values[levels+1] , signals[1].values[levels+1] , params );
				_advanceState = opticalFlow.getStats();
			}
			if( Verbose.set )
			{
				printf( "Alignment: L1-Difference: %g (%g); L2-Difference: %g (%g)\n" , advanceState.l1Difference , advanceState.l1Difference / advanceState.l1Norm , advanceState.l2Difference , advanceState.l2Difference / advanceState.l2Norm );
				printf( "Optical flow time: %.2f(s)\n" , timer.elapsed() );
			}
			if( Out.set )
				if( OutputAllRotations.set ) OutputAdvectedMeshes( -1 , sMeshes[0].vertices , sMeshes[1].vertices );
				else
				{
					OutputAdvectedMesh( In.values[0] , Out.values[0] , sMeshes[0].vertices );
					OutputAdvectedMesh( In.values[1] , Out.values[1] , sMeshes[1].vertices );
				}
		}

	}
#ifdef NO_OPEN_GL
#else // !NO_OPEN_GL
	else
	{
		OpticalFlowVisualization< double , DivFree > *&ofv = OpticalFlowVisualizationWrapper< DivFree >::ofv;
		if( Rotate.set )
			if( Rotate.value>=0 && Rotate.value<candidateRotations.size() ) ApplyRotation( sMeshes[0] , sMeshes[1] , candidateRotations[ Rotate.value ] );
			else fprintf( stderr , "[ERROR] Rotation index out of bounds: %d <= %d < %d\n" , 0 , Rotate.value , (int)candidateRotations.size() ) , exit( 0 );

		ofv = new OpticalFlowVisualization< double , DivFree >( sMeshes[0] , colors[0] , &signals[0].values[1] , sMeshes[1] , colors[1] , &signals[1].values[1] , (unsigned int)signals[0].values.size()-2 , BandLimitDampening.value , SmoothingWeight.value , BaseResolution.value , IntegrationSamples.value , StepSize.value , AdvectionSubSteps.value , RescaleFlow.set , SolvesPerLevel.value , FullVerbose.set );
		ofv->showFPS = false;
		ofv->showHelp = false;
		ofv->edgeMode = ShowEdges.set;
		ofv->showColor = ShowColors.set;
		if( InCamera.set ) printf( "Reading camera\n" ) , ofv->readCamera( InCamera.value );
		if( Animate.set ) ofv->animationSteps = -1;
		glutInitDisplayMode( GLUT_DEPTH | GLUT_RGB | GLUT_DOUBLE );
		glutInitWindowPosition( 0 , 0 );
		glutInitWindowSize( ofv->screenWidth , ofv->screenHeight );
		glutInit( &argc , argv );
		char windowName[1024];
		sprintf( windowName , "Optical Flow Viewer ( %s %s )" , In.values[0] , In.values[1] );
		glutCreateWindow( windowName );

		if( glewInit()!=GLEW_OK ) fprintf( stderr , "[ERROR] glewInit failed\n" ) , exit( 0 );
		glutIdleFunc         ( OpticalFlowVisualizationWrapper< DivFree >::Idle );
		glutDisplayFunc      ( OpticalFlowVisualizationWrapper< DivFree >::Display );
		glutReshapeFunc      ( OpticalFlowVisualizationWrapper< DivFree >::Reshape );
		glutMouseFunc        ( OpticalFlowVisualizationWrapper< DivFree >::MouseFunc );
		glutMotionFunc       ( OpticalFlowVisualizationWrapper< DivFree >::MotionFunc );
		glutPassiveMotionFunc( OpticalFlowVisualizationWrapper< DivFree >::PassiveMotionFunc );
		glutKeyboardFunc     ( OpticalFlowVisualizationWrapper< DivFree >::KeyboardFunc );
		glutSpecialFunc      ( OpticalFlowVisualizationWrapper< DivFree >::SpecialFunc );
		glutMainLoop();
		delete ofv;
	}
#endif // NO_OPEN_GL
	return EXIT_SUCCESS;    
}

int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set || !InSignal.set )
	{
		Usage( argv[0] );
		return EXIT_FAILURE;
	}
	omp_set_num_threads( std::max< int >( Threads.value , 1 ) );

	Verbose.set |= FullVerbose.set;

	if( NoDivergenceFree.set ) return _main< double , false >( argc , argv );
	else                       return _main< double , true  >( argc , argv );
}