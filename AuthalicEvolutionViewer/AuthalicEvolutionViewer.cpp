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

#include <cstdlib>
#include <vector>
#include <omp.h>

#include "Misha/Geometry.h"
#include "Misha/Algebra.h"
#include "Misha/Ply.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Timer.h"
#ifdef NO_OPEN_GL
#include "Misha/AuthalicEvolution.h"
#else // !NO_OPEN_GL
#include "AuthalicEvolutionViewer.inl"
#endif // NO_OPEN_GL

cmdLineParameter< char* >
	In( "in" ) ,									// Input spherical parameterization
	Out( "out" );									// Output authalic spherical parameterization
cmdLineParameter< int >
	Resolution( "res" , 256 ) ,						// Spherical grid resolution
	SubSteps( "subSteps" , 1 ) ,					// Number of geodesic steps to take per advection
	EvolutionSteps( "steps" , 100 ) ,				// Number of authalic evolution steps to perform
	Threads( "threads" , omp_get_num_procs() ) ,	// Number of threads to use for parallelization
	CenteringIterations( "cIters" , 3 );			// Number of centering iterations to perform after each evolution step
cmdLineParameter< float >
	StepSize( "stepSize" , 0.005f ) ,				// The size of an advection step
	Smoothness( "smooth" , 2.5e-3f );				// That diffusion factor for smoothing the signal prior to computing the gradient
cmdLineReadable
	UseGoldenSectionSearch( "useGSS" ) ,			// Should a golden section search be used for centering?
	ShowMeshOnly( "meshOnly" ) ,					// Should only the mesh be shown?
	Verbose( "verbose" );							// Should verbose output be provided?

void Usage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input spherical mesh>\n" , In.name );
	printf( "\t[--%s <output spherical mesh>]\n" , Out.name );
	printf( "\t[--%s <flow step size>=%f]\n" , StepSize.name , StepSize.value );
	printf( "\t[--%s <smoothness>=%f]\n" , Smoothness.name , Smoothness.value );
	printf( "\t[--%s <Moebius centering iterations>=%d]\n" , CenteringIterations.name , CenteringIterations.value );
	printf( "\t[--%s <spherical grid resolution>=%d]\n" , Resolution.name , Resolution.value );
	printf( "\t[--%s <advection sub-steps>=%d]\n" , SubSteps.name , SubSteps.value );
	printf( "\t[--%s <evolution steps>=%d]\n" , EvolutionSteps.name , EvolutionSteps.value );
	printf( "\t[--%s]\n" , UseGoldenSectionSearch.name );
	printf( "\t[--%s]\n" , ShowMeshOnly.name );
	printf( "\t[--%s]\n" , Verbose.name );
}
cmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&StepSize ,
	&Smoothness ,
	&Threads ,
	&EvolutionSteps ,
	&Resolution ,
	&SubSteps ,
	&CenteringIterations ,
	&UseGoldenSectionSearch ,
	&ShowMeshOnly ,
	&Verbose ,
	NULL
};

#ifdef NO_OPEN_GL
#else // !NO_OPEN_GL
struct AuthalicEvolutionVisualizationWrapper
{
	static AuthalicEvolutionVisualization< double > aev;
	static void Idle             ( void )                                  { return aev.Idle(); }
	static void Display          ( void )                                  { return aev.Display(); }
	static void KeyboardFunc     ( unsigned char key , int x , int y )     { return aev.KeyboardFunc( key , x , y ); }
	static void SpecialFunc      ( int key , int x , int y )               { return aev.SpecialFunc( key , x , y ); }
	static void Reshape          ( int w , int h )                         { return aev.Reshape( w , h ); }
	static void MouseFunc        ( int button , int state , int x , int y ){ return aev.MouseFunc( button , state , x , y ); }
	static void MotionFunc       ( int x , int y )                         { return aev.MotionFunc( x , y ); }
	static void PassiveMotionFunc( int x , int y )                         { return aev.PassiveMotionFunc( x , y ); }
};
AuthalicEvolutionVisualization< double > AuthalicEvolutionVisualizationWrapper::aev;
#endif // NO_OPEN_GL

int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		Usage( argv[0] );
		return EXIT_FAILURE;
	}
	omp_set_num_threads( std::max< int >( Threads.value , 1 ) );

#ifdef NO_OPEN_GL
#else // !NO_OPEN_GL
	if( Out.set )
#endif // NO_OPEN_GL
	{
		Timer timer;

		double originalFlippedTriangleFraction , currentFlippedTriangleFraction;
		typename SphericalGeometry::AuthalicEvolution< double >::Stats originalStats , currentStats;

		typename SphericalGeometry::AuthalicEvolution< double >::Parameters parameters;
		parameters.resolution = Resolution.value;
		parameters.smoothness = Smoothness.value;
		parameters.stepSize = StepSize.value;
		parameters.subSteps = SubSteps.value;

		SphericalGeometry::Mesh< double > sMesh;
		std::vector< Point3D< double > > meshVertices , colors;

		bool hasColor = sMesh.read( In.value , meshVertices , colors );

		typename SphericalGeometry::AuthalicEvolution< double > authalicEvolution( sMesh , meshVertices , parameters );

		originalStats = authalicEvolution.getStats();
		originalFlippedTriangleFraction = authalicEvolution.flippedTriangleFraction();
		for( int i=0 ; i<EvolutionSteps.value ; i++ ) authalicEvolution.advance();
		currentStats = authalicEvolution.getStats();
		currentFlippedTriangleFraction = authalicEvolution.flippedTriangleFraction();

		if( Out.set )
			if( hasColor ) sMesh.write( Out.value , meshVertices , colors );
			else           sMesh.write( Out.value , meshVertices );

		printf( "Authalic evolution: %.2f(s)\n" , timer.elapsed() );
		printf( "\tStats: [ %g , %g ] +/- %g -> [ %g , %g ] +/- %g\n" , originalStats.min , originalStats.max , originalStats.dev , currentStats.min , currentStats.max , currentStats.dev );
		printf( "\tFlipped triangles: %g -> %g\n" , originalFlippedTriangleFraction , currentFlippedTriangleFraction );
	}
#ifdef NO_OPEN_GL
#else // !NO_OPEN_GL
	else
	{
		AuthalicEvolutionVisualization< double >& aev = AuthalicEvolutionVisualizationWrapper::aev;
		if( EvolutionSteps.set ) aev.animationSteps = EvolutionSteps.value;
		aev.verbose = Verbose.set;
		aev.parameters.centeringIterations = CenteringIterations.value;
		aev.parameters.useGoldenSectionSearch = UseGoldenSectionSearch.set;
		aev.init( In.value , Smoothness.value , StepSize.value , Resolution.value , SubSteps.value , ShowMeshOnly.set );

		glutInitDisplayMode( GLUT_DEPTH | GLUT_RGB | GLUT_DOUBLE );
		glutInitWindowSize( aev.screenWidth , aev.screenHeight );
		glutInit( &argc , argv );
		char windowName[1024];
		sprintf( windowName , "Authalic Evolution Viewer ( %s )" , In.value );
		glutCreateWindow( windowName );

		if( glewInit()!=GLEW_OK ) fprintf( stderr , "[ERROR] glewInit failed\n" ) , exit( 0 );
		glutIdleFunc         ( AuthalicEvolutionVisualizationWrapper::Idle );
		glutDisplayFunc      ( AuthalicEvolutionVisualizationWrapper::Display );
		glutReshapeFunc      ( AuthalicEvolutionVisualizationWrapper::Reshape );
		glutMouseFunc        ( AuthalicEvolutionVisualizationWrapper::MouseFunc );
		glutMotionFunc       ( AuthalicEvolutionVisualizationWrapper::MotionFunc );
		glutPassiveMotionFunc( AuthalicEvolutionVisualizationWrapper::PassiveMotionFunc );
		glutKeyboardFunc     ( AuthalicEvolutionVisualizationWrapper::KeyboardFunc );
		glutSpecialFunc      ( AuthalicEvolutionVisualizationWrapper::SpecialFunc );

		glutMainLoop();
	}
#endif // NO_OPEN_GL
	return EXIT_SUCCESS;    
}