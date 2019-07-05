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
#include "Misha/MoebiusCentering.h"
#ifdef NO_OPEN_GL
#include "Misha/CMCF.h"
#else // !NO_OPEN_GL
#include "CMCFViewer.inl"
#endif // NO_OPEN_GL

static const int DefaultCenterToInversion = SphericalGeometry::MoebiusCentering< double >::Parameters::PARAMETERS_POINCARE;

cmdLineParameter< char* >
	In( "in" ) ,											// Input mesh
	Out( "out" );											// Output spherical mesh
cmdLineParameter< int >
	Threads( "threads" , omp_get_num_procs() ) ,			// Number of parallelization threads
	EvolutionSteps( "steps" , 100 ) ,						// Number of cMCF evolution steps to perform
	CenterToInversion( "c2i" , DefaultCenterToInversion );	// Center to inversion conversion type
cmdLineParameter< float >
	StepSize( "stepSize" , 0.1f ) ,							// cMCF step size
	CenteringCutOff( "cutOff" , 1e-10f );					// Centering cut-off
cmdLineReadable
	FullVerbose( "fullVerbose" ) , 							// Should very verbose output be provided?
	Verbose( "verbose" );									// Should verbose output be provided?

void Usage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name );
	printf( "\t[--%s <output spherical mesh>]\n" , Out.name );
	printf( "\t[--%s <evolution steps>=%d]\n" , EvolutionSteps.name , EvolutionSteps.value );
	printf( "\t[--%s <flow step size>=%f]\n" , StepSize.name , StepSize.value );
	printf( "\t[--%s <Moebius centering cut-off>=%g]\n" , CenteringCutOff.name , CenteringCutOff.value );
	printf( "\t[--%s <center to inversion type>=%d]\n" , CenterToInversion.name , CenterToInversion.value );
	printf( "\t\t0] Trivial\n" );
	printf( "\t\t1] Golden section search\n" );
	printf( "\t\t2] Poincare\n" );
	printf( "\t[--%s]\n" , FullVerbose.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

cmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&EvolutionSteps ,
	&StepSize ,
	&CenteringCutOff ,
	&CenterToInversion ,
	&FullVerbose ,
	&Verbose ,
	NULL
};


#ifdef NO_OPEN_GL
#else // !NO_OPEN_GL
struct CMCFVisualizationWrapper
{
	static CMCFVisualization< double > cmcfv;
	static void Idle             ( void )                                  { return cmcfv.Idle(); }
	static void Display          ( void )                                  { return cmcfv.Display(); }
	static void KeyboardFunc     ( unsigned char key , int x , int y )     { return cmcfv.KeyboardFunc( key , x , y ); }
	static void SpecialFunc      ( int key , int x , int y )               { return cmcfv.SpecialFunc( key , x , y ); }
	static void Reshape          ( int w , int h )                         { return cmcfv.Reshape( w , h ); }
	static void MouseFunc        ( int button , int state , int x , int y ){ return cmcfv.MouseFunc( button , state , x , y ); }
	static void MotionFunc       ( int x , int y )                         { return cmcfv.MotionFunc( x , y ); }
	static void PassiveMotionFunc( int x , int y )                         { return cmcfv.PassiveMotionFunc( x , y ); }
};
CMCFVisualization< double > CMCFVisualizationWrapper::cmcfv;
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
	Verbose.set |= FullVerbose.set;

#ifdef NO_OPEN_GL
#else // !NO_OPEN_GL
	if( Out.set )
#endif // NO_OPEN_GL
	{
		Timer timer;

		typename SphericalGeometry::CMCF< double >::Stats originalStats , currentStats;
		typename SphericalGeometry::CMCF< double >::Parameters parameters;
		parameters.stepSize = StepSize.value;

		std::vector< Point3D< double > > vertices , colors;
		std::vector< TriangleIndex > triangles;
		bool hasColor;
		{
			std::vector< PlyColorVertex< float > > _vertices;
			int fileType;
			bool readFlags[ PlyColorVertex< float >::ReadComponents ];
			PlyReadTriangles( In.value , _vertices , triangles , PlyColorVertex< float >::ReadProperties , readFlags , PlyColorVertex< float >::ReadComponents , fileType );
			hasColor = ( readFlags[3] || readFlags[6] ) && (readFlags[4] && readFlags[7] ) && ( readFlags[5] && readFlags[8] );

			vertices.resize( _vertices.size() );
#pragma omp parallel for
			for( int i=0 ; i<_vertices.size() ; i++ ) vertices[i] = Point3D< double >( _vertices[i].point );
			colors.resize( _vertices.size() );
			if( hasColor )
#pragma omp parallel for
				for( int i=0 ; i<_vertices.size() ; i++ ) colors[i] = Point3D< double >( _vertices[i].color );
			else
			{
				for( int i=0 ; i<triangles.size() ; i++ )
				{
					Point3D< double > v[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
					Point3D< double > n = Point3D< double >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
					for( int j=0 ; j<3 ; j++ ) colors[ triangles[i][j] ] += n;
				}
				auto NormalColor = []( Point3D< double > n )
				{
					n /= -Length( n );
					Point3D< double > c = ( -n + Point3D< double >( 1. , 1. , 1. ) ) * 128.;
					for( int d=0 ; d<3 ; d++ )
					{
						if( c[d]>255. ) c[d] = 255.;
						if( c[d]<  0. ) c[d] =   0.;
					}
					return c;
				};

				for( int i=0 ; i<colors.size() ; i++ ) colors[i] = NormalColor( colors[i] );
			}
		}

		SphericalGeometry::Mesh< double > sMesh;

		typename SphericalGeometry::CMCF< double > cmcf( vertices , triangles , sMesh , parameters );

		originalStats = cmcf.getStats();
		for( int i=0 ; i<EvolutionSteps.value ; i++ )
		{
			cmcf.advance();
			if( FullVerbose.set )
			{
				currentStats = cmcf.getStats();
				printf( "cMCF[%d] Delta / log(QC) / Radial: %.4f / %.4f / %.4f    \r" , i+1 , currentStats.deformationScale , log(currentStats.quasiConformalRatio) , currentStats.radialDeviation );
			}
		}
		if( FullVerbose.set ) printf( "\n" );
		sMesh.setMasses( vertices );

#pragma omp parallel for
		for( int i=0 ; i<sMesh.vertices.size() ; i++ ) sMesh.vertices[i] /= Length( sMesh.vertices[i] );

		if( CenterToInversion.value==SphericalGeometry::MoebiusCentering< double >::Parameters::PARAMETERS_CENTER_TO_INVERSION || CenterToInversion.value==SphericalGeometry::MoebiusCentering< double >::Parameters::PARAMETERS_GOLDEN_SECTION_SEARCH || CenterToInversion.value==SphericalGeometry::MoebiusCentering< double >::Parameters::PARAMETERS_POINCARE )
		{
			typename SphericalGeometry::MoebiusCentering< double >::Parameters parameters( CenterToInversion.value );
			typename SphericalGeometry::MoebiusCentering< double > mc( sMesh , parameters );

			static const int MAX_ITERATIONS = 50;
			int iters;

			for( iters=0 ; iters<MAX_ITERATIONS ; iters++ )
			{
				mc.advance();
				typename SphericalGeometry::MoebiusCentering< double >::Stats stats = mc.getStats();
				if( FullVerbose.set ) printf( "%d] %g (%g %g %g)       \r" , iters+1 , Length( stats.center ) , stats.center[0] , stats.center[1] , stats.center[2] );
				if( stats.center.squareNorm()<=CenteringCutOff.value*CenteringCutOff.value ) break;
			}
			if( FullVerbose.set ) printf( "\n" );
			if( iters==MAX_ITERATIONS && CenterToInversion.value!=SphericalGeometry::MoebiusCentering< double >::Parameters::PARAMETERS_GOLDEN_SECTION_SEARCH )
			{
				fprintf( stderr , "[WARNING] Failed to meet centering threshold after %d iterations, trying golden-section search\n" , MAX_ITERATIONS );
				typename SphericalGeometry::MoebiusCentering< double >::Parameters parameters( SphericalGeometry::MoebiusCentering< double >::Parameters::PARAMETERS_GOLDEN_SECTION_SEARCH );
				typename SphericalGeometry::MoebiusCentering< double > mc( sMesh , parameters );

				for( iters=0 ; iters<MAX_ITERATIONS ; iters++ )
				{
					mc.advance();
					typename SphericalGeometry::MoebiusCentering< double >::Stats stats = mc.getStats();
					if( FullVerbose.set ) printf( "%d] %g (%g %g %g)       \r" , iters+1 , Length( stats.center ) , stats.center[0] , stats.center[1] , stats.center[2] );
					if( stats.center.squareNorm()<=CenteringCutOff.value*CenteringCutOff.value ) break;
				}
				if( FullVerbose.set ) printf( "\n" );
			}
			if( iters==MAX_ITERATIONS ) fprintf( stderr , "[WARNING] Failed to meet centering threshold after %d iterations\n" , MAX_ITERATIONS );
		}

		if( Out.set ) sMesh.write( Out.value , vertices , colors );

		if( Verbose.set ) printf( "cMCF: %.2f(s)\n" , timer.elapsed() );
	}
#ifdef NO_OPEN_GL
#else // !NO_OPEN_GL
	else
	{
		CMCFVisualization< double >& cmcfv = CMCFVisualizationWrapper::cmcfv;
		if( EvolutionSteps.set ) cmcfv.animationSteps = EvolutionSteps.value;
		cmcfv.init( In.value , StepSize.value );

		glutInitDisplayMode( GLUT_DEPTH | GLUT_RGB | GLUT_DOUBLE );
		glutInitWindowSize( cmcfv.screenWidth , cmcfv.screenHeight );
		glutInit( &argc , argv );
		char windowName[1024];
		sprintf( windowName , "cMCF Viewer ( %s )" , In.value );
		glutCreateWindow( windowName );

		if( glewInit()!=GLEW_OK ) fprintf( stderr , "[ERROR] glewInit failed\n" ) , exit( 0 );
		glutIdleFunc         ( CMCFVisualizationWrapper::Idle );
		glutDisplayFunc      ( CMCFVisualizationWrapper::Display );
		glutReshapeFunc      ( CMCFVisualizationWrapper::Reshape );
		glutMouseFunc        ( CMCFVisualizationWrapper::MouseFunc );
		glutMotionFunc       ( CMCFVisualizationWrapper::MotionFunc );
		glutPassiveMotionFunc( CMCFVisualizationWrapper::PassiveMotionFunc );
		glutKeyboardFunc     ( CMCFVisualizationWrapper::KeyboardFunc );
		glutSpecialFunc      ( CMCFVisualizationWrapper::SpecialFunc );

		glutMainLoop();
	}
#endif // NO_OPEN_GL
	return EXIT_SUCCESS;    
}