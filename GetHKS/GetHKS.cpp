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
#include <algorithm>
#include <omp.h>
#include "Misha/CmdLineParser.h"
#include "Misha/Ply.h"
#include "Misha/Signals.h"

const float DefaultTimeScales[] = { 1.f , 0.5f , 0.25f , 0.125f , 0.0625f , 0.01f };
const unsigned int DefaulTimeScalesCount = sizeof( DefaultTimeScales ) / sizeof( float );

cmdLineParameter< char * >
	InMesh( "in" ) ,															// Input mesh
	InSpectrum( "spectrum" ) ,													// Input spectrum
	Out( "out" );																// Output HKS functions
cmdLineParameters< float >
	HKSTimeScales( "hksTimes" , DefaultTimeScales , DefaulTimeScalesCount );	// Time-scales for the HKS functions
cmdLineParameter< float >
	SpectrumOffset( "off" , 100.f );											// Offset used for the "shift" part of "invert-and-shift"
cmdLineParameter< int >
	SpectrumDimension( "dim" , 200 );											// Maximum dimension of the spectrum
cmdLineReadable
	SpectrumLump( "lump" ) ,													// Should the mass matrix be lumped?
	Verbose( "verbose" );														// Should verbose output be provided?


void Usage( const char* ex ) 
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , InMesh.name );
	printf( "\t[--%s <input spectrum>]\n" , InSpectrum.name );
	printf( "\t[--%s <output signal>]\n" , Out.name );
	printf( "\t[--%s <# hks time scales> <hks time scales>]\n" , HKSTimeScales.name );
	printf( "\t[--%s <maximum specturm dimensions>=%d]\n" , SpectrumDimension.name , SpectrumDimension.value );
	printf( "\t[--%s <spectral offset>=%f]\n" , SpectrumOffset.name , SpectrumOffset.value );
	printf( "\t[--%s]\n" , SpectrumLump.name );
	printf( "\t[--%s]\n" , Verbose.name );
}
cmdLineReadable* params[] = { &InMesh , &Out , &InSpectrum , &HKSTimeScales , &Verbose , &SpectrumDimension , &SpectrumOffset , &SpectrumLump , NULL };

template< typename Real >
void _main( int argc , char *argv[] )
{
	Timer tmr;
	Signal< Real > signal;
	signal.values.resize( HKSTimeScales.count );

	Spectrum< Real > spectrum;

	if( InSpectrum.set ) spectrum.read( InSpectrum.value );
	else if( InMesh.set )
	{
		Timer tmr;
		std::vector< TriangleIndex > triangles;
		std::vector< Point3D< Real > > vertices;

		//////////////////////
		// Read in the data //
		{
			int file_type;
			std::vector< PlyVertex< float > > _vertices;
			PlyReadTriangles( InMesh.value , _vertices ,  triangles , PlyVertex< float >::ReadProperties , NULL , PlyVertex< float >::ReadComponents , file_type );
			vertices.resize( _vertices.size() );
			for( int i=0 ; i<_vertices.size() ; i++ ) vertices[i] = Point3D< Real >( _vertices[i].point );
		}
		// Read in the data //
		//////////////////////

		spectrum.set( vertices , triangles , SpectrumDimension.value , SpectrumOffset.value , SpectrumLump.set );
		if( Verbose.set ) printf( "Got spectrum: %.2f(s)\n" , tmr.elapsed() );
	}
	else fprintf( stderr , "[ERROR] Either an input mesh or a spectrum needs to be provided\n" ) , exit( 0 );
	tmr.reset();
	for( int t=0 ; t<HKSTimeScales.count ; t++ ) Signal< Real >::SetHKS( spectrum , HKSTimeScales.values[t] , signal.values[t] );
	if( Verbose.set ) printf( "Got HKS: %.2f(s)\n" , tmr.elapsed() );
	if( Out.set ) signal.write( Out.value );
}

int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !InMesh.set && !InSpectrum.set )
	{
		Usage( argv[0] );
		return EXIT_FAILURE;
	}
	if( !HKSTimeScales.count )
	{
		fprintf( stderr , "[ERROR] No time scales specified\n" );
		Usage( argv[0] );
		return EXIT_FAILURE;
	}
	_main< double >( argc , argv );
	return EXIT_SUCCESS;
}