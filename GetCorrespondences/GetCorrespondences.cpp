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

#include <cstdlib>
#include <vector>
#include <omp.h>

#include "Misha/CmdLineParser.h"
#include "Misha/Timer.h"
#include "Misha/SourceToTargetCorrespondence.h"

cmdLineParameterArray< char* , 2 >
	In( "in" ) ,				// Input source and target spherical meshes
	Out( "out" );				// Output source-to-target and target-to-source correspondences
cmdLineParameter< int >
	Resolution( "res" , 256 );	// Spherical mesh resolution
cmdLineReadable
	Verbose( "verbose" );		// Should verbose output be provided?

void Usage( const char* ex ) 
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input source/target spherical meshes>\n" , In.name );
	printf( "\t[--%s <output source-to-target/target-to-source mapping>]\n" , Out.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

cmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&Verbose ,
	NULL
};

int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		Usage( argv[0] );
		return EXIT_FAILURE;
	}

	Timer timer;
	SourceTargetCorrespondences< double > correspondences( In.values[0] , In.values[1] , Resolution.value );
	if( Verbose.set ) printf( "Got correspondences: %.2f(s)\n" , timer.elapsed() );

	if( Out.set )
	{
		SourceTargetCorrespondences< double >::WriteCorrespondences( Out.values[0] , correspondences.correspondences[0] );
		SourceTargetCorrespondences< double >::WriteCorrespondences( Out.values[1] , correspondences.correspondences[1] );
	}

	return EXIT_SUCCESS;    
}