/*
Copyright (c) 2018, Michael Kazhdan
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

#ifndef VISUALIZATION_INCLUDED
#define VISUALIZATION_INCLUDED

#include <algorithm>
#include <GL/glew.h>
#include <vector>
#include <sys/timeb.h>
#include <algorithm>
#include "Misha/Image.h"
#include "Misha/Array.h"


namespace Visualization
{
	static const int KEY_UPARROW    = 101;
	static const int KEY_DOWNARROW	= 103;
	static const int KEY_LEFTARROW	= 100;
	static const int KEY_RIGHTARROW	= 102;
	static const int KEY_PGUP		= 104;
	static const int KEY_PGDN		= 105;
	static const int KEY_CTRL_C     =   3;
	static const int KEY_BACK_SPACE =   8;
	static const int KEY_ENTER      =  13;
	static const int KEY_ESC        =  27;

	double Time( void )
	{
#ifdef WIN32
		struct _timeb t;
		_ftime(&t);
		return double(t.time)+double(t.millitm)/1000.0;
#else // WIN32
		struct timeval t;
		gettimeofday(&t,NULL);
		return t.tv_sec+(double)t.tv_usec/1000000;
#endif // WIN32
	}

	struct KeyboardCallBack
	{
		char key;
		char prompt[1024];
		char description[1024];
		void (*callBackFunction)( struct Viewable* , const char* );
		struct Viewable* viewable;
		KeyboardCallBack( struct Viewable* viewable , char key , const char* description , void (*callBackFunction)( Viewable* , const char* ) );
		KeyboardCallBack( struct Viewable* viewable , char key , const char* description , const char* prompt , void ( *callBackFunction )( Viewable* , const char* ) );
	};

	struct Viewable
	{
	protected:
		const double _MIN_FPS_TIME = 0.5;
		double _lastFPSTime;
		int _lastFPSCount;
		double _fps;
		int _currentFrame , _totalFrames;
		bool _exitAfterSnapshot , _exitAfterVideo , _fullScreen;
		int _screenWidth , _screenHeight;
	public:
		int screenWidth , screenHeight;
		void *font , *promptFont;
		int fontHeight , promptFontHeight;
		bool showHelp , showInfo , showFPS;
		void (*promptCallBack)( Viewable* , const char* );
		char promptString[1024];
		int promptLength;
		char* snapshotName;
		char* videoHeader;
		bool flushImage;

		std::vector< KeyboardCallBack > callBacks;
		std::vector< char* > info;
		Viewable( void );
		void setSnapshot( const char* sName , bool exitAfterSnapshot=true );
		void setVideo( const char* sName , int frames , bool exitAfterVideo=true );
		virtual void display( void ) {}
		virtual void idle( void ) {}
		virtual void keyboardFunc( unsigned char key , int x , int y ) {}
		virtual void specialFunc( int key, int x, int y ) {}
		virtual void mouseFunc( int button , int state , int x , int y ) {}
		virtual void motionFunc( int x , int y ) {}
		virtual void passiveMotionFunc( int x , int y ) {}

		void Idle             ( void );
		void KeyboardFunc     ( unsigned char key , int x , int y );
		void SpecialFunc      ( int key, int x, int y );
		void Display          ( void );
		void Reshape          ( int w , int h );
		void MouseFunc        ( int button , int state , int x , int y );
		void MotionFunc       ( int x , int y );
		void PassiveMotionFunc( int x , int y );

		static void             QuitCallBack( Viewable*   , const char* ){ exit( 0 ); }
		static void        ToggleFPSCallBack( Viewable* v , const char* ){ v->showFPS  = !v->showFPS ; }
		static void       ToggleHelpCallBack( Viewable* v , const char* ){ v->showHelp = !v->showHelp; }
		static void       ToggleInfoCallBack( Viewable* v , const char* ){ v->showInfo = !v->showInfo; }
		static void ToggleFullScreenCallBack( Viewable* v , const char* )
		{
			v->_fullScreen = !v->_fullScreen;
			if( !v->_fullScreen ) glutReshapeWindow( v->_screenWidth , v->_screenHeight );
			else                  glutFullScreen();
		}
		static void SetFrameBufferCallBack( Viewable* v , const char* prompt );

		static void WriteLeftString( int x , int y , void* font , const char* format , ... );
		static int StringWidth( void* font , const char* format , ... );
		void writeLeftString( int x , int y , const char* format , ... ) const;
		void writeRightString( int x , int y , const char* format , ... ) const;
		int stringWidth( const char* format , ... ) const;

		void saveFrameBuffer( const char* fileName , int whichBuffer=GL_BACK );
	};

	struct Viewer
	{
		static Viewable* viewable;
		static void Run( Viewable* viewable , int argc , char* argv[] , const char* windowName="" );
		static void Idle             ( void );
		static void KeyboardFunc     ( unsigned char key , int x , int y );
		static void SpecialFunc      ( int key, int x, int y );
		static void Display          ( void );
		static void Reshape          ( int w , int h );
		static void MouseFunc        ( int button , int state , int x , int y );
		static void MotionFunc       ( int x , int y );
		static void PassiveMotionFunc( int x , int y );
	};

	//////////////////////
	// KeyboardCallBack //
	//////////////////////
	KeyboardCallBack::KeyboardCallBack( Viewable* viewable , char key , const char* description , void (*callBackFunction)( Viewable* , const char* ) )
	{
		this->viewable = viewable;
		this->key = key;
		strcpy( this->description , description );
		prompt[0] = 0;
		this->callBackFunction = callBackFunction;
	}

	KeyboardCallBack::KeyboardCallBack( Viewable* viewable , char key , const char* description , const char* prompt , void ( *callBackFunction )( Viewable* , const char* ) )
	{
		this->viewable = viewable;
		this->key = key;
		strcpy( this->description , description );
		strcpy( this->prompt , prompt );
		this->callBackFunction = callBackFunction;
	}

	//////////////
	// Viewable //
	//////////////
	void Viewable::Reshape( int w , int h )
	{
		screenWidth = w , screenHeight = h;
		if( !_fullScreen ) _screenWidth = w , _screenHeight = h;
		glViewport( 0 , 0 , screenWidth , screenHeight );
	}
	void Viewable::MouseFunc( int button , int state , int x , int y ){ mouseFunc( button , state , x , y ); }
	void Viewable::MotionFunc( int x , int y ){ motionFunc( x , y );}
	void Viewable::PassiveMotionFunc( int x , int y ){ passiveMotionFunc( x , y );}
	void Viewable::Idle( void )
	{
		if( snapshotName )
		{
			if( flushImage )
			{
				flushImage = false;
				glutPostRedisplay();
				return;
			}
			else
			{
				saveFrameBuffer( snapshotName , GL_FRONT );
				delete[] snapshotName;
				snapshotName = NULL;
				if( _exitAfterSnapshot ) exit( 0 );
			}
		}
		else if( videoHeader && _currentFrame<_totalFrames )
		{
			char snapshotName[512];
			sprintf( snapshotName , "%s.%04d.jpg" , videoHeader , _currentFrame );
			saveFrameBuffer( snapshotName , GL_FRONT );
			printf( "Writing frame %d / %d to %s\r" , _currentFrame , _totalFrames , snapshotName );
			_currentFrame++;
			if( _currentFrame==_totalFrames )
			{
				printf( "\n" );
				if( _exitAfterVideo ) exit( 0 );
			}
		}
		idle();
	}
	void Viewable::KeyboardFunc( unsigned char key , int x , int y )
	{
		if( promptCallBack )
		{
			size_t len = strlen( promptString );
			if( key==KEY_BACK_SPACE )
			{
				if( len>promptLength ) promptString[len-1] = 0;
			}
			else if( key==KEY_ENTER )
			{
				promptCallBack( this , promptString+promptLength );
				promptString[0] = 0;
				promptLength = 0;
				promptCallBack = NULL;
			}
			else if( key==KEY_CTRL_C )
			{
				promptString[0] = 0;
				promptLength = 0;
				promptCallBack = NULL;
			}
			else if( key>=32 && key<=126 ) // ' ' to '~'
			{
				promptString[ len ] = key;
				promptString[ len+1 ] = 0;
			}
			glutPostRedisplay();
			return;
		}
		switch( key )
		{
		case KEY_CTRL_C:
			exit( 0 );
			break;
		default:
			for( int i=0 ; i<callBacks.size() ; i++ ) if( callBacks[i].key==key )
			{
				if( strlen( callBacks[i].prompt ) )
				{
					sprintf( promptString , "%s: " , callBacks[i].prompt );
					promptLength = int( strlen( promptString ) );
					promptCallBack = callBacks[i].callBackFunction;
				}
				else (*callBacks[i].callBackFunction)( this , NULL );
				break;
			}
		}
		keyboardFunc( key , x , y );
		glutPostRedisplay();
	}

	void Viewable::SpecialFunc( int key , int x , int y ){ specialFunc( key , x , y );}
	void Viewable::Display( void )
	{
		glClearColor( 1 , 1 , 1 , 1 );
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

		display();

		_lastFPSCount++;
		double t = Time();
		if( t-_lastFPSTime > _MIN_FPS_TIME )
		{
			_fps = (double)_lastFPSCount / (t-_lastFPSTime);
			_lastFPSCount = 0;
			_lastFPSTime = t;
		}
		if( showFPS ) writeRightString( 5 , screenHeight - fontHeight - 5 , "%d x %d @ %.2f" , screenWidth , screenHeight , _fps );

		glDisable( GL_LIGHTING );
		int offset = fontHeight/2;
		if( showHelp )
		{
			{
				GLint vp[4];
				glGetIntegerv( GL_VIEWPORT , vp );

				glMatrixMode( GL_PROJECTION );
				glPushMatrix();
				glLoadIdentity();
				glOrtho( vp[0] , vp[2] , vp[1] , vp[3] , 0 , 1 );

				glMatrixMode( GL_MODELVIEW );
				glPushMatrix();
				glLoadIdentity();

				int x=0 , y = offset;
				for( int i=0 ; i<callBacks.size() ; i++ ) if( strlen( callBacks[i].description ) ) x = std::max< int >( x , stringWidth( "\'%c\': %s" , callBacks[i].key , callBacks[i].description ) ) , y += fontHeight + offset;

				GLint srcAlpha , dstAlpha;
				glGetIntegerv( GL_BLEND_SRC_ALPHA , &srcAlpha );
				glGetIntegerv( GL_BLEND_DST_ALPHA , &dstAlpha );

				glEnable( GL_BLEND );
				glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
				glBegin( GL_QUADS );
				glColor4f( 1.f , 1.f , 1.f , 0.5f );
				glVertex2i( screenWidth-5 , 0 ) , glVertex2i( screenWidth-(x+15) , 0 ) , glVertex2i( screenWidth-(x+15) , y ) , glVertex2i( screenWidth-5 , y );
				glEnd();
				glBlendFunc( srcAlpha , dstAlpha );
				glDisable( GL_BLEND );
				glDisable( GL_DEPTH_TEST );
				glLineWidth( 2.f );
				glBegin( GL_LINE_LOOP );
				glColor4f( 0.f , 0.f , 0.f , 1.f );
				glVertex2i( screenWidth-5 , 0 ) , glVertex2i( screenWidth-(x+15) , 0 ) , glVertex2i( screenWidth-(x+15) , y ) , glVertex2i( screenWidth-5 , y );
				glEnd();

				glMatrixMode( GL_PROJECTION );
				glPopMatrix();

				glMatrixMode( GL_MODELVIEW );
				glPopMatrix();
			}

			{
				int y = offset , width = 0;

				for( int i=0 ; i<callBacks.size() ; i++ ) if( strlen( callBacks[i].description ) )
					width = std::max< int >( width , stringWidth( "\'%c\': %s" , callBacks[i].key , callBacks[i].description ) );
				for( int i=0 ; i<callBacks.size() ; i++ ) if( strlen( callBacks[i].description ) )
					writeLeftString( screenWidth - 10 - width , y , "\'%c\': %s" , callBacks[i].key , callBacks[i].description ) , y += fontHeight + offset;
			}
		}
		if( showInfo && info.size() )
		{
			{
				GLint vp[4];
				glGetIntegerv( GL_VIEWPORT , vp );

				glMatrixMode( GL_MODELVIEW );
				glPushMatrix();
				glLoadIdentity();
				glMatrixMode( GL_PROJECTION );
				glPushMatrix();
				glLoadIdentity();
				glOrtho( vp[0] , vp[2] , vp[1] , vp[3] , 0 , 1 );
				int x=0 , y = offset;
				for( int i=0 ; i<info.size() ; i++ ) if( strlen( info[i] ) ) x = std::max< int >( x , glutBitmapLength( font , (unsigned char*) info[i] ) ) , y += fontHeight + offset;
				glEnable( GL_BLEND );
				GLint srcAlpha , dstAlpha;
				glGetIntegerv( GL_BLEND_SRC_ALPHA , &srcAlpha );
				glGetIntegerv( GL_BLEND_DST_ALPHA , &dstAlpha );
				glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
				glBegin( GL_QUADS );
				glColor4f( 1.f , 1.f , 1.f , 0.5f );
				glVertex2i( 5 , 0 ) , glVertex2i( x+15 , 0 ) , glVertex2i( x+15 , y ) , glVertex2i( 5 , y );
				glEnd();
				glBlendFunc( srcAlpha , dstAlpha );
				glDisable( GL_BLEND );
				glDisable( GL_DEPTH_TEST );
				glLineWidth( 2.f );
				glBegin( GL_LINE_LOOP );
				glColor4f( 0.f , 0.f , 0.f , 1.f );
				glVertex2i( 5 , 0 ) , glVertex2i( x+15 , 0 ) , glVertex2i( x+15 , y ) , glVertex2i( 5 , y );
				glEnd();

				glMatrixMode( GL_PROJECTION );
				glPopMatrix();

				glMatrixMode( GL_MODELVIEW );
				glPopMatrix();
			}
			{
				int y = offset;
				for( int i=0 ; i<info.size() ; i++ ) if( strlen( info[i] ) )
					writeLeftString( 10 , y , "%s" , info[i] ) , y += fontHeight + offset;
			}
		}
		if( strlen( promptString ) )
		{
			void* _font = font;
			int _fontHeight = fontHeight;
			font = promptFont;
			fontHeight = promptFontHeight;

			int sw = StringWidth ( font , promptString );
			glColor4f( 1.f , 1.f , 1.f , 0.5 );
			glEnable( GL_BLEND );
			GLint srcAlpha , dstAlpha;
			glGetIntegerv( GL_BLEND_SRC_ALPHA , &srcAlpha );
			glGetIntegerv( GL_BLEND_DST_ALPHA , &dstAlpha );
			glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
			glBegin( GL_QUADS );
			{
				glVertex2i(     0 , screenHeight              );
				glVertex2i( sw+20 , screenHeight              );
				glVertex2i( sw+20 , screenHeight-fontHeight*2 );
				glVertex2i(     0 , screenHeight-fontHeight*2 );
			}
			glEnd();
			glBlendFunc( srcAlpha , dstAlpha );
			glDisable( GL_BLEND );
			glColor4f( 0.f , 0.f , 0.f , 1.f );
			glLineWidth( 2.f );
			glBegin( GL_LINE_LOOP );
			{
				glVertex2i(     0 , screenHeight              );
				glVertex2i( sw+20 , screenHeight              );
				glVertex2i( sw+20 , screenHeight-fontHeight*2 );
				glVertex2i(     0 , screenHeight-fontHeight*2 );
			}
			glEnd();
			writeLeftString( 10 , screenHeight-fontHeight-fontHeight/2 , promptString );
			font = _font;
			fontHeight = _fontHeight;
		}
		glutSwapBuffers();
	}

	void Viewable::WriteLeftString( int x , int y , void* font , const char* format , ... )
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}

		GLint vp[4];
		glGetIntegerv( GL_VIEWPORT , vp );

		glMatrixMode( GL_PROJECTION );
		glPushMatrix();
		glLoadIdentity();
		glOrtho( vp[0] , vp[2] , vp[1] , vp[3] , 0 , 1 );

		glMatrixMode( GL_MODELVIEW );
		glPushMatrix();
		glLoadIdentity();

		GLint matrixMode;
		glGetIntegerv( GL_MATRIX_MODE , &matrixMode );
		int depth = glIsEnabled( GL_DEPTH_TEST );
		int lighting = glIsEnabled( GL_LIGHTING );
		glDisable( GL_DEPTH_TEST );
		glDisable( GL_LIGHTING );
		glColor4f( 0 , 0 , 0 , 1 );
		glRasterPos2i( x , y );
		int len = int( strlen( str ) );
		for( int i=0 ; i<len ; i++ ) glutBitmapCharacter( font , str[i] );
		if( depth ) glEnable( GL_DEPTH_TEST );
		if( lighting ) glEnable( GL_LIGHTING );

		glMatrixMode( GL_PROJECTION );
		glPopMatrix();

		glMatrixMode( GL_MODELVIEW );
		glPopMatrix();

		glMatrixMode( matrixMode );
	}
	int Viewable::StringWidth( void* font , const char* format , ... )
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		return glutBitmapLength( font , (unsigned char*) str );
	}
	int Viewable::stringWidth( const char* format , ... ) const
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		return glutBitmapLength( font , (unsigned char*) str );
	}
	void Viewable::writeLeftString( int x , int y , const char* format , ... ) const
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		WriteLeftString( x , y , font , str );
	}
	void Viewable::writeRightString( int x , int y , const char* format , ... ) const
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		WriteLeftString( screenWidth-x-glutBitmapLength( font , (unsigned char*) str ) , y , font  ,str );
	}
	void Viewable::saveFrameBuffer( const char* fileName , int whichBuffer )
	{
		Pointer( float ) pixels = AllocPointer< float >( sizeof(float) * 3 * screenWidth * screenHeight );
		Pointer( unsigned char ) _pixels = AllocPointer< unsigned char >( sizeof(unsigned char) * 3 * screenWidth * screenHeight );
		glReadBuffer( whichBuffer );
		glReadPixels( 0 , 0 , screenWidth , screenHeight , GL_RGB , GL_FLOAT , pixels );
		for( int j=0 ; j<screenHeight ; j++ ) for( int i=0 ; i<screenWidth ; i++ ) for( int c=0 ; c<3 ; c++ )
		{
			int ii = int( pixels[ c + i * 3 + ( screenHeight - 1 - j ) * screenWidth * 3 ]*256 );
			if( ii<  0 ) ii =   0;
			if( ii>255 ) ii = 255;
			_pixels[ c + i * 3 + j * screenWidth * 3 ] = (unsigned char)ii;
		}
		FreePointer( pixels );
		ImageWriter::Write( fileName , _pixels , screenWidth , screenHeight , 3 , 95 );
		FreePointer( _pixels );
	}

	Viewable::Viewable( void )
	{
		callBacks.push_back( KeyboardCallBack( this , KEY_ESC    , "" , QuitCallBack ) );
		callBacks.push_back( KeyboardCallBack( this , KEY_CTRL_C , "" , QuitCallBack ) );
		callBacks.push_back( KeyboardCallBack( this , 'F' , "toggle fps"  , ToggleFPSCallBack ) );
		callBacks.push_back( KeyboardCallBack( this , 'H' , "toggle help" , ToggleHelpCallBack ) );
		callBacks.push_back( KeyboardCallBack( this , 'I' , "toggle info" , ToggleInfoCallBack ) );
		callBacks.push_back( KeyboardCallBack( this , 'i' , "save frame buffer" , "Ouput image" , SetFrameBufferCallBack ) );
		callBacks.push_back( KeyboardCallBack( this , 'f' , "toggle full screen" , ToggleFullScreenCallBack ) );
		snapshotName = videoHeader = NULL;
		flushImage = false;
		showHelp = showInfo = showFPS = true;
		_exitAfterSnapshot = _exitAfterVideo = false;
		_fullScreen = false;
		_screenWidth = _screenHeight = screenWidth = screenHeight = 512;
		font = GLUT_BITMAP_HELVETICA_12;
		fontHeight = 12;
		promptFont = GLUT_BITMAP_TIMES_ROMAN_24;
		promptFontHeight = 24;
		promptCallBack = NULL;
		strcpy( promptString , "" );
		promptLength = 0;

		_lastFPSTime = Time();
		_lastFPSCount = 0;
		_fps = 0;
	}
	void Viewable::setSnapshot( const char* sName , bool exitAfterSnapshot )
	{
		_exitAfterSnapshot = exitAfterSnapshot;
		snapshotName = new char[ strlen( sName ) + 1 ];
		strcpy( snapshotName , sName );
		showHelp = showInfo = showFPS = false;
		flushImage = true;
	}
	void Viewable::setVideo( const char* vHeader , int frames , bool exitAfterVideo )
	{
		_exitAfterVideo = exitAfterVideo;
		videoHeader = new char[ strlen( vHeader ) + 1 ];
		strcpy( videoHeader , vHeader );
		showHelp = showInfo = showFPS = false;
		_currentFrame = 0;
		_totalFrames = frames;
	}
	void Viewable::SetFrameBufferCallBack( Viewable* v , const char* prompt )
	{
		if( prompt )
		{
			v->snapshotName = new char[ strlen(prompt)+1 ];
			strcpy( v->snapshotName , prompt );
			v->flushImage = true;
		}
	}

	/////////////////////////
	// VisualizationViewer //
	/////////////////////////

	Viewable* Viewer::viewable = NULL;
	void Viewer::Run( Viewable* v , int argc , char* argv[] , const char* windowName )
	{
		viewable = v;
		glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
		glutInitWindowSize( viewable->screenWidth , viewable->screenHeight );
		glutInit( &argc , argv );
		glutCreateWindow( windowName );

		if( glewInit()!=GLEW_OK ) fprintf( stderr , "glewInit failed. Exiting...\n" ) , exit(0);

		glutIdleFunc      ( Idle );
		glutDisplayFunc   ( Display );
		glutReshapeFunc   ( Reshape );
		glutMouseFunc     ( MouseFunc );
		glutMotionFunc    ( MotionFunc );
		glutKeyboardFunc  ( KeyboardFunc );
		glutSpecialFunc   ( SpecialFunc );

		glutMainLoop();
	}
	void Viewer::Idle( void ){ viewable->Idle(); }
	void Viewer::KeyboardFunc( unsigned char key , int x , int y ){ viewable->KeyboardFunc( key , x , y ); }
	void Viewer::SpecialFunc( int key , int x , int y ){ viewable->SpecialFunc( key , x ,  y ); }
	void Viewer::Display( void ){ viewable->Display(); }
	void Viewer::Reshape( int w , int h ){ viewable->Reshape( w , h ); }
	void Viewer::MouseFunc( int button , int state , int x , int y ){ viewable->MouseFunc( button , state , x , y ); }
	void Viewer::MotionFunc( int x , int y ){ viewable->MotionFunc( x , y ); }
}
#endif // VISUALIZATION_INCLUDED