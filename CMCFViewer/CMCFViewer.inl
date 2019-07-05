#include <string>
#include "GL/glew.h"
#include "GL/glut.h"
#include <Eigen/Dense>
#include "Misha/SphericalGeometry.h"
#include "Misha/FEM.h"
#include "Misha/Timer.h"
#include "Misha/Visualization.h"
#include "Misha/Camera.h"
#include "Misha/Ply.h"
#include "Misha/SparseMatrix.h"
#include "Misha/SphericalSignals.h"
#include "Misha/CMCF.h"

template< typename Real >
struct CMCFVisualization : public Visualization::Viewable
{
protected:
	static void _SetNormals( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , std::vector< Point3D< Real > > &normals );
	typename SphericalGeometry::CMCF< Real >::Stats _currentStats;
	SphericalGeometry::Mesh< Real > _sMesh;
	std::vector< Point3D< Real > > _meshVertices , _normals , _colors;
	std::vector< TriangleIndex > _triangles;
	typename SphericalGeometry::CMCF< Real > *_cmcf;

	GLuint _vbo , _ebo;
	int _count;
	std::vector< Real > _vertexWeights , _vertexFaceWeights , _faceWeights;
	std::vector< std::vector< int > > _faces;

	void _setElementBufferObject( void );
	void _setVertexBufferObject( void );
	void _updateVertexBufferObject( void );
public:
	typename SphericalGeometry::CMCF< Real >::Parameters parameters;

	Camera camera;
	float zoom;
	int oldX , oldY , newX , newY;
	bool scaling , rotating , panning , useLight , edgeMode , showVectors , showColor;
	int animationSteps;

	GLfloat lightAmbient[4] , lightDiffuse[4] , lightSpecular[4] , shapeSpecular[4] , shapeSpecularShininess;
	char resolutionStr[1024];
	char statsStr[1024];
	char timeStr[1024];

	CMCFVisualization( void );
	~CMCFVisualization( void );

	void init( const char *sourceFile , Real stepSize );
	void advance( void );
	void outputSphericalMesh( const char *fileName ) const;

	void idle( void );
	void keyboardFunc( unsigned char key , int x , int y );
	void specialFunc( int key, int x, int y );
	void setProjectionMatrix( void );
	void setModelviewMatrix( void );
	void display( void );
	bool select( int x , int  y , Point3D< float >& out );
	void mouseFunc( int button , int state , int x , int y );
	void motionFunc( int x , int y );
	void passiveMotionFunc( int x , int y );

	static void TakeVideoCallBack( Visualization::Viewable* v , const char* str )
	{
		char _str[1024];
		strcpy( _str , str );
		char *header = NULL;
		int frames = 0;
		char *t = strtok( _str , "|" );
		int count = 0;
		while( t!=NULL )
		{
			if( count==0 ) header = t;
			else if( count==1 ) frames = atoi(t);
			t = strtok( NULL , "|" );
			count++;
		}
		if( count!=2 ) fprintf( stderr , "[WARNING] Failed to tokenize: %s (%d)\n" , str , count );
		else
		{
			( (CMCFVisualization*)v)->animationSteps = frames;
			v->setVideo( header , frames , false );
		}
	}
	static void ToggleColorsCallBack( Visualization::Viewable* v , const char* ){ ( (CMCFVisualization*)v)->showColor = !( (CMCFVisualization*)v)->showColor; }
	static void ToggleEdgesCallBack( Visualization::Viewable* v , const char* ){ ( (CMCFVisualization*)v)->edgeMode = !( (CMCFVisualization*)v)->edgeMode; }
	static void ToggleLightCallBack( Visualization::Viewable* v , const char* ){ ( (CMCFVisualization*)v)->useLight = !( (CMCFVisualization*)v)->useLight; }
	static void ToggleAnimationCallBack( Visualization::Viewable* v , const char* ){ ( (CMCFVisualization*)v)->animationSteps = ( (CMCFVisualization*)v)->animationSteps ? 0 : -1; }
	static void OutputSphericalMeshCallBack( Visualization::Viewable* v , const char* str ){ ( (CMCFVisualization*)v)->outputSphericalMesh( str ); }
	static void AdvanceCallBack( Visualization::Viewable* v , const char* ){ ( (CMCFVisualization*)v)->advance(); }
};

template< typename Real >
void CMCFVisualization< Real >::_SetNormals( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , std::vector< Point3D< Real > > &normals )
{
	normals.resize( vertices.size() );
#pragma omp parallel for
	for( int i=0 ; i<normals.size() ; i++ ) normals[i] = Point3D< Real >( normals[i] );

	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< double > v[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
		Point3D< double > n = Point3D< double >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
		for( int j=0 ; j<3 ; j++ ) normals[ triangles[i][j] ] += n;
	}

#pragma omp parallel for
	for( int i=0 ; i<normals.size() ; i++ ) normals[i] /= (Real)Length( normals[i] );
}

template< typename Real >
void CMCFVisualization< Real >::advance( void )
{
	Timer timer;

	_count++;

	double time;
	{
		timer.reset();
		_cmcf->advance();
		time = timer.elapsed();
	}
	_currentStats = _cmcf->getStats();

	_updateVertexBufferObject();

	sprintf( timeStr , "Time[%d]: %.2f(s)" , _count , time );
	sprintf( statsStr , "Delta / log(QC) / Radial: %.4f / %.4f / %.4f" , _currentStats.deformationScale , log(_currentStats.quasiConformalRatio) , _currentStats.radialDeviation );
}

template< typename Real >
CMCFVisualization< Real >::CMCFVisualization( void )
{
	_cmcf = NULL;
	_count = 0;
	edgeMode = false;
	screenWidth = screenHeight = 512;
	showColor = true;
	animationSteps = 0;
	useLight = true;
}

template< typename Real >
CMCFVisualization< Real >::~CMCFVisualization( void )
{
	if( _cmcf ) delete _cmcf;
	_cmcf = NULL;
}

template< typename Real >
void CMCFVisualization< Real >::init( const char *meshFile , Real stepSize )
{
	parameters.stepSize = stepSize;

	bool hasColor;
	{
		std::vector< PlyColorVertex< float > > _vertices;
		int fileType;
		bool readFlags[ PlyColorVertex< float >::ReadComponents ];
		PlyReadTriangles( In.value , _vertices , _triangles , PlyColorVertex< float >::ReadProperties , readFlags , PlyColorVertex< float >::ReadComponents , fileType );
		hasColor = ( readFlags[3] || readFlags[6] ) && (readFlags[4] && readFlags[7] ) && ( readFlags[5] && readFlags[8] );

		_meshVertices.resize( _vertices.size() );
#pragma omp parallel for
		for( int i=0 ; i<_vertices.size() ; i++ ) _meshVertices[i] = Point3D< double >( _vertices[i].point );
		_colors.resize( _vertices.size() );
		if( hasColor )
#pragma omp parallel for
			for( int i=0 ; i<_vertices.size() ; i++ ) _colors[i] = Point3D< double >( _vertices[i].color );
		else
		{
			_SetNormals( _meshVertices , _triangles , _colors );
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

			for( int i=0 ; i<_colors.size() ; i++ ) _colors[i] = NormalColor( _colors[i] );
		}
	}

	if( _cmcf ) delete _cmcf;
	_cmcf = new SphericalGeometry::CMCF< Real >( _meshVertices , _triangles , _sMesh , parameters );
	_currentStats = _cmcf->getStats();

	///////////////////////////////////////
	camera.position = Point3D< double >( 0 , 0 , 1. ) , zoom = 1.f / 0.95f;
	_vbo = _ebo = 0;

	lightAmbient [0] = lightAmbient [1] = lightAmbient [2] = 0.25f , lightAmbient [3] = 1.f;
	lightDiffuse [0] = lightDiffuse [1] = lightDiffuse [2] = 0.70f , lightDiffuse [3] = 1.f;
	lightSpecular[0] = lightSpecular[1] = lightSpecular[2] = 0.25f , lightSpecular[3] = 1.f;
	shapeSpecular[0] = shapeSpecular[1] = shapeSpecular[2] = 1.00f , shapeSpecular[3] = 1.f;
	shapeSpecularShininess = 128;

	rotating = scaling = panning = false;
	callBacks.push_back( Visualization::KeyboardCallBack( this , 'V' , "video" , "Header|Frames" , TakeVideoCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , 'l' , "toggle light" , ToggleLightCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , 'c' , "toggle colors" , ToggleColorsCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , 'o' , "output mesh" , "Name" , OutputSphericalMeshCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , 'e' , "toggle edges" , ToggleEdgesCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , '+' , "advance" , AdvanceCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , ' ' , "advance" , ToggleAnimationCallBack ) );
	info.push_back( statsStr );
	info.push_back( timeStr );
	info.push_back( resolutionStr );

	sprintf( resolutionStr , "Vertices / Triangles: %d / %d" , (int)_sMesh.vertices.size() , (int)_triangles.size() );
	sprintf( timeStr , "Time[%d]: ---" , _count );
	sprintf( statsStr , "Delta / log(QC) / Radial: --- / %.4f / %.4f" , log(_currentStats.quasiConformalRatio) , _currentStats.radialDeviation );
}

template< typename Real >
void CMCFVisualization< Real >::outputSphericalMesh( const char *fileName ) const { _sMesh.write( fileName , _meshVertices , _colors ); }

template< typename Real >
void CMCFVisualization< Real >::_setElementBufferObject( void )
{
	glGenBuffers( 1 , &_ebo );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , _ebo );
	glBufferData( GL_ELEMENT_ARRAY_BUFFER , _triangles.size() * sizeof( int ) * 3 ,  &_triangles[0] , GL_STATIC_DRAW );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , 0 );
}

template< typename Real >
void CMCFVisualization< Real >::_setVertexBufferObject( void )
{
	glGenBuffers( 1 , &_vbo );

	Real r=0;
	for( int i=0 ; i<_sMesh.vertices.size() ; i++ ) r = std::max< Real >( r , Point3D< Real >::SquareNorm( _sMesh.vertices[i] ) );
	r = (Real)sqrt(r);

	_SetNormals( _sMesh.vertices , _triangles , _normals );
	Point3D< float > *vertexData = new Point3D< float >[ _sMesh.vertices.size()*3 ];

	glBindBuffer( GL_ARRAY_BUFFER , _vbo );
	for( int i=0 ; i<_sMesh.vertices.size() ; i++ ) vertexData[ _sMesh.vertices.size()*0 + i ] = Point3D< float >( _sMesh.vertices[i] / r );
	for( int i=0 ; i<_sMesh.vertices.size() ; i++ ) vertexData[ _sMesh.vertices.size()*1 + i ] = Point3D< float >( _normals[i] );
	for( int i=0 ; i<_sMesh.vertices.size() ; i++ ) vertexData[ _sMesh.vertices.size()*2 + i ] = Point3D< float >( _colors[i]/255.f );
	glBufferData( GL_ARRAY_BUFFER , 9 * _sMesh.vertices.size() * sizeof( float ) , (float*)vertexData , GL_DYNAMIC_DRAW );
	glBindBuffer( GL_ARRAY_BUFFER , 0 );

	delete[] vertexData;
}

template< typename Real >
void CMCFVisualization< Real >::_updateVertexBufferObject( void )
{
	Real r=0;
	for( int i=0 ; i<_sMesh.vertices.size() ; i++ ) r = std::max< Real >( r , Point3D< Real >::SquareNorm( _sMesh.vertices[i] ) );
	r = (Real)sqrt(r);

	_SetNormals( _sMesh.vertices , _triangles , _normals );
	Point3D< float > *vertexData = new Point3D< float >[ _sMesh.vertices.size()*2 ];

	glBindBuffer( GL_ARRAY_BUFFER , _vbo );
	for( int i=0 ; i<_sMesh.vertices.size() ; i++ ) vertexData[ _sMesh.vertices.size()*0 + i ] = Point3D< float >( _sMesh.vertices[i] / r );
	for( int i=0 ; i<_sMesh.vertices.size() ; i++ ) vertexData[ _sMesh.vertices.size()*1 + i ] = Point3D< float >( _normals[i] );
	glBufferSubData( GL_ARRAY_BUFFER , 0 , 6 * _sMesh.vertices.size() * sizeof( float ) , (float*)vertexData );

	delete[] vertexData;
}

template< typename Real >
void CMCFVisualization< Real >::setProjectionMatrix( void )
{
	int width = screenWidth , height = screenHeight;

	glViewport( 0 , 0 , width , height );

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();

	float ar = (float)width/(float)height , ar_r = 1.f/ar;
	if( width>height ) glOrtho( -ar*zoom , ar*zoom , -zoom , zoom , -2.f , 2.f );
	else               glOrtho( -zoom , zoom , -ar_r*zoom , ar_r*zoom , -2.f , 2.f );
}

template< typename Real >
void CMCFVisualization< Real >::setModelviewMatrix( void )
{
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	camera.draw();
	GLfloat scale[16];
	for( int i=0 ; i<16 ; i++ ) scale[i] = 0;
	scale[4*0+0] = scale[4*1+1] = scale[4*2+2] = 1.f;
	scale[4*3+3] = 1.f;
	glMultMatrixf( scale );
}

template< typename Real >
void CMCFVisualization< Real >::display( void )
{
	if( edgeMode )
	{
		glPolygonMode( GL_FRONT_AND_BACK , GL_LINE );
		// 	polygonOffsetFactor = polygonOffsetUnits = 1.f;
		int smoothState = glIsEnabled( GL_LINE_SMOOTH );
		int blendState = glIsEnabled( GL_BLEND );
		int lightState = glIsEnabled( GL_LIGHTING );
		int multisampleState = glIsEnabled( GL_MULTISAMPLE );
		if( true )
		{
			glDisable( GL_MULTISAMPLE );
			glEnable( GL_BLEND );
			glEnable( GL_LINE_SMOOTH );
			glHint( GL_LINE_SMOOTH_HINT , GL_NICEST );
			glDepthMask( GL_FALSE );
			glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
		}
		glDisable( GL_LIGHTING );
		glColor3f( 0.25f , 0.25f , 0.25f );
		glLineWidth( 0.125f );
	}
	else
	{
		glPolygonMode( GL_FRONT_AND_BACK , GL_FILL );
		glDepthMask( GL_TRUE );
	}

	if( !_ebo ) _setElementBufferObject();
	if( !_vbo ) _setVertexBufferObject();
	setProjectionMatrix();
	setModelviewMatrix();

	glEnable( GL_NORMALIZE );
	glEnable( GL_DEPTH_TEST );
	glDisable( GL_CULL_FACE );

	GLfloat lPosition[4];
	{
		Point3D< float > d = camera.up + camera.right - camera.forward*5;
		lPosition[0] = d[0] , lPosition[1] = d[1] , lPosition[2] = d[2];
	}
	lPosition[3] = 0.0;	
	glLightModeli( GL_LIGHT_MODEL_LOCAL_VIEWER , GL_FALSE );
	glLightModeli( GL_LIGHT_MODEL_TWO_SIDE , GL_TRUE );
	glLightfv( GL_LIGHT0 , GL_AMBIENT , lightAmbient );
	glLightfv( GL_LIGHT0 , GL_DIFFUSE , lightDiffuse );
	glLightfv( GL_LIGHT0 , GL_SPECULAR , lightSpecular );
	glLightfv( GL_LIGHT0 , GL_POSITION , lPosition );
	glEnable( GL_LIGHT0 );
	if( useLight ) glEnable ( GL_LIGHTING );
	else           glDisable( GL_LIGHTING );

	glColorMaterial( GL_FRONT_AND_BACK , GL_AMBIENT_AND_DIFFUSE );
	glEnable( GL_COLOR_MATERIAL );

	glMaterialfv( GL_FRONT_AND_BACK , GL_SPECULAR  , shapeSpecular );
	glMaterialf ( GL_FRONT_AND_BACK , GL_SHININESS , shapeSpecularShininess );

	glBindBuffer( GL_ARRAY_BUFFER , _vbo );
	glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_NORMAL_ARRAY );
	if( showColor ) glEnableClientState( GL_COLOR_ARRAY );
	else            glDisableClientState( GL_COLOR_ARRAY );
	const std::vector< Point3D< Real > > &sphereVertices = _sMesh.vertices;
	glVertexPointer  ( 3 , GL_FLOAT , 0 , (GLubyte*)NULL + sizeof( float ) * _sMesh.vertices.size() * 0 );
	glNormalPointer  (     GL_FLOAT , 0 , (GLubyte*)NULL + sizeof( float ) * _sMesh.vertices.size() * 3 );
	glColorPointer   ( 3 , GL_FLOAT , 0 , (GLubyte*)NULL + sizeof( float ) * _sMesh.vertices.size() * 6 );
	glColor3f( 0.75f , 0.75f , 0.75f );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , _ebo );
	glDrawElements( GL_TRIANGLES , (GLsizei)( _triangles.size() * 3 ) , GL_UNSIGNED_INT , NULL );
	glDisableClientState( GL_NORMAL_ARRAY );
	glDisableClientState( GL_COLOR_ARRAY );
	glBindBuffer( GL_ARRAY_BUFFER , 0 );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , 0 );
}

template< typename Real >
void CMCFVisualization< Real >::mouseFunc( int button , int state , int x , int y )
{
	newX = x , newY = y;

	rotating = scaling = panning = false;
	if( button==GLUT_LEFT_BUTTON  )
		if( glutGetModifiers() & GLUT_ACTIVE_CTRL ) panning = true;
		else                                        rotating = true;
	else if( button==GLUT_RIGHT_BUTTON )            scaling = true;
}
template< typename Real >
void CMCFVisualization< Real >::motionFunc( int x , int y )
{
	oldX = newX , oldY = newY , newX = x , newY = y;

	int imageSize = std::min< int >( screenWidth , screenHeight );
	float rel_x = (newX - oldX) / (float)imageSize * 2;
	float rel_y = (newY - oldY) / (float)imageSize * 2;
	float pRight = -rel_x * zoom , pUp = rel_y * zoom;
	float sForward = rel_y*4;
	float rRight = rel_y , rUp = rel_x;

	float scale = zoom;

	if     ( scaling ) zoom *= (float)pow( 0.9 , (double)sForward );
	else if( panning ) camera.translate( camera.right * pRight + camera.up * pUp );
	else if( rotating ) camera.rotateUp( rUp*scale ) , camera.rotateRight( rRight*scale );

	glutPostRedisplay();
}
template< typename Real >
void CMCFVisualization< Real >::passiveMotionFunc( int x , int y )
{
}
template< typename Real >
void CMCFVisualization< Real >::idle( void )
{
	if( !promptCallBack )
	{
		if( animationSteps )
		{
			advance();
			if( animationSteps>0 ) animationSteps--;
			glutPostRedisplay();
		}
	}
}

template< typename Real >
void CMCFVisualization< Real >::keyboardFunc( unsigned char key , int x , int y )
{
	float scale = zoom;
	switch( key )
	{
	case 'q': camera.rotateUp( -M_PI/128 * scale ) ; break;
	case 'Q': camera.rotateUp( -M_PI/  4 ) ; break;
	case 'w': camera.rotateUp(  M_PI/128 * scale ) ; break;
	case 'W': camera.rotateUp(  M_PI/  4 ) ; break;
	case 'a': camera.rotateRight( -M_PI/128 * scale ) ; break;
	case 'A': camera.rotateRight( -M_PI/  4 ) ; break;
	case 'z': camera.rotateRight(  M_PI/128 * scale ) ; break;
	case 'Z': camera.rotateRight(  M_PI/  4 ) ; break;
	case 's': camera.rotateForward( -M_PI/128 * scale ) ; break;
	case 'S': camera.rotateForward( -M_PI/  4 ) ; break;
	case 'x': camera.rotateForward(  M_PI/128 * scale ) ; break;
	case 'X': camera.rotateForward(  M_PI/  4 ) ; break;
	}
}

template< typename Real >
void CMCFVisualization< Real >::specialFunc( int key , int x , int y )
{
	float stepSize = 10.f / ( screenWidth + screenHeight );
	if( glutGetModifiers()&GLUT_ACTIVE_CTRL ) stepSize /= 16;
	float panSize = stepSize*2 , scaleSize = stepSize*2;

	switch( key )
	{
	case Visualization::KEY_UPARROW:    zoom *= (float)pow( 0.9 , 1 ) ; break;
	case Visualization::KEY_DOWNARROW:  zoom *= (float)pow( 0.9 , -1 ) ; break;
	case Visualization::KEY_LEFTARROW:  camera.translate(  camera.right * panSize ) ; break;
	case Visualization::KEY_RIGHTARROW: camera.translate( -camera.right * panSize ) ; break;
	case Visualization::KEY_PGUP:       camera.translate( -camera.up * panSize ) ; break;
	case Visualization::KEY_PGDN:       camera.translate(  camera.up * panSize ) ; break;
	}
	glutPostRedisplay();
}
