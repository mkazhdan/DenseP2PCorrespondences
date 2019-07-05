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
#include "Misha/AuthalicEvolution.h"

template< typename Real > Point3D< Real > NormalColor( Point3D< Real > n );

template< typename Real >
struct AuthalicEvolutionVisualization : public Visualization::Viewable
{
protected:
	static const unsigned int _Windows = 2;
	unsigned int _windows;
	Real _originalFlippedTriangleFraction;
	typename SphericalGeometry::AuthalicEvolution< Real >::Stats _originalStats , _currentStats;
	std::vector< TriangleIndex > _triangles;

	SphericalGeometry::Mesh< Real > _sMesh;
	std::vector< Point3D< Real > > _meshVertices;
	typename SphericalGeometry::AuthalicEvolution< Real > *_authalicEvolution;

	std::vector< Real > _meshLogScaleFactors;
	std::vector< Point3D< Real > > _sphereVertices , _colors[2];
	GLuint _vbo[_Windows] , _ebo[_Windows];
	Real _vectorScale;
	int _count;
	std::vector< std::pair< Point3D< Real > , Point3D< Real > > > _flowFieldSamples;
	std::vector< Real > _vertexWeights , _vertexFaceWeights , _faceWeights;
	std::vector< std::vector< int > > _faces;

	void _setElementBufferObject( void );
	void _setVertexBufferObject( void );
	void _updateVertexBufferObject( void );
	void _setScaleFactorColors( void );

	unsigned int _window( int x ) const;
public:
	typename SphericalGeometry::AuthalicEvolution< Real >::Parameters parameters;

	Camera camera;
	float zoom;
	int oldX , oldY , newX , newY;
	bool scaling , rotating , panning , useLight , showGrid , edgeMode , showVectors , remapScaleFactors , showColor;
	int animationSteps;
	bool verbose;

	GLfloat lightAmbient[4] , lightDiffuse[4] , lightSpecular[4] , shapeSpecular[4] , shapeSpecularShininess;
	char resolutionStr[1024];
	char statsStr[1024];
	char timeStr[1024];

	AuthalicEvolutionVisualization( void );
	~AuthalicEvolutionVisualization( void );

	void init( const char *sourceFile , Real smooth , Real stepSize , int resolution , int subSteps , bool showMeshOnly );
	void advance( void );
	void setFibonacci( unsigned int samples );
	void setFibonacciSamples( void );
	void outputSphericalMesh( const char *fileName ) const;

	void idle( void );
	void keyboardFunc( unsigned char key , int x , int y );
	void specialFunc( int key, int x, int y );
	void setProjectionMatrix( unsigned int modelIndex );
	void setModelviewMatrix( unsigned int modelIndex );
	void display( void );
	bool select( int x , int  y , Point3D< float >& out );
	void mouseFunc( int button , int state , int x , int y );
	void motionFunc( int x , int y );
	void passiveMotionFunc( int x , int y );

	static void VectorScaleUpCallBack( Visualization::Viewable* v , const char* ){ ( (AuthalicEvolutionVisualization*)v)->_vectorScale *= 1.1f; }
	static void VectorScaleDownCallBack( Visualization::Viewable* v , const char* ){ ( (AuthalicEvolutionVisualization*)v)->_vectorScale /= 1.1f; }
	static void SampleScaleUpCallBack( Visualization::Viewable* v , const char* ){ ( (AuthalicEvolutionVisualization*)v)->setFibonacci( (int)ceil( ( (AuthalicEvolutionVisualization*)v)->_flowFieldSamples.size()*1.1 ) ); }
	static void SampleScaleDownCallBack( Visualization::Viewable* v , const char* ){ ( (AuthalicEvolutionVisualization*)v)->setFibonacci( (int)ceil( ( (AuthalicEvolutionVisualization*)v)->_flowFieldSamples.size()/1.1 ) ); }
	static void ToggleRemapScaleFactorsCallBack( Visualization::Viewable* v , const char* )
	{
		( (AuthalicEvolutionVisualization*)v)->remapScaleFactors = !( (AuthalicEvolutionVisualization*)v)->remapScaleFactors;
		( (AuthalicEvolutionVisualization*)v)->_updateVertexBufferObject();
		glutPostRedisplay();
	}
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
			( (AuthalicEvolutionVisualization*)v)->animationSteps = frames;
			v->setVideo( header , frames , false );
		}
	}
	static void ToggleVectorsCallBack( Visualization::Viewable* v , const char* ){ ( (AuthalicEvolutionVisualization*)v)->showVectors = !( (AuthalicEvolutionVisualization*)v)->showVectors; }
	static void ToggleColorsCallBack( Visualization::Viewable* v , const char* ){ ( (AuthalicEvolutionVisualization*)v)->showColor = !( (AuthalicEvolutionVisualization*)v)->showColor; }
	static void ToggleGridCallBack( Visualization::Viewable* v , const char* ){ ( (AuthalicEvolutionVisualization*)v)->showGrid = !( (AuthalicEvolutionVisualization*)v)->showGrid; }
	static void ToggleEdgesCallBack( Visualization::Viewable* v , const char* ){ ( (AuthalicEvolutionVisualization*)v)->edgeMode = !( (AuthalicEvolutionVisualization*)v)->edgeMode; }
	static void ToggleLightCallBack( Visualization::Viewable* v , const char* ){ ( (AuthalicEvolutionVisualization*)v)->useLight = !( (AuthalicEvolutionVisualization*)v)->useLight; }
	static void ToggleAnimationCallBack( Visualization::Viewable* v , const char* ){ ( (AuthalicEvolutionVisualization*)v)->animationSteps = ( (AuthalicEvolutionVisualization*)v)->animationSteps ? 0 : -1; }
	static void OutputSphericalMeshCallBack( Visualization::Viewable* v , const char* str ){ ( (AuthalicEvolutionVisualization*)v)->outputSphericalMesh( str ); }
	static void AdvanceCallBack( Visualization::Viewable* v , const char* ){ ( (AuthalicEvolutionVisualization*)v)->advance(); }
};

template< typename Real >
unsigned int AuthalicEvolutionVisualization< Real >::_window( int x ) const
{
	return x / ( screenWidth/_windows );
}

template< typename Real >
void AuthalicEvolutionVisualization< Real >::setFibonacci( unsigned int samples )
{
	static const Real PHI = (Real)( ( 1.+sqrt(5.) ) / 2. );
	static const Real PHI_INV = (Real)( 1. / PHI );
	_flowFieldSamples.resize( samples );
#pragma omp parallel for
	for( int i=0 ; i<(int)samples ; i++ )
	{
		Real theta = (Real)acos( 1 - ( (Real)(2*i+1) ) / samples );
		Real phi = PHI_INV * 2 * i * M_PI;
		_flowFieldSamples[i].first = Point3D< Real >( cos(phi)*sin(theta) , sin(phi)*sin(theta) , cos(theta) );
	}
	setFibonacciSamples();
}
template< typename Real >
void AuthalicEvolutionVisualization< Real >::setFibonacciSamples( void )
{
#pragma omp parallel for
	for( int i=0 ; i<_flowFieldSamples.size() ; i++ ) _flowFieldSamples[i].second = _authalicEvolution->flowField( _flowFieldSamples[i].first );
}

template< typename Real >
void AuthalicEvolutionVisualization< Real >::advance( void )
{
	Timer timer;

	_count++;

	double time;
	{
		timer.reset();
		_authalicEvolution->advance();
		time = timer.elapsed();
	}
	_currentStats = _authalicEvolution->getStats( _meshLogScaleFactors );
	setFibonacciSamples();

	_updateVertexBufferObject();

	Real flippedTriangles = _authalicEvolution->flippedTriangleFraction();
	sprintf( timeStr ,  "Time[%d]: %.2f(s)" , _count , time );
	sprintf( statsStr , "[min max] +/- dev (flipped): [%.3f %.3f] +/- %.3f (%.4f)" , _currentStats.min , _currentStats.max , _currentStats.dev , flippedTriangles );
//	sprintf( statsStr , "[%.3f %.3f] +/- %.3f (%.4f) -> [%.3f %.3f] +/- %.3f (%.4f)" , _originalStats.min ,  _originalStats.max , _originalStats.dev , _originalFlippedTriangleFraction , _currentStats.min , _currentStats.max , _currentStats.dev , flippedTriangles );
	if( verbose )
	{
		if( _count==1 ) printf( "%d] [%.3f %.3f] +/- %.3f (%.4f)\n" , 0 , _originalStats.min , _originalStats.max , _originalStats.dev , _originalFlippedTriangleFraction );
		printf( "%d] [%.3f %.3f] +/- %.3f (%.4f)\n" , _count , _currentStats.min , _currentStats.max , _currentStats.dev , flippedTriangles );
	}
}

template< typename Real >
AuthalicEvolutionVisualization< Real >::AuthalicEvolutionVisualization( void )
{
	_authalicEvolution = NULL;
	_vectorScale = 1;
	_count = 0;
	edgeMode = false;
	remapScaleFactors = false;
	screenHeight = 512;
	showColor = true;
	animationSteps = 0;
	showGrid = true;
	useLight = true;
	showVectors = true;
}

template< typename Real >
AuthalicEvolutionVisualization< Real >::~AuthalicEvolutionVisualization( void )
{
	if( _authalicEvolution ) delete _authalicEvolution;
	_authalicEvolution = NULL;
}

template< typename Real >
void AuthalicEvolutionVisualization< Real >::init( const char *sphereMeshFile , Real smooth , Real stepSize , int resolution , int subSteps , bool showMeshOnly )
{
	parameters.resolution = resolution;
	parameters.smoothness = smooth;
	parameters.stepSize = stepSize;
	parameters.subSteps = subSteps;
	_windows = showMeshOnly ? 1 : _Windows;
	screenWidth = screenHeight * _windows;

	_sMesh.read( sphereMeshFile , _meshVertices , _colors[0] );
	if( _authalicEvolution ) delete _authalicEvolution;
	_authalicEvolution = new SphericalGeometry::AuthalicEvolution< Real >( _sMesh , _meshVertices , parameters );
	_originalStats = _authalicEvolution->getStats( _meshLogScaleFactors );
	_originalFlippedTriangleFraction = _authalicEvolution->flippedTriangleFraction();

	if( !_colors[0].size() )
	{
		_colors[0].resize( _meshVertices.size() );
		for( int t=0 ; t<_sMesh.polygons.size() ; t++ )
		{
			Point3D< Real > v[] = { _meshVertices[ _sMesh.polygons[t][0] ] , _meshVertices[ _sMesh.polygons[t][1] ] , _meshVertices[ _sMesh.polygons[t][2] ] };
			Point3D< Real > n = Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
			for( int j=0 ; j<3 ; j++ ) _colors[0][ _sMesh.polygons[t][j] ] += n;
		}
#pragma omp parallel for
		for( int i=0 ; i<_colors[0].size() ; i++ ) _colors[0][i] = NormalColor( -_colors[0][i]/Point3D< Real >::Length( _colors[0][i] ) );
	}

	// Triangulate the mesh
	_triangles = Triangulate( _meshVertices , _sMesh.polygons );
	_sphereVertices = _sMesh.vertices;

	// Make the original mesh have area 4\pi
	{
		Real area = 0;
		for( int i=0 ; i<_triangles.size() ; i++ )
		{
			Point3D< Real > v[] = { _meshVertices[ _triangles[i][0] ] , _meshVertices[ _triangles[i][1] ] , _meshVertices[ _triangles[i][2] ] };
			area += Point3D< Real >::Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2;
		}
		Real scl = (Real)sqrt( 4. * M_PI / area );
		for( int i=0 ; i<_meshVertices.size() ; i++ ) _meshVertices[i] *= scl;
	}

	_faces.resize( SphericalGeometry::Tessellation< Real >::GridFaceNum( resolution ) );
	_sphereVertices.resize( SphericalGeometry::Tessellation< Real >::GridVertexNum( resolution ) );
#pragma omp parallel for
	for( int i=0 ; i<_sphereVertices.size() ; i++ ) _sphereVertices[i] = SphericalGeometry::Tessellation< Real >::GridVertex( i , resolution );
#pragma omp parallel for
	for( int i=0 ; i<_faces.size() ; i++ )
	{
		std::vector< unsigned int > face = SphericalGeometry::Tessellation< Real >::GridFace( i , resolution );
		_faces[i].resize( face.size() );
		for( int j=0 ; j<face.size() ; j++ ) _faces[i][j] = (int)face[j];
	}

	_vertexWeights.resize( _sphereVertices.size() );
	_vertexFaceWeights.resize( _sphereVertices.size() );
	_faceWeights.resize( _faces.size() );
	Real area = 0;
#pragma omp parallel for reduction( + : area )
	for( int i=0 ; i<_faces.size() ; i++ )
	{
		Real a = SphericalGeometry::Tessellation< Real >::GridFaceArea( i , resolution );
		_faceWeights[i] = a;
		area += a * _faces[i].size();
		for( int j=0 ; j<_faces[i].size() ; j++ )
#pragma omp atomic
			_vertexFaceWeights[ _faces[i][j] ] += a;
	}
#pragma omp parallel for
	for( int i=0 ; i<_sphereVertices.size() ; i++ ) _vertexWeights[i] = _vertexFaceWeights[i] / area;
	setFibonacci( 1024 );

	///////////////////////////////////////
	camera.position = Point3D< double >( 0 , 0 , 2. ) , zoom = 1.f / 0.95f;
	for( int i=0 ; i<_Windows ; i++ ) _vbo[i] = _ebo[i] = 0;

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
	callBacks.push_back( Visualization::KeyboardCallBack( this , 'v' , "toggle vectors" , ToggleVectorsCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , 'r' , "remap scale factors" , ToggleRemapScaleFactorsCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , '+' , "advance" , AdvanceCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , '{' , "sample down" , SampleScaleDownCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , '}' , "sample up" , SampleScaleUpCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , '[' , "scale down" , VectorScaleDownCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , ']' , "scale up" , VectorScaleUpCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , ' ' , "advance" , ToggleAnimationCallBack ) );
	info.push_back( statsStr );
	info.push_back( timeStr );
	info.push_back( resolutionStr );

	sprintf( resolutionStr , "Vertices / Triangles / Spherical Resolution: %d / %d / %d" , (int)_sMesh.vertices.size() , (int)_triangles.size() , (int)resolution );
	sprintf( timeStr , "Time[%d]: ---" , _count );
	sprintf( statsStr , "[min max] +/- dev (flipped): [%.3f %.3f] +/- %.3f (%.4f)" , _originalStats.min , _originalStats.max , _originalStats.dev , _originalFlippedTriangleFraction );
}

template< typename Real >
void AuthalicEvolutionVisualization< Real >::outputSphericalMesh( const char *fileName ) const { _sMesh.write( fileName , _meshVertices , _colors[0] ); }

template< typename Real >
void AuthalicEvolutionVisualization< Real >::_setElementBufferObject( void )
{
	glGenBuffers( _windows , _ebo );
	for( int i=0 ; i<2 ; i++ )
	{
		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , _ebo[i] );
		glBufferData( GL_ELEMENT_ARRAY_BUFFER , _triangles.size() * sizeof( int ) * 3 ,  &_triangles[0] , GL_STATIC_DRAW );
		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , 0 );
	}
}

template< typename Real >
void AuthalicEvolutionVisualization< Real >::_setScaleFactorColors( void )
{
	auto HeatMapFunction = [&]( Real v )
	{
		if( remapScaleFactors ) v /= std::max< Real >( (Real)fabs( _currentStats.min ) , (Real)fabs( _currentStats.max ) );
		else                    v /= std::max< Real >( (Real)fabs( _originalStats.min ) , (Real)fabs( _originalStats.max ) );
		return Point3D< float >( HeatMap( v ) );
	};
	_colors[1].resize( _sMesh.vertices.size() );

	for( int i=0 ; i<_colors[1].size() ; i++ ) _colors[1][i] = HeatMapFunction( _meshLogScaleFactors[i] ) * 255.f;
}

template< typename Real >
void AuthalicEvolutionVisualization< Real >::_setVertexBufferObject( void )
{
	_setScaleFactorColors();

	glGenBuffers( _windows , _vbo );

	Point3D< float > *vertexData = new Point3D< float >[ std::max< size_t >( _sMesh.vertices.size() , _sphereVertices.size() )*3 ];

	for( unsigned int w=0 ; w<_windows ; w++ )
	{
		const std::vector< Point3D< Real > > &sphereVertices = _sMesh.vertices;
		glBindBuffer( GL_ARRAY_BUFFER , _vbo[w] );
		for( int i=0 ; i<sphereVertices.size() ; i++ ) vertexData[i] = vertexData[ sphereVertices.size() + i ] = Point3D< float >( sphereVertices[i] );
		for( int i=0 ; i<sphereVertices.size() ; i++ ) vertexData[ sphereVertices.size()*2 + i ] = Point3D< float >( _colors[w][i]/255.f );
		glBufferData( GL_ARRAY_BUFFER , 9 * sphereVertices.size() * sizeof( float ) , (float*)vertexData , GL_DYNAMIC_DRAW );
		glBindBuffer( GL_ARRAY_BUFFER , 0 );
	}

	delete[] vertexData;
}

template< typename Real >
void AuthalicEvolutionVisualization< Real >::_updateVertexBufferObject( void )
{
	_setScaleFactorColors();

	Point3D< float > *vertexData = new Point3D< float >[ std::max< size_t >( _sMesh.vertices.size() , _sphereVertices.size() )*3 ];


	for( unsigned int w=0 ; w<_windows ; w++ )
	{
		const std::vector< Point3D< Real > > &sphereVertices = _sMesh.vertices;
		glBindBuffer( GL_ARRAY_BUFFER , _vbo[w] );
		for( int i=0 ; i<sphereVertices.size() ; i++ ) vertexData[i] = vertexData[ sphereVertices.size() + i ] = Point3D< float >( sphereVertices[i] );
		for( int i=0 ; i<sphereVertices.size() ; i++ ) vertexData[ sphereVertices.size()*2 + i ] = Point3D< float >( _colors[w][i]/255.f );
		glBufferSubData( GL_ARRAY_BUFFER , 0 , 9 * sphereVertices.size() * sizeof( float ) , (float*)vertexData );
	}

	delete[] vertexData;
}

template< typename Real >
void AuthalicEvolutionVisualization< Real >::setProjectionMatrix( unsigned int window )
{
	int width = screenWidth/_windows , height = screenHeight;

	glViewport( (window*screenWidth)/_windows , 0 , width , height );

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();

	float ar = (float)width/(float)height , ar_r = 1.f/ar;
	if( width>height ) glOrtho( -ar*zoom , ar*zoom , -zoom , zoom , -2.f , 2.f );
	else               glOrtho( -zoom , zoom , -ar_r*zoom , ar_r*zoom , -2.f , 2.f );
}

template< typename Real >
void AuthalicEvolutionVisualization< Real >::setModelviewMatrix( unsigned int window )
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
void AuthalicEvolutionVisualization< Real >::display( void )
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

	for( unsigned int window=0 ; window<_windows ; window++ )
	{
		if( !_ebo[window] ) _setElementBufferObject();
		if( !_vbo[window] ) _setVertexBufferObject();
		setProjectionMatrix( window );
		setModelviewMatrix( window );

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

		glBindBuffer( GL_ARRAY_BUFFER , _vbo[window] );
		glEnableClientState( GL_VERTEX_ARRAY );
		glEnableClientState( GL_NORMAL_ARRAY );
		if( showColor ) glEnableClientState( GL_COLOR_ARRAY );
		else            glDisableClientState( GL_COLOR_ARRAY );
		const std::vector< Point3D< Real > > &sphereVertices = _sMesh.vertices;
		glVertexPointer  ( 3 , GL_FLOAT , 0 , (GLubyte*)NULL + sizeof( float ) * sphereVertices.size() * 0 );
		glNormalPointer  (     GL_FLOAT , 0 , (GLubyte*)NULL + sizeof( float ) * sphereVertices.size() * 3 );
		glColorPointer   ( 3 , GL_FLOAT , 0 , (GLubyte*)NULL + sizeof( float ) * sphereVertices.size() * 6 );
		glColor3f( 0.75f , 0.75f , 0.75f );
		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , _ebo[window] );
		glDrawElements( GL_TRIANGLES , (GLsizei)( _triangles.size() * 3 ) , GL_UNSIGNED_INT , NULL );
		glDisableClientState( GL_NORMAL_ARRAY );
		glDisableClientState( GL_COLOR_ARRAY );
		glBindBuffer( GL_ARRAY_BUFFER , 0 );
		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , 0 );

		if( window==1 && showVectors )
		{
			Point3D< float > f = camera.forward / 256;
			glDisable( GL_LIGHTING );
			glBegin( GL_TRIANGLES );
			for( int i=0 ; i<(int)_flowFieldSamples.size() ; i++ )
			{
				Point3D< float > p = Point3D< float >( _flowFieldSamples[i].first );
				Point3D< float > d = Point3D< float >( _flowFieldSamples[i].second * parameters.stepSize * _vectorScale );
				Point3D< float > _n = Point3D< float >::CrossProduct( d , p );
				_n /= 20;
				glColor3f( 1.0f , 1.0f , 1.0f ) , glVertex3f( p[0]+d[0] , p[1]+d[1] , p[2]+d[2] );
				glColor3f( 0.0f , 0.0f , 0.0f ) , glVertex3f( p[0]-_n[0] , p[1]-_n[1] , p[2]-_n[2] ) , glVertex3f( p[0]+_n[0] , p[1]+_n[1] , p[2]+_n[2] );
			}
			glEnd();
		}
	}

	glViewport( 0 , 0 , screenWidth , screenHeight );

	if( showGrid )
	{
		glMatrixMode( GL_PROJECTION );
		glPushMatrix();
		glLoadIdentity();
		glOrtho( 0 , screenWidth , 0 , screenHeight , 0.f , 1.f );
		glMatrixMode( GL_MODELVIEW );
		glPushMatrix();
		glLoadIdentity();

		glDisable( GL_DEPTH_TEST );
		glDisable( GL_LIGHTING );
		glColor3f( 0.f , 0.f , 0.f );
		glLineWidth( 2 );
		glBegin( GL_LINES );
		for( unsigned int i=1 ; i<_windows ; i++ )
		{
			glVertex2i( (i*screenWidth)/_windows , 0 );
			glVertex2i( (i*screenWidth)/_windows , screenHeight );
		}
		glEnd();
		glEnable( GL_DEPTH_TEST );
		glMatrixMode( GL_PROJECTION ) ; glPopMatrix();
		glMatrixMode( GL_MODELVIEW  ) ; glPopMatrix();
	}
}

template< typename Real >
void AuthalicEvolutionVisualization< Real >::mouseFunc( int button , int state , int x , int y )
{
	unsigned int window = _window(x);

	newX = x , newY = y;

	rotating = scaling = panning = false;
	if( button==GLUT_LEFT_BUTTON  )
		if( glutGetModifiers() & GLUT_ACTIVE_CTRL ) panning = true;
		else                                        rotating = true;
	else if( button==GLUT_RIGHT_BUTTON )            scaling = true;
}
template< typename Real >
void AuthalicEvolutionVisualization< Real >::motionFunc( int x , int y )
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
void AuthalicEvolutionVisualization< Real >::passiveMotionFunc( int x , int y )
{
}
template< typename Real >
void AuthalicEvolutionVisualization< Real >::idle( void )
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
void AuthalicEvolutionVisualization< Real >::keyboardFunc( unsigned char key , int x , int y )
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
void AuthalicEvolutionVisualization< Real >::specialFunc( int key , int x , int y )
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

template< typename Real >
Point3D< Real > NormalColor( Point3D< Real > n )
{
	Point3D< Real > c = ( -n + Point3D< Real >( 1 , 1 , 1 ) ) * Real(128);
	for( int d=0 ; d<3 ; d++ )
	{
		if( c[d]>Real(255) ) c[d] = Real(255);
		if( c[d]<Real(  0) ) c[d] = Real(  0);
	}
	return c;
}

template< typename Real >
void SetScaleFactors( const std::vector< Point3D< Real > > &paramVertices , const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , std::vector< Real > &scaleFactors )
{
	scaleFactors.resize( vertices.size() );
	static std::vector< Real > paramAreas( vertices.size() );
#pragma omp parallel for
	for( int i=0 ; i<scaleFactors.size() ; i++ ) scaleFactors[i] = paramAreas[i] = 0;
	std::vector< Real > &areas = scaleFactors;
#pragma omp parallel for
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< Real > pv[] = { paramVertices[ triangles[i][0] ] , paramVertices[ triangles[i][1] ] , paramVertices[ triangles[i][2] ] };
		Point3D< Real >  v[] = {      vertices[ triangles[i][0] ] ,      vertices[ triangles[i][1] ] ,      vertices[ triangles[i][2] ] };
		Real pArea = Point3D< Real >::Length( Point3D< Real >::CrossProduct( pv[1]-pv[0] , pv[2]-pv[0] ) ) / 6;
		Real  area = Point3D< Real >::Length( Point3D< Real >::CrossProduct(  v[1]- v[0] ,  v[2]- v[0] ) ) / 6;
		for( int j=0 ; j<3 ; j++ )
		{
#pragma omp atomic
			paramAreas[ triangles[i][j] ] += pArea;
#pragma omp atomic
			areas[ triangles[i][j] ] += area;
		}
	}
#pragma omp parallel for
	for( int i=0 ; i<vertices.size() ; i++ ) scaleFactors[i] = areas[i] / paramAreas[i];
}

template< typename Real >
void SetVertexWeights( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , std::vector< Real > &vertexWeights , std::vector< Real > &vertexTriangleWeights , std::vector< Real > &triangleWeights )
{
	vertexWeights.resize( vertices.size() );
	vertexTriangleWeights.resize( vertices.size() );
	triangleWeights.resize( triangles.size() );
#pragma omp parallel for
	for( int i=0 ; i<vertices.size() ; i++ ) vertexWeights[i] = vertexTriangleWeights[i] = 0;

	Real area = 0;
#pragma omp parallel for
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< Real > v[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
		Real a = Point3D< Real >::Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2;
		triangleWeights[i] = a;
#pragma omp atomic
		area += a;
		for( int j=0 ; j<3 ; j++ )
#pragma omp atomic
			vertexTriangleWeights[ triangles[i][j] ] += a/3;
	}
#pragma omp parallel for
	for( int i=0 ; i<vertices.size() ; i++ ) vertexWeights[i] = vertexTriangleWeights[i] / area;
}

template< typename Real >
Point3D< Real > HSVtoRGB( Point3D< Real > hsv )
{
	// From FvD
	if( hsv[1]<=0.0 ) return Point3D< Real >( hsv[2] , hsv[2] , hsv[2] );

	hsv[0] = (Real)fmod( hsv[0] , 2.0 * M_PI );
	if( hsv[0]<0.0 ) hsv[0] += (Real)( 2.0 * M_PI );
	hsv[0] /= M_PI / 3.0;
	int i = (int)floor(hsv[0]);
	Real f = hsv[0] - i;
	Real p = hsv[2] * (Real)(1. - hsv[1]);
	Real q = hsv[2] * (Real)(1. - (hsv[1]*f));
	Real t = hsv[2] * (Real)(1. - (hsv[1]*(1.-f)));
	Point3D< Real > rgb;
	switch(i) 
	{
	case 0:  rgb[0] = hsv[2] ; rgb[1] = t      ; rgb[2] = p      ; break;
	case 1:  rgb[0] = q      ; rgb[1] = hsv[2] ; rgb[2] = p      ; break;
	case 2:  rgb[0] = p      ; rgb[1] = hsv[2] ; rgb[2] = t      ; break;
	case 3:  rgb[0] = p      ; rgb[1] = q      ; rgb[2] = hsv[2] ; break;
	case 4:  rgb[0] = t      ; rgb[1] = p      ; rgb[2] = hsv[2] ; break;
	default: rgb[0] = hsv[2] ; rgb[1] = p      ; rgb[2] = q      ; break;
	}
	return rgb;
}
template< typename Real >
Point3D< Real > HeatMap( Real value )
{
	// Clamp to the range [-1,1]
	value = std::min< Real >( std::max< Real >( value , -1.f ) , 1.f );
	// Normalize to the range [0,1]
	value = (value+1)/2;

	Point3D< Real > cold( 2. * M_PI , 1. , 1. ) , hot( 4. * M_PI/3 , 1. , 1. );
	return HSVtoRGB( ( cold*(1.-value) + hot*value ) );
}

template< typename Real >
std::vector< TriangleIndex > Triangulate( const std::vector< Point3D< Real > > &vertices , const std::vector< std::vector< int > > &polygons )
{
	std::vector< TriangleIndex > triangles;

	int count = 0;
	for( int i=0 ; i<polygons[i].size() ; i++ ) count += (int)polygons[i].size()-2;
	triangles.reserve( count );
#pragma omp parallel for
	for( int i=0 ; i<polygons.size() ; i++ )
	{
		MinimalAreaTriangulation< Real > MAT;
		const std::vector< int >& poly = polygons[i];
		std::vector< Point3D< Real > > _vertices( poly.size() );
		std::vector< TriangleIndex > _triangles;
		for( int j=0 ; j<poly.size() ; j++ ) _vertices[j] = vertices[ poly[j] ];
		MAT.GetTriangulation( _vertices , _triangles );
		for( int i=0 ; i<_triangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) _triangles[i][j] = poly[ _triangles[i][j] ];
#pragma omp critical
		triangles.insert( triangles.end() , _triangles.begin() , _triangles.end() );
	}
	return triangles;
}
