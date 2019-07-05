#define NEW_CODE
//#define W_CYCLE

#include "GL/glew.h"
#include "GL/glut.h"
#include <string>
#include <functional>
#include <limits>
#include <Eigen/Dense>
#include "Misha/Geometry.h"
#include "Misha/Algebra.h"
#include "Misha/Ply.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Timer.h"
#include "Misha/SphericalGeometry.h"
#include "Misha/SphericalSignals.h"
#include "Misha/Visualization.h"
#include "Misha/Camera.h"
#include "Misha/SparseMatrix.h"
#include "Misha/Solver.h"
#include "Misha/SO3Alignment.h"
#include "Misha/OpticalFlow.h"

template< typename Real , bool DivergenceFree >
struct OpticalFlowVisualization : public Visualization::Viewable
{
	struct Stats
	{
		Real min , max , weightSum , valueSum , varianceSum;
		Stats( void ) : weightSum(0) , valueSum(0) , varianceSum(0) , min( std::numeric_limits< Real >::infinity() ) , max( -std::numeric_limits< Real >::infinity() ){}
		Stats &operator += ( const Stats &stats )
		{
			min = std::min< Real >( min , stats.min );
			max = std::max< Real >( max , stats.max );
			weightSum += stats.weightSum;
			valueSum += stats.valueSum;
			varianceSum += stats.varianceSum;
			return *this;
		}
		Stats operator + ( const Stats &stats )
		{
			Stats s = *this;
			s += stats;
			return s;
		}
		Real average( void ) const { return valueSum / weightSum; }
		Real deviation( void ) const { return (Real)sqrt( varianceSum / weightSum ); }
	};
	Stats GetStats( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , const std::vector< Real > &signal );
	Point3D< Real > HSVtoRGB( Point3D< Real > hsv );
	Point3D< Real > HeatMap( Real value );
	Real NormalizeValue( Stats stats , Real v ){ return ( v-stats.average() ) / ( 2*stats.deviation() ); }
	Point3D< float > HeatMapFunction( Stats stats , Real v ){ return Point3D< float >( HeatMap( NormalizeValue( stats , v ) ) ); }

	struct OpticalFlowState
	{
	protected:
		unsigned int _iter;
#ifdef W_CYCLE
		unsigned int _maxLevel;
#endif // W_CYCLE
	public:
		unsigned int maxLevel , solvesPerLevel;
		unsigned int iter( void ) const{ return solvesPerLevel==0 ? 0 : (_iter % solvesPerLevel); }
		unsigned int level( void ) const{ return std::min< unsigned int >( solvesPerLevel==0 ? maxLevel : _iter / solvesPerLevel , maxLevel ); }
		unsigned int levels( void ) const { return maxLevel + 1; }
		// V-Cycle
		void reset( void )
		{
			_iter = 0;
#ifdef W_CYCLE
			_maxLevel = 0;
#endif // W_CYCLE
		}
		bool advance( void )
		{
			_iter++;
#ifdef W_CYCLE
			if( _iter<(_maxLevel+1)*solvesPerLevel ) return true;
			else if( _maxLevel<maxLevel )
			{
				_iter = 0;
				_maxLevel++;
			}
			else return false;
#else // !W_CYCLE
			return _iter<(maxLevel+1)*solvesPerLevel;
#endif // W_CYCLE
		};
	};

protected:
	enum
	{
		MESH_SOURCE ,
		MESH_TARGET ,
		MESH_COUNT
	};
	typename SphericalGeometry::OpticalFlow< Real , DivergenceFree > *_opticalFlow;

	std::vector< TriangleIndex > _triangles[MESH_COUNT];
	std::vector< Point3D< Real > > _colors[MESH_COUNT];
	std::vector< std::vector< Point3D< Real > > > _signalColors[MESH_COUNT];
	std::vector< std::vector< Real > > _vertexData[MESH_COUNT];

	SphericalGeometry::Mesh< Real > _sMeshes[MESH_COUNT];
	std::vector< Point3D< Real > > _vfPositions;
	std::vector< Point3D< Real > > _flowFieldSamples;

	GLuint _VBO[MESH_COUNT] , _EBO[MESH_COUNT];
	Real _vectorScale;
	typename SphericalGeometry::OpticalFlow< Real , DivergenceFree >::Parameters _params;

	int _startResolution;

	OpticalFlowState _opticalFlowState;

	void _setElementBufferObject( void );
	void _setVertexBufferObject( void );
	void _updateVertexBufferObject( void );

	void _setFibonacci( unsigned int samples );
	void _initFibonacci( unsigned int samples );
	void _setFibonacci( void );

	unsigned int _mesh( int x ) const;

public:
	bool verbose;
	int animationSteps;
	Camera camera;
	float zoom;
	int oldX , oldY , newX , newY;
	bool scaling , rotating , panning , useLight , edgeMode , showVectors , showColor , showGrid;

	GLfloat lightAmbient[4] , lightDiffuse[4] , lightSpecular[4] , shapeSpecular[4] , shapeSpecularShininess;
	char resolutionStr[1024];
	char statsStr[1024];

	void readCamera( const char *fileName );
	void writeCamera( const char *fileName ) const;

	OpticalFlowVisualization( const SphericalGeometry::Mesh< Real > &sourceSMesh , const std::vector< Point3D< Real > > &sourceColors , const std::vector< Real > *sourceValues , const SphericalGeometry::Mesh< Real > &targetSMesh , const std::vector< Point3D< Real > > &targetColors , const std::vector< Real > *targetValues , unsigned int levels , Real bandLimitDampening , Real sWeight , int startResolution , int integrationSamples , Real stepSize , unsigned int subSteps , bool rescaleFlow , unsigned int solvesPerLevel , bool verbose );
	~OpticalFlowVisualization( void );

	bool advance( typename SphericalGeometry::OpticalFlow< Real , DivergenceFree >::Stats &stats , bool updateVBO );
	const std::vector< Point3D< Real > > &advectedSourceVertices( void ) const { return _sMeshes[ MESH_SOURCE ].vertices; }
	const std::vector< Point3D< Real > > &advectedTargetVertices( void ) const { return _sMeshes[ MESH_TARGET ].vertices; }

	void idle( void );
	void keyboardFunc( unsigned char key , int x , int y );
	void specialFunc( int key, int x, int y );
	void setProjectionMatrix( unsigned int modelIndex );
	void setModelviewMatrix( unsigned int modelIndex );
	void display( void );
	void mouseFunc( int button , int state , int x , int y );
	void motionFunc( int x , int y );
	void passiveMotionFunc( int x , int y );

	static void OutputCamera( Visualization::Viewable* v , const char *prompt )
	{ 
		OpticalFlowVisualization *ofv = (OpticalFlowVisualization*)v;
		ofv->writeCamera( prompt );
	}
	static void SampleScaleUpCallBack( Visualization::Viewable* v , const char* ){ ( (OpticalFlowVisualization*)v)->_setFibonacci( (int)ceil( ( (OpticalFlowVisualization*)v)->_vfPositions.size()*1.1 ) ); }
	static void SampleScaleDownCallBack( Visualization::Viewable* v , const char* ){ ( (OpticalFlowVisualization*)v)->_setFibonacci( (int)ceil( ( (OpticalFlowVisualization*)v)->_vfPositions.size()/1.1 ) ); }
	static void ToggleVectorsCallBack( Visualization::Viewable* v , const char* ){ ( (OpticalFlowVisualization*)v)->showVectors = !( (OpticalFlowVisualization*)v)->showVectors; }
	static void ToggleColorsCallBack( Visualization::Viewable* v , const char* )
	{
		OpticalFlowVisualization *ofv = (OpticalFlowVisualization*)v;
		ofv->showColor = !ofv->showColor;
		ofv->_updateVertexBufferObject();
	}
	static void ToggleGridCallBack( Visualization::Viewable* v , const char* ){ ( (OpticalFlowVisualization*)v)->showGrid = !( (OpticalFlowVisualization*)v)->showGrid; }
	static void ToggleEdgesCallBack( Visualization::Viewable* v , const char* ){ ( (OpticalFlowVisualization*)v)->edgeMode = !( (OpticalFlowVisualization*)v)->edgeMode; }
	static void ToggleLightCallBack( Visualization::Viewable* v , const char* ){ ( (OpticalFlowVisualization*)v)->useLight = !( (OpticalFlowVisualization*)v)->useLight; }
	static void AdvanceCallBack( Visualization::Viewable* v , const char* ){ ( (OpticalFlowVisualization*)v)->animationSteps++; }
	static void ToggleAnimationCallBack( Visualization::Viewable* v , const char* ){ ( (OpticalFlowVisualization*)v)->animationSteps = ( (OpticalFlowVisualization*)v)->animationSteps ? 0 : -1; }
	static void VectorScaleUpCallBack( Visualization::Viewable* v , const char* ){ ( (OpticalFlowVisualization*)v)->_vectorScale *= 1.1f; }
	static void VectorScaleDownCallBack( Visualization::Viewable* v , const char* ){ ( (OpticalFlowVisualization*)v)->_vectorScale /= 1.1f; }
};

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::readCamera( const char *fileName )
{
	FILE *fp = fopen( fileName , "r" );
	if( !fp ) fprintf( stderr , "[WARNING] Could not open camera file for reading: %s\n" , fileName );
	else
	{
		camera.read( fp );
		if( fscanf( fp , " %f" , &zoom )!=1 ) fprintf( stderr , "[ERROR] Could not read camera zoom from: %s\n" , fileName );
		fclose( fp );
	}
}
template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::writeCamera( const char *fileName ) const
{
	FILE *fp = fopen( fileName , "w" );
	if( !fp ) fprintf( stderr , "[WARNING] Could not open camera file for writing: %s\n" , fileName );
	else
	{
		camera.write( fp );
		fprintf( fp , " %f\n" , zoom );
		fclose( fp );
	}
}

template< typename Real , bool DivergenceFree >
unsigned int OpticalFlowVisualization< Real , DivergenceFree >::_mesh( int x ) const
{
	return x / ( screenWidth/MESH_COUNT );
}

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::_initFibonacci( unsigned int samples )
{
	static const Real PHI = (Real)( ( 1.+sqrt(5.) ) / 2. );
	static const Real PHI_INV = (Real)( 1. / PHI );
	_vfPositions.resize( samples );
	_flowFieldSamples.resize( samples );
#pragma omp parallel for
	for( int i=0 ; i<(int)samples ; i++ )
	{
		Real theta = (Real)acos( 1 - ( (Real)(2*i+1) ) / samples );
		Real phi = PHI_INV * 2 * i * M_PI;
		_vfPositions[i] = Point3D< Real >( cos(phi)*sin(theta) , sin(phi)*sin(theta) , cos(theta) );
	}
}

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::_setFibonacci( unsigned int samples )
{
	_initFibonacci( samples );
	_setFibonacci();
}

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::_setFibonacci( void )
{
#pragma omp parallel for
	for( int i=0 ; i<(int)_vfPositions.size() ; i++ ) _flowFieldSamples[i] = _opticalFlow->flowField( _vfPositions[i] );
}

template< typename Real , bool DivergenceFree >
bool OpticalFlowVisualization< Real , DivergenceFree >::advance( typename SphericalGeometry::OpticalFlow< Real , DivergenceFree >::Stats &stats , bool updateVBO )
{
	if( _opticalFlowState.solvesPerLevel==0 ) return false;
	double time = 0;
	Timer timer;

	if( _opticalFlowState.iter()==0 )
	{
		if( _opticalFlow ) delete _opticalFlow;
		_params.resolution = _startResolution<<_opticalFlowState.level();
		_opticalFlow = new typename SphericalGeometry::OpticalFlow< Real , DivergenceFree >( _sMeshes[0] , _sMeshes[1] , _vertexData[0][ _opticalFlowState.level() ] , _vertexData[1][ _opticalFlowState.level() ] , _params );
	}
	_opticalFlow->advance();
	time = timer.elapsed();

	OpticalFlowState oldState = _opticalFlowState;
	bool moreToDo = _opticalFlowState.advance();

	// Rasterize the advected signals
	stats = _opticalFlow->getStats();
	if( updateVBO )
	{
		_setFibonacci();
		_updateVertexBufferObject();
		sprintf( statsStr , "Set flow field / Advected / Set Signals[%d](%d): %.2f(s)\n" , _startResolution<<oldState.level() , oldState.iter()+1 , time );
	}
	if( verbose ) printf( "Set flow field / Advected / Set Signals[%d](%d): %.2f(s)\n" , _startResolution<<oldState.level() , oldState.iter()+1 , time );

	if( !moreToDo ) _opticalFlowState.reset() , animationSteps = 0;
	return moreToDo;
}

template< typename Real , bool DivergenceFree >
OpticalFlowVisualization< Real , DivergenceFree >::OpticalFlowVisualization( const SphericalGeometry::Mesh< Real > &sourceSMesh , const std::vector< Point3D< Real > > &sourceColors , const std::vector< Real > *sourceValues , const SphericalGeometry::Mesh< Real > &targetSMesh , const std::vector< Point3D< Real > > &targetColors , const std::vector< Real > *targetValues , unsigned int levels , Real bandLimitDampening , Real sWeight , int startResolution , int integrationSamples , Real stepSize , unsigned int subSteps , bool rescaleFlow , unsigned int solvesPerLevel , bool verbose )
{
	this->verbose = verbose;
	_opticalFlow = NULL;
	animationSteps = 0;
	_vectorScale = 1;
	edgeMode = false;
	screenHeight = 512;
	showColor = false;
	useLight = true;
	showVectors = true;
	showGrid = true;
	screenWidth = 2*screenHeight;

	_opticalFlowState.maxLevel = levels-1;
	_opticalFlowState.solvesPerLevel = solvesPerLevel;
	_opticalFlowState.reset();

	_params.subSteps = subSteps;
	_params.rescaleFlow = false;
	_params.bandLimitDampening = bandLimitDampening;
	_params.smoothingWeight = sWeight;
	_params.integrationSamples = integrationSamples;
	_params.stepSize = stepSize;
	_params.useSemiImplicit = true;
	_params.rescaleFlow = rescaleFlow;
	_startResolution = startResolution;

	_sMeshes[MESH_SOURCE] = sourceSMesh;
	_sMeshes[MESH_TARGET] = targetSMesh;
	_triangles[MESH_SOURCE] = _sMeshes[MESH_SOURCE].triangulate();
	_triangles[MESH_TARGET] = _sMeshes[MESH_TARGET].triangulate();

	// Copy the per-vertex colors
	_colors[MESH_SOURCE] = sourceColors;
	_colors[MESH_TARGET] = targetColors;

	for( int m=0 ; m<MESH_COUNT ; m++ ) _signalColors[m].resize( levels ) , _vertexData[m].resize( levels );

	// Copy the per-vertex signals
	for( unsigned int d=0 ; d<levels ; d++ )
	{
		_vertexData[MESH_SOURCE][d] = sourceValues[d];
		_vertexData[MESH_TARGET][d] = targetValues[d];
	}

	// Transform the per-vertex signals into per-vertex colors
	for( unsigned int l=0 ; l<levels ; l++ )
	{
		Stats stats;
		for( int m=0 ; m<MESH_COUNT ; m++ ) stats = GetStats( _sMeshes[m].vertices , _triangles[m] , _vertexData[m][l] );
		for( int m=0 ; m<MESH_COUNT ; m++ ) 
		{
			_signalColors[m][l].resize( _vertexData[m][l].size() );
			for( int j=0 ; j<_signalColors[m][l].size() ; j++ ) _signalColors[m][l][j] = HeatMapFunction( stats , _vertexData[m][l][j] ) * 255.f;
		}
	}

	_initFibonacci( 1024 );

	////////////////////////////

	camera.position = Point3D< double >( 0 , 0 , 2. ) , zoom = 1.f / 0.95f;
	for( int i=0 ; i<MESH_COUNT ; i++ ) _VBO[i] = _EBO[i] = 0;

	lightAmbient [0] = lightAmbient [1] = lightAmbient [2] = 0.25f , lightAmbient [3] = 1.f;
	lightDiffuse [0] = lightDiffuse [1] = lightDiffuse [2] = 0.70f , lightDiffuse [3] = 1.f;
	lightSpecular[0] = lightSpecular[1] = lightSpecular[2] = 0.25f , lightSpecular[3] = 1.f;
	shapeSpecular[0] = shapeSpecular[1] = shapeSpecular[2] = 1.00f , shapeSpecular[3] = 1.f;
	shapeSpecularShininess = 128;

	rotating = scaling = panning = false;
	callBacks.push_back( Visualization::KeyboardCallBack( this , 'l' , "toggle light" , ToggleLightCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , 'c' , "toggle colors" , ToggleColorsCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , 'e' , "toggle edges" , ToggleEdgesCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , 'v' , "toggle vectors" , ToggleVectorsCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , 'g' , "toggle grid" , ToggleGridCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , '+' , "advance" , AdvanceCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , ' ' , "advance" , ToggleAnimationCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , '{' , "sample down" , SampleScaleDownCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , '}' , "sample up" , SampleScaleUpCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , '[' , "scale down" , VectorScaleDownCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , ']' , "scale up" , VectorScaleUpCallBack ) );
	callBacks.push_back( Visualization::KeyboardCallBack( this , 'C' , "output camera" , "Filename" , OutputCamera ) );
	info.push_back( resolutionStr );
	info.push_back( statsStr );

	sprintf( resolutionStr , "Vertices / Polygons / Spherical Resolution: %d , %d / %d , %d / %d" , (int)_sMeshes[0].vertices.size() , (int)_sMeshes[1].vertices.size() , (int)_sMeshes[0].polygons.size() , (int)_sMeshes[1].polygons.size() , (int)( _startResolution<<(_opticalFlowState.levels()-1) ) );
	sprintf( statsStr , "" );
}

template< typename Real , bool DivergenceFree >
OpticalFlowVisualization< Real , DivergenceFree >::~OpticalFlowVisualization( void )
{
	if( _opticalFlow ) delete _opticalFlow;
	_opticalFlow = NULL;
}

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::_setElementBufferObject( void )
{
	glGenBuffers( MESH_COUNT , _EBO );
	for( int i=0 ; i<MESH_COUNT ; i++ )
	{
		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , _EBO[i] );
		glBufferData( GL_ELEMENT_ARRAY_BUFFER , _triangles[i].size() * sizeof( int ) * 3 ,  &_triangles[i][0] , GL_STATIC_DRAW );
		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , 0 );
	}
}

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::_setVertexBufferObject( void )
{
	glGenBuffers( MESH_COUNT , _VBO );
	for( unsigned int m=0 ; m<MESH_COUNT ; m++ )
	{
		Point3D< float > *vertexData = new Point3D< float >[ _sMeshes[m].vertices.size()*3 ];
		glBindBuffer( GL_ARRAY_BUFFER , _VBO[m] );
		for( int i=0 ; i<_sMeshes[m].vertices.size() ; i++ ) vertexData[i] = vertexData[ _sMeshes[m].vertices.size() + i ] = Point3D< float >( _sMeshes[m].vertices[i] );
		if( showColor ) for( int i=0 ; i<_sMeshes[m].vertices.size() ; i++ ) vertexData[ _sMeshes[m].vertices.size()*2 + i ] = Point3D< float >( _colors[m][i]/255.f );
		else            for( int i=0 ; i<_sMeshes[m].vertices.size() ; i++ ) vertexData[ _sMeshes[m].vertices.size()*2 + i ] = Point3D< float >( _signalColors[m][ _opticalFlowState.level() ][i]/255.f );
		glBufferData( GL_ARRAY_BUFFER , 9 * _sMeshes[m].vertices.size() * sizeof( float ) , (float*)vertexData , GL_DYNAMIC_DRAW );
		delete[] vertexData;
	}
}

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::_updateVertexBufferObject( void )
{
	for( unsigned int m=0 ; m<MESH_COUNT ; m++ )
	{
		Point3D< float > *vertexData = new Point3D< float >[ _sMeshes[m].vertices.size()*3 ];
		glBindBuffer( GL_ARRAY_BUFFER , _VBO[m] );
		for( int i=0 ; i<_sMeshes[m].vertices.size() ; i++ ) vertexData[i] = vertexData[ _sMeshes[m].vertices.size() + i ] = Point3D< float >( _sMeshes[m].vertices[i] );
		if( showColor ) for( int i=0 ; i<_sMeshes[m].vertices.size() ; i++ ) vertexData[ _sMeshes[m].vertices.size()*2 + i ] = Point3D< float >( _colors[m][i]/255.f );
		else            for( int i=0 ; i<_sMeshes[m].vertices.size() ; i++ ) vertexData[ _sMeshes[m].vertices.size()*2 + i ] = Point3D< float >( _signalColors[m][ _opticalFlowState.level() ][i]/255.f );
		glBufferSubData( GL_ARRAY_BUFFER , 0 , 9 * _sMeshes[m].vertices.size() * sizeof( float ) , (float*)vertexData );
		delete[] vertexData;
	}
}

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::setProjectionMatrix( unsigned int mesh )
{
	int width = screenWidth/2 , height = screenHeight;
	switch( mesh )
	{
	case MESH_SOURCE: glViewport( 0 , 0 , width , height ) ; break;
	case MESH_TARGET: glViewport( width , 0 , width , height ) ; break;
	}

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();

	float ar = (float)width/(float)height , ar_r = 1.f/ar;
	if( width>height ) glOrtho( -ar*zoom , ar*zoom , -zoom , zoom , -2.f , 2.f );
	else               glOrtho( -zoom , zoom , -ar_r*zoom , ar_r*zoom , -2.f , 2.f );
}

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::setModelviewMatrix( unsigned int )
{
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	camera.draw();
}

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::display( void )
{
	for( unsigned int mesh=0 ; mesh<MESH_COUNT ; mesh++ )
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
		else glPolygonMode( GL_FRONT_AND_BACK , GL_FILL );
		if( !_EBO[mesh] || !_EBO[mesh] ) _setElementBufferObject();
		if( !_VBO[mesh] || !_VBO[mesh] ) _setVertexBufferObject();
		setProjectionMatrix( mesh );
		setModelviewMatrix( mesh );

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

		glBindBuffer( GL_ARRAY_BUFFER , _VBO[mesh] );
		glEnableClientState( GL_VERTEX_ARRAY );
		glEnableClientState( GL_NORMAL_ARRAY );
		glEnableClientState( GL_COLOR_ARRAY );
		glVertexPointer  ( 3 , GL_FLOAT , 0 , (GLubyte*)NULL + sizeof( float ) * _sMeshes[mesh].vertices.size() * 0 );
		glNormalPointer  (     GL_FLOAT , 0 , (GLubyte*)NULL + sizeof( float ) * _sMeshes[mesh].vertices.size() * 3 );
		glColorPointer   ( 3 , GL_FLOAT , 0 , (GLubyte*)NULL + sizeof( float ) * _sMeshes[mesh].vertices.size() * 6 );
		glColor3f( 0.75f , 0.75f , 0.75f );
		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , _EBO[mesh] );
		glDrawElements( GL_TRIANGLES , (GLsizei)( _triangles[mesh].size() * 3 ) , GL_UNSIGNED_INT , NULL );
		glDisableClientState( GL_NORMAL_ARRAY );
		glDisableClientState( GL_COLOR_ARRAY );
		glBindBuffer( GL_ARRAY_BUFFER , 0 );
		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , 0 );
		glDepthMask( GL_TRUE );


		auto DrawVectorField = [&]( const std::vector< Point3D< Real > > &vf , Real scale )
		{
			glDisable( GL_LIGHTING );
			glBegin( GL_TRIANGLES );
			for( int i=0 ; i<(int)_vfPositions.size() ; i++ )
			{
				Point3D< float > p = Point3D< float >( _vfPositions[i] );
				Point3D< float > d = Point3D< float >( vf[i] * scale );
				Point3D< float > _n = Point3D< float >::CrossProduct( d , p );
				_n /= 20;
				glColor3f( 1.0f , 1.0f , 1.0f ) , glVertex3f( p[0]+d[0] , p[1]+d[1] , p[2]+d[2] );
				glColor3f( 0.0f , 0.0f , 0.0f ) , glVertex3f( p[0]-_n[0] , p[1]-_n[1] , p[2]-_n[2] ) , glVertex3f( p[0]+_n[0] , p[1]+_n[1] , p[2]+_n[2] );
			}
			glEnd();
		};

		if( showVectors )
		{
			glPolygonMode( GL_FRONT_AND_BACK , GL_FILL );
			if( true )
			{
				if     ( mesh==MESH_SOURCE ) DrawVectorField( _flowFieldSamples ,  _vectorScale );
				else if( mesh==MESH_TARGET ) DrawVectorField( _flowFieldSamples , -_vectorScale );
			}
		}

		if( showGrid )
		{
			const int GridLines = 4;
			unsigned int width = screenWidth/2 , height = screenHeight;
			glMatrixMode( GL_PROJECTION );
			glPushMatrix();
			glLoadIdentity();
			if( width>height ) glOrtho( -1.f , 1.f , -(float)height/width , (float)height/width , 0.f , 1.f );
			else               glOrtho( -(float)width/height , (float)width/height , -1.f , 1.f , 0.f , 1.f );
			glMatrixMode( GL_MODELVIEW );
			glPushMatrix();
			glLoadIdentity();

			glDisable( GL_DEPTH_TEST );
			glDisable( GL_LIGHTING );
			glColor3f( 0.25f , 0.25f , 0.25f );
			glLineWidth( 1 );
			glBegin( GL_LINES );
			for( int i=-GridLines ; i<=GridLines ; i++ )
			{
				glVertex2f( (float)i/GridLines , -1.f ) , glVertex2f( (float)i/GridLines , 1.f );
				glVertex2f( -1.f , (float)i/GridLines ) , glVertex2f( 1.f , (float)i/GridLines );
			}
			glEnd();
			glEnable( GL_DEPTH_TEST );
			glMatrixMode( GL_PROJECTION ) ; glPopMatrix();
			glMatrixMode( GL_MODELVIEW  ) ; glPopMatrix();
		}
	}

	glViewport( 0 , 0 , screenWidth , screenHeight );

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
		glLineWidth( 4 );
		glBegin( GL_LINES );
		glVertex2i( screenWidth/2 , 0 ) , glVertex2i( screenWidth/2 , screenHeight );
		glEnd();
		glEnable( GL_DEPTH_TEST );
		glMatrixMode( GL_PROJECTION ) ; glPopMatrix();
		glMatrixMode( GL_MODELVIEW  ) ; glPopMatrix();
	}
}

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::mouseFunc( int button , int state , int x , int y )
{
	newX = x , newY = y;

	rotating = scaling = panning = false;
	if( button==GLUT_LEFT_BUTTON  )
		if( glutGetModifiers() & GLUT_ACTIVE_CTRL ) panning = true;
		else                                        rotating = true;
	else if( button==GLUT_RIGHT_BUTTON )            scaling = true;
}

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::motionFunc( int x , int y )
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

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::passiveMotionFunc( int x , int y )
{
}

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::idle( void )
{
	if( !promptCallBack )
	{
		if( animationSteps )
		{
			typename SphericalGeometry::template OpticalFlow< Real , DivergenceFree >::Stats stats;
			advance( stats , true );
			if( animationSteps>0 ) animationSteps--;
			glutPostRedisplay();
		}
	}
}

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::keyboardFunc( unsigned char key , int x , int y )
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

template< typename Real , bool DivergenceFree >
void OpticalFlowVisualization< Real , DivergenceFree >::specialFunc( int key , int x , int y )
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
///////////////////////////////////////////////////////////////////////////////

template< typename Real , bool DivergenceFree >
typename OpticalFlowVisualization< Real , DivergenceFree >::Stats OpticalFlowVisualization< Real , DivergenceFree >::GetStats( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , const std::vector< Real > &signal )
{
	std::vector< Real > vWeights , vtWeights , tWeights;
	{
		vWeights.resize( vertices.size() );
		vtWeights.resize( vertices.size() );
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ ) vWeights[i] = 0;

		Real area = 0;
#pragma omp parallel for reduction( + : area )
		for( int i=0 ; i<triangles.size() ; i++ )
		{
			Point3D< Real > v[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
			Real a = Point3D< Real >::Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2;
			area += a;
			for( int j=0 ; j<3 ; j++ )
#pragma omp atomic
				vtWeights[ triangles[i][j] ] += a/3;
		}
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ ) vWeights[i] = vtWeights[i] / area;
	}

	Stats stats;

	stats.min = stats.max = signal[0];
	for( int i=0 ; i<signal.size() ; i++ )
	{
		stats.min = std::min< Real >( stats.min , signal[i] ) , stats.max = std::max< Real >( stats.max , signal[i] );
		stats.weightSum += vWeights[i];
		stats.valueSum += signal[i] * vWeights[i];
	}
	for( int i=0 ; i<signal.size() ; i++ ) stats.varianceSum += ( signal[i] - stats.average() ) * ( signal[i] - stats.average() ) * vWeights[i];

	return stats;
}

template< typename Real , bool DivergenceFree >
Point3D< Real > OpticalFlowVisualization< Real , DivergenceFree >::HSVtoRGB( Point3D< Real > hsv )
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

template< typename Real , bool DivergenceFree >
Point3D< Real > OpticalFlowVisualization< Real , DivergenceFree >::HeatMap( Real value )
{
	// Clamp to the range [-1,1]
	value = std::min< Real >( std::max< Real >( value , -1.f ) , 1.f );
	// Normalize to the range [0,1]
	value = (value+1)/2;

	Point3D< Real > cold( 2. * M_PI , 1. , 1. ) , hot( 4. * M_PI/3 , 1. , 1. );
	return HSVtoRGB( ( cold*(1.-value) + hot*value ) );
}
