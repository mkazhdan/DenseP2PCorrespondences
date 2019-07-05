#include <stdio.h>
#include <stdlib.h>


/* constants for the biCompression field */
#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

/* Some magic numbers */

#define BMP_BF_TYPE 0x4D42
/* word BM */

#define BMP_BF_OFF_BITS 54
/* 14 for file header + 40 for info header (not sizeof(), but packed size) */

#define BMP_BI_SIZE 40
/* packed size of info header */

#ifndef _WIN32
typedef struct tagBITMAPFILEHEADER
{
	unsigned short int bfType;
	unsigned int bfSize;
	unsigned short int bfReserved1;
	unsigned short int bfReserved2;
	unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
	unsigned int biSize;
	int biWidth;
	int biHeight;
	unsigned short int biPlanes;
	unsigned short int biBitCount;
	unsigned int biCompression;
	unsigned int biSizeImage;
	int biXPelsPerMeter;
	int biYPelsPerMeter;
	unsigned int biClrUsed;
	unsigned int biClrImportant;
} BITMAPINFOHEADER;
#endif // !_WIN32


inline bool BMPReader::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	BITMAPFILEHEADER bmfh;
	BITMAPINFOHEADER bmih;

	FILE* fp = fopen( fileName , "rb" );
	if( !fp ) fprintf( stderr , "Failed to open: %s\n" , fileName ) , exit(0);

	fread( &bmfh , sizeof( BITMAPFILEHEADER ) , 1 , fp );
	fread( &bmih , sizeof( BITMAPINFOHEADER ) , 1 , fp );

	if( bmfh.bfType!=BMP_BF_TYPE || bmfh.bfOffBits!=BMP_BF_OFF_BITS ) fprintf( stderr , "[ERROR] BMPReader::GetInfo: Bad bitmap file header\n" ) , fclose( fp ) , exit( 0 );
	if( bmih.biSize!=BMP_BI_SIZE || bmih.biWidth<=0 || bmih.biHeight<=0 || bmih.biPlanes!=1 || bmih.biBitCount!=24 || bmih.biCompression!=BI_RGB ) fprintf( stderr , "[ERROR] BMPReader::GetInfo: Bad bitmap file info\n" ) , fclose( fp ) , exit( 0 );
	width           = bmih.biWidth;
	height          = bmih.biHeight;
	channels        = 3;
	int lineLength = width * channels;
	if( lineLength % 4 ) lineLength = (lineLength / 4 + 1) * 4;
	if( bmih.biSizeImage!=lineLength*height ) fprintf( stderr , "[ERROR] BMPReader::GetInfo: Bad bitmap image size\n" ) , fclose( fp ) , exit( 0 );
	fclose( fp );
	return true;
}

inline BMPReader::BMPReader( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	_currentRow = 0;
	BITMAPFILEHEADER bmfh;
	BITMAPINFOHEADER bmih;

	_info.fp = fopen( fileName , "rb" );
	_info.data = NULL;
	if( !_info.fp ) fprintf( stderr , "[ERROR] BMPInitRead: Failed to open: %s\n" , fileName ) , exit(0);

	fread( &bmfh , sizeof( BITMAPFILEHEADER ) , 1 , _info.fp );
	fread( &bmih , sizeof( BITMAPINFOHEADER ) , 1 , _info.fp );

	if( bmfh.bfType!=BMP_BF_TYPE || bmfh.bfOffBits!=BMP_BF_OFF_BITS ) fprintf( stderr , "[ERROR] BMPReader::BMPReader: Bad bitmap file header\n" ) , exit( 0 );
	if( bmih.biSize!=BMP_BI_SIZE || bmih.biWidth<=0 || bmih.biHeight<=0 || bmih.biPlanes!=1 || bmih.biBitCount!=24 || bmih.biCompression!=BI_RGB ) fprintf( stderr , "[ERROR] BMPInitReadColor: Bad bitmap file info\n" ) , exit( 0 );

	_info.width = width = bmih.biWidth;
	height = bmih.biHeight;
	channels = 3;
	_info.lineLength = width * 3;
	if( _info.lineLength % 4 ) _info.lineLength = (_info.lineLength / 4 + 1) * 4;
	if( bmih.biSizeImage!=_info.lineLength*height ) fprintf( stderr , "[ERROR] BMPReader::BMPReader: Bad bitmap image size\n" ) , exit( 0 );

	fseek( _info.fp , (long) bmfh.bfOffBits , SEEK_SET );
	fseek( _info.fp , (long) _info.lineLength * height , SEEK_CUR );
}
inline BMPReader::~BMPReader( void ){ fclose( _info.fp ); }
inline unsigned int BMPReader::nextRow( unsigned char* row )
{
	fseek( _info.fp , -_info.lineLength , SEEK_CUR );
	fread( row , 1 , _info.lineLength , _info.fp );
	fseek( _info.fp , -_info.lineLength , SEEK_CUR );
	if( ferror(_info.fp) ) fprintf( stderr , "[ERROR] BMPReader::nextRow: Error reading bitmap row\n" ) , exit( 0 );
	for( int i=0 ; i<_info.width ; i++ ) { unsigned char temp = row[i*3] ; row[i*3] = row[i*3+2] ; row[i*3+2] = temp; }
	return _currentRow++;
}

inline BMPWriter::BMPWriter( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int )
{
	_currentRow = 0;

	BITMAPFILEHEADER bmfh;
	BITMAPINFOHEADER bmih;

	_info.fp = fopen( fileName , "wb" );
	if( !_info.fp ) fprintf( stderr , "BMPWriter::BMPWriter: Failed to open: %s\n" , fileName ) , exit(0);
	_info.width = width;

	_info.lineLength = width * 3;	/* RGB */
	if( _info.lineLength % 4 ) _info.lineLength = (_info.lineLength / 4 + 1) * 4;
	/* Write file header */

	bmfh.bfType = BMP_BF_TYPE;
	bmfh.bfSize = BMP_BF_OFF_BITS + _info.lineLength * height;
	bmfh.bfReserved1 = 0;
	bmfh.bfReserved2 = 0;
	bmfh.bfOffBits = BMP_BF_OFF_BITS;

	fwrite( &bmfh , sizeof(BITMAPFILEHEADER) , 1 , _info.fp );

	bmih.biSize = BMP_BI_SIZE;
	bmih.biWidth = width;
	bmih.biHeight = -(int)height;
	bmih.biPlanes = 1;
	bmih.biBitCount = 24;			/* RGB */
	bmih.biCompression = BI_RGB;	/* RGB */
	bmih.biSizeImage = _info.lineLength * (unsigned int) bmih.biHeight;	/* RGB */
	bmih.biXPelsPerMeter = 2925;
	bmih.biYPelsPerMeter = 2925;
	bmih.biClrUsed = 0;
	bmih.biClrImportant = 0;

	fwrite( &bmih , sizeof(BITMAPINFOHEADER) , 1 , _info.fp );

	_info.data = (unsigned char*)malloc( _info.lineLength * sizeof( unsigned char ) );
	if( !_info.data ) fprintf( stderr , "[ERROR] BMPWriter::BMPWriter: Could not allocate memory for bitmap data\n" ) , exit( 0 );
}
inline BMPWriter::~BMPWriter( void )
{
	free( _info.data );
	fclose( _info.fp );
}
inline unsigned int BMPWriter::nextRow( const unsigned char* row )
{
	for( int i=0 ; i<_info.width ; i++ )
	{
		_info.data[i*3+0] = row[i*3+2];
		_info.data[i*3+1] = row[i*3+1];
		_info.data[i*3+2] = row[i*3+0];
	}
	fwrite( _info.data , sizeof(unsigned char) , _info.width*3 , _info.fp );
	int nbytes = _info.width*3;
	while( nbytes % 4 ) putc( 0 , _info.fp ) , nbytes++;
	return ++_currentRow;
}
inline unsigned int BMPWriter::nextRows( const unsigned char* rows , unsigned int rowNum )
{
	for( unsigned int i=0 ; i<rowNum ; i++ ) nextRow( rows + _info.width*i );
	return _currentRow;
}
