#ifndef BMP_INCLUDED
#define BMP_INCLUDED


struct BMPInfo
{
	unsigned char* data;
	FILE* fp;
	int width , lineLength;
};

struct BMPReader : public ImageReader
{
	BMPReader( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
	~BMPReader( void );
	unsigned int nextRow( unsigned char* row );
	static bool GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
protected:
	unsigned int _currentRow;
	BMPInfo _info;
};

struct BMPWriter : public ImageWriter
{
	BMPWriter( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality=100 );
	~BMPWriter( void );
	unsigned int nextRow( const unsigned char* row );
	unsigned int nextRows( const unsigned char* row , unsigned int rowNum );
protected:
	BMPInfo _info;
	unsigned int _currentRow;
};

#include "BMP.inl"
#endif //BMP_INCLUDED
