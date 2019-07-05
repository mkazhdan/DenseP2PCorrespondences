#ifndef SPARSE_MATRIX_INTERFACE_INCLUDED
#define SPARSE_MATRIX_INTERFACE_INCLUDED

#define FORCE_TWO_BYTE_ALIGNMENT 1
#include "Array.h"
#include <vector>


#if FORCE_TWO_BYTE_ALIGNMENT
#pragma pack(push)
#pragma pack(2)
#endif // FORCE_TWO_BYTE_ALIGNMENT
template< class T , class IndexType >
struct MatrixEntry
{
	MatrixEntry( void )             { N =-1 , Value = 0; }
	MatrixEntry( IndexType i )      { N = i , Value = 0; }
	MatrixEntry( IndexType n , T v ){ N = n , Value = v; }
	IndexType N;
	T Value;
};
#if FORCE_TWO_BYTE_ALIGNMENT
#pragma pack(pop)
#endif // FORCE_TWO_BYTE_ALIGNMENT


enum
{
	MULTIPLY_ADD = 1 ,
	MULTIPLY_NEGATE = 2
};

template< class T , class const_iterator > class SparseMatrixInterface
{
public:
	virtual const_iterator begin( size_t row ) const = 0;
	virtual const_iterator end  ( size_t row ) const = 0;
	virtual size_t Rows   ( void )             const = 0;
	virtual size_t RowSize( size_t idx )       const = 0;

	size_t Entries( void ) const;

	double SquareNorm( void ) const;
	double SquareASymmetricNorm( void ) const;
	double SquareASymmetricNorm( int& idx1 , int& idx2 ) const;

	template< class T2 > void Multiply      (           ConstPointer( T2 )  In , Pointer( T2 ) Out , int multiplyFlag=0 ) const;
	template< class T2 > void MultiplyScaled( T scale , ConstPointer( T2 )  In , Pointer( T2 ) Out , int multiplyFlag=0 ) const;
	template< class T2 > void Multiply      (                Pointer( T2 )  In , Pointer( T2 ) Out , int multiplyFlag=0 ) const { Multiply      (         ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }
	template< class T2 > void MultiplyScaled( T scale ,      Pointer( T2 )  In , Pointer( T2 ) Out , int multiplyFlag=0 ) const { MultiplyScaled( scale , ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }

	template< class T2 > void SetDiagonal( Pointer( T2 ) diagonal ) const;
	template< class T2 > void JacobiIteration( ConstPointer( T2 ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , Pointer( T2 ) Mx , T2 sor ) const;
	template< class T2 > void JacobiIteration( ConstPointer( T2 ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , T2 sor ) const { Pointer( T2 ) Mx = AllocPointer< T2 >( Rows() ) ; JacobiIteration( diagonal , b , x , Mx , sor ) ; FreePointer( Mx ); }
#if 1
	template< class T2 , bool StripDiagonal=false > void GSIteration(                                                        ConstPointer( T2 ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward ) const;
	template< class T2 , bool StripDiagonal=false > void GSIteration( std::vector< std::vector< int > >& multiColorIndices , ConstPointer( T2 ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward ) const;
#else
	template< class T2 > void GSIteration(                                                        ConstPointer( T2 ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward ) const;
	template< class T2 > void GSIteration( std::vector< std::vector< int > >& multiColorIndices , ConstPointer( T2 ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward ) const;

#endif
};
// Assuming that the SPDOperator class defines:
//		int SPDOperator::Rows( void ) const
//		auto SPDOperator::Multiply( ConstPointer( T ) , Pointer( T ) ) const
template< class SPDOperator , class T > int SolveCG( const SPDOperator& M , ConstPointer( T ) b , int iters , Pointer( T ) x , T eps=1e-8 , bool solveNormal=false );

#include "SparseMatrixInterface.inl"
#endif // SPARSE_MATRIX_INTERFACE_INCLUDED
