/*
Copyright (c) 2013, Michael Kazhdan
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

#ifndef ALGEBRA_INCLUDED
#define ALGEBRA_INCLUDED

template<class Element>
class Group
{
public:
	// For this to work, need to define:
	// void Element::SetIdentity	(void);
	// void Element::Multiply		(const Element& e)
	// void Element::Invert			(void)

	static Element Identity(void)
	{
		Element out;
		out.SetIdentity();
		return out;
	}
	Element inverse(void) const
	{
		Element out=*(Element*)this;
		out.Invert();
		return out;
	}
	friend Element  operator *  (const Element& e1,const Element& e2)
	{
		Element out=e1;
		out.Multiply(e2);
		return out;
	}
	friend Element operator /  (const Element& e1,const Element& e2)
	{
		Element out=e1;
		Element inv=e2.Invert();
		out.Multiply(inv);
		return out;
	}
	friend Element& operator *= (Element& e1,const Element& e2)
	{
		e1.Multiply(e2);
		return e1;
	}
	friend Element& operator /= (Element& e1,const Element& e2)
	{
		Element inv=e2;
		inv.Invert();
		e1.Multiply(inv);
		return e1;
	}
};

template<class Real,class Element>
class VectorSpace
{
public:
	typedef Real R;
	// For this to work, need to define:
	// void Element::Add(const Element& e)
	// void Element::Scale(Real s)
	friend Element  operator -  ( const Element& e )
	{
		Element out=e;
		out.Scale(-1);
		return out;
	}

	friend Element operator + (const Element& e1,const Element& e2)
	{
		Element out=e1;
		out.Add(e2);
		return out;
	};
	friend Element  operator -  (const Element& e1,const Element& e2)
	{
		Element out=e2;
		out.Scale(-1);
		out.Add(e1);
		return out;
	}
	friend Element& operator += (Element& e1,const Element& e2)
	{
		e1.Add(e2);
		return e1;
	}
	friend Element& operator -= (Element& e1,const Element& e2)
	{
		Element neg=e2;
		neg.Scale(-1);
		return e1+=neg;
	}

	friend Element  operator * (const Element& e,const Real& s)
	{
		Element out=e;
		out.Scale(s);
		return out;
	}
	friend Element  operator * (const Real& s,const Element& e)
	{
		return e*s;
	}
	friend Element  operator /  (const Element& e,const Real& s)
	{
		return e*(Real(1)/s);
	}

	friend Element& operator *= (Element& e,const Real& s)
	{
		e.Scale(s);
		return e;
	}
	friend Element& operator /= (Element& e,const Real& s)
	{
		return e*=Real(1)/s;
	}
};

template< class Real , class ... Elements > class VectorSpaceDirectSumElement;

template< class Real , class ... Elements >
class VectorSpaceDirectSum
{
	typedef VectorSpaceDirectSumElement< Real , Elements ... > _Element;
public:
	struct Element : public std::tuple< Elements ... >{};
	typedef Real R;
	// For this to work, each element needs to define:
	// void Element::Add(const Element& e)
	// void Element::Scale(Real s)
	friend _Element  operator -  ( const _Element& e )
	{
		_Element out = e;
		_Scale( out , -1 );
		return out;
	}

	friend _Element operator + (const _Element& e1,const _Element& e2)
	{
		_Element out=e1;
		_Add( out , e2 );
		return out;
	};
	friend _Element  operator -  (const _Element& e1,const _Element& e2)
	{
		_Element out=e2;
		_Scale( out , -1 );
		_Add( out , e1 );
		return out;
	}
	friend _Element& operator += (_Element& e1,const _Element& e2)
	{
		_Add( e1 , e2 );
		return e1;
	}
	friend _Element& operator -= (_Element& e1,const _Element& e2)
	{
		_Element neg=e2;
		_Scale( neg , -1 );
		return e1+=neg;
	}

	friend _Element  operator * (const _Element& e,const Real& s)
	{
		_Element out=e;
		_Scale( out , s );
		return out;
	}
	friend _Element  operator * (const Real& s,const _Element& e)
	{
		return e*s;
	}
	friend _Element  operator /  (const _Element& e,const Real& s)
	{
		return e*(Real(1)/s);
	}

	friend _Element& operator *= (_Element& e,const Real& s)
	{
		_Scale( e , s );
		return e;
	}
	friend _Element& operator /= (_Element& e,const Real& s)
	{
		return e*=Real(1)/s;
	}
protected:
	template< int Idx=0 > static typename std::enable_if< ( Idx< sizeof...(Elements) ) >::type _Add  ( _Element &e1 , const _Element &e2 ){ std::get< Idx >( e1 ) += std::get< Idx >( e2 ) ; _Add< Idx+1 >( e1 , e2 ); }
	template< int Idx=0 > static typename std::enable_if< ( Idx< sizeof...(Elements) ) >::type _Scale( _Element &e , Real s ){ std::get< Idx >( e ) *= s ; _Scale< Idx+1 >( e , s ); }
	template< int Idx=0 > static typename std::enable_if< ( Idx==sizeof...(Elements) ) >::type _Add  ( _Element &e1 , const _Element &e2 ){}
	template< int Idx=0 > static typename std::enable_if< ( Idx==sizeof...(Elements) ) >::type _Scale( _Element &e , Real s ){}
};
template< class Real , class ... Elements >
class VectorSpaceDirectSumElement : public VectorSpaceDirectSum< Real , Elements ... >::Element , VectorSpaceDirectSum< Real , Elements ... >{};

template<class Real,class Element>
class InnerProductSpace : public VectorSpace< Real , Element >
{
public:
	// For this to work, need to define:
	// Real Element::InnerProduct	(const Element& e) const

	static Real SquareNorm		(const Element& e)						{ return e.InnerProduct(e); }
	static Real Dot				(const Element& e1,const Element& e2)	{ return e1.InnerProduct(e2); }
	static Real SquareDistance	(const Element& e1,const Element& e2)	{ return SquareNorm(e1-e2); }
	static Real Length          (const Element& e)                      { return Real( sqrt( e.InnerProduct(e) ) ); }
	Real squareNorm( void ) const { return SquareNorm( *( ( Element* )this ) ); }
};

template<class Real,class Element>
class Algebra : public VectorSpace<Real,Element>
{
public:
	virtual void SetIdentity	(void)				= 0;
	virtual void Multiply		(const Element& e)	= 0;

	static Element Identity(void)
	{
		Element out;
		out.SetIdentity();
		return out;
	}

	friend Element  operator *  (const Element& e1,const Element& e2)
	{
		Element out=e1;
		out.Multiply(e2);
		return out;
	}

	friend Element& operator *= (Element& e1,const Element& e2)
	{
		e1.Multiply(e2);
		return e1;
	}
};
template< class Element >
class Field
{
public:
	typedef Element R;
	// For this to work, need to define:
	// void Element::SetAdditiveIdentity		(void);
	// void Element::SetMultiplicativeIdentity	(void);
	// void Element::Add						(const Element& e);
	// void Element::Multiply					(const Element& e);
	// void Element::Negate						(void);
	// void Element::Invert						(void);

	static Element AdditiveIdentity(void)
	{
		Element out;
		out.SetAdditiveIdentity();
		return out;
	}
	static Element MultiplicativeIdentity(void)
	{
		Element out;
		out.SetMultiplicativeIdentity();
		return out;
	}
	Element additiveInverse(void) const
	{
		Element out=*(Element*)this;
		out.Negate();
		return out;
	}
	Element multiplicativeInverse(void) const
	{
		Element out=*(Element*)this;
		out.Invert();
		return out;
	}
	friend Element operator + ( const Element& e1 , const Element& e2 )
	{
		Element out=e1;
		out.Add(e2);
		return out;
	}
	friend Element operator * ( const Element& e1 , const Element& e2 )
	{
		Element out=e1;
		out.Multiply(e2);
		return out;
	}
	friend Element operator - (const Element& e1,const Element& e2)
	{
		Element out=e1;
		Element inv=e2.Negate();
		out.Add(inv);
		return out;
	}
	friend Element operator / (const Element& e1,const Element& e2)
	{
		Element out=e1;
		Element inv=e2.Invert();
		out.Multiply(inv);
		return out;
	}
	friend Element& operator += (Element& e1,const Element& e2)
	{
		e1.Add(e2);
		return e1;
	}
	friend Element& operator *= (Element& e1,const Element& e2)
	{
		e1.Multiply(e2);
		return e1;
	}
	friend Element& operator -= (Element& e1,const Element& e2)
	{
		Element inv=e2;
		inv.Negate();
		e1.Add(inv);
		return e1;
	}
	friend Element& operator /= (Element& e1,const Element& e2)
	{
		Element inv=e2;
		inv.Invert();
		e1.Multiply(inv);
		return e1;
	}
};
#endif // ALGEBRA_INCLUDED
