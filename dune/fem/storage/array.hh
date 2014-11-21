#ifndef DUNE_FEM_ARRAY_HH
#define DUNE_FEM_ARRAY_HH

#include <cassert>

#include <dune/fem/misc/bartonnackmaninterface.hh>
#include <dune/fem/misc/metaprogramming.hh>

#include <dune/fem/storage/arrayallocator.hh>

namespace Dune
{
  namespace Fem
  {

    /** \class ArrayInterface
     *  \ingroup VectorClasses
     *  \brief abstract array interface
     */
    template< class AT >
    class ArrayInterface
    : public BartonNackmanInterface< ArrayInterface< AT >, typename AT::ArrayType >
    {
      typedef ArrayInterface< AT > ThisType;
      typedef BartonNackmanInterface< ThisType, typename AT::ArrayType > BaseType;

    public:
      //! type of the traits
      typedef AT Traits;

      //! type of the implementation (Barton-Nackman)
      typedef typename Traits::ArrayType ArrayType;

      //! type of this interface
      typedef ThisType ArrayInterfaceType;

      //! type of the array elements
      typedef typename Traits::ElementType ElementType;

      //! make consistent with std::vector
      typedef ElementType value_type ;

      //! type of constant iterator
      typedef typename Traits::ConstIteratorType ConstIteratorType;

      //! type of (non-constant) iterator
      typedef typename Traits::IteratorType IteratorType;

    protected:
      using BaseType::asImp;

    public:
      /** \brief access an array element
       *
       *  \param[in]  index  index of the array element to access
       *
       *  \returns a const reference to the array element
       */
      const ElementType &operator[] ( unsigned int index ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
        return asImp()[ index ];
      }

      /** \brief access an array element
       *
       *  \param[in]  index  index of the array element to access
       *
       *  \returns a reference to the array element
       */
      ElementType &operator[] ( unsigned int index )
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp()[ index ] );
        return asImp()[ index ];
      }

      /** \brief fill the array with copies of an element
       *
       *  \param[in]  element  element wich shall be copied into every array
       *                       entry
       */
      void assign ( const ElementType &element )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().assign( element ) );
      }

      /** \brief copy another array to this one
       *
       *  Copies the data from another array to this one. Both arrays must be of
       *  the same size.
       *
       *  \param[in]  other  array to copy
       */
      template< class T >
      void assign( const ArrayInterface< T > &other )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().assign( other ) );
      }

      /** \brief obtain begin iterator
       *
       *  \returns an iterator pointing to the first array element
       */
      ConstIteratorType begin () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
        return asImp().begin();
      }

      /** \brief obtain begin iterator
       *
       *  \returns an iterator pointing to the first array element
       */
      IteratorType begin ()
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
        return asImp().begin();
      }

      /** \brief obtain end iterator
       *
       *  \returns an iterator pointing behind the last array element
       */
      ConstIteratorType end () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
        return asImp().end();
      }

      /** \brief obtain end iterator
       *
       *  \returns an iterator pointing behind the last array element
       */
      IteratorType end ()
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
        return asImp().end();
      }

      /** obtain the size of the array
       *
       *  \returns the size of the array
       */
      unsigned int size () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().size() );
        return asImp().size();
      }
    };


    template< class Array >
    struct SupportsArrayInterface
    {
      typedef ArrayInterface< typename Array::Traits > ArrayInterfaceType;
      static const bool v = Conversion< Array, ArrayInterfaceType >::exists;
    };


    template< class Element, class Array >
    class ArrayDefaultIterator
    {
      typedef ArrayDefaultIterator< Element, Array > ThisType;

    public:
      typedef Element ElementType;

      typedef Array ArrayType;

    protected:
      ArrayType &array_;
      unsigned int index_;

    public:
      ArrayDefaultIterator ( ArrayType &array, unsigned int index )
      : array_( array ),
        index_( index )
      {
        assert( index <= array.size() );
      }

      ArrayDefaultIterator( const ThisType &other )
      : array_( other.array_ ),
        index_( other.index_ )
      {}

      ThisType &operator= ( const ThisType &other )
      {
        assert( &(other.array_) == &array_ );
        index_ = other.index_;
      }

      ElementType &operator* ()
      {
        assert( index_ < array_.size() );
        return array_[ index_ ];
      }

      ThisType &operator++ ()
      {
        assert( index_ < array_.size() );
        ++index_;
        return *this;
      }

      bool operator== ( const ThisType &other ) const
      {
        assert( &(other.array_) == &array_ );
        return index_ == other.index_;
      }

      bool operator!= ( const ThisType &other ) const
      {
        return !(*this == other);
      }

      unsigned int index () const
      {
        return index_;
      }
    };


    template< class ElementImp, class ArrayImp >
    struct ArrayDefaultTraits
    {
      typedef ElementImp ElementType;

      typedef ArrayImp ArrayType;

      typedef ArrayDefaultIterator< ElementType, ArrayType > IteratorType;

      typedef ArrayDefaultIterator< const ElementType, const ArrayType > ConstIteratorType;
    };


    /** \class ArrayDefault
     *  \ingroup VectorClasses
     *  \brief default implementation of the ArrayInterface
     */
    template< class ElementImp, class ArrayImp >
    class ArrayDefault
    : public ArrayInterface< ArrayDefaultTraits< ElementImp, ArrayImp > >
    {
    public:
      //! type of the array elements
      typedef ElementImp ElementType;

      //! type of the implementation (Barton-Nackman)
      typedef ArrayImp ArrayType;

      //! type of the traits
      typedef ArrayDefaultTraits< ElementType, ArrayType > Traits;

    private:
      typedef ArrayDefault< ElementType, ArrayType > ThisType;
      typedef ArrayInterface< Traits > BaseType;

    public:
      using BaseType :: size;

    protected:
      using BaseType :: asImp;

    public:
      //! type of constant iterator
      typedef typename Traits :: ConstIteratorType ConstIteratorType;

      //! type of (non-constant) iterator
      typedef typename Traits :: IteratorType IteratorType;

    public:
      /** \copydoc Dune::Fem::ArrayInterface::assign(const ElementType &element) */
      inline void assign ( const ElementType &element )
      {
        ArrayType &imp = asImp();
        const unsigned int size = imp.size();
        for( unsigned int i = 0; i < size; ++i )
          imp[ i ] = element;
      }

      /** \copydoc Dune::Fem::ArrayInterface::assign(const ArrayInterface<T> &other) */
      template< class T >
      inline void assign( const ArrayInterface< T > &other )
      {
        ArrayType &imp = asImp();
        const unsigned int size = imp.size();
        assert( size == other.size() );
        for( unsigned int i = 0; i < size; ++i )
          imp[ i ] = other[ i ];
      }

      /** \copydoc Dune::Fem::ArrayInterface::begin() const */
      inline ConstIteratorType begin () const
      {
        return ConstIteratorType( asImp(), 0 );
      }

      /** \copydoc Dune::Fem::ArrayInterface::begin() */
      inline IteratorType begin ()
      {
        return IteratorType( asImp(), 0 );
      }

      /** \copydoc Dune::Fem::ArrayInterface::end() const */
      inline ConstIteratorType end () const
      {
        return ConstIteratorType( asImp(), size() );
      }

      /** \copydoc Dune::Fem::ArrayInterface::end() */
      inline IteratorType end ()
      {
        return IteratorType( asImp(), size() );
      }
    };


    /** \class ArrayWrapper
     *  \ingroup VectorClasses
     *  \brief implementation of the ArrayInterface wrapping a pointer to an
     *         array of elements
     */
    template< class ElementImp >
    class ArrayWrapper
    : public ArrayDefault< ElementImp, ArrayWrapper< ElementImp > >
    {
    public:
      //! type of the array elements
      typedef ElementImp ElementType;

    private:
      typedef ArrayWrapper< ElementType > ThisType;
      typedef ArrayDefault< ElementType, ThisType > BaseType;

    private:
      const unsigned int size_;
      ElementType *elements_;

    public:
      /** \brief create an ArrayWrapper from a size and a pointer */
      inline ArrayWrapper ( unsigned int size, ElementType *elements )
      : size_( size ),
        elements_( elements )
      {
        assert( elements_ != NULL );
      }

      /** \copydoc Dune::Fem::ArrayInterface::operator[](unsigned int index) const */
      inline const ElementType &operator[] ( unsigned int index ) const
      {
        assert( index < size_ );
        return elements_[ index ];
      }

      /** \copydoc Dune::Fem::ArrayInterface::operator[](unsigned int index) */
      inline ElementType &operator[] ( unsigned int index )
      {
        assert( index < size_ );
        return elements_[ index ];
      }

      /** \copydoc Dune::Fem::ArrayInterface::size */
      inline unsigned int size () const
      {
        return size_;
      }
    };


    /** \class FixedSizeArray
     *  \ingroup VectorClasses
     *  \brief standard array with fixed size
     *
     *  This array's size is a compile time constant. Hence, the memory needed
     *  to store the array elements can be provided within the object. There
     *  is no need to allocate it dynamically.
     *
     *  Basically, FixedSizeArray is just a standard C++ array. It's strength
     *  lies in the following two properties:
     *  - It can be passed by value or reference. A standard C++ array is always
     *    passed as a pointer.
     *  - It satisfies the array interface. It can be used whenever this
     *    interface is explicitly required.
     *  .
     *
     *  \param  Element  type of the array elements (must be default
     *                   constructable)
     *  \param  Size     number of elements in the array
     */
    template< class Element, unsigned int Size >
    class FixedSizeArray
    : public ArrayDefault< Element, FixedSizeArray< Element, Size > >
    {
      typedef FixedSizeArray< Element, Size > ThisType;
      typedef ArrayDefault< Element, ThisType > BaseType;

    public:
      //! type of the array elements
      typedef Element ElementType;

      //! compile time constant size of the array
      static const unsigned int fixedSize = Size;

    protected:
      ElementType elements_[ fixedSize ];

    public:
      /** \brief default constructor
       *
       *  The array elements are not initialized with this constructor
       */
      inline FixedSizeArray ()
      {}

      /** \brief initializing constructor
       *
       *  Initializes the entire array with a default value
       *
       *  \param[in]  element  default value
       */
      inline explicit FixedSizeArray ( const ElementType &element )
      {
        assign( element );
      }

      /** \brief copy constructor
       *
       *  \param[in]  other  array to copy
       */
      inline FixedSizeArray ( const ThisType &other )
      {
        assign( other );
      }

      /** \copydoc Dune::Fem::ArrayInterface::operator[](unsigned int index) const */
      inline const ElementType &operator[] ( unsigned int index ) const
      {
        assert( index < fixedSize );
        return elements_[ index ];
      }

      /** \copydoc Dune::Fem::ArrayInterface::operator[](unsigned int index) */
      inline ElementType &operator[] ( unsigned int index )
      {
        assert( index < fixedSize );
        return elements_[ index ];
      }

      /** \copydoc Dune::Fem::ArrayInterface::size() const */
      inline unsigned int size () const
      {
        return fixedSize;
      }
    };


    template< class Element,
              template< class > class ArrayAllocator = DefaultArrayAllocator >
    class DynamicArray
    : public ArrayDefault< Element, DynamicArray< Element, ArrayAllocator > >
    {
      typedef DynamicArray< Element, ArrayAllocator > ThisType;
      typedef ArrayDefault< Element, ThisType > BaseType;

    public:
      typedef Element ElementType;
      typedef ElementType FieldType ;
      typedef ElementType value_type ;

    protected:
      typedef ArrayAllocator< ElementType > ArrayAllocatorType;

      typedef typename ArrayAllocatorType :: ElementPtrType ElementPtrType;

    protected:
      ArrayAllocatorType allocator_;

      unsigned int size_;
      ElementPtrType elements_;

    public:
      using BaseType :: assign;

    public:
      inline explicit DynamicArray ( unsigned int size = 0 )
      : allocator_()
      {
        size_ = size;
        allocator_.allocate( size_, elements_ );
      }

      inline explicit DynamicArray ( const ArrayAllocatorType &arrayAllocator,
                                     unsigned int size = 0 )
      : allocator_( arrayAllocator )
      {
        size_ = size;
        allocator_.allocate( size_, elements_ );
      }

      inline DynamicArray ( unsigned int size,
                            const ElementType &element )
      : allocator_()
      {
        size_ = size;
        allocator_.allocate( size_, elements_ );
        assign( element );
      }

      inline DynamicArray ( const ArrayAllocatorType &arrayAllocator,
                            unsigned int size,
                            const ElementType defaultElement )
      : allocator_( arrayAllocator )
      {
        size_ = size;
        allocator_.allocate( size_, elements_ );
        assign( defaultElement );
      }

      inline DynamicArray ( const ThisType &other )
      : allocator_( other.allocator_ )
      {
        size_ = other.size_;
        allocator_.allocate( size_, elements_ );
        assign( other );
      }

      inline ~DynamicArray ()
      {
        allocator_.free( elements_ );
      }

      inline const ElementType& operator[] ( unsigned int index ) const
      {
        assert( index < size_ );
        return elements_[ index ];
      }

      inline ElementType& operator[] ( unsigned int index )
      {
        assert( index < size_ );
        return elements_[ index ];
      }

      inline void append ( const ElementType &element )
      {
        const unsigned int oldSize = size_;
        resize( oldSize + 1 );
        elements_[ oldSize ] = element;
      }

      template< class T >
      inline void append ( const ArrayInterface< T > &array )
      {
        const unsigned int arraySize = array.size();

        const unsigned int oldSize = size_;
        resize( oldSize + arraySize );

        for( unsigned int i = 0; i < arraySize; ++i )
          elements_[ oldSize + i ] = array[ i ];
      }

      /** \copydoc Dune::Fem::ArrayInterface::assign(const ArrayInterface<T> &other) */
      template< class T >
      inline void assign( const ArrayInterface< T > &other )
      {
        resize( other.size() );
        for( unsigned int i = 0; i < size_; ++i )
          elements_[ i ] = other[ i ];
      }

      inline ElementType *leakPointer ()
      {
        return (ElementType *)elements_;
      }

      inline const ElementType *leakPointer () const
      {
        return (const ElementType *)elements_;
      }

      inline void reserve ( unsigned int newSize )
      {
        allocator_.reserve( newSize, elements_ );
      }

      inline void resize ( unsigned int newSize )
      {
        const unsigned int oldSize = size_;
        if( newSize == oldSize )
          return;

        allocator_.reallocate( oldSize, newSize, elements_ );
        size_ = newSize;
      }

      inline void resize ( unsigned int newSize,
                           const ElementType &defaultElement )
      {
        const unsigned int oldSize = size_;
        if( newSize == oldSize )
          return;

        allocator_.reallocate( oldSize, newSize, elements_ );
        size_ = newSize;
        for( unsigned int i = oldSize; i < newSize; ++i )
          elements_[ i ] = defaultElement;
      }

      inline unsigned int size () const
      {
        return size_;
      }
    };

    // Capabilities
    // ------------

    namespace Capabilities
    {

      template< class Array >
      struct HasLeakPointer
      : public Fem::MetaBool< false >
      {};

      template< class Element, template< class > class ArrayAllocator >
      struct HasLeakPointer< Fem :: DynamicArray< Element, ArrayAllocator > >
      : public Fem::MetaBool< true >
      {};

    }

  } // namespace Fem

} // namespace Dune

#include "array_inline.hh"

#endif // #ifndef DUNE_FEM_ARRAY_HH
