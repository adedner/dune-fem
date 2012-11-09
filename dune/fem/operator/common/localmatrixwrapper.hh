#ifndef DUNE_FEM_LOCALMATRIXWRAPPER_HH
#define DUNE_FEM_LOCALMATRIXWRAPPER_HH

#include <dune/fem/operator/common/localmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    template< class LocalMatrixStackImp >
    class LocalMatrixWrapper;


    template< class LocalMatrixStackImp >
    struct LocalMatrixWrapperTraits
    {
      typedef LocalMatrixStackImp LocalMatrixStackType;

      typedef typename LocalMatrixStackImp :: ObjectType WrappedLocalMatrixType;

      typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;

      typedef typename WrappedLocalMatrixType :: RangeFieldType RangeFieldType;
      
      typedef typename WrappedLocalMatrixType :: DomainSpaceType DomainSpaceType;
      typedef typename WrappedLocalMatrixType :: RangeSpaceType RangeSpaceType;

      typedef typename WrappedLocalMatrixType :: LittleBlockType LittleBlockType;
    };


    template< class LocalMatrixStackImp >
    class LocalMatrixWrapper
    : public LocalMatrixInterface< LocalMatrixWrapperTraits< LocalMatrixStackImp > >
    {
    public:
      //! type of the local matrix stack
      typedef LocalMatrixStackImp LocalMatrixStackType;

      //! type of the traits
      typedef LocalMatrixWrapperTraits< LocalMatrixStackType > Traits;

    private:
      typedef LocalMatrixWrapper< LocalMatrixStackType > ThisType;
      typedef LocalMatrixInterface< Traits > BaseType;

    public:
      //! type of the wrapped local matrix
      typedef typename Traits :: WrappedLocalMatrixType WrappedLocalMatrixType;

      typedef typename Traits :: RangeFieldType RangeFieldType;

      typedef typename BaseType :: DomainSpaceType DomainSpaceType;
      typedef typename BaseType :: RangeSpaceType RangeSpaceType;

      typedef typename BaseType :: DomainBasisFunctionSetType
        DomainBasisFunctionSetType;
      typedef typename BaseType :: RangeBasisFunctionSetType
        RangeBasisFunctionSetType;

    private:
      typedef typename LocalMatrixStackType :: PointerType
        WrappedLocalMatrixPtrType;

    private:
      // ObjectPointer to the actual local matrix
      // (the pointer is required to keep the reference alive)
      WrappedLocalMatrixPtrType localMatrixPtr_;

      // reference to the actual local matrix
      WrappedLocalMatrixType &localMatrix_;

    public:
      //! constructor creating an uninitialized local matrix
      inline explicit LocalMatrixWrapper ( LocalMatrixStackType &stack )
      : localMatrixPtr_( stack.getObject() ),
        localMatrix_( *localMatrixPtr_ )
      {
      }
      
      //! constructor initializing the wrapped local matrix
      template< class DomainEntityType, class RangeEntityType >
      inline LocalMatrixWrapper( LocalMatrixStackType &stack,
                                 const DomainEntityType &domainEntity,
                                 const RangeEntityType &rangeEntity )
      : localMatrixPtr_( stack.getObject() ),
        localMatrix_( *localMatrixPtr_ )
      {
        // initialize the wrapped local matrix with the entities
        localMatrix().init( domainEntity, rangeEntity );
      }

      /** \brief copy constructor
       *
       *  \param[in]  other  LocalMatrixWrapper to copy
       */
      inline LocalMatrixWrapper ( const ThisType &other )
      : localMatrixPtr_( other.localMatrixPtr_ ),
        localMatrix_( *localMatrixPtr_ )
      {
      }

    private:
      // prohibit assignment
      inline ThisType &operator= ( const ThisType & );

    public:
      /** \copydoc Dune::Fem::LocalMatrixInterface::init */
      template< class DomainEntityType, class RangeEntityType >
      inline void init ( const DomainEntityType &domainEntity,
                         const RangeEntityType &rangeEntity )
      {
        localMatrix().init( domainEntity, rangeEntity );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::add */
      inline void add ( const int localRow,
                        const int localCol,
                        const RangeFieldType &value )
      {
        localMatrix().add( localRow, localCol, value );
      }
      
      /** \copydoc Dune::Fem::LocalMatrixInterface::set */
      inline void set ( const int localRow,
                        const int localCol,
                        const RangeFieldType &value )
      {
        localMatrix().set( localRow, localCol, value );
      }
      
      /** \copydoc Dune::Fem::LocalMatrixInterface::unitRow */
      inline void unitRow ( const int localRow ) DUNE_DEPRECATED
      {
        localMatrix().unitRow( localRow ); 
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::clearRow */
      inline void clearRow ( const int localRow )
      {
        localMatrix().clearRow( localRow );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::clearRow */
      inline void clearCol ( const int localCol )
      {
        localMatrix().clearCol( localCol );
      }
      
      /** \copydoc Dune::Fem::LocalMatrixInterface::get */
      inline const RangeFieldType get ( const int localRow,
                                        const int localCol ) const
      {
        return localMatrix().get( localRow, localCol );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::scale */
      inline void scale ( const RangeFieldType& scalar ) 
      {
        return localMatrix().scale( scalar );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::clear */
      inline void clear ()
      {
        return localMatrix().clear();
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::resort */
      inline void resort ()
      {
        return localMatrix().resort();
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::rows */
      inline int rows () const
      {
        return localMatrix().rows();
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::columns */
      inline int columns () const
      {
        return localMatrix().columns();
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::multiplyAdd */
      template <class DomainLocalFunctionImp, 
                class RangeLocalFunctionImp>
      inline void multiplyAdd(const DomainLocalFunctionImp& dLf,
                              RangeLocalFunctionImp& rLf)
      {
        localMatrix().multiplyAdd( dLf, rLf); 
      }
      
      /** \copydoc Dune::Fem::LocalMatrixInterface::domainSpace */
      inline const DomainSpaceType &domainSpace () const
      {
        return localMatrix().domainSpace();
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::rangeSpace */
      inline const RangeSpaceType &rangeSpace () const
      {
        return localMatrix().rangeSpace();
      }
      
      /** \copydoc Dune::Fem::LocalMatrixInterface::domainBasisFunctionSet */
      inline const DomainBasisFunctionSetType &domainBasisFunctionSet () const
      {
        return localMatrix().domainBasisFunctionSet();
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::rangeBasisFunctionSet */
      inline const RangeBasisFunctionSetType &rangeBasisFunctionSet () const
      {
        return localMatrix().rangeBasisFunctionSet();
      }

    protected:
      inline const WrappedLocalMatrixType &localMatrix () const
      {
        return localMatrix_;
      }

      inline WrappedLocalMatrixType &localMatrix ()
      {
        return localMatrix_;
      }
    };

  } // namespace Fem

#if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: LocalMatrixWrapper ;

#endif // DUNE_FEM_COMPATIBILITY
  
} // namespace Dune

#endif // #ifndef DUNE_FEM_LOCALMATRIXWRAPPER_HH
