#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_MUTABLE_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_MUTABLE_HH

//-s system includes 
#include <cassert>

//- Dune includes 
#include <dune/fem/storage/objectstack.hh>
#include <dune/fem/function/localfunction/localfunction.hh>

namespace Dune
{

  namespace Fem 
  {

    template< class >
    struct DiscreteFunctionTraits;

    template< class >
    class ConstLocalFunction;

    //**************************************************************************
    //
    //  --MutableLocalFunction 
    //
    //**************************************************************************
    //! Manages the getting and deleting of local function pointers and 
    //! acts like a local functions 
    template < class DiscreteFunction >
    class MutableLocalFunction
    : public LocalFunction< typename DiscreteFunctionTraits< DiscreteFunction > :: DiscreteFunctionSpaceType ::  BasisFunctionSetType,
      typename DiscreteFunctionTraits< DiscreteFunction > :: LocalDofVectorType >
    {
      typedef MutableLocalFunction< DiscreteFunction > ThisType;
      typedef LocalFunction< typename DiscreteFunctionTraits< DiscreteFunction > :: DiscreteFunctionSpaceType ::  BasisFunctionSetType,
      typename DiscreteFunctionTraits< DiscreteFunction > :: LocalDofVectorType > BaseType;

    public:
      //! type of DiscreteFunction 
      typedef DiscreteFunction DiscreteFunctionType;

      //! type of the entity, the local function lives on is given by the space
      typedef typename BaseType::EntityType EntityType;

      //! type of local Dof vector object
      typedef typename BaseType::LocalDofVectorType LocalDofVectorType;

      //! type of BasisFunctionSet
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;


      //! cast from ConstLocalFunction
      MutableLocalFunction( const ConstLocalFunction< DiscreteFunctionType > &constLocalFunction ) DUNE_DEPRECATED
      : BaseType( constLocalFunction.basisFunctionSet(), LocalDofVectorType( constLocalFunction.discreteFunction_.localDofVectorAllocator() ) ),
        discreteFunction_( const_cast< DiscreteFunctionType& >( constLocalFunction.discreteFunction_ ) )
      {
        discreteFunction_.getLocalDofs( constLocalFunction.entity(), localDofVector() );
      }


      //! Constructor creating empty local function from given discrete function 
      explicit MutableLocalFunction ( DiscreteFunctionType &discreteFunction ) 
      : BaseType( LocalDofVectorType( discreteFunction.localDofVectorAllocator() ) ),
        discreteFunction_( discreteFunction )
      {
      }

      //! Constructor creating local function from given discrete function and entity, not empty
      explicit MutableLocalFunction ( DiscreteFunctionType &discreteFunction, const EntityType &entity ) 
      : BaseType( discreteFunction.space().basisFunctionSet( entity ), LocalDofVectorType( discreteFunction.localDofVectorAllocator() ) ), 
        discreteFunction_( discreteFunction )
      {
        discreteFunction_.getLocalDofs( entity, localDofVector() );
      }

      MutableLocalFunction ( const ThisType & ) = default;
      MutableLocalFunction ( ThisType && ) = default;

      void init ( const EntityType &entity )
      {
        BaseType::init( discreteFunction_.space().basisFunctionSet( entity ) );
        discreteFunction_.getLocalDofs( entity, localDofVector() );
      }

      using BaseType::localDofVector;

    private:
      DiscreteFunctionType &discreteFunction_;
      // prohibit assignment 
      ThisType &operator= ( const ThisType & );
    }; // end MutableLocalFunction  

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_MUTABLE_HH