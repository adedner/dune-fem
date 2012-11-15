#ifndef DUNE_FEM_FUNCTION_HH
#define DUNE_FEM_FUNCTION_HH

// dune-common includes 
#include <dune/common/fvector.hh>

// dune-fem includes 
#include <dune/fem/misc/bartonnackmaninterface.hh>
#include <dune/fem/operator/common/mapping.hh>
#include <dune/fem/version.hh>


namespace Dune
{

  namespace Fem
  {

    /** @addtogroup Functions
        Functions are mapping from one finite dimensional
        vector space into another, e.g.,
        \f$K^n\f$ into \f$L^m\f$.
        They are element of a 
        \ref FunctionSpaceInterface "function space".

        \remark The interface is given by Function.
        @{
    **/
  

    /*! \brief 
        Abstract class representing a function


        Template parameters are:
        -  FunctionSpaceImp      type of the function space where the function 
                                 belongs to.
        -  FunctionImp           type of the implemented function (Barton-Nackman)
       
        @interfaceclass
    **/
    template< class FunctionSpaceImp, class FunctionImp >
    class Function
    : public BartonNackmanInterface< Function< FunctionSpaceImp, FunctionImp >,
                                     FunctionImp >,
      public Mapping < typename FunctionSpaceImp :: DomainFieldType,
                       typename FunctionSpaceImp :: RangeFieldType,
                       typename FunctionSpaceImp :: DomainType, 
                       typename FunctionSpaceImp :: RangeType > 
    {
      typedef Function< FunctionSpaceImp, FunctionImp > ThisType;
      typedef BartonNackmanInterface< ThisType, FunctionImp > BaseType;

    public:
      //! type of function space this function belongs to 
      typedef FunctionSpaceImp FunctionSpaceType;

      //! type of the implementation (Barton-Nackman)
      typedef FunctionImp FunctionType;

      //! field type of domain
      typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
      //! field type of range
      typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
      //! domain type
      typedef typename FunctionSpaceType :: DomainType DomainType;
      //! range type
      typedef typename FunctionSpaceType :: RangeType RangeType;
      //! jacobian type
      typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;
      //! hessian type
      typedef typename FunctionSpaceType :: HessianRangeType HessianRangeType;

      //! type of mapping base class
      typedef Mapping< DomainFieldType, RangeFieldType, DomainType, RangeType >
        MappingType;
   
    protected:
      using BaseType::asImp;

    protected:
      /** \brief default constructor */
      Function ()
      {}

      /** \brief constructor storing the function space
       *
       *  \param[in]  fSpace  function space for this function
       */
      explicit Function ( const FunctionSpaceType &fSpace ) DUNE_DEPRECATED
      {}

    protected:
      // Disallow copying of function, but allow copying of derived classes 
      Function ( const ThisType &other ) {} 
    private:  
      // Disallow assignment
      ThisType &operator= ( const ThisType & );
      
    public:
      //! destructor 
      virtual ~Function ()
      {}

      /** \brief application operator call evaluate 
          \param[in] arg argument 
          \param[out] dest destination, i.e. f(arg) 
      */
      virtual void operator()(const DomainType & arg, RangeType & dest) const 
      {
        evaluate(arg,dest);
      }

      /** \brief evaluate the function
       *
       *  \param[in]  x    evaluation point 
       *  \param[out] ret  value of the function in x
       */
      void evaluate ( const DomainType &x, RangeType &ret ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().evaluate( x, ret ) );
      }

      /** \brief evaluate the Jacobian of the function
       *
       *  \param[in]  x    evaluation point
       *  \param[out] ret  value of the Jacobian in x
       */
      void jacobian ( const DomainType &x, JacobianRangeType &ret ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().jacobian( x, ret ) );
      }

      /** \brief evaluate the hessian of the function
       *
       *  \param[in]  x    evaluation point
       *  \param[out] ret  value of the hessian in x
       */
      void hessian ( const DomainType &x, HessianRangeType &ret ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().hessian( x, ret ) );
      }

      /** \brief evaluate a derivative of the function
       * 
       *  \param[in]  diffVariable  vector describing the partial derivative to
       *                            evaluate
       *  \param[in]  x             evaluation point
       *  \param[out] ret           value of the derivative in x 
       */    
      template< int diffOrder >
      DUNE_VERSION_DEPRECATED(1,4,remove)
      void evaluate ( const FieldVector< int, diffOrder > &diffVariable,
                      const DomainType &x,
                      RangeType &ret ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().evaluate( diffVariable, x, ret ) );
      }

      /** \brief Obtain the related function space
       * 
       *  \returns a reference to the function space 
       */ 
      const FunctionSpaceType &space () const DUNE_DEPRECATED
      {
        static FunctionSpaceType functionSpace;
        return functionSpace;
      }

    private:
      //! Helper function for Mapping
      //! With this function, a combined mapping can choose the right application
      //! operator (i.e. the one from Mapping itself, or from Function/Operator)
      //! \note: Do not override this definition
      virtual void apply (const DomainType& arg, RangeType& dest) const 
      {
        operator()(arg, dest);
      }
    };

    ///@} 

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_HH
