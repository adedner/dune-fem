#ifndef DUNE_FEM_DIFFUSIONMODEL_HH
#define DUNE_FEM_DIFFUSIONMODEL_HH

#include <dune/fem/misc/bartonnackmaninterface.hh>
#include <dune/fem/operator/common/mapping.hh>

#ifndef DUNE_FEM_COMPATIBILITY
#define DUNE_FEM_COMPATIBILITY 1
#endif

namespace Dune
{

  /*! \ingroup FEM
   *  \class DiffusionModelInterface
   *  \brief Interface for a mathematical model of a DiffusionOperator
   *
   *  The DiffusionModelInterface specifies the way mathematical data can
   *  be incorporated into a DiffusionOperator, i.e. a linear operator
   *  performing
   *  \f[
   *  -\nabla \cdot \bigl( a( x ) \nabla u \bigr).
   *  \f]
   *  Here the matematical data is exactly the function \f$a\f$, which we call
   *  diffusive flux.
   *
   *  \param  FunctionSpaceImp   function space to work in
   *  \param  DiffusionModelImp  actual implementation (Barton-Nackman)
   */
  template< class FunctionSpaceImp, class DiffusionModelImp >
  class DiffusionModelInterface
  : public BartonNackmanInterface
    < DiffusionModelInterface< FunctionSpaceImp, DiffusionModelImp >,
      DiffusionModelImp >
  {
  public:
    //! type of the function space we are using
    typedef FunctionSpaceImp FunctionSpaceType;

    //! type of the implementation (Barton-Nackman)
    typedef DiffusionModelImp DiffusionModelType;

  private:
    typedef DiffusionModelInterface< FunctionSpaceType, DiffusionModelType >
      ThisType;
    typedef BartonNackmanInterface< ThisType, DiffusionModelType > BaseType;

  public:
    //! type of the interface
    typedef ThisType DiffusionModelInterfaceType;

    //! type of points within the domain
    typedef typename FunctionSpaceType :: DomainType DomainType;
    //! type of points within the range
    typedef typename FunctionSpaceType :: RangeType RangeType;
    //! type of the Jacobian (evaluated in some point)
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    //! field type of the domain
    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    //! field type of the range
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  protected:
    using BaseType :: asImp;

  public:
    /** \brief calculate the diffusive flux \f$a( x ) \nabla u\f$ in a point
     *
     *  \param[in]  diffVariable  vector describin the partial derivative to
     *                            evaluate
     *  \param[in]  entity        entity to evaluate the flux on
     *  \param[in]  x             evaluation point (in local coordinates)
     *  \param[in]  gradient      \f$\nabla u\f$ in the evaluation point
     *  \param[out] flux          variable to receive the evaluated flux
     */
    template< int diffOrder, class EntityType, class PointType >
    inline void diffusiveFlux ( const FieldVector< deriType, diffOrder > &diffVariable,
                                const EntityType &entity,
                                const PointType &x,
                                const JacobianRangeType &gradient,
                                JacobianRangeType &flux ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().diffusiveFlux( diffVariable, entity, x, gradient, flux ) );
    }

    /** \brief calculate the diffusive flux \f$a( x ) \nabla u\f$ in a point
     *
     *  \param[in]  entity      entity to evaluate the flux on
     *  \param[in]  x           evaluation point (in local coordinates)
     *  \param[in]  gradient    \f$\nabla u\f$ in the evaluation point
     *  \param[out] flux        variable to receive the evaluated flux
     */
    template< class EntityType, class PointType >
    inline void diffusiveFlux ( const EntityType &entity,
                                const PointType &x,
                                const JacobianRangeType &gradient,
                                JacobianRangeType &flux ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().diffusiveFlux( entity, x, gradient, flux ) );
    }

#if DUNE_FEM_COMPATIBILITY
    /** \brief calculate the diffusive flux \f$a( x ) \nabla u\f$ in a quadrature
     *         point
     *
     *  \param[in]  entity      entity to evaluate the flux on
     *  \param[in]  quadrature  quadrature to use
     *  \param[in]  quadPoint   number of the point within the quadrature
     *  \param[in]  gradient    \f$\nabla u\f$ in the evaluation point
     *  \param[out] flux        variable to receive the evaluated flux
     */
    template< class EntityType, class QuadratureType >
    inline void diffusiveFlux ( const EntityType &entity,
                                const QuadratureType &quadrature,
                                const int quadPoint,
                                const JacobianRangeType &gradient,
                                JacobianRangeType &flux ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().diffusiveFlux( entity, quadrature, quadPoint, gradient, flux ) );
    }
#endif
  };


  
  template< class FunctionSpaceImp, class DiffusionModelImp >
  class DiffusionModelDefault
  : public DiffusionModelInterface< FunctionSpaceImp, DiffusionModelImp >
  {
  public:
    //! type of the function space we are using
    typedef FunctionSpaceImp FunctionSpaceType;

    //! type of the implementation (Barton-Nackman)
    typedef DiffusionModelImp DiffusionModelType;

  private:
    typedef DiffusionModelDefault< FunctionSpaceType, DiffusionModelType >
      ThisType;
    typedef DiffusionModelInterface< FunctionSpaceType, DiffusionModelType >
      BaseType;

  public:
    //! type of points within the domain
    typedef typename FunctionSpaceType :: DomainType DomainType;
    //! type of points within the range
    typedef typename FunctionSpaceType :: RangeType RangeType;
    //! type of the Jacobian (evaluated in some point)
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    //! field type of the domain
    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    //! field type of the range
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    using BaseType :: diffusiveFlux;

  protected:
    using BaseType :: asImp;

  public:
    /** \copydoc Dune::DiffusionModelInterface::diffusiveFlux(const EntityType &entity,const PointType &x,const JacobianRangeType &gradient,JacobianRangeType &flux) const
     *
     *  The default implementation calls
     *  \code
     *  FieldVector< deriType, 0 > diffVar;
     *  diffusiveFlux( diffVar, entity, x, gradient, flux );
     *  \endcode
     */
    template< class EntityType, class PointType >
    inline void diffusiveFlux ( const EntityType &entity,
                                const PointType &x,
                                const JacobianRangeType &gradient,
                                JacobianRangeType &flux ) const
    {
      FieldVector< deriType, 0 > diffVar;
      asImp().diffusiveFlux( diffVar, entity, x, gradient, flux );
    }

#if DUNE_FEM_COMPATIBILITY
    /** \copydoc Dune::DiffusionModelInterface::diffusiveFlux(const EntityType &entity,const QuadratureType &quadrature,const int quadPoint,const JacobianRangeType &gradient,JacobianRangeType &flux) const
     *
     *  The default implementation calls
     *  \code
     *  diffusiveFlux( entity, quadrature[ quadpoint ], gradient, flux );
     *  \endcode
     */
    template< class EntityType, class QuadratureType >
    inline void diffusiveFlux ( const EntityType &entity,
                                const QuadratureType &quadrature,
                                const int quadPoint,
                                const JacobianRangeType &gradient,
                                JacobianRangeType &flux ) const
    {
      asImp().diffusiveFlux( entity, quadrature[ quadPoint ], gradient, flux );
    }
#endif
  };

}

#endif
