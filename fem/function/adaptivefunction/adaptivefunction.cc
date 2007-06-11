
namespace Dune {

  //- AdaptiveDiscreteFunction
  // nothing here, everything in adaptivefunction.hh
  
  
  //- AdaptiveLocalFunction
  template <class DiscreteFunctionSpaceImp>
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  AdaptiveLocalFunction(const DiscreteFunctionSpaceType& spc,
                        DofStorageType& dofVec) :
    spc_(spc),
    dofVec_(dofVec),
    values_(),
    numDofs_(0),
    tmp_(0),
    tmpGrad_(0),
    init_(false),
    multipleGeometryTypes_(spc_.multipleGeometryTypes()),
    baseSet_(0),
    en_(0),
    geoType_(0) // init as Vertex 
  {}
  
  template <class DiscreteFunctionSpaceImp>
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  AdaptiveLocalFunction(const ThisType& other) :
    spc_(other.spc_),
    dofVec_(other.dofVec_),
    values_(),
    numDofs_(0),
    tmp_(0.0),
    tmpGrad_(0.0),
    init_(false),
    multipleGeometryTypes_(spc_.multipleGeometryTypes()),
    baseSet_(0),
    en_(0),
    geoType_(0) // init as Vertex 
  {}

  template <class DiscreteFunctionSpaceImp>
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  ~AdaptiveLocalFunction() {}  

  template <class DiscreteFunctionSpaceImp>
  typename AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::DofType&
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  operator[] (const int num) 
  {
    assert(init_);
    assert(num >= 0 && num < numDofs());
    // check that storage (dofVec_) and mapper are in sync:
    assert(dofVec_.size() >= spc_.size());
    return (* (values_[num]));
  }
  
  template <class DiscreteFunctionSpaceImp>
  const typename AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::DofType&
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  operator[] (const int num) const 
  {
    assert(init_);
    assert(num >= 0 && num < numDofs());
    // check that storage (dofVec_) and mapper are in sync:
    assert(dofVec_.size() >= spc_.size());
    return (* (values_[num]));
  }

  template <class DiscreteFunctionSpaceImp>
  int AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  numDofs() const
  {
    assert( numDofs_ == values_.size() );
    return numDofs_;
  }

  template <class DiscreteFunctionSpaceImp>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  evaluate(const DomainType& x, RangeType& ret) const 
  {
    assert(init_);
    ret = 0.0;

    const int numDof = this->numDofs();
    for (int i = 0; i < numDof; ++i) 
    {
      this->baseFunctionSet().evaluate(i, x, tmp_);
      ret.axpy( (*values_[i]) , tmp_ );
    }
  }

  template <class DiscreteFunctionSpaceImp>
  template <class QuadratureType>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  evaluate(const QuadratureType& quad, 
           const int quadPoint, 
           RangeType& ret) const 
  {
    assert(init_);
    ret = 0.0;

    const int numDof = this->numDofs();
    for (int i = 0; i < numDof; ++i) 
    {
      this->baseFunctionSet().evaluate(i, quad,quadPoint, tmp_);
      ret.axpy( (*values_[i]) , tmp_ );
    }
  }

  template <class DiscreteFunctionSpaceImp>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  jacobian(EntityType& en, 
     const DomainType& x, 
     JacobianRangeType& ret) const
  {
    jacobian(x,ret);
  }

  template <class DiscreteFunctionSpaceImp>
  template <class QuadratureType>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  jacobian(EntityType& en, 
           QuadratureType& quad, 
           int quadPoint, 
           JacobianRangeType& ret) const
  {
    jacobian(quad,quadPoint,ret);
  }

  template <class DiscreteFunctionSpaceImp>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  jacobian(const DomainType& x, JacobianRangeType& ret) const
  {
    assert(init_);
    enum { dim = EntityType::dimension };
    typedef typename DiscreteFunctionSpaceImp::GridType::ctype ctype;

    // get jacobian inverse 
    typedef FieldMatrix<ctype, dim, dim> JacobianInverseType;
    const JacobianInverseType& jti = 
          en().geometry().jacobianInverseTransposed(x);

    ret = 0.0;
    const BaseFunctionSetType& bSet = this->baseFunctionSet();
    
    const int numDof = this->numDofs();
    for (int i = 0; i < numDof; ++i) 
    {
      // evaluate gradient on reference element
      bSet.jacobian(i, x, tmpGrad_);

      // apply element specific values 
      for (int l = 0; l < dimRange; ++l) 
      {
        tmpGrad_[l] *= *values_[i];
        jti.umv(tmpGrad_[l], ret[l]);
      }
    }    
  }

  template <class DiscreteFunctionSpaceImp>
  template <class QuadratureType>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::
  jacobian(const QuadratureType& quad, 
           const int quadPoint, 
           JacobianRangeType& ret) const
  {
    assert(init_);
    enum { dim = EntityType::dimension };

    ret = 0.0;
    const BaseFunctionSetType& bSet = this->baseFunctionSet();
    typedef typename DiscreteFunctionSpaceImp::GridType::ctype ctype;
    
    typedef FieldMatrix<ctype, dim, dim> JacobianInverseType;
    const JacobianInverseType& jti = 
      en().geometry().jacobianInverseTransposed(quad.point(quadPoint));

    const int numDof = this->numDofs();
    for (int i = 0; i < numDof; ++i) 
    {
      // evaluate gradient on reference element
      bSet.jacobian(i, quad,quadPoint, tmpGrad_);

      // apply element specific values 
      for (int l = 0; l < dimRange; ++l) 
      {
        tmpGrad_[l] *= *(values_[i]);
        jti.umv(tmpGrad_[l], ret[l]);
      }
    }
  }

  template <class DiscreteFunctionSpaceImp>
  const typename
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::BaseFunctionSetType& 
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::baseFunctionSet() const {
    assert(baseSet_);
    return *baseSet_;
  }
  
  template <class DiscreteFunctionSpaceImp>
  const typename
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::EntityType & 
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp >::en() const 
  {
    assert( en_ );
    return *en_;
  }
  

  // --init
  template <class DiscreteFunctionSpaceImp>
  void AdaptiveLocalFunction<DiscreteFunctionSpaceImp>::
  init(const EntityType& en) 
  {
    // NOTE: if init is false, then LocalFunction has been create before. 
    // if multipleGeometryTypes_ is true, then grid has elements 
    // of different geometry type (hybrid grid) and we have to check geometry
    // type again, if not we skip this part, because calling the entity's
    // geometry method is not a cheep call 
    
    if( !init_ || multipleGeometryTypes_ )
    {
      if( geoType_ != en.geometry().type() )
      {
        baseSet_ = &spc_.baseFunctionSet(en);

        numDofs_ = baseSet_->numBaseFunctions();
        values_.resize(numDofs_);

        init_ = true;
        geoType_ = en.geometry().type();
      }
    }

    // cache entity
    en_ = &en;

    assert( geoType_ == en.geometry().type() );
    const int numOfDof = numDofs(); 
    for (int i = 0; i < numOfDof; ++i) 
    {
      values_[i] = &(this->dofVec_[spc_.mapToGlobal(en, i)]);
    }

    return ;
  }

  template <class DiscreteFunctionSpaceImp>
  template <class QuadratureType>
  inline void 
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp>::
  axpy(const QuadratureType& quad, 
       const int quadPoint, 
       const RangeType& factor)
  {
    const int numDof = this->numDofs();
    for(int i=0; i<numDof; ++i) 
    {
      this->baseFunctionSet().evaluate( i , quad, quadPoint, tmp_ );
      (*values_[i]) += tmp_ * factor;
    }
  }

  template <class DiscreteFunctionSpaceImp>
  template <class QuadratureType>
  inline void 
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp>::
  axpy(const QuadratureType& quad, 
       const int quadPoint, 
       const RangeType& factor1,
       const JacobianRangeType& factor2)
  {
    const int numDof = this->numDofs();

    const JacobianInverseType& jti = 
      en().geometry().jacobianInverseTransposed(quad.point(quadPoint));
    rightMultiply( factor2, jti, factorInv_ );

    for(int i=0; i<numDof; ++i)
    {
      // evaluate gradient on reference element
      this->baseFunctionSet().evaluate(i, quad, quadPoint, tmp_ );
      (*values_[i]) += tmp_ * factor1;
      this->baseFunctionSet().jacobian(i, quad, quadPoint, tmpGrad_);
      for (int l = 0; l < dimRange; ++l) 
      {
        (*values_[i]) += tmpGrad_[l] * factorInv_[l];
      }
    }
  }
  template <class DiscreteFunctionSpaceImp>
  template <class QuadratureType>
  inline void 
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp>::
  axpy(const QuadratureType& quad, 
       const int quadPoint, 
       const JacobianRangeType& factor)
  {
    const int numDof = this->numDofs();

    const JacobianInverseType& jti = 
      en().geometry().jacobianInverseTransposed(quad.point(quadPoint));
    rightMultiply( factor, jti, factorInv_ );

    for(int i=0; i<numDof; ++i)
    {
      // evaluate gradient on reference element
      this->baseFunctionSet().jacobian(i, quad, quadPoint , tmpGrad_);
      for (int l = 0; l < dimRange; ++l) 
      {
        (*values_[i]) += tmpGrad_[l] * factorInv_[l];
      }
    }
  }
  template <class DiscreteFunctionSpaceImp>
  inline void 
  AdaptiveLocalFunction<DiscreteFunctionSpaceImp>::
  rightMultiply(const JacobianRangeType& factor,
                const JacobianInverseType& jInv,
                JacobianRangeType& result) const 
  {
    enum { rows = JacobianRangeType :: rows };
    enum { cols = JacobianInverseType :: rows };
    for (int i=0; i<rows; ++i)
    {
      for (int j=0; j<cols; ++j)  
      {
        result[i][j] = 0;
        for (int k=0; k<cols; ++k)
        {
          result[i][j] += factor[i][k] * jInv[k][j];
        }
      }
    }
  }


  ////////////////////////////////////////////////////////
  //- AdaptiveDiscreteFunction (specialisation)
  ////////////////////////////////////////////////////////
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  AdaptiveDiscreteFunction<
     CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  ~AdaptiveDiscreteFunction() 
  {
  }
  
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  typename AdaptiveDiscreteFunction<
    CombinedSpace<ContainedFunctionSpaceImp, N, p> >::SubDiscreteFunctionType
  AdaptiveDiscreteFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  subFunction(int component) 
  {
    SubSpaceType& subSpace = this->spc_.subSpace(component);
    return SubDiscreteFunctionType(std::string("Subfunction of ")+this->name(),
                                   subSpace,
                                   this->dofStorage());
  }

  //- AdaptiveLocalFunction (Specialisation for CombinedSpace)
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  AdaptiveLocalFunction(const DiscreteFunctionSpaceType& spc,
                        DofStorageType& dofVec) :
    spc_(spc),
    dofVec_(dofVec),
    values_(),
    numDofs_(0),
    cTmp_(0.0),
    cTmpGradRef_(0.0),
    cTmpGradReal_(0.0),
    tmp_(0.0),
    init_(false),
    multipleGeometryTypes_(spc_.multipleGeometryTypes()),
    baseSet_(0),
    en_(0),
    geoType_(0) // init as Vertex 
  {}

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  AdaptiveLocalFunction(const ThisType& other) :
    spc_(other.spc_),
    dofVec_(other.dofVec_),
    values_(),
    numDofs_(0),
    cTmp_(0.0),
    cTmpGradRef_(0.0),
    cTmpGradReal_(0.0),
    tmp_(0.0),
    init_(false),
    multipleGeometryTypes_(spc_.multipleGeometryTypes()),
    baseSet_(0),
    en_(0),
    geoType_(0) // init as Vertex 
  {}

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  ~AdaptiveLocalFunction() {}

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  typename AdaptiveLocalFunction<
    CombinedSpace<ContainedFunctionSpaceImp, N, p> >::DofType&
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  operator[] (const int num) 
  {
    assert(num >= 0 && num < numDofs());
    return *values_[num/N][static_cast<SizeType>(num%N)];
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  const typename AdaptiveLocalFunction<
    CombinedSpace<ContainedFunctionSpaceImp, N, p> >::DofType&
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  operator[] (const int num) const 
  {
    assert(num >= 0 && num < numDofs());
    return *values_[num/N][static_cast<SizeType>(num%N)];
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  int AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  numDofs() const 
  {
    assert( numDofs_ == (values_.size()*N) );
    return numDofs_;
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  evaluate(const DomainType& x, 
           RangeType& result) const
  {
    result = 0.0;
    assert((values_.size()) == this->baseFunctionSet().numDifferentBaseFunctions());
    const int valSize = values_.size();
    for (int i = 0; i < valSize; ++i) 
    {
      // Assumption: scalar contained base functions
      this->baseFunctionSet().evaluateScalar(i, x, cTmp_);
      for (SizeType j = 0; j < N; ++j) 
      {
        result[j] += cTmp_[0]*(*values_[i][j]);
      }
    }
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  template <class QuadratureType>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  evaluate(const QuadratureType& quad, 
           const int quadPoint, 
           RangeType& ret) const
  {
    ret = 0;
    assert((values_.size()) == this->baseFunctionSet().numDifferentBaseFunctions());
    const int valSize = values_.size();
    for (int i = 0; i < valSize; ++i) 
    {
      // Assumption: scalar contained base functions
      this->baseFunctionSet().evaluateScalar(i, quad, quadPoint, cTmp_);
      for (SizeType j = 0; j < N; ++j) 
      {
        ret[j] += cTmp_[0]*(*values_[i][j]);
      }
    }
  }
  
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  jacobian(EntityType& en, 
           const DomainType& x, 
           JacobianRangeType& result) const
  {
    jacobian(x,result);
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  jacobian(const DomainType& x, JacobianRangeType& result) const
  {
    result = 0.0;

    const JacobianInverseType& jInv = 
      en().geometry().jacobianInverseTransposed(x);

    const int numDiffBaseFct = this->baseFunctionSet().numDifferentBaseFunctions();
    for (int i = 0; i < numDiffBaseFct; ++i) 
    {
      this->baseFunctionSet().jacobianScalar(i, x, cTmpGradRef_);
      cTmpGradReal_ = 0.0;
      jInv.umv(cTmpGradRef_[0], cTmpGradReal_[0]);

      for (SizeType j = 0; j < N; ++j) 
      {
        // Assumption: ContainedDimRange == 1
        result[j].axpy(*values_[i][j], cTmpGradReal_[0]);
      }
    }

  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  template<class QuadratureType>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  jacobian(EntityType& en, 
           QuadratureType& quad, 
           int quadPoint, 
           JacobianRangeType& result) const 
  {
    result = 0.0;

    const JacobianInverseType& jInv = 
      en().geometry().jacobianInverseTransposed(quad.point(quadPoint));

    const int numDiffBaseFct = 
      this->baseFunctionSet().numDifferentBaseFunctions();
    for (int i = 0; i < numDiffBaseFct; ++i) 
    {
      this->baseFunctionSet().jacobianScalar(i, quad, quadPoint , cTmpGradRef_);
      cTmpGradReal_ = 0.0;
      jInv.umv(cTmpGradRef_[0], cTmpGradReal_[0]);

      for (SizeType j = 0; j < N; ++j) 
      {
        // Assumption: ContainedDimRange == 1
        result[j].axpy(*values_[i][j], cTmpGradReal_[0]);
      }
    }
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  template<class QuadratureType>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  jacobian(const QuadratureType& quad, 
           const int quadPoint, 
           JacobianRangeType& ret) const 
  {
    jacobian(quad.point(quadPoint), ret);
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  int AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  numDifferentBaseFunctions() const 
  {
    return values_.size();
  }

  // --init
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  void AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  init(const EntityType& en) 
  {
    // NOTE: if init is false, then LocalFunction has been create before. 
    // if fSpace_.multipleGeometryTypes() is true, then grid has elements 
    // of different geometry type (hybrid grid) and we have to check geometry
    // type again, if not we skip this part, because calling the entity's
    // geometry method is not a cheep call 
    
    if( !init_ || multipleGeometryTypes_ )
    {
      if( geoType_ != en.geometry().type() )
      {
        baseSet_ = &spc_.baseFunctionSet(en);

        numDofs_ = baseSet_->numDifferentBaseFunctions();
        values_.resize(numDofs_);
        
        // real dof number is larger 
        numDofs_ *= N;

        init_ = true;
        geoType_ = en.geometry().type();
      }
    }

    assert( geoType_ == en.geometry().type() );

    // cache entity
    en_ = &en;
    
    const int numDDof = values_.size();
    assert( values_.size()*N == numDofs_);
    for (int i = 0; i < numDDof; ++i) 
    {
      // apply local mapping (see adaptivefunction.hh)
      MapLocalDofs<GridType,p>::map(spc_,en,i,dofVec_,values_); 
    } // end for i
  }
  
  // --baseFunctionSet 
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  const typename AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >:: BaseFunctionSetType& 
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  baseFunctionSet() const 
  {
    assert( baseSet_ );
    return *baseSet_;
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  const typename AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >:: EntityType& 
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  en() const 
  {
    assert( en_ );
    return *en_;
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  template <class QuadratureType>
  inline void 
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  axpy(const QuadratureType& quad, 
       const int quadPoint, 
       const RangeType& factor)
  {
    const int numDDof = values_.size();
    for(int i=0; i<numDDof; ++i) 
    {
      this->baseFunctionSet().evaluateScalar(i , quad, quadPoint, cTmp_ );
      for(int j=0; j<N; ++j)
      {
        (*values_[i][j]) += cTmp_[0] * factor[j];
      }
    }
  }

  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  template <class QuadratureType>
  inline void 
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  axpy(const QuadratureType& quad, 
       const int quadPoint, 
       const JacobianRangeType& factor)
  {
    // get jacobian inverse for given point 
    const JacobianInverseType& jInv = 
      en().geometry().jacobianInverseTransposed(quad.point(quadPoint));

    // apply jacobian inverse 
    rightMultiply( factor, jInv, factorInv_ );
  
    const int numDDof = values_.size();
    for(int i=0; i<numDDof; ++i) 
    {
      // evaluate gradient on reference element
      this->baseFunctionSet().jacobianScalar( i, quad, quadPoint , cTmpGradRef_ );
      for (SizeType j = 0; j < N; ++j) 
      {
        (*(values_[i][j])) += cTmpGradRef_[0] * factorInv_[j]; 
      }
    }
  }
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  template <class QuadratureType>
  inline void 
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  axpy(const QuadratureType& quad, 
       const int quadPoint, 
       const RangeType& factor1,
       const JacobianRangeType& factor2)
  {
    // get jacobian inverse for given point 
    const JacobianInverseType& jInv = 
      en().geometry().jacobianInverseTransposed(quad.point(quadPoint));

    // apply jacobian inverse 
    rightMultiply( factor2, jInv, factorInv_ );
  
    const int numDDof = values_.size();
    for(int i=0; i<numDDof; ++i) 
    {
      // evaluate gradient on reference element
      this->baseFunctionSet().evaluateScalar(i , quad, quadPoint, cTmp_ );
      this->baseFunctionSet().jacobianScalar(i , quad, quadPoint, cTmpGradRef_ );
      for (SizeType j = 0; j < N; ++j) 
      {
        (*(values_[i][j])) += cTmp_[0] * factor1[j] +
          cTmpGradRef_[0] * factorInv_[j]; 
      }
    }
  }
  
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  inline void 
  AdaptiveLocalFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  rightMultiply(const JacobianRangeType& factor,
                const JacobianInverseType& jInv,
                JacobianRangeType& result) const 
  {
    enum { rows = JacobianRangeType :: rows };
    enum { cols = JacobianInverseType :: rows };
    for (int i=0; i<rows; ++i)
    {
      for (int j=0; j<cols; ++j)  
      {
        result[i][j] = 0;
        for (int k=0; k<cols; ++k)
        {
          result[i][j] += factor[i][k] * jInv[k][j];
        }
      }
    }
  }
  
} // end namespace Dune
