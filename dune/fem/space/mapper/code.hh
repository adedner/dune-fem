#ifndef DUNE_FEM_DOFMAPPER_CODE_HH
#define DUNE_FEM_DOFMAPPER_CODE_HH

#include <algorithm>
#include <cassert>
#include <iostream>

namespace Dune
{

  namespace Fem
  {

    // DofMapperCode
    // -------------

    class DofMapperCode
    {
    protected:
      typedef const unsigned int *ConstIterator;
      typedef unsigned int *Iterator;

      DofMapperCode ( unsigned int numBlocks, unsigned int numDofs )
      {
        code_ = new unsigned int[ size( numBlocks, numDofs ) ];
        code_[ 0 ] = numBlocks;
        code_[ 1 ] = numDofs;
      }

    public:
      DofMapperCode ()
      {
        code_ = new unsigned int[ size( 0, 0 ) ];
        code_[ 1 ] = code_[ 0 ] = 0;
      }

      DofMapperCode ( const DofMapperCode &other )
      {
        code_ = new unsigned int[ other.size() ];
        std::copy( ConstIterator( other.code_ ), other.end(), code_ );
      }

      ~DofMapperCode ()
      {
        delete[] code_;
      }

      const DofMapperCode &operator= ( const DofMapperCode &other )
      {
        if( size() != other.size() )
        {
          delete[] code_;
          code_ = new unsigned int[ other.size() ];
        }
        std::copy( ConstIterator( other.code_ ), other.end(), code_ );
        return *this;
      }

      template< class Functor >
      void operator() ( Functor f ) const
      {
        for( ConstIterator it = begin(); it != end(); )
        {
          const unsigned int gtIndex = *(it++);
          const unsigned int subEntity = *(it++);
          unsigned int nDofs = *(it++);
          f( gtIndex, subEntity, ConstIterator( it ), ConstIterator( it + nDofs ) );
          it += nDofs;
        }
      }

      unsigned int numBlocks () const { return code_[ 0 ]; }
      unsigned int numDofs () const { return code_[ 1 ]; }

      friend std::ostream &operator<< ( std::ostream &out, const DofMapperCode &code )
      {
        out << "{ " << code.numBlocks() << ", " << code.numDofs();
        for( DofMapperCode::ConstIterator it = code.begin(); it != code.end(); ++it )
          out << ", " << *it;
        return out << " }";
      }

    protected:
      ConstIterator begin () const { return code_ + 2; }
      Iterator begin () { return code_ + 2; }
      ConstIterator end () const { return code_ + size(); }
      Iterator end () { return code_ + size(); }

      std::size_t size () const { return size( numBlocks(), numDofs() ); }

      /**@internal
       *
       * @param[in] numBlocks The number of sub-entities which carry DoFs
       *
       * @param[in] numDofs The total number of DoFs
       */
      static std::size_t size ( unsigned int numBlocks, unsigned int numDofs )
      {
        // Format of the code_ array:
        //
        // code_[0]: number of sub-entities with DoFs
        // code_[1]: total number of DoFs
        //
        // For all k = 0 ... (numBlocks-1)
        // (NB: k corresponds ot a sub-entity with DoFs)
        // It follows a variable size block:
        //
        // code_[offset_k + 0]: geometry of the sub-entity
        // code_[offset_k + 1]: local number for given codim (0 ... refElem.size(cd))
        // code_[offset_k + 2]: #DoFs attached to  this sub-entity
        //
        // code_[offset_k + 2 + j]:
        // for all (j = 0 ... #subEntityDofs) the local index of the
        // given DoF, where "local" now means the number of the
        // corresponding local basis function of the bulk element,
        // i.e. not the numbering inside the entity.
        //
        // offset_(k+1) is then just the start of the next block ...
        //
        // Oh boy :)
        return 2 + 3*numBlocks + numDofs;
      }

      unsigned int *code_;
    };



    // DofMapperCodeWriter
    // -------------------

    class DofMapperCodeWriter
    : public DofMapperCode
    {
    public:
      DofMapperCodeWriter ( unsigned int numBlocks, unsigned int numDofs )
      : DofMapperCode( numBlocks, numDofs )
      {}

      const unsigned int &operator[] ( unsigned int i ) const
      {
        assert( (std::ptrdiff_t)i < end() - begin() );
        return begin()[ i ];
      }

      unsigned int &operator[] ( unsigned int i )
      {
        assert( (std::ptrdiff_t)i < end() - begin() );
        return begin()[ i ];
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DOFMAPPER_CODE_HH
