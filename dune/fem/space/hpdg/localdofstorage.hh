#ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_LOCALDOFSTORAGE_HH
#define DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_LOCALDOFSTORAGE_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <ostream>
#include <vector>

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // LocalDofStorage
      // ---------------

      template< class GlobalKey >
      class LocalDofStorage
      {
        using ThisType = LocalDofStorage< GlobalKey >;

        using container = std::vector< GlobalKey >;

        typedef uint16_t SizeType;
      public:
        /** \brief global key type */
        using value_type = GlobalKey;

        /** \brief iterator type */
        using iterator = typename container::iterator;
        /** \brief iterator type */
        using const_iterator = typename container::const_iterator;

        /** \name Construction
         *  \{
         */

        /** \brief default constructor */
        LocalDofStorage () : size_( 0 ), new_size_( 0 ), active_( true )
        {
        }

        /** \} */

        /** \name Copying and assignment
         *  \{
         */

        /** \brief copy constructor */
        LocalDofStorage ( const ThisType & ) = default;

        /** \brief move constructor */
        LocalDofStorage ( ThisType && ) = default;

        /** \brief assignment operator */
        ThisType &operator= ( const ThisType & ) = default;

        /** \brief move assignment operator */
        ThisType &operator= ( ThisType && ) = default;

        /** \name Iterators
         *  \{
         */

        /** \brief return iterator to beginning */
        iterator begin () { return dofs_.begin(); }

        /** \brief return iterator to beginning */
        const_iterator begin () const { return dofs_.cbegin(); }

        /** \brief return iterator to end */
        iterator end () { return begin() + size(); }

        /** \brief return iterator to end */
        const_iterator end () const { return begin() + size(); }

        /** \name Capacity
         *  \{
         */

        /** \brief return number of dofs */
        SizeType size () const { return size_; }

        /** \brief enlarge dof vector
         *
         *  \param[in]  new_size  minimum new size for dof vector
         *  \param[in]  dof  running dof number (will be incremented if dof vector is enlarged)
         *
         *  \returns (incremented) dof number
         */
        template< class Function >
        Function reserve ( SizeType new_size, Function function )
        {
          // remember new size
          new_size_ = new_size;

          // enlarge dof vector if needed
          const SizeType old_size = dofs_.size();
          if( old_size < new_size_ )
          {
            dofs_.resize( new_size );
            for( SizeType i = old_size; i < new_size; ++i )
              dofs_[ i ] = function();
          }

          return std::move( function );
        }

        /** \brief remove marked dofs from dof vector */
        template< class Function >
        void resize ( Function function )
        {
          if( size_ != new_size_ )
          {
            assert( new_size_ <= dofs_.size() );
            size_ = new_size_;

            const SizeType holes = dofs_.size() - size_;
            std::for_each( dofs_.rbegin(), dofs_.rbegin() + holes, function );

            dofs_.resize( size_ );
            // dofs_.shrink_to_fit();
          }
        }

        /** \brief remove marked dofs from dof vector */
        template< class Function >
        void resize ( Function function, const bool activate )
        {
          resize( function );

          // mark dofs as active/inactive for dof compression
          active_ = activate;
        }

        /** \brief remove marked dofs from dof vector */
        void clear ()
        {
          size_ = 0;
          new_size_ = 0;

          dofs_.resize( size_ );

          // mark dofs as inactive
          active_ = false;
        }

        /** \} */

        /** \name Element access
         *  \{
         */

        /** \brief access element */
        GlobalKey &operator[] ( std::size_t n ) { return dofs_[ n ]; }

        /** \brief access element */
        const GlobalKey &operator[] ( std::size_t n ) const { return dofs_[ n ]; }

        /** \} */

#ifndef DOXYGEN

        friend std::ostream &operator<< ( std::ostream &ostream, const LocalDofStorage &storage )
        {
          const auto &dofs = storage.dofs_;
          const SizeType size = dofs.size();
          for( SizeType i = 0; i < size; ++i )
            ostream << dofs[ i ] << " ";
          ostream << "[" << storage.size_  << "; " << storage.new_size_ << "]";
          return ostream;
        }

#endif // #ifndef DOXYGEN

        bool active () const { return active_; }
        void deactivate ()   { active_ = false; }
        void activate ()     { active_ = true; }

      private:
        // 65536 dofs per element should be sufficient, otherwise increase
        SizeType size_, new_size_;
        bool active_; // flag for dof compression
        std::vector< GlobalKey > dofs_;
      };

    } // namespace hpDG

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_LOCALDOFSTORAGE_HH
