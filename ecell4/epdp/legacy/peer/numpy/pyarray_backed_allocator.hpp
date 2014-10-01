#ifndef OBJECTMATRIX_PEER_NUMPY_PYARRAY_BACKED_ALLOCATOR_HPP
#define OBJECTMATRIX_PEER_NUMPY_PYARRAY_BACKED_ALLOCATOR_HPP

#include <cstddef>
#include <limits>
#include <numpy/arrayobject.h>

namespace peer {

namespace util
{
    namespace detail
    {
        class pyarray_backed_allocator_base
        {
        public:
            typedef std::size_t size_type;
            typedef std::ptrdiff_t difference_type;

        protected:
            struct allocator_state
            {
                bool giveup_ownership;

            protected:
                std::size_t refcount;

            public:
                allocator_state(bool _giveup_ownership = false)
                    : giveup_ownership(_giveup_ownership), refcount(1) {} 

                bool release()
                {
                    return 0 == --refcount;
                }

                allocator_state& add_ref()
                {
                    ++refcount;
                    return *this;
                }
            };

        public:
            pyarray_backed_allocator_base(bool giveup_ownership)
                : state_(new allocator_state(giveup_ownership)) {}

            pyarray_backed_allocator_base(const pyarray_backed_allocator_base& that)
                : state_(&that.state_->add_ref()) {}

            ~pyarray_backed_allocator_base()
            {
                if (state_->release())
                {
                    delete state_;
                }
            }

            void giveup_ownership()
            {
                
                state_->giveup_ownership = true;
            }

        protected:
            void* _nalloc(const size_type sz, const size_type n) const
            {
                if (static_cast<size_type>(static_cast<double>(sz) * n) !=
                        sz * n)
                {
                    throw std::bad_alloc();
                }

                void* retval = PyDataMem_NEW(sz * n);
                if (!retval)
                {
                    throw std::bad_alloc();
                }
                return retval;
            }

            // it's possible that "free" is previously defined as a
            // preprocessor macro.
            void _free(void* ptr) const
            {   
                if (ptr && !state_->giveup_ownership)
                {
                    PyDataMem_FREE(ptr);
                }
            }

        protected:
            allocator_state* state_;
        };
    }

    template<typename T_>
    class pyarray_backed_allocator
        : public detail::pyarray_backed_allocator_base
    {
    public:
        typedef T_ value_type;
        typedef T_* pointer;
        typedef const T_* const_pointer;
        typedef T_& reference;
        typedef const T_& const_reference;

        template<typename Tother_>
        struct rebind
        {
            typedef pyarray_backed_allocator<Tother_> other;
        };

    public:
        pyarray_backed_allocator(bool giveup_ownership = false)
            : pyarray_backed_allocator_base(giveup_ownership) {}

        pyarray_backed_allocator(const pyarray_backed_allocator_base& that)
            : pyarray_backed_allocator_base(that) {}

        pointer address(reference r) const
        {
            return &r;
        }

        const_pointer address(const_reference r) const
        {
            return &r;
        }

        value_type* allocate(size_type n, const void* hint = 0) const
        {
            return reinterpret_cast<value_type*>(this->_nalloc(sizeof(T_), n));
        }

        void construct(T_* p, const T_& src) const
        {
            new(p) T_(src);
        }

        void destroy(T_* p) const
        {
            p->~T_(); // XXX: does this work for PODs?
        }

        size_type max_size() const
        {
            return std::numeric_limits<size_type>::max() / sizeof(T_);
        }

        void deallocate(T_* p, size_type n) const
        {
            this->_free(p);
        }

        bool operator==(const pyarray_backed_allocator_base& rhs)
        {
            return state_->giveup_ownership == rhs.state_->giveup_ownership;
        }

        bool operator!=(const pyarray_backed_allocator_base& rhs)
        {
            return !operator==(rhs);
        }
    };

} // namespace util

} // namespace peer

#endif /* OBJECTMATRIX_PEER_NUMPY_PYARRAY_BACKED_ALLOCATOR_HPP */
