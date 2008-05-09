#ifndef OBJECTMATRIX_PEER_GENERATOR_SUPPORT_HPP
#define OBJECTMATRIX_PEER_GENERATOR_SUPPORT_HPP

#include <boost/coroutine/generator.hpp>
#include "peer/enumerator.hpp"

namespace peer {

namespace util {
    template<typename Tgen_>
    struct GeneratorEnumeratorAdapter
        : public Enumerator<typename boost::remove_pointer<
                typename Tgen_::result_type>::type&>
    {
        typedef typename boost::remove_pointer<
                typename Tgen_::result_type>::type& result_type;

        GeneratorEnumeratorAdapter(const Tgen_& gen): gen_(gen) {}

        virtual result_type next()
        {
            if (!*gen_)
                throw StopIteration();
            return *gen_();
        }

    private:
        Tgen_ gen_;
    };
} // namespace util

} // namespace peer


#endif /* OBJECTMATRIX_PEER_GENERATOR_SUPPORT_HPP */
