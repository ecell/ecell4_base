#ifndef UTILS_HPP
#define UTILS_HPP

template<template<typename, typename> class TTmapper_>
struct make_get_mapper_mf
{
    template<typename Tkey_, typename Tval_>
    struct meta_type {
        typedef TTmapper_<Tkey_, Tval_> type;
    };
};

#endif /* UTILS_HPP */
