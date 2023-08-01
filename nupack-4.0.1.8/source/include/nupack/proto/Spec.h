#pragma once
#include <nupack/proto/public.pb.h>
#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/cat.hpp>
#include <variant>

#include <nupack/reflect/Serialize.h>

namespace nupack {
    NUPACK_NAMESPACE(spec);
}

namespace nupack::spec {

struct Dummy {
    template <class T>
    T& operator()(T &t) {return t;}
};

/******************************************************************************/

template <class T>
struct Mapping {
    using type = T;
    static T load(T t) {return t;}

    static void save(T &out, T t) {out = std::move(t);}
};

template <class T>
T to_proto(typename Mapping<T>::type const &t) {
    T out;
    Mapping<T>::save(out, t);
    return out;
}

template <class T>
T to_proto(typename Mapping<T>::type &&t) {
    T out;
    Mapping<T>::save(out, std::move(t));
    return out;
}

template <class T>
auto from_proto(T const &t) {return Mapping<T>::load(t);}

/******************************************************************************/

template <class K, class V>
struct Mapping<google::protobuf::Map<K, V>> {
    using type = std::map<typename Mapping<K>::type, typename Mapping<V>::type>;

    static type load(google::protobuf::Map<K, V> const &m) {
        type out;
        for (auto const &[k, v] : m)
            out.emplace(Mapping<K>::load(k), Mapping<V>::load(v));
        return out;
    }

    static void save(google::protobuf::Map<K, V> &out, type const &t) {
        for (auto const &[k, v] : t)
            out.insert({to_proto<K>(k), to_proto<V>(v)});
    }
};

/******************************************************************************/

template <class T>
struct Mapping<google::protobuf::RepeatedField<T>> {
    using type = vec<typename Mapping<T>::type>;

    static type load(google::protobuf::RepeatedField<T> const &v) {
        type out;
        out.reserve(v.size());
        for (auto const &x : v) out.emplace_back(Mapping<T>::load(x));
        return out;
    }

    static void save(google::protobuf::RepeatedField<T> &out, type const &t) {
        out.Add(t.begin(), t.end());
    }
};


template <class T>
struct Mapping<google::protobuf::RepeatedPtrField<T>> {
    using type = vec<typename Mapping<T>::type>;

    static type load(google::protobuf::RepeatedPtrField<T> const &v) {
        type out;
        out.reserve(v.size());
        for (auto const &x : v) out.emplace_back(Mapping<T>::load(x));
        return out;
    }

    static void save(google::protobuf::RepeatedPtrField<T> &out, type const &t) {
        for (auto const &x : t) out.Add(to_proto<T>(x));
    }
};

template <class P, class T>
struct Custom {
    using type = T;
    static T load(P const &p) {return T::load(p);}
    static void save(P &p, T t) {t.save(p);}
};

template <class T>
static constexpr bool is_allocated = std::is_class_v<T> && !std::is_same_v<T, std::string>;

/******************************************************************************/

#define NUPACK_WRAP_MEMBER(r, T, index, name) typename Mapping<std::decay_t<decltype( std::declval<proto::T>().name() )>>::type name;

#define NUPACK_WRAP_GET(r, data, index, name) BOOST_PP_COMMA_IF(index) Mapping<std::decay_t<decltype(t.name())>>::load(t.name())

#define NUPACK_WRAP_SET(r, T, index, name) \
    {  using U = std::decay_t< decltype(std::declval<proto::T>().name()) >; \
    if constexpr(is_allocated<U>) { \
        Mapping<U>::save( *(C()(out)).BOOST_PP_CAT(mutable_, name)(), name); \
    } else { \
        (C()(out)).BOOST_PP_CAT(set_, name)(name); \
    }}

#define NUPACK_WRAP(T, ...) \
struct T { \
    BOOST_PP_SEQ_FOR_EACH_I(NUPACK_WRAP_MEMBER, T, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
    NUPACK_REFLECT(T, __VA_ARGS__) \
    static T load(proto::T const &t) { \
        return { BOOST_PP_SEQ_FOR_EACH_I(NUPACK_WRAP_GET,, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) }; \
    } \
    template <class C=Dummy> \
    void save(proto::T &out) const { \
        BOOST_PP_SEQ_FOR_EACH_I(NUPACK_WRAP_SET, T, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)); \
    } \
    proto::T dump() const {proto::T p; save(p); return p;} \
}; \
template <> struct Mapping<proto::T> : Custom<proto::T, T> {}

/******************************************************************************/

#define WRAP_VARIANT_ITEM(r, T, index, name) BOOST_PP_COMMA_IF(index) typename Mapping<std::decay_t<decltype(std::declval<proto::T>().name())>>::type

#define WRAP_VARIANT_SAVE(r, T, index, name) { \
    using U = std::decay_t<decltype(std::declval<proto::T>().name())>; \
    if (auto p = std::get_if<typename Mapping<U>::type>(&t)) { \
        Mapping<U>::save(*out.BOOST_PP_CAT(mutable_, name)(), *p); return; \
    }}

#define WRAP_VARIANT_LOAD(r, T, index, name) if (t.BOOST_PP_CAT(has_, name)()) return Mapping<std::decay_t<decltype(std::declval<proto::T>().name())>>::load(t.name());

#define NUPACK_WRAP_VARIANT(T, ...) \
using T = std::variant< BOOST_PP_SEQ_FOR_EACH_I(WRAP_VARIANT_ITEM, T, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) >; \
template <> \
struct Mapping<proto::T> { \
    using type = T; \
    static T load(proto::T const &t) { \
        BOOST_PP_SEQ_FOR_EACH_I(WRAP_VARIANT_LOAD, T, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)); \
        throw std::runtime_error("bad variant"); \
    } \
    static void save(proto::T &out, T const &t) { \
        BOOST_PP_SEQ_FOR_EACH_I(WRAP_VARIANT_SAVE, T, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
        throw std::runtime_error("bad variant"); \
    } \
};

#define NUPACK_WRAP_ENUM(cls) using proto::cls;

#include "Macros.h"

#undef NUPACK_WRAP
#undef NUPACK_WRAP_VARIANT
#undef NUPACK_WRAP_ENUM

/******************************************************************************/

template <class T>
typename Mapping<T>::type load_stream(std::istream &is) {
    T raw;
    if (!raw.ParseFromIstream(&is)) throw std::runtime_error("parse failure");
    return Mapping<T>::load(std::move(raw));
}

template <class T>
typename Mapping<T>::type load_string(std::string const &s) {
    T raw;
    if (!raw.ParseFromString(s)) throw std::runtime_error("parse failure");
    return Mapping<T>::load(std::move(raw));
}

/******************************************************************************/

}
