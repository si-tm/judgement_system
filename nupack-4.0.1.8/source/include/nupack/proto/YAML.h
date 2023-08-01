/**
Provides conversion between protobuf objects/strings and YAML objects/strings
*/

#pragma once
#include <nupack/proto/public.pb.h>
#include <yaml-cpp/yaml.h>
#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/cat.hpp>

namespace nupack::proto {

std::string proto_to_yaml(std::string const &s);

std::string yaml_to_proto(std::string s);

template <class X>
std::string to_proto_string(X const &x) {
	std::string s;
	if (!x.SerializeToString(&s)) throw std::runtime_error("Invalid protobuf object");
	return s;
}

template <class T>
T from_proto_string(std::string const &s) {
	T raw;
    if (!raw.ParseFromString(s)) throw std::runtime_error("Invalid protobuf string");
	return raw;
}

template <class T>
std::string to_yaml_string(T const &t) {
	std::stringstream ss;
	ss << YAML::Node(t);
	return ss.str();
}

template <class T>
T from_yaml_string(std::string const &s) {
	return YAML::Load(s).as<T>();
} 

struct Dummy {
    template <class T>
    T& operator()(T &t) {return t;}
};

constexpr char const *nice_case(std::string_view s) {
  	if (s == "h_split")   return "H_split";
  	if (s == "n_split")   return "N_split";
  	if (s == "dg_clamp")  return "dG_clamp";
  	if (s == "m_bad")     return "M_bad";
  	if (s == "m_reseed")  return "M_reseed";
  	if (s == "m_reopt")   return "M_reopt";
	return s.data();
}

template <class T>
constexpr bool skip_yaml(T const &t) {return false;}

inline bool skip_yaml(std::string const &t) {return t.empty();}

/******************************************************************************/

template <class T>
static constexpr bool is_allocated = std::is_class_v<T> && !std::is_same_v<T, std::string>;

// YAML::Node encode_map(bool block, std::vector<std::pair<std::string &&, YAML::Node &&>> v) {
// 	YAML::Node node;
// 	if (!block) node.SetStyle(YAML::EmitterStyle::Flow);
// 	for (auto &[k, n] : v) node[std::move(k)] = std::move(n);
// 	return node;
// }

// template <class V>
// YAML::Node encode_list(bool block, V const &v) {
// 	YAML::Node node;
// 	if (!block) node.SetStyle(YAML::EmitterStyle::Flow);
// 	for (auto const &x : v) node.push_back(x);
// 	return node;
// }

// template <class T>
// bool decode_list(YAML::Node const &node, T *t) {
// 	if (!node.IsSequence()) return false;
// 	for (auto const &i : node) t->Add(i.as<typename T::value_type>());
// 	return true;
// }


}

namespace YAML {

/******************************************************************************/

template <class T>
struct convert<google::protobuf::RepeatedField<T>> {
    static Node encode(google::protobuf::RepeatedField<T> const &v) {
        Node out;
        out.SetStyle(YAML::EmitterStyle::Flow);
        for (auto const &x : v) out.push_back(x);
        return out;
    }

    static bool decode(Node const &t, google::protobuf::RepeatedField<T> &out) {
        for (auto const &x : t)	out.Add(x.as<T>());
        return true;
    }
};

/******************************************************************************/

template <class T>
struct convert<google::protobuf::RepeatedPtrField<T>> {
    static Node encode(google::protobuf::RepeatedPtrField<T> const &v) {
        Node out;
        if constexpr(std::is_same_v<T, std::string>) {
            std::size_t n = 0;
            for (auto const &s : v) n = std::max(n, s.size());
            if (n < 40) out.SetStyle(YAML::EmitterStyle::Flow);
        }
        for (auto const &x : v) out.push_back(x);
        return out;
    }

    static bool decode(Node const &t, google::protobuf::RepeatedPtrField<T> &out) {
        for (auto const &x : t) out.Add(x.as<T>());
        return true;
    }
};

/******************************************************************************/

#define NUPACK_YAML_ENUM(T) template<> struct convert<T> { \
	static Node encode(T t) { \
		std::string s = BOOST_PP_CAT(T, _Name)(t); \
		for (char &c : s) c = std::tolower(c); \
		return Node(s); \
	} \
	static bool decode(Node const &n, T &t) { \
		if (!n.IsScalar()) return false; \
		std::string s = n.Scalar(); \
		for (char &c : s) c = std::toupper(c); \
		return BOOST_PP_CAT(T, _Parse)(s, &t); \
	} \
};

/******************************************************************************/

// Only put YAML nodes in that are not 1) Null 2) Empty list 3) Empty string
#define NUPACK_YAML_GET(r, data, index, name) { \
	YAML::Node x(t.name()); \
	if (!(nupack::proto::skip_yaml(t.name()) || x.IsNull() || (x.IsSequence() && x.size() == 0))) \
		node[nupack::proto::nice_case(BOOST_PP_STRINGIZE(name))] = std::move(x); \
}

#define NUPACK_YAML_SET(r, T, index, name) { \
	using U = std::decay_t< decltype(std::declval<T>().name()) >; \
    if (auto const &p = n[nupack::proto::nice_case(BOOST_PP_STRINGIZE(name))]) { \
        if constexpr(nupack::proto::is_allocated<U>) { \
            *C()(out).BOOST_PP_CAT(mutable_, name)() = p.as<U>(); \
        } else { \
            C()(out).BOOST_PP_CAT(set_, name)(p.as<U>()); \
        } \
    } else { \
        if constexpr(std::is_same_v<U, std::string>) C()(out).BOOST_PP_CAT(set_, name)(""); \
    } \
}

#define NUPACK_YAML_WRAP(T, ...) \
template<> struct convert<T> { \
    static Node encode(T const &t) { \
    	Node node; \
        BOOST_PP_SEQ_FOR_EACH_I(NUPACK_YAML_GET,, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)); \
        return node; \
    } \
    template <class C=nupack::proto::Dummy> \
    static bool decode(Node const &n, T &out) { \
    	if (!n.IsMap()) return false; \
        BOOST_PP_SEQ_FOR_EACH_I(NUPACK_YAML_SET, T, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)); \
        return true; \
    } \
};

/******************************************************************************/

#define NUPACK_YAML_VARIANT_SAVE(r, T, index, name) { \
    using U = std::decay_t<decltype(std::declval<T>().name())>; \
    if (type == nupack::proto::nice_case(BOOST_PP_STRINGIZE(name))) { \
        *out.BOOST_PP_CAT(mutable_, name)() = node.as<U>(); \
        return true; \
    } \
}

#define NUPACK_YAML_VARIANT_LOAD(r, T, index, name) if (t.BOOST_PP_CAT(has_, name)()) {Node o(t.name()); o["type"] = BOOST_PP_STRINGIZE(name); return o;}

#define NUPACK_YAML_VARIANT(T, ...) \
template <> \
struct convert<T> { \
    static Node encode(T const &t) { \
        BOOST_PP_SEQ_FOR_EACH_I(NUPACK_YAML_VARIANT_LOAD, T, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)); \
        throw std::runtime_error("Invalid variant"); \
    } \
    static bool decode(Node const &node, T &out) { \
    	std::string const &type = node["type"].Scalar(); \
        BOOST_PP_SEQ_FOR_EACH_I(NUPACK_YAML_VARIANT_SAVE, T, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
        return false; \
    } \
};

/******************************************************************************/

#define NUPACK_WRAP(cls, ...) NUPACK_YAML_WRAP(nupack::proto::cls, __VA_ARGS__)
#define NUPACK_WRAP_VARIANT(cls, ...) NUPACK_YAML_VARIANT(nupack::proto::cls, __VA_ARGS__)
#define NUPACK_WRAP_ENUM(cls) NUPACK_YAML_ENUM(nupack::proto::cls)

#include "Macros.h"

#undef NUPACK_WRAP
#undef NUPACK_WRAP_VARIANT
#undef NUPACK_WRAP_ENUM

}

