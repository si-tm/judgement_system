/** \file IO.cc
 * @brief Defines IO functions and stream operators for many basic objects, and parameter read-in
 */

#include <nupack/types/PairList.h>
#include <nupack/types/Sequence.h>
#include <nupack/types/Alphabet.h>
#include <nupack/types/IO.h>
#include <nupack/types/LRU.h>
#include <nupack/types/Complex.h>
#include <nupack/types/Named.h>

#include <nupack/reflect/Serialize.h>
#include <nupack/reflect/Print.h>
#include <nupack/iteration/Patterns.h>
#include <nupack/algorithms/Utility.h>

#include <stack>
#include <boost/functional/hash.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <charconv>
#include <unordered_set>

#include <ctype.h>
#include <cctype>

namespace nupack {

/******************************************************************************************/

auto sequence_separator() {return boost::is_any_of(",+ \n\t");}

namespace io {

/******************************************************************************************/

std::string run_length_encoding(std::string s, std::size_t min) noexcept {
    if (min < 2) return s;
    char *out = s.data(), *end = out + s.size();

    for (auto b = s.begin(); b != s.end();) {
        char const c = *(out++) = *(b++);

        std::size_t n = 0;
        for (; b != s.end() && *b == c; ++b) ++n; // count run length

        if (n) {
            if (n < min-1) {
                std::fill(out, out + n, c);
                out += n;
            } else {
                out = std::to_chars(out, end, n+1).ptr;
            }
        }
    }

    s.erase(out - s.data() + s.begin(), s.end());
    s.shrink_to_fit();
    return s;
}

bool iequals(std::string_view s, std::string_view t) noexcept {return boost::iequals(s, t);}

std::shared_ptr<std::ostream> default_out = std::shared_ptr<std::ostream>(&std::cout, NoOp());
std::mutex default_out_guard;

bool is_on_next_line(std::istream &is, string const &s) {
    string line;
    auto cur = is.tellg();
    if (!is.good()) return false;
    std::getline(is, line);
    bool ret = (line.find(s) != line.npos);
    is.seekg(cur);
    return ret;
}

/******************************************************************************************/

std::istream & go_to_number(std::istream &is) {
    auto cur = is.tellg();
    string line;
    while (std::getline(is, line)) {
        auto number_pos = line.find_first_of("+-.1234567890");
        if (number_pos != line.npos)
            if (line.find_first_not_of(" \n\r\t") == number_pos) return is.seekg(cur);
        cur = is.tellg();
    }
    throw std::runtime_error("Reached end of file while looking for numbers");
}

/******************************************************************************************/

std::istream & goto_line_after(std::istream &is, string const &s){
    NUPACK_ASSERT(is.good());
    string line;
    while (std::getline(is, line) && is.good()) {
        if (line.find(s) != line.npos) return is;
    }
    throw std::runtime_error("Reached end of file while looking for: " + s);
}

/******************************************************************************************/

std::istream & skip_comments(std::istream &is, string const &comment_start) {
    NUPACK_ASSERT(is.good());
    string line;
    while (is_on_next_line(is, comment_start)) {
        if (!is.good()) throw std::runtime_error("Reached end of file while skipping comments");
        std::getline(is, line);
    }
    return is;
}

/******************************************************************************************/

string peek(std::istream &is){
    string line;
    auto cur = is.tellg();
    std::getline(is, line);
    is.seekg(cur);
    return line;
}

/******************************************************************************************/

std::size_t decoded_length(std::string_view s) noexcept {
    std::size_t out = 0;
    run_length_decode(s, [&](Ignore, std::size_t n) {out += n;});
    return out;
}

/******************************************************************************************/

std::size_t decoded_length_without_separators(std::string_view s) noexcept {
    std::size_t out = 0;
    run_length_decode(s, [&, sep=sequence_separator()](char c, std::size_t n) {if (!sep(c)) out += n;});
    return out;
}

/******************************************************************************************/

/// Convert dot-parens to pair array
void to_pairs(iseq *v, iseq const n, char const *dp) {
    auto const sep = sequence_separator();
    iseq i = 0, s = n - 1;
    // Number of times to repeat a character, 1 if not given
    auto repeat = [](auto &c) {
        auto d = std::strtoul(++c, const_cast<char **>(&c), 10);
        return d ? d : 1u;
    };
    char const *c = dp;
    while (*c) {
        if (*c == '(') {
            for (auto d = repeat(c); d--; ++i) v[s--] = i; // push on stack
        } else if (*c == ')') {
            auto d = repeat(c);
            if (s + d >= n) NUPACK_ERROR("unmatched ) parenthesis", dp, i);
            for (; d--; ++i) {
                v[i] = v[++s]; // pop off stack
                v[v[i]] = i;
            }
        } else if (*c == '.') {
            for (auto d = repeat(c); d--; ++i) v[i] = i;
        } else if (sep(*c)) {
            ++c;
        } else {
            auto index = c - dp;
            NUPACK_ERROR("bad dot-parens character", dp, index, *c, int(*c));
        }
    };
    if (s != n - 1) NUPACK_ERROR("unmatched ( parenthesis");
    if (i != n) NUPACK_ERROR("dot-parens-plus parsing failed");
}

}

/******************************************************************************************/

string Alphabet::to_string(View<Base const *> v, std::size_t n) const {
    string out(v.size(), '_');
    if (data) zip(out, v, [&](char &c, Base b) {c = data->letters[b.value];});
    return io::run_length_encoding(std::move(out), n);
}

string Alphabet::to_string(View<Wildcard const *> v, std::size_t n) const {
    string out(v.size(), '_');
    if (data) zip(out, v, [&](char &c, Wildcard b) {
        if (auto it = data->wildcard_letters.find(b); it != data->wildcard_letters.end())
            c = it->second;;
    });
    return io::run_length_encoding(std::move(out), n);
}

Wildcard Alphabet::complement(Wildcard w) const noexcept {
    WildcardIndex out = 0;
    for (BaseIndex i = 0; w.value; ++i, w.value >>= 1) if (w.value & 1)
        out |= 1 << +complement(Base::from_index(i));
    return Wildcard{out};
}

/******************************************************************************************/

void disable_noncomplement_closing(BasePairing &p, Alphabet const &a) {
    for (auto i : a.all()) for (auto j : a.all())
        if (i != a.complement(j)) p.closing.at(p.length * +i + +j) = false;
}

/******************************************************************************************/

// s = ''.join(chr(a) for a in range(ord('a'), ord('z')+1)) + ''.join(chr(a) for a in range(ord('A'), ord('Z')+1)) + "#$&~|!?/=:;`"
// s = 'acgtACGT' + ''.join(c for c in s if c not in 'acgtACGT') + '_'
// print(s)
char Base::raw_char() const noexcept {
    static_assert(capacity == 64);
    return "acgtACGTbdefhijklmnopqrsuvwxyzBDEFHIJKLMNOPQRSUVWXYZ#$&~|!?/=:;`_"[value];
}

// print([s.index(chr(i)) if chr(i) in s else 64 for i in range(128)])
Base Base::from_raw_char(char c) noexcept {
    static_assert(null().value == 64);
    static constexpr BaseIndex const lookup[128] =
        {64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
         64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 57, 64, 52, 53, 64, 54, 64, 64, 64,
         64, 64, 64, 64, 64, 59, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 61, 62, 64, 60, 64,
         58, 64, 4, 30, 5, 31, 32, 33, 6, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 7,
         46, 47, 48, 49, 50, 51, 64, 64, 64, 64, 64, 63, 0, 8, 1, 9, 10, 11, 2, 12, 13, 14,
         15, 16, 17, 18, 19, 20, 21, 22, 23, 3, 24, 25, 26, 27, 28, 29, 64, 56, 64, 55, 64};
    int const i = c;
    return from_index(0 <= i && i < 128 ? lookup[i] : 64);
}

/******************************************************************************************/

void Sequence::load_repr(std::string_view s) {
    *this = Sequence(io::decoded_length(s), [&](auto p, Ignore) {
        io::run_length_decode(s, [&](char c, std::size_t n) {
            p = std::fill_n(p, n, Base::from_raw_char(c));
        });
    });
}

string raw_sequence_string(View<Base const *> v, std::size_t n) {
    string out;
    out.reserve(v.size());
    for (Base b : v) out.push_back(b.raw_char());
    if (n > 2) out = io::run_length_encoding(std::move(out), n);
    return out;
}

/******************************************************************************************/

// s = ''.join(chr(a) for a in range(ord('a'), ord('z')+1)) + ''.join(chr(a) for a in range(ord('A'), ord('Z')+1)) + "#$&~|!?/=:;`"
// p = '_ACMGRSVUWYHKDBN'
// s = p + ''.join(c for c in s if c not in p)
char Wildcard::raw_char() const noexcept {
    return value < 65 ? "_ACMGRSVUWYHKDBNabcdefghijklmnopqrstuvwxyzEFIJLOPQTXZ#$&~|!?/=:;`"[value] : '\0';
}

// print([s.index(chr(i)) if chr(i) in s else 0 for i in range(128)])
Wildcard Wildcard::from_raw_char(char c) noexcept {
    static constexpr WildcardIndex const lookup[128] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 58, 0, 53, 54, 0, 55, 0, 0, 0, 0, 0, 0, 0, 0, 60, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 62, 63, 0, 61, 0, 59, 0, 1, 14, 2, 13, 42, 43, 4, 11, 44, 45, 12, 46,
        3, 15, 47, 48, 49, 5, 6, 50, 8, 7, 9, 51, 10, 52, 0, 0, 0, 0, 0, 64, 16, 17, 18, 19,
        20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
        41, 0, 57, 0, 56, 0};
    return from_index(lookup[reinterpret_cast<unsigned char const &>(c)]);
}

std::ostream &operator<<(std::ostream &os, Wildcard b) {return os << b.raw_char();}

/******************************************************************************************/

void Domain::load_repr(std::string_view s) {
    *this = Domain(io::decoded_length(s), [&](auto p, Ignore) {
        io::run_length_decode(s, [&](char c, std::size_t n) {
            p = std::fill_n(p, n, Wildcard::from_raw_char(c));
        });
    });
}

string raw_domain_string(Subdomain v, std::size_t n) {
    std::stringstream ss;
    bool simple = true;
    for (Wildcard w : v) {
        if (char c = w.raw_char()) ss << w.raw_char();
        else {simple = false; ss << '[' << +w << ']';}
    }
    string out = ss.str();
    if (simple && n > 2) out = io::run_length_encoding(std::move(out), n);
    return out;
}

/******************************************************************************************/

bool Alphabet::is_palindromic(Subsequence v) const noexcept {
    return std::equal(v.begin(), v.end(), std::make_reverse_iterator(v.end()), std::make_reverse_iterator(v.begin()),
        [&](Base a, Base b) {return a == complement(b);});
}

bool Alphabet::is_palindromic(Subdomain v) const noexcept {
    return std::equal(v.begin(), v.end(), std::make_reverse_iterator(v.end()), std::make_reverse_iterator(v.begin()),
        [&](Wildcard a, Wildcard b) {return a == complement(b);});
}

/******************************************************************************************/

struct BaseDefinition {
    string letter, input, name, complement, material;
    NUPACK_REFLECT(BaseDefinition, letter, input, name, complement, material);
};

void AlphabetData::load_repr(json const &j) {
    materials = vmap(j.at("materials"), [](auto const &j) {
        MaterialData m;
        auto const &v = j.at("letter").template get<string>();
        NUPACK_REQUIRE(len(v), ==, 1, "material string must be exactly one letter", j);
        m.prefix = v[0];
        j.at("name").get_to(m.name);
        fill(m.bases, Base::invalid());
        fill(m.wildcards, Wildcard::null());
        m.bases['_'] = Base::null();
        m.wildcards['_'] = Wildcard::null();
        return m;
    });
    NUPACK_REQUIRE(len(materials), >, 0, "At least one material must be specified");

    auto const defs = j.at("bases").get<vec<BaseDefinition>>();
    NUPACK_REQUIRE(len(defs), <=, Base::capacity, "requested alphabet too long");

    fill(letters, '_');
    fill(complements, Base::invalid());
    fill(bases, Base::invalid());

    for (auto const &def : defs) {
        NUPACK_REQUIRE(len(def.letter), ==, 1, "exactly one letter for each base must be specified", def);
        NUPACK_REQUIRE(len(def.complement), ==, 1, "exactly one complement for each base must be specified", def);
    }

    izip(defs, [&](auto i, auto const &def) {
        letters[i] = def.letter[0];
        auto mat = find_if(materials, [&](auto const &m) {return m.name == def.material;});
        NUPACK_ASSERT(mat != materials.end(), "material not found in list of materials", def);

        auto c = find_if(defs, [&](auto const &d) {return d.letter[0] == def.complement[0];});
        NUPACK_ASSERT(c != defs.end(), "complement not found in list of bases", def);
        complements[i] = Base::from_index(c - defs.begin());

        for (char l : def.input) mat->bases[l] = Base::from_index(i);
        Wildcard const w(WildcardIndex(1) << i);
        // for (char b : def.input) mat->wildcards[b] = w;
        wildcard_letters.emplace(w, def.letter[0]);

        NUPACK_ASSERT(bases[def.letter[0]] == Base::invalid(), "Base letter has already been specified", def.letter[0]);
        bases[def.letter[0]] = Base::from_index(i);
        base_names.emplace_back(def.name);
    });

    for (auto const &[k, v] : j.at("wildcards").items()) {
        NUPACK_REQUIRE(k.size(), ==, 1, "Expected one character in wildcard key");
        WildcardIndex i = 0;
        for (char c : v.get<string>()) i |= WildcardIndex(1) << +bases[c];
        wildcard_letters.emplace(Wildcard{i}, k[0]);
        for (auto &m : materials) m.wildcards[k[0]] = Wildcard{i};
    }
    zip(materials, j.at("materials"), [&](auto &m, auto const &j) {
        for (auto const &[k, v] : j.at("wildcards").items()) {
            NUPACK_REQUIRE(k.size(), ==, 1, "Expected one character in wildcard key");
            WildcardIndex i = 0;
            for (char c : v.template get<string>()) i |= WildcardIndex(1) << +bases[c];
            wildcard_letters.emplace(Wildcard{i}, k[0]);
            m.wildcards[k[0]] = Wildcard{i};
        }
    });
}

/******************************************************************************************/

json AlphabetData::save_repr() const {
    json j;

    j["bases"] = vmap(range(length()), [&](auto b) -> json {
        auto material = find_if(materials, [&](auto const &m) {return contains(m.bases, Base::from_index(b));});
        NUPACK_ASSERT(material != materials.end(), "missing material");
        string input;
        izip(material->bases, [&](auto i, Base c) {if (Base::from_index(b) == c) input.push_back(i);});
        return {
            {"letter", string(1, letters[b])},
            {"name", base_names.at(b)},
            {"complement", string(1, letters[+complements[b]])},
            {"material", material->name},
            {"input", std::move(input)}
        };
    });

    auto global_wildcards = json::object();
    j["materials"] = vmap(materials, [&](auto const &m) -> json {
        auto wildcards = json::object();
        izip(m.wildcards, [&](auto c, Wildcard w) {
            if (w == Wildcard::null()) return;
            string s;
            for (auto i : w.indices()) s.push_back(letters[i]);
            bool global = all_of(materials, [&](auto const &m) {return m.wildcards[c] == w;});
            (global ? global_wildcards : wildcards)[string(1, char(c))] = std::move(s);
        });
        return {
            {"name", m.name},
            {"letter", string(1, m.prefix)},
            {"wildcards", std::move(wildcards)}
        };
    });
    j["wildcards"] = std::move(global_wildcards);

    return j;
}

/******************************************************************************************/

BasePairing load_pairing(Alphabet const &a, json const &j) {
    BasePairing bp(a.length());
    for (auto const &k : j) {
        auto s = k.get<std::string>();
        NUPACK_REQUIRE(len(s), ==, 2, "Incorrect number of bases given in JSON for base pairs");
        Base i = a.get().bases[s[0]], j = a.get().bases[s[1]];
        bp.pairing[bp.length * +j + +i] = bp.pairing[bp.length * +i + +j] = true;
        bp.closing[bp.length * +j + +i] = bp.closing[bp.length * +i + +j] = true;
        bp.possible_pairs[+j].emplace_back(i);
        bp.possible_pairs[+i].emplace_back(j);
    }
    for (auto &p : bp.possible_pairs) {
        p = unique_sorted(std::move(p));
        p.shrink_to_fit();
    }
    return bp;
}

/******************************************************************************************/

json Alphabet::save_repr() const {
    return data ? data->save_repr() : json();
}

void Alphabet::load_repr(json const &j) {
    if (j.is_null()) return;
    data = std::make_shared<AlphabetData>();
    data->load_repr(j);
    m_length = data->length();
}

/******************************************************************************************/

Base Alphabet::first_base(char c) const {return material(0).bases[c];}
Wildcard Alphabet::first_wildcard(char c) const {return material(0).wildcards[c];}

Sequence Alphabet::sequence(std::string_view input) const {
    std::size_t length = 0;
    auto const &materials = get().materials;
    auto material = materials.begin();
    io::run_length_decode(input, [&](char character, std::size_t n) {
        if (Base b = material->bases[character]; b != Base::invalid()) { // look for base within current material
            length += n;
        } else { // look for new material
            material = find_if(materials, [character](auto const &m) {return m.prefix == character;});
            if (material == materials.end()) NUPACK_ERROR("Invalid input string for sequence", input, character, int(character));
        }
    });
    return Sequence(length, [&](Base *p, Ignore) {
        auto material = materials.begin();
        io::run_length_decode(input, [&](char character, std::size_t n) {
            if (Base b = material->bases[character]; b != Base::invalid()) { // look for base within current material
                p = std::fill_n(p, n, b);
            } else { // look for new material
                material = find_if(materials, [character](auto const &m) {return m.prefix == character;});
            }
        });
    });
}

Domain Alphabet::domain(std::string_view input) const {
    std::size_t length = 0;
    io::run_length_decode(input, [&](char character, std::size_t n) {
        auto const &materials = get().materials;
        auto material = materials.begin();
        if (Wildcard b = material->wildcards[character]; b != Wildcard::null() || character == '_') { // look for base within current material
            length += n;
        } else { // look for new material
            material = find_if(materials, [character](auto const &m) {return m.prefix == character;});
            if (material == materials.end()) NUPACK_ERROR("Invalid input string for domain", input, character, int(character));
        }
    });
    return Domain(length, [&](Wildcard *p, Ignore) {
        auto const &materials = get().materials;
        auto material = materials.begin();
        io::run_length_decode(input, [&](char character, std::size_t n) { // look for base within current material
            if (Wildcard b = material->wildcards[character]; b != Wildcard::null() || character == '_') {
                p = std::fill_n(p, n, b);
            } else { // look for new material
                material = find_if(materials, [character](auto const &m) {return m.prefix == character;});
            }
        });
    });
}

char Alphabet::operator()(Wildcard b) const noexcept {
    if (!data) return '\0';
    if (auto it = data->wildcard_letters.find(b); it != data->wildcard_letters.end())
        return it->second;
    return '\0';
}

// DNA default alphabet
Alphabet const DNA(json::parse(R"xxx({
    "bases": [
        {"letter": "A", "name": "DNA-A", "complement": "T", "material": "DNA", "input": "Aa"},
        {"letter": "C", "name": "DNA-C", "complement": "G", "material": "DNA", "input": "Cc"},
        {"letter": "G", "name": "DNA-G", "complement": "C", "material": "DNA", "input": "Gg"},
        {"letter": "T", "name": "DNA-T", "complement": "A", "material": "DNA", "input": "UuTt"}
    ],
    "materials": [
        {"name": "DNA", "letter": "d", "wildcards": {}}
    ],
    "wildcards": {"A": "A", "C": "C", "G": "G", "U": "T", "T": "T", "R": "AG", "M": "AC", "S": "CG", "W": "AT", "K": "GT", "Y": "CT", "V": "ACG", "H": "ACT", "D": "AGT", "B": "CGT", "N": "ACGT",
                  "a": "A", "c": "C", "g": "G", "u": "T", "t": "T",            "m": "AC", "s": "CG", "w": "AT", "k": "GT", "y": "CT", "v": "ACG", "h": "ACT",             "b": "CGT", "n": "ACGT"}
})xxx"));

// RNA default alphabet
Alphabet const RNA(json::parse(R"xxx({
    "bases": [
        {"letter": "A", "name": "RNA-A", "complement": "U", "material": "RNA", "input": "Aa"},
        {"letter": "C", "name": "RNA-C", "complement": "G", "material": "RNA", "input": "Cc"},
        {"letter": "G", "name": "RNA-G", "complement": "C", "material": "RNA", "input": "Gg"},
        {"letter": "U", "name": "RNA-U", "complement": "A", "material": "RNA", "input": "UuTt"}
    ],
    "materials": [
        {"name": "RNA", "letter": "r", "wildcards": {}}
    ],
    "wildcards": {"A": "A", "C": "C", "G": "G", "U": "U", "T": "U", "A": "A", "C": "C", "G": "G", "U": "U", "T": "U", "R": "AG", "M": "AC", "S": "CG", "W": "AU", "K": "GU", "Y": "CU", "V": "ACG", "H": "ACU", "D": "AGU", "B": "CGU", "N": "ACGU",
                  "a": "A", "c": "C", "g": "G", "u": "U", "t": "U", "a": "A", "c": "C", "g": "G", "u": "U", "t": "U",            "m": "AC", "s": "CG", "w": "AU", "k": "GU", "y": "CU", "v": "ACG", "h": "ACU",             "b": "CGU", "n": "ACGU"}
})xxx"));

// Mixed RNA DNA alphabet
Alphabet const RNADNA(json::parse(R"xxx({
    "bases": [
        {"letter": "A", "name": "DNA-A", "complement": "T", "material": "DNA", "input": "Aa"},
        {"letter": "C", "name": "DNA-C", "complement": "G", "material": "DNA", "input": "Cc"},
        {"letter": "G", "name": "DNA-G", "complement": "C", "material": "DNA", "input": "Gg"},
        {"letter": "T", "name": "DNA-T", "complement": "A", "material": "DNA", "input": "Tt"},
        {"letter": "a", "name": "RNA-A", "complement": "u", "material": "RNA", "input": "Aa"},
        {"letter": "c", "name": "RNA-C", "complement": "g", "material": "RNA", "input": "Cc"},
        {"letter": "g", "name": "RNA-G", "complement": "c", "material": "RNA", "input": "Gg"},
        {"letter": "u", "name": "RNA-U", "complement": "a", "material": "RNA", "input": "Uu"}
    ],
    "materials": [
        {"name": "DNA", "letter": "d", "wildcards": {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "M": "AC", "S": "CG", "W": "AT", "K": "GT", "Y": "CT", "V": "ACG", "H": "ACT", "D": "AGT", "B": "CGT", "N": "ACGT"}},
        {"name": "RNA", "letter": "r", "wildcards": {"A": "a", "C": "c", "G": "g", "U": "u", "R": "ag", "M": "ac", "S": "cg", "W": "au", "K": "gu", "Y": "cu", "V": "acg", "H": "acu", "D": "agu", "B": "cgu", "N": "acgu"}}
    ],
    "wildcards": {}
})xxx"));

// Mixed RNA 2OMe-RNA alphabet
Alphabet const RNA2OMeRNA(json::parse(R"xxx({
    "bases": [
        {"letter": "A", "name": "RNA-A", "complement": "U", "material": "RNA", "input": "Aa"},
        {"letter": "C", "name": "RNA-C", "complement": "G", "material": "RNA", "input": "Cc"},
        {"letter": "G", "name": "RNA-G", "complement": "C", "material": "RNA", "input": "Gg"},
        {"letter": "U", "name": "RNA-U", "complement": "A", "material": "RNA", "input": "TtUu"},
        {"letter": "a", "name": "2OMe-RNA-A", "complement": "u", "material": "2OMe-RNA", "input": "Aa"},
        {"letter": "c", "name": "2OMe-RNA-C", "complement": "g", "material": "2OMe-RNA", "input": "Cc"},
        {"letter": "g", "name": "2OMe-RNA-G", "complement": "c", "material": "2OMe-RNA", "input": "Gg"},
        {"letter": "u", "name": "2OMe-RNA-U", "complement": "a", "material": "2OMe-RNA", "input": "Uu"}
    ],
    "materials": [
        {"name": "RNA", "letter": "r", "wildcards": {"A": "A", "C": "C", "G": "G", "U": "U", "T": "U", "R": "AG", "M": "AC", "S": "CG", "W": "AU", "K": "GU", "Y": "CU", "V": "ACG", "H": "ACU", "D": "AGU", "B": "CGU", "N": "ACGU"}},
        {"name": "2OMe-RNA", "letter": "m", "wildcards": {"A": "a", "C": "c", "G": "g", "U": "u", "T": "u", "R": "ag", "M": "ac", "S": "cg", "W": "au", "K": "gu", "Y": "cu", "V": "acg", "H": "acu", "D": "agu", "B": "cgu", "N": "acgu"}}
    ],
    "wildcards": {}
})xxx"));

Sequence operator"" _4(char const *s, std::size_t n) {
    return DNA.sequence(std::string_view(s, n));
}

/******************************************************************************************/

Sequence& Sequence::operator+=(Subsequence const &s) {
    return *this = Sequence(size() + s.size(), [&](Base *p, auto n) {
        p = std::copy(begin(), end(), p);
        std::copy(s.begin(), s.end(), p);
    });
}

Domain& Domain::operator+=(Subdomain const &s) {
    return *this = Domain(size() + s.size(), [&](Wildcard *p, auto n) {
        p = std::copy(begin(), end(), p);
        std::copy(s.begin(), s.end(), p);
    });
}

/******************************************************************************************/

vec<string> split_sequence_string(string_view s2) {
    vec<string> strs;
    string s(s2);
    boost::trim_if(s, sequence_separator());
    boost::split(strs, s, sequence_separator(), boost::token_compress_on);
    return strs;
}

/******************************************************************************************/

vec<std::array<typename PairList::value_type, 4>> PairList::pseudoknots() const {
    vec<std::array<value_type, 4>> out;
    for_pseudoknots(*this, [&](auto ...is) {out.push_back({static_cast<value_type>(is)...});});
    return out;
}

/******************************************************************************************/

PairList PairList::with_null_bases(vec<iseq> const &strand_lengths) const {
    constexpr auto null = std::numeric_limits<iseq>::max();
    data_type out(len(values) + 2 * len(strand_lengths), null);
    // Copy in old values, skip null bases
    auto d = values.begin();
    auto o = out.begin() + 1;
    for (auto l : strand_lengths) {
        std::copy_n(d, l, o);
        d += l;
        o += l + 2;
    }
    // Offset values
    iseq offset = 2 * len(strand_lengths) - 1, end = len(values);
    for (auto l : reversed(strand_lengths)) {
        auto next_end = end - l;
        for (auto &i : out) if (i >= next_end && i < end) i += offset;
        end = next_end;
        offset -= 2;
    }
    // Fix up null bases
    izip(out, [](auto i, auto &j) {if (j == null) j = i;});
    if (!Release) izip(out, [&](auto i, auto &j) {NUPACK_REQUIRE(out[j], ==, i, values, out);});
    return PairList{std::move(out)};
}

/******************************************************************************************/

std::size_t PairList::symmetry() const {
    auto v = values;
    // vector of circular offsets
    izip(v, [s=v.size()](auto i, auto &j) {j = i < j ? j - i : s + j - i;});
    return rotational_symmetry(v);
}

/******************************************************************************************/

NamedSequence NamedSequence::operator~() const {
    return NamedSequence(alphabet.reverse_complement(static_cast<Sequence const &>(*this)), complement_name(name), alphabet);
}

NamedDomain NamedDomain::operator~() const {
    return NamedDomain(alphabet.reverse_complement(static_cast<Domain const &>(*this)), complement_name(name), alphabet);
}

/******************************************************************************************/

template <class V>
string generate_complex_name(V const &v) {
    std::stringstream os;
    char c = '(';
    for (auto const &s : v) {
        os << c << s.name;
        c = '+';
    }
    os << ')';
    return os.str();
}

NamedComplex::NamedComplex(vec<NamedSequence> v, string_view n, real b) :
    strands(std::move(v)), name(n), bonus(b) {
    if (name.empty()) name = generate_complex_name(strands);
}

NamedComplex NamedComplex::lowest_rotation() const {
    NamedComplex out = *this;
    out.strands = ::nupack::lowest_rotation(std::move(out.strands));
    return out;
}

TargetComplex TargetComplex::lowest_rotation() const {
    TargetComplex out = *this;
    auto i = lowest_rotational_order(out.strands);
    std::rotate(out.strands.begin(), out.strands.begin()+i, out.strands.end());
    out.structure.rotate(i);
    return out;
}

ComplexSet ComplexSet::join(View<ComplexSet const *> v) {
    ComplexSet out;
    for (auto const &c : v) {
        cat(out.strands, c.strands);
        cat(out.complexes, c.complexes);
    }
    out.strands = unique_sorted(std::move(out.strands));
    out.complexes = unique_sorted(std::move(out.complexes));
    return out;
}

Tube Tube::join(View<Tube const *> v) {
    Tube out;
    for (auto const &t : v) {
        cat(out.strands, t.strands);
        cat(out.complexes, t.complexes);
    }
    out.strands = unique_sorted(std::move(out.strands));
    out.complexes = unique_sorted(std::move(out.complexes));
    out.concentrations.resize(out.strands.size());
    for (auto const &t : v) {
        zip(t.strands, t.concentrations, [&](auto const &s, real x) {
            if (auto i = binary_search_index(out.strands, s); i != out.strands.size())
                out.concentrations[i] += x;
            else NUPACK_ERROR("Complex strand is not in list of tube strands", s, out.strands);
        });
    }
    return out;
}

/******************************************************************************************/

TargetComplex::TargetComplex(vec<TargetStrand> v, PairList p, string_view n, real b)
    : strands(std::move(v)), name(n), bonus(b) {
    if (!p.empty()) {
        auto const length = sum(strands, [](auto const &s) {return s.nt();});
        NUPACK_REQUIRE(len(p), ==, length, "TargetComplex: structure length must match strand lengths", name, strands, p);
        structure = Structure(std::move(p), prefixes<Nicks>(false,
            indirect_view(strands, [](auto const &s) {return s.nt();})));
    }
    if (name.empty()) name = generate_complex_name(strands);
}

void TargetComplex::resolve(TargetComplex const &b) {
    if (b.bonus) {
        if (bonus) {NUPACK_REQUIRE(bonus, ==, b.bonus);}
        else bonus = b.bonus;
    }
    if (!b.structure.empty()) {
        if (!structure.empty()) {NUPACK_REQUIRE(structure, ==, b.structure, "target structures do not match or were given in multiple alignments");}
        else {
            structure = b.structure;
            strands = b.strands; // copy strands so in same order as structure
        }
    }
    if (b.name != name && name == generate_complex_name(strands))
        name = b.name;
}

/******************************************************************************************/

ComplexSet::ComplexSet(vec<NamedSequence> s, vec<NamedComplex> c)
    : strands(std::move(s)), complexes(std::move(c)) {
    std::unordered_set<NamedSequence> const unique(strands.begin(), strands.end());
    NUPACK_REQUIRE(len(unique), ==, len(strands), "Duplicate strands were given");

    complexes.erase(std::unique(complexes.begin(), complexes.end(), [](auto &a, auto &b) {
        if (lowest_rotation(a.strands) != lowest_rotation(b.strands)) return false;

        if (a.bonus && b.bonus) {NUPACK_REQUIRE(a.bonus, ==, b.bonus, "Equivalent complexes have different bonuses");}
        else if (a.bonus && !b.bonus) b.bonus = a.bonus;
        else if (!a.bonus && b.bonus) a.bonus = b.bonus;

        if (!a.name.empty() && !b.name.empty()) {} //{NUPACK_REQUIRE(a.name, ==, b.name, "Equivalent complexes have different names");}
        else if (!a.name.empty() && b.name.empty()) b.name = a.name;
        else if (a.name.empty() && !b.name.empty()) a.name = b.name;

        return true;
    }), complexes.end());

    for (auto const &c : complexes) for (auto const &s : c.strands)
        NUPACK_ASSERT(unique.count(s), "Complex contains a strand that was not input");
}

/******************************************************************************************/

Tube::Tube(ComplexSet c, vec<real> x, string_view n)
    : ComplexSet(std::move(c)), concentrations(std::move(x)), name(n) {
    NUPACK_REQUIRE(strands.size(), ==, concentrations.size());
    for (auto c : concentrations) NUPACK_REQUIRE(c, >=, 0, "concentration should be non-negative");
}

/******************************************************************************************/

TargetTube::TargetTube(vec<TargetComplex> v, vec<real> c, string_view s) : name(s) {

    for (auto c : concentrations)
        NUPACK_REQUIRE(c, >=, 0, "Concentration should be non-negative");

    std::map<vec<TargetStrand>, std::pair<TargetComplex, real>> lookup;
    izip(v, [&](auto i, auto &x) {
        auto [it, inserted] = lookup.try_emplace(lowest_rotation(x.strands));
        if (inserted) {
            it->second.first = std::move(x);
            it->second.second = i < c.size() ? c[i] : 0;
        } else {
            it->second.first.resolve(std::move(x));
            it->second.second += i < c.size() ? c[i] : 0;
        }
    });

    auto pairs = vmap(lookup, [](auto &p) {return std::move(p.second);});
    n_on_targets = std::stable_partition(pairs.begin(), pairs.end(),
        [](auto const &p) {return p.second != 0;}) - pairs.begin();

    complexes = indirect_view(pairs, [](auto &p) {return std::move(p.first);});
    concentrations = indirect_view(pairs, [](auto const &p) {return p.second;});

    for (auto const &i : on_targets())
        NUPACK_ASSERT(!i.structure.empty(), "On-targets must be given a structure");
}

/******************************************************************************************/

NamedDomainList operator+(NamedDomainList a, NamedDomainList b) {
    a.insert(a.end(), std::make_move_iterator(b.begin()), std::make_move_iterator(b.end()));
    return a;
}

NamedDomainList operator+(NamedDomainList a, NamedDomain b) {a.emplace_back(std::move(b)); return a;}

NamedDomainList operator+(NamedDomain a, NamedDomainList b) {
    b.insert(b.begin(), std::move(a));
    return b;
}

NamedDomainList operator+(NamedDomain a, NamedDomain b) {
    return {std::move(a), std::move(b)};
}

NamedDomainList operator~(NamedDomainList const &v) {
    return indirect_view(reversed(v), std::bit_not<NamedDomain>());
}

/******************************************************************************************/

PairList reduce_pairs(PairList const &p, SequenceList const &v) {
    vec<int> strand;
    izip(v, [&](auto i, auto const &s) {
        strand.emplace_back(-1);
        for (auto x : indices(s)) strand.emplace_back(i);
        strand.emplace_back(-1);
    });
    PairList o;
    for (auto x : p) {
        if (strand[x] != -1) o.values.emplace_back(x - 1 - 2 * strand[x]);
    }
    return o;
}

// If a single pair separates X and Y, return its left index, right index, and +1 or -1 if it is addition or deletion
// Otherwise return 0 0 0
std::tuple<uint, uint, int> single_pair_mismatch(PairList const &X, PairList const &Y) {
    auto const [x, y] = std::mismatch(X.begin(), X.end(), Y.begin());
    if (x != X.end()) {
        auto const a = x - X.begin(); // index of left base
        if (*x == a || *y == a) {
            auto const b = std::max(*x, *y); // index of right base
            if (std::equal(x+1, next(X, b), y+1) && std::equal(next(X, b)+1, X.end(), next(Y, b)+1)) {
                NUPACK_QUICK_ASSERT((X ^ Y) == 2);
                std::tuple<uint, uint, int> out{a, b, *x == a ? 1 : -1};
                if (std::get<2>(out) > 0) {
                    NUPACK_QUICK_ASSERT(X.n_pairs()+1 == Y.n_pairs());
                    NUPACK_QUICK_ASSERT(X[a] == a);
                    NUPACK_QUICK_ASSERT(X[b] == b);
                    NUPACK_QUICK_ASSERT(Y[a] == b);
                    NUPACK_QUICK_ASSERT(Y[b] == a);
                }
                if (std::get<2>(out) < 0) {
                    NUPACK_QUICK_ASSERT(X.n_pairs() == Y.n_pairs()+1);
                    NUPACK_QUICK_ASSERT(Y[a] == a);
                    NUPACK_QUICK_ASSERT(Y[b] == b);
                    NUPACK_QUICK_ASSERT(X[a] == b);
                    NUPACK_QUICK_ASSERT(X[b] == a);
                }
                return out;
            }
        }
    }
    NUPACK_DEBUG_ASSERT((X ^ Y) != 2);
    return {0, 0, 0};
}

/******************************************************************************************/

vec<std::size_t> loop_indices(PairList const &p) {
    vec<std::size_t> base_to_loop(len(p), len(p));
    for (auto i : indices(p)) if (p[i] > i) {
        base_to_loop[i] = i;
        for (auto k = p[i+1]; k != i; k = p[k + 1]) {base_to_loop[k] = i;}
    }
    return base_to_loop;
}

}

namespace std {
    size_t hash<nupack::Sequence>::operator()(nupack::Sequence const &s) const noexcept {
        return nupack::range_hash(nupack::view(s));
    }

    size_t hash<nupack::Domain>::operator()(nupack::Domain const &s) const noexcept {
        return nupack::range_hash(nupack::view(s));
    }

    size_t hash<nupack::Complex>::operator()(nupack::Complex const &m) const noexcept {
        return nupack::range_hash(nupack::view(m));
    }

    size_t hash<nupack::NamedSequence>::operator()(nupack::NamedSequence const &m) const noexcept {
        return hash<string>()(m.name);
    }

    size_t hash<nupack::NamedDomain>::operator()(nupack::NamedDomain const &m) const noexcept {
        return hash<string>()(m.name);
    }

    size_t hash<nupack::AlphabetData>::operator()(nupack::AlphabetData const &a) const noexcept {
        return nupack::range_hash(a.base_names);
    }
}

