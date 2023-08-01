#pragma once
#include "../iteration/Patterns.h"
#include "../common/Error.h"
#include "../reflect/Print.h"
#include <boost/config.hpp>

namespace nupack::thermo {

// The priorities are designated to reflect logical dependency.
// If A depends on B (for a given i,j element), the priority of A must be higher than B
// Multiple orders would probably work. but this one is fine
using Priority = std::uint_fast8_t;
using Index = int;
using Material = uint;

struct D0 {
    static constexpr Index value = 0; 
    explicit constexpr operator bool() const {return false;}
    template <class I, NUPACK_IF(std::is_integral_v<I>)>
    constexpr operator I() const {return 0;}
    template <class I, NUPACK_IF(std::is_integral_v<I>)>
    friend constexpr I operator+(D0, I i) {return i;}
    template <class I, NUPACK_IF(std::is_integral_v<I>)>
    friend constexpr I operator+(I i, D0) {return i;}
    template <class I, NUPACK_IF(std::is_integral_v<I>)>
    friend constexpr I operator-(I i, D0) {return i;}
};

struct D1 {
    static constexpr Index value = 1; 
    explicit constexpr operator bool() const {return true;}
    
    template <class I, NUPACK_IF(std::is_integral_v<I>)>
    constexpr operator I() const {return 1;}

    template <class I, NUPACK_IF(std::is_integral_v<I>)>
    friend constexpr I operator+(D1, I i) {return ++i;}
    template <class I, NUPACK_IF(std::is_integral_v<I>)>
    friend constexpr I operator+(I i, D1) {return ++i;}
    template <class I, NUPACK_IF(std::is_integral_v<I>)>
    friend constexpr I operator-(I i, D1) {return --i;}
};

static constexpr D0 const d0{};
static constexpr D0 const All{};
static constexpr D1 const d1{};
static constexpr D1 const On{};

inline constexpr D0 operator+(D0, D0) {return {};}
inline constexpr D1 operator+(D1, D0) {return {};}
inline constexpr D1 operator+(D0, D1) {return {};}
inline constexpr Index operator+(D1, D1) {return 2;}

template <class T>
static constexpr bool can_backtrack = decltype(std::decay_t<T>::can_backtrack())::value;

template <std::size_t I>
struct RecursionName;

#define NUPACK_TMP(N, NAME) struct NAME##_t { \
        constexpr auto recursion() const {return *this;} \
        constexpr Priority priority() const {return N;} \
    }; \
    static constexpr auto NAME = NAME##_t();

NUPACK_TMP(0,  X);   // fast interior
NUPACK_TMP(0,  E);   // region bonus


NUPACK_TMP(10,  M2);  // multi pair
NUPACK_TMP(20,  B);   // paired PF
NUPACK_TMP(30,  Z);   // paired plus terminal penalty
NUPACK_TMP(40,  D);   // paired plus terminal penalty if can close helix
NUPACK_TMP(50,  YA);  // mismatch
NUPACK_TMP(60,  YB);  // other mismatch
NUPACK_TMP(70,  MD);  // was CM
NUPACK_TMP(80,  MC);  // was CS
NUPACK_TMP(90,  MCS); // was CSS
NUPACK_TMP(110, MS);  // multi complement. was CBS
NUPACK_TMP(120, CD);  // was C
NUPACK_TMP(130, S);   // paired complement
NUPACK_TMP(140, M);   // multi
NUPACK_TMP(150, MA);   // multi + alpha
NUPACK_TMP(160, Q);   // total PF
NUPACK_TMP(170, N);   // exterior loop

// Mixed Materials
NUPACK_TMP(1,    X_single);
NUPACK_TMP(2,    X_mixed);
NUPACK_TMP(3,    X_singlemix);
NUPACK_TMP(11,   MM_single);
NUPACK_TMP(12,   MM_mixed);
NUPACK_TMP(41,   D_single);
NUPACK_TMP(42,   D_mixed);
NUPACK_TMP(90,   L_single);
NUPACK_TMP(100,  L_mixed);
NUPACK_TMP(102,  MA_single);
NUPACK_TMP(104,  MA_mixed);

NUPACK_TMP(71,  MD_single);  
NUPACK_TMP(72,  MD_mixed);  

NUPACK_TMP(81,  MC_single);  
NUPACK_TMP(82,  MC_mixed);  
NUPACK_TMP(91,  MCS_single); 
NUPACK_TMP(92,  MCS_mixed);
NUPACK_TMP(111,  MS_single);
NUPACK_TMP(112,  MS_mixed);
NUPACK_TMP(121,  CD_new); 
NUPACK_TMP(141,  M_single);
NUPACK_TMP(142,  M_mixed);




// Kinetic stuff
NUPACK_TMP(190, HM); 
NUPACK_TMP(200, MH); 
NUPACK_TMP(210, QH); 
NUPACK_TMP(220, HQ); 
NUPACK_TMP(230, HN); 
NUPACK_TMP(240, NH); 

NUPACK_TMP(170, AH); 
NUPACK_TMP(180, HA); 
NUPACK_TMP(240, LH);
NUPACK_TMP(240, HL);
NUPACK_TMP(240, HLC);
NUPACK_TMP(240, LHC);

NUPACK_TMP(240, HQC);
NUPACK_TMP(240, QHC);

NUPACK_TMP(240, HNC);
NUPACK_TMP(240, NHC);

NUPACK_TMP(240, HQH);
NUPACK_TMP(240, HNH);
NUPACK_TMP(240, HLH);

NUPACK_TMP(240, CHQH);
NUPACK_TMP(240, CHNH);
NUPACK_TMP(240, CHLH);

NUPACK_TMP(240, M1);
NUPACK_TMP(240, M3);
NUPACK_TMP(240, MO);
NUPACK_TMP(240, ND);


// NUPACK_TMP(3,  B);   // paired PF
// NUPACK_TMP(4,  Z);   // paired plus terminal penalty
// NUPACK_TMP(7,  YA);  // mismatch
// NUPACK_TMP(8,  YB);  // other mismatch
// NUPACK_TMP(13, S);   // paired complement
// NUPACK_TMP(16, Q);   // total PF
// NUPACK_TMP(17, N);   // exterior loop

#undef NUPACK_TMP

std::string_view recursion_name(Priority i);

struct Addressed {};

/******************************************************************************************/

// Code holding the current status of the DP calculation
// Negative numbers refer to specific error codes
// Non-negative numbers refer to the index of the first diagonal that has not yet been completed successfully.
struct Stat {
    int value;
    constexpr explicit Stat(int i) : value(i) {}
    NUPACK_REFLECT(Stat, value);

    constexpr bool operator==(Stat const &s) const {return value == s.value;}
    constexpr bool operator!=(Stat const &s) const {return value != s.value;}
    
    static constexpr Stat interrupt() {return Stat(-1);}
    static constexpr Stat bug() {return Stat(-2);}

    uint start_diagonal() const noexcept {return value;}

    friend std::ostream &operator<<(std::ostream &os, Stat const &s) {
        if (s.value == -1) return os << "interrupt";
        if (s.value == -2) return os << "bug";
        return os << s.value;
    }
};

/******************************************************************************************/

class AtomicStat {
    std::atomic<int> value{-3};
public:
    void store_failure(Stat const &s) noexcept {value.store(s.value, std::memory_order_relaxed);}
    Stat load(int desired) const noexcept {
        auto v = value.load(std::memory_order_relaxed);
        return Stat(v == -3 ? desired : v);
    }
    bool should_stop() const noexcept {return value.load(std::memory_order_relaxed) != -3;}
};

/******************************************************************************************/

}
