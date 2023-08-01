#pragma once
#include <nupack/types/Sequence.h>

namespace nupack {

/******************************************************************************************/

enum class MoveType {Add, Del, Join};

/******************************************************************************************/

struct BasePairAddition : TotallyOrdered {
    BasePairAddition() = default;
    BasePairAddition(BaseIter b1_, BaseIter b2_, iseq s1_, iseq s2_, real dlogp_, real rate_=0):
        s1(s1_), s2(s2_), b1(b1_), b2(b2_), dlogp(dlogp_), rate(rate_) {
        NUPACK_ASSERT(s1 <= s2);
        if (s1 == s2) NUPACK_ASSERT(b1 <= b2);
    };

    iseq s1, s2;
    BaseIter b1, b2;
    real dlogp, rate;

    NUPACK_REFLECT(BasePairAddition, s1, s2, b1, b2, dlogp, rate);

    friend bool operator< (BasePairAddition const &m1, BasePairAddition const &m2) {
        return std::make_tuple(m1.s1, m1.s2, m1.b1, m1.b2)
             < std::make_tuple(m2.s1, m2.s2, m2.b1, m2.b2);
    }

    friend bool operator== (BasePairAddition const &m1, BasePairAddition const &m2) {
        return std::make_tuple(m1.s1, m1.s2, m1.b1, m1.b2)
            == std::make_tuple(m2.s1, m2.s2, m2.b1, m2.b2);
    }
};

/******************************************************************************************/

struct JoinMove {
    usize o1, o2; // loop 1, loop 2
    uint index; // Index of type of interaction.
    real margin, scale; // cumulative information for choosing a rate

    NUPACK_REFLECT(JoinMove, o1, o2, index, margin, scale);
};

/******************************************************************************************/

/// Represents a base pair addition between two bases in two different complexes
/// Does not include information on which loops the bases are in
struct JoinLoc {
    iseq s;         // index of the sequence this base is on (I think)
    BaseIter b;     // iterator to the base (dereference to get the base)
    
    // free energy change contributed by this base, partial rate constant contributed by this base
    std::optional<std::pair<real, real>> dlogp_propensity;

    JoinLoc() {};
    JoinLoc(iseq s_, BaseIter b_) : s(s_), b(b_) {};
    JoinLoc(iseq s_, BaseIter b_, real dlogp, real h) : s(s_), b(b_), dlogp_propensity(std::pair<real, real>(dlogp, h)) {};

    NUPACK_REFLECT(JoinLoc, s, b, dlogp_propensity);
};

/******************************************************************************************/

struct LoopJoinMove {
    uint seq1, seq2;
    BaseIter base1, base2;

    NUPACK_REFLECT(LoopJoinMove, seq1, seq2, base1, base2);
};

/******************************************************************************************/

struct ComplexJoinMove {
    usize loop1, loop2;
    uint seq1, seq2;
    BaseIter base1, base2;
    real dlogp, rate;

    operator LoopJoinMove() const {return {.seq1=seq1, .seq2=seq2, .base1=base1, .base2=base2};}

    NUPACK_REFLECT(ComplexJoinMove, loop1, loop2, seq1, seq2, base1, base2, dlogp, rate);
};

/******************************************************************************************/

struct BasePairDeletion {
    using is_member_ordered = True;
    real dlogp = 0;
    real rate = 0;

    NUPACK_REFLECT(BasePairDeletion, dlogp, rate);
};

/******************************************************************************************/

/// Small container to hold info for last move a state took
struct StateMove {
    using is_member_ordered = True;
    std::pair<iseq, iseq> base_pair; /// Last base pair formed or deleted
    real rate; /// Rate of the last step
    NUPACK_REFLECT(StateMove, base_pair, rate);
};

/******************************************************************************************/

}
