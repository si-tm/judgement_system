#pragma once
#include <nupack/model/Move.h>
#include <nupack/model/ModelVariants.h>

namespace nupack::kmc {

/******************************************************************************************/

using BasePairList = small_vec<std::pair<Base, Base>, 8>;

/******************************************************************************************/

// From a nick index progressing to a base pair index (inclusive, exclusive), return the column factor of stacking
template <class T, class V, class M>
auto stack_left(V const &v, int const nick, int pair, M const &model) {
    using A = StackMatrix<T, std::plus<T>, std::multiplies<T>>;
    StackColumn<T, std::plus<T>, std::multiplies<T>> x({1, 1}); // initialization that takes care of the nick
    
    if (pair == 0) pair = len(v);

    if (pair < nick) {
        for (auto s : range(nick + 1, len(v))) 
            x = x * factor<A>(len(v[s-1]), len(v[s]), false, stack_model<true>(v, s, model));
        for (auto s : range(pair)) 
            x = x * factor<A>(len(v[s-1]), len(v[s]), false, stack_model<true>(v, s, model));
    } else {
        for (auto s : lrange(nick + 1, pair))
            x = x * factor<A>(len(v[s-1]), len(v[s]), false, stack_model<true>(v, s, model));
    }
    return x;
}

/******************************************************************************************/

// From a base pair index progressing to a nick index (inclusive, exclusive), return the column factor of stacking
template <class T, class V, class M>
auto stack_right(V const &v, int const pair, int nick, M const &model) {
    using A = StackMatrix<T, std::plus<T>, std::multiplies<T>>;
    StackColumn<T, std::plus<T>, std::multiplies<T>> x({1, 0});

    if (nick == 0) nick = len(v);

    if (nick < pair) {
        for (auto s : ~range(nick)) 
            x = factor<A>(len(v[s-1]), len(v[s]), false, stack_model<true>(v, s, model)) * x;
        for (auto s : ~range(pair, len(v))) 
            x = factor<A>(len(v[s-1]), len(v[s]), false, stack_model<true>(v, s, model)) * x;
    } else {
        for (auto s : ~range(pair, nick)) 
            x = factor<A>(len(v[s-1]), len(v[s]), false, stack_model<true>(v, s, model)) * x;
    }
    return x;
}

/******************************************************************************************/

// Same as StackModel but artificially gives 1 for the coaxial stack
template <class M>
struct EndStackModel : StackModel<M, true> {
    EndStackModel(StackModel<M, true> const &o) : StackModel<M, true>(o) {}

    real stack() const  {return this->one();}
};

/******************************************************************************************/

template <class T>
struct FullStackAccumulator {
    using A = StackMatrix<T, std::plus<T>, std::multiplies<T>>;

    // For a base pair attaching to given seq and base index, get the stacking column associated with the BPs BEFORE that base pair
    template <class V, class M>
    auto begin_stack(V const &v, int const nick, int seq, int base, M const &model) const {
        StackColumn<T, std::plus<T>, std::multiplies<T>> x({1, 1}); // initialization that takes care of the nick
        // stack_model only looks backwards, so we can just use it while modifying the sequence length appropriately

        // every pair after the nick but before seq
        if (seq != nick) {
            if (nick < seq) {
                for (auto s : range(nick+1, seq))
                    x = x * factor<A>(len(v[s-1]), len(v[s]), false, stack_model<true>(v, s, model));
            } else {
                for (auto s : range(nick+1, len(v)))
                    x = x * factor<A>(len(v[s-1]), len(v[s]), false, stack_model<true>(v, s, model));
                for (auto s : range(seq))
                    x = x * factor<A>(len(v[s ? s-1 : len(v)-1]), len(v[s]), false, stack_model<true>(v, s, model));
            }
            x = x * factor<A>(len(v[seq ? seq-1 : len(v)-1]), base+1, false, stack_model<true>(v, seq, model));
        }
        return x;
    }

    // For a base pair attaching to given seq and base index, get the stacking column associated with the BPs AFTER that base pair
    template <class V, class M>
    auto end_stack(V const &v, int nick, int seq, int base, M const &model) const {
        StackColumn<T, std::plus<T>, std::multiplies<T>> x({1, 0});
        // if (nick == 0) nick = len(v);
        auto const seq1 = (seq+1) % len(v);

        if (seq1 != nick) {
            if (1 < (len(v) + nick - seq1) % len(v)) {
                auto seq2 = (seq1+1) % len(v);
                if (seq2 < nick) {
                    for (auto s : ~range(seq2, nick))
                        x = factor<A>(len(v[s ? s-1 : len(v)-1]), len(v[s]), false, stack_model<true>(v, s, model)) * x;
                } else {
                    for (auto s : ~range(nick)) 
                        x = factor<A>(len(v[s ? s-1 : len(v)-1]), len(v[s]), false, stack_model<true>(v, s, model)) * x;
                    for (auto s : ~range(seq2, len(v))) 
                        x = factor<A>(len(v[s ? s-1 : len(v)-1]), len(v[s]), false, stack_model<true>(v, s, model)) * x;
                }
            }
            x = factor<A>(len(v[seq])-base, len(v[seq1]), false, EndStackModel<M>(stack_model<true>(v, seq1, model))) * x;
        }

        return x;
    }
};

/******************************************************************************************/

template <class T>
struct StackAccumulator {
    using A = StackMatrix<T, std::plus<T>, std::multiplies<T>>;
    using B = StackColumn<T, std::plus<T>, std::multiplies<T>>;
    vec<B> bcache, ecache;

    template <class V, class M>
    StackAccumulator(V const &v, int const nick, M const &model) {
        bcache.reserve(len(v));
        ecache.reserve(len(v));
        bcache.push_back({1, 1});
        ecache.push_back({1, 0});
        // every pair after the nick but before seq
        for (auto s : range(nick+1, len(v)))
            bcache.emplace_back(bcache.back() * factor<A>(len(v[s-1]), len(v[s]), false, stack_model<true>(v, s, model)));
        for (auto s : range(nick))
            bcache.emplace_back(bcache.back() * factor<A>(len(v[s ? s-1 : len(v)-1]), len(v[s]), false, stack_model<true>(v, s, model)));

        for (auto s : ~range(nick)) // progress backwards from the nick
            ecache.emplace_back(factor<A>(len(v[s ? s-1 : len(v)-1]), len(v[s]), false, stack_model<true>(v, s, model)) * ecache.back());
        for (auto s : ~range(nick + 1, len(v))) // progress backwards from the end until the nick
            ecache.emplace_back(factor<A>(len(v[s ? s-1 : len(v)-1]), len(v[s]), false, stack_model<true>(v, s, model)) * ecache.back());
    }

    // For a base pair attaching to given seq and base index, get the stacking column associated with the BPs BEFORE that base pair
    template <class V, class M>
    B begin_stack(V const &v, int const nick, int seq, int base, M const &model) const {
        if (seq == nick) return bcache[0];
        auto const &prefix = bcache[(seq - nick - 1) % len(v)];
        return prefix * factor<A>(len(v[seq ? seq-1 : len(v)-1]), base+1, false, stack_model<true>(v, seq, model));
    }

    // For a base pair attaching to given seq and base index, get the stacking column associated with the BPs AFTER that base pair
    template <class V, class M>
    B end_stack(V const &v, int nick, int seq, int base, M const &model) const {
        int const N  = len(v);
        auto const seq1 = (seq+1) % N;
        if (seq1 == nick) return ecache[0];

        auto const &suffix = ecache[1 < (N + nick - seq1) % N ?  (N + nick - seq1 - 1) % N : 0];
        return factor<A>(len(v[seq])-base, len(v[seq1]), false, EndStackModel<M>(stack_model<true>(v, seq1, model))) * suffix;
    }
};


/******************************************************************************************/

inline constexpr bool dangle_match(Base a, Base b) {return a.is_null() || a == b;}

/******************************************************************************************/

template <class Model>
auto coaxial_factor(Model const &model, bool which, Base a, Base b, Base c, Base d) {
    return model.boltz(model.coaxial_stack_energy(a, b, c, d) + (which ? model.terminal_penalty(c, d) : model.terminal_penalty(a, b)));
}

/******************************************************************************************/

template <class V, class F>
void for_dangles(V const &v, F &&f) {
    for (Base b : v) f(Base(b));
    f(Base::null());
}

/******************************************************************************************/

// Left column for base i of sequence s to bind from exterior loop v
template <class O, class A, class V, class M>
void accumulate_left_stacks(O &&o, A const &acc, M const &model, BasePairList const &pairs, V const &v, uint nick, uint s, uint i) {
    auto const bases = model.alphabet().all();
    Base const bi = v[s][i];
    // Dangling bases on either side of the new base pair, or null if there is none there
    Base const bi1 = i+2 < len(v[s]) ? v[s][i+1] : Base::null();
    Base const bi0 = 1 < i ? v[s][i-1] : Base::null();

    // This gives the pfunc of every BP up to the new one, depending if the new one stacks or not to the left
    // it looks at every base up to BP s and also the one possibly dangling after
    // This gives the pfunc of every BP up to the new one, depending if the one after the new one stacks or not to the left
    auto const x = acc.begin_stack(v, nick, s, i, model);
    auto const y = acc.end_stack(v, nick, s, i, model);

    auto const dangle = [&model](Base id, Base i, Base j, Base jd) {
        return model.boltz(
            model.terminal_penalty(i, j) +
            (id.is_null() ? (jd.is_null() ? 0 : model.dangle5(i, j, jd)) : 
                            (jd.is_null() ? model.dangle3(id, i, j) : model.terminal_mismatch(id, i, j, jd)))
        );
    };

    // No coaxial stacks
    for (auto j : bases) {
        if (!model.pairing().can_close(bi, j)) {
            o += sq(len(bases) + 1);
        } else for_dangles(bases, [&](Base j1) {
            for_dangles(bases, [&](Base j0) {
                // partition function given the sets of dangles on the left side.
                // PF of the left side given pair to j, j1 and j0 dangling
                
                // x says the additional PF given that j0 is dangle or not.
                // y says the additional PF given that next BP stacks left or not.
                real lo = x(0) * dangle(Base::null(), bi, j, j1) // no dangle to the left of the new BP
                        + (bi0.is_null() ? 0 : x(1) * dangle(bi0, bi, j, j1)); // dangle to the left of the new BP if it's allowed
                real hi = y(0) * (dangle(j0, j, bi, Base::null()) + (bi1.is_null() ? 0 : dangle(j0, j, bi, bi1))) // next pair is not stacking left. then the new one may or may not be able to danlgle
                        + y(1) * (bi1.is_null() ? 0 : (dangle(j0, j, bi, Base::null()) + (i + 3 == len(v[s]) ? 0 : dangle(j0, j, bi, bi1)))); // all possibilities if next pair is stacking left
                o << lo * hi;
            });
        });
    }
    // Coaxial stacks on left, top and bottom
    if (len(v[s]) == 3 && !front(v[s]).is_null() && !back(v[s]).is_null()) {
        for (auto j : bases) {
            if (model.pairing().can_close(bi, j))
                o << x(1) * y(1) * coaxial_factor(model, true, back(v[s ? s-1 : len(v)-1]), front(v[s]), bi, j) * coaxial_factor(model, false, j, bi, back(v[s]), front(v[(s+1) % len(v)])); // BP must stack left, BP after must stack left
            else ++o;
        }
    } else o += len(bases);

    // Coaxial stacks on right, top and bottom
    for (auto bi_ : bases) {
        if (bi_ == bi) {
            real lo = x(0); // BP must not stack left
            real hi = y(0) + (bi1.is_null() ? 0 : y(1)); // BP after may only stack left if a dangle exists
            o << lo * hi; 
        } else ++o;
    }

    // Coaxial stack on top left only
    for (auto j : bases) {
        if (model.pairing().can_close(bi, j) && bi1.is_null() && !v[s][i+1].is_null()) {
            for_dangles(bases, [&](Base j1) {
                o << (x(0) * dangle(Base::null(), bi, j, j1) + x(1) * (bi0.is_null() ? 0 : dangle(bi0, bi, j, j1))) 
                     * y(1) * coaxial_factor(model, false, j, bi, back(v[s]), front(v[(s+1) % len(v)])); // BP must not stack left, BP after must stack left
            });
        } else o += len(bases) + 1;
    }

    // Coaxial stack on bottom left only
    for (auto j : bases) {
        if (model.pairing().can_close(bi, j) && bi0.is_null() && !v[s][i-1].is_null()) {
            for_dangles(bases, [&](Base j0) {
                o << x(1) * coaxial_factor(model, true, back(v[s ? s-1 : len(v)-1]), front(v[s]), bi, j)
                     * (y(0) * (dangle(j0, j, bi, Base::null()) + (bi1.is_null() ? 0 : dangle(j0, j, bi, bi1))) 
                        + (bi1.is_null() ? 0 : (y(1) * (dangle(j0, j, bi, Base::null()) + (i+3 == len(v[s]) ? 0 : dangle(j0, j, bi, bi1)))))); // BP must stack left, BP after can stack left only if a dangle exists.
            });
        } else o += len(bases) + 1;
    }
    // Coaxial stack on top right only
    for (auto bi_ : bases) {
        if (bi != bi_) o += len(bases) + 1;
        else for_dangles(bases, [&](Base bi0_) {
            if (dangle_match(bi0_, bi0))
                o << x(!bi0_.is_null()) * (y(0) + (bi1.is_null() ? 0 : y(1))); // BP must not stack left, BP after may only stack left if a dangle exists
            else ++o;
        });
    }
    // Coaxial stack on bottom right only
    for (auto bi_ : bases) {
        if (bi != bi_) o += len(bases) + 1;
        else for_dangles(bases, [&](Base bi1_) {
            if (dangle_match(bi1_, bi1)) 
                o << x(0) * (y(0) + y(1) * (!bi1.is_null() && (bi1_.is_null() || i+3 != len(v[s])) ? 1 : 0)); // BP must not stack left, BP after may or may not stack left
            else ++o;
        });
    }
    // Coaxial stack top left and bottom right
    if (bi1.is_null() && !v[s][i+1].is_null()) {
        for (auto const &[bi_, bj] : pairs) {
            if (bi == bi_)
                o << x(0) * y(1) * coaxial_factor(model, false, bj, bi, back(v[s]), front(v[(s+1) % len(v)])); // BP must not stack left. BP after must stack left
            else ++o;
        }
    } else o += len(pairs);

    // Coaxial stack bottom left and top right
    if (bi0.is_null() && !v[s][i-1].is_null()) {
        for (auto const &[bi_, bj] : pairs) {
            if (bi == bi_)
                o << x(1) * (y(0) + (bi1.is_null() ? 0 : y(1))) * coaxial_factor(model, true, back(v[s ? s-1 : len(v)-1]), front(v[s]), bi, bj); // BP must stack. BP after may or may not stack left
            else ++o;
        }
    } else o += len(pairs);
}

/******************************************************************************************/

// Right column for base j of sequence s to bind from exterior loop v
template <class O, class A, class V, class M>
void accumulate_right_stacks(O &&o, A &&acc, M const &model, BasePairList const &pairs, V const &v, uint nick, uint s, uint j) {
    auto const bases = model.alphabet().all();
    Base const bj = v[s][j];
    Base const bj1 = j+2 < len(v[s]) ? v[s][j+1] : Base::null();
    Base const bj0 = 1 < j ? v[s][j-1] : Base::null();
    
    auto x = acc.begin_stack(v, nick, s, j, model);
    // This gives the pfunc of every BP up to the new one, depending if the one after the new one stacks or not to the left
    auto y = acc.end_stack(v, nick, s, j, model);

    // No coaxial stacks
    for (auto bj_ : bases) {
        if (bj != bj_) o += sq(len(bases) + 1);
        else for_dangles(bases, [&](Base bj1_) {
            if (!dangle_match(bj1_, bj1)) o += len(bases) + 1;
            else for_dangles(bases, [&](Base bj0_) {
                if (dangle_match(bj0_, bj0)) {
                    real lo = bj0_.is_null() ? x(0) : x(1);
                    real hi = bj1_.is_null() ? (y(0) + (bj1.is_null() ? 0 : y(1))) // possibilities if no dangle to the right
                                                   : (y(0) + (j+3 == len(v[s]) ? 0 : y(1)));
                    o << lo * hi;
                } else ++o;
            });
        }); 
    }

    // Coaxial stacks on left, top and bottom
    for (auto bj_ : bases) {
        if (bj_ == bj) o << x(0) * (y(0) + (bj1.is_null() ? 0 : y(1))); // BP cannot stack left. BP after can stack left only if dangle exists.
        else ++o;
    }

    // Coaxial stacks on right, top and bottom
    for (auto bi : bases) {
        if (model.pairing().can_close(bi, bj) && len(v[s]) == 3 && !front(v[s]).is_null() && !back(v[s]).is_null()) {
            o << x(1) * y(1) * coaxial_factor(model, true, back(v[s ? s-1 : len(v)-1]), front(v[s]), bj, bi) 
                             * coaxial_factor(model, false, bi, bj, back(v[s]), front(v[(s+1) % len(v)])); // BP stack left, BP after must stack left
        } else ++o;
    }

    // Coaxial stack on top left only
    for (auto bj_ : bases) {
        if (bj != bj_) o += len(bases) + 1;
        else for_dangles(bases, [&](Base j1) {
            if (dangle_match(j1, bj1)) {
                o << x(0) * (y(0) + y(1) * (bj1.is_null() || (j+3 == len(v[s]) && !j1.is_null()) ? 0 : 1));
            } else ++o;
        });
    }

    // Coaxial stack on bottom left only
    for (auto bj_ : bases) {
        if (bj != bj_) o += len(bases) + 1;
        else for_dangles(bases, [&](Base j0) {
            if (dangle_match(j0, bj0)) {
                o << x(j0.is_null() ? 0 : 1)
                   * (y(0) + (bj1.is_null() ? 0 : y(1)));
            } else ++o;
        });
    }

    auto const dangle = [&model](Base id, Base i, Base j, Base jd) {
        return model.boltz(
            model.terminal_penalty(i, j) +
            (id.is_null() ? (jd.is_null() ? 0 : model.dangle5(i, j, jd)) : 
                           (jd.is_null() ? model.dangle3(id, i, j) : model.terminal_mismatch(id, i, j, jd)))
        );
    };

    // Coaxial stack on top right only
    if (j == 1 && !v[s][0].is_null()) {
        for (auto bi : bases)
            if (model.pairing().can_close(bj, bi)) {
                for_dangles(bases, [&](Base bi0_) {
                    o << x(1) * coaxial_factor(model, true, back(v[s ? s-1 : len(v)-1]), front(v[s]), bj, bi)
                    * (y(0) * (dangle(bi0_, bi, bj, Base::null()) + (bj1.is_null() ? 0 : dangle(bi0_, bi, bj, bj1))) 
                        + (bj1.is_null() ? 0 : (y(1) * (dangle(bi0_, bi, bj, Base::null()) + (j+3 == len(v[s]) ? 0 : dangle(bi0_, bi, bj, bj1))))));
                });
            } else o += len(bases) + 1;
    } else o += len(bases) * (len(bases) + 1);
            
    // Coaxial stack on bottom right only
    if (bj1.is_null() && !v[s][j+1].is_null()) {
        for (auto bi : bases)
            if (model.pairing().can_close(bj, bi)) {
                for_dangles(bases, [&](Base bi1_) {
                    o << (x(0) * dangle(Base::null(), bj, bi, bi1_) + (bj0.is_null() ? 0 : x(1) * dangle(bj0, bj, bi, bi1_))) 
                        * y(1) * coaxial_factor(model, false, bi, bj, back(v[s]), front(v[(s+1) % len(v)]));
                });
            } else o += len(bases) + 1;
    } else o += len(bases) * (len(bases) + 1);

    // Coaxial stack top left and bottom right
    for (auto const &[bi, bj_] : pairs)
        if (bj == bj_ && bj1.is_null() && !back(v[s]).is_null()) {
            o << x(0) * y(1) * coaxial_factor(model, false, bi, bj, back(v[s]), front(v[(s+1) % len(v)]));
        } else ++o;
    
    // Coaxial stack bottom left and top right
    for (auto const &[bi, bj_] : pairs)
        if (bj == bj_ && j == 1 && !front(v[s]).is_null()) {
            o << x(1) * coaxial_factor(model, true, back(v[s ? s-1 : len(v)-1]), front(v[s]), bj, bi) 
                 * (y(0) + (bj1.is_null() ? 0 : y(1)));
        } else ++o;

}

/******************************************************************************************/

}