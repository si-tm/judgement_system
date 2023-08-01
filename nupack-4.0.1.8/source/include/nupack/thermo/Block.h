#pragma once
#include "../model/ModelVariants.h"
#include "BasicPF.h"
#include "CoaxialPF.h"
#include "Tensor.h"
#include "FastInterior.h"

namespace nupack::thermo {

/******************************************************************************************/

template <class T>
struct RegionBonus : E_t {
    native_type<T> one;

    template <class Constants>
    void allocate(Constants const &c, std::size_t i, std::size_t j) {one = c.one();}

    void next_diagonal() const noexcept {}

    bool calculate(Ignore, Ignore, Ignore, Ignore) const noexcept {return true;}

    auto operator()(Index i, Index j) const noexcept {return one;}

    template <class U, class Converter>
    void assign_and_clear(RegionBonus<U> const &r, Converter const &c) {c(r.one, one);}

    static auto can_backtrack() {return False();}

    void clear() {}
};

/******************************************************************************************/

template <class T, class Algorithm>
struct Block;

/******************************************************************************************/

template <class T>
struct Block<T, NoStacking> {
    using value_type = T;
    using StrandData = DangleData<native_type<T>>;

    RegionBonus<T>    E;
    FastInterior<T> X;
    RecursionMatrix<BasicN,  T,  SpanSecond> N;
    RecursionMatrix<BasicM2, T,  SpanSecond> M2;
    RecursionMatrix<BasicM3, T,  SpanSecond> M3;
    RecursionMatrix<BasicB,  T,  SpanFirst>  B;
    RecursionMatrix<BasicZ,  T,  SpanBoth>   Z;
    RecursionMatrix<BasicD,  T,  SpanSecond> D;
    RecursionMatrix<BasicYA, T,  SpanBoth>   YA;
    RecursionMatrix<BasicYB, T,  SpanBoth>   YB;
    RecursionMatrix<BasicMS, T,  SpanFirst>  MS;
    RecursionMatrix<BasicM1, T,  SpanSecond> M1;
    RecursionMatrix<BasicM,  T,  SpanSecond> M;
    RecursionMatrix<BasicS,  T,  SpanFirst>  S;
    RecursionMatrix<BasicQ,  T,  SpanSecond> Q;

    NUPACK_REFLECT(Block, E, X, N, M2, M3, B, Z, D, YA, YB, MS, M1, M, S, Q);
};

/******************************************************************************************/

template <class T>
struct AllOrMinDangleBlock {
    using value_type = T;
    using StrandData = DangleData<native_type<T>>;

    RegionBonus<T>    E;
    FastInterior<T> X;
    RecursionMatrix<BasicN,   T, SpanSecond> N;
    RecursionMatrix<BasicM2,  T, SpanSecond> M2;
    RecursionMatrix<BasicB,   T, SpanFirst>  B;
    RecursionMatrix<BasicZ,   T, SpanBoth>   Z;
    RecursionMatrix<BasicD,   T, SpanSecond> D;
    RecursionMatrix<BasicYA,  T, SpanBoth>   YA;
    RecursionMatrix<BasicYB,  T, SpanBoth>   YB;
    RecursionMatrix<DangleMS, T, SpanFirst>  MS;
    RecursionMatrix<DangleM,  T, SpanSecond> M;
    RecursionMatrix<DangleS,  T, SpanFirst>  S;
    RecursionMatrix<BasicQ,   T, SpanSecond> Q;

    NUPACK_REFLECT(AllOrMinDangleBlock, E, X, N, M2, B, Z, D, YA, YB, MS, M, S, Q);
};

template <class T>
struct Block<T, AllDangles> : AllOrMinDangleBlock<T> {};

template <class T>
struct Block<T, MinDangles> : AllOrMinDangleBlock<T> {};
    
/******************************************************************************************/

template <class R, class T, class L>
struct Table : R {
    vec<RecursionMatrix<R, T, L>> table;
    NUPACK_REFLECT(Table, table);

    template <class Constants>
    void allocate(Constants const &c, std::size_t i, std::size_t j) {
        table.resize(c.model.energy_model.alphabet().length());
        for (auto &t : table) t.allocate(c, i, j);
    }

    template <class S, class Constants, class ...Ts>
    bool calculate(S, Index i, Index j, Constants const &c, Ts const &...args) {
        BaseIndex x = 0;
        for (auto &m : table)
            if (!m.calculate(S(), i, j, c, args..., Base(x++))) return false;
        return true;
    }

    template <class I, class J>
    decltype(auto) operator()(I const &i, J const &j, Base x) const {return table[+x](i, j);} 

    void next_diagonal() const noexcept {}

    template <class U, class F>
    void assign_and_clear(Table<R, U, L> &&u, F &&f)  {table.resize(u.table.size()); zip(table, u.table, [&](auto &x, auto &y) {x.assign_and_clear(std::move(y), f);});}

    static constexpr auto can_backtrack() {return False();}
};

/******************************************************************************************/

template <class R, class T, class L>
struct Table2 : Table<R, T, L> {
    BaseIndex n;

    template <class Constants>
    void allocate(Constants const &c, std::size_t i, std::size_t j) {
        n = c.model.energy_model.alphabet().length();
        this->table.resize(n * n);
        for (auto &t : this->table) t.allocate(c, i, j);
    }

    template <class S, class Constants>
    bool calculate(S, Index i, Index j, Constants const &c) {
        auto m = this->table.begin();
        for (auto x : range(n)) for (auto y : range(n)) {
            if (!m->calculate(S(), i, j, c, Base(x), Base(y))) return false;
            ++m;
        }
        return true;
    }

    template <class I, class J>
    decltype(auto) operator()(I const &i, J const &j, Base x, Base y) const {return this->table[+x * n + +y](i, j);} 
};

/******************************************************************************************/

template <class T>
struct Toggle {
    std::pair<T, T> ts;

    NUPACK_REFLECT(Toggle, ts);

    template <class Constants>
    void allocate(Constants const &c, std::size_t i, std::size_t j) {
        for_each(ts, [&](auto &m) {m.allocate(c, i, j);});
    }

    template <class S, class Constants, class ...Ts>
    bool calculate(S, Index i, Index j, Constants const &c, Ts const &...args) {
        return std::get<0>(ts).calculate(S(), i, j, c, args..., d0)
            && std::get<1>(ts).calculate(S(), i, j, c, args..., d1);
    }

    template <class I, class J, class N, class ...Ts>
    decltype(auto) operator()(I const &i, J const &j, N, Ts &&...ts) const {return std::get<N::value>(this->ts)(i, j, std::forward<Ts>(ts)...);} 
    
    template <class I, class J>
    decltype(auto) operator()(I const &i, J const &j) const {return std::get<0>(this->ts)(i, j);} 

    decltype(auto) corner(bool i, bool j) const {return std::get<0>(ts).corner(i, j);}

    void next_diagonal() const noexcept {}
    static constexpr auto can_backtrack() {return False();}

    template <class ...Ts>
    decltype(auto) addressed(Ts &&...args) const {return std::get<0>(ts).addressed(std::forward<Ts>(args)...);}
};

/******************************************************************************************/

template <class R, class T, class L>
struct Toggle2 : R, Toggle<Toggle<RecursionMatrix<R, T, L>>> {};

template <class R, class T, class L>
struct Toggle1 : R, Toggle<RecursionMatrix<R, T, L>> {};

template <class R, class T, class L>
struct ToggleTable : R, Toggle<Table<R, T, L>> {};

/******************************************************************************************/

template <class T>
struct Block<T, Stacking> : MemberComparable {
    using value_type = T;
    using StrandData = CoaxialRows<native_type<T>>;

    RegionBonus<T>    E;
    FastInterior<T> X;
    RecursionMatrix<CoaxialN,     T, SpanBoth>    N;
    RecursionMatrix<CoaxialM2,  T, SpanSecond>  M2;
    RecursionMatrix<CoaxialB,   T, SpanBoth>    B;
    RecursionMatrix<BasicZ,     T, SpanBoth>    Z;
    RecursionMatrix<BasicD,     T, SpanBoth>    D;
    RecursionMatrix<BasicYA,    T, SpanBoth>    YA;
    RecursionMatrix<BasicYB,    T, SpanBoth>    YB;
    RecursionMatrix<CoaxialMD,  T, SpanSecond>  MD;
    RecursionMatrix<CoaxialMC,  T, SpanSecond>  MC;
    RecursionMatrix<CoaxialMCS, T, SpanFirst>   MCS;
    RecursionMatrix<CoaxialMS,  T, SpanFirst>   MS; // MS == CMS + MCS; CMS never referenced alone
    RecursionMatrix<CoaxialCD,  T, SpanSecond>  CD;
    RecursionMatrix<CoaxialS,   T, SpanFirst>   S;
    RecursionMatrix<CoaxialM,   T, SpanBoth>    M;
    RecursionMatrix<CoaxialQ,   T, SpanBoth>    Q;
    NUPACK_REFLECT(Block, E, X, N, M2, B, Z, D, YA, YB, MD, MC, MCS, MS, CD, S, M, Q);
};

/******************************************************************************************/

template <class Block, class Constants>
void allocate_block(Block &b, Constants const &c) {
    for_each(members_of(b), [&](auto &m) {m.allocate(c, c.lsize(), c.rsize());});
}



}
