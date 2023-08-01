#pragma once
// #include <nupack/thermo/Graph.h>
#include <nupack/common/Constants.h>
#include <nupack/types/Matrix.h>
#include <nupack/standard/Vec.h>
#include <nupack/thermo/BasicPF.h>
#include <nupack/thermo/Block.h>
#include <nupack/model/Model.h>
#include <nupack/reflect/Serialize.h>
#include <nupack/thermo/Engine.h>

namespace nupack::thermo {

template <class R, class T, class L>
struct MaterialTable : R {
    vec<RecursionMatrix<R, T, L>> table;
    NUPACK_REFLECT(MaterialTable, table);

    template <class Constants>
    void allocate(Constants const &c, std::size_t i, std::size_t j) {
        table.resize(c.single_material_models.size());
        for (auto &t : table) t.allocate(c, i, j);
    }

    template <class S, class Constants, class ...Ts>
    bool calculate(S, Index i, Index j, Constants const &c, Ts const &...args) {
        Material material_index = 0;
        for (auto &m : table) {
            if (!m.calculate(S(), i, j, c, args..., c.single_material_models[material_index], material_index)) return false;
            ++material_index;
        }
        return true;
    }

    template <class I, class J>
    decltype(auto) operator()(I const &i, J const &j, Material x) const {return table[x](i, j);} 

    void next_diagonal() const noexcept {}

    template <class U>
    void assign_and_clear(MaterialTable<R, U, L> &&u)  {table.resize(u.table.size()); zip(table, u.table, [](auto &x, auto &y) {x.assign_and_clear(std::move(y));});}

    static constexpr auto can_backtrack() {return False();}
};

}
