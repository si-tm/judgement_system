#pragma once
#include <taskflow/taskflow.hpp>
#include <nupack/reflect/Reflection.h>
#include <nupack/reflect/Print.h>
#include <nupack/common/Error.h>

namespace nupack {

template <>
struct io::Printer<tf::Task> {
    void operator()(std::ostream &os, tf::Task const &t, Indent) const {
        if (t.empty()) os << "Task()";
        else os << "Task(\"" << t.name() << "\")";
    }
};

/******************************************************************************************/

class SharedError {
    struct Impl {
        std::atomic<bool> set = false; // set before ptr is stored
        std::atomic<bool> stored = false; // set after ptr is stored
        std::exception_ptr exception;
    };
    std::shared_ptr<Impl> ptr = std::make_shared<Impl>();

public:

    NUPACK_REFLECT(SharedError, ptr);

    void rethrow_if_set() const;
    void set_to_current_exception() noexcept;
    bool is_set() const noexcept;
    bool clear() noexcept;

    template <class Exc>
    void set_exception(Exc exc) noexcept {
        try {throw std::move(exc);}
        catch (Exc const &) {set_to_current_exception();}
    }

    template <class F, class ...Args>
    void invoke_noexcept(F &&f, Args &&...args) noexcept {
        try {
            f(fw<Args>(args)...);
        } catch (...) {
            set_to_current_exception();
        }
    }
};

/******************************************************************************************/

struct SharedExecutor {
    std::shared_ptr<tf::Executor> impl;

    SharedExecutor() = default;

    explicit SharedExecutor(uint n) 
        : impl(n ? std::make_shared<tf::Executor>(n) : std::make_shared<tf::Executor>()) {
            NUPACK_REQUIRE(n, <=, 10000, "maximum threads is 10000 for safety in locking mechanisms");
        }

    explicit operator bool() const {return bool(impl);}

    auto run(tf::Taskflow &flow) {NUPACK_ASSERT(bool(impl)); return impl->run(flow);}

    NUPACK_REFLECT(SharedExecutor, impl);
};

/******************************************************************************************/

}