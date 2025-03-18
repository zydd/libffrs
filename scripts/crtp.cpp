#include <cstdio>
#include <tuple>


template<typename...Fs>
class CRTP : public Fs... {
public:
    template<typename F>
    static constexpr CRTP const& cast(F _this) {
        return static_cast<CRTP const&>(*_this);
    }

    template<typename F>
    static constexpr CRTP& cast_mut(F _this) {
        return static_cast<CRTP&>(*_this);
    }

    template<typename...F0>
    inline CRTP(F0&&... f0):
        F0(std::move(f0))...
    { }
};

template<typename T0>
struct s_1 {
    template<typename T, typename C>
    struct types1 {
        using S1 = T0;
        int s1=666;
        types1(int a): s1(a) { printf("s_1()\n"); };
    };
};

template<typename T, typename C>
struct s_2 { int s2; s_2(): s2(26) { printf("s_2()\n"); }; };

template<typename T, typename C>
struct s_3 {
    int s2; s_3(): s2(39) {
        printf("s_3(), s1:%d\n", T(C::cast(this).s1));
    };

    void m_3(T t) {
        printf("m_3(%d)\n", t);
    }
};

template<typename T, template<class, class>typename...Fs>
class GF : public CRTP<Fs<T, GF<T, Fs...>>...>  {
public:
    using GFT = T;
    using CRTP<Fs<T, GF>...>::CRTP;
};

using typ = GF<char, s_1<char>::types1, s_2, s_3>;

int main() {
    typ a = typ(typ::types1(256 + 72));

    printf("%d %d %d\n", a.s1, a.s_2::s2, a.s_3::s2);
    a.m_3(33);

    return 0;
}