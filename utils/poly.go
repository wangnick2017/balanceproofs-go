package utils

import (
    "github.com/alinush/go-mcl"
    "github.com/wangnick2017/balanceproofs-go/fft"
    "math/bits"
)

func IsPolyZero(a []mcl.Fr) bool {

    n := len(a)
    if n == 0 {
        panic("IsPolyZero Error")
    }
    var flag bool
    flag = true
    for i := 0; i < n && flag; i++ {
        flag = flag && a[i].IsZero()
    }
    return flag
}

func PolyCondense(a []mcl.Fr) []mcl.Fr {
    n := len(a)
    if n == 0 {
        panic("PolyCondense Error")
    }

    i := n
    for i > 1 {
        if !a[i-1].IsZero() {
            break
        }
        i--
    }
    return a[:i]
}

func nextPowOf2(v uint64) uint64 {
    if v == 0 {
        return 1
    }
    return uint64(1) << bits.Len64(v-1)
}

func Max(x, y int) int {
    if x < y {
        return y
    }
    return x
}

func MulVecFr(a, b []mcl.Fr) []mcl.Fr {

    n := len(a)
    if n == len(b) && n > 0 {
        result := make([]mcl.Fr, n)
        for i := 0; i < n; i++ {
            mcl.FrMul(&result[i], &a[i], &b[i])
        }
        return result
    }
    result := make([]mcl.Fr, 0)
    return result
}

func PolyAdd(a []mcl.Fr, b []mcl.Fr) []mcl.Fr {

    if IsPolyZero(a) {
        return PolyCondense(b)
    }

    if IsPolyZero(b) {
        return PolyCondense(a)
    }

    aLen := len(a)
    bLen := len(b)
    n := Max(aLen, bLen)
    c := make([]mcl.Fr, n)

    for i := 0; i < n; i++ {
        if i < aLen {
            mcl.FrAdd(&c[i], &c[i], &a[i])
        }
        if i < bLen {
            mcl.FrAdd(&c[i], &c[i], &b[i])
        }
    }
    c = PolyCondense(c)
    return c
}

func PolyDiv(A []mcl.Fr, B []mcl.Fr) ([]mcl.Fr, []mcl.Fr) {
    if IsPolyZero(B) {
        var zero mcl.Fr
        zero.Clear()
        return A, []mcl.Fr{zero}
        //panic("PolyDiv: Cannot divide by zero polynomial.")
    }

    if len(B) > len(A) {
        panic("PolyDiv: Deg(B) should be <= Ded(A)")
    }

    a := make([]mcl.Fr, len(A))
    for i := 0; i < len(a); i++ {
        a[i] = A[i]
    }
    aPos := len(a) - 1
    bPos := len(B) - 1
    diff := aPos - bPos
    out := make([]mcl.Fr, diff+1)
    for diff >= 0 {
        quot := &out[diff]
        mcl.FrDiv(quot, &a[aPos], &B[bPos])
        var tmp, tmp2 mcl.Fr
        for i := bPos; i >= 0; i-- {
            // In steps: a[diff + i] -= b[i] * quot
            // tmp =  b[i] * quot
            mcl.FrMul(&tmp, quot, &B[i])
            // tmp2 = a[diff + i] - tmp
            mcl.FrSub(&tmp2, &a[diff+i], &tmp)
            // a[diff + i] = tmp2
            a[diff+i] = tmp2
        }
        aPos -= 1
        diff -= 1
    }
    a = PolyCondense(a)
    return out, a
}

func PolyMul(a []mcl.Fr, b []mcl.Fr) []mcl.Fr {
    if IsPolyZero(a) || IsPolyZero(b) {
        var zero mcl.Fr
        zero.Clear()
        return []mcl.Fr{zero}
    }

    aLen := len(a)
    bLen := len(b)
    if aLen == bLen && aLen == 1 {
        c := make([]mcl.Fr, 1)
        mcl.FrMul(&c[0], &a[0], &b[0])
        return c
    }
    n := uint64(2 * Max(aLen, bLen))
    n = nextPowOf2(n)

    var padding []mcl.Fr

    padding = make([]mcl.Fr, n-uint64(aLen))
    a = append(a, padding...)

    padding = make([]mcl.Fr, n-uint64(bLen))
    b = append(b, padding...)

    l := uint8(bits.Len64(n)) - 1 // n = 8 => 3 or 4?
    fs := fft.NewSettings(l)

    evalsA := fs.FFT_Fr(a, false)
    evalsB := fs.FFT_Fr(b, false)

    res := fs.FFT_Fr(MulVecFr(evalsA, evalsB), true)
    res = res[:(aLen + bLen - 1)]
    //res = PolyCondense(res)
    return res
}

type SubProduct struct {
    Poly  []mcl.Fr
    left  *SubProduct
    right *SubProduct
}

func SubProductTree(index []mcl.Fr) *SubProduct {
    n := uint64(len(index))
    if n == 1 {
        ret := make([]mcl.Fr, 2)
        mcl.FrNeg(&ret[0], &index[0])
        ret[1].SetInt64(1)
        return &SubProduct{
            Poly:  ret,
            left:  nil,
            right: nil,
        }
    }
    n2 := uint64(n / 2)
    ls := SubProductTree(index[:n2])
    rs := SubProductTree(index[n2:])
    //fmt.Printf("n=%d leftpoly=%d ri=%d\n", n, len(ls.poly), len(rs.poly))
    return &SubProduct{
        Poly:  PolyMul(ls.Poly, rs.Poly),
        left:  ls,
        right: rs,
    }
}

func PolyInterpolation(index []mcl.Fr, eval_with_prime []mcl.Fr, sp *SubProduct) []mcl.Fr {
    n := uint64(len(index))
    if n == 1 {
        return []mcl.Fr{eval_with_prime[0]}
    }
    n2 := n / 2
    l := PolyInterpolation(index[:n2], eval_with_prime[:n2], sp.left)
    r := PolyInterpolation(index[n2:], eval_with_prime[n2:], sp.right)
    return PolyAdd(PolyMul(l, sp.right.Poly), PolyMul(r, sp.left.Poly))
}

func PolyMultiEvaluate(f []mcl.Fr, sp *SubProduct) []mcl.Fr {
    n := uint64(len(f))
    if n == 0 {
        return f
    }
    if n == 1 {
        l := len(sp.Poly) - 1
        ret := make([]mcl.Fr, l)
        for i := 0; i < l; i++ {
            ret[i] = f[0]
        }
        return ret
    }
    if len(sp.Poly) <= 2 {
        ret := make([]mcl.Fr, 1)
        mcl.FrEvaluatePolynomial(&ret[0], f, &(sp.Poly[0]))
        return ret
    }
    _, aL := PolyDiv(f, sp.left.Poly)
    _, aR := PolyDiv(f, sp.right.Poly)

    l := PolyMultiEvaluate(aL, sp.left)
    r := PolyMultiEvaluate(aR, sp.right)

    return append(l, r...)
}

func PolyDifferentiate(a []mcl.Fr) []mcl.Fr {
    n := uint64(len(a))
    if n == 0 {
        panic("PolyDifferentiate: Input is empty")
    }
    if n == 1 {
        ret := make([]mcl.Fr, 1)
        ret[0].Clear()
        return ret
    }
    c := make([]mcl.Fr, n)
    var temp mcl.Fr
    for i := uint64(1); i < n; i++ {
        temp.SetInt64(int64(i))
        mcl.FrMul(&c[i], &a[i], &temp)
    }
    return c[1:]
}
