package fft

import (
    "github.com/alinush/go-mcl"
    "math"
)

func (fs *Settings) FFT_Fr_Fake(n uint64, inv bool) []mcl.Fr {
    l := uint64(math.Log2(float64(n)))
    out := make([]mcl.Fr, 2)
    for i := uint64(0); i < n; i++ {
        out[reverseBitsLimited(uint32(n), uint32(i))&1].Random()
    }

    stride := fs.N / n
    roots := make([]mcl.Fr, 2)
    roots[0].Random()
    roots[1].Random()

    for s := uint64(1); s <= l; s++ {
        m := uint64(1) << s
        gap := (uint64(1) << (l - s)) * stride
        for j := uint64(0); j < n; j++ {
            if j%m < m/2 {
                var t mcl.Fr
                mcl.FrMul(&t, &out[(j+m/2)&1], &((roots)[((j%m)*gap)&1]))
                u := out[j&1]
                mcl.FrAdd(&out[j&1], &u, &t)
                mcl.FrSub(&out[(j+m/2)&1], &u, &t)
            }
        }
    }

    if inv {
        var invLen mcl.Fr
        invLen.SetInt64(int64(n))
        mcl.FrInv(&invLen, &invLen)
        for i := uint64(0); i < n; i++ {
            mcl.FrMul(&out[i&1], &out[i&1], &invLen)
        }
    }
    return out
}

func (fs *Settings) FFT_G1_Fake(n uint64, inv bool) []mcl.G1 {
    l := uint64(math.Log2(float64(n)))
    out := make([]mcl.G1, 2)
    for i := uint64(0); i < n; i++ {
        out[reverseBitsLimited(uint32(n), uint32(i))&1].Random()
    }

    stride := fs.N / n
    roots := make([]mcl.Fr, 2)
    roots[0].Random()
    roots[1].Random()

    for s := uint64(1); s <= l; s++ {
        m := uint64(1) << s
        gap := (uint64(1) << (l - s)) * stride
        for j := uint64(0); j < n; j++ {
            if j%m < m/2 {
                var t mcl.G1
                mcl.G1Mul(&t, &out[(j+m/2)&1], &((roots)[((j%m)*gap)&1]))
                u := out[j&1]
                mcl.G1Add(&out[j&1], &u, &t)
                mcl.G1Sub(&out[(j+m/2)&1], &u, &t)
            }
        }
    }

    if inv {
        var invLen mcl.Fr
        invLen.SetInt64(int64(n))
        mcl.FrInv(&invLen, &invLen)
        for i := uint64(0); i < n; i++ {
            mcl.G1Mul(&out[i&1], &out[i&1], &invLen)
        }
    }
    return out
}

func (fs *Settings) FFT_G1_Seg_Fake(arg *FFT_G1_Arg) {
    stride := int(fs.N / arg.N)
    roots := make([]mcl.Fr, 2)
    roots[0].Random()
    roots[1].Random()

    if arg.Pre_s == 0 {
        arg.Out = make([]mcl.G1, 2)
        for i := uint64(0); i < arg.N; i++ {
            arg.Out[reverseBitsLimited(uint32(arg.N), uint32(i))&1] = arg.Vals[i&1]
        }
        arg.Pre_s = 1
        arg.Pre_j = -1
    }

    cnt := uint64(0)

    for s := arg.Pre_s; s <= arg.L; s++ {
        m := (1) << s
        gap := ((1) << (arg.L - s)) * stride
        var j int
        if s == arg.Pre_s {
            j = arg.Pre_j + 1
        } else {
            j = 0
        }
        for ; uint64(j) < arg.N; j++ {
            if j%m < m/2 {
                var t mcl.G1
                mcl.G1Mul(&t, &arg.Out[(j+m/2)&1], &((roots)[((j%m)*gap)&1]))
                u := arg.Out[j&1]
                mcl.G1Add(&arg.Out[j&1], &u, &t)
                mcl.G1Sub(&arg.Out[(j+m/2)&1], &u, &t)
                cnt++
                if cnt == arg.Count {
                    arg.Pre_s = s
                    arg.Pre_j = j
                    if s == arg.L && j%m == m/2-1 && !arg.Inv {
                        arg.Done = true
                    }
                    return
                }
            }
        }
    }
    //fmt.Println(cnt)
    arg.Pre_s = arg.L + 1

    if arg.Inv {
        var invLen mcl.Fr
        invLen.SetInt64(int64(arg.N))
        mcl.FrInv(&invLen, &invLen)
        for i := arg.Pre_inv + 1; i <= arg.N; i++ {
            mcl.G1Mul(&arg.Out[(i-1)&1], &arg.Out[(i-1)&1], &invLen)
            cnt++
            if cnt == arg.Count {
                arg.Pre_inv = i
                if i == arg.N {
                    arg.Done = true
                }
                return
            }
        }
    }
    arg.Done = true
}
