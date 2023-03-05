package vcs

import (
    "fmt"
    "github.com/alinush/go-mcl"
    "github.com/wangnick2017/balanceproofs-go/asvc"
    "github.com/wangnick2017/balanceproofs-go/utils"
    "math/rand"
    "time"
)

func BucketCommitOpenAll() {
    fmt.Println("Bucketing Commit time and OpenAll time:")
    for L := 20; L <= 30; L = L + 2 {
        fmt.Printf("L=%d\t", L)
        vc := VBUC{}

        l := uint8(L)
        n := 1 << l
        vc.Init(l)

        vec := make([]mcl.Fr, n, n)
        for i := 0; i < n; i++ {
            vec[i].Random()
        }
        dt := time.Now()
        vc.Commit(vec)
        duration := time.Since(dt)
        fmt.Printf("Commit time: %f sec\t", duration.Seconds())

        dt = time.Now()
        vc.OpenAll(vec)
        duration = time.Since(dt)
        fmt.Printf("OpenAll time: %f sec\n", duration.Seconds())
    }

}

func GetXYZ(i uint64, n1 uint64, n2 uint64, n3 uint64) (x uint64, y uint64, z uint64) {
    x = i / (n2 * n3)
    y = (i % (n2 * n3)) / n3
    z = i % n3
    return x, y, z
}

func BucketProofSize() {
    fmt.Println("Bucketing Proof size:")
    m := 1024
    for L := 20; L <= 30; L = L + 2 {
        lis := utils.GenerateIndices(uint64(1<<L), m)
        cnt := 1
        c := 1
        x, y, _ := GetXYZ(lis[0], uint64(1<<(L/4)), uint64(1<<L/4), uint64(1<<(L-L/4-L/4)))
        for i := 1; i < m; i++ {
            x_, y_, _ := GetXYZ(lis[i], uint64(1<<(L/4)), uint64(1<<(L/4)), uint64(1<<(L-L/4-L/4)))
            if x != x_ || y != y_ {
                cnt++
            }
            if x != x_ {
                c++
            }
            x, y = x_, y_
        }
        fmt.Printf("L=%d\t", L)
        fmt.Printf("Halving: %f KB\t", float64(cnt+c+1)/1024.0*48.0)
        fmt.Printf("Without halving: %f KB\n", float64(cnt+cnt+c+1)/1024.0*48.0)
    }
}

func BucketIndvidual_Fake() {
    fmt.Println("Bucketing Verify individual time:")
    for L := 20; L <= 30; L = L + 2 {
        fmt.Printf("L=%d\t", L)
        vc := VBUC{}
        vc.Init_Fake(uint8(L))
        t := 0.0
        for i := 0; i < 1024; i++ {
            var digest mcl.G1
            digest.Random()
            var v mcl.Fr
            v.Random()
            var Pii, Psi, Pi mcl.G1
            Pii.Random()
            Psi.Random()
            Pi.Random()
            dt := time.Now()
            vc.VerifySingle_Fake(digest, BucProofSingle{
                Pii: Pii,
                Psi: Psi,
                Pi:  Pi,
            }, asvc.Val{
                Index: 0,
                Y:     v,
            })
            duration := time.Since(dt)
            t = t + duration.Seconds()
        }

        fmt.Printf("Verify ind time: %f sec\n", t/1024.0)
    }
}

func BucketAgg_Fake() {
    fmt.Println("Bucketing Agg and Verify agg time:")
    m := 1024
    for L := 20; L <= 30; L = L + 2 {
        fmt.Printf("L=%d\t", L)
        vc := VBUC{}
        vc.Init_Fake(uint8(L))
        lis := utils.GenerateIndices(uint64(1<<L), m)
        lp := make([]BucProofSingle, m)
        for i := 0; i < m; i++ {
            lp[i] = BucProofSingle{}
            lp[i].Pii.Random()
            lp[i].Psi.Random()
            lp[i].Pii.Random()
        }
        dt := time.Now()
        agg := vc.Aggregate_Fake(lis, lp)
        duration := time.Since(dt)
        fmt.Printf("Agg time: %f sec\t", duration.Seconds())
        vals := make([][][]mcl.Fr, len(agg.X))
        for i := 0; i < len(agg.X); i++ {
            vals[i] = make([][]mcl.Fr, len(agg.Y[i]))
            for j := 0; j < len(agg.Y[i]); j++ {
                vals[i][j] = make([]mcl.Fr, len(agg.Z[i][j]))
                for k := 0; k < len(agg.Z[i][j]); k++ {
                    vals[i][j][k].Random()
                }
            }
        }
        var digest mcl.G1
        digest.Random()
        dt = time.Now()
        vc.VerifyAggregation_Fake(digest, agg, vals)
        duration = time.Since(dt)
        fmt.Printf("Verify agg time: %f sec\n", duration.Seconds())
    }
}

func BucketQuery_Fake() {
    fmt.Println("Bucketing Query time:")
    for L := 20; L <= 30; L = L + 2 {
        fmt.Printf("L=%d\t", L)
        vc := VBUC{}
        vc.Init_Fake(uint8(L))
        vec := make([]mcl.Fr, 2)
        vec[0].Random()
        vec[1].Random()
        proofs := vc.OpenAll_Fake(vec)
        aux := vc.InitAux()
        for i := 0; i < (1<<(vc.L3/2))*(int(vc.L3))/2; i++ {
            var delta mcl.Fr
            delta.Random()
            id := rand.Uint64() % vc.N3
            r := asvc.UpdateReq{
                Index: id,
                Delta: delta,
            }
            aux[0][0] = append(aux[0][0], r)
        }
        dt := time.Now()
        vc.Query_Fake(0, proofs, aux)
        duration := time.Since(dt)
        fmt.Printf("Query time: %f sec\n", duration.Seconds())
    }
}

func BucketUpdate_Fake() {
    fmt.Println("Bucketing Update all time:")
    for L := 20; L <= 30; L = L + 2 {
        fmt.Printf("L=%d\t", L)
        vc := VBUC{}
        vc.Init_Fake(uint8(L))
        vec := make([]mcl.Fr, 2)
        vec[0].Random()
        vec[1].Random()
        proofs := vc.OpenAll_Fake(vec)
        aux := vc.InitAux()
        r := make([]asvc.UpdateReq, 1024)
        for i := 0; i < 1024; i++ {
            var delta mcl.Fr
            delta.Random()
            id := rand.Uint64() % vc.N3
            r[i] = asvc.UpdateReq{
                Index: id,
                Delta: delta,
            }
        }
        dt := time.Now()
        for i := 0; i < 1024; i++ {
            vec, aux = vc.UpdateAll_Fake(&proofs, vec, r[i], aux)
        }
        duration := time.Since(dt)
        tt := duration.Seconds() / 1024
        fmt.Printf("Update time: %f sec \n", tt)
    }
}
