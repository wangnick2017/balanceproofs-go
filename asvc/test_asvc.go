package asvc

import (
    "fmt"
    "github.com/alinush/go-mcl"
    "github.com/wangnick2017/balanceproofs-go/utils"
    "time"
)

func AsvcCommitOpenAll() {
    fmt.Println("BalanceProofs/asvc Commit time and OpenAll time:")
    for L := 20; L <= 30; L = L + 2 {
        fmt.Printf("L=%d\t", L)
        vc := ASVC{}

        l := uint8(L)
        n := 1 << l
        vc.Init(l, []mcl.G1{}, []mcl.G2{})

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

func AsvcAgg_Fake() {
    fmt.Println("BalanceProofs/asvc Agg and Verify agg time:")
    for L := 20; L <= 30; L = L + 2 {
        fmt.Printf("L=%d\t", L)
        a := ASVC{}

        l := uint8(L)
        a.Init_Fake(l)
        vec := make([]mcl.Fr, 2)
        vec[0].Random()
        vec[1].Random()
        var proof mcl.G1
        proof.Random()
        aggs := make([]Inp, 1024)
        aggvs := make([]Val, 1024)
        lis := utils.GenerateIndices(uint64(1<<a.L), 1024)
        for j := 0; j < 1024; j++ {
            id := lis[j]
            aggs[j] = Inp{Index: id, Proof: proof}
            aggvs[j] = Val{Index: id, Y: vec[id&1]}
        }
        dt := time.Now()
        a.Aggregate_Fake(aggs)
        duration := time.Since(dt)
        var digest mcl.G1
        digest.Random()
        fmt.Printf("Agg time: %f sec\t", duration.Seconds())

        dt = time.Now()
        a.VerifyAggregation_Fake(digest, proof, aggvs)
        duration = time.Since(dt)
        fmt.Printf("Verify agg time: %f sec\n", duration.Seconds())
    }
}

func AsvcIndvidual_Fake() {
    fmt.Println("BalanceProofs/asvc Verify individual time:")
    for L := 20; L <= 30; L = L + 2 {
        fmt.Printf("L=%d\t", L)
        vc := ASVC{}
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
            vc.VerifySingle_Fake(digest, Pi, Val{
                Index: 0,
                Y:     v,
            })
            duration := time.Since(dt)
            t = t + duration.Seconds()
        }

        fmt.Printf("Verify ind time: %f sec\n", t/1024.0)
    }
}
