package vcs

import (
    "fmt"
    "github.com/alinush/go-mcl"
    "github.com/wangnick2017/balanceproofs-go/asvc"
    "math/rand"
    "time"
)

func BasicUpdate_Fake() {
    fmt.Println("BalanceProofs Update all time:")
    for L := 20; L <= 30; L = L + 2 {
        fmt.Printf("L=%d\t", L)
        vc := VBAS{}
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
            id := rand.Uint64() % vc.N
            r[i] = asvc.UpdateReq{
                Index: id,
                Delta: delta,
            }
        }
        dt := time.Now()
        for i := 0; i < 1024; i++ {
            proofs, vec, aux = vc.UpdateAll_Fake(proofs, vec, r[i], aux)
        }
        duration := time.Since(dt)
        tt := duration.Seconds() / 1024
        fmt.Printf("Update time: %f sec \n", tt)
    }
}

func BasicQuery_Fake() {
    fmt.Println("BalanceProofs Query time:")
    for L := 20; L <= 30; L = L + 2 {
        fmt.Printf("L=%d\t", L)
        a := VBAS{}
        a.Init_Fake(uint8(L))
        n := 1 << L

        aux := a.InitAux()

        Lp := 1
        for ; Lp < L; Lp = Lp + Lp {
        }
        for i := 0; i < (1<<(L/2))*L/Lp; i++ {
            var delta mcl.Fr
            delta.Random()
            id := rand.Uint64() % uint64(n)
            r := asvc.UpdateReq{
                Index: id,
                Delta: delta,
            }
            aux = append(aux, r)
        }
        proofs := a.OpenAll_Fake([]mcl.Fr{})
        dt := time.Now()
        a.Query_Fake(uint64(0), proofs, aux)
        duration := time.Since(dt)
        fmt.Printf("Query time: %f sec\n", duration.Seconds())
    }
}
