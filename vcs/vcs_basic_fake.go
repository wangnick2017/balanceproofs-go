package vcs

import (
    "github.com/alinush/go-mcl"
    "github.com/wangnick2017/balanceproofs-go/asvc"
)

func (vbas *VBAS) Init_Fake(L uint8) {
    vbas.L = L
    vbas.N = 1 << L
    vbas.asvc = &asvc.ASVC{}
    vbas.asvc.Init_Fake(L)
    vbas.args = nil
}

func (vbas *VBAS) Commit_Fake(vector []mcl.Fr) mcl.G1 {
    return vbas.asvc.Commit_Fake(vector)
}

func (vbas *VBAS) OpenAll_Fake(vector []mcl.Fr) []mcl.G1 {
    return vbas.asvc.OpenAll_Fake(vector)
}

func (vbas *VBAS) updateVector_Fake(vector []mcl.Fr, list []asvc.UpdateReq) []mcl.Fr {
    for i := 0; i < len(list); i++ {
        temp := vector[list[i].Index&1]
        mcl.FrAdd(&vector[list[i].Index&1], &temp, &list[i].Delta)
    }
    return vector
}

func (vbas *VBAS) Query_Fake(index uint64, proofs []mcl.G1, aux []asvc.UpdateReq) mcl.G1 {
    var roots mcl.Fr
    roots.Random()

    g := make([]mcl.G1, len(aux))
    fr := make([]mcl.Fr, len(aux))
    for i := 0; i < len(aux); i++ {
        if index == aux[i].Index {
            g[i] = vbas.asvc.U_i[index&1]
            fr[i] = aux[i].Delta
        } else {
            var omega_i mcl.Fr
            mcl.FrAdd(&omega_i, &roots, &roots)

            var w_i_j mcl.G1
            mcl.G1Sub(&w_i_j, &vbas.asvc.A_i[aux[i].Index&1], &vbas.asvc.A_i[index&1])

            var n, ao mcl.Fr
            n.SetInt64(int64(vbas.asvc.N))
            mcl.FrDiv(&ao, &roots, &n)
            mcl.FrMul(&ao, &ao, &aux[i].Delta)
            mcl.FrDiv(&omega_i, &ao, &omega_i)

            g[i] = w_i_j
            fr[i] = omega_i
        }
    }
    var temp mcl.G1
    mcl.G1MulVec(&temp, g, fr)
    mcl.G1Add(&temp, &temp, &proofs[index&1])
    return temp
}

func (vbas *VBAS) UpdateCommitment_Fake(digest mcl.G1, req asvc.UpdateReq) mcl.G1 {
    return vbas.asvc.UpdateCommitment_Fake(digest, req)
}

func (vbas *VBAS) UpdateAll_Fake(proofs []mcl.G1, vector []mcl.Fr, req asvc.UpdateReq, aux []asvc.UpdateReq) ([]mcl.G1, []mcl.Fr, []asvc.UpdateReq) {
    aux = append(aux, req)
    l := uint64(len(aux))
    if l*l < vbas.asvc.N && l < 1024/2 {
        if vbas.args == nil {
            return proofs, vector, aux
        } else {
            panic("error: update all")
        }
    }

    if vbas.args == nil || vbas.args.Done {
        vector = vbas.updateVector_Fake(vector, aux)
        vbas.P = len(aux)
        Lp := uint8(1)
        for ; Lp < vbas.asvc.L; Lp = Lp + Lp {
        }
        vbas.args = &asvc.OpenAllStepArg{
            Vector: vector,
            Step:   0,
            Proofs: proofs,
            Count:  (1 << (vbas.asvc.L / 2)) * uint64(Lp) / 2,
            Done:   false,
        }
        vbas.asvc.OpenAllStep_Fake(vbas.args)
        return proofs, vector, aux
    }
    vbas.args.Step++
    vbas.asvc.OpenAllStep_Fake(vbas.args)
    if vbas.args.Done {
        aux = aux[vbas.P:]
    }
    return vbas.args.Proofs, vbas.args.Vector, aux
}

func (vbas *VBAS) UpdateProof_Fake(proof mcl.G1, index uint64, req asvc.UpdateReq) mcl.G1 {
    return vbas.asvc.UpdateProof_Fake(proof, index, req)
}

func (vbas *VBAS) Aggregate_Fake(aggs []asvc.Inp) mcl.G1 {
    return vbas.asvc.Aggregate_Fake(aggs)
}

func (vbas *VBAS) VerifySingle_Fake(digest mcl.G1, proof mcl.G1, v asvc.Val) bool {
    return vbas.asvc.VerifySingle_Fake(digest, proof, v)
}

func (vbas *VBAS) VerifyAggregation_Fake(digest mcl.G1, proof mcl.G1, aggvs []asvc.Val) bool {
    return vbas.asvc.VerifyAggregation_Fake(digest, proof, aggvs)
}
