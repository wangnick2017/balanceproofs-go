package vcs

import (
    "github.com/alinush/go-mcl"
    "github.com/wangnick2017/balanceproofs-go/asvc"
)

type VBAS struct {
    N uint64
    L uint8
    P int

    args *asvc.OpenAllStepArg

    asvc *asvc.ASVC
}

func (vbas *VBAS) Init(L uint8) {
    vbas.L = L
    vbas.N = 1 << L
    vbas.asvc = &asvc.ASVC{}
    vbas.asvc.Init(L, []mcl.G1{}, []mcl.G2{})
    vbas.args = nil
}

func (vbas *VBAS) InitAux() []asvc.UpdateReq {
    return make([]asvc.UpdateReq, 0)
}

func (vbas *VBAS) Commit(vector []mcl.Fr) mcl.G1 {
    return vbas.asvc.Commit(vector)
}

func (vbas *VBAS) Open(index uint64, vector []mcl.Fr, aux []asvc.UpdateReq) mcl.G1 {
    temp := vbas.asvc.Open(index, vector)
    for i := 0; i < len(aux); i++ {
        temp = vbas.UpdateProof(temp, index, aux[i])
    }
    return temp
}

func (vbas *VBAS) OpenAll(vector []mcl.Fr) []mcl.G1 {
    return vbas.asvc.OpenAll(vector)
}

func (vbas *VBAS) updateVector(vector []mcl.Fr, list []asvc.UpdateReq) []mcl.Fr {
    for i := 0; i < len(list); i++ {
        temp := vector[list[i].Index]
        mcl.FrAdd(&vector[list[i].Index], &temp, &list[i].Delta)
    }
    return vector
}

func (vbas *VBAS) Query(index uint64, proofs []mcl.G1, aux []asvc.UpdateReq) mcl.G1 {
    g := make([]mcl.G1, len(aux))
    fr := make([]mcl.Fr, len(aux))
    for i := 0; i < len(aux); i++ {
        if index == aux[i].Index {
            g[i] = vbas.asvc.U_i[index]
            fr[i] = aux[i].Delta
        } else {
            var omega_i mcl.Fr
            mcl.FrSub(&omega_i, &vbas.asvc.FFT.RootsOfUnity[aux[i].Index*2], &vbas.asvc.FFT.RootsOfUnity[index*2])

            var w_i_j mcl.G1
            mcl.G1Sub(&w_i_j, &vbas.asvc.A_i[aux[i].Index], &vbas.asvc.A_i[index])

            var n, ao mcl.Fr
            n.SetInt64(int64(vbas.asvc.N))
            mcl.FrDiv(&ao, &vbas.asvc.FFT.RootsOfUnity[aux[i].Index*2], &n)
            mcl.FrMul(&ao, &ao, &aux[i].Delta)
            mcl.FrDiv(&omega_i, &ao, &omega_i)

            g[i] = w_i_j
            fr[i] = omega_i
        }
    }
    var temp mcl.G1
    mcl.G1MulVec(&temp, g, fr)
    mcl.G1Add(&temp, &temp, &proofs[index])
    return temp
}

func (vbas *VBAS) UpdateCommitment(digest mcl.G1, req asvc.UpdateReq) mcl.G1 {
    return vbas.asvc.UpdateCommitment(digest, req)
}

func (vbas *VBAS) UpdateAll(proofs []mcl.G1, vector []mcl.Fr, req asvc.UpdateReq, aux []asvc.UpdateReq) ([]mcl.G1, []mcl.Fr, []asvc.UpdateReq) {
    aux = append(aux, req)
    l := uint64(len(aux))
    if l*l < vbas.asvc.N {
        if vbas.args == nil {
            return proofs, vector, aux
        } else {
            panic("error: update all")
        }
    }

    if vbas.args == nil || vbas.args.Done {
        vector = vbas.updateVector(vector, aux)
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
        vbas.asvc.OpenAllStep(vbas.args)
        return proofs, vector, aux
    }
    vbas.args.Step++
    vbas.asvc.OpenAllStep(vbas.args)
    if vbas.args.Done {
        aux = aux[vbas.P:]
    }
    return vbas.args.Proofs, vbas.args.Vector, aux
}

func (vbas *VBAS) UpdateProof(proof mcl.G1, index uint64, req asvc.UpdateReq) mcl.G1 {
    return vbas.asvc.UpdateProof(proof, index, req)
}

func (vbas *VBAS) Aggregate(aggs []asvc.Inp) mcl.G1 {
    return vbas.asvc.Aggregate(aggs)
}

func (vbas *VBAS) VerifySingle(digest mcl.G1, proof mcl.G1, v asvc.Val) bool {
    return vbas.asvc.VerifySingle(digest, proof, v)
}

func (vbas *VBAS) VerifyAggregation(digest mcl.G1, proof mcl.G1, aggvs []asvc.Val) bool {
    return vbas.asvc.VerifyAggregation(digest, proof, aggvs)
}
