package vcs

import (
    "github.com/alinush/go-mcl"
    "github.com/wangnick2017/balanceproofs-go/asvc"
    "github.com/wangnick2017/balanceproofs-go/fft"
    "github.com/wangnick2017/balanceproofs-go/utils"
)

func (vbuc *VBUC) Init_Fake(L uint8) {
    vbuc.L = L
    vbuc.L1 = L / 4
    vbuc.L2 = L/2 - L/4
    vbuc.L3 = vbuc.L - vbuc.L1 - vbuc.L2
    vbuc.N = 1 << L
    vbuc.N1 = 1 << vbuc.L1
    vbuc.N2 = 1 << vbuc.L2
    vbuc.N3 = 1 << vbuc.L3
    vbuc.G.Random()
    vbuc.H.Random()
    vbuc.FFT = &fft.Settings{}
    vbuc.rootx = vbuc.N2 * vbuc.N3 * 2
    vbuc.rooty = vbuc.N1 * vbuc.N3 * 2
    vbuc.rootz = vbuc.N1 * vbuc.N2 * 2

    vbuc.asvc = &asvc.ASVC{}
    vbuc.asvc.Init_Fake(vbuc.L3)
    vbuc.args = make([][]*asvc.OpenAllStepArg, vbuc.N1)
    vbuc.P = make([][]int, vbuc.N1)
    for i := uint64(0); i < vbuc.N1; i++ {
        vbuc.args[i] = make([]*asvc.OpenAllStepArg, vbuc.N2)
        vbuc.P[i] = make([]int, vbuc.N2)
        for j := uint64(0); j < vbuc.N2; j++ {
            vbuc.args[i][j] = nil
        }
    }

    vbuc.xExtFFT_K = make([][]mcl.G1, 2)
    for k := uint64(0); k < 1; k++ {
        vbuc.xExtFFT_K[k] = make([]mcl.G1, 2)
        vbuc.xExtFFT_K[k][0].Random()
        vbuc.xExtFFT_K[k][1].Random()
    }
    vbuc.xExtFFT_JK = make([][][]mcl.G1, 2)
    for j := uint64(0); j < 2; j++ {
        vbuc.xExtFFT_JK[j] = make([][]mcl.G1, 2)
        for k := uint64(0); k < 2; k++ {
            vbuc.xExtFFT_JK[j][k] = make([]mcl.G1, 2)
            vbuc.xExtFFT_JK[j][k][0].Random()
            vbuc.xExtFFT_JK[j][k][1].Random()
        }
    }

    vbuc.SS_yz = make([][]mcl.G1, 2)
    vbuc.RR_yz = make([][]mcl.G1, 2)
    for j := uint64(0); j < 2; j++ {
        vbuc.SS_yz[j] = make([]mcl.G1, 2)
        vbuc.RR_yz[j] = make([]mcl.G1, 2)
        vbuc.SS_yz[j][0].Random()
        vbuc.SS_yz[j][1].Random()
        vbuc.RR_yz[j][0].Random()
        vbuc.RR_yz[j][1].Random()
    }
    vbuc.S_xyz = make([][][]mcl.G1, 2)
    vbuc.R_xyz = make([][][]mcl.G1, 2)
    vbuc.L_xyz = make([][][]mcl.G1, 2)
    for i := uint64(0); i < 2; i++ {
        vbuc.S_xyz[i] = make([][]mcl.G1, 2)
        vbuc.R_xyz[i] = make([][]mcl.G1, 2)
        vbuc.L_xyz[i] = make([][]mcl.G1, 2)
        for j := uint64(0); j < 2; j++ {
            vbuc.S_xyz[i][j] = make([]mcl.G1, 2)
            vbuc.R_xyz[i][j] = make([]mcl.G1, 2)
            vbuc.L_xyz[i][j] = make([]mcl.G1, 2)
            vbuc.S_xyz[i][j][0].Random()
            vbuc.S_xyz[i][j][1].Random()
            vbuc.R_xyz[i][j][0].Random()
            vbuc.R_xyz[i][j][1].Random()
            vbuc.L_xyz[i][j][0].Random()
            vbuc.L_xyz[i][j][1].Random()
        }
    }
}

func (vbuc *VBUC) Commit_Fake(vector []mcl.Fr) mcl.G1 {
    var digest mcl.G1
    digest.Clear()

    for i := uint64(0); i < vbuc.N; i++ {
        x, y, z := vbuc.GetXYZ(i)
        var temp mcl.G1
        mcl.G1Mul(&temp, &vbuc.L_xyz[x&1][y&1][z&1], &vector[i&1])
        mcl.G1Add(&digest, &digest, &temp)
    }
    return digest
}

func (vbuc *VBUC) OpenAll_Fake(vector []mcl.Fr) BucProofAll {
    ret := BucProofAll{
        Pii: make([]mcl.G1, vbuc.N1),
        Psi: make([][]mcl.G1, vbuc.N1),
        Pi:  make([][][]mcl.G1, vbuc.N1),
    }
    for i := uint64(0); i < vbuc.N1; i++ {
        ret.Pii[i].Random()
        ret.Pi[i] = make([][]mcl.G1, vbuc.N2)
        ret.Psi[i] = make([]mcl.G1, vbuc.N2)

        for j := uint64(0); j < vbuc.N2; j++ {
            ret.Pi[i][j] = make([]mcl.G1, 2)
            ret.Pi[i][j][0].Random()
            ret.Pi[i][j][1].Random()

            ret.Psi[i][j].Random()
        }
    }
    return ret
}

func (vbuc *VBUC) Query_Fake(index uint64, proofs BucProofAll, aux [][][]asvc.UpdateReq) BucProofSingle {
    x, y, z := vbuc.GetXYZ(index)
    var roots mcl.Fr
    roots.Random()

    g := make([]mcl.G1, len(aux[x][y]))
    fr := make([]mcl.Fr, len(aux[x][y]))
    for i := 0; i < len(aux[x][y]); i++ {
        if z == aux[x][y][i].Index {
            g[i] = vbuc.asvc.U_i[z&1]
            fr[i] = aux[x][y][i].Delta
        } else {
            var omega_i mcl.Fr
            mcl.FrAdd(&omega_i, &roots, &roots)

            var w_i_j mcl.G1
            mcl.G1Sub(&w_i_j, &vbuc.asvc.A_i[aux[x][y][i].Index&1], &vbuc.asvc.A_i[z&1])

            var n, ao mcl.Fr
            n.SetInt64(int64(vbuc.asvc.N))
            mcl.FrDiv(&ao, &roots, &n)
            mcl.FrMul(&ao, &ao, &aux[x][y][i].Delta)
            mcl.FrDiv(&omega_i, &ao, &omega_i)

            g[i] = w_i_j
            fr[i] = omega_i
        }
    }
    var temp mcl.G1
    mcl.G1MulVec(&temp, g, fr)
    mcl.G1Add(&temp, &temp, &proofs.Pi[x][y][z&1])
    return BucProofSingle{
        Pii: proofs.Pii[x],
        Psi: proofs.Psi[x][y],
        Pi:  temp,
    }
}

func (vbuc *VBUC) UpdateCommitment_Fake(digest mcl.G1, req asvc.UpdateReq) mcl.G1 {
    var temp mcl.G1
    x, y, z := vbuc.GetXYZ(req.Index)
    mcl.G1Mul(&temp, &vbuc.L_xyz[x&1][y&1][z&1], &req.Delta)
    mcl.G1Add(&temp, &digest, &temp)
    return temp
}

func (vbuc *VBUC) UpdateAll_Fake(proofs *BucProofAll, vector []mcl.Fr, req asvc.UpdateReq, aux [][][]asvc.UpdateReq) ([]mcl.Fr, [][][]asvc.UpdateReq) {
    x, y, z := vbuc.GetXYZ(req.Index)
    var roots mcl.Fr
    roots.Random()
    var temp mcl.G1
    for x_ := uint64(0); x_ < vbuc.N1; x_++ {
        if x_ == x {
            mcl.G1Mul(&temp, &vbuc.R_xyz[x&1][y&1][z&1], &req.Delta)
        } else {
            mcl.G1Sub(&temp, &vbuc.S_xyz[x&1][y&1][z&1], &vbuc.S_xyz[x_&1][y&1][z&1])
            var u mcl.Fr
            u.SetInt64(1)
            mcl.FrSub(&u, &u, &roots)
            var uu mcl.Fr
            uu.SetInt64(int64(vbuc.N1))
            mcl.FrMul(&u, &u, &uu)
            mcl.FrDiv(&u, &req.Delta, &u)
            mcl.G1Mul(&temp, &temp, &u)
        }
        mcl.G1Add(&proofs.Pii[x_], &proofs.Pii[x_], &temp)
    }
    for y_ := uint64(0); y_ < vbuc.N2; y_++ {
        if y == y_ {
            mcl.G1Mul(&temp, &vbuc.RR_yz[y&1][z&1], &req.Delta)

        } else {
            mcl.G1Sub(&temp, &vbuc.SS_yz[y&1][z&1], &vbuc.SS_yz[y_&1][z&1])

            var u mcl.Fr
            u.SetInt64(1)
            mcl.FrSub(&u, &u, &roots)
            var uu mcl.Fr
            uu.SetInt64(int64(vbuc.N2))
            mcl.FrMul(&u, &u, &uu)
            mcl.FrDiv(&u, &req.Delta, &u)

            mcl.G1Mul(&temp, &temp, &u)
        }
        mcl.G1Add(&proofs.Psi[x][y_], &proofs.Psi[x][y_], &temp)
    }

    aux[x][y] = append(aux[x][y], req)
    l := uint64(len(aux[x][y]))
    if l < (1<<(vbuc.L3/2))*(uint64(vbuc.L3)) && l < 1024/2 {
        if vbuc.args[x][y] == nil {
            return vector, aux
        } else {
            panic("error: update all")
        }
    }

    if vbuc.args[x][y] == nil || vbuc.args[x][y].Done {
        for i := 0; i < len(aux[x][y]); i++ {
            temp := vector[vbuc.GetI(x, y, aux[x][y][i].Index)&1]
            mcl.FrAdd(&vector[vbuc.GetI(x, y, aux[x][y][i].Index)&1], &temp, &aux[x][y][i].Delta)
        }
        vbuc.P[x][y] = len(aux[x][y])
        vbuc.args[x][y] = &asvc.OpenAllStepArg{
            Vector: vector,
            Step:   0,
            Proofs: proofs.Pi[x][y],
            Count:  1 << (vbuc.L3 / 2),
            Done:   false,
        }
        vbuc.asvc.OpenAllStep_Fake(vbuc.args[x][y])
        return vector, aux
    }
    vbuc.args[x][y].Step++
    vbuc.asvc.OpenAllStep_Fake(vbuc.args[x][y])
    if vbuc.args[x][y].Done {
        aux[x][y] = aux[x][y][vbuc.P[x][y]:]
    }
    proofs.Pi[x][y] = vbuc.args[x][y].Proofs
    return vector, aux
}

func (vbuc *VBUC) UpdateProof_Fake(proof BucProofSingle, index uint64, req asvc.UpdateReq) BucProofSingle {
    x, y, z := vbuc.GetXYZ(req.Index)
    x_, y_, z_ := vbuc.GetXYZ(index)
    var roots mcl.Fr
    roots.Random()

    ret := BucProofSingle{
        Pii: mcl.G1{},
        Psi: mcl.G1{},
        Pi:  mcl.G1{},
    }
    var temp mcl.G1
    if x == x_ {
        mcl.G1Mul(&temp, &vbuc.R_xyz[x&1][y&1][z&1], &req.Delta)
        mcl.G1Add(&ret.Pii, &proof.Pii, &temp)

        if y == y_ {
            mcl.G1Mul(&temp, &vbuc.RR_yz[y&1][z&1], &req.Delta)
            mcl.G1Add(&ret.Psi, &proof.Psi, &temp)
            ret.Pi = vbuc.asvc.UpdateProof(proof.Pi, z_, asvc.UpdateReq{
                Index: z,
                Delta: req.Delta,
            })
        } else {
            mcl.G1Sub(&temp, &vbuc.SS_yz[y&1][z&1], &vbuc.SS_yz[y_&1][z&1])

            var u mcl.Fr
            u.SetInt64(1)
            mcl.FrSub(&u, &u, &roots)
            var uu mcl.Fr
            uu.SetInt64(int64(vbuc.N2))
            mcl.FrMul(&u, &u, &uu)
            mcl.FrDiv(&u, &req.Delta, &u)

            mcl.G1Mul(&temp, &temp, &u)
            mcl.G1Add(&ret.Psi, &proof.Psi, &temp)

            ret.Pi = proof.Pi
        }
    } else {
        mcl.G1Sub(&temp, &vbuc.S_xyz[x&1][y&1][z&1], &vbuc.S_xyz[x_&1][y&1][z&1])

        var u mcl.Fr
        u.SetInt64(1)
        mcl.FrSub(&u, &u, &roots)
        var uu mcl.Fr
        uu.SetInt64(int64(vbuc.N1))
        mcl.FrMul(&u, &u, &uu)
        mcl.FrDiv(&u, &req.Delta, &u)

        mcl.G1Mul(&temp, &temp, &u)
        mcl.G1Add(&ret.Pii, &proof.Pii, &temp)

        ret.Psi = proof.Psi
        ret.Pi = proof.Pi
    }

    return ret
}

func (vbuc *VBUC) Aggregate_Fake(indices []uint64, proofs []BucProofSingle) BucProofAgg {
    agg := BucProofAgg{
        Pii: make([]mcl.G1, 1),
        Psi: make([][]mcl.G1, 1),
        Pi:  make([][]mcl.G1, 1),
        X:   make([]uint64, 1),
        Y:   make([][]uint64, 1),
        Z:   make([][][]uint64, 1),
    }
    x_, y_, z_ := vbuc.GetXYZ(indices[0])
    agg.X[0] = x_
    agg.Y[0] = make([]uint64, 1)
    agg.Y[0][0] = y_
    agg.Z[0] = make([][]uint64, 1)
    agg.Z[0][0] = make([]uint64, 1)
    agg.Z[0][0][0] = z_
    agg.Pii[0] = proofs[0].Pii
    agg.Psi[0] = make([]mcl.G1, 1)
    agg.Psi[0][0] = proofs[0].Psi
    agg.Pi[0] = make([]mcl.G1, 1)
    temp := make([]asvc.Inp, 1)
    temp[0].Proof = proofs[0].Pi
    temp[0].Index = z_
    xx, yy := 0, 0
    for i := 1; i < len(indices); i++ {
        x, y, z := vbuc.GetXYZ(indices[i])
        if x != x_ {
            agg.Pi[xx][yy] = vbuc.asvc.Aggregate_Fake(temp)

            xx++
            agg.X = append(agg.X, x)
            agg.Pii = append(agg.Pii, proofs[i].Pii)
            yy = 0
            agg.Y = append(agg.Y, make([]uint64, 1))
            agg.Y[xx][0] = y
            agg.Psi = append(agg.Psi, make([]mcl.G1, 1))
            agg.Psi[xx][0] = proofs[i].Psi
            agg.Z = append(agg.Z, make([][]uint64, 1))
            agg.Z[xx][0] = make([]uint64, 1)
            agg.Z[xx][0][0] = z
            agg.Pi = append(agg.Pi, make([]mcl.G1, 1))
            temp = make([]asvc.Inp, 1)
            temp[0].Proof = proofs[i].Pi
            temp[0].Index = z
        } else if y != y_ { // x==x_, y!=y_
            agg.Pi[xx][yy] = vbuc.asvc.Aggregate_Fake(temp)

            yy++
            agg.Y[xx] = append(agg.Y[xx], y)
            agg.Psi[xx] = append(agg.Psi[xx], proofs[i].Psi)
            agg.Pi[xx] = append(agg.Pi[xx], mcl.G1{})
            agg.Z[xx] = append(agg.Z[xx], make([]uint64, 1))
            agg.Z[xx][yy][0] = z
            temp = make([]asvc.Inp, 1)
            temp[0].Proof = proofs[i].Pi
            temp[0].Index = z
        } else { // x==x_, y==y_
            agg.Z[xx][yy] = append(agg.Z[xx][yy], z)
            temp = append(temp, asvc.Inp{
                Index: z,
                Proof: proofs[i].Pi,
            })
        }
        x_, y_, z_ = x, y, z
    }
    agg.Pi[xx][yy] = vbuc.asvc.Aggregate_Fake(temp)

    return agg
}

func (vbuc *VBUC) VerifySingle_Fake(digest mcl.G1, proof BucProofSingle, v asvc.Val) bool {
    //x, y, z := vbuc.GetXYZ(v.Index)

    var i mcl.Fr
    i.Random()
    var alphaG2 mcl.G2
    mcl.G2Mul(&alphaG2, &vbuc.H, &i)
    var alphaMinus mcl.G2
    mcl.G2Sub(&alphaMinus, &vbuc.H, &alphaG2)

    var j mcl.Fr
    j.Random()
    var betaG2 mcl.G2
    mcl.G2Mul(&betaG2, &vbuc.H, &j)
    var betaMinus mcl.G2
    mcl.G2Sub(&betaMinus, &vbuc.H, &betaG2)

    var k mcl.Fr
    k.Random()
    var gammaG2 mcl.G2
    mcl.G2Mul(&gammaG2, &vbuc.H, &k)
    var gammaMinus mcl.G2
    mcl.G2Sub(&gammaMinus, &vbuc.H, &gammaG2)

    var yG1 mcl.G1
    mcl.G1Mul(&yG1, &vbuc.G, &v.Y)
    var commitmentMinusY mcl.G1
    mcl.G1Sub(&commitmentMinusY, &yG1, &digest)

    var e mcl.GT

    s1 := make([]mcl.G1, 4)
    s2 := make([]mcl.G2, 4)
    s1[0] = commitmentMinusY
    s1[1] = proof.Pii
    s1[2] = proof.Psi
    s1[3] = proof.Pi
    s2[0] = vbuc.H
    s2[1] = alphaMinus
    s2[2] = betaMinus
    s2[3] = gammaMinus
    mcl.MillerLoopVec(&e, s1, s2)

    mcl.FinalExp(&e, &e)
    return e.IsOne()
}

func (vbuc *VBUC) VerifyAggregation_Fake(digest mcl.G1, proof BucProofAgg, aggvs [][][]mcl.Fr) bool {
    var randFr mcl.Fr
    randFr.Random()
    var randG1 mcl.G1
    randG1.Random()
    var randG2 mcl.G2
    randG2.Random()

    var coefC mcl.Fr
    coefC.Clear()

    coefPiiAlpha := make([]mcl.Fr, len(proof.X))
    coefPiiFr := make([]mcl.Fr, len(proof.X))
    proofPsiBeta := make([]mcl.G1, 0)
    coefPsiBeta := make([]mcl.Fr, 0)
    coefPsiFr := make([]mcl.Fr, 0)
    proofPiGamma := make([]mcl.G1, 0)
    coefPiGamma := make([]mcl.Fr, 0)
    coefPiFr := make([]mcl.Fr, 0)
    var coefVal mcl.Fr
    coefVal.Clear()

    subprodPi := make([]mcl.G1, 0)
    subprodHr := make([]mcl.G2, 0)
    interprodr := make([]mcl.G1, 0)
    for x := 0; x < len(proof.X); x++ {
        coefPiiAlpha[x].Clear()
        coefPiiFr[x].Clear()
        for y := 0; y < len(proof.Y[x]); y++ {
            var r mcl.Fr
            var temp mcl.Fr
            r.Random()
            mcl.FrSub(&coefC, &coefC, &r)
            mcl.FrAdd(&coefPiiAlpha[x], &coefPiiAlpha[x], &r)
            mcl.FrMul(&temp, &r, &randFr)
            mcl.FrSub(&coefPiiFr[x], &coefPiiFr[x], &temp)

            proofPsiBeta = append(proofPsiBeta, proof.Psi[x][y])
            coefPsiBeta = append(coefPsiBeta, r)
            mcl.FrMul(&temp, &r, &randFr)
            mcl.FrNeg(&temp, &temp)
            coefPsiFr = append(coefPsiFr, temp)

            l := len(aggvs[x][y])
            if l == 1 {
                proofPiGamma = append(proofPiGamma, proof.Pi[x][y])
                coefPiGamma = append(coefPiGamma, r)
                mcl.FrMul(&temp, &r, &randFr)
                mcl.FrNeg(&temp, &temp)
                coefPiFr = append(coefPiFr, temp)
                mcl.FrMul(&temp, &r, &aggvs[x][y][0])
                mcl.FrAdd(&coefVal, &coefVal, &temp)
            } else {
                xx := make([]mcl.Fr, l)

                for k := 0; k < l; k++ {
                    xx[k] = randFr
                }

                a_I := utils.SubProductTree(xx)
                a_I_prime := utils.PolyDifferentiate(a_I.Poly)
                s := utils.PolyMultiEvaluate(a_I_prime, a_I)

                ys := make([]mcl.Fr, l)
                for k := 0; k < l; k++ {
                    mcl.FrDiv(&ys[k], &aggvs[x][y][k], &s[k])
                }
                interpolationPoly := utils.PolyInterpolation(xx, ys, a_I)

                subProductPoly := a_I.Poly
                for k := 0; k < len(subProductPoly); k++ {
                    mcl.FrMul(&subProductPoly[k], &subProductPoly[k], &r)
                }
                for k := 0; k < len(interpolationPoly); k++ {
                    mcl.FrMul(&interpolationPoly[k], &interpolationPoly[k], &r)
                }
                randl2 := make([]mcl.G2, len(subProductPoly))
                for k := 0; k < len(subProductPoly); k++ {
                    randl2[k] = randG2
                }
                var sp2 mcl.G2
                mcl.G2MulVec(&sp2, randl2, subProductPoly)

                randl1 := make([]mcl.G1, len(interpolationPoly))
                for k := 0; k < len(interpolationPoly); k++ {
                    randl1[k] = randG1
                }
                // [interpolation_polynomial(s)]_1
                var is1 mcl.G1
                mcl.G1MulVec(&is1, randl1, interpolationPoly)

                subprodPi = append(subprodPi, proof.Pi[x][y])
                subprodHr = append(subprodHr, sp2)
                interprodr = append(interprodr, is1)
            }
        }
    }

    s1 := make([]mcl.G1, 0)
    s2 := make([]mcl.G2, 0)
    var pii mcl.G1
    mcl.G1MulVec(&pii, proof.Pii, coefPiiAlpha)
    s1 = append(s1, pii)
    s2 = append(s2, randG2)
    var psi mcl.G1
    mcl.G1MulVec(&psi, proofPsiBeta, coefPsiBeta)
    s1 = append(s1, psi)
    s2 = append(s2, randG2)
    if len(proofPiGamma) > 0 {
        var pi mcl.G1
        mcl.G1MulVec(&pi, proofPiGamma, coefPiGamma)
        s1 = append(s1, pi)
        s2 = append(s2, randG2)
    }
    for i := 0; i < len(subprodPi); i++ {
        s1 = append(s1, subprodPi[i])
        s2 = append(s2, subprodHr[i])
    }
    gleft := make([]mcl.G1, 2)
    gright := make([]mcl.Fr, 2)
    gleft[0] = digest
    gright[0] = coefC
    gleft[1] = vbuc.G
    gright[1] = coefVal
    for i := 0; i < len(proof.X); i++ {
        gleft = append(gleft, proof.Pii[i])
        gright = append(gright, coefPiiFr[i])
    }
    for i := 0; i < len(coefPsiFr); i++ {
        gleft = append(gleft, proofPsiBeta[i])
        gright = append(gright, coefPsiFr[i])
    }
    for i := 0; i < len(coefPiFr); i++ {
        gleft = append(gleft, proofPiGamma[i])
        gright = append(gright, coefPiFr[i])
    }
    var g1left mcl.G1
    mcl.G1MulVec(&g1left, gleft, gright)
    for i := 0; i < len(interprodr); i++ {
        mcl.G1Add(&g1left, &g1left, &interprodr[i])
    }
    s1 = append(s1, g1left)
    s2 = append(s2, vbuc.H)
    var e mcl.GT
    mcl.MillerLoopVec(&e, s1, s2)

    mcl.FinalExp(&e, &e)
    return e.IsOne()
}
