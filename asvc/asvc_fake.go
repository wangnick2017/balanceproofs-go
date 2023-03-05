package asvc

import (
    "github.com/alinush/go-mcl"
    "github.com/wangnick2017/balanceproofs-go/fft"
    "github.com/wangnick2017/balanceproofs-go/utils"
)

func (asvc *ASVC) Init_Fake(L uint8) {
    asvc.L = L
    asvc.N = 1 << L

    asvc.G.Random()
    asvc.H.Random()
    asvc.FFT = &fft.Settings{}
    asvc.FFT.RootsOfUnity = make([]mcl.Fr, 1500)

    //ToeplitzPart1
    asvc.xExtFFT = make([]mcl.G1, 2)
    asvc.xExtFFT[0].Random()
    asvc.xExtFFT[1].Random()

    //compute asvc.L_i, a_i, u_i, a
    asvc.L_i = make([]mcl.G1, 2)
    asvc.L_i[0].Random()
    asvc.L_i[1].Random()
    asvc.A_i = make([]mcl.G1, 2)
    asvc.A_i[0].Random()
    asvc.A_i[1].Random()

    asvc.U_i = make([]mcl.G1, 2)
    asvc.U_i[0].Random()
    asvc.U_i[1].Random()

    asvc.SecretG1 = make([]mcl.G1, 1500)
    asvc.SecretG2 = make([]mcl.G2, 1500)
    for i := 0; i < 1500; i++ {
        asvc.SecretG1[i].Random()
        asvc.SecretG2[i].Random()
        asvc.FFT.RootsOfUnity[i].Random()
    }
}

func (asvc *ASVC) Commit_Fake(vector []mcl.Fr) mcl.G1 {
    var digest mcl.G1
    digest.Clear()

    for i := uint64(0); i < asvc.N; i++ {
        var temp mcl.G1
        mcl.G1Mul(&temp, &asvc.L_i[i&1], &vector[i&1])
        mcl.G1Add(&digest, &digest, &temp)
    }

    return digest
}

func (asvc *ASVC) OpenAll_Fake(vector []mcl.Fr) []mcl.G1 {
    ret := make([]mcl.G1, 10)
    for i := uint64(0); i < 10; i++ {
        ret[i].Random()
    }

    return ret
}

func (asvc *ASVC) OpenAllStep_Fake(openArg *OpenAllStepArg) {
    n := asvc.N
    n2 := n * 2

    switch {
    case openArg.Step == 0:
        openArg.firstFFTDone = false
        poly := asvc.FFT.FFT_Fr_Fake(n, true)

        //ToeplitzcoeffsStep

        // [last poly item] + [0]*(n+1) + [poly items except first and last]
        openArg.toeplitzCoeffs = make([]mcl.Fr, 2)
        openArg.toeplitzCoeffs[0] = poly[(n-1)&1]
        for i := uint64(1); i <= n+1; i++ {
            openArg.toeplitzCoeffs[i&1].Clear()
        }
        for i, j := n+2, 1; i < n2; i, j = i+1, j+1 {
            openArg.toeplitzCoeffs[i&1] = poly[j&1]
        }

    case openArg.Step == 1:
        //ToeplitzPart2
        openArg.toeplitzCoeffsFFT = asvc.FFT.FFT_Fr_Fake(n2, false)

    case openArg.Step > 1 && openArg.Step <= 1+n2/openArg.Count:
        if openArg.Step == 2 {
            openArg.hExtFFT = make([]mcl.G1, 2)
        }
        for i, j := (openArg.Step-2)*openArg.Count, uint64(0); j < openArg.Count; i, j = i+1, j+1 {
            mcl.G1Mul(&openArg.hExtFFT[i&1], &asvc.xExtFFT[i&1], &openArg.toeplitzCoeffsFFT[i&1])
        }

    default:
        if openArg.Step == 2+n2/openArg.Count {
            openArg.fftArg = &fft.FFT_G1_Arg{
                N:       n2,
                L:       uint64(asvc.L) + 1,
                Vals:    openArg.hExtFFT,
                Out:     nil,
                Inv:     true,
                Pre_s:   0,
                Pre_j:   0,
                Pre_inv: 0,
                Count:   openArg.Count,
                Done:    false,
            }
        }
        if !openArg.firstFFTDone {
            asvc.FFT.FFT_G1_Seg_Fake(openArg.fftArg)
            if openArg.fftArg.Done {
                openArg.h = openArg.fftArg.Out[:2]
                openArg.firstFFTDone = true
                openArg.fftArg = nil
            }
        } else {
            if openArg.fftArg == nil {
                openArg.fftArg = &fft.FFT_G1_Arg{
                    N:       n,
                    L:       uint64(asvc.L),
                    Vals:    openArg.h,
                    Out:     nil,
                    Inv:     false,
                    Pre_s:   0,
                    Pre_j:   0,
                    Pre_inv: 0,
                    Count:   openArg.Count,
                    Done:    false,
                }
            }
            asvc.FFT.FFT_G1_Seg_Fake(openArg.fftArg)
            if openArg.fftArg.Done {
                openArg.Proofs = openArg.fftArg.Out
                openArg.Done = true
            }
        }
    }
}

func (asvc *ASVC) UpdateCommitment_Fake(digest mcl.G1, req UpdateReq) mcl.G1 {
    var temp mcl.G1
    mcl.G1Mul(&temp, &asvc.L_i[req.Index&1], &req.Delta)
    mcl.G1Add(&temp, &digest, &temp)
    return temp
}

func (asvc *ASVC) UpdateProof_Fake(proof mcl.G1, index uint64, req UpdateReq) mcl.G1 {
    var roots mcl.Fr
    roots.Random()
    if index == req.Index {
        var temp mcl.G1
        mcl.G1Mul(&temp, &asvc.U_i[index&1], &req.Delta)
        mcl.G1Add(&temp, &proof, &temp)
        return temp
    }

    var omega_i mcl.Fr
    mcl.FrAdd(&omega_i, &roots, &roots)

    var w_i_j mcl.G1
    mcl.G1Sub(&w_i_j, &asvc.A_i[req.Index&1], &asvc.A_i[index&1])

    var n, ao mcl.Fr
    n.SetInt64(int64(asvc.N))
    mcl.FrDiv(&ao, &roots, &n)
    mcl.FrMul(&ao, &ao, &req.Delta)
    mcl.FrDiv(&omega_i, &ao, &omega_i)

    mcl.G1Mul(&w_i_j, &w_i_j, &omega_i)

    var temp mcl.G1
    mcl.G1Add(&temp, &proof, &w_i_j)
    return temp
}

func (asvc *ASVC) Aggregate_Fake(aggs []Inp) mcl.G1 {
    l := len(aggs)
    if l == 1 {
        return aggs[0].Proof
    }
    coef := make([]mcl.Fr, l)
    for i := 0; i < l; i++ {
        coef[i].Random()
    }
    a_I := utils.SubProductTree(coef)
    a_I_prime := utils.PolyDifferentiate(a_I.Poly)
    c := utils.PolyMultiEvaluate(a_I_prime, a_I)

    var res mcl.G1
    res.Clear()
    for i := 0; i < l; i++ {
        var c_pi mcl.G1
        mcl.FrInv(&c[i], &c[i])
        p := aggs[i].Proof
        p.Random()
        mcl.G1Mul(&c_pi, &p, &c[i])
        mcl.G1Add(&res, &res, &c_pi)
    }
    return res
}

func (asvc *ASVC) VerifySingle_Fake(digest mcl.G1, proof mcl.G1, v Val) bool {
    var i mcl.Fr
    i.Random()
    var xG2 mcl.G2
    mcl.G2Mul(&xG2, &asvc.H, &i)
    var sMinuxX mcl.G2
    mcl.G2Sub(&sMinuxX, &asvc.H, &xG2)

    var yG1 mcl.G1
    mcl.G1Mul(&yG1, &asvc.G, &v.Y)
    var commitmentMinusY mcl.G1
    mcl.G1Sub(&commitmentMinusY, &digest, &yG1)

    // e([commitment - y], [1]) = e([proof],  [s - x])
    var e1, e2 mcl.GT
    mcl.Pairing(&e1, &commitmentMinusY, &asvc.H)
    mcl.Pairing(&e2, &proof, &sMinuxX)
    return e1.IsEqual(&e2)
}

func (asvc *ASVC) VerifyAggregation_Fake(digest mcl.G1, proof mcl.G1, aggvs []Val) bool {
    l := len(aggvs)
    x := make([]mcl.Fr, l)

    for i := 0; i < l; i++ {
        x[i] = asvc.FFT.RootsOfUnity[i]
    }

    //dt := time.Now()
    a_I := utils.SubProductTree(x)
    //duration := time.Since(dt)
    //fmt.Println("haha", duration.Seconds())
    a_I_prime := utils.PolyDifferentiate(a_I.Poly)
    s := utils.PolyMultiEvaluate(a_I_prime, a_I)

    ys := make([]mcl.Fr, l)
    for i := 0; i < l; i++ {
        mcl.FrDiv(&ys[i], &aggvs[i].Y, &s[i])
    }
    interpolationPoly := utils.PolyInterpolation(x, ys, a_I)

    subProductPoly := a_I.Poly

    var sp2 mcl.G2
    mcl.G2MulVec(&sp2, asvc.SecretG2[:len(subProductPoly)], subProductPoly)

    // [interpolation_polynomial(s)]_1
    var is1 mcl.G1
    mcl.G1MulVec(&is1, asvc.SecretG1[:len(interpolationPoly)], interpolationPoly)
    // [commitment - interpolation_polynomial(s)]_1 = [commit]_1 - [interpolation_polynomial(s)]_1
    var commitMinusInterpolation mcl.G1
    mcl.G1Sub(&commitMinusInterpolation, &digest, &is1)

    // Verify the pairing equation
    //
    // e([commitment - interpolation_polynomial(s)], [1]) = e([proof],  [s^n - x^n])
    var e1, e2 mcl.GT
    mcl.Pairing(&e1, &commitMinusInterpolation, &asvc.H)
    mcl.Pairing(&e2, &proof, &sp2)
    return e1.IsEqual(&e2)
}
