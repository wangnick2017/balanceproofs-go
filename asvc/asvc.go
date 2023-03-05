package asvc

import (
    "github.com/alinush/go-mcl"
    "github.com/wangnick2017/balanceproofs-go/fft"
    "github.com/wangnick2017/balanceproofs-go/utils"
)

type UpdateReq struct {
    Index uint64
    Delta mcl.Fr
}

type Inp struct {
    Index uint64
    Proof mcl.G1
}

type Val struct {
    Index uint64
    Y     mcl.Fr
}

type ASVC struct {
    N        uint64
    L        uint8
    SecretG1 []mcl.G1
    SecretG2 []mcl.G2
    G        mcl.G1
    H        mcl.G2
    FFT      *fft.Settings
    L_i      []mcl.G1
    U_i      []mcl.G1
    A_i      []mcl.G1
    a        mcl.G1
    xExtFFT  []mcl.G1
}

func (asvc *ASVC) Init(L uint8, s1 []mcl.G1, s2 []mcl.G2) {
    asvc.L = L
    asvc.N = 1 << L
    if len(s1) == 0 {
        asvc.SecretG1, asvc.SecretG2 = utils.GenerateTestingSetup(asvc.N)
    } else {
        asvc.SecretG1 = s1
        asvc.SecretG2 = s2
    }

    asvc.G = asvc.SecretG1[0]
    asvc.H = asvc.SecretG2[0]
    asvc.FFT = fft.NewSettings(L + 1)

    //ToeplitzPart1
    x := make([]mcl.G1, asvc.N, asvc.N)
    for i, j := uint64(0), asvc.N-2; i < asvc.N-1; i, j = i+1, j-1 {
        x[i] = asvc.SecretG1[j]
    }
    x[asvc.N-1].Clear()
    n2 := asvc.N * 2
    xExt := make([]mcl.G1, n2, n2)
    for i := uint64(0); i < asvc.N; i++ {
        xExt[i] = x[i]
    }
    for i := asvc.N; i < n2; i++ {
        xExt[i].Clear()
    }
    asvc.xExtFFT = asvc.FFT.FFT_G1(xExt, false)

    //compute asvc.L_i, a_i, u_i, a
    mcl.G1Sub(&asvc.a, &asvc.SecretG1[asvc.N], &asvc.SecretG1[0])
    asvc.L_i = asvc.FFT.FFT_G1(asvc.SecretG1[:asvc.N], true)
    asvc.A_i = asvc.OpenOnes()

    h := make([]mcl.G1, asvc.N)
    for i := uint64(0); i <= asvc.N-2; i++ {
        var j mcl.Fr
        j.SetInt64(int64(i + 1))
        mcl.G1Mul(&h[i], &asvc.SecretG1[asvc.N-2-i], &j)
    }
    h[asvc.N-1].Clear()

    asvc.U_i = asvc.FFT.FFT_G1(h, false)
    for i := uint64(0); i < asvc.N; i++ {
        var nFr mcl.Fr
        nFr.SetInt64(int64(asvc.N))
        mcl.FrDiv(&nFr, &asvc.FFT.RootsOfUnity[i*2], &nFr)
        mcl.G1Mul(&asvc.U_i[i], &asvc.U_i[i], &nFr)
    }
}

func (asvc *ASVC) Commit(vector []mcl.Fr) mcl.G1 {
    var digest mcl.G1
    mcl.G1MulVec(&digest, asvc.L_i, vector)

    return digest
}

// Open could also be implemented in O(n) time through n calls of UpdateProof()
func (asvc *ASVC) Open(index uint64, vector []mcl.Fr) mcl.G1 {
    poly := asvc.FFT.FFT_Fr(vector, true)

    // divisor = [-index, 1]
    divisor := [2]mcl.Fr{}
    divisor[1].SetInt64(1)
    divisor[0] = asvc.FFT.RootsOfUnity[index*2]
    mcl.FrNeg(&divisor[0], &divisor[0])
    // quot = poly / divisor
    quotientPolynomial, _ := utils.PolyDiv(poly, divisor[:])

    // evaluate quotient poly at shared secret, in G1
    var proof mcl.G1
    mcl.G1MulVec(&proof, asvc.SecretG1[:len(quotientPolynomial)], quotientPolynomial)
    return proof
}

func (asvc *ASVC) OpenOnes() []mcl.G1 {
    //ToeplitzPart1
    x := make([]mcl.G1, asvc.N, asvc.N)
    for i, j := uint64(0), asvc.N-1; i < asvc.N; i, j = i+1, j-1 {
        x[i] = asvc.SecretG1[j]
    }
    n2 := asvc.N * 2
    xExt := make([]mcl.G1, n2, n2)
    for i := uint64(0); i < asvc.N; i++ {
        xExt[i] = x[i]
    }
    for i := asvc.N; i < n2; i++ {
        xExt[i].Clear()
    }
    xExtFFT := asvc.FFT.FFT_G1(xExt, false)

    //ToeplitzcoeffsStep
    // [last poly item] + [0]*(n+1) + [poly items except first and last]
    toeplitzCoeffs := make([]mcl.Fr, n2, n2)
    toeplitzCoeffs[0].SetInt64(1)
    for i := uint64(1); i < n2; i++ {
        toeplitzCoeffs[i].Clear()
    }

    //ToeplitzPart2
    toeplitzCoeffsFFT := asvc.FFT.FFT_Fr(toeplitzCoeffs, false)
    hExtFFT := make([]mcl.G1, n2, n2)
    for i := uint64(0); i < n2; i++ {
        mcl.G1Mul(&hExtFFT[i], &xExtFFT[i], &toeplitzCoeffsFFT[i])
    }

    //h := fk.ToeplitzPart3(hExtFFT)
    out := asvc.FFT.FFT_G1(hExtFFT, true)
    // Only the top half is the Toeplitz product, the rest is padding
    h := out[:len(out)/2]

    return asvc.FFT.FFT_G1(h, false)
}

func (asvc *ASVC) OpenAll(vector []mcl.Fr) []mcl.G1 {
    poly := asvc.FFT.FFT_Fr(vector, true)

    //ToeplitzcoeffsStep
    n := asvc.N
    n2 := n * 2
    // [last poly item] + [0]*(n+1) + [poly items except first and last]
    toeplitzCoeffs := make([]mcl.Fr, n2, n2)
    toeplitzCoeffs[0] = poly[n-1]
    for i := uint64(1); i <= n+1; i++ {
        toeplitzCoeffs[i].Clear()
    }
    for i, j := n+2, 1; i < n2; i, j = i+1, j+1 {
        toeplitzCoeffs[i] = poly[j]
    }

    //ToeplitzPart2
    toeplitzCoeffsFFT := asvc.FFT.FFT_Fr(toeplitzCoeffs, false)
    hExtFFT := make([]mcl.G1, n2, n2)
    for i := uint64(0); i < n2; i++ {
        mcl.G1Mul(&hExtFFT[i], &asvc.xExtFFT[i], &toeplitzCoeffsFFT[i])
    }

    //h := fk.ToeplitzPart3(hExtFFT)
    out := asvc.FFT.FFT_G1(hExtFFT, true)
    // Only the top half is the Toeplitz product, the rest is padding
    h := out[:len(out)/2]

    return asvc.FFT.FFT_G1(h, false)
}

type OpenAllStepArg struct {
    Vector []mcl.Fr
    Step   uint64
    Proofs []mcl.G1
    Done   bool
    Count  uint64

    toeplitzCoeffs    []mcl.Fr
    toeplitzCoeffsFFT []mcl.Fr
    hExtFFT           []mcl.G1
    h                 []mcl.G1
    firstFFTDone      bool
    fftArg            *fft.FFT_G1_Arg
}

func (asvc *ASVC) OpenAllStep(openArg *OpenAllStepArg) {
    n := asvc.N
    n2 := n * 2

    switch {
    case openArg.Step == 0:
        openArg.firstFFTDone = false
        poly := asvc.FFT.FFT_Fr(openArg.Vector, true)

        //ToeplitzcoeffsStep

        // [last poly item] + [0]*(n+1) + [poly items except first and last]
        openArg.toeplitzCoeffs = make([]mcl.Fr, n2, n2)
        openArg.toeplitzCoeffs[0] = poly[n-1]
        for i := uint64(1); i <= n+1; i++ {
            openArg.toeplitzCoeffs[i].Clear()
        }
        for i, j := n+2, 1; i < n2; i, j = i+1, j+1 {
            openArg.toeplitzCoeffs[i] = poly[j]
        }

    case openArg.Step == 1:
        //ToeplitzPart2
        openArg.toeplitzCoeffsFFT = asvc.FFT.FFT_Fr(openArg.toeplitzCoeffs, false)

    case openArg.Step > 1 && openArg.Step <= 1+n2/openArg.Count:
        if openArg.Step == 2 {
            openArg.hExtFFT = make([]mcl.G1, n2, n2)
        }
        for i, j := (openArg.Step-2)*openArg.Count, uint64(0); j < openArg.Count; i, j = i+1, j+1 {
            mcl.G1Mul(&openArg.hExtFFT[i], &asvc.xExtFFT[i], &openArg.toeplitzCoeffsFFT[i])
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
            asvc.FFT.FFT_G1_Seg(openArg.fftArg)
            if openArg.fftArg.Done {
                openArg.h = openArg.fftArg.Out[:n]
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
            asvc.FFT.FFT_G1_Seg(openArg.fftArg)
            if openArg.fftArg.Done {
                openArg.Proofs = openArg.fftArg.Out
                openArg.Done = true
            }
        }
    }
}

func (asvc *ASVC) UpdateCommitment(digest mcl.G1, req UpdateReq) mcl.G1 {
    var temp mcl.G1
    mcl.G1Mul(&temp, &asvc.L_i[req.Index], &req.Delta)
    mcl.G1Add(&temp, &digest, &temp)
    return temp
}

func (asvc *ASVC) UpdateProof(proof mcl.G1, index uint64, req UpdateReq) mcl.G1 {
    if index == req.Index {
        var temp mcl.G1
        mcl.G1Mul(&temp, &asvc.U_i[index], &req.Delta)
        mcl.G1Add(&temp, &proof, &temp)
        return temp
    }

    var omega_i mcl.Fr
    mcl.FrSub(&omega_i, &asvc.FFT.RootsOfUnity[req.Index*2], &asvc.FFT.RootsOfUnity[index*2])

    var w_i_j mcl.G1
    mcl.G1Sub(&w_i_j, &asvc.A_i[req.Index], &asvc.A_i[index])

    var n, ao mcl.Fr
    n.SetInt64(int64(asvc.N))
    mcl.FrDiv(&ao, &asvc.FFT.RootsOfUnity[req.Index*2], &n)
    mcl.FrMul(&ao, &ao, &req.Delta)
    mcl.FrDiv(&omega_i, &ao, &omega_i)

    mcl.G1Mul(&w_i_j, &w_i_j, &omega_i)

    var temp mcl.G1
    mcl.G1Add(&temp, &proof, &w_i_j)
    return temp
}

func (asvc *ASVC) Aggregate(aggs []Inp) mcl.G1 {
    l := len(aggs)
    if l == 1 {
        return aggs[0].Proof
    }
    coef := make([]mcl.Fr, l)
    for i := 0; i < l; i++ {
        coef[i] = asvc.FFT.RootsOfUnity[aggs[i].Index*2]
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
        mcl.G1Mul(&c_pi, &p, &c[i])
        mcl.G1Add(&res, &res, &c_pi)
    }
    return res
}

func (asvc *ASVC) VerifySingle(digest mcl.G1, proof mcl.G1, v Val) bool {
    var i mcl.Fr
    i = asvc.FFT.RootsOfUnity[v.Index*2]
    var xG2 mcl.G2
    mcl.G2Mul(&xG2, &asvc.H, &i)
    var sMinuxX mcl.G2
    mcl.G2Sub(&sMinuxX, &asvc.SecretG2[1], &xG2)

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

func (asvc *ASVC) VerifyAggregation(digest mcl.G1, proof mcl.G1, aggvs []Val) bool {
    l := len(aggvs)
    x := make([]mcl.Fr, l)

    for i := 0; i < l; i++ {
        x[i] = asvc.FFT.RootsOfUnity[aggvs[i].Index*2]
    }

    a_I := utils.SubProductTree(x)
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
