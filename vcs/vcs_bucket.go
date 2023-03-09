package vcs

import (
    "github.com/alinush/go-mcl"
    "github.com/wangnick2017/balanceproofs-go/asvc"
    "github.com/wangnick2017/balanceproofs-go/fft"
    "github.com/wangnick2017/balanceproofs-go/utils"
)

type VBUC struct {
    N        uint64
    N1       uint64
    N2       uint64
    N3       uint64
    L        uint8
    L1       uint8
    L2       uint8
    L3       uint8
    P        [][]int
    SecretG1 [][][]mcl.G1
    SecretG2 [][][]mcl.G2
    FFT      *fft.Settings
    rootx    uint64
    rooty    uint64
    rootz    uint64
    G        mcl.G1
    H        mcl.G2
    L_xyz    [][][]mcl.G1
    S_xyz    [][][]mcl.G1
    R_xyz    [][][]mcl.G1
    SS_yz    [][]mcl.G1
    RR_yz    [][]mcl.G1

    args       [][]*asvc.OpenAllStepArg
    asvc       *asvc.ASVC
    xExtFFT_K  [][]mcl.G1
    xExtFFT_JK [][][]mcl.G1
}

func (vbuc *VBUC) toeplitzPart1(s []mcl.G1, n uint64) []mcl.G1 {
    x := make([]mcl.G1, n, n)
    for i, j := uint64(0), n-2; i < n-1; i, j = i+1, j-1 {
        x[i] = s[j]
    }
    x[n-1].Clear()
    n2 := n * 2
    xExt := make([]mcl.G1, n2, n2)
    for i := uint64(0); i < n; i++ {
        xExt[i] = x[i]
    }
    for i := n; i < n2; i++ {
        xExt[i].Clear()
    }
    return vbuc.FFT.FFT_G1(xExt, false)
}

func (vbuc *VBUC) Init(L uint8) {
    vbuc.L = L
    vbuc.L1 = L / 4
    vbuc.L2 = L/2 - vbuc.L1
    vbuc.L3 = vbuc.L - vbuc.L1 - vbuc.L2
    vbuc.N = 1 << L
    vbuc.N1 = 1 << vbuc.L1
    vbuc.N2 = 1 << vbuc.L2
    vbuc.N3 = 1 << vbuc.L3
    vbuc.SecretG1, vbuc.SecretG2 = utils.GenerateTestingSetup3(vbuc.N1, vbuc.N2, vbuc.N3)
    vbuc.G = vbuc.SecretG1[0][0][0]
    vbuc.H = vbuc.SecretG2[0][0][0]
    vbuc.FFT = fft.NewSettings(L + 1)
    vbuc.rootx = vbuc.N2 * vbuc.N3 * 2
    vbuc.rooty = vbuc.N1 * vbuc.N3 * 2
    vbuc.rootz = vbuc.N1 * vbuc.N2 * 2

    vbuc.asvc = &asvc.ASVC{}
    vbuc.asvc.Init(vbuc.L3, vbuc.SecretG1[0][0], vbuc.SecretG2[0][0])
    vbuc.args = make([][]*asvc.OpenAllStepArg, vbuc.N1)
    vbuc.P = make([][]int, vbuc.N1)
    for i := uint64(0); i < vbuc.N1; i++ {
        vbuc.args[i] = make([]*asvc.OpenAllStepArg, vbuc.N2)
        vbuc.P[i] = make([]int, vbuc.N2)
        for j := uint64(0); j < vbuc.N2; j++ {
            vbuc.args[i][j] = nil
        }
    }

    vbuc.xExtFFT_K = make([][]mcl.G1, vbuc.N3)
    for k := uint64(0); k < vbuc.N3; k++ {
        s := make([]mcl.G1, vbuc.N2)
        for j := uint64(0); j < vbuc.N2; j++ {
            s[j] = vbuc.SecretG1[0][j][k]
        }
        vbuc.xExtFFT_K[k] = vbuc.toeplitzPart1(s, vbuc.N2)
    }
    vbuc.xExtFFT_JK = make([][][]mcl.G1, vbuc.N2)
    for j := uint64(0); j < vbuc.N2; j++ {
        vbuc.xExtFFT_JK[j] = make([][]mcl.G1, vbuc.N3)
        for k := uint64(0); k < vbuc.N3; k++ {
            s := make([]mcl.G1, vbuc.N1)
            for i := uint64(0); i < vbuc.N1; i++ {
                s[i] = vbuc.SecretG1[i][j][k]
            }
            vbuc.xExtFFT_JK[j][k] = vbuc.toeplitzPart1(s, vbuc.N1)
        }
    }

    L_z := make([][][]mcl.G1, vbuc.N1+1)
    L_yz := make([][][]mcl.G1, vbuc.N1+1)

    for i := uint64(0); i <= vbuc.N1; i++ {
        L_z[i] = make([][]mcl.G1, vbuc.N2+1)
        L_yz[i] = make([][]mcl.G1, vbuc.N2)
        for j := uint64(0); j <= vbuc.N2; j++ {
            L_z[i][j] = vbuc.FFT.FFT_G1(vbuc.SecretG1[i][j][:vbuc.N3], true)
        }
        for j := uint64(0); j < vbuc.N2; j++ {
            L_yz[i][j] = make([]mcl.G1, vbuc.N3)
        }
    }
    vbuc.SS_yz = make([][]mcl.G1, vbuc.N2)
    vbuc.RR_yz = make([][]mcl.G1, vbuc.N2)
    for j := uint64(0); j < vbuc.N2; j++ {
        vbuc.SS_yz[j] = make([]mcl.G1, vbuc.N3)
        vbuc.RR_yz[j] = make([]mcl.G1, vbuc.N3)
    }
    vbuc.S_xyz = make([][][]mcl.G1, vbuc.N1)
    vbuc.R_xyz = make([][][]mcl.G1, vbuc.N1)
    vbuc.L_xyz = make([][][]mcl.G1, vbuc.N1)
    for i := uint64(0); i < vbuc.N1; i++ {
        vbuc.S_xyz[i] = make([][]mcl.G1, vbuc.N2)
        vbuc.R_xyz[i] = make([][]mcl.G1, vbuc.N2)
        vbuc.L_xyz[i] = make([][]mcl.G1, vbuc.N2)
        for j := uint64(0); j < vbuc.N2; j++ {
            vbuc.S_xyz[i][j] = make([]mcl.G1, vbuc.N3)
            vbuc.R_xyz[i][j] = make([]mcl.G1, vbuc.N3)
            vbuc.L_xyz[i][j] = make([]mcl.G1, vbuc.N3)
        }
    }
    for k := uint64(0); k < vbuc.N3; k++ {
        for i := uint64(0); i <= vbuc.N1; i++ {
            l := make([]mcl.G1, vbuc.N2)
            for j := uint64(0); j < vbuc.N2; j++ {
                l[j] = L_z[i][j][k]
            }
            l_prime := vbuc.FFT.FFT_G1(l, true)
            for j := uint64(0); j < vbuc.N2; j++ {
                L_yz[i][j][k] = l_prime[j]
            }
            if i == 0 {
                ss := vbuc.OpenOnes(vbuc.N2, l)
                h := make([]mcl.G1, vbuc.N2)
                for j := uint64(0); j <= vbuc.N2-2; j++ {
                    var jj mcl.Fr
                    jj.SetInt64(int64(j + 1))
                    mcl.G1Mul(&h[j], &l[vbuc.N2-2-j], &jj)
                }
                h[vbuc.N2-1].Clear()
                u_j := vbuc.FFT.FFT_G1(h, false)
                for j := uint64(0); j < vbuc.N2; j++ {
                    var nFr mcl.Fr
                    nFr.SetInt64(int64(vbuc.N2))
                    mcl.FrDiv(&nFr, &vbuc.FFT.RootsOfUnity[j*vbuc.rooty], &nFr)
                    mcl.G1Mul(&u_j[j], &u_j[j], &nFr)
                }

                for j := uint64(0); j < vbuc.N2; j++ {
                    vbuc.SS_yz[j][k] = ss[j]
                    vbuc.RR_yz[j][k] = u_j[j]
                }
            }
        }

    }
    for k := uint64(0); k < vbuc.N3; k++ {
        for j := uint64(0); j < vbuc.N2; j++ {
            l := make([]mcl.G1, vbuc.N1)
            for i := uint64(0); i < vbuc.N1; i++ {
                l[i] = L_yz[i][j][k]
            }
            l_prime := vbuc.FFT.FFT_G1(l, true)
            ss := vbuc.OpenOnes(vbuc.N1, l)

            h := make([]mcl.G1, vbuc.N1)
            for i := uint64(0); i <= vbuc.N1-2; i++ {
                var jj mcl.Fr
                jj.SetInt64(int64(i + 1))
                mcl.G1Mul(&h[i], &l[vbuc.N1-2-i], &jj)
            }
            h[vbuc.N1-1].Clear()
            u_i := vbuc.FFT.FFT_G1(h, false)
            for i := uint64(0); i < vbuc.N1; i++ {
                var nFr mcl.Fr
                nFr.SetInt64(int64(vbuc.N1))
                mcl.FrDiv(&nFr, &vbuc.FFT.RootsOfUnity[i*vbuc.rootx], &nFr)
                mcl.G1Mul(&u_i[i], &u_i[i], &nFr)
            }

            for i := uint64(0); i < vbuc.N1; i++ {
                vbuc.L_xyz[i][j][k] = l_prime[i]
                vbuc.S_xyz[i][j][k] = ss[i]
                vbuc.R_xyz[i][j][k] = u_i[i]
            }
        }
    }

}

func (vbuc *VBUC) OpenOnes(n uint64, s []mcl.G1) []mcl.G1 {
    //ToeplitzPart1
    x := make([]mcl.G1, n, n)
    for i, j := uint64(0), n-1; i < n; i, j = i+1, j-1 {
        x[i] = s[j]
    }
    n2 := n * 2
    xExt := make([]mcl.G1, n2, n2)
    for i := uint64(0); i < n; i++ {
        xExt[i] = x[i]
    }
    for i := n; i < n2; i++ {
        xExt[i].Clear()
    }
    xExtFFT := vbuc.FFT.FFT_G1(xExt, false)

    //ToeplitzcoeffsStep
    // [last poly item] + [0]*(n+1) + [poly items except first and last]
    toeplitzCoeffs := make([]mcl.Fr, n2, n2)
    toeplitzCoeffs[0].SetInt64(1)
    for i := uint64(1); i < n2; i++ {
        toeplitzCoeffs[i].Clear()
    }

    //ToeplitzPart2
    toeplitzCoeffsFFT := vbuc.FFT.FFT_Fr(toeplitzCoeffs, false)
    hExtFFT := make([]mcl.G1, n2, n2)
    for i := uint64(0); i < n2; i++ {
        mcl.G1Mul(&hExtFFT[i], &xExtFFT[i], &toeplitzCoeffsFFT[i])
    }

    //h := fk.ToeplitzPart3(hExtFFT)
    out := vbuc.FFT.FFT_G1(hExtFFT, true)
    // Only the top half is the Toeplitz product, the rest is padding
    h := out[:len(out)/2]

    return vbuc.FFT.FFT_G1(h, false)
}

func (vbuc *VBUC) InitAux() [][][]asvc.UpdateReq {
    ret := make([][][]asvc.UpdateReq, vbuc.N1)
    for i := uint64(0); i < vbuc.N1; i++ {
        ret[i] = make([][]asvc.UpdateReq, vbuc.N2)
        for j := uint64(0); j < vbuc.N2; j++ {
            ret[i][j] = make([]asvc.UpdateReq, 0)
        }
    }
    return ret
}

func (vbuc *VBUC) GetXYZ(i uint64) (x uint64, y uint64, z uint64) {
    x = i / (vbuc.N2 * vbuc.N3)
    y = (i % (vbuc.N2 * vbuc.N3)) / vbuc.N3
    z = i % vbuc.N3
    return x, y, z
}

func (vbuc *VBUC) GetI(x uint64, y uint64, z uint64) uint64 {
    return x*(vbuc.N2*vbuc.N3) + y*vbuc.N3 + z
}

func (vbuc *VBUC) Commit(vector []mcl.Fr) mcl.G1 {
    var digest mcl.G1
    digest.Clear()

    for i := uint64(0); i < vbuc.N; i++ {
        x, y, z := vbuc.GetXYZ(i)
        var temp mcl.G1
        mcl.G1Mul(&temp, &vbuc.L_xyz[x][y][z], &vector[i])
        mcl.G1Add(&digest, &digest, &temp)
    }
    return digest
}

type BucProofAll struct {
    Pii []mcl.G1
    Psi [][]mcl.G1
    Pi  [][][]mcl.G1
}

type BucProofSingle struct {
    Pii mcl.G1
    Psi mcl.G1
    Pi  mcl.G1
}

type BucProofAgg struct {
    Pii []mcl.G1
    Psi [][]mcl.G1
    Pi  [][]mcl.G1
    X   []uint64
    Y   [][]uint64
    Z   [][][]uint64
}

func (vbuc *VBUC) OpenAll(vector []mcl.Fr) BucProofAll {
    ret := BucProofAll{
        Pii: make([]mcl.G1, vbuc.N1),
        Psi: make([][]mcl.G1, vbuc.N1),
        Pi:  make([][][]mcl.G1, vbuc.N1),
    }
    Poly := make([][][]mcl.Fr, vbuc.N1)
    for i := uint64(0); i < vbuc.N1; i++ {
        ret.Pi[i] = make([][]mcl.G1, vbuc.N2)
        poly := make([][]mcl.Fr, vbuc.N2)
        Poly[i] = make([][]mcl.Fr, vbuc.N2)
        ret.Psi[i] = make([]mcl.G1, vbuc.N2)

        for j := uint64(0); j < vbuc.N2; j++ {
            vec := vector[vbuc.GetI(i, j, 0):vbuc.GetI(i, j, vbuc.N3)]
            ret.Pi[i][j] = vbuc.asvc.OpenAll(vec)

            poly[j] = vbuc.FFT.FFT_Fr(vec, true)

            Poly[i][j] = make([]mcl.Fr, vbuc.N3)
            ret.Psi[i][j].Clear()
        }

        for k := uint64(0); k < vbuc.N3; k++ {
            coef := make([]mcl.Fr, vbuc.N2)
            for j := uint64(0); j < vbuc.N2; j++ {
                coef[j] = poly[j][k]
            }
            coef_prime := vbuc.FFT.FFT_Fr(coef, true)
            for j := uint64(0); j < vbuc.N2; j++ {
                Poly[i][j][k] = coef_prime[j]
            }

            /*Psik := make([]mcl.G1, vbuc.N2)
              //manually compute PsiK
              haha := make([]mcl.G1, len(coef_prime))
              for jj := 0; jj < len(coef_prime); jj++ {
              	haha[jj] = vbuc.SecretG1[0][jj][k]
              }
              for j := uint64(0); j < vbuc.N2; j++ {
              	// divisor = [-index, 1]
              	divisor := [2]mcl.Fr{}
              	divisor[1].SetInt64(1)
              	divisor[0] = vbuc.FFT.RootsOfUnity[j*vbuc.rooty]
              	mcl.FrNeg(&divisor[0], &divisor[0])
              	// quot = poly / divisor
              	quotientPolynomial, _ := utils.PolyDiv(coef_prime, divisor[:])
              	// evaluate quotient poly at shared secret, in G1
              	mcl.G1MulVec(&Psik[j], haha[:len(quotientPolynomial)], quotientPolynomial)

              }*/

            PsiK := vbuc.OpenAllBucket(coef_prime, vbuc.xExtFFT_K[k])
            for j := uint64(0); j < vbuc.N2; j++ {
                mcl.G1Add(&ret.Psi[i][j], &ret.Psi[i][j], &PsiK[j])
            }
        }
        ret.Pii[i].Clear()

        ////check psi
        //for j := uint64(0); j < vbuc.N2; j++ {
        //	var jj mcl.Fr
        //	jj = vbuc.FFT.RootsOfUnity[j*vbuc.rooty]
        //	var betaG2 mcl.G2
        //	mcl.G2Mul(&betaG2, &vbuc.H, &jj)
        //	var betaMinus mcl.G2
        //	mcl.G2Sub(&betaMinus, &vbuc.SecretG2[0][1][0], &betaG2)
        //
        //	var commitmentMinusY mcl.G1
        //	//mcl.G1Sub(&commitmentMinusY, &testtest, &inter[j])
        //	commitmentMinusY = testtest
        //
        //	var e1, e3, e4 mcl.GT
        //	mcl.Pairing(&e1, &commitmentMinusY, &vbuc.H)
        //	mcl.MillerLoop(&e3, &ret.Psi[i][j], &betaMinus)
        //	mcl.MillerLoop(&e4, &inter[j], &vbuc.H)
        //	mcl.GTAdd(&e3, &e3, &e4)
        //
        //	//fmt.Println("psi check:(", i, ",", j, ") ", e1.IsEqual(&e3))
        //}
        //var digest mcl.G1
        //digest.Clear()
        //
        //for j := uint64(0); j < vbuc.N2; j++ {
        //	for k := uint64(0); k < vbuc.N3; k++ {
        //		var temp mcl.G1
        //		mcl.G1Mul(&temp, &vbuc.l_yz[0][j][k], &vector[vbuc.GetI(i, j, k)])
        //		mcl.G1Add(&digest, &digest, &temp)
        //	}
        //}
        //fmt.Println("testtest ", digest.IsEqual(&testtest))

        //var testDigest mcl.G1
        //testDigest.Clear()
        //for j := uint64(0); j < vbuc.N2; j++ {
        //	for k := uint64(0); k < vbuc.N3; k++ {
        //		var temp mcl.G1
        //		mcl.G1Mul(&temp, &vbuc.SecretG1[0][j][k], &Poly[i][j][k])
        //		mcl.G1Add(&testDigest, &testDigest, &temp)
        //	}
        //}
        //fmt.Println("testdigest ", digest.IsEqual(&testDigest))

        //for j := uint64(0); j < vbuc.N2; j++ {
        //	for k := uint64(0); k < vbuc.N3; k++ {
        //		var jj mcl.Fr
        //		jj = vbuc.FFT.RootsOfUnity[j*vbuc.rooty]
        //		var betaG2 mcl.G2
        //		mcl.G2Mul(&betaG2, &vbuc.H, &jj)
        //		var betaMinus mcl.G2
        //		mcl.G2Sub(&betaMinus, &vbuc.SecretG2[0][1][0], &betaG2)
        //
        //		var kk mcl.Fr
        //		kk = vbuc.FFT.RootsOfUnity[k*vbuc.rootz]
        //		var gammaG2 mcl.G2
        //		mcl.G2Mul(&gammaG2, &vbuc.H, &kk)
        //		var gammaMinus mcl.G2
        //		mcl.G2Sub(&gammaMinus, &vbuc.SecretG2[0][0][1], &gammaG2)
        //
        //		var yG1 mcl.G1
        //		mcl.G1Mul(&yG1, &vbuc.G, &vector[vbuc.GetI(i, j, k)])
        //		var commitmentMinusY mcl.G1
        //		mcl.G1Sub(&commitmentMinusY, &yG1, &digest)
        //
        //		var e mcl.GT
        //
        //		s1 := make([]mcl.G1, 3)
        //		s2 := make([]mcl.G2, 3)
        //		s1[0] = commitmentMinusY
        //		s1[1] = ret.Psi[i][j]
        //		s1[2] = ret.Pi[i][j][k]
        //		s2[0] = vbuc.H
        //		s2[1] = betaMinus
        //		s2[2] = gammaMinus
        //		mcl.MillerLoopVec(&e, s1, s2)
        //
        //		mcl.FinalExp(&e, &e)
        //		if !e.IsOne() {
        //			fmt.Println("psi check:(", i, ",", j, ",", k, ") ", e.IsOne())
        //		}
        //
        //	}
        //
        //}
    }
    for j := uint64(0); j < vbuc.N2; j++ {
        for k := uint64(0); k < vbuc.N3; k++ {
            coef := make([]mcl.Fr, vbuc.N1)
            for i := uint64(0); i < vbuc.N1; i++ {
                coef[i] = Poly[i][j][k]
            }
            coef_prime := vbuc.FFT.FFT_Fr(coef, true)

            PiiJK := vbuc.OpenAllBucket(coef_prime, vbuc.xExtFFT_JK[j][k])
            for i := uint64(0); i < vbuc.N1; i++ {
                mcl.G1Add(&ret.Pii[i], &ret.Pii[i], &PiiJK[i])
            }
        }
    }
    return ret
}

func (vbuc *VBUC) OpenAllBucket(poly []mcl.Fr, xExtFFT []mcl.G1) []mcl.G1 {
    //poly := vbuc.FFT.FFT_Fr(vector, true)

    //ToeplitzcoeffsStep
    n := uint64(len(poly))
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
    toeplitzCoeffsFFT := vbuc.FFT.FFT_Fr(toeplitzCoeffs, false)
    hExtFFT := make([]mcl.G1, n2, n2)
    for i := uint64(0); i < n2; i++ {
        mcl.G1Mul(&hExtFFT[i], &xExtFFT[i], &toeplitzCoeffsFFT[i])
    }

    //h := fk.ToeplitzPart3(hExtFFT)
    out := vbuc.FFT.FFT_G1(hExtFFT, true)
    // Only the top half is the Toeplitz product, the rest is padding
    h := out[:len(out)/2]

    return vbuc.FFT.FFT_G1(h, false)
}

func (vbuc *VBUC) Query(index uint64, proofs BucProofAll, aux [][][]asvc.UpdateReq) BucProofSingle {
    x, y, z := vbuc.GetXYZ(index)
    if len(aux[x][y]) == 0 {
        return BucProofSingle{
            Pii: proofs.Pii[x],
            Psi: proofs.Psi[x][y],
            Pi:  proofs.Pi[x][y][z],
        }
    }

    g := make([]mcl.G1, len(aux[x][y]))
    fr := make([]mcl.Fr, len(aux[x][y]))
    for i := 0; i < len(aux[x][y]); i++ {
        if z == aux[x][y][i].Index {
            g[i] = vbuc.asvc.U_i[z]
            fr[i] = aux[x][y][i].Delta
        } else {
            var omega_i mcl.Fr
            mcl.FrSub(&omega_i, &vbuc.asvc.FFT.RootsOfUnity[aux[x][y][i].Index*2], &vbuc.asvc.FFT.RootsOfUnity[z*2])

            var w_i_j mcl.G1
            mcl.G1Sub(&w_i_j, &vbuc.asvc.A_i[aux[x][y][i].Index], &vbuc.asvc.A_i[z])

            var n, ao mcl.Fr
            n.SetInt64(int64(vbuc.asvc.N))
            mcl.FrDiv(&ao, &vbuc.asvc.FFT.RootsOfUnity[aux[x][y][i].Index*2], &n)
            mcl.FrMul(&ao, &ao, &aux[x][y][i].Delta)
            mcl.FrDiv(&omega_i, &ao, &omega_i)

            g[i] = w_i_j
            fr[i] = omega_i
        }
    }
    var temp mcl.G1
    mcl.G1MulVec(&temp, g, fr)
    mcl.G1Add(&temp, &temp, &proofs.Pi[x][y][z])
    return BucProofSingle{
        Pii: proofs.Pii[x],
        Psi: proofs.Psi[x][y],
        Pi:  temp,
    }
}

func (vbuc *VBUC) UpdateCommitment(digest mcl.G1, req asvc.UpdateReq) mcl.G1 {
    var temp mcl.G1
    x, y, z := vbuc.GetXYZ(req.Index)
    mcl.G1Mul(&temp, &vbuc.L_xyz[x][y][z], &req.Delta)
    mcl.G1Add(&temp, &digest, &temp)
    return temp
}

func (vbuc *VBUC) UpdateAll(proofs *BucProofAll, vector []mcl.Fr, req asvc.UpdateReq, aux [][][]asvc.UpdateReq) ([]mcl.Fr, [][][]asvc.UpdateReq) {
    x, y, z := vbuc.GetXYZ(req.Index)
    var temp mcl.G1
    for x_ := uint64(0); x_ < vbuc.N1; x_++ {
        if x_ == x {
            mcl.G1Mul(&temp, &vbuc.R_xyz[x][y][z], &req.Delta)
        } else {
            mcl.G1Sub(&temp, &vbuc.S_xyz[x][y][z], &vbuc.S_xyz[x_][y][z])
            var u mcl.Fr
            u.SetInt64(1)
            mcl.FrSub(&u, &u, &vbuc.FFT.RootsOfUnity[(x_-x)*vbuc.rootx])
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
            mcl.G1Mul(&temp, &vbuc.RR_yz[y][z], &req.Delta)

        } else {
            mcl.G1Sub(&temp, &vbuc.SS_yz[y][z], &vbuc.SS_yz[y_][z])

            var u mcl.Fr
            u.SetInt64(1)
            mcl.FrSub(&u, &u, &vbuc.FFT.RootsOfUnity[(y_-y)*vbuc.rooty])
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
    if l < (1<<(vbuc.L3/2))*(uint64(vbuc.L3)) {
        if vbuc.args[x][y] == nil {
            return vector, aux
        } else {
            panic("error: update all")
        }
    }

    if vbuc.args[x][y] == nil || vbuc.args[x][y].Done {
        for i := 0; i < len(aux[x][y]); i++ {
            temp := vector[vbuc.GetI(x, y, aux[x][y][i].Index)]
            mcl.FrAdd(&vector[vbuc.GetI(x, y, aux[x][y][i].Index)], &temp, &aux[x][y][i].Delta)
        }
        vbuc.P[x][y] = len(aux[x][y])
        vbuc.args[x][y] = &asvc.OpenAllStepArg{
            Vector: vector[vbuc.GetI(x, y, 0):vbuc.GetI(x, y, vbuc.N3)],
            Step:   0,
            Proofs: proofs.Pi[x][y],
            Count:  1 << (vbuc.L3 / 2),
            Done:   false,
        }
        vbuc.asvc.OpenAllStep(vbuc.args[x][y])
        return vector, aux
    }
    vbuc.args[x][y].Step++
    vbuc.asvc.OpenAllStep(vbuc.args[x][y])
    if vbuc.args[x][y].Done {
        aux[x][y] = aux[x][y][vbuc.P[x][y]:]
    }
    proofs.Pi[x][y] = vbuc.args[x][y].Proofs
    return vector, aux
}

func (vbuc *VBUC) UpdateProof(proof BucProofSingle, index uint64, req asvc.UpdateReq) BucProofSingle {
    x, y, z := vbuc.GetXYZ(req.Index)
    x_, y_, z_ := vbuc.GetXYZ(index)

    ret := BucProofSingle{
        Pii: mcl.G1{},
        Psi: mcl.G1{},
        Pi:  mcl.G1{},
    }
    var temp mcl.G1
    if x == x_ {
        mcl.G1Mul(&temp, &vbuc.R_xyz[x][y][z], &req.Delta)
        mcl.G1Add(&ret.Pii, &proof.Pii, &temp)

        if y == y_ {
            mcl.G1Mul(&temp, &vbuc.RR_yz[y][z], &req.Delta)
            mcl.G1Add(&ret.Psi, &proof.Psi, &temp)
            ret.Pi = vbuc.asvc.UpdateProof(proof.Pi, z_, asvc.UpdateReq{
                Index: z,
                Delta: req.Delta,
            })
        } else {
            mcl.G1Sub(&temp, &vbuc.SS_yz[y][z], &vbuc.SS_yz[y_][z])

            var u mcl.Fr
            u.SetInt64(1)
            mcl.FrSub(&u, &u, &vbuc.FFT.RootsOfUnity[(y_-y)*vbuc.rooty])
            var uu mcl.Fr
            uu.SetInt64(int64(vbuc.N2))
            mcl.FrMul(&u, &u, &uu)
            mcl.FrDiv(&u, &req.Delta, &u)

            mcl.G1Mul(&temp, &temp, &u)
            mcl.G1Add(&ret.Psi, &proof.Psi, &temp)

            ret.Pi = proof.Pi
        }
    } else {
        mcl.G1Sub(&temp, &vbuc.S_xyz[x][y][z], &vbuc.S_xyz[x_][y][z])

        var u mcl.Fr
        u.SetInt64(1)
        mcl.FrSub(&u, &u, &vbuc.FFT.RootsOfUnity[(x_-x)*vbuc.rootx])
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

func (vbuc *VBUC) Aggregate(indices []uint64, proofs []BucProofSingle) BucProofAgg {
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
            agg.Pi[xx][yy] = vbuc.asvc.Aggregate(temp)

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
            agg.Pi[xx][yy] = vbuc.asvc.Aggregate(temp)

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
    agg.Pi[xx][yy] = vbuc.asvc.Aggregate(temp)

    return agg
}

func (vbuc *VBUC) VerifySingle(digest mcl.G1, proof BucProofSingle, v asvc.Val) bool {
    x, y, z := vbuc.GetXYZ(v.Index)

    var i mcl.Fr
    i = vbuc.FFT.RootsOfUnity[x*vbuc.rootx]
    var alphaG2 mcl.G2
    mcl.G2Mul(&alphaG2, &vbuc.H, &i)
    var alphaMinus mcl.G2
    mcl.G2Sub(&alphaMinus, &vbuc.SecretG2[1][0][0], &alphaG2)

    var j mcl.Fr
    j = vbuc.FFT.RootsOfUnity[y*vbuc.rooty]
    var betaG2 mcl.G2
    mcl.G2Mul(&betaG2, &vbuc.H, &j)
    var betaMinus mcl.G2
    mcl.G2Sub(&betaMinus, &vbuc.SecretG2[0][1][0], &betaG2)

    var k mcl.Fr
    k = vbuc.FFT.RootsOfUnity[z*vbuc.rootz]
    var gammaG2 mcl.G2
    mcl.G2Mul(&gammaG2, &vbuc.H, &k)
    var gammaMinus mcl.G2
    mcl.G2Sub(&gammaMinus, &vbuc.SecretG2[0][0][1], &gammaG2)

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

func (vbuc *VBUC) VerifyAggregation(digest mcl.G1, proof BucProofAgg, aggvs [][][]mcl.Fr) bool {
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
            mcl.FrMul(&temp, &r, &vbuc.FFT.RootsOfUnity[proof.X[x]*vbuc.rootx])
            mcl.FrSub(&coefPiiFr[x], &coefPiiFr[x], &temp)

            proofPsiBeta = append(proofPsiBeta, proof.Psi[x][y])
            coefPsiBeta = append(coefPsiBeta, r)
            mcl.FrMul(&temp, &r, &vbuc.FFT.RootsOfUnity[proof.Y[x][y]*vbuc.rooty])
            mcl.FrNeg(&temp, &temp)
            coefPsiFr = append(coefPsiFr, temp)

            l := len(aggvs[x][y])
            if l == 1 {
                proofPiGamma = append(proofPiGamma, proof.Pi[x][y])
                coefPiGamma = append(coefPiGamma, r)
                mcl.FrMul(&temp, &r, &vbuc.FFT.RootsOfUnity[proof.Z[x][y][0]*vbuc.rootz])
                mcl.FrNeg(&temp, &temp)
                coefPiFr = append(coefPiFr, temp)
                mcl.FrMul(&temp, &r, &aggvs[x][y][0])
                mcl.FrAdd(&coefVal, &coefVal, &temp)
            } else {
                xx := make([]mcl.Fr, l)

                for k := 0; k < l; k++ {
                    xx[k] = vbuc.FFT.RootsOfUnity[proof.Z[x][y][k]*vbuc.rootz]
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
                var sp2 mcl.G2
                mcl.G2MulVec(&sp2, vbuc.SecretG2[0][0][:len(subProductPoly)], subProductPoly)

                // [interpolation_polynomial(s)]_1
                var is1 mcl.G1
                mcl.G1MulVec(&is1, vbuc.SecretG1[0][0][:len(interpolationPoly)], interpolationPoly)

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
    s2 = append(s2, vbuc.SecretG2[1][0][0])
    var psi mcl.G1
    mcl.G1MulVec(&psi, proofPsiBeta, coefPsiBeta)
    s1 = append(s1, psi)
    s2 = append(s2, vbuc.SecretG2[0][1][0])
    if len(proofPiGamma) > 0 {
        var pi mcl.G1
        mcl.G1MulVec(&pi, proofPiGamma, coefPiGamma)
        s1 = append(s1, pi)
        s2 = append(s2, vbuc.SecretG2[0][0][1])
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
