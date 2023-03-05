package utils

import (
    "github.com/alinush/go-mcl"
)

// GenerateTestingSetup creates a setup of n values from the given secret. **for testing purposes only**
func GenerateTestingSetup(n uint64) ([]mcl.G1, []mcl.G2) {
    var s mcl.Fr
    s.Random()

    var sPow mcl.Fr
    sPow.SetInt64(1)

    var G mcl.G1
    var H mcl.G2
    G.Random()
    H.Random()

    s1Out := make([]mcl.G1, n+1, n+1)
    s2Out := make([]mcl.G2, n+1, n+1)
    for i := uint64(0); i <= n; i++ {
        mcl.G1Mul(&s1Out[i], &G, &sPow)
        mcl.G2Mul(&s2Out[i], &H, &sPow)
        tmp := sPow
        mcl.FrMul(&sPow, &tmp, &s)
    }
    return s1Out, s2Out
}

// GenerateTestingSetup3 creates a setup of n values from the given secret. **for testing purposes only**
func GenerateTestingSetup3(n1 uint64, n2 uint64, n3 uint64) ([][][]mcl.G1, [][][]mcl.G2) {
    var alpha, beta, gamma mcl.Fr
    alpha.Random()
    beta.Random()
    gamma.Random()

    var alphaPow, betaPow mcl.Fr
    alphaPow.SetInt64(1)

    var sPow mcl.Fr
    sPow.SetInt64(1)

    var G mcl.G1
    var H mcl.G2
    G.Random()
    H.Random()

    s1 := make([][][]mcl.G1, n1+1)
    s2 := make([][][]mcl.G2, n1+1)
    for i := uint64(0); i <= n1; i++ {
        s1[i] = make([][]mcl.G1, n2+1)
        s2[i] = make([][]mcl.G2, n2+1)
        betaPow.SetInt64(1)
        for j := uint64(0); j <= n2; j++ {
            s1[i][j] = make([]mcl.G1, n3+1)
            s2[i][j] = make([]mcl.G2, n3+1)
            var pow mcl.Fr
            mcl.FrMul(&pow, &alphaPow, &betaPow)
            for k := uint64(0); k <= n3; k++ {
                mcl.G1Mul(&s1[i][j][k], &G, &pow)
                mcl.G2Mul(&s2[i][j][k], &H, &pow)
                mcl.FrMul(&pow, &pow, &gamma)
            }
            mcl.FrMul(&betaPow, &betaPow, &beta)
        }
        mcl.FrMul(&alphaPow, &alphaPow, &alpha)
    }
    return s1, s2
}
