package fft

import (
    "github.com/alinush/go-mcl"
)

type Settings struct {
    N uint64
    // the generator used to get all roots of unity
    RootOfUnity *mcl.Fr
    // domain, starting and ending with 1 (duplicate!)
    RootsOfUnity        []mcl.Fr
    ReverseRootsOfUnity []mcl.Fr
}

func NewSettings(L uint8) *Settings {
    width := uint64(1) << L
    root := &Scale2RootOfUnity[L]
    rootz := make([]mcl.Fr, 2)
    rootz[0].SetInt64(1)
    rootz[1] = Scale2RootOfUnity[L]
    for i := 1; !rootz[i].IsOne(); {
        rootz = append(rootz, mcl.Fr{})
        this := &rootz[i]
        i++
        mcl.FrMul(&rootz[i], this, &Scale2RootOfUnity[L])
    }

    rootzReverse := make([]mcl.Fr, len(rootz), len(rootz))
    copy(rootzReverse, rootz)
    for i, j := uint64(0), uint64(len(rootz)-1); i < j; i, j = i+1, j-1 {
        rootzReverse[i], rootzReverse[j] = rootzReverse[j], rootzReverse[i]
    }

    return &Settings{
        N:                   width,
        RootOfUnity:         root,
        RootsOfUnity:        rootz,
        ReverseRootsOfUnity: rootzReverse,
    }
}
