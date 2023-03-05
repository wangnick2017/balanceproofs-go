package main

import (
    "github.com/wangnick2017/balanceproofs-go/asvc"
    "github.com/wangnick2017/balanceproofs-go/fft"
    "github.com/wangnick2017/balanceproofs-go/vcs"
    "math/rand"
    "time"
)

func Init() {
    fft.InitGlobals()
    rand.Seed(time.Now().Unix())
}

func main() {
    Init()

    //Careful: asvc.AsvcCommitOpenAll() and vcs.BasicUpdate_Fake() are time-consuming

    //asvc.AsvcCommitOpenAll()
    //vcs.BasicUpdate_Fake()
    vcs.BasicQuery_Fake()
    asvc.AsvcIndvidual_Fake()
    asvc.AsvcAgg_Fake()

    //Careful: vcs.BucketCommitOpenAll() is time-consuming

    //vcs.BucketCommitOpenAll()
    vcs.BucketUpdate_Fake()
    vcs.BucketQuery_Fake()
    vcs.BucketIndvidual_Fake()
    vcs.BucketAgg_Fake()
    vcs.BucketProofSize()
}
