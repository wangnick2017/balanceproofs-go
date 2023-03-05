package utils

import (
    "math/rand"
    "sort"
)

func GenerateIndices(n uint64, m int) []uint64 {
    ret := make([]uint64, m)
    for i := 0; i < m; i = i + 1 {

        temp := rand.Uint64() % (n - 3)
        for j := 0; j < i; j++ {
            if ret[j] == temp {
                temp = temp + 1
            }
        }
        ret[i] = temp

        sort.Slice(ret[0:i+1], func(i, j int) bool { return ret[i] < ret[j] })
    }
    return ret
}
