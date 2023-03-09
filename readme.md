# Balanceproofs-go

The implementation of [Balanceproofs](https://eprint.iacr.org/2022/864) in Go.

This repo depends on:
- [go-mcl](https://github.com/alinush/go-mcl) for elliptic curve operations.
- [go-kzg](https://github.com/protolambda/go-kzg) as a reference to implement KZG proofs.


## Instructions

- Install golang
- Install ```mcl```
   ```bash
    $ git clone https://github.com/herumi/mcl
    $ cd mcl/
    $ git checkout 35a39d27 #herumi/mcl v1.35
    $ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    $ cmake --build build
    $ sudo cmake --build build --target install
    $ sudo ldconfig
   ```
- Run ```main.go```
   ```bash
    $ go run main.go
   ```
