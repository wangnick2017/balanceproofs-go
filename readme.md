# Balanceproofs-go

The implementation of Balanceproofs in go.

This repo depends on:
- [go-mcl](https://github.com/alinush/go-mcl/) for elliptic curve operations.

[Balanceproofs]: https://eprint.iacr.org/2022/864

## Instructions

### Software requirements
- Install golang, python
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
