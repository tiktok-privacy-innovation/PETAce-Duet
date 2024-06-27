# PETAce-Duet

## Introduction

PETAce-Duet is a collection of general-purpose two-party secure computing protocols.
It is one of the crucial building blocks of the framework [PETAce](https://github.com/tiktok-privacy-innovation/PETAce).

As implied by the name "Duet", PETAce-Duet implements various two-party MPC protocols as follows:

- [ABY](https://www.ndss-symposium.org/wp-content/uploads/2017/09/08_2_1.pdf) protocols: arithmetic sharing, boolean sharing, share conversion, and operations over shares (e.g., addition and multiplication over arithmetic shares and bit-wise AND operation over boolean shares).

- Secure comparison protocols from ABY and [Cheetah](https://www.usenix.org/system/files/sec22-huang-zhicong.pdf).

- [Secret shared shuffle](https://link.springer.com/chapter/10.1007/978-3-030-64840-4_12). The permutation network will be available in a future release.

- Homomorphic Encryption (HE): Duet integrates the Paillier cryptosystem to improve two-party secure computation.
For example, HE enables efficient alternative solutions to problems such as secure shuffling and linear matrix operations.
A generic implementation of the Paillier cryptosystem based on the GMP library is provided in PETAce-Solo.
On supported processors, users can switch to the Intel Paillier Cryptosystem Library (IPCL) for extreme performance.
We support the conversion between HE ciphertexts and arithmetic shares. As a result, we achieve an "ABH" framework where arithmetic shares, Boolean shares, and HE shares (secret shares encrypted by HE) can be converted to one another efficiently.

- Fully Homomorphic Encryption (FHE): Duet integrates the BFV scheme to improve two-party Beaver triple generation.
The underlying BFV functionalities depend on Microsoft SEAL.

The core functionalities of Duet are written in C++ to provide the best performance.
A Python interface is provided in [PETAce](https://github.com/tiktok-privacy-innovation/PETAce).

## Requirements

| System | Toolchain                                             |
|--------|-------------------------------------------------------|
| Linux  | Clang++ (>= 5.0) or GNU G++ (>= 5.5), CMake (>= 3.15) |

| Required dependency                                                            | Tested version | Use                               |
|--------------------------------------------------------------------------------|----------------|-----------------------------------|
| [Eigen](https://gitlab.com/libeigen/eigen)                                     | 3.4.0          | Matrix and vector templates       |
| [Intel Paillier Cryptosystem Library (IPCL)](https://github.com/intel/pailliercryptolib)                                     | 2.0.0          | Paillier Encryption       |
| [Microsoft SEAL](https://github.com/microsoft/SEAL)                                     | 4.1.0          | Fully homomorphic encryption       |
| [PETAce-Solo](https://github.com/tiktok-privacy-innovation/PETAce-Solo)       | 0.3.0          | Cryptography primitives           |
| [PETAce-Network](https://github.com/tiktok-privacy-innovation/PETAce-Network) | 0.3.0          | Network communication protocols   |
| [PETAce-Verse](https://github.com/tiktok-privacy-innovation/PETAce-Verse)     | 0.3.0          | Primitive cryptographic protocols |

| Optional dependency                                    | Tested version | Use                    |
|--------------------------------------------------------|----------------|------------------------|
| [GoogleTest](https://github.com/google/googletest)     | 1.12.1         | For running tests      |
| [Google Logging](https://github.com/google/glog)   | 0.4.0          | For running benchmarks |
| [TCLAP](https://github.com/mirror/tclap)           | 1.2.2          | For running examples and benchmarks |

## Building PETAce-Duet

We assume that all commands presented below are executed in the root directory of PETAce-Duet.

To build PETAce-Duet library (optionally with test, benchmark, and example):

```bash
cmake -S . -B build -DDUET_BUILD_TEST=ON -DDUET_BUILD_BENCH=ON -DDUET_BUILD_EXAMPLE=ON
cmake --build build
```

Output binaries can be found in `build/lib/` and `build/bin/` directories.

| Compile Options          | Values        | Default | Description                                         |
|--------------------------|---------------|---------|-----------------------------------------------------|
| `CMAKE_BUILD_TYPE`       | Release/Debug | Release | Debug mode decreases run-time performance.          |
| `DUET_BUILD_SHARED_LIBS` | ON/OFF        | OFF     | Build a shared library if set to ON.                |
| `DUET_BUILD_BENCH`       | ON/OFF        | ON      | Build C++ benchmark if set to ON.                   |
| `DUET_BUILD_EXAMPLE`     | ON/OFF        | ON      | Build C++ example if set to ON.                     |
| `DUET_BUILD_TEST`        | ON/OFF        | ON      | Build C++ test if set to ON.                        |
| `DUET_BUILD_DEPS`        | ON/OFF        | ON      | Download and build unmet dependencies if set to ON. |

For more compilation options, such as enabling IPCL and network agents, please refer to [PETAce-Solo](https://github.com/tiktok-privacy-innovation/PETAce-Solo) and [PETAce-Network](https://github.com/tiktok-privacy-innovation/PETAce-Network).

We have prepared three code examples to show you how to use PETAce-Duet.
You can find more details in [Duet Examples](example/README.md).

To run c++ examples, execute the following in commands in separate terminal sessions.

```bash
./build/bin/duet_example -p 0
./build/bin/duet_example -p 1
```

You can also use `./build/bin/duet_example -h` to learn more details.

We have also prepared a series of performance test cases for PETAce-Duet.

To run the benchmark, execute the following in commands in separate terminal sessions.

```bash
./build/bin/duet_bench -p 0 --log_path ./duet0.log
./build/bin/duet_bench -p 1 --log_path ./duet1.log
```

You can also use `./build/bin/duet_bench -h` to learn more details.

## Contribution

Please check [Contributing](CONTRIBUTING.md) for more details.

## Code of Conduct

Please check [Code of Conduct](CODE_OF_CONDUCT.md) for more details.

## License

This project is licensed under the [Apache-2.0 License](LICENSE).

## Citing PETAce

To cite PETAce in academic papers, please use the following BibTeX entries.

### Version 0.3.0

```tex
    @misc{petace,
        title = {PETAce (release 0.3.0)},
        howpublished = {\url{https://github.com/tiktok-privacy-innovation/PETAce}},
        month = Jun,
        year = 2024,
        note = {TikTok Pte. Ltd.},
        key = {PETAce}
    }
```
