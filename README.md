# PETAce-Duet

## Introduction

PETAce-Duet is a collection of general-purpose two-party secure computing protocols.
It is one of the many components in [the framework PETAce](https://github.com/tiktok-privacy-innovation/PETAce).

PETAce-Duet implements the following protocols that involve two parties, as implied by the name "Duet".

- [ABY](https://www.ndss-symposium.org/wp-content/uploads/2017/09/08_2_1.pdf) protocols: additive sharing, binary sharing, share conversion, and operations over shares (e.g., addition and multiplication over additive shares and bit-wise AND operation over binary shares).
-  Secure comparison protocols from ABY and [Cheetah](https://www.usenix.org/system/files/sec22-huang-zhicong.pdf).
- [Secret shared shuffle](https://link.springer.com/chapter/10.1007/978-3-030-64840-4_12). Permutation network will be available in a future release.

## Requirements

| System | Toolchain                                             |
|--------|-------------------------------------------------------|
| Linux  | Clang++ (>= 5.0) or GNU G++ (>= 5.5), CMake (>= 3.15) |

| Required dependency                                                            | Tested version | Use                               |
|--------------------------------------------------------------------------------|----------------|-----------------------------------|
| [Eigen](https://gitlab.com/libeigen/eigen)                                     | 3.4.0          | Matrix and vector templates       |
| [PETAce-Solo](https://github.com/tiktok-privacy-innovation/PETAce-Solo)       | 0.1.0          | Cryptography primitives           |
| [PETAce-Network](https://github.com/tiktok-privacy-innovation/PETAce-Network) | 0.1.0          | Network communication protocols   |
| [PETAce-Verse](https://github.com/tiktok-privacy-innovation/PETAce-Verse)     | 0.1.0          | Primitive cryptographic protocols |

| Optional dependency                                    | Tested version | Use                    |
|--------------------------------------------------------|----------------|------------------------|
| [GoogleTest](https://github.com/google/googletest)     | 1.12.1         | For running tests      |
| [GoogleBenchmark](https://github.com/google/benchmark) | 1.6.1          | For running benchmarks |
| [TCLAP](https://github.com/mirror/tclap)     | 1.22.2         | For running examples      |

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

We have prepared three code examples to show you how to use PETAce-Duet.
You can find more details in [Duet Examples](example/README.md).

To run examples, execute the following in commands in separate terminal sessions.

```bash
./build/bin/duet_example -p 0
./build/bin/duet_example -p 1
```

You can also use `./build/bin/duet_example -h` to learn more details

## Contribution

Please check [Contributing](CONTRIBUTING.md) for more details.

## Code of Conduct

Please check [Code of Conduct](CODE_OF_CONDUCT.md) for more details.

## License

This project is licensed under the [Apache-2.0 License](LICENSE).

## Citing PETAce

To cite PETAce in academic papers, please use the following BibTeX entries.

### Version 0.1.0

```tex
    @misc{petace,
        title = {PETAce (release 0.1.0)},
        howpublished = {\url{https://github.com/tiktok-privacy-innovation/PETAce}},
        month = July,
        year = 2023,
        note = {TikTok Pte. Ltd.},
        key = {PETAce}
    }
```
