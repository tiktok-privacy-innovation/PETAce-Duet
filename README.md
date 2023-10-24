# PETAce-Duet

## Introduction

PETAce-Duet is a collection of general-purpose two-party secure computing protocols.
It is one of the crucial building blocks of the framework [PETAce](https://github.com/tiktok-privacy-innovation/PETAce).

As implied by the name "Duet", PETAce-Duet implements various two-party MPC protocols as follows:

- [ABY](https://www.ndss-symposium.org/wp-content/uploads/2017/09/08_2_1.pdf) protocols: arithmetic sharing, boolean sharing, share conversion, and operations over shares (e.g., addition and multiplication over arithmetic shares and bit-wise AND operation over boolean shares).

- Secure comparison protocols from ABY and [Cheetah](https://www.usenix.org/system/files/sec22-huang-zhicong.pdf).

- [Secret shared shuffle](https://link.springer.com/chapter/10.1007/978-3-030-64840-4_12). The permutation network will be available in a future release.

- Homomorphic Encryption(HE) Capability: Duet implements the Paillier encryption to support homomorphic encryption(HE) for two-party MPC, the underlying Paillier functionalities are supported through Intel Paillier Cryptosystem Library(IPCL).
We also support the conversion between HE ciphers and arithmetic sharings, as a result, we achieve an "ABH" framework where arithmetic sharings, boolean sharings, and HE sharings can be converted to each other efficiently.
HE provides an efficient alternative for problems such as secure shuffling and linear matrix operations.

The core functionalities of Duet are written in C++ to provide the best performance.
Meanwhile, as one of the significant updates in version v0.2.0, we provide a Python interface "SecureNumpy" so that users can directly write MPC programs with convenient Python interfaces.
Similar to other MPC frameworks, we build the Duet virtual machine to connect Python APIs with C++ functionalities.
Register-based Duet virtual machines translate Python codes to their corresponding bytecode, and then execute the bytecode instructions with C++ functions.
We provide an example in duet_example.py to demonstrate how to build your MPC applications with the Python APIs.

## Requirements

| System | Toolchain                                             |
|--------|-------------------------------------------------------|
| Linux  | Clang++ (>= 5.0) or GNU G++ (>= 5.5), CMake (>= 3.15) |

| Required dependency                                                            | Tested version | Use                               |
|--------------------------------------------------------------------------------|----------------|-----------------------------------|
| [Eigen](https://gitlab.com/libeigen/eigen)                                     | 3.4.0          | Matrix and vector templates       |
| [Intel Paillier Cryptosystem Library (IPCL)](https://github.com/intel/pailliercryptolib)                                     | 2.0.0          | Paillier Encryption       |
| [PETAce-Solo](https://github.com/tiktok-privacy-innovation/PETAce-Solo)       | 0.2.0          | Cryptography primitives           |
| [PETAce-Network](https://github.com/tiktok-privacy-innovation/PETAce-Network) | 0.2.0          | Network communication protocols   |
| [PETAce-Verse](https://github.com/tiktok-privacy-innovation/PETAce-Verse)     | 0.2.0          | Primitive cryptographic protocols |

| Optional dependency                                    | Tested version | Use                    |
|--------------------------------------------------------|----------------|------------------------|
| [GoogleTest](https://github.com/google/googletest)     | 1.12.1         | For running tests      |
| [GoogleBenchmark](https://github.com/google/benchmark) | 1.6.1          | For running benchmarks |
| [TCLAP](https://github.com/mirror/tclap)     | 1.22.2         | For running examples      |

## Building PETAce-Duet

PETAce-Duet depends on IPCL, so you should follow the commands below to install IPCL before building PETAce-Duet:

First, install NASM with `apt install nasm`, or build and install NASM from source code downloaded from [NASM's official page](https://www.nasm.us/). Assume that you are working in the root directory of NASM source code.
```shell
./configure
make -j
make install
```

Second, build and install IPCL using the following scripts.
Assume that IPCL is cloned into the directory `${IPCL}` and will be installed to the directory `${IPCL_INSTALL_DIR}`.
```shell
cmake -B ${IPCL}/build -S ${IPCL} -DCMAKE_INSTALL_PREFIX=${IPCL_INSTALL_DIR} -DCMAKE_BUILD_TYPE=Release -DIPCL_TEST=OFF -DIPCL_BENCHMARK=OFF
cmake --build ${IPCL}/build -j
cmake --build ${IPCL}/build --target install
```

Finnally, you should set the environment varialble `IPCL_DIR` to be the directory where ipcl is installed. Refer to [this link](https://github.com/intel/pailliercryptolib/blob/development/example/README.md) for more information.

Now we can build PETAce-Duet with the instructions below. We assume that all commands presented below are executed in the root directory of PETAce-Duet.

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
### Building python
To build PETAce-Duet python package:

```bash
cmake -S . -B build -DDUET_BUILD_PYTHON=ON
cmake --build build
cd build
make wheel
cd python/wheel
pip3 install petace-0.2.0-py3-none-any.whl
```

We have prepared three code examples to show you how to use PETAce-Duet.
You can find more details in [Duet Examples](example/README.md).

To run c++ examples, execute the following in commands in separate terminal sessions.

```bash
./build/bin/duet_example -p 0
./build/bin/duet_example -p 1
```

To run python examples, execute the following in commands in separate terminal sessions.

```bash
python3 ./example/duet_example.py -p 0
python3 ./example/duet_example.py -p 1
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

### Version 0.2.0

```tex
    @misc{petace,
        title = {PETAce (release 0.2.0)},
        howpublished = {\url{https://github.com/tiktok-privacy-innovation/PETAce}},
        month = Oct,
        year = 2023,
        note = {TikTok Pte. Ltd.},
        key = {PETAce}
    }
```
