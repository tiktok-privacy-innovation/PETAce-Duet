// Copyright 2023 TikTok Pte. Ltd.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <cmath>

#include <algorithm>

#include "gmp.h"
#include "gmpxx.h"

#include "solo/ahe_paillier.h"
#include "solo/prng.h"

#include "duet/util/consts.h"
#include "duet/util/prng.h"
namespace petace {
namespace duet {

inline block read_block_from_dev_urandom() {
    block ret;
    solo::PRNG::get_random_byte_array(sizeof(ret), reinterpret_cast<solo::Byte*>(&ret));
    return ret;
}

inline std::int64_t double_to_fixed(double input) {
    return static_cast<std::int64_t>(input * ((std::int64_t)1 << kFixedPointPrecision));
}

inline double fixed_to_double(std::int64_t input) {
    return static_cast<double>(input) / ((std::int64_t)1 << kFixedPointPrecision);
}

inline std::size_t ceil_log2(std::size_t in) {
    return static_cast<std::size_t>(std::ceil(std::log2(in)));
}

inline mpz_class get_random_mpz(std::shared_ptr<PRNG> prng, std::size_t bits) {
    mpz_class out = 0;
    std::size_t byte_count = (bits + 7 / 8);
    std::vector<solo::Byte> in(byte_count);
    do {
        prng->get_unique_rand_gen()->generate(byte_count, in.data());
        mpz_import(out.get_mpz_t(), in.size(), 1, 1, 0, 0, in.data());
        out >>= (byte_count * 8 - bits);
    } while (out == 0 || mpz_sizeinbase(out.get_mpz_t(), 2) != bits);
    return out;
}

inline void mpz_bn_from_bytes(const petace::solo::Byte* in, std::size_t in_byte_count, mpz_class& out) {
    if (in == nullptr) {
        throw std::invalid_argument("in is nullptr");
    }
    mpz_import(out.get_mpz_t(), in_byte_count, -1, sizeof(petace::solo::Byte), -1, 0, in);
}

inline void mpz_bn_to_bytes(const mpz_class& in, petace::solo::Byte* out, std::size_t out_byte_count) {
    if (out == nullptr) {
        throw std::invalid_argument("out is nullptr");
    }
    std::size_t length = (mpz_sizeinbase(in.get_mpz_t(), 2) + 7) / 8;
    mpz_export(out, nullptr, -1, sizeof(petace::solo::Byte), -1, 0, in.get_mpz_t());
    if (length < out_byte_count) {
        std::fill_n(out + length, out_byte_count - length, petace::solo::Byte('\x00'));
    }
}

}  // namespace duet
}  // namespace petace
