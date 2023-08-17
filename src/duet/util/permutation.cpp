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

#include "duet/util/permutation.h"

#include <algorithm>

#include "duet/util/common.h"
#include "duet/util/prng.h"

namespace petace {
namespace duet {

Permutation::Permutation(const std::vector<std::size_t>& data) {
    std::vector<std::size_t> sort_tmp = data;
    std::sort(sort_tmp.begin(), sort_tmp.end());
    for (std::size_t i = 0; i < sort_tmp.size(); ++i) {
        if (sort_tmp[i] != i) {
            throw std::invalid_argument("index should like 0, 1, 2, ..., n-1");
        }
    }
    data_ = data;
}

Permutation::Permutation(std::size_t size) {
    gen_random_permutation(read_block_from_dev_urandom(), size);
}

Permutation::Permutation(const PRGSeed& seed, std::size_t size) {
    gen_random_permutation(seed, size);
}

const std::vector<std::size_t>& Permutation::data() const {
    return data_;
}

std::size_t Permutation::size() const {
    return data_.size();
}

std::size_t Permutation::operator[](std::size_t index) const {
    if (index > (size() - 1)) {
        throw std::invalid_argument("index shuld less than permutation size ");
    }
    return data_[index];
}

Permutation Permutation::inverse() const {
    std::vector<std::size_t> out(size());
    for (std::size_t i = 0; i < size(); ++i) {
        out[data_[i]] = i;
    }
    return std::move(Permutation(out));
}

void Permutation::gen_random_permutation(const PRGSeed& seed, std::size_t size) {
    data_.resize(size);
    for (std::size_t i = 0; i < size; ++i) {
        data_[i] = i;
    }

    solo::PRNGFactory prng_factory(solo::PRNGScheme::AES_ECB_CTR);
    std::vector<solo::Byte> bytes_seed(sizeof(block));
    memcpy(bytes_seed.data(), reinterpret_cast<solo::Byte*>(const_cast<block*>(&seed)), sizeof(block));
    auto prng = prng_factory.create(bytes_seed, sizeof(std::size_t) * size);

    // maybe update to other interface
    std::shuffle(data_.begin(), data_.end(), solo::PRNGStandard(prng));
}

ArithMatrix Permutation::permute(const ArithMatrix& in) const {
    ArithMatrix out(in.rows(), in.cols());
    std::size_t rows = in.rows();
    for (std::size_t i = 0; i < rows; ++i) {
        out.shares.row(i) = in.shares.row(data_[i]);
    }
    return out;
}

Permutation Permutation::combine(const Permutation& other) const {
    std::vector<std::size_t> tmp(size() + other.size());
    std::size_t my_size = size();
    for (std::size_t i = 0; i < my_size; ++i) {
        tmp[i] = data_[i];
    }

    for (std::size_t i = 0; i < other.size(); ++i) {
        tmp[i + my_size] = other[i] + my_size;
    }
    return std::move(Permutation(tmp));
}

}  // namespace duet
}  // namespace petace
