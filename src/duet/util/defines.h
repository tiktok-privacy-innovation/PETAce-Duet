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

#include <emmintrin.h>

#include <vector>

#include "Eigen/Dense"

namespace petace {
namespace duet {

// some declarations
using block = __m128i;

// parameters for cheetah secure comparison
// note: when setting up these parameters, we need to guarantee that the number of the blocks is a power of 2.
const std::size_t kBlockBitLength = 1;
const std::size_t kOTSize = std::size_t(1) << kBlockBitLength;

using PRGSeed = block;
using OTLabel = block;
using OTChoice = std::int8_t;
using GGMTreeNode = block;
using GGMTreePuncChoice = std::int8_t;
using Buffer = std::vector<std::int8_t>;

template <typename T>
using PlainMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

class BoolMatrix {
public:
    BoolMatrix() {
    }

    BoolMatrix(std::size_t x_size, std::size_t y_size) {
        resize(x_size, y_size);
    }

    PlainMatrix<std::int64_t> shares;

    std::size_t rows() const {
        return shares.rows();
    }
    std::size_t cols() const {
        return shares.cols();
    }
    std::size_t size() const {
        return shares.size();
    }

    void resize(std::size_t x_size, std::size_t y_size) {
        shares.resize(x_size, y_size);
    }

    BoolMatrix operator^(const BoolMatrix& b) const {
        BoolMatrix c;
        c.shares.resize(b.rows(), b.cols());
        for (std::size_t i = 0; i < b.size(); i++) {
            c.shares(i) = shares(i) ^ b.shares(i);
        }
        return c;
    }

    BoolMatrix operator&(const BoolMatrix& b) const {
        BoolMatrix c;
        c.shares.resize(b.rows(), b.cols());
        for (std::size_t i = 0; i < b.size(); i++) {
            c.shares(i) = shares(i) & b.shares(i);
        }
        return c;
    }
};

class ArithMatrix {
public:
    ArithMatrix() {
    }

    ArithMatrix(std::size_t x_size, std::size_t y_size) {
        resize(x_size, y_size);
    }

    PlainMatrix<std::int64_t> shares;

    std::size_t rows() const {
        return shares.rows();
    }
    std::size_t cols() const {
        return shares.cols();
    }
    std::size_t size() const {
        return shares.size();
    }

    void resize(std::size_t x_size, std::size_t y_size) {
        shares.resize(x_size, y_size);
    }

    ArithMatrix operator+(const ArithMatrix& b) const {
        ArithMatrix c;
        c.shares.resize(b.rows(), b.cols());
        c.shares = shares + b.shares;
        return c;
    }

    ArithMatrix operator-(const ArithMatrix& b) const {
        ArithMatrix c;
        c.shares.resize(b.rows(), b.cols());
        c.shares = shares - b.shares;
        return c;
    }

    ArithMatrix operator*(const ArithMatrix& b) const {
        ArithMatrix c;
        c.shares.resize(b.rows(), b.cols());
        c.shares = shares.cwiseProduct(b.shares);
        return c;
    }
};

}  // namespace duet
}  // namespace petace
