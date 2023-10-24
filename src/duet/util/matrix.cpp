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

#include "duet/util/matrix.h"

namespace petace {
namespace duet {

const Matrix<std::int64_t>& BoolMatrix::shares() const {
    return matrix();
}

Matrix<std::int64_t>& BoolMatrix::shares() {
    return matrix();
}

BoolMatrix BoolMatrix::operator^(const BoolMatrix& b) const {
    BoolMatrix c;
    c.resize(b.rows(), b.cols());
    for (std::size_t i = 0; i < b.size(); i++) {
        c(i) = (*this)(i) ^ b(i);
    }
    return c;
}

BoolMatrix BoolMatrix::operator&(const BoolMatrix& b) const {
    BoolMatrix c;
    c.resize(b.rows(), b.cols());
    for (std::size_t i = 0; i < b.size(); i++) {
        c(i) = (*this)(i)&b(i);
    }
    return c;
}

const Matrix<std::int64_t>& ArithMatrix::shares() const {
    return matrix();
}

Matrix<std::int64_t>& ArithMatrix::shares() {
    return matrix();
}

ArithMatrix ArithMatrix::operator+(const ArithMatrix& b) const {
    ArithMatrix c;
    c.resize(b.rows(), b.cols());
    c.shares() = (*this).shares() + b.shares();
    return c;
}

ArithMatrix ArithMatrix::operator-(const ArithMatrix& b) const {
    ArithMatrix c;
    c.resize(b.rows(), b.cols());
    c.shares() = (*this).shares() - b.shares();
    return c;
}

ArithMatrix ArithMatrix::operator*(const ArithMatrix& b) const {
    ArithMatrix c;
    c.resize(b.rows(), b.cols());
    c.shares() = (*this).shares().cwiseProduct(b.shares());
    return c;
}

}  // namespace duet
}  // namespace petace
