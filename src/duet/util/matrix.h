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

#include <stdexcept>
#include <vector>

#include "solo/ahe_paillier.h"

#include "duet/util/defines.h"

namespace petace {
namespace duet {

template <typename T>
class MatrixTypeBase {
public:
    MatrixTypeBase() = default;

    ~MatrixTypeBase() = default;

    MatrixTypeBase(std::size_t rows, std::size_t cols) {
        resize(rows, cols);
    }

    std::size_t rows() const {
        return static_cast<std::size_t>(matrix_.rows());
    }

    std::size_t cols() const {
        return static_cast<std::size_t>(matrix_.cols());
    }

    std::size_t size() const {
        return static_cast<std::size_t>(matrix_.size());
    }

    void resize(std::size_t rows, std::size_t cols) {
        matrix_.resize(static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(cols));
    }

    std::vector<std::size_t> shape() {
        std::vector<std::size_t> ret;
        ret.push_back(rows());
        ret.push_back(cols());
        return ret;
    }

    const Matrix<T>& matrix() const {
        return matrix_;
    }

    Matrix<T>& matrix() {
        return matrix_;
    }

    T& operator()(std::size_t index) {
        return matrix_(static_cast<Eigen::Index>(index));
    }

    const T& operator()(std::size_t index) const {
        return matrix_(static_cast<Eigen::Index>(index));
    }

    T& operator()(std::size_t row, std::size_t col) {
        return matrix_(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
    }

    const T& operator()(std::size_t row, std::size_t col) const {
        return matrix_(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
    }

private:
    Matrix<T> matrix_;
};

template <typename T>
class PublicMatrix : public MatrixTypeBase<T> {
public:
    PublicMatrix() = default;

    ~PublicMatrix() = default;

    PublicMatrix(std::size_t rows, std::size_t cols) : MatrixTypeBase<T>(rows, cols) {
    }
};

using PublicMatrixBool = PublicMatrix<std::int64_t>;

template <typename T>
class PrivateMatrix : public MatrixTypeBase<T> {
public:
    PrivateMatrix() = default;

    ~PrivateMatrix() = default;

    explicit PrivateMatrix(std::size_t party_id) : party_id_(party_id) {
    }

    PrivateMatrix(std::size_t rows, std::size_t cols, std::size_t party_id)
            : MatrixTypeBase<T>(rows, cols), party_id_(party_id) {
    }

    void set_party_id(std::size_t party_id) {
        party_id_ = party_id;
    }

    std::size_t party_id() const {
        return party_id_;
    }

    void index_like(std::size_t rows, std::size_t cols, std::size_t axis, std::size_t party_id) {
        this->resize(rows, cols);
        if (this->party_id() == party_id) {
            if (axis == 0) {
                for (size_t i = 0; i < rows; ++i) {
                    for (size_t j = 0; j < cols; ++j) {
                        (*this)(i * cols + j) = static_cast<T>(i);
                    }
                }
            } else if (axis == 1) {
                for (size_t i = 0; i < rows; ++i) {
                    for (size_t j = 0; j < cols; ++j) {
                        (*this)(i * cols + j) = static_cast<T>(j);
                    }
                }
            }
        }
    }

    void index_like(std::size_t rows, std::size_t cols, std::size_t party_id) {
        this->resize(rows, cols);
        if (this->party_id() == party_id) {
            for (size_t i = 0; i < this->size(); ++i) {
                (*this)(i) = static_cast<T>(i);
            }
        }
    }

private:
    std::size_t party_id_ = 0;
};

using PrivateMatrixBool = PrivateMatrix<std::int64_t>;
using PrivateMatrixInt64 = PrivateMatrix<std::int64_t>;

class BoolMatrix : public MatrixTypeBase<std::int64_t> {
public:
    BoolMatrix() = default;

    ~BoolMatrix() = default;

    BoolMatrix(std::size_t rows, std::size_t cols) : MatrixTypeBase<std::int64_t>(rows, cols) {
    }

    const Matrix<std::int64_t>& shares() const;

    Matrix<std::int64_t>& shares();

    BoolMatrix operator^(const BoolMatrix& b) const;

    BoolMatrix operator&(const BoolMatrix& b) const;
};

class ArithMatrix : public MatrixTypeBase<std::int64_t> {
public:
    ArithMatrix() = default;

    ~ArithMatrix() = default;

    ArithMatrix(std::size_t rows, std::size_t cols) : MatrixTypeBase<std::int64_t>(rows, cols) {
    }

    const Matrix<std::int64_t>& shares() const;

    Matrix<std::int64_t>& shares();

    ArithMatrix operator+(const ArithMatrix& b) const;

    ArithMatrix operator-(const ArithMatrix& b) const;

    ArithMatrix operator*(const ArithMatrix& b) const;
};

template <typename T, typename R>
class HEMatrix {
public:
    HEMatrix() {
    }

    HEMatrix(std::size_t x_size, std::size_t y_size) {
        resize(x_size, y_size);
    }

    HEMatrix(std::size_t x_size, std::size_t y_size, std::size_t party) : party_(party) {
        resize(x_size, y_size);
    }

    std::size_t rows() const {
        return rows_;
    }

    std::size_t cols() const {
        return cols_;
    }

    std::size_t party() const {
        return party_;
    }

    void set_party(std::size_t party) {
        party_ = party;
    }

    std::size_t size() const {
        return rows_ * cols_;
    }

    std::size_t num_ciphers() const {
        return cipher_.slot_count();
    }

    const T& ciphers() const {
        return cipher_;
    }

    T& ciphers() {
        return cipher_;
    }

    R& operator()(std::size_t index) {
        return cipher_[index];
    }

    const R& operator()(std::size_t index) const {
        return cipher_[index];
    }

    R& operator()(std::size_t row, std::size_t col) {
        return cipher_[row * cols_ + col];
    }

    const R& operator()(std::size_t row, std::size_t col) const {
        return cipher_[row * cols_ + col];
    }

    void set_ciphers(T& cipher) {
        cipher_ = cipher;
    }

    void resize(std::size_t x_size, std::size_t y_size) {
        rows_ = x_size;
        cols_ = y_size;
        cipher_.resize(rows_ * cols_);
    }

private:
    T cipher_;
    std::size_t rows_;
    std::size_t cols_;
    std::size_t party_;
};

using PaillierMatrix = HEMatrix<std::vector<solo::ahepaillier::Ciphertext>, solo::ahepaillier::Ciphertext>;

template <typename T>
void matrix_block(
        const T& src, T& dst, std::size_t begin_row, std::size_t begin_col, std::size_t row_num, std::size_t col_num) {
    dst.matrix() = src.matrix().block(static_cast<Eigen::Index>(begin_row), static_cast<Eigen::Index>(begin_col),
            static_cast<Eigen::Index>(row_num), static_cast<Eigen::Index>(col_num));
}

template <typename T>
void matrix_block(const PrivateMatrix<T>& src, PrivateMatrix<T>& dst, std::size_t begin_row, std::size_t begin_col,
        std::size_t row_num, std::size_t col_num, std::size_t party_id) {
    if (src.party_id() == party_id) {
        dst.matrix() = src.matrix().block(static_cast<Eigen::Index>(begin_row), static_cast<Eigen::Index>(begin_col),
                static_cast<Eigen::Index>(row_num), static_cast<Eigen::Index>(col_num));
        dst.set_party_id(src.party_id());
    }
}

template <typename T>
void vstack(const T& src_0, const T& src_1, T& dst) {
    if (src_0.cols() != src_1.cols()) {
        throw std::invalid_argument("not support broadcast.");
    }
    dst.resize(src_0.rows() + src_1.rows(), src_0.cols());
    dst.matrix() << src_0.matrix(), src_1.matrix();
}

template <typename T>
void vstack(const PrivateMatrix<T>& src_0, const PrivateMatrix<T>& src_1, PrivateMatrix<T>& dst, std::size_t party_id) {
    if (src_0.party_id() != src_1.party_id()) {
        throw std::invalid_argument("two private matrix must have same party.");
    }

    if (src_0.party_id() == party_id) {
        if (src_0.cols() != src_1.cols()) {
            throw std::invalid_argument("not support broadcast.");
        }
        dst.resize(src_0.rows() + src_1.rows(), src_0.cols());
        dst.matrix() << src_0.matrix(), src_1.matrix();
    }
}

template <typename T>
void hstack(const T& src_0, const T& src_1, T& dst) {
    if (src_0.rows() != src_1.rows()) {
        throw std::invalid_argument("not support broadcast.");
    }
    dst.resize(src_0.rows(), src_0.cols() + src_1.cols());
    dst.matrix() << src_0.matrix(), src_1.matrix();
}

template <typename T>
void hstack(const PrivateMatrix<T>& src_0, const PrivateMatrix<T>& src_1, PrivateMatrix<T>& dst, std::size_t party_id) {
    if (src_0.party_id() != src_1.party_id()) {
        throw std::invalid_argument("two private matrix must have same party.");
    }

    if (src_0.party_id() == party_id) {
        if (src_0.rows() != src_1.rows()) {
            throw std::invalid_argument("not support broadcast.");
        }
        dst.resize(src_0.rows(), src_0.cols() + src_1.cols());
        dst.matrix() << src_0.matrix(), src_1.matrix();
    }
}

}  // namespace duet
}  // namespace petace
