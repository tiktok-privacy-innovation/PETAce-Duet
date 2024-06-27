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

#include <stdexcept>
#include <vector>

#include "duet/duet.h"
#include "duet/util/common.h"
#include "duet/util/consts.h"
#include "duet/util/secret_shared_shuffle.h"

namespace petace {
namespace duet {

void Duet::groupby_sum(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, const ArithMatrix& encoding,
        ArithMatrix& out) const {
    if ((encoding.size() <= 0) || (in.size() <= 0)) {
        throw std::invalid_argument("the encoding size or data size can not less 0");
    }

    std::size_t row = encoding.rows();
    std::size_t col = encoding.cols();

    std::size_t sum_size = in.cols();
    ArithMatrix mul_result;
    ArithMatrix sum_result;
    out.resize(sum_size, col);
    for (std::size_t j = 0; j < sum_size; j++) {
        std::vector<std::size_t> ind(col, j);
        ArithMatrix sum_col(row, col);
        sum_col.shares() = in.shares()(Eigen::all, ind);
        elementwise_mul(net, sum_col, encoding, mul_result);
        sum(mul_result, sum_result);
        out.shares().row(j) = sum_result.shares().row(0);
    }

    return;
}

void Duet::groupby_max(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, const ArithMatrix& encoding,
        ArithMatrix& out) const {
    if ((encoding.size() <= 0) || (in.size() <= 0)) {
        throw std::invalid_argument("the encoding size or data size can not less 0");
    }

    std::size_t row = encoding.rows();
    std::size_t col = encoding.cols();

    std::size_t max_size = in.cols();
    ArithMatrix mul_result(row, col);
    ArithMatrix min_mul_result(row, col);
    ArithMatrix max_result, max_index;

    ArithMatrix min_value(row, col);
    if (party_id_ == 0) {
        min_value.shares().setConstant(0);
    } else {
        min_value.shares().setConstant(kMinValue);
    }

    ArithMatrix one_matrix(row, col);
    if (party_id_ == 0) {
        one_matrix.shares().setZero();
    } else {
        one_matrix.shares().setConstant(1UL << (kFixedPointPrecision));
    }

    ArithMatrix one_sub_encoding(row, col);
    sub(one_matrix, encoding, one_sub_encoding);

    out.resize(max_size, col);
    for (std::size_t j = 0; j < max_size; j++) {
        std::vector<std::size_t> ind(col, j);
        ArithMatrix max_col(row, col);
        max_col.shares() = in.shares()(Eigen::all, ind);
        elementwise_mul(net, max_col, encoding, mul_result);
        elementwise_mul(net, min_value, one_sub_encoding, min_mul_result);
        add(mul_result, min_mul_result, mul_result);
        argmax_and_max(net, mul_result, max_index, max_result);
        out.shares().row(j) = max_result.shares().row(0);
    }

    return;
}

void Duet::groupby_min(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, const ArithMatrix& encoding,
        ArithMatrix& out) const {
    if ((encoding.size() <= 0) || (in.size() <= 0)) {
        throw std::invalid_argument("the encoding size or data size can not less 0");
    }

    std::size_t row = encoding.rows();
    std::size_t col = encoding.cols();

    std::size_t min_size = in.cols();
    ArithMatrix mul_result(row, col);
    ArithMatrix max_mul_result(row, col);
    ArithMatrix min_result, min_index;

    ArithMatrix max_value(row, col);
    if (party_id_ == 0) {
        max_value.shares().setConstant(0);
    } else {
        max_value.shares().setConstant(kMaxValue);
    }

    ArithMatrix one_matrix(row, col);
    if (party_id_ == 0) {
        one_matrix.shares().setZero();
    } else {
        one_matrix.shares().setConstant(1UL << (kFixedPointPrecision));
    }

    ArithMatrix one_sub_encoding(row, col);
    sub(one_matrix, encoding, one_sub_encoding);

    out.resize(min_size, col);
    for (std::size_t j = 0; j < min_size; j++) {
        std::vector<std::size_t> ind(col, j);
        ArithMatrix min_col(row, col);
        min_col.shares() = in.shares()(Eigen::all, ind);
        elementwise_mul(net, min_col, encoding, mul_result);
        elementwise_mul(net, max_value, one_sub_encoding, max_mul_result);
        add(mul_result, max_mul_result, mul_result);
        argmin_and_min(net, mul_result, min_index, min_result);
        out.shares().row(j) = min_result.shares().row(0);
    }

    return;
}

void Duet::groupby_count(const ArithMatrix& in, const ArithMatrix& encoding, ArithMatrix& out) const {
    if ((encoding.size() <= 0) || (in.size() <= 0)) {
        throw std::invalid_argument("the encoding size or data size can not less 0");
    }

    std::size_t col = encoding.cols();

    std::size_t count_size = in.cols();
    ArithMatrix count_result;
    sum(encoding, count_result);

    out.resize(count_size, col);
    for (std::size_t j = 0; j < count_size; j++) {
        out.shares().row(j) = count_result.shares().row(0);
    }

    return;
}

void Duet::group_then_sum_by_grouped_count(
        const PublicMatrix<double>& grouped_count, const ArithMatrix& in, ArithMatrix& out) const {
    if (static_cast<std::size_t>(grouped_count.matrix().sum()) != in.rows()) {
        throw std::invalid_argument("matrix size is not match to grouped count.");
    }
    if (in.cols() != 1) {
        throw std::invalid_argument("now only support the matrix which is one cols.");
    }
    out.resize(grouped_count.rows(), in.cols());
    std::size_t index = 0;
    for (std::size_t i = 0; i < grouped_count.rows(); ++i) {
        out(i) = 0;
        for (std::size_t j = 0; j < static_cast<std::size_t>(grouped_count(i)); ++j) {
            out(i) += in(index);
            index++;
        }
    }
}

}  // namespace duet
}  // namespace petace
