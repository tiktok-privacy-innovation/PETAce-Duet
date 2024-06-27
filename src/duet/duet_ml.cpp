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

void Duet::sigmoid(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, ArithMatrix& out) const {
    ArithMatrix x_2;
    ArithMatrix x_3;
    ArithMatrix ax_3;

    elementwise_mul(net, in, in, x_2);
    elementwise_mul(net, in, x_2, x_3);
    elementwise_mul(kSigmoidParmas[1], in, out);
    elementwise_mul(kSigmoidParmas[0], x_3, ax_3);
    add(out, kSigmoidParmas[2], out);
    add(out, ax_3, out);
}

void Duet::split_by_condition(const std::shared_ptr<network::Network>& net, const BoolMatrix& cond,
        const ArithMatrix& in, ArithMatrix& out0, ArithMatrix& out1) const {
    if (cond.cols() != 1) {
        throw std::invalid_argument("only support the condition matrix in (x, 1) shape");
    }
    PublicMatrixBool tmp;
    reveal(net, cond, tmp);
    std::int64_t count = tmp.matrix().sum();
    std::size_t index_0 = 0;
    std::size_t index_1 = 0;
    out0.resize(in.rows() - count, in.cols());
    out1.resize(count, in.cols());
    for (std::size_t i = 0; i < cond.rows(); ++i) {
        if (tmp(i) == 0) {
            out0.shares().row(index_0) = in.shares().row(i);
            index_0++;
        } else {
            out1.shares().row(index_1) = in.shares().row(i);
            index_1++;
        }
    }
}

}  // namespace duet
}  // namespace petace
