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

Duet::Duet(const std::shared_ptr<network::Network>& net, std::size_t party_id) : party_id_(party_id) {
    block data_from_recv;
    block data_to_send = read_block_from_dev_urandom();

    if (party_id_ == 0) {
        send_block(net, &data_to_send, 1);
        recv_block(net, &data_from_recv, 1);
    } else {
        recv_block(net, &data_from_recv, 1);
        send_block(net, &data_to_send, 1);
    }

    // set ot generator
    rand_generator_ = std::make_shared<PRNG>(data_to_send ^ data_from_recv);
    block seed = _mm_set_epi64x(
            rand_generator_->get_unique_rand<std::int64_t>(), rand_generator_->get_unique_rand<std::int64_t>());
    ot_generator_ = std::make_shared<ObliviousTransfer>(party_id_, seed);
    ot_generator_->initialize(net);
    st_generator_ = std::make_shared<STGenerator>(party_id_, ot_generator_);

    // set triple
    rand_bool_triple_generator_ = std::make_shared<BooleanTriple>(net, party_id_, ot_generator_);
    rand_arith_triple_generator_ =
            TripleFactory::get_instance().build(TripleScheme::FHE, net, party_id_, rand_generator_, ot_generator_);

    paillier_engine_ = std::make_shared<Paillier>(party_id_);

    ByteVector pk_byte(paillier_engine_->get_public_key_byte_count());
    ByteVector other_pk_byte(paillier_engine_->get_public_key_byte_count());
    paillier_engine_->serialize_public_key_to_bytes(paillier_engine_->get_pk(), pk_byte.data());
    if (party_id_ == 0) {
        net->send_data(pk_byte.data(), paillier_engine_->get_public_key_byte_count() * sizeof(solo::Byte));
        net->recv_data(other_pk_byte.data(), paillier_engine_->get_public_key_byte_count() * sizeof(solo::Byte));
    } else {
        net->recv_data(other_pk_byte.data(), paillier_engine_->get_public_key_byte_count() * sizeof(solo::Byte));
        net->send_data(pk_byte.data(), paillier_engine_->get_public_key_byte_count() * sizeof(solo::Byte));
    }
    auto pk_other = std::make_shared<petace::solo::ahepaillier::PublicKey>(kPaillierKeySize, true);
    paillier_engine_->deserialize_public_key_from_bytes(other_pk_byte.data(), pk_other);
    paillier_engine_->set_pk_other(pk_other);

    return;
}

void Duet::share(
        const std::shared_ptr<network::Network>& net, const PrivateMatrix<double>& in, ArithMatrix& out) const {
    if ((in.size() == 0) && (in.party_id() == party_id_)) {
        throw std::invalid_argument("share(): matrix size is equal to zero.");
    }

    std::vector<std::size_t> row_and_col{in.rows(), in.cols()};
    sync_shape(net, in.party_id(), row_and_col, row_and_col);

    out.resize(row_and_col[0], row_and_col[1]);
    if (in.party_id() == party_id_) {
        for (std::size_t i = 0; i < out.size(); i++) {
            out(i) = double_to_fixed(in(i)) - rand_generator_->get_common_rand();
        }
    } else {
        for (std::size_t i = 0; i < out.size(); i++) {
            out(i) = rand_generator_->get_common_rand();
        }
    }
    return;
}

void Duet::share(const std::shared_ptr<network::Network>& net, const PrivateMatrixInt64& in, ArithMatrix& out) const {
    if ((in.size() == 0) && (in.party_id() == party_id_)) {
        throw std::invalid_argument("share(): matrix size is equal to zero.");
    }

    std::vector<std::size_t> row_and_col{in.rows(), in.cols()};
    sync_shape(net, in.party_id(), row_and_col, row_and_col);

    out.resize(row_and_col[0], row_and_col[1]);
    if (in.party_id() == party_id_) {
        for (std::size_t i = 0; i < out.size(); i++) {
            out(i) = in(i)-rand_generator_->get_common_rand();
        }
    } else {
        for (std::size_t i = 0; i < out.size(); i++) {
            out(i) = rand_generator_->get_common_rand();
        }
    }
    return;
}

void Duet::share(const std::shared_ptr<network::Network>& net, const PrivateMatrixBool& in, BoolMatrix& out) {
    if ((in.size() == 0) && (in.party_id() == party_id_)) {
        throw std::invalid_argument("share(): matrix size is equal to zero.");
    }

    std::vector<std::size_t> row_and_col{in.rows(), in.cols()};
    sync_shape(net, in.party_id(), row_and_col, row_and_col);

    out.resize(row_and_col[0], row_and_col[1]);
    if (in.party_id() == party_id_) {
        for (std::size_t i = 0; i < out.size(); i++) {
            out(i) = in(i) ^ rand_generator_->get_common_rand();
        }
    } else {
        for (std::size_t i = 0; i < out.size(); i++) {
            out(i) = rand_generator_->get_common_rand();
        }
    }
    return;
}

void Duet::share(const std::shared_ptr<network::Network>& net, const PrivateMatrixBool& in, BoolMatrix& out) const {
    if ((in.size() == 0) && (in.party_id() == party_id_)) {
        throw std::invalid_argument("share(): matrix size is equal to zero.");
    }

    std::size_t row;
    std::size_t col;
    if (in.party_id() == party_id_) {
        row = in.rows();
        col = in.cols();
        net->send_data(&row, sizeof(std::size_t));
        net->send_data(&col, sizeof(std::size_t));
    } else {
        net->recv_data(&row, sizeof(std::size_t));
        net->recv_data(&col, sizeof(std::size_t));
    }

    out.resize(row, col);
    if (in.party_id() == party_id_) {
        for (std::size_t i = 0; i < out.size(); i++) {
            out(i) = in(i) ^ rand_generator_->get_common_rand();
        }
    } else {
        for (std::size_t i = 0; i < out.size(); i++) {
            out(i) = rand_generator_->get_common_rand();
        }
    }
    return;
}

void Duet::reveal(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& in, PrivateMatrix<double>& out) const {
    if (in.size() == 0) {
        throw std::invalid_argument("reveal(): matrix size is equal to zero.");
    }
    ArithMatrix fixed_matrix(in.rows(), in.cols());
    if (party_id_ != out.party_id()) {
        send_matrix(net, &in.shares(), 1);
    } else {
        recv_matrix(net, &fixed_matrix.shares(), 1);
    }

    out.resize(in.rows(), in.cols());

    if (party_id_ == out.party_id()) {
        fixed_matrix = fixed_matrix + in;

        for (std::size_t i = 0; i < in.size(); i++) {
            out(i) = fixed_to_double(fixed_matrix(i));
        }
    }
    return;
}

void Duet::reveal(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& in, PublicMatrix<double>& out) const {
    PrivateMatrix<double> p0(0);
    PrivateMatrix<double> p1(1);
    reveal(net, in, p0);
    reveal(net, in, p1);
    if (party_id_ == 0) {
        out.matrix() = p0.matrix();
    } else {
        out.matrix() = p1.matrix();
    }
    return;
}

void Duet::reveal(const std::shared_ptr<network::Network>& net, const BoolMatrix& in, PrivateMatrixBool& out) const {
    if (in.size() == 0) {
        throw std::invalid_argument("reveal(): matrix size is equal to zero.");
    }
    ArithMatrix fixed_matrix(in.rows(), in.cols());
    if (party_id_ != out.party_id()) {
        send_matrix(net, &in.shares(), 1);
    } else {
        recv_matrix(net, &fixed_matrix.shares(), 1);
    }

    out.resize(in.rows(), in.cols());

    if (party_id_ == out.party_id()) {
        for (std::size_t i = 0; i < in.size(); i++) {
            out(i) = fixed_matrix(i) ^ in(i);
        }
    }
    return;
}

void Duet::reveal(const std::shared_ptr<network::Network>& net, const BoolMatrix& in, PublicMatrixBool& out) const {
    PrivateMatrixBool p0(0);
    PrivateMatrixBool p1(1);
    reveal(net, in, p0);
    reveal(net, in, p1);
    if (party_id_ == 0) {
        out.matrix() = p0.matrix();
    } else {
        out.matrix() = p1.matrix();
    }
    return;
}

void Duet::add(const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.shares() = x.shares() + y.shares();
    return;
}

void Duet::add(const ArithMatrix& x, const PublicMatrix<double>& y, ArithMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), x.cols());
    if (party_id_ == 0) {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z(i) = x(i) + double_to_fixed(y(i));
        }
    } else {
        z.shares() = x.shares();
    }
    return;
}

void Duet::add(const ArithMatrix& x, PublicDouble y, ArithMatrix& z) const {
    if (party_id_ == 0) {
        z.shares() = x.shares().array() + double_to_fixed(y);
    } else {
        z.shares() = x.shares();
    }
    return;
}

void Duet::sub(const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.shares() = x.shares() - y.shares();
    return;
}

void Duet::sub(const ArithMatrix& x, const PublicMatrix<double>& y, ArithMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), x.cols());
    if (party_id_ == 0) {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z(i) = x(i) - double_to_fixed(y(i));
        }
    } else {
        z.shares() = x.shares();
    }
    return;
}

void Duet::sub(const ArithMatrix& x, PublicDouble y, ArithMatrix& z) const {
    if (party_id_ == 0) {
        z.shares() = x.shares().array() - double_to_fixed(y);
    } else {
        z.shares() = x.shares();
    }
    return;
}

void Duet::elementwise_bool_and(
        const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const BoolMatrix& y, BoolMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t row = x.rows();
    std::size_t col = x.cols();
    z.resize(row, col);

    BoolMatrix triple_a(row, col);
    BoolMatrix triple_b(row, col);
    BoolMatrix triple_c(row, col);
    for (std::size_t i = 0; i < x.size(); ++i) {
        auto triple = rand_bool_triple_generator_->get_rand_triple(net, party_id_);
        triple_a(i) = triple[0];
        triple_b(i) = triple[1];
        triple_c(i) = triple[2];
    }

    BoolMatrix e(row, col);
    BoolMatrix f(row, col);
    for (std::size_t i = 0; i < e.size(); ++i) {
        e(i) = x(i) ^ triple_a(i);
        f(i) = y(i) ^ triple_b(i);
    }

    BoolMatrix reveal_e(row, col);
    BoolMatrix reveal_f(row, col);
    if (party_id_ == 0) {
        send_matrix(net, &e.shares(), 1);
        send_matrix(net, &f.shares(), 1);
        recv_matrix(net, &reveal_e.shares(), 1);
        recv_matrix(net, &reveal_f.shares(), 1);
    } else {
        recv_matrix(net, &reveal_e.shares(), 1);
        recv_matrix(net, &reveal_f.shares(), 1);
        send_matrix(net, &e.shares(), 1);
        send_matrix(net, &f.shares(), 1);
    }

    for (std::size_t i = 0; i < reveal_e.size(); ++i) {
        reveal_e(i) = reveal_e(i) ^ e(i);
        reveal_f(i) = reveal_f(i) ^ f(i);
    }

    if (party_id_ == 0) {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z(i) = (reveal_f(i) & triple_a(i)) ^ (reveal_e(i) & triple_b(i)) ^ triple_c(i);
        }
    } else {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z(i) = (reveal_e(i) & reveal_f(i)) ^ (reveal_f(i) & triple_a(i)) ^ (reveal_e(i) & triple_b(i)) ^
                   triple_c(i);
        }
    }
    return;
}

void Duet::elementwise_bool_and(const PublicMatrixBool& x, const BoolMatrix& y, BoolMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), x.cols());
    for (std::size_t i = 0; i < z.size(); ++i) {
        z(i) = y(i) & x(i);
    }
    return;
}

void Duet::elementwise_bool_or(const PublicMatrixBool& x, const BoolMatrix& y, BoolMatrix& z) const {
    if ((x.rows() != y.rows()) || (x.cols() != y.cols())) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), x.cols());

    if (party_id_ == 0) {
        for (std::size_t i = 0; i < x.size(); ++i) {
            z(i) = x(i) | y(i);
        }
    } else {
        for (std::size_t i = 0; i < x.size(); ++i) {
            z(i) = (~x(i)) & y(i);
        }
    }
    return;
}

void Duet::elementwise_bool_or(
        const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const BoolMatrix& y, BoolMatrix& z) const {
    BoolMatrix tmp0;
    BoolMatrix tmp1;
    elementwise_bool_xor(x, y, tmp0);
    elementwise_bool_and(net, x, y, tmp1);
    elementwise_bool_xor(tmp0, tmp1, z);
    return;
}

void Duet::elementwise_bool_xor(const BoolMatrix& x, const BoolMatrix& y, BoolMatrix& z) const {
    if ((x.rows() != y.rows()) || (x.cols() != y.cols())) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), x.cols());
    for (std::size_t i = 0; i < x.size(); ++i) {
        z(i) = x(i) ^ y(i);
    }
    return;
}

void Duet::elementwise_bool_xor(const PublicMatrixBool& x, const BoolMatrix& y, BoolMatrix& z) const {
    if ((x.rows() != y.rows()) || (x.cols() != y.cols())) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), x.cols());
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (party_id_ == 0) {
            z(i) = x(i) ^ y(i);
        } else {
            z(i) = y(i);
        }
    }
    return;
}

void Duet::elementwise_bool_not(const BoolMatrix& x, BoolMatrix& z) const {
    if (x.size() == 0) {
        throw std::invalid_argument("matrix size is zero.");
    }
    z.resize(x.rows(), x.cols());

    if (party_id_ == 0) {
        z.shares() = x.shares();
    } else {
        for (std::size_t i = 0; i < x.size(); ++i) {
            z(i) = ~x(i);
        }
    }
    return;
}

void Duet::elementwise_mul(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y,
        ArithMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t row = x.rows();
    std::size_t col = x.cols();
    z.resize(row, col);

    ArithMatrix triple_a(row, col);
    ArithMatrix triple_b(row, col);
    ArithMatrix triple_c(row, col);
    for (std::size_t i = 0; i < x.size(); ++i) {
        auto triple = rand_arith_triple_generator_->get_rand_triple(net);
        triple_a(i) = triple[0];
        triple_b(i) = triple[1];
        triple_c(i) = triple[2];
    }

    ArithMatrix e(row, col);
    ArithMatrix f(row, col);
    for (std::size_t i = 0; i < e.size(); ++i) {
        e(i) = x(i) - triple_a(i);
        f(i) = y(i) - triple_b(i);
    }

    ArithMatrix reveal_e(row, col);
    ArithMatrix reveal_f(row, col);
    if (party_id_ == 0) {
        send_matrix(net, &e.shares(), 1);
        send_matrix(net, &f.shares(), 1);
        recv_matrix(net, &reveal_e.shares(), 1);
        recv_matrix(net, &reveal_f.shares(), 1);
    } else {
        recv_matrix(net, &reveal_e.shares(), 1);
        recv_matrix(net, &reveal_f.shares(), 1);
        send_matrix(net, &e.shares(), 1);
        send_matrix(net, &f.shares(), 1);
    }

    for (std::size_t i = 0; i < reveal_e.size(); ++i) {
        reveal_e(i) = reveal_e(i) + e(i);
        reveal_f(i) = reveal_f(i) + f(i);
    }

    if (party_id_ == 0) {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z(i) = ((reveal_f(i) * triple_a(i)) + (reveal_e(i) * triple_b(i)) + triple_c(i)) >> kFixedPointPrecision;
            z(i) += 1;
        }
    } else {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z(i) = ((reveal_e(i) * reveal_f(i)) + (reveal_f(i) * triple_a(i)) + (reveal_e(i) * triple_b(i)) +
                           triple_c(i)) >>
                   kFixedPointPrecision;
        }
    }
    return;
}

void Duet::elementwise_mul(const std::shared_ptr<network::Network>& net, const PublicMatrix<std::int64_t>& precision,
        const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t row = x.rows();
    std::size_t col = x.cols();
    z.resize(row, col);

    ArithMatrix triple_a(row, col);
    ArithMatrix triple_b(row, col);
    ArithMatrix triple_c(row, col);
    for (std::size_t i = 0; i < x.size(); ++i) {
        auto triple = rand_arith_triple_generator_->get_rand_triple(net);
        triple_a(i) = triple[0];
        triple_b(i) = triple[1];
        triple_c(i) = triple[2];
    }

    ArithMatrix e(row, col);
    ArithMatrix f(row, col);
    for (std::size_t i = 0; i < e.size(); ++i) {
        e(i) = x(i) - triple_a(i);
        f(i) = y(i) - triple_b(i);
    }

    ArithMatrix reveal_e(row, col);
    ArithMatrix reveal_f(row, col);
    if (party_id_ == 0) {
        send_matrix(net, &e.shares(), 1);
        send_matrix(net, &f.shares(), 1);
        recv_matrix(net, &reveal_e.shares(), 1);
        recv_matrix(net, &reveal_f.shares(), 1);
    } else {
        recv_matrix(net, &reveal_e.shares(), 1);
        recv_matrix(net, &reveal_f.shares(), 1);
        send_matrix(net, &e.shares(), 1);
        send_matrix(net, &f.shares(), 1);
    }

    for (std::size_t i = 0; i < reveal_e.size(); ++i) {
        reveal_e(i) = reveal_e(i) + e(i);
        reveal_f(i) = reveal_f(i) + f(i);
    }

    if (party_id_ == 0) {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z(i) = ((reveal_f(i) * triple_a(i)) + (reveal_e(i) * triple_b(i)) + triple_c(i)) >> precision(i);
        }
    } else {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z(i) = ((reveal_e(i) * reveal_f(i)) + (reveal_f(i) * triple_a(i)) + (reveal_e(i) * triple_b(i)) +
                           triple_c(i)) >>
                   precision(i);
        }
    }
    return;
}

void Duet::elementwise_mul(const PublicMatrix<double>& x, const ArithMatrix& y, ArithMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), x.cols());
    for (std::size_t i = 0; i < z.size(); ++i) {
        z(i) = (y(i) * double_to_fixed(x(i))) >> kFixedPointPrecision;
    }
    return;
}

void Duet::elementwise_mul(PublicDouble x, const ArithMatrix& y, ArithMatrix& z) const {
    z.shares() = y.shares() * double_to_fixed(x);
    for (std::size_t i = 0; i < z.size(); ++i) {
        z(i) = z(i) >> kFixedPointPrecision;
    }
    return;
}

void Duet::elementwise_div(const ArithMatrix& x, const PublicMatrix<double>& y, ArithMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), x.cols());
    PublicMatrix<double> tmp(y.rows(), y.cols());
    for (std::size_t i = 0; i < y.size(); ++i) {
        tmp(i) = 1 / y(i);
    }
    elementwise_mul(tmp, x, z);
    return;
}

void Duet::elementwise_div(const ArithMatrix& x, PublicDouble y, ArithMatrix& z) const {
    elementwise_mul(1.0 / y, x, z);
}

void Duet::elementwise_div(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y,
        ArithMatrix& z) const {
    std::size_t row = x.rows();
    std::size_t col = x.cols();

    ArithMatrix positive_ones(row, col);
    ArithMatrix minus_ones(row, col);
    ArithMatrix sign_matrix(row, col);
    ArithMatrix positive_y(row, col);
    BoolMatrix sign_y(row, col);

    if (party_id_ == 0) {
        positive_ones.shares().setConstant(double_to_fixed(1));
        minus_ones.shares().setConstant(double_to_fixed(-1));
    } else {
        positive_ones.shares().setZero();
        minus_ones.shares().setZero();
    }

    a2b(net, y, sign_y);
    multiplexer(net, sign_y, positive_ones, minus_ones, sign_matrix);
    elementwise_mul(net, sign_matrix, y, positive_y);

    PublicMatrix<int64_t> alpha(row, col);
    alpha.matrix().setZero();

    PublicMatrix<int64_t> shift_const(row, col);
    ArithMatrix sub_matrix(row, col);
    BoolMatrix cmp_matrix(row, col);
    PrivateMatrix<int64_t> cmp_plain_p0(row, col, 0);
    PrivateMatrix<int64_t> cmp_plain_p1(row, col, 1);
    for (std::size_t i = kPowDepth - 1; i--;) {
        shift_const.matrix().setConstant(static_cast<std::uint64_t>(1) << i);
        shift_const.matrix() = shift_const.matrix() + alpha.matrix();
        for (std::size_t j = 0; j < shift_const.size(); j++) {
            shift_const(j) = static_cast<std::uint64_t>(1) << shift_const(j);
        }
        if (party_id_ == 0) {
            sub_matrix.shares() = shift_const.matrix() - positive_y.shares();
        } else {
            sub_matrix.shares() = -positive_y.shares();
        }
        a2b(net, sub_matrix, cmp_matrix);

        reveal(net, cmp_matrix, cmp_plain_p0);
        reveal(net, cmp_matrix, cmp_plain_p1);

        if (party_id_ == 0) {
            alpha.matrix() = alpha.matrix() + (static_cast<std::uint64_t>(1) << i) * cmp_plain_p0.matrix();
        } else {
            alpha.matrix() = alpha.matrix() + (static_cast<std::uint64_t>(1) << i) * cmp_plain_p1.matrix();
        }
    }

    PublicMatrix<int64_t> one_matrix(row, col);
    one_matrix.matrix().setOnes();

    PublicMatrix<int64_t> precision(row, col);
    precision.matrix() = alpha.matrix() + one_matrix.matrix();

    for (std::size_t i = 0; i < one_matrix.size(); i++) {
        one_matrix(i) = static_cast<std::uint64_t>(1) << precision(i);
    }

    ArithMatrix refactor_y(row, col);
    refactor_y = positive_y;

    PublicMatrix<int64_t> two_point_nine(row, col);
    for (std::size_t i = 0; i < two_point_nine.size(); i++) {
        two_point_nine(i) = static_cast<int64_t>(
                kTwoPointNine * (static_cast<double>(static_cast<std::uint64_t>(1) << precision(i))));
    }
    ArithMatrix two_y(row, col);
    two_y.shares() = 2 * refactor_y.shares();

    ArithMatrix w0(row, col);
    ArithMatrix w1(row, col);
    if (party_id_ == 0) {
        w0.shares() = two_point_nine.matrix() - two_y.shares();
    } else {
        w0.shares() = -two_y.shares();
    }

    ArithMatrix b_w0(row, col);
    elementwise_mul(net, precision, refactor_y, w0, b_w0);

    ArithMatrix epsilon0(row, col);
    ArithMatrix epsilon1(row, col);
    if (party_id_ == 0) {
        epsilon0.shares() = one_matrix.matrix() - b_w0.shares();
    } else {
        epsilon0.shares() = -b_w0.shares();
    }

    elementwise_mul(net, precision, epsilon0, epsilon0, epsilon1);

    ArithMatrix one_add_epsilon0(row, col);
    ArithMatrix one_add_epsilon1(row, col);
    if (party_id_ == 0) {
        one_add_epsilon0.shares() = one_matrix.matrix() + epsilon0.shares();
        one_add_epsilon1.shares() = one_matrix.matrix() + epsilon1.shares();
    } else {
        one_add_epsilon0.shares() = epsilon0.shares();
        one_add_epsilon1.shares() = epsilon1.shares();
    }

    elementwise_mul(net, precision, one_add_epsilon0, one_add_epsilon1, epsilon0);
    elementwise_mul(net, precision, epsilon0, w0, epsilon1);

    PublicMatrix<std::int64_t> refactor_precision(row, col);
    one_matrix.matrix().setOnes();
    refactor_precision.matrix() = 2 * precision.matrix() - kFixedPointPrecision * one_matrix.matrix();

    elementwise_mul(net, refactor_precision, epsilon1, sign_matrix, epsilon0);

    elementwise_mul(net, epsilon0, x, z);

    return;
}

void Duet::mat_mul(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y,
        ArithMatrix& z) const {
    if (x.cols() != y.rows()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), y.cols());
    ArithMatrix row_expand;
    ArithMatrix row_expand_mul;
    ArithMatrix row_expand_mul_sum;
    for (std::size_t i = 0; i < x.rows(); i++) {
        auto tmp = x.shares().row(i).transpose().replicate(1, y.cols());
        row_expand.shares() = tmp;
        elementwise_mul(net, row_expand, y, row_expand_mul);
        sum(row_expand_mul, row_expand_mul_sum);
        z.shares().row(i) = row_expand_mul_sum.shares().row(0);
    }
    return;
}

void Duet::mat_mul(const PublicMatrix<double>& x, const ArithMatrix& y, ArithMatrix& z) const {
    if (x.cols() != y.rows()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), y.cols());
    PublicMatrix<int64_t> public_matrix(x.rows(), x.cols());
    for (std::size_t i = 0; i < x.size(); ++i) {
        public_matrix(i) = (double_to_fixed(x(i)));
    }
    z.shares() = public_matrix.matrix() * y.shares();
    for (std::size_t i = 0; i < z.size(); i++) {
        z(i) = z(i) >> kFixedPointPrecision;
    }
    return;
}

void Duet::mat_mul(const ArithMatrix& x, const PublicMatrix<double>& y, ArithMatrix& z) const {
    if (x.cols() != y.rows()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), y.cols());
    PublicMatrix<int64_t> public_matrix(y.rows(), y.cols());
    for (std::size_t i = 0; i < y.size(); ++i) {
        public_matrix(i) = (double_to_fixed(y(i)));
    }
    z.shares() = x.shares() * public_matrix.matrix();
    for (std::size_t i = 0; i < z.size(); i++) {
        z(i) = z(i) >> kFixedPointPrecision;
    }
    return;
}

void Duet::greater(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    less_than_zero(net, y - x, z);
    return;
}

void Duet::greater(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const PublicMatrix<double>& y,
        BoolMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t size = x.size();
    ArithMatrix c(x.rows(), x.cols());
    if (party_id_ == 0) {
        for (std::size_t j = 0; j < size; j++) {
            c(j) = double_to_fixed(y(j)) - x(j);
        }
    } else {
        c.shares() = -x.shares();
    }
    less_than_zero(net, c, z);
    return;
}

void Duet::greater(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, PublicDouble y, BoolMatrix& z) const {
    PublicMatrix<double> tmp(x.rows(), x.cols());
    tmp.matrix().setConstant(y);
    greater(net, x, tmp, z);
    return;
}

void Duet::less(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z) const {
    if (x.size() != static_cast<std::size_t>(y.size())) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    less_than_zero(net, x - y, z);
    return;
}

void Duet::less(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const PublicMatrix<double>& y,
        BoolMatrix& z) const {
    if (x.size() != static_cast<std::size_t>(y.size())) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t size = x.size();
    ArithMatrix c(x.rows(), x.cols());
    if (party_id_ == 0) {
        for (std::size_t j = 0; j < size; j++) {
            c(j) = x(j) - double_to_fixed(y(j));
        }
    } else {
        c.shares() = x.shares();
    }
    less_than_zero(net, c, z);
    return;
}

void Duet::less(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, PublicDouble y, BoolMatrix& z) const {
    PublicMatrix<double> tmp(x.rows(), x.cols());
    tmp.matrix().setConstant(y);
    less(net, x, tmp, z);
    return;
}

void Duet::greater_equal(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    ArithMatrix c(x.rows(), x.cols());
    c.shares() = x.shares() - y.shares();

    less_than_zero(net, c, z);

    for (std::size_t j = 0; j < x.size(); j++) {
        if (party_id_ == 0) {
            z(j) ^= 1;
        }
    }
    return;
}

void Duet::greater_equal(const std::shared_ptr<network::Network>& net, const ArithMatrix& x,
        const PublicMatrix<double>& y, BoolMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    ArithMatrix c(x.rows(), x.cols());
    if (party_id_ == 0) {
        for (std::size_t j = 0; j < x.size(); j++) {
            c(j) = x(j) - double_to_fixed(y(j));
        }
    } else {
        c.shares() = x.shares();
    }

    less_than_zero(net, c, z);

    for (std::size_t j = 0; j < x.size(); j++) {
        if (party_id_ == 0) {
            z(j) ^= 1;
        }
    }

    return;
}

void Duet::equal(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), y.cols());

    // check if all inputs are good
    if (kBlockBitLength < 1 || kBlockBitLength > 8) {
        throw std::invalid_argument("kBlockBitLength has to be from 1 to 8");
    }

    ArithMatrix xx(x.rows(), x.cols());
    if (party_id_ == 0) {
        xx = x - y;
    } else {
        xx = y - x;
    }

    std::size_t matrix_size = xx.size();
    std::size_t num_blocks = (static_cast<std::uint64_t>(1) << kKoggeStonePpaDepth) / kBlockBitLength;
    if (num_blocks * kBlockBitLength != (static_cast<std::uint64_t>(1) << kKoggeStonePpaDepth)) {
        num_blocks += 1;
    }
    std::size_t num_total_blocks = num_blocks * matrix_size;

    // step 1: break the input to digits
    std::vector<std::int8_t> digits(num_total_blocks, 0);
    const std::int8_t mask_block_bit_length = (static_cast<std::int8_t>(1) << kBlockBitLength) - 1;

    for (std::size_t i = 0; i < matrix_size; ++i) {
        for (std::size_t j = 0; j < num_blocks; ++j) {
            std::size_t shft = j * kBlockBitLength;
            digits[i * num_blocks + j] = static_cast<int8_t>((xx(i) >> shft) & mask_block_bit_length);
        }
    }

    std::vector<std::int64_t> leaf_eq;

    // step 2: generate OT messages
    if (party_id_ == 0) {
        std::vector<std::vector<std::int64_t>> s_ot_msg(num_total_blocks);
        std::vector<std::vector<std::int64_t>> t_ot_msg(num_total_blocks);

        // party_0 is the sender
        for (std::size_t i = 0; i < num_total_blocks; ++i) {
            // generate random elements
            leaf_eq.push_back(rand_generator_->get_unique_rand<std::int8_t>() & 0x1);
            for (std::size_t j = 0; j < kOTSize; j++) {
                t_ot_msg[i].push_back((digits[i] == static_cast<std::int8_t>(j)) ^ leaf_eq[i]);
            }
        }
        ot_generator_->get_nch1_ot(net, kOTSize, t_ot_msg);
    } else {
        // party_1 is the receiver
        ot_generator_->get_nch1_ot(net, kOTSize, digits, leaf_eq);
    }

    // step 3: compute the circuit
    std::size_t log_block_number = ceil_log2(num_blocks);
    for (std::size_t i = 1; i < log_block_number + 1; ++i) {
        std::size_t cur_num_blocks = num_blocks / (static_cast<std::uint64_t>(1) << i);
        std::size_t num_AND = matrix_size * cur_num_blocks;

        BoolMatrix temp_and_input_1(1, num_AND);
        BoolMatrix temp_and_input_2(1, num_AND);
        BoolMatrix temp_and_res(1, num_AND);

        for (std::size_t j = 0; j < cur_num_blocks; ++j) {
            for (std::size_t k = 0; k < matrix_size; k++) {
                temp_and_input_1(k * cur_num_blocks + j) = leaf_eq[k * cur_num_blocks * 2 + 2 * j];
                temp_and_input_2(k * cur_num_blocks + j) = leaf_eq[k * cur_num_blocks * 2 + 2 * j + 1];
            }
        }

        pack_and_evaluate_and(net, temp_and_input_1, temp_and_input_2, temp_and_res);
        for (std::size_t j = 0; j < cur_num_blocks; ++j) {
            for (std::size_t k = 0; k < matrix_size; k++) {
                leaf_eq[k * cur_num_blocks + j] = temp_and_res(k * cur_num_blocks + j) & 0x1;
            }
        }
    }
    for (std::size_t i = 0; i < matrix_size; ++i) {
        z(i) = leaf_eq[i] & 0x1;
    }
    return;
}

void Duet::equal(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const PublicMatrix<double>& y,
        BoolMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), y.cols());

    // check if all inputs are good
    if (kBlockBitLength < 1 || kBlockBitLength > 8) {
        throw std::invalid_argument("kBlockBitLength has to be from 1 to 8");
    }

    ArithMatrix xx(x.rows(), x.cols());
    if (party_id_ == 0) {
        xx.shares() = x.shares();
    } else {
        for (std::size_t j = 0; j < x.size(); j++) {
            xx(j) = double_to_fixed(y(j)) - x(j);
        }
    }

    std::size_t matrix_size = xx.size();
    std::size_t num_blocks = (static_cast<std::uint64_t>(1) << kKoggeStonePpaDepth) / kBlockBitLength;
    if (num_blocks * kBlockBitLength != (static_cast<std::uint64_t>(1) << kKoggeStonePpaDepth)) {
        num_blocks += 1;
    }
    std::size_t num_total_blocks = num_blocks * matrix_size;

    // step 1: break the input to digits
    std::vector<std::int8_t> digits(num_total_blocks, 0);
    const std::int8_t mask_block_bit_length = (static_cast<std::int8_t>(1) << kBlockBitLength) - 1;

    for (std::size_t i = 0; i < matrix_size; ++i) {
        for (std::size_t j = 0; j < num_blocks; ++j) {
            std::size_t shft = j * kBlockBitLength;
            digits[i * num_blocks + j] = static_cast<int8_t>((xx(i) >> shft) & mask_block_bit_length);
        }
    }

    std::vector<std::int64_t> leaf_eq;

    // step 2: generate OT messages
    if (party_id_ == 0) {
        std::vector<std::vector<std::int64_t>> s_ot_msg(num_total_blocks);
        std::vector<std::vector<std::int64_t>> t_ot_msg(num_total_blocks);

        // party_0 is the sender
        for (std::size_t i = 0; i < num_total_blocks; ++i) {
            // generate random elements
            leaf_eq.push_back(rand_generator_->get_unique_rand<std::int8_t>() & 0x1);
            for (std::size_t j = 0; j < kOTSize; j++) {
                t_ot_msg[i].push_back((digits[i] == static_cast<std::int8_t>(j)) ^ leaf_eq[i]);
            }
        }
        ot_generator_->get_nch1_ot(net, kOTSize, t_ot_msg);
    } else {
        // party_1 is the receiver
        ot_generator_->get_nch1_ot(net, kOTSize, digits, leaf_eq);
    }

    // step 3: compute the circuit
    std::size_t log_block_number = ceil_log2(num_blocks);
    for (std::size_t i = 1; i < log_block_number + 1; ++i) {
        std::size_t cur_num_blocks = num_blocks / (static_cast<std::uint64_t>(1) << i);
        std::size_t num_AND = matrix_size * cur_num_blocks;

        BoolMatrix temp_and_input_1(1, num_AND);
        BoolMatrix temp_and_input_2(1, num_AND);
        BoolMatrix temp_and_res(1, num_AND);

        for (std::size_t j = 0; j < cur_num_blocks; ++j) {
            for (std::size_t k = 0; k < matrix_size; k++) {
                temp_and_input_1(k * cur_num_blocks + j) = leaf_eq[k * cur_num_blocks * 2 + 2 * j];
                temp_and_input_2(k * cur_num_blocks + j) = leaf_eq[k * cur_num_blocks * 2 + 2 * j + 1];
            }
        }

        pack_and_evaluate_and(net, temp_and_input_1, temp_and_input_2, temp_and_res);
        for (std::size_t j = 0; j < cur_num_blocks; ++j) {
            for (std::size_t k = 0; k < matrix_size; k++) {
                leaf_eq[k * cur_num_blocks + j] = temp_and_res(k * cur_num_blocks + j) & 0x1;
            }
        }
    }
    for (std::size_t i = 0; i < matrix_size; ++i) {
        z(i) = leaf_eq[i] & 0x1;
    }
    return;
}

void Duet::not_equal(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z) const {
    equal(net, x, y, z);
    for (std::size_t j = 0; j < z.size(); j++) {
        if (party_id_ == 0) {
            z(j) ^= 1;
        }
    }
    return;
}

void Duet::multiplexer(
        const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const ArithMatrix& y, ArithMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t size = x.size();
    std::size_t row = x.rows();
    std::size_t col = x.cols();

    ArithMatrix r(row, col);
    ArithMatrix s0(row, col);
    ArithMatrix s1(row, col);
    for (std::size_t i = 0; i < size; i++) {
        r(i) = rand_generator_->get_unique_rand<std::int64_t>();
    }

    for (std::size_t i = 0; i < size; i++) {
        if (x(i) == 0) {
            s0(i) = -r(i);
            s1(i) = y(i) - r(i);
        } else {
            s0(i) = y(i) - r(i);
            s1(i) = -r(i);
        }
    }

    ArithMatrix y0(row, col);
    ArithMatrix y1(row, col);
    ArithMatrix rb(row, col);
    bool* k = new bool[size];
    if (party_id_ == 0) {
        for (std::size_t i = 0; i < size; i++) {
            std::int8_t choice;
            std::int64_t msg;
            ot_generator_->get_random_ot(net, choice, msg);
            k[i] = choice ^ x(i);
            rb(i) = msg;
        }

        send_bool(net, k, size);
        recv_matrix(net, &y0.shares(), 1);
        recv_matrix(net, &y1.shares(), 1);
        for (std::size_t i = 0; i < size; i++) {
            if (x(i) == 0) {
                z(i) = y0(i) ^ rb(i);
            } else {
                z(i) = y1(i) ^ rb(i);
            }
        }
    } else {
        for (std::size_t i = 0; i < size; i++) {
            std::vector<std::int64_t> msg;
            ot_generator_->get_random_ot(net, msg);
            y0(i) = msg[0];
            y1(i) = msg[1];
        }
        recv_bool(net, k, size);

        for (std::size_t i = 0; i < size; i++) {
            if (k[i] == 0) {
                y0(i) ^= s0(i);
                y1(i) ^= s1(i);

            } else {
                std::int64_t t = s0(i) ^ y1(i);
                y1(i) = s1(i) ^ y0(i);
                y0(i) = t;
            }
        }
        send_matrix(net, &y0.shares(), 1);
        send_matrix(net, &y1.shares(), 1);
    }

    if (party_id_ == 1) {
        for (std::size_t i = 0; i < size; i++) {
            std::int8_t choice;
            std::int64_t msg;
            ot_generator_->get_random_ot(net, choice, msg);
            k[i] = choice ^ x(i);
            rb(i) = msg;
        }
        send_bool(net, k, size);
        recv_matrix(net, &y0.shares(), 1);
        recv_matrix(net, &y1.shares(), 1);
        for (std::size_t i = 0; i < size; i++) {
            if (x(i) == 0) {
                z(i) = y0(i) ^ rb(i);
            } else {
                z(i) = y1(i) ^ rb(i);
            }
        }
    } else {
        for (std::size_t i = 0; i < size; i++) {
            std::vector<std::int64_t> msg;
            ot_generator_->get_random_ot(net, msg);
            y0(i) = msg[0];
            y1(i) = msg[1];
        }
        recv_bool(net, k, size);
        for (std::size_t i = 0; i < size; i++) {
            if (k[i] == 0) {
                y0(i) ^= s0(i);
                y1(i) ^= s1(i);

            } else {
                std::int64_t t = s0(i) ^ y1(i);
                y1(i) = s1(i) ^ y0(i);
                y0(i) = t;
            }
        }
        send_matrix(net, &y0.shares(), 1);
        send_matrix(net, &y1.shares(), 1);
    }

    for (std::size_t i = 0; i < size; i++) {
        z(i) += r(i);
    }
    delete[] k;

    return;
}

void Duet::multiplexer(const std::shared_ptr<network::Network>& net, const BoolMatrix& alpha, const ArithMatrix& x,
        const ArithMatrix& y, ArithMatrix& z) const {
    if ((alpha.size() != x.size()) || (alpha.size() != y.size())) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t size = alpha.size();
    std::size_t row = alpha.rows();
    std::size_t col = alpha.cols();

    z.resize(row, col);

    ArithMatrix r(row, col);
    ArithMatrix s0(row, col);
    ArithMatrix s1(row, col);
    for (std::size_t i = 0; i < size; i++) {
        r(i) = rand_generator_->get_unique_rand<std::int64_t>();
    }

    for (std::size_t i = 0; i < size; i++) {
        if (alpha(i) == 0) {
            s0(i) = x(i) - r(i);
            s1(i) = y(i) - r(i);
        } else {
            s0(i) = y(i) - r(i);
            s1(i) = x(i) - r(i);
        }
    }

    ArithMatrix y0(row, col);
    ArithMatrix y1(row, col);
    ArithMatrix rb(row, col);
    bool* k = new bool[size];
    if (party_id_ == 0) {
        for (std::size_t i = 0; i < size; i++) {
            std::int8_t choice;
            std::int64_t msg;
            ot_generator_->get_random_ot(net, choice, msg);
            k[i] = choice ^ alpha(i);
            rb(i) = msg;
        }

        send_bool(net, k, size);
        recv_matrix(net, &y0.shares(), 1);
        recv_matrix(net, &y1.shares(), 1);
        for (std::size_t i = 0; i < size; i++) {
            if (alpha(i) == 0) {
                z(i) = y0(i) ^ rb(i);
            } else {
                z(i) = y1(i) ^ rb(i);
            }
        }
    } else {
        for (std::size_t i = 0; i < size; i++) {
            std::vector<std::int64_t> msg;
            ot_generator_->get_random_ot(net, msg);
            y0(i) = msg[0];
            y1(i) = msg[1];
        }
        recv_bool(net, k, size);

        for (std::size_t i = 0; i < size; i++) {
            if (k[i] == 0) {
                y0(i) ^= s0(i);
                y1(i) ^= s1(i);

            } else {
                std::int64_t t = s0(i) ^ y1(i);
                y1(i) = s1(i) ^ y0(i);
                y0(i) = t;
            }
        }
        send_matrix(net, &y0.shares(), 1);
        send_matrix(net, &y1.shares(), 1);
    }

    if (party_id_ == 1) {
        for (std::size_t i = 0; i < size; i++) {
            std::int8_t choice;
            std::int64_t msg;
            ot_generator_->get_random_ot(net, choice, msg);
            k[i] = choice ^ alpha(i);
            rb(i) = msg;
        }
        send_bool(net, k, size);
        recv_matrix(net, &y0.shares(), 1);
        recv_matrix(net, &y1.shares(), 1);
        for (std::size_t i = 0; i < size; i++) {
            if (alpha(i) == 0) {
                z(i) = y0(i) ^ rb(i);
            } else {
                z(i) = y1(i) ^ rb(i);
            }
        }
    } else {
        for (std::size_t i = 0; i < size; i++) {
            std::vector<std::int64_t> msg;
            ot_generator_->get_random_ot(net, msg);
            y0(i) = msg[0];
            y1(i) = msg[1];
        }
        recv_bool(net, k, size);
        for (std::size_t i = 0; i < size; i++) {
            if (k[i] == 0) {
                y0(i) ^= s0(i);
                y1(i) ^= s1(i);

            } else {
                std::int64_t t = s0(i) ^ y1(i);
                y1(i) = s1(i) ^ y0(i);
                y0(i) = t;
            }
        }
        send_matrix(net, &y0.shares(), 1);
        send_matrix(net, &y1.shares(), 1);
    }

    for (std::size_t i = 0; i < size; i++) {
        z(i) += r(i);
    }
    delete[] k;

    return;
}

void Duet::kogge_stone_ppa(
        const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const BoolMatrix& y, BoolMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t row = x.rows();
    std::size_t col = x.cols();
    std::size_t depth = kKoggeStonePpaDepth;
    BoolMatrix g1(row, col);
    BoolMatrix p1(row, col);
    BoolMatrix g(row, col);
    BoolMatrix p(row, col);
    std::int64_t keep_masks[kKoggeStonePpaDepth] = {0x0000000000000001, 0x0000000000000003, 0x000000000000000f,
            0x00000000000000ff, 0x000000000000ffff, 0x00000000ffffffff};

    elementwise_bool_and(net, x, y, g);
    for (std::size_t i = 0; i < x.size(); i++) {
        p(i) = x(i) ^ y(i);
    }

    for (std::size_t i = 0; i < depth; i++) {
        std::int64_t shift = static_cast<std::uint64_t>(1) << i;
        for (std::size_t k = 0; k < p.size(); k++) {
            p1(k) = p(k) << shift;
        }
        for (std::size_t k = 0; k < g.size(); k++) {
            g1(k) = g(k) << shift;
        }

        if (party_id_ == 0) {
            for (std::size_t k = 0; k < p.size(); k++) {
                p1(k) ^= keep_masks[i];
            }
        }
        elementwise_bool_and(net, p, g1, g1);
        for (std::size_t k = 0; k < g.size(); k++) {
            g(k) ^= g1(k);
        }
        elementwise_bool_and(net, p, p1, p);
    }

    for (std::size_t k = 0; k < g.size(); k++) {
        g1(k) = g(k) << 1;
    }

    for (std::size_t k = 0; k < g.size(); k++) {
        z(k) = g1(k) ^ x(k) ^ y(k);
    }

    return;
}

void Duet::a2b(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, BoolMatrix& z) const {
    if (x.size() == 0) {
        throw std::invalid_argument("a2b(): matrix size is equal to zero.");
    }
    std::size_t size = x.size();
    std::size_t row = x.rows();
    std::size_t col = x.cols();
    z.resize(row, col);

    BoolMatrix input_0(row, col);
    BoolMatrix input_1(row, col);
    if (party_id_ == 0) {
        for (std::size_t j = 0; j < size; j++) {
            input_0(j) = x(j) ^ rand_generator_->get_common_rand();
            input_1(j) = rand_generator_->get_common_rand();
        }
    } else {
        for (std::size_t j = 0; j < size; j++) {
            input_0(j) = rand_generator_->get_common_rand();
            input_1(j) = x(j) ^ rand_generator_->get_common_rand();
        }
    }
    kogge_stone_ppa(net, input_0, input_1, z);
    for (std::size_t j = 0; j < size; j++) {
        z(j) = (z(j) >> 63) & 0x1;
    }

    return;
}

void Duet::millionaire(const std::shared_ptr<network::Network>& net, Matrix<std::int64_t>& x, BoolMatrix& y,
        std::size_t bit_length) const {
    // check if all inputs are good
    if (kBlockBitLength < 1 || kBlockBitLength > 8) {
        throw std::invalid_argument("kBlockBitLength has to be from 1 to 8");
    }

    std::size_t matrix_size = x.size();
    std::size_t num_blocks = bit_length / kBlockBitLength;
    if (num_blocks * kBlockBitLength != bit_length) {
        num_blocks += 1;
    }
    std::size_t num_total_blocks = num_blocks * matrix_size;

    // step 1: break the input to digits
    std::vector<std::int8_t> digits(num_total_blocks, 0);
    const std::int8_t mask_block_bit_length = (static_cast<std::int8_t>(1) << kBlockBitLength) - 1;

    for (std::size_t i = 0; i < matrix_size; ++i) {
        for (std::size_t j = 0; j < num_blocks; ++j) {
            std::size_t shft = j * kBlockBitLength;
            digits[i * num_blocks + j] = static_cast<int8_t>((x(i) >> shft) & mask_block_bit_length);
        }
    }

    std::vector<std::int64_t> leaf(2 * num_total_blocks);
    std::vector<std::vector<std::int64_t>> ot_msg(2 * num_total_blocks);

    // step 2: generate OT messages
    if (party_id_ == 0) {
        // party_0 is the sender
        for (std::size_t i = 0; i < num_total_blocks; ++i) {
            // generate random elements
            std::int8_t rand_data = rand_generator_->get_unique_rand<std::int8_t>();
            leaf[i] = rand_data & 0x1;
            leaf[num_total_blocks + i] = (rand_data >> 1) & 0x1;
            for (std::size_t j = 0; j < kOTSize; j++) {
                ot_msg[i].push_back((digits[i] < static_cast<std::int8_t>(j)) ^ leaf[i]);
            }
            for (std::size_t j = 0; j < kOTSize; j++) {
                ot_msg[num_total_blocks + i].push_back(
                        (digits[i] == static_cast<std::int8_t>(j)) ^ leaf[num_total_blocks + i]);
            }
        }
        ot_generator_->get_nch1_ot(net, kOTSize, ot_msg);
    } else {
        // party_1 is the receiver
        std::vector<int8_t> double_digits;
        std::copy(digits.begin(), digits.end(), std::back_inserter(double_digits));
        std::copy(digits.begin(), digits.end(), std::back_inserter(double_digits));
        ot_generator_->get_nch1_ot(net, kOTSize, double_digits, leaf);
    }

    // step 3: compute the circuit
    std::size_t log_block_number = ceil_log2(num_blocks);
    for (std::size_t i = 1; i < log_block_number + 1; ++i) {
        std::size_t cur_num_blocks = num_blocks / (static_cast<std::uint64_t>(1) << i);
        std::size_t num_AND = matrix_size * cur_num_blocks;

        BoolMatrix temp_and_input_1(1, 2 * num_AND);
        BoolMatrix temp_and_input_2(1, 2 * num_AND);
        BoolMatrix temp_and_res(1, 2 * num_AND);
        std::vector<std::int64_t> leaf_lt_temp(num_AND, 0);

        for (std::size_t j = 0; j < cur_num_blocks; ++j) {
            for (std::size_t k = 0; k < matrix_size; k++) {
                temp_and_input_1(k * cur_num_blocks + j) = leaf[k * cur_num_blocks * 2 + 2 * j];
                temp_and_input_2(k * cur_num_blocks + j) = leaf[num_total_blocks + k * cur_num_blocks * 2 + 2 * j + 1];
            }
        }

        for (std::size_t j = 0; j < cur_num_blocks; ++j) {
            for (std::size_t k = 0; k < matrix_size; k++) {
                temp_and_input_1(num_AND + k * cur_num_blocks + j) =
                        leaf[num_total_blocks + k * cur_num_blocks * 2 + 2 * j];
                temp_and_input_2(num_AND + k * cur_num_blocks + j) =
                        leaf[num_total_blocks + k * cur_num_blocks * 2 + 2 * j + 1];
            }
        }

        pack_and_evaluate_and(net, temp_and_input_1, temp_and_input_2, temp_and_res);

        for (std::size_t j = 0; j < cur_num_blocks; ++j) {
            for (std::size_t k = 0; k < matrix_size; k++) {
                leaf_lt_temp[k * cur_num_blocks + j] =
                        (leaf[k * cur_num_blocks * 2 + 2 * j + 1] ^ temp_and_res(k * cur_num_blocks + j)) & 0x1;
            }
        }

        leaf.resize(2 * leaf_lt_temp.size());
        leaf.insert(leaf.begin(), leaf_lt_temp.begin(), leaf_lt_temp.end());
        for (std::size_t j = 0; j < cur_num_blocks; ++j) {
            for (std::size_t k = 0; k < matrix_size; k++) {
                leaf[leaf_lt_temp.size() + k * cur_num_blocks + j] =
                        temp_and_res(num_AND + k * cur_num_blocks + j) & 0x1;
            }
        }
        num_total_blocks = num_total_blocks / 2;
    }
    for (std::size_t i = 0; i < matrix_size; ++i) {
        y(i) = leaf[i] & 0x1;
    }
}

void Duet::pack_and_evaluate_and(
        const std::shared_ptr<network::Network>& net, BoolMatrix& x, BoolMatrix& y, BoolMatrix& z) const {
    // step 1: packing
    std::size_t num_elements = x.size();
    std::size_t num_pack = static_cast<std::size_t>(ceil(static_cast<double>(num_elements) / 64));
    BoolMatrix x_packed(1, num_pack);
    BoolMatrix y_packed(1, num_pack);
    BoolMatrix z_packed(1, num_pack);
    x_packed.shares().setZero();
    y_packed.shares().setZero();
    for (size_t i = 0; i < num_pack; i++) {
        for (size_t j = 0; j < 64; j++) {
            if (i * 64 + j < num_elements) {
                x_packed(i) |= (x(i * 64 + j) << j);
                y_packed(i) |= (y(i * 64 + j) << j);
            }
        }
    }
    elementwise_bool_and(net, x_packed, y_packed, z_packed);

    // step 2: unpacking
    for (size_t i = 0; i < num_pack; i++) {
        for (size_t j = 0; j < 64; j++) {
            if (i * 64 + j < num_elements) {
                z(i * 64 + j) = (z_packed(i) >> j) & 0x1;
            }
        }
    }
    return;
}

void Duet::less_than_zero(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, BoolMatrix& y,
        std::size_t bit_length) const {
    std::size_t num_elements = x.size();

    // step 1: split the inputs to be msb and the rest
    Matrix<std::int64_t> msb(1, num_elements);
    Matrix<std::int64_t> rest(1, num_elements);
    for (std::size_t i = 0; i < num_elements; i++) {
        msb(i) = static_cast<std::uint64_t>(x(i)) >> (bit_length - 1);
        rest(i) = x(i) & ((static_cast<std::uint64_t>(1) << (bit_length - 1)) - 1);
    }
    // step 2: compute the carry bit of rest
    BoolMatrix carry(1, num_elements);
    if (party_id_ == 0) {
        // party 0 set the input of millionare to be 2^{l - 1} - 1 - rest
        Matrix<std::int64_t> millionare_input(1, num_elements);
        for (std::size_t i = 0; i < num_elements; i++) {
            millionare_input(i) = (static_cast<std::uint64_t>(1) << (bit_length - 1)) - 1 - rest(i);
        }
        millionaire(net, millionare_input, carry);
    } else {
        // party 1 set the input of millionare to be rest
        millionaire(net, rest, carry);
    }
    // step 3: compute the result
    y.resize(x.rows(), x.cols());
    for (std::size_t i = 0; i < num_elements; i++) {
        y(i) = msb(i) ^ carry(i);
    }
    return;
}

void Duet::sync_shape(const std::shared_ptr<network::Network>& net, const std::size_t party,
        const std::vector<std::size_t>& in, std::vector<std::size_t>& out) const {
    if (party == party_id_) {
        net->send_data(in.data(), sizeof(std::size_t) * 2);
    } else {
        out.resize(2);
        net->recv_data(out.data(), sizeof(std::size_t) * 2);
    }
    return;
}

}  // namespace duet
}  // namespace petace
