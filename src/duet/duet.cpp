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

#include "duet/duet.h"

#include <stdexcept>
#include <vector>

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
    // DoublePrg::fixed_key_aes.set_key(seed);
    ot_generator_ = std::make_shared<OTGenerator>(party_id_, seed);
    ot_generator_->initialize(net);
    st_generator_ = std::make_shared<STGenerator>(party_id_, ot_generator_);

    // set triplet
    rand_bool_triplet_generator_ = std::make_shared<BooleanTriplet>(net, party_id_, ot_generator_);
    rand_arith_triplet_generator_ = std::make_shared<ArithmeticTriplet>(net, party_id_, rand_generator_, ot_generator_);

    paillier_engine_ = std::make_shared<PaillierEngine>(party_id_);

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
    std::shared_ptr<petace::solo::ahepaillier::PublicKey> pk_other = nullptr;
    paillier_engine_->deserialize_public_key_from_bytes(other_pk_byte.data(), pk_other);
    paillier_engine_->set_pk_other(pk_other);

    return;
}

void Duet::share(const std::shared_ptr<network::Network>& net, const PrivateMatrix<double>& in, ArithMatrix& out) {
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
            out(i) = double_to_fixed(in(i)) - rand_generator_->get_common_rand();
        }
    } else {
        for (std::size_t i = 0; i < out.size(); i++) {
            out(i) = rand_generator_->get_common_rand();
        }
    }
    return;
}

void Duet::share(const std::shared_ptr<network::Network>& net, const PrivateMatrixBool& in, ArithMatrix& out) {
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
            out(i) = in(i)-rand_generator_->get_common_rand();
        }
    } else {
        for (std::size_t i = 0; i < out.size(); i++) {
            out(i) = rand_generator_->get_common_rand();
        }
    }
    return;
}

void Duet::reveal(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, PrivateMatrix<double>& out) {
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

void Duet::reveal(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, PublicMatrix<double>& out) {
    PrivateMatrix<double> p0(0);
    PrivateMatrix<double> p1(1);
    reveal(net, in, p0);
    reveal(net, in, p1);
    if (party_id_ == 0) {
        out.matrix() = p0.matrix();
    } else {
        out.matrix() = p1.matrix();
    }
}

void Duet::reveal_bool(const std::shared_ptr<network::Network>& net, const BoolMatrix& in, PrivateMatrixBool& out) {
    if (in.size() == 0) {
        throw std::invalid_argument("reveal_bool(): matrix size is equal to zero.");
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

void Duet::add(const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.shares() = x.shares() + y.shares();
    return;
}

void Duet::sub(const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.shares() = x.shares() - y.shares();
    return;
}

void Duet::elementwise_bool_mul(
        const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const BoolMatrix& y, BoolMatrix& z) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t row = x.rows();
    std::size_t col = x.cols();
    z.resize(row, col);

    BoolMatrix triplet_a(row, col);
    BoolMatrix triplet_b(row, col);
    BoolMatrix triplet_c(row, col);
    for (std::size_t i = 0; i < x.size(); ++i) {
        auto triplet = rand_bool_triplet_generator_->get_rand_triplet(net, party_id_);
        triplet_a(i) = triplet[0];
        triplet_b(i) = triplet[1];
        triplet_c(i) = triplet[2];
    }

    BoolMatrix e(row, col);
    BoolMatrix f(row, col);
    for (std::size_t i = 0; i < e.size(); ++i) {
        e(i) = x(i) ^ triplet_a(i);
        f(i) = y(i) ^ triplet_b(i);
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
            z(i) = (reveal_f(i) & triplet_a(i)) ^ (reveal_e(i) & triplet_b(i)) ^ triplet_c(i);
        }
    } else {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z(i) = (reveal_e(i) & reveal_f(i)) ^ (reveal_f(i) & triplet_a(i)) ^ (reveal_e(i) & triplet_b(i)) ^
                   triplet_c(i);
        }
    }
    return;
}

void Duet::elementwise_bool_or(const PublicMatrixBool& x, const BoolMatrix& y, BoolMatrix& z) {
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

void Duet::elementwise_mul(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t row = x.rows();
    std::size_t col = x.cols();
    z.resize(row, col);

    ArithMatrix triplet_a(row, col);
    ArithMatrix triplet_b(row, col);
    ArithMatrix triplet_c(row, col);
    for (std::size_t i = 0; i < x.size(); ++i) {
        auto triplet = rand_arith_triplet_generator_->get_rand_triplet(net);
        triplet_a(i) = triplet[0];
        triplet_b(i) = triplet[1];
        triplet_c(i) = triplet[2];
    }

    ArithMatrix triplet_a_1(row, col);
    ArithMatrix triplet_b_1(row, col);
    ArithMatrix triplet_c_1(row, col);

    if (party_id_ == 0) {
        send_matrix(net, &triplet_a.shares(), 1);
        send_matrix(net, &triplet_b.shares(), 1);
        send_matrix(net, &triplet_c.shares(), 1);
    } else {
        recv_matrix(net, &triplet_a_1.shares(), 1);
        recv_matrix(net, &triplet_b_1.shares(), 1);
        recv_matrix(net, &triplet_c_1.shares(), 1);
    }

    ArithMatrix e(row, col);
    ArithMatrix f(row, col);
    for (std::size_t i = 0; i < e.size(); ++i) {
        e(i) = x(i) - triplet_a(i);
        f(i) = y(i) - triplet_b(i);
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
            z(i) = ((reveal_f(i) * triplet_a(i)) + (reveal_e(i) * triplet_b(i)) + triplet_c(i)) >> kFixedPointPrecision;
        }
    } else {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z(i) = ((reveal_e(i) * reveal_f(i)) + (reveal_f(i) * triplet_a(i)) + (reveal_e(i) * triplet_b(i)) +
                           triplet_c(i)) >>
                   kFixedPointPrecision;
        }
    }
    return;
}

void Duet::elementwise_mul(const std::shared_ptr<network::Network>& net, const PublicMatrix<std::int64_t>& precision,
        const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t row = x.rows();
    std::size_t col = x.cols();
    z.resize(row, col);

    ArithMatrix triplet_a(row, col);
    ArithMatrix triplet_b(row, col);
    ArithMatrix triplet_c(row, col);
    for (std::size_t i = 0; i < x.size(); ++i) {
        auto triplet = rand_arith_triplet_generator_->get_rand_triplet(net);
        triplet_a(i) = triplet[0];
        triplet_b(i) = triplet[1];
        triplet_c(i) = triplet[2];
    }

    ArithMatrix triplet_a_1(row, col);
    ArithMatrix triplet_b_1(row, col);
    ArithMatrix triplet_c_1(row, col);

    if (party_id_ == 0) {
        send_matrix(net, &triplet_a.shares(), 1);
        send_matrix(net, &triplet_b.shares(), 1);
        send_matrix(net, &triplet_c.shares(), 1);
    } else {
        recv_matrix(net, &triplet_a_1.shares(), 1);
        recv_matrix(net, &triplet_b_1.shares(), 1);
        recv_matrix(net, &triplet_c_1.shares(), 1);
    }

    ArithMatrix e(row, col);
    ArithMatrix f(row, col);
    for (std::size_t i = 0; i < e.size(); ++i) {
        e(i) = x(i) - triplet_a(i);
        f(i) = y(i) - triplet_b(i);
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
            z(i) = ((reveal_f(i) * triplet_a(i)) + (reveal_e(i) * triplet_b(i)) + triplet_c(i)) >> precision(i);
        }
    } else {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z(i) = ((reveal_e(i) * reveal_f(i)) + (reveal_f(i) * triplet_a(i)) + (reveal_e(i) * triplet_b(i)) +
                           triplet_c(i)) >>
                   precision(i);
        }
    }
    return;
}

void Duet::elementwise_mul(const PublicMatrix<double>& x, const ArithMatrix& y, ArithMatrix& z) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), x.cols());
    for (std::size_t i = 0; i < z.size(); ++i) {
        z(i) = (y(i) * double_to_fixed(x(i))) >> kFixedPointPrecision;
    }

    return;
}

void Duet::elementwise_share_div_public(const ArithMatrix& x, const PublicMatrix<double>& y, ArithMatrix& z) {
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

void Duet::scalar_mul(PublicDouble x, const ArithMatrix& y, ArithMatrix& z) {
    z.shares() = y.shares() * double_to_fixed(x);
    for (std::size_t i = 0; i < z.size(); ++i) {
        z(i) = z(i) >> kFixedPointPrecision;
    }
    return;
}

void Duet::share_sub_public_double(const ArithMatrix& x, PublicDouble y, ArithMatrix& z) {
    if (party_id_ == 0) {
        z.shares() = x.shares().array() - double_to_fixed(y);
    } else {
        z.shares() = x.shares();
    }
}

void Duet::elementwise_div(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z) {
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
        shift_const.matrix().setConstant(1UL << i);
        shift_const.matrix() = shift_const.matrix() + alpha.matrix();
        for (std::size_t j = 0; j < shift_const.size(); j++) {
            shift_const(j) = 1UL << shift_const(j);
        }
        if (party_id_ == 0) {
            sub_matrix.shares() = shift_const.matrix() - positive_y.shares();
        } else {
            sub_matrix.shares() = -positive_y.shares();
        }
        a2b(net, sub_matrix, cmp_matrix);

        reveal_bool(net, cmp_matrix, cmp_plain_p0);
        reveal_bool(net, cmp_matrix, cmp_plain_p1);

        if (party_id_ == 0) {
            alpha.matrix() = alpha.matrix() + (1UL << i) * cmp_plain_p0.matrix();
        } else {
            alpha.matrix() = alpha.matrix() + (1UL << i) * cmp_plain_p1.matrix();
        }
    }

    PublicMatrix<int64_t> one_matrix(row, col);
    one_matrix.matrix().setOnes();

    PublicMatrix<int64_t> precision(row, col);
    precision.matrix() = alpha.matrix() + one_matrix.matrix();

    for (std::size_t i = 0; i < one_matrix.size(); i++) {
        one_matrix(i) = 1UL << precision(i);
    }

    ArithMatrix refactor_y(row, col);
    refactor_y = positive_y;

    PublicMatrix<int64_t> two_point_nine(row, col);
    for (std::size_t i = 0; i < two_point_nine.size(); i++) {
        two_point_nine(i) = static_cast<int64_t>(2.9142 * (static_cast<double>(1UL << precision(i))));
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

void Duet::kogge_stone_ppa(
        const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const BoolMatrix& y, BoolMatrix& z) {
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

    elementwise_bool_mul(net, x, y, g);
    for (std::size_t i = 0; i < x.size(); i++) {
        p(i) = x(i) ^ y(i);
    }

    for (std::size_t i = 0; i < depth; i++) {
        std::int64_t shift = 1L << i;
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
        elementwise_bool_mul(net, p, g1, g1);
        for (std::size_t k = 0; k < g.size(); k++) {
            g(k) ^= g1(k);
        }
        elementwise_bool_mul(net, p, p1, p);
    }

    for (std::size_t k = 0; k < g.size(); k++) {
        g1(k) = g(k) << 1;
    }

    for (std::size_t k = 0; k < g.size(); k++) {
        z(k) = g1(k) ^ x(k) ^ y(k);
    }

    return;
}

void Duet::a2b(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, BoolMatrix& z) {
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

void Duet::greater(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    a2b(net, y - x, z);
    return;
}

void Duet::greater(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const PublicMatrix<double>& y,
        BoolMatrix& z) {
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
    a2b(net, c, z);
    return;
}

void Duet::less(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z) {
    if (x.size() != static_cast<std::size_t>(y.size())) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    a2b(net, x - y, z);
    return;
}

void Duet::less(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const PublicMatrix<double>& y,
        BoolMatrix& z) {
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
    a2b(net, c, z);
    return;
}

void Duet::greater_equal(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    ArithMatrix c(x.rows(), x.cols());
    c.shares() = x.shares() - y.shares();

    a2b(net, c, z);

    for (std::size_t j = 0; j < x.size(); j++) {
        if (party_id_ == 0) {
            z(j) ^= 1;
        }
    }

    return;
}

void Duet::equal(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.resize(x.rows(), x.cols());
    BoolMatrix c0(x.rows(), x.cols());
    BoolMatrix c1(x.rows(), x.cols());

    greater_equal(net, x, y, c0);
    greater_equal(net, y, x, c1);

    for (std::size_t j = 0; j < x.size(); j++) {
        z(j) = (c0(j)) ^ (c1(j));
        if (party_id_ == 0) {
            z(j) ^= 1;
        }
    }

    return;
}

void Duet::sum(const ArithMatrix& x, ArithMatrix& z) const {
    if (x.size() == 0) {
        throw std::invalid_argument("sum(): matrix size is equal to zero.");
    }
    z.resize(1, x.cols());
    z.shares().row(0) = x.shares().colwise().sum();
    return;
}

void Duet::millionaire(
        const std::shared_ptr<network::Network>& net, Matrix<std::int64_t>& x, BoolMatrix& y, std::size_t bit_length) {
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
    const std::int8_t mask_block_bit_length = (1 << kBlockBitLength) - 1;

    for (std::size_t i = 0; i < matrix_size; ++i) {
        for (std::size_t j = 0; j < num_blocks; ++j) {
            std::size_t shft = j * kBlockBitLength;
            digits[i * num_blocks + j] = static_cast<int8_t>((x(i) >> shft) & mask_block_bit_length);
        }
    }

    std::vector<std::int64_t> leaf_lt;
    std::vector<std::int64_t> leaf_eq;

    // step 2: generate OT messages
    if (party_id_ == 0) {
        std::vector<std::vector<std::int64_t>> s_ot_msg(num_total_blocks);
        std::vector<std::vector<std::int64_t>> t_ot_msg(num_total_blocks);

        // party_0 is the sender
        for (std::size_t i = 0; i < num_total_blocks; ++i) {
            // generate random elements
            leaf_lt.push_back(rand_generator_->get_unique_rand<std::int8_t>() & 0x1);
            leaf_eq.push_back(rand_generator_->get_unique_rand<std::int8_t>() & 0x1);
            for (std::size_t j = 0; j < kOTSize; j++) {
                s_ot_msg[i].push_back((digits[i] < static_cast<std::int8_t>(j)) ^ leaf_lt[i]);
                t_ot_msg[i].push_back((digits[i] == static_cast<std::int8_t>(j)) ^ leaf_eq[i]);
            }
        }
        ot_generator_->get_nch1_ot(net, kOTSize, s_ot_msg);
        ot_generator_->get_nch1_ot(net, kOTSize, t_ot_msg);
    } else {
        // party_1 is the receiver

        ot_generator_->get_nch1_ot(net, kOTSize, digits, leaf_lt);
        ot_generator_->get_nch1_ot(net, kOTSize, digits, leaf_eq);
    }

    // step 3: compute the circuit
    std::size_t log_block_number = ceil_log2(num_blocks);
    for (std::size_t i = 1; i < log_block_number + 1; ++i) {
        std::size_t cur_num_blocks = num_blocks / (std::size_t(1) << i);
        std::size_t num_AND = matrix_size * cur_num_blocks;

        BoolMatrix temp_and_input_1(1, num_AND);
        BoolMatrix temp_and_input_2(1, num_AND);
        BoolMatrix temp_and_res(1, num_AND);
        std::vector<std::int64_t> leaf_lt_temp(num_AND, 0);

        for (std::size_t j = 0; j < cur_num_blocks; ++j) {
            for (std::size_t k = 0; k < matrix_size; k++) {
                temp_and_input_1(k * cur_num_blocks + j) = leaf_lt[k * cur_num_blocks * 2 + 2 * j];
                temp_and_input_2(k * cur_num_blocks + j) = leaf_eq[k * cur_num_blocks * 2 + 2 * j + 1];
            }
        }

        pack_and_evaluate_and(net, temp_and_input_1, temp_and_input_2, temp_and_res);

        for (std::size_t j = 0; j < cur_num_blocks; ++j) {
            for (std::size_t k = 0; k < matrix_size; k++) {
                leaf_lt_temp[k * cur_num_blocks + j] =
                        (leaf_lt[k * cur_num_blocks * 2 + 2 * j + 1] ^ temp_and_res(k * cur_num_blocks + j)) & 0x1;
                temp_and_input_1(k * cur_num_blocks + j) = leaf_eq[k * cur_num_blocks * 2 + 2 * j];
                temp_and_input_2(k * cur_num_blocks + j) = leaf_eq[k * cur_num_blocks * 2 + 2 * j + 1];
            }
        }
        leaf_lt.assign(leaf_lt_temp.begin(), leaf_lt_temp.end());

        pack_and_evaluate_and(net, temp_and_input_1, temp_and_input_2, temp_and_res);
        for (std::size_t j = 0; j < cur_num_blocks; ++j) {
            for (std::size_t k = 0; k < matrix_size; k++) {
                leaf_eq[k * cur_num_blocks + j] = temp_and_res(k * cur_num_blocks + j) & 0x1;
            }
        }
    }
    for (std::size_t i = 0; i < matrix_size; ++i) {
        y(i) = leaf_lt[i] & 0x1;
    }
}

void Duet::pack_and_evaluate_and(
        const std::shared_ptr<network::Network>& net, BoolMatrix& x, BoolMatrix& y, BoolMatrix& z) {
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
    elementwise_bool_mul(net, x_packed, y_packed, z_packed);

    // step 2: unpacking
    for (size_t i = 0; i < num_pack; i++) {
        for (size_t j = 0; j < 64; j++) {
            if (i * 64 + j < num_elements) {
                z(i * 64 + j) = (z_packed(i) >> j) & 0x1;
            }
        }
    }
}

void Duet::less_than_zero(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, BoolMatrix& y, std::size_t bit_length) {
    std::size_t num_elements = x.size();

    // step 1: split the inputs to be msb and the rest
    Matrix<std::int64_t> msb(1, num_elements);
    Matrix<std::int64_t> rest(1, num_elements);
    for (std::size_t i = 0; i < num_elements; i++) {
        msb(i) = static_cast<std::uint64_t>(x(i)) >> (bit_length - 1);
        rest(i) = x(i) & ((1 << (bit_length - 1)) - 1);
    }
    // step 2: compute the carry bit of rest
    BoolMatrix carry(1, num_elements);
    if (party_id_ == 0) {
        // party 0 set the input of millionare to be 2^{l - 1} - 1 - rest
        Matrix<std::int64_t> millionare_input(1, num_elements);
        for (std::size_t i = 0; i < num_elements; i++) {
            millionare_input(i) = (1 << (bit_length - 1)) - 1 - rest(i);
        }
        millionaire(net, millionare_input, carry);
    } else {
        // party 1 set the input of millionare to be rest
        millionaire(net, rest, carry);
    }
    // step 3: compute the result
    for (std::size_t i = 0; i < num_elements; i++) {
        y(i) = msb(i) ^ carry(i);
    }
}

void Duet::multiplexer(
        const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const ArithMatrix& y, ArithMatrix& z) {
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
        const ArithMatrix& y, ArithMatrix& z) {
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

void Duet::shuffle(const std::shared_ptr<network::Network>& net, const PrivateMatrix<double>& in, ArithMatrix& out) {
    std::size_t row;
    std::size_t col;
    std::size_t matrix_size;
    SecretSharedShuffle ss_shuffle;
    if (in.party_id() == party_id_) {
        row = in.rows();
        col = in.cols();
        matrix_size = in.size();
        net->send_data(&row, sizeof(std::size_t));
        net->send_data(&col, sizeof(std::size_t));
        Matrix<std::int64_t> x(row, col);
        Matrix<std::int64_t> a;
        Matrix<std::int64_t> b;
        Matrix<std::int64_t> x_sub_a;
        for (std::size_t i = 0; i < matrix_size; ++i) {
            x(i) = double_to_fixed(in(i));
        }
        st_generator_->get_st(row, col, net, a, b);
        ss_shuffle.passive_phase_1(x, a, x_sub_a);
        out.shares() = b;
        send_matrix(net, &x_sub_a, 1);

    } else {
        net->recv_data(&row, sizeof(std::size_t));
        net->recv_data(&col, sizeof(std::size_t));
        Matrix<std::int64_t> delta;
        Matrix<std::int64_t> x_sub_a(row, col);
        Matrix<std::int64_t> share;
        Permutation p(row);
        st_generator_->get_st(row, col, net, p, delta);
        recv_matrix(net, &x_sub_a, 1);
        ss_shuffle.active_phase_1(p, x_sub_a, delta, share);
        out.shares() = share;
    }
}

void Duet::shuffle(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, ArithMatrix& out) {
    std::size_t row;
    std::size_t col;
    SecretSharedShuffle ss_shuffle;
    row = in.rows();
    col = in.cols();
    if (party_id_ == 0) {
        Matrix<std::int64_t> a0;
        Matrix<std::int64_t> b0;
        Matrix<std::int64_t> x_sub_a0;
        Matrix<std::int64_t> share0;
        st_generator_->get_st(row, col, net, a0, b0);
        ss_shuffle.passive_phase_1(in.shares(), a0, x_sub_a0);
        share0 = b0;
        send_matrix(net, &x_sub_a0, 1);

        Matrix<std::int64_t> delta1;
        Matrix<std::int64_t> x_sub_a1(row, col);
        Matrix<std::int64_t> share1;
        Permutation p1(row);
        st_generator_->get_st(row, col, net, p1, delta1);
        recv_matrix(net, &x_sub_a1, 1);
        ss_shuffle.active_phase_1(p1, x_sub_a1, delta1, share1);
        out.shares() = share1 + p1.permute(share0);
    } else {
        Matrix<std::int64_t> delta0;
        Matrix<std::int64_t> x_sub_a0(row, col);
        Matrix<std::int64_t> share0;
        Permutation p0(row);
        st_generator_->get_st(row, col, net, p0, delta0);
        recv_matrix(net, &x_sub_a0, 1);
        ss_shuffle.active_phase_1(p0, x_sub_a0, delta0, share0);
        share0 = share0 + p0.permute(in.shares());

        Matrix<std::int64_t> a1;
        Matrix<std::int64_t> b1;
        Matrix<std::int64_t> x_sub_a1;
        Matrix<std::int64_t> share1;
        st_generator_->get_st(row, col, net, a1, b1);
        ss_shuffle.passive_phase_1(share0, a1, x_sub_a1);
        out.shares() = b1;
        send_matrix(net, &x_sub_a1, 1);
    }
}

void Duet::row_major_argmax_and_max(const std::shared_ptr<network::Network>& net, const ArithMatrix& in,
        ArithMatrix& max_index, ArithMatrix& max_value) {
    std::size_t rows = in.rows();
    std::size_t cols = in.cols();
    PrivateMatrix<double> index(0);
    index.index_like(rows, cols, 0, party_id_);
    ArithMatrix index_share;
    share(net, index, index_share);
    ArithMatrix now = in;
    ArithMatrix up;
    ArithMatrix down;
    ArithMatrix rem;
    BoolMatrix cmp;
    ArithMatrix index_now = index_share;
    ArithMatrix index_up;
    ArithMatrix index_down;
    ArithMatrix index_rem;
    ArithMatrix tmp;
    ArithMatrix index_tmp;
    while (now.rows() > 1) {
        std::size_t half_num = now.rows() / 2;
        std::size_t rem_num = now.rows() % 2;
        matrix_block(now, up, 0, 0, half_num, cols);
        matrix_block(now, down, half_num, 0, half_num, cols);
        matrix_block(index_now, index_up, 0, 0, half_num, cols);
        matrix_block(index_now, index_down, half_num, 0, half_num, cols);
        if (rem_num == 1) {
            matrix_block(now, rem, half_num * 2, 0, 1, cols);
            matrix_block(index_now, index_rem, half_num * 2, 0, 1, cols);
        }
        greater_equal(net, up, down, cmp);
        multiplexer(net, cmp, down, up, now);
        multiplexer(net, cmp, index_down, index_up, index_now);

        if (rem_num == 1) {
            vstack(now, rem, tmp);
            now = tmp;
            vstack(index_now, index_rem, index_tmp);
            index_now = index_tmp;
        }
    }
    max_value = now;
    max_index = index_now;
}

void Duet::encrypt(const PrivateMatrix<std::int64_t>& input, PaillierMatrix& output, size_t using_self_pk) {
    paillier_engine_->encrypt(input, output, using_self_pk);
    output.set_party(party_id_ * using_self_pk + (1 - party_id_) * (1 - using_self_pk));
}

void Duet::decrypt(const PaillierMatrix& input, PrivateMatrix<std::int64_t>& output) {
    paillier_engine_->decrypt(input, output);
}

void Duet::encrypt(const PrivateMatrix<double>& input, PaillierMatrix& output, size_t using_self_pk) {
    paillier_engine_->encrypt_double(input, output, using_self_pk);
    output.set_party(party_id_ * using_self_pk + (1 - party_id_) * (1 - using_self_pk));
}

void Duet::decrypt(const PaillierMatrix& input, PrivateMatrix<double>& output) {
    paillier_engine_->decrypt_double(input, output);
}

void Duet::add(const PrivateMatrix<std::int64_t>& x, const PaillierMatrix& y, PaillierMatrix& z) {
    if (static_cast<std::size_t>(x.size()) != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    solo::ahepaillier::Plaintext pt;
    paillier_engine_->encode(x.matrix(), pt);
    paillier_engine_->add(y.ciphers(), pt, z.ciphers());
    z.resize(y.rows(), y.cols());
    z.set_party(y.party());
}

void Duet::add(const PaillierMatrix& x, const PaillierMatrix& y, PaillierMatrix& z) {
    if (static_cast<std::size_t>(x.size()) != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    paillier_engine_->add(x.ciphers(), y.ciphers(), z.ciphers());
    z.resize(y.rows(), y.cols());
    z.set_party(y.party());
}

void Duet::add(const PrivateMatrix<double>& x, const PaillierMatrix& y, PaillierMatrix& z) {
    if (static_cast<std::size_t>(x.size()) != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    solo::ahepaillier::Plaintext pt;
    std::vector<std::uint64_t> vec;
    vec.resize(x.size());
    for (std::size_t i = 0; i < static_cast<std::size_t>(x.size()); ++i) {
        vec[i] = double_to_fixed(x(i));
    }
    paillier_engine_->encode(vec, pt);
    paillier_engine_->add(y.ciphers(), pt, z.ciphers());
    z.resize(y.rows(), y.cols());
    z.set_party(y.party());
}

void Duet::mul(const PrivateMatrix<std::int64_t>& x, const PaillierMatrix& y, PaillierMatrix& z) {
    if (static_cast<std::size_t>(x.size()) != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    solo::ahepaillier::Plaintext pt;
    paillier_engine_->encode(x.matrix(), pt);
    paillier_engine_->mul(y.ciphers(), pt, z.ciphers());
    z.resize(y.rows(), y.cols());
    z.set_party(y.party());
}

void Duet::h2a(const std::shared_ptr<network::Network>& net, const PaillierMatrix& in, ArithMatrix& out) {
    std::size_t row;
    std::size_t col;
    std::size_t party_he_stored_owner = 1 - in.party();
    std::vector<solo::ahepaillier::BigNum> random_r_buffer;
    PaillierMatrix x_plus_r;
    PaillierMatrix x_plus_r_receiver;
    solo::ahepaillier::BigNum two_power_64(solo::ahepaillier::BigNum::One());
    solo::ahepaillier::utils::bn_lshift(two_power_64, 64);
    solo::ahepaillier::BigNum two_power_64_plus_lambda(solo::ahepaillier::BigNum::One());
    solo::ahepaillier::utils::bn_lshift(two_power_64_plus_lambda, 64 + kStatisticalLambda);
    x_plus_r.set_party(1 - party_he_stored_owner);
    if (party_id_ == party_he_stored_owner) {
        for (std::size_t i = 0; i < in.size(); i++) {
            solo::ahepaillier::BigNum r =
                    two_power_64_plus_lambda +
                    (solo::ahepaillier::utils::get_random_bn(static_cast<int>(kPaillierKeySize)) %
                            (*(paillier_engine_->get_pk_other()->getN()) - two_power_64_plus_lambda));
            random_r_buffer.emplace_back(r);
            out(i) = ipcl_bn_to_int64(((-1) * r) % two_power_64);
        }
        row = in.rows();
        col = in.cols();
        out.resize(row, col);
        x_plus_r.resize(row, col);
        solo::ahepaillier::Plaintext pt;
        paillier_engine_->bn_to_pt(random_r_buffer, pt);
        paillier_engine_->add(in.ciphers(), pt, x_plus_r.ciphers());
        net->send_data(&row, sizeof(std::size_t));
        net->send_data(&col, sizeof(std::size_t));
        send_cipher(net, x_plus_r, kPaillierKeySize);
    } else {
        net->recv_data(&row, sizeof(std::size_t));
        net->recv_data(&col, sizeof(std::size_t));
        out.resize(row, col);
        x_plus_r.resize(row, col);
        recv_cipher(net, paillier_engine_->get_pk(), x_plus_r, kPaillierKeySize);
        paillier_engine_->decrypt(x_plus_r, out.shares());
    }
}

void Duet::a2h(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, PaillierMatrix& out) {
    PaillierMatrix enc_share;
    std::size_t party_he_receiver = 1 - out.party();

    if (party_id_ != party_he_receiver) {
        enc_share.resize(in.rows(), in.cols());
        paillier_engine_->encrypt(in.shares(), enc_share);
        enc_share.set_party(party_id_);
        send_cipher(net, enc_share, kPaillierKeySize);
    } else {
        recv_cipher(net, paillier_engine_->get_pk_other(), enc_share, kPaillierKeySize);
        PaillierMatrix self_share;
        paillier_engine_->encrypt(in.shares(), self_share, 0);
        self_share.set_party(1 - party_id_);
        paillier_engine_->add(self_share.ciphers(), enc_share.ciphers(), out.ciphers());
        out.resize(in.rows(), in.cols());
    }
}

}  // namespace duet
}  // namespace petace
