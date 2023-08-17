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
    block seed = _mm_set_epi64x(rand_generator_->get_unique_rand(), rand_generator_->get_unique_rand());
    // DoublePrg::fixed_key_aes.set_key(seed);
    ot_generator_ = std::make_shared<OTGenerator>(party_id_, seed);
    ot_generator_->initialize(net);
    st_generator_ = std::make_shared<STGenerator>(party_id_, ot_generator_);

    // set triplet
    rand_bool_triplet_generator_ = std::make_shared<BooleanTriplet>();
    rand_bool_triplet_generator_->initialize(net, party_id_, ot_generator_);

    rand_arith_triplet_generator_ = std::make_shared<ArithmeticTriplet>();
    rand_arith_triplet_generator_->initialize(net, party_id_, rand_generator_, ot_generator_);

    return;
}

void Duet::share(const std::shared_ptr<network::Network>& net, std::size_t party, const PlainMatrix<double>& in,
        ArithMatrix& out) {
    if ((in.size() == 0) && (party == party_id_)) {
        throw std::invalid_argument("matrix size is equal to zero.");
    }

    std::size_t row;
    std::size_t col;
    if (party == party_id_) {
        row = in.rows();
        col = in.cols();
        net->send_data(&row, sizeof(std::size_t));
        net->send_data(&col, sizeof(std::size_t));
    } else {
        net->recv_data(&row, sizeof(std::size_t));
        net->recv_data(&col, sizeof(std::size_t));
    }

    out.resize(row, col);
    if (party == party_id_) {
        for (std::size_t i = 0; i < out.size(); i++) {
            out.shares(i) = float_to_fixed(in(i)) - rand_generator_->get_common_rand();
        }
    } else {
        for (std::size_t i = 0; i < out.size(); i++) {
            out.shares(i) = rand_generator_->get_common_rand();
        }
    }
    return;
}

void Duet::share(const std::shared_ptr<network::Network>& net, std::size_t party, const PlainMatrix<std::int64_t>& in,
        ArithMatrix& out) {
    if ((in.size() == 0) && (party == party_id_)) {
        throw std::invalid_argument("matrix size is equal to zero.");
    }

    std::size_t row;
    std::size_t col;
    if (party == party_id_) {
        row = in.rows();
        col = in.cols();
        net->send_data(&row, sizeof(std::size_t));
        net->send_data(&col, sizeof(std::size_t));
    } else {
        net->recv_data(&row, sizeof(std::size_t));
        net->recv_data(&col, sizeof(std::size_t));
    }

    out.resize(row, col);
    if (party == party_id_) {
        for (std::size_t i = 0; i < out.size(); i++) {
            out.shares(i) = in(i)-rand_generator_->get_common_rand();
        }
    } else {
        for (std::size_t i = 0; i < out.size(); i++) {
            out.shares(i) = rand_generator_->get_common_rand();
        }
    }
    return;
}

void Duet::reveal(const std::shared_ptr<network::Network>& net, std::size_t party, const ArithMatrix& in,
        PlainMatrix<double>& out) {
    if (in.size() == 0) {
        throw std::invalid_argument("matrix size is equal to zero.");
    }
    ArithMatrix fixed_matrix(in.rows(), in.cols());
    if (party_id_ != party) {
        send_matrix(net, &in.shares, 1);
    } else {
        recv_matrix(net, &fixed_matrix.shares, 1);
    }

    out.resize(in.rows(), in.cols());

    if (party_id_ == party) {
        fixed_matrix = fixed_matrix + in;

        for (std::size_t i = 0; i < in.size(); i++) {
            out(i) = fixed_to_float(fixed_matrix.shares(i));
        }
    }
    return;
}

void Duet::reveal_bool(const std::shared_ptr<network::Network>& net, std::size_t party, const BoolMatrix& in,
        PlainMatrix<std::int64_t>& out) {
    if (in.size() == 0) {
        throw std::invalid_argument("matrix size is equal to zero.");
    }
    ArithMatrix fixed_matrix(in.rows(), in.cols());
    if (party_id_ != party) {
        send_matrix(net, &in.shares, 1);
    } else {
        recv_matrix(net, &fixed_matrix.shares, 1);
    }

    out.resize(in.rows(), in.cols());

    if (party_id_ == party) {
        for (std::size_t i = 0; i < in.size(); i++) {
            out(i) = fixed_matrix.shares(i) ^ in.shares(i);
        }
    }
    return;
}

void Duet::add(const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.shares = x.shares + y.shares;
    return;
}

void Duet::sub(const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z) const {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    z.shares = x.shares - y.shares;
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
        triplet_a.shares(i) = triplet[0];
        triplet_b.shares(i) = triplet[1];
        triplet_c.shares(i) = triplet[2];
    }

    BoolMatrix e(row, col);
    BoolMatrix f(row, col);
    for (std::size_t i = 0; i < e.size(); ++i) {
        e.shares(i) = x.shares(i) ^ triplet_a.shares(i);
        f.shares(i) = y.shares(i) ^ triplet_b.shares(i);
    }

    BoolMatrix reveal_e(row, col);
    BoolMatrix reveal_f(row, col);
    if (party_id_ == 0) {
        send_matrix(net, &e.shares, 1);
        send_matrix(net, &f.shares, 1);
        recv_matrix(net, &reveal_e.shares, 1);
        recv_matrix(net, &reveal_f.shares, 1);
    } else {
        recv_matrix(net, &reveal_e.shares, 1);
        recv_matrix(net, &reveal_f.shares, 1);
        send_matrix(net, &e.shares, 1);
        send_matrix(net, &f.shares, 1);
    }

    for (std::size_t i = 0; i < reveal_e.size(); ++i) {
        reveal_e.shares(i) = reveal_e.shares(i) ^ e.shares(i);
        reveal_f.shares(i) = reveal_f.shares(i) ^ f.shares(i);
    }

    if (party_id_ == 0) {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z.shares(i) = (reveal_f.shares(i) & triplet_a.shares(i)) ^ (reveal_e.shares(i) & triplet_b.shares(i)) ^
                          triplet_c.shares(i);
        }
    } else {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z.shares(i) = (reveal_e.shares(i) & reveal_f.shares(i)) ^ (reveal_f.shares(i) & triplet_a.shares(i)) ^
                          (reveal_e.shares(i) & triplet_b.shares(i)) ^ triplet_c.shares(i);
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
        triplet_a.shares(i) = triplet[0];
        triplet_b.shares(i) = triplet[1];
        triplet_c.shares(i) = triplet[2];
    }

    ArithMatrix triplet_a_1(row, col);
    ArithMatrix triplet_b_1(row, col);
    ArithMatrix triplet_c_1(row, col);

    if (party_id_ == 0) {
        send_matrix(net, &triplet_a.shares, 1);
        send_matrix(net, &triplet_b.shares, 1);
        send_matrix(net, &triplet_c.shares, 1);
    } else {
        recv_matrix(net, &triplet_a_1.shares, 1);
        recv_matrix(net, &triplet_b_1.shares, 1);
        recv_matrix(net, &triplet_c_1.shares, 1);
    }

    ArithMatrix e(row, col);
    ArithMatrix f(row, col);
    for (std::size_t i = 0; i < e.size(); ++i) {
        e.shares(i) = x.shares(i) - triplet_a.shares(i);
        f.shares(i) = y.shares(i) - triplet_b.shares(i);
    }

    ArithMatrix reveal_e(row, col);
    ArithMatrix reveal_f(row, col);
    if (party_id_ == 0) {
        send_matrix(net, &e.shares, 1);
        send_matrix(net, &f.shares, 1);
        recv_matrix(net, &reveal_e.shares, 1);
        recv_matrix(net, &reveal_f.shares, 1);
    } else {
        recv_matrix(net, &reveal_e.shares, 1);
        recv_matrix(net, &reveal_f.shares, 1);
        send_matrix(net, &e.shares, 1);
        send_matrix(net, &f.shares, 1);
    }

    for (std::size_t i = 0; i < reveal_e.size(); ++i) {
        reveal_e.shares(i) = reveal_e.shares(i) + e.shares(i);
        reveal_f.shares(i) = reveal_f.shares(i) + f.shares(i);
    }

    if (party_id_ == 0) {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z.shares(i) = ((reveal_f.shares(i) * triplet_a.shares(i)) + (reveal_e.shares(i) * triplet_b.shares(i)) +
                                  triplet_c.shares(i)) >>
                          16;
        }
    } else {
        for (std::size_t i = 0; i < z.size(); ++i) {
            z.shares(i) = ((reveal_e.shares(i) * reveal_f.shares(i)) + (reveal_f.shares(i) * triplet_a.shares(i)) +
                                  (reveal_e.shares(i) * triplet_b.shares(i)) + triplet_c.shares(i)) >>
                          16;
        }
    }
    return;
}

void Duet::kogge_stone_ppa(
        const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const BoolMatrix& y, BoolMatrix& z) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t row = x.rows();
    std::size_t col = x.cols();
    std::size_t depth = 6;
    BoolMatrix g1(row, col);
    BoolMatrix p1(row, col);
    BoolMatrix g(row, col);
    BoolMatrix p(row, col);
    std::int64_t keep_masks[6] = {0x0000000000000001, 0x0000000000000003, 0x000000000000000f, 0x00000000000000ff,
            0x000000000000ffff, 0x00000000ffffffff};

    elementwise_bool_mul(net, x, y, g);
    for (std::size_t i = 0; i < x.size(); i++) {
        p.shares(i) = x.shares(i) ^ y.shares(i);
    }

    for (std::size_t i = 0; i < depth; i++) {
        std::int64_t shift = 1L << i;
        for (std::size_t k = 0; k < p.size(); k++) {
            p1.shares(k) = p.shares(k) << shift;
        }
        for (std::size_t k = 0; k < g.size(); k++) {
            g1.shares(k) = g.shares(k) << shift;
        }

        if (party_id_ == 0) {
            for (std::size_t k = 0; k < p.size(); k++) {
                p1.shares(k) ^= keep_masks[i];
            }
        }
        elementwise_bool_mul(net, p, g1, g1);
        for (std::size_t k = 0; k < g.size(); k++) {
            g.shares(k) ^= g1.shares(k);
        }
        elementwise_bool_mul(net, p, p1, p);
    }

    for (std::size_t k = 0; k < g.size(); k++) {
        g1.shares(k) = g.shares(k) << 1;
    }

    for (std::size_t k = 0; k < g.size(); k++) {
        z.shares(k) = g1.shares(k) ^ x.shares(k) ^ y.shares(k);
    }

    return;
}

void Duet::a2b(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, BoolMatrix& z) {
    if (x.size() == 0) {
        throw std::invalid_argument("matrix size is equal to zero.");
    }
    std::size_t size = x.size();
    std::size_t row = x.rows();
    std::size_t col = x.cols();

    BoolMatrix input_0(row, col);
    BoolMatrix input_1(row, col);
    if (party_id_ == 0) {
        for (std::size_t j = 0; j < size; j++) {
            input_0.shares(j) = x.shares(j) ^ rand_generator_->get_common_rand();
            input_1.shares(j) = rand_generator_->get_common_rand();
        }
    } else {
        for (std::size_t j = 0; j < size; j++) {
            input_0.shares(j) = rand_generator_->get_common_rand();
            input_1.shares(j) = x.shares(j) ^ rand_generator_->get_common_rand();
        }
    }
    kogge_stone_ppa(net, input_0, input_1, z);
    for (std::size_t j = 0; j < size; j++) {
        z.shares(j) = (z.shares(j) >> 63) & 0x1;
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

void Duet::greater(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const PlainMatrix<double>& y,
        BoolMatrix& z) {
    if (x.size() != static_cast<std::size_t>(y.size())) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t size = x.size();
    ArithMatrix c(x.rows(), x.cols());
    if (party_id_ == 0) {
        for (std::size_t j = 0; j < size; j++) {
            c.shares(j) = float_to_fixed(y(j)) - x.shares(j);
        }
    } else {
        c.shares = -x.shares;
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

void Duet::less(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const PlainMatrix<double>& y,
        BoolMatrix& z) {
    if (x.size() != static_cast<std::size_t>(y.size())) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::size_t size = x.size();
    ArithMatrix c(x.rows(), x.cols());
    if (party_id_ == 0) {
        for (std::size_t j = 0; j < size; j++) {
            c.shares(j) = x.shares(j) - float_to_fixed(y(j));
        }
    } else {
        c.shares = x.shares;
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
    c.shares = x.shares - y.shares;

    a2b(net, c, z);

    for (std::size_t j = 0; j < x.size(); j++) {
        if (party_id_ == 0) {
            z.shares(j) ^= 1;
        }
    }

    return;
}

void Duet::equal(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    BoolMatrix c0(x.rows(), x.cols());
    BoolMatrix c1(x.rows(), x.cols());

    greater_equal(net, x, y, c0);
    greater_equal(net, y, x, c1);

    for (std::size_t j = 0; j < x.size(); j++) {
        z.shares(j) = (c0.shares(j)) ^ (c1.shares(j));
        if (party_id_ == 0) {
            z.shares(j) ^= 1;
        }
    }

    return;
}

void Duet::sum(const ArithMatrix& x, ArithMatrix& z) const {
    if (x.size() == 0) {
        throw std::invalid_argument("matrix size is equal to zero.");
    }
    z.resize(1, x.cols());
    z.shares.row(0) = x.shares.colwise().sum();
    return;
}

void Duet::millionaire(const std::shared_ptr<network::Network>& net, PlainMatrix<std::int64_t>& x, BoolMatrix& y,
        std::size_t bit_length) {
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
            leaf_lt.push_back(rand_generator_->get_unique_rand_int8() & 0x1);
            leaf_eq.push_back(rand_generator_->get_unique_rand_int8() & 0x1);
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
                temp_and_input_1.shares(k * cur_num_blocks + j) = leaf_lt[k * cur_num_blocks * 2 + 2 * j];
                temp_and_input_2.shares(k * cur_num_blocks + j) = leaf_eq[k * cur_num_blocks * 2 + 2 * j + 1];
            }
        }

        pack_and_evaluate_and(net, temp_and_input_1, temp_and_input_2, temp_and_res);

        for (std::size_t j = 0; j < cur_num_blocks; ++j) {
            for (std::size_t k = 0; k < matrix_size; k++) {
                leaf_lt_temp[k * cur_num_blocks + j] =
                        (leaf_lt[k * cur_num_blocks * 2 + 2 * j + 1] ^ temp_and_res.shares(k * cur_num_blocks + j)) &
                        0x1;
                temp_and_input_1.shares(k * cur_num_blocks + j) = leaf_eq[k * cur_num_blocks * 2 + 2 * j];
                temp_and_input_2.shares(k * cur_num_blocks + j) = leaf_eq[k * cur_num_blocks * 2 + 2 * j + 1];
            }
        }
        leaf_lt.assign(leaf_lt_temp.begin(), leaf_lt_temp.end());

        pack_and_evaluate_and(net, temp_and_input_1, temp_and_input_2, temp_and_res);
        for (std::size_t j = 0; j < cur_num_blocks; ++j) {
            for (std::size_t k = 0; k < matrix_size; k++) {
                leaf_eq[k * cur_num_blocks + j] = temp_and_res.shares(k * cur_num_blocks + j) & 0x1;
            }
        }
    }
    for (std::size_t i = 0; i < matrix_size; ++i) {
        y.shares(i) = leaf_lt[i] & 0x1;
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
    x_packed.shares.setZero();
    y_packed.shares.setZero();
    for (size_t i = 0; i < num_pack; i++) {
        for (size_t j = 0; j < 64; j++) {
            if (i * 64 + j < num_elements) {
                x_packed.shares(i) |= (x.shares(i * 64 + j) << j);
                y_packed.shares(i) |= (y.shares(i * 64 + j) << j);
            }
        }
    }
    elementwise_bool_mul(net, x_packed, y_packed, z_packed);

    // step 2: unpacking
    for (size_t i = 0; i < num_pack; i++) {
        for (size_t j = 0; j < 64; j++) {
            if (i * 64 + j < num_elements) {
                z.shares(i * 64 + j) = (z_packed.shares(i) >> j) & 0x1;
            }
        }
    }
}

void Duet::less_than_zero(
        const std::shared_ptr<network::Network>& net, const ArithMatrix& x, BoolMatrix& y, std::size_t bit_length) {
    std::size_t num_elements = x.size();

    // step 1: split the inputs to be msb and the rest
    PlainMatrix<std::int64_t> msb(1, num_elements);
    PlainMatrix<std::int64_t> rest(1, num_elements);
    for (std::size_t i = 0; i < num_elements; i++) {
        msb(i) = static_cast<std::uint64_t>(x.shares(i)) >> (bit_length - 1);
        rest(i) = x.shares(i) & ((1 << (bit_length - 1)) - 1);
    }
    // step 2: compute the carry bit of rest
    BoolMatrix carry(1, num_elements);
    if (party_id_ == 0) {
        // party 0 set the input of millionare to be 2^{l - 1} - 1 - rest
        PlainMatrix<std::int64_t> millionare_input(1, num_elements);
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
        y.shares(i) = msb(i) ^ carry.shares(i);
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
        r.shares(i) = rand_generator_->get_unique_rand();
    }

    for (std::size_t i = 0; i < size; i++) {
        if (x.shares(i) == 0) {
            s0.shares(i) = -r.shares(i);
            s1.shares(i) = y.shares(i) - r.shares(i);
        } else {
            s0.shares(i) = y.shares(i) - r.shares(i);
            s1.shares(i) = -r.shares(i);
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
            k[i] = choice ^ x.shares(i);
            rb.shares(i) = msg;
        }

        send_bool(net, k, size);
        recv_matrix(net, &y0.shares, 1);
        recv_matrix(net, &y1.shares, 1);
        for (std::size_t i = 0; i < size; i++) {
            if (x.shares(i) == 0) {
                z.shares(i) = y0.shares(i) ^ rb.shares(i);
            } else {
                z.shares(i) = y1.shares(i) ^ rb.shares(i);
            }
        }
    } else {
        for (std::size_t i = 0; i < size; i++) {
            std::vector<std::int64_t> msg;
            ot_generator_->get_random_ot(net, msg);
            y0.shares(i) = msg[0];
            y1.shares(i) = msg[1];
        }
        recv_bool(net, k, size);

        for (std::size_t i = 0; i < size; i++) {
            if (k[i] == 0) {
                y0.shares(i) ^= s0.shares(i);
                y1.shares(i) ^= s1.shares(i);

            } else {
                std::int64_t t = s0.shares(i) ^ y1.shares(i);
                y1.shares(i) = s1.shares(i) ^ y0.shares(i);
                y0.shares(i) = t;
            }
        }
        send_matrix(net, &y0.shares, 1);
        send_matrix(net, &y1.shares, 1);
    }

    if (party_id_ == 1) {
        for (std::size_t i = 0; i < size; i++) {
            std::int8_t choice;
            std::int64_t msg;
            ot_generator_->get_random_ot(net, choice, msg);
            k[i] = choice ^ x.shares(i);
            rb.shares(i) = msg;
        }
        send_bool(net, k, size);
        recv_matrix(net, &y0.shares, 1);
        recv_matrix(net, &y1.shares, 1);
        for (std::size_t i = 0; i < size; i++) {
            if (x.shares(i) == 0) {
                z.shares(i) = y0.shares(i) ^ rb.shares(i);
            } else {
                z.shares(i) = y1.shares(i) ^ rb.shares(i);
            }
        }
    } else {
        for (std::size_t i = 0; i < size; i++) {
            std::vector<std::int64_t> msg;
            ot_generator_->get_random_ot(net, msg);
            y0.shares(i) = msg[0];
            y1.shares(i) = msg[1];
        }
        recv_bool(net, k, size);
        for (std::size_t i = 0; i < size; i++) {
            if (k[i] == 0) {
                y0.shares(i) ^= s0.shares(i);
                y1.shares(i) ^= s1.shares(i);

            } else {
                std::int64_t t = s0.shares(i) ^ y1.shares(i);
                y1.shares(i) = s1.shares(i) ^ y0.shares(i);
                y0.shares(i) = t;
            }
        }
        send_matrix(net, &y0.shares, 1);
        send_matrix(net, &y1.shares, 1);
    }

    for (std::size_t i = 0; i < size; i++) {
        z.shares(i) += r.shares(i);
    }
    delete[] k;

    return;
}

void Duet::shuffle(const std::shared_ptr<network::Network>& net, std::size_t party, const PlainMatrix<double>& in,
        ArithMatrix& out) {
    std::size_t row;
    std::size_t col;
    std::size_t matrix_size;
    SecretSharedShuffle ss_shuffle;
    if (party == party_id_) {
        row = in.rows();
        col = in.cols();
        matrix_size = in.size();
        net->send_data(&row, sizeof(std::size_t));
        net->send_data(&col, sizeof(std::size_t));
        PlainMatrix<std::int64_t> x(row, col);
        PlainMatrix<std::int64_t> a;
        PlainMatrix<std::int64_t> b;
        PlainMatrix<std::int64_t> x_sub_a;
        for (std::size_t i = 0; i < matrix_size; ++i) {
            x(i) = float_to_fixed(in(i));
        }
        st_generator_->get_st(row, col, net, a, b);
        ss_shuffle.passive_phase_1(x, a, x_sub_a);
        out.shares = b;
        send_matrix(net, &x_sub_a, 1);

    } else {
        net->recv_data(&row, sizeof(std::size_t));
        net->recv_data(&col, sizeof(std::size_t));
        PlainMatrix<std::int64_t> delta;
        PlainMatrix<std::int64_t> x_sub_a(row, col);
        PlainMatrix<std::int64_t> share;
        Permutation p(row);
        st_generator_->get_st(row, col, net, p, delta);
        recv_matrix(net, &x_sub_a, 1);
        ss_shuffle.active_phase_1(p, x_sub_a, delta, share);
        out.shares = share;
    }
}

void Duet::shuffle(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, ArithMatrix& out) {
    std::size_t row;
    std::size_t col;
    SecretSharedShuffle ss_shuffle;
    row = in.rows();
    col = in.cols();
    if (party_id_ == 0) {
        PlainMatrix<std::int64_t> a0;
        PlainMatrix<std::int64_t> b0;
        PlainMatrix<std::int64_t> x_sub_a0;
        PlainMatrix<std::int64_t> share0;
        st_generator_->get_st(row, col, net, a0, b0);
        ss_shuffle.passive_phase_1(in.shares, a0, x_sub_a0);
        share0 = b0;
        send_matrix(net, &x_sub_a0, 1);

        PlainMatrix<std::int64_t> delta1;
        PlainMatrix<std::int64_t> x_sub_a1(row, col);
        PlainMatrix<std::int64_t> share1;
        Permutation p1(row);
        st_generator_->get_st(row, col, net, p1, delta1);
        recv_matrix(net, &x_sub_a1, 1);
        ss_shuffle.active_phase_1(p1, x_sub_a1, delta1, share1);
        out.shares = share1 + p1.permute(share0);
    } else {
        PlainMatrix<std::int64_t> delta0;
        PlainMatrix<std::int64_t> x_sub_a0(row, col);
        PlainMatrix<std::int64_t> share0;
        Permutation p0(row);
        st_generator_->get_st(row, col, net, p0, delta0);
        recv_matrix(net, &x_sub_a0, 1);
        ss_shuffle.active_phase_1(p0, x_sub_a0, delta0, share0);
        share0 = share0 + p0.permute(in.shares);

        PlainMatrix<std::int64_t> a1;
        PlainMatrix<std::int64_t> b1;
        PlainMatrix<std::int64_t> x_sub_a1;
        PlainMatrix<std::int64_t> share1;
        st_generator_->get_st(row, col, net, a1, b1);
        ss_shuffle.passive_phase_1(share0, a1, x_sub_a1);
        out.shares = b1;
        send_matrix(net, &x_sub_a1, 1);
    }
}

}  // namespace duet
}  // namespace petace
