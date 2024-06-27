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

void Duet::encrypt(const PrivateMatrix<std::int64_t>& input, PaillierMatrix& output, size_t using_self_pk) const {
    paillier_engine_->encrypt(input, output, using_self_pk);
    output.set_party(party_id_ * using_self_pk + (1 - party_id_) * (1 - using_self_pk));
}

void Duet::decrypt(const PaillierMatrix& input, PrivateMatrix<std::int64_t>& output) const {
    paillier_engine_->decrypt(input, output);
}

void Duet::encrypt(const PrivateMatrix<double>& input, PaillierMatrix& output, size_t using_self_pk) const {
    paillier_engine_->encrypt_double(input, output, using_self_pk);
    output.set_party(party_id_ * using_self_pk + (1 - party_id_) * (1 - using_self_pk));
}

void Duet::decrypt(const PaillierMatrix& input, PrivateMatrix<double>& output) const {
    paillier_engine_->decrypt_double(input, output);
}

void Duet::add(
        const PrivateMatrix<std::int64_t>& x, const PaillierMatrix& y, PaillierMatrix& z, const bool self_pk) const {
    if (static_cast<std::size_t>(x.size()) != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::vector<solo::ahepaillier::Plaintext> pt(x.size());
    for (std::size_t i = 0; i < static_cast<std::size_t>(x.size()); ++i) {
        paillier_engine_->encode(x(i), pt[i]);
    }
    paillier_engine_->add(y.ciphers(), pt, z.ciphers(), self_pk);

    z.resize(y.rows(), y.cols());
    z.set_party(y.party());
}

void Duet::add(const PaillierMatrix& x, const PaillierMatrix& y, PaillierMatrix& z, const bool self_pk) const {
    if (static_cast<std::size_t>(x.size()) != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    paillier_engine_->add(x.ciphers(), y.ciphers(), z.ciphers(), self_pk);
    z.resize(y.rows(), y.cols());
    z.set_party(y.party());
}

void Duet::add(const PrivateMatrix<double>& x, const PaillierMatrix& y, PaillierMatrix& z, const bool self_pk) const {
    if (static_cast<std::size_t>(x.size()) != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::vector<solo::ahepaillier::Plaintext> pt;
    pt.resize(y.size());

    for (std::size_t i = 0; i < static_cast<std::size_t>(x.size()); ++i) {
        paillier_engine_->encode(double_to_fixed(x(i)), pt[i]);
    }
    paillier_engine_->add(y.ciphers(), pt, z.ciphers(), self_pk);
    z.resize(y.rows(), y.cols());
    z.set_party(y.party());
}

void Duet::mul(
        const PrivateMatrix<std::int64_t>& x, const PaillierMatrix& y, PaillierMatrix& z, const bool self_pk) const {
    if (static_cast<std::size_t>(x.size()) != y.size()) {
        throw std::invalid_argument("matrix size is not equal.");
    }
    std::vector<solo::ahepaillier::Plaintext> pt(x.size());
    for (std::size_t i = 0; i < static_cast<std::size_t>(x.size()); ++i) {
        paillier_engine_->encode(x(i), pt[i]);
    }

    z.resize(y.rows(), y.cols());
    paillier_engine_->mul(y.ciphers(), pt, z.ciphers(), self_pk);
    z.set_party(y.party());
}

void Duet::h2a(const std::shared_ptr<network::Network>& net, const PaillierMatrix& in, ArithMatrix& out) const {
    std::size_t row;
    std::size_t col;
    std::size_t party_he_stored_owner = 1 - in.party();
    std::vector<mpz_class> random_r_buffer;
    PaillierMatrix x_plus_r;
    PaillierMatrix x_plus_r_receiver;
    mpz_class two_power_64 = mpz_class(kTwoPowerSixtyFour);
    mpz_class two_power_64_plus_lambda = mpz_class(kTwoPowerSixtyFour);
    two_power_64_plus_lambda = two_power_64 * mpz_class(pow(2, kStatisticalLambda));
    x_plus_r.set_party(1 - party_he_stored_owner);
    if (party_id_ == party_he_stored_owner) {
        for (std::size_t i = 0; i < in.size(); i++) {
            mpz_class r = get_random_mpz(rand_generator_, kPaillierKeySize / 2);
            while (r < two_power_64_plus_lambda) {
                r = get_random_mpz(rand_generator_, kPaillierKeySize / 2);
            }
            random_r_buffer.emplace_back(r);
            mpz_class out_r;
            mpz_class r_minus = r * mpz_class(-1);
            mpz_mod(out_r.get_mpz_t(), r_minus.get_mpz_t(), mpz_class(kTwoPowerSixtyFour).get_mpz_t());
            out(i) = static_cast<std::int64_t>(out_r.get_ui());
        }
        row = in.rows();
        col = in.cols();
        out.resize(row, col);
        x_plus_r.resize(row, col);
        std::vector<solo::ahepaillier::Plaintext> pt;
        pt.resize(row * col);
        std::vector<petace::solo::Byte> temp((kPaillierKeySize / 2 + 7) / 8);
        for (std::size_t i = 0; i < row * col; i++) {
            mpz_bn_to_bytes(random_r_buffer[i], temp.data(), temp.size());
            pt[i].deserialize_from_bytes(temp.data(), temp.size());
        }
        paillier_engine_->add(in.ciphers(), pt, x_plus_r.ciphers(), false);
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

void Duet::a2h(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, PaillierMatrix& out) const {
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
        out.resize(in.rows(), in.cols());
        add(self_share, enc_share, out, false);
    }
}

}  // namespace duet
}  // namespace petace
