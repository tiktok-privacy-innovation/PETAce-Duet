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

void Duet::sum(const ArithMatrix& x, ArithMatrix& z) const {
    if (x.size() == 0) {
        throw std::invalid_argument("sum(): matrix size is equal to zero.");
    }
    z.resize(1, x.cols());
    z.shares().row(0) = x.shares().colwise().sum();
    return;
}

void Duet::shuffle(
        const std::shared_ptr<network::Network>& net, const PrivateMatrix<double>& in, ArithMatrix& out) const {
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

void Duet::shuffle(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, ArithMatrix& out) const {
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

void Duet::shuffle(const std::shared_ptr<network::Network>& net, const PrivatePermutation& perm, const ArithMatrix& in,
        ArithMatrix& out) const {
    std::size_t row, col;
    row = in.rows();
    col = in.cols();

    PaillierMatrix cipher_a2h(row, col, 1 - perm.party_id());
    a2h(net, in, cipher_a2h);

    PaillierMatrix shuffle_cipher_a2h(row, col, 1 - perm.party_id());
    if (party_id_ == perm.party_id()) {
        std::vector<mpz_class> decode_cipher(row * col);
        ByteVector temp(kPaillierCipherSize);
        for (std::size_t i = 0; i < row * col; i++) {
            cipher_a2h(i).serialize_to_bytes(temp.data(), temp.size());
            mpz_bn_from_bytes(temp.data(), temp.size(), decode_cipher[i]);
        }

        std::vector<std::vector<mpz_class>> decode_cipher_dim_2(row, std::vector<mpz_class>(col));
        std::vector<std::vector<mpz_class>> shuffle_decode_cipher_dim_2(row, std::vector<mpz_class>(col));
        for (std::size_t i = 0; i < row; i++) {
            for (std::size_t j = 0; j < col; j++) {
                decode_cipher_dim_2[i][j] = decode_cipher[i * col + j];
            }
        }
        for (std::size_t i = 0; i < row; i++) {
            shuffle_decode_cipher_dim_2[i] = decode_cipher_dim_2[perm[i]];
        }
        std::vector<mpz_class> shuffle_decode_cipher(row * col);
        for (std::size_t i = 0; i < row; i++) {
            for (std::size_t j = 0; j < col; j++) {
                shuffle_decode_cipher[i * col + j] = shuffle_decode_cipher_dim_2[i][j];
            }
        }
        std::vector<solo::ahepaillier::Ciphertext> shuffle_ciphers(row * col);
        for (std::size_t i = 0; i < shuffle_decode_cipher.size(); i++) {
            mpz_bn_to_bytes(shuffle_decode_cipher[i], temp.data(), temp.size());
            shuffle_ciphers[i].deserialize_from_bytes(paillier_engine_->get_pk_other(), temp.data(), temp.size());
        }
        shuffle_cipher_a2h.set_ciphers(shuffle_ciphers);
    }

    out.resize(row, col);
    h2a(net, shuffle_cipher_a2h, out);
}

void Duet::quick_sort(const std::shared_ptr<network::Network>& net, const std::vector<std::size_t>& index,
        ArithMatrix& matrix) const {
    // shuffled before sort
    ArithMatrix cmp_matrix(matrix.rows(), matrix.cols());
    std::vector<std::size_t> next_index;
    BoolMatrix cmp_share;
    PublicMatrixBool cmp_ret;

    bool needs_do = false;
    for (std::size_t i = 0; i < index.size(); i += 2) {
        if (index[i] < index[i + 1]) {
            needs_do = true;
        }
    }
    if (!needs_do) {
        return;
    }

    // let cmp parallel
    for (std::size_t i = 0; i < index.size(); i += 2) {
        if (index[i] < index[i + 1]) {
            for (std::size_t j = index[i] + 1; j <= index[i + 1]; ++j) {
                cmp_matrix(j) = matrix(index[i]);
            }
        }
    }

    greater(net, cmp_matrix, matrix, cmp_share);

    reveal(net, cmp_share, cmp_ret);

    for (std::size_t i = 0; i < index.size(); i += 2) {
        if (index[i + 1] > index[i]) {
            std::size_t end = index[i + 1];
            std::size_t begin = index[i] + 1;
            while (begin <= end) {
                if (cmp_ret(begin) == 0) {
                    while (cmp_ret(end) == 0 && end >= begin) {
                        end--;
                    }
                    if (end >= begin) {
                        std::swap(matrix(begin), matrix(end));
                        end--;
                        begin++;
                    }
                } else {
                    begin++;
                }
            }
            std::swap(matrix(end), matrix(index[i]));
            if (end != 0) {
                next_index.push_back(index[i]);
                next_index.push_back(end - 1);
            }

            next_index.push_back(end + 1);
            next_index.push_back(index[i + 1]);
        }
    }
    quick_sort(net, next_index, matrix);
}

void Duet::quick_sort(const std::shared_ptr<network::Network>& net, const std::vector<std::size_t>& index,
        std::size_t col_index, ArithMatrix& matrix) const {
    // must shuffle
    ArithMatrix cmp_matrix(matrix.rows(), 1);
    ArithMatrix sort_column(matrix.rows(), 1);
    std::vector<std::size_t> next_index;
    BoolMatrix cmp_share;
    PublicMatrixBool cmp_ret;
    ArithMatrix tmp(1, matrix.cols());

    bool needs_do = false;
    for (std::size_t i = 0; i < index.size(); i += 2) {
        if (index[i] < index[i + 1]) {
            needs_do = true;
        }
    }
    if (!needs_do) {
        return;
    }

    for (std::size_t i = 0; i < matrix.rows(); i++) {
        sort_column(i) = matrix(i * matrix.cols() + col_index);
    }

    // let cmp parallel
    for (std::size_t i = 0; i < index.size(); i += 2) {
        if (index[i] < index[i + 1]) {
            for (std::size_t j = index[i] + 1; j <= index[i + 1]; ++j) {
                cmp_matrix(j) = sort_column(index[i]);
            }
        }
    }

    greater(net, cmp_matrix, sort_column, cmp_share);

    reveal(net, cmp_share, cmp_ret);

    for (std::size_t i = 0; i < index.size(); i += 2) {
        if (index[i + 1] > index[i]) {
            std::size_t end = index[i + 1];
            std::size_t begin = index[i] + 1;
            while (begin <= end) {
                if (cmp_ret(begin) == 0) {
                    while (cmp_ret(end) == 0 && end >= begin) {
                        end--;
                    }
                    if (end >= begin) {
                        tmp.shares().row(0) = matrix.shares().row(begin);
                        matrix.shares().row(begin) = matrix.shares().row(end);
                        matrix.shares().row(end) = tmp.shares().row(0);

                        end--;
                        begin++;
                    }
                } else {
                    begin++;
                }
            }
            tmp.shares().row(0) = matrix.shares().row(end);
            matrix.shares().row(end) = matrix.shares().row(index[i]);
            matrix.shares().row(index[i]) = tmp.shares().row(0);

            if (end != 0) {
                next_index.push_back(index[i]);
                next_index.push_back(end - 1);
            }

            next_index.push_back(end + 1);
            next_index.push_back(index[i + 1]);
        }
    }
    quick_sort(net, next_index, col_index, matrix);
}

void Duet::quick_sort(const std::shared_ptr<network::Network>& net, ArithMatrix& matrix) const {
    PrivatePermutation perm0(matrix.rows());
    perm0.set_party_id(0);
    shuffle(net, perm0, matrix, matrix);

    PrivatePermutation perm1(matrix.rows());
    perm1.set_party_id(1);
    shuffle(net, perm1, matrix, matrix);

    std::vector<std::size_t> index = {0, matrix.size() - 1};
    quick_sort(net, index, matrix);
}

void Duet::quick_sort(const std::shared_ptr<network::Network>& net, std::size_t col_index, const ArithMatrix& in,
        ArithMatrix& out) const {
    out.shares() = in.shares();

    PrivatePermutation perm0(out.rows());
    perm0.set_party_id(0);
    shuffle(net, perm0, out, out);

    PrivatePermutation perm1(out.rows());
    perm1.set_party_id(1);
    shuffle(net, perm1, out, out);

    std::vector<std::size_t> index = {0, in.rows() - 1};
    quick_sort(net, index, col_index, out);
}

void Duet::argmax_and_max(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, ArithMatrix& max_index,
        ArithMatrix& max_value) const {
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

void Duet::argmin_and_min(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, ArithMatrix& min_index,
        ArithMatrix& min_value) const {
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
        less(net, up, down, cmp);
        multiplexer(net, cmp, down, up, now);
        multiplexer(net, cmp, index_down, index_up, index_now);

        if (rem_num == 1) {
            vstack(now, rem, tmp);
            now = tmp;
            vstack(index_now, index_rem, index_tmp);
            index_now = index_tmp;
        }
    }
    min_value = now;
    min_index = index_now;
}

}  // namespace duet
}  // namespace petace
