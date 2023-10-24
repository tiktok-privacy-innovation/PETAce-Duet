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

#include "test_utils.h"

#include <stdlib.h>

#include "duet/util/consts.h"

namespace petace {
namespace duet {

bool is_equal_plain_matrix(const Matrix<double>& m1, const Matrix<double>& m2, double e) {
    if ((m1.rows() != m2.rows()) || (m1.cols() != m2.cols()) || m1.size() == 0) {
        return false;
    }
    for (std::size_t i = 0; i < static_cast<std::size_t>(m1.size()); i++) {
        double d = m1(i) - m2(i);
        if ((d < -e) || (d > e)) {
            return false;
        }
    }

    return true;
}

bool is_equal_private_matrix(const PrivateMatrix<double>& m1, const PrivateMatrix<double>& m2, double e) {
    if (m1.party_id() != m2.party_id()) {
        return false;
    } else {
        return is_equal_plain_matrix(m1.matrix(), m2.matrix(), e);
    }
}

bool is_equal_plain_matrix(const Matrix<std::int64_t>& m1, const Matrix<std::int64_t>& m2) {
    if ((m1.rows() != m2.rows()) || (m1.cols() != m2.cols())) {
        return false;
    }
    for (std::size_t i = 0; i < static_cast<std::size_t>(m1.size()); i++) {
        int64_t d = m1(i) - m2(i);
        if (d != 0) {
            return false;
        }
    }

    return true;
}

bool is_equal_private_matrix(const PrivateMatrixBool& m1, const PrivateMatrixBool& m2) {
    if (m1.party_id() != m2.party_id()) {
        return false;
    } else {
        return is_equal_plain_matrix(m1.matrix(), m2.matrix());
    }
}

double get_rand_double() {
    double res;
    double threshold = 1.0 * (1 << kFixedPointPrecision);

    do {
        res = 1.0 * static_cast<double>(random()) / static_cast<double>(random());
        if (random() & 1) {
            res = -res;
        }
    } while ((res > threshold) || (res < -threshold));

    return res;
}

void get_rand_plain_matrix(Matrix<double>& plain) {
    for (std::size_t i = 0; i < static_cast<std::size_t>(plain.size()); i++) {
        plain(i) = get_rand_double();
    }
}

void get_rand_private_matrix(PrivateMatrix<double>& private_m, std::size_t party_id) {
    if (private_m.party_id() == party_id) {
        get_rand_plain_matrix(private_m.matrix());
    }
}

std::int64_t get_rand_int64() {
    std::int64_t res;
    std::int64_t threshold = (1 << kFixedPointPrecision);
    do {
        res = static_cast<std::int64_t>(random());
        if (random() & 1) {
            res = -res;
        }
    } while ((res > threshold) || (res < -threshold));
    return res;
}

std::int64_t get_rand_bool() {
    return static_cast<std::int64_t>(random()) % 2;
}

void get_rand_plain_bool_matrix(Matrix<std::int64_t>& plain) {
    for (std::size_t i = 0; i < static_cast<std::size_t>(plain.size()); i++) {
        plain(i) = get_rand_bool();
    }
}

void get_rand_private_bool_matrix(PrivateMatrixBool& private_m, std::size_t party_id) {
    if (private_m.party_id() == party_id) {
        get_rand_plain_bool_matrix(private_m.matrix());
    }
}

void get_rand_public_bool_matrix(PublicMatrixBool& private_m) {
    get_rand_plain_bool_matrix(private_m.matrix());
}

void get_rand_plain_matrix(Matrix<std::int64_t>& plain) {
    for (std::size_t i = 0; i < static_cast<std::size_t>(plain.size()); i++) {
        plain(i) = get_rand_int64();
    }
}

void get_rand_private_matrix(PrivateMatrix<int64_t>& private_m, std::size_t party_id) {
    if (private_m.party_id() == party_id) {
        get_rand_plain_matrix(private_m.matrix());
    }
}

}  // namespace duet
}  // namespace petace
