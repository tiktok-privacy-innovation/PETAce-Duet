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

#include "duet/util/defines.h"
#include "duet/util/matrix.h"

namespace petace {
namespace duet {

bool is_equal_plain_matrix(const Matrix<double>& m1, const Matrix<double>& m2, double e);

bool is_equal_private_matrix(const PrivateMatrix<double>& m1, const PrivateMatrix<double>& m2, double e);

bool is_equal_plain_matrix(const Matrix<std::int64_t>& m1, const Matrix<std::int64_t>& m2);

bool is_equal_private_matrix(const PrivateMatrixBool& m1, const PrivateMatrixBool& m2);

double get_rand_double();

void get_rand_plain_matrix(Matrix<double>& plain);

void get_rand_private_matrix(PrivateMatrix<double>& private_m, std::size_t party_id);

std::int64_t get_rand_int64();

std::int64_t get_rand_bool();

void get_rand_plain_bool_matrix(Matrix<std::int64_t>& plain);

void get_rand_private_bool_matrix(PrivateMatrixBool& private_m, std::size_t party_id);

void get_rand_public_bool_matrix(PublicMatrixBool& private_m);

void get_rand_plain_matrix(Matrix<std::int64_t>& plain);

void get_rand_private_matrix(PrivateMatrix<int64_t>& private_m, std::size_t party_id);

}  // namespace duet
}  // namespace petace
