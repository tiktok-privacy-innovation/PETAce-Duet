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

#include "gtest/gtest.h"
#include "print_utils.h"

#include "duet/util/matrix.h"
#include "duet/util/permutation.h"

namespace petace {
namespace duet {

TEST(UtilTest, permute_and_inv) {
    std::vector<std::size_t> test_vector = {0, 1, 2, 3, 4};
    Permutation p(5);

    auto ret0 = p.permute(test_vector);
    EXPECT_EQ(ret0, p.data());

    Permutation inv_p = p.inverse();
    auto ret1 = inv_p.permute(ret0);
    EXPECT_EQ(ret1, test_vector);

    auto ret2 = p.inverse(ret0);
    EXPECT_EQ(ret2, test_vector);
}

TEST(UtilTest, permute_index) {
    Permutation p(5);
    EXPECT_EQ(p[0], p.data()[0]);
    EXPECT_THROW(p[5], std::invalid_argument);
}

TEST(UtilTest, permute_matrix) {
    Permutation p0(5);
    Matrix<std::size_t> test_matrix0(5, 1);
    std::size_t matrix_size = 5;
    for (std::size_t i = 0; i < matrix_size; ++i) {
        test_matrix0(i) = i;
    }
    auto ret0 = p0.permute(test_matrix0);
    std::vector<std::size_t> ret_vector0(5);
    for (std::size_t i = 0; i < ret_vector0.size(); ++i) {
        ret_vector0[i] = ret0(i);
    }
    EXPECT_EQ(ret_vector0, p0.data());

    ArithMatrix test_matrix1(5, 1);
    for (std::size_t i = 0; i < matrix_size; ++i) {
        test_matrix1(i) = i;
    }
    auto ret1 = p0.permute(test_matrix1);
    std::vector<std::size_t> ret_vector1(5);
    for (std::size_t i = 0; i < ret_vector1.size(); ++i) {
        ret_vector1[i] = ret1(i);
    }
    EXPECT_EQ(ret_vector1, p0.data());
}

TEST(UtilTest, permute_creat) {
    std::vector<std::size_t> test_p0 = {2, 1, 0, 3, 4};
    Permutation p0(test_p0);
    EXPECT_EQ(test_p0, p0.data());

    std::vector<std::size_t> test_p1 = {2, 1, 0, 3, 5};
    EXPECT_THROW(Permutation p1(test_p1), std::invalid_argument);
}

TEST(UtilTest, permute_combine) {
    std::vector<std::size_t> test_p0 = {0, 1, 2, 3};
    std::vector<std::size_t> test_p1 = {0, 1, 2, 3, 4, 5, 6, 7};
    Permutation p0(test_p0);
    Permutation p1(test_p0);
    Permutation p2 = p0.combine(p1);
    EXPECT_EQ(p2.size(), 8);
    EXPECT_EQ(p2.data(), test_p1);
}

TEST(UtilTest, private_matrix) {
    PrivateMatrix<std::int64_t> pm(2, 2, 1);
    EXPECT_EQ(1, pm.party_id());
    Matrix<std::int64_t> m(2, 2);
    pm(0) = 233;
    m(0) = 456;
    EXPECT_EQ(pm(0), 233);
    EXPECT_EQ(pm.matrix()(0), 233);
    EXPECT_EQ(pm(0, 0), 233);
    EXPECT_EQ(pm.matrix()(0, 0), 233);
    pm.matrix()(0) = 789;
    EXPECT_EQ(pm(0), 789);
    EXPECT_EQ(pm.matrix()(0), 789);
    EXPECT_EQ(pm(0, 0), 789);
    EXPECT_EQ(pm.matrix()(0, 0), 789);
    pm.matrix() = m;
    EXPECT_EQ(pm(0), 456);
    EXPECT_EQ(pm.matrix(), m);
    PrivateMatrix<std::int64_t> sub_pm;
    matrix_block(pm, sub_pm, 0, 0, 1, 1, 0);
    EXPECT_EQ(sub_pm.size(), 0);
    matrix_block(pm, sub_pm, 0, 0, 1, 1, 1);
    EXPECT_EQ(sub_pm.size(), 1);
    EXPECT_EQ(sub_pm.party_id(), 1);
    EXPECT_EQ(sub_pm(0), 456);
    Matrix<std::int64_t> vm(4, 2);
    vm << m, m;
    PrivateMatrix<std::int64_t> vpm;
    vstack(pm, pm, vpm, 0);
    EXPECT_EQ(vpm.size(), 0);
    vstack(pm, pm, vpm, 1);
    EXPECT_EQ(vpm.size(), 8);
    EXPECT_EQ(vpm.matrix(), vm);
    Matrix<std::int64_t> hm(2, 4);
    hm << m, m;
    PrivateMatrix<std::int64_t> hpm;
    hstack(pm, pm, hpm, 0);
    EXPECT_EQ(hpm.size(), 0);
    hstack(pm, pm, hpm, 1);
    EXPECT_EQ(hpm.size(), 8);
    EXPECT_EQ(hpm.matrix(), hm);
}

TEST(UtilTest, private_matrix_index) {
    PrivateMatrix<std::int64_t> pm_0(0);
    Matrix<std::int64_t> m_0(2, 2);
    m_0 << 0, 0, 1, 1;
    pm_0.index_like(2, 2, 0, 0);
    EXPECT_EQ(pm_0.matrix(), m_0);

    PrivateMatrix<std::int64_t> pm_1(0);
    Matrix<std::int64_t> m_1(2, 2);
    m_1 << 0, 1, 0, 1;
    pm_1.index_like(2, 2, 1, 0);
    EXPECT_EQ(pm_1.matrix(), m_1);

    PrivateMatrix<std::int64_t> pm_2(0);
    Matrix<std::int64_t> m_2(2, 2);
    m_2 << 0, 1, 2, 3;
    pm_2.index_like(2, 2, 0);
    EXPECT_EQ(pm_2.matrix(), m_2);
}

TEST(UtilTest, public_matrix) {
    PublicMatrix<std::int64_t> cm(2, 2);
    Matrix<std::int64_t> m(2, 2);
    cm(0) = 233;
    m(0) = 456;
    EXPECT_EQ(cm(0), 233);
    EXPECT_EQ(cm.matrix()(0), 233);
    cm.matrix()(0) = 789;
    EXPECT_EQ(cm(0), 789);
    EXPECT_EQ(cm.matrix()(0), 789);
    cm.matrix() = m;
    EXPECT_EQ(cm(0), 456);
    EXPECT_EQ(cm.matrix(), m);
    Matrix<std::int64_t> vm(4, 2);
    vm << m, m;
    PublicMatrix<std::int64_t> vcm;
    vstack(cm, cm, vcm);
    EXPECT_EQ(vcm.size(), 8);
    EXPECT_EQ(vcm.matrix(), vm);
    Matrix<std::int64_t> hm(2, 4);
    hm << m, m;
    PublicMatrix<std::int64_t> hcm;
    hstack(cm, cm, hcm);
    EXPECT_EQ(hcm.size(), 8);
    EXPECT_EQ(hcm.matrix(), hm);
}

}  // namespace duet
}  // namespace petace
