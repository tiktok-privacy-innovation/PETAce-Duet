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

#include <thread>

#include "gtest/gtest.h"
#include "print_utils.h"

#include "network/net_factory.h"
#include "network/net_socket.h"
#include "network/network.h"

#include "duet/st_generator/st_generator.h"
#include "duet/util/common.h"
#include "duet/util/secret_shared_shuffle.h"

namespace petace {
namespace duet {

class ShuffleTest : public ::testing::Test {
public:
    bool is_block_equal(block& first, block& second) {
        return (first[0] ^ first[1] ^ second[0] ^ second[1]) == 0;
    }

    void st_generator(std::size_t party_id, Permutation& p, Matrix<std::int64_t>& delta, Matrix<std::int64_t>& a,
            Matrix<std::int64_t>& b) {
        std::size_t rows = 12;
        std::size_t cols = 3;
        petace::network::NetParams net_params;
        if (party_id == 0) {
            net_params.remote_addr = "127.0.0.1";
            net_params.remote_port = 8890;
            net_params.local_port = 8891;
        } else {
            net_params.remote_addr = "127.0.0.1";
            net_params.remote_port = 8891;
            net_params.local_port = 8890;
        }

        auto net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params);

        block seed = read_block_from_dev_urandom();
        std::shared_ptr<OTGenerator> ot = std::make_shared<OTGenerator>(party_id, seed);
        ot->initialize(net);
        STGenerator st(party_id, ot);
        if (party_id == 0) {
            st.get_st(rows, cols, net, p, delta);
        } else {
            st.get_st(rows, cols, net, a, b);
        }
    }

private:
    block fixed_key_ = {0, 1};
};

TEST_F(ShuffleTest, double_prg) {
    // KAT test
    PRGSeed a = {0, 1};
    PRGSeed b;
    PRGSeed c;
    DoublePrg::gen(a, b, c);
    PRGSeed aim_b = {-3054159662337734857LL, 2941533170686271127LL};
    PRGSeed aim_c = {1203591567212739478LL, -9079467822864142152LL};
    EXPECT_EQ(is_block_equal(b, aim_b), true);
    EXPECT_EQ(is_block_equal(c, aim_c), true);
}

TEST_F(ShuffleTest, opv_test) {
    GGMTreeNode zero = {0, 0};
    std::size_t n = 12;
    std::size_t aim_layer = 4;
    std::size_t pos = 9;
    ObliviousPuncturedVector test(n);
    std::size_t ot_num = test.get_ot_num();
    EXPECT_EQ(aim_layer, ot_num);

    std::vector<OTChoice> choice;
    std::vector<OTChoice> aim_choice = {0, 1, 1, 0};
    test.active_phase_1(pos, choice);
    EXPECT_EQ(choice, aim_choice);

    std::vector<std::vector<GGMTreeNode>> all_level_sums;
    std::vector<GGMTreeNode> leaves_passive;
    test.passive_phase_1(all_level_sums, leaves_passive);
    // mock ot
    std::vector<GGMTreeNode> needs_level_sums(ot_num);
    for (std::size_t i = 0; i < ot_num; ++i) {
        if (choice[i] == 0) {
            needs_level_sums[i] = all_level_sums[i][0];
        } else {
            needs_level_sums[i] = all_level_sums[i][1];
        }
    }

    std::vector<GGMTreeNode> leaves_active;
    test.active_phase_2(pos, needs_level_sums, leaves_active);

    EXPECT_EQ(is_block_equal(leaves_active[9], zero), true);
    for (std::size_t i = 0; i < n; ++i) {
        if (i != pos) {
            EXPECT_EQ(is_block_equal(leaves_active[i], leaves_passive[i]), true);
        } else {
            EXPECT_EQ(is_block_equal(leaves_active[i], leaves_passive[i]), false);
        }
    }
}

TEST_F(ShuffleTest, secert_shared_shuffle_matrix_test) {
    std::size_t rows = 128;
    GGMTreeNode zero = {0, 0};
    Permutation p(rows);
    ShareTranslation test(rows);
    std::vector<std::vector<OTChoice>> all_choice;
    test.active_phase_1(p, all_choice);
    EXPECT_EQ(all_choice.size(), rows);
    EXPECT_EQ(all_choice[0].size(), 7);

    std::vector<std::vector<std::vector<GGMTreeNode>>> all_levels_sums;
    std::vector<std::vector<GGMTreeNode>> passive_all_leaves_matrix;
    test.passive_phase_1(all_levels_sums, passive_all_leaves_matrix);

    // mock_ot
    std::vector<std::vector<GGMTreeNode>> all_need_levels_sums(rows, std::vector<GGMTreeNode>(all_choice[0].size()));
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < all_choice[0].size(); ++j) {
            if (all_choice[i][j] == 0) {
                all_need_levels_sums[i][j] = all_levels_sums[i][j][0];
            } else {
                all_need_levels_sums[i][j] = all_levels_sums[i][j][1];
            }
        }
    }

    std::vector<std::vector<GGMTreeNode>> active_all_leaves_matrix;
    test.active_phase_2(p, all_need_levels_sums, active_all_leaves_matrix);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < rows; ++j) {
            if (p[i] != j) {
                EXPECT_EQ(is_block_equal(active_all_leaves_matrix[i][j], passive_all_leaves_matrix[i][j]), true);
            } else {
                EXPECT_EQ(is_block_equal(active_all_leaves_matrix[i][j], zero), true);
            }
        }
    }
}

TEST_F(ShuffleTest, secert_shared_shuffle_arith_test) {
    // gen st
    std::size_t rows = 12;
    std::size_t cols = 3;
    std::size_t matrix_size = rows * cols;
    Permutation p(rows);
    ShareTranslation test(rows);
    std::vector<std::vector<OTChoice>> all_choice;
    test.active_phase_1(p, all_choice);
    EXPECT_EQ(all_choice.size(), rows);
    EXPECT_EQ(all_choice[0].size(), 4);

    std::vector<std::vector<std::vector<GGMTreeNode>>> all_levels_sums;
    Matrix<std::int64_t> a;
    Matrix<std::int64_t> b;
    test.passive_phase_1<std::int64_t>(cols, all_levels_sums, a, b);

    // mock_ot
    std::vector<std::vector<GGMTreeNode>> all_need_levels_sums(rows, std::vector<GGMTreeNode>(all_choice[0].size()));
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < all_choice[0].size(); ++j) {
            if (all_choice[i][j] == 0) {
                all_need_levels_sums[i][j] = all_levels_sums[i][j][0];
            } else {
                all_need_levels_sums[i][j] = all_levels_sums[i][j][1];
            }
        }
    }

    Matrix<std::int64_t> delta;
    test.active_phase_2<std::int64_t>(cols, p, all_need_levels_sums, delta);

    Matrix<std::int64_t> pa = p.permute(a);
    EXPECT_EQ(pa - b, delta);

    // test for secert shared shuffle;
    Matrix<std::int64_t> x;
    x.resize(rows, cols);
    for (std::size_t i = 0; i < matrix_size; ++i) {
        x(i) = i;
    }
    SecretSharedShuffle ss_shuffle_test;
    Matrix<std::int64_t> x_sub_a;
    ss_shuffle_test.passive_phase_1(x, a, x_sub_a);
    Matrix<std::int64_t> active_share;
    Matrix<std::int64_t> passive_share = b;
    ss_shuffle_test.active_phase_1(p, x_sub_a, delta, active_share);
    Matrix<std::int64_t> px = p.permute(x);
    EXPECT_EQ(px, active_share + passive_share);
}

TEST_F(ShuffleTest, secert_shared_shuffle_bool_test) {
    // gen st
    std::size_t rows = 12;
    std::size_t cols = 3;
    Permutation p(rows);
    ShareTranslation test(rows);
    std::vector<std::vector<OTChoice>> all_choice;
    test.active_phase_1(p, all_choice);
    EXPECT_EQ(all_choice.size(), rows);
    EXPECT_EQ(all_choice[0].size(), 4);

    std::vector<std::vector<std::vector<GGMTreeNode>>> all_levels_sums;
    Matrix<std::int64_t> a;
    Matrix<std::int64_t> b;
    test.passive_phase_1<std::int64_t>(cols, all_levels_sums, a, b, true);

    // mock_ot
    std::vector<std::vector<GGMTreeNode>> all_need_levels_sums(rows, std::vector<GGMTreeNode>(all_choice[0].size()));
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < all_choice[0].size(); ++j) {
            if (all_choice[i][j] == 0) {
                all_need_levels_sums[i][j] = all_levels_sums[i][j][0];
            } else {
                all_need_levels_sums[i][j] = all_levels_sums[i][j][1];
            }
        }
    }

    Matrix<std::int64_t> delta;
    test.active_phase_2<std::int64_t>(cols, p, all_need_levels_sums, delta, true);

    Matrix<std::int64_t> pa = p.permute(a);

    Matrix<std::int64_t> pa_xor_b;
    pa_xor_b.resize(rows, cols);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            pa_xor_b(i * cols + j) = pa(i * cols + j) ^ b(i * cols + j);
        }
    }

    EXPECT_EQ(pa_xor_b, delta);

    // test for secert shared shuffle;
}

TEST_F(ShuffleTest, st_generator_test) {
    std::vector<std::thread> threads;
    Permutation p(12);
    Matrix<std::int64_t> delta;
    Matrix<std::int64_t> a;
    Matrix<std::int64_t> b;

    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { st_generator(i, p, delta, a, b); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    Matrix<std::int64_t> pa = p.permute(a);
    EXPECT_EQ(pa - b, delta);
}

}  // namespace duet
}  // namespace petace
