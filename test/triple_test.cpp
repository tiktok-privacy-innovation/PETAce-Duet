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

#include <stdlib.h>
#include <unistd.h>

#include <memory>
#include <thread>
#include <utility>

#include "gtest/gtest.h"

#include "network/net_factory.h"
#include "network/net_socket.h"
#include "network/network.h"
#include "solo/prng.h"

#include "duet/duet.h"
#include "duet/util/common.h"

namespace petace {
namespace duet {

class TripleTest : public ::testing::Test {
public:
    void SetUp() {
        net_params0.remote_addr = "127.0.0.1";
        net_params0.remote_port = 8890;
        net_params0.local_port = 8891;

        net_params1.remote_addr = "127.0.0.1";
        net_params1.remote_port = 8891;
        net_params1.local_port = 8890;
    }

    void test_triple(bool is_sender, block& common_seed, petace::duet::TripleScheme triple_scheme,
            std::vector<std::vector<std::int64_t>>& triples) {
        std::shared_ptr<petace::network::Network> net;
        if (is_sender) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        // set ot generator

        std::shared_ptr<PRNG> rand_generator = std::make_shared<PRNG>(common_seed);
        block seed = _mm_set_epi64x(
                rand_generator->get_unique_rand<std::int64_t>(), rand_generator->get_unique_rand<std::int64_t>());
        std::shared_ptr<ObliviousTransfer> ot_generator = std::make_shared<ObliviousTransfer>(is_sender, seed);
        ot_generator->initialize(net);
        std::shared_ptr<ArithmeticTriple> rand_arith_triple_generator =
                TripleFactory::get_instance().build(triple_scheme, net, is_sender, rand_generator, ot_generator);

        for (std::size_t i = 0; i < 100; i++) {
            auto triple = rand_arith_triple_generator->get_rand_triple(net);
            triples.emplace_back(triple);
        }
        return;
    }

public:
    std::thread t_[2];
    petace::network::NetParams net_params0;
    petace::network::NetParams net_params1;
};

TEST_F(TripleTest, fhe_triple_test) {
    std::vector<std::thread> threads;
    std::vector<std::vector<std::vector<std::int64_t>>> triples(2);
    block common_seed = read_block_from_dev_urandom();
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_triple(i, common_seed, TripleScheme::FHE, triples[i]); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    for (std::size_t i = 0; i < 100; ++i) {
        EXPECT_EQ((triples[0][i][0] + triples[1][i][0]) * (triples[0][i][1] + triples[1][i][1]),
                (triples[0][i][2] + triples[1][i][2]));
    }
}

TEST_F(TripleTest, ot_triple_test) {
    std::vector<std::thread> threads;
    std::vector<std::vector<std::vector<std::int64_t>>> triples(2);
    block common_seed = read_block_from_dev_urandom();
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_triple(i, common_seed, TripleScheme::OT, triples[i]); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    for (std::size_t i = 0; i < 100; ++i) {
        EXPECT_EQ((triples[0][i][0] + triples[1][i][0]) * (triples[0][i][1] + triples[1][i][1]),
                (triples[0][i][2] + triples[1][i][2]));
    }
}

}  // namespace duet
}  // namespace petace
