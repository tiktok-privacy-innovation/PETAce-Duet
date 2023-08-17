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

#include "duet/ot_generator/ot_generator.h"

class OtTest : public ::testing::Test {
public:
    void SetUp() {
        net_params0.remote_addr = "127.0.0.1";
        net_params0.remote_port = 8890;
        net_params0.local_port = 8891;

        net_params1.remote_addr = "127.0.0.1";
        net_params1.remote_port = 8891;
        net_params1.local_port = 8890;
    }

    void rand_ot(bool is_sender) {
        std::shared_ptr<petace::network::Network> net;
        if (is_sender) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }

        srandom(static_cast<std::uint32_t>(time(nullptr)));
        petace::duet::OTGenerator ot(is_sender, _mm_set_epi64x(random(), random()));
        ot.initialize(net);

        if (is_sender) {
            ot.get_random_ot(net, send_msg_);
            net->send_data(send_msg_.data(), send_msg_.size() * sizeof(std::int64_t));
        } else {
            ot.get_random_ot(net, choice_, recv_msg_);
            send_msg_.resize(2);
            net->recv_data(send_msg_.data(), send_msg_.size() * sizeof(std::int64_t));
        }
    }

    void rand_2ot(bool is_sender) {
        std::shared_ptr<petace::network::Network> net;
        if (is_sender) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }

        srandom(static_cast<std::uint32_t>(time(nullptr)));
        petace::duet::OTGenerator ot(is_sender, _mm_set_epi64x(random(), random()));
        ot.initialize(net);

        if (is_sender) {
            ot.get_random_ot(net, send_msg_);
        } else {
            ot.get_random_ot(net, choice_, recv_msg_);
        }

        if (is_sender) {
            ot.get_random_ot(net, choice_, recv_msg_);
            net->send_data(&recv_msg_, sizeof(std::int64_t));
            net->send_data(&choice_, sizeof(std::int8_t));
        } else {
            ot.get_random_ot(net, send_msg_);
            net->recv_data(&recv_msg_, sizeof(std::int64_t));
            net->recv_data(&choice_, sizeof(std::int8_t));
        }
    }

    void choices_ot(bool is_sender) {
        std::shared_ptr<petace::network::Network> net;
        if (is_sender) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }

        srandom(static_cast<std::uint32_t>(time(nullptr)));
        petace::duet::OTGenerator ot(is_sender, _mm_set_epi64x(random(), random()));
        ot.initialize(net);

        std::size_t ot_len = 16;
        choices_.resize(ot_len);
        send_msgs_.resize(ot_len);
        for (std::size_t i = 0; i < ot_len; i++) {
            send_msgs_[i].resize(2);
            choices_[i] = random() & 0x1;
            send_msgs_[i][0] = random();
            send_msgs_[i][1] = random();
        }

        if (is_sender) {
            ot.get_standard_ot(net, send_msgs_);
            for (std::size_t i = 0; i < ot_len; i++) {
                net->send_data(&send_msgs_[i][0], 2 * sizeof(std::int64_t));
            }
        } else {
            ot.get_standard_ot(net, choices_, recv_msgs_);
            for (std::size_t i = 0; i < ot_len; i++) {
                net->recv_data(&send_msgs_[i][0], 2 * sizeof(std::int64_t));
            }
        }
    }

    void correlated_ot(bool is_sender) {
        std::shared_ptr<petace::network::Network> net;
        if (is_sender) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }

        srandom(static_cast<std::uint32_t>(time(nullptr)));
        petace::duet::OTGenerator ot(is_sender, _mm_set_epi64x(random(), random()));
        ot.initialize(net);

        std::size_t ot_len = 16;
        choices_.resize(ot_len);
        send_msgs_.resize(ot_len);
        for (std::size_t i = 0; i < ot_len; i++) {
            send_msgs_[i].resize(2);
        }

        if (is_sender) {
            ot.get_correlated_ot(net, ot_len, delta_, send_msgs_);
            for (std::size_t i = 0; i < ot_len; i++) {
                net->send_data(&send_msgs_[i][0], 2 * sizeof(std::int64_t));
            }
        } else {
            ot.get_correlated_ot(net, ot_len, choices_, recv_msgs_);
            for (std::size_t i = 0; i < ot_len; i++) {
                net->recv_data(&send_msgs_[i][0], 2 * sizeof(std::int64_t));
            }
        }
    }

    void one_of_n_ot(bool is_sender) {
        std::shared_ptr<petace::network::Network> net;
        if (is_sender) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }

        srandom(static_cast<std::uint32_t>(time(nullptr)));
        petace::duet::OTGenerator ot(is_sender, _mm_set_epi64x(random(), random()));
        ot.initialize(net);

        std::size_t ot_len = 32;
        std::size_t select_size = 8;
        choices_.resize(ot_len);
        send_msgs_.resize(ot_len);
        for (std::size_t i = 0; i < ot_len; i++) {
            send_msgs_[i].resize(select_size);
            choices_[i] = random() & 0x7;
            for (std::size_t j = 0; j < select_size; j++) {
                send_msgs_[i][j] = random();
            }
        }

        if (is_sender) {
            ot.get_nch1_ot(net, select_size, send_msgs_);
            for (std::size_t i = 0; i < ot_len; i++) {
                net->send_data(&send_msgs_[i][0], select_size * sizeof(std::int64_t));
            }
        } else {
            ot.get_nch1_ot(net, select_size, choices_, recv_msgs_);
            for (std::size_t i = 0; i < ot_len; i++) {
                net->recv_data(&send_msgs_[i][0], select_size * sizeof(std::int64_t));
            }
        }
    }

public:
    std::thread t_[2];
    std::vector<std::int64_t> send_msg_;
    std::int64_t recv_msg_;
    std::int8_t choice_;
    std::vector<std::int8_t> choices_;
    std::vector<std::vector<std::int64_t>> send_msgs_;
    std::vector<std::int64_t> recv_msgs_;
    std::int64_t delta_ = 0x12345678;
    petace::network::NetParams net_params0;
    petace::network::NetParams net_params1;
};

TEST_F(OtTest, rand_ot) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        rand_ot(true);
        exit(EXIT_SUCCESS);
    } else {
        rand_ot(false);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        EXPECT_EQ(send_msg_[choice_], recv_msg_);
    }
}

TEST_F(OtTest, rand_2ot) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        rand_2ot(true);
        exit(EXIT_SUCCESS);
    } else {
        rand_2ot(false);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        EXPECT_EQ(send_msg_[choice_], recv_msg_);
    }
}

TEST_F(OtTest, choices_ot) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        choices_ot(true);
        exit(EXIT_SUCCESS);
    } else {
        choices_ot(false);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        for (std::size_t i = 0; i < 16; i++) {
            EXPECT_EQ(send_msgs_[i][choices_[i]], recv_msgs_[i]);
        }
    }
}

TEST_F(OtTest, correlated_ot) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        correlated_ot(true);
        exit(EXIT_SUCCESS);
    } else {
        correlated_ot(false);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        for (std::size_t i = 0; i < 16; i++) {
            EXPECT_EQ(send_msgs_[i][choices_[i]], recv_msgs_[i]);
            EXPECT_EQ(send_msgs_[i][1] - send_msgs_[i][0], delta_);
        }
    }
}

TEST_F(OtTest, one_of_n_ot) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        one_of_n_ot(true);
        exit(EXIT_SUCCESS);
    } else {
        one_of_n_ot(false);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        for (std::size_t i = 0; i < 32; i++) {
            EXPECT_EQ(send_msgs_[i][choices_[i]], recv_msgs_[i]);
        }
    }
}
