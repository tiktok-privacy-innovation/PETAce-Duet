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
#include "test_utils.h"

#include "duet/vm/duet_vm.h"
#include "duet/vm/instruction.h"
#include "duet/vm/operator.h"
#include "duet/vm/register.h"

namespace petace {
namespace duet {

class DuetVMTest : public ::testing::Test {
public:
    void SetUp() {
        net_params0.remote_addr = "127.0.0.1";
        net_params0.remote_port = 8890;
        net_params0.local_port = 8891;

        net_params1.remote_addr = "127.0.0.1";
        net_params1.remote_port = 8891;
        net_params1.local_port = 8890;
    }

    DuetVM init_vm(std::size_t party_id) {
        petace::network::NetParams net_params;
        if (party_id == 0) {
            net_params = net_params0;
        } else {
            net_params = net_params1;
        }
        DuetVM vm(net_params, petace::network::NetScheme::SOCKET, party_id);
        return vm;
    }

    void prepare_private_data(std::size_t party_id, DuetVM& vm, std::vector<RegisterAddress>& reg_addr,
            std::size_t rows, std::size_t cols) {
        std::shared_ptr<PrivateMatrix<double>> a = std::make_shared<PrivateMatrix<double>>(rows, cols, 0);
        std::shared_ptr<PrivateMatrix<double>> b = std::make_shared<PrivateMatrix<double>>(rows, cols, 0);
        if (party_id == 0) {
            for (std::size_t i = 0; i < 5; ++i) {
                (*a)(i) = static_cast<double>(i);
            }
            for (std::size_t i = 0; i < 5; ++i) {
                (*b)(i) = static_cast<double>(i * 10);
            }
        }
        reg_addr.push_back(vm.set_data(a));
        reg_addr.push_back(vm.set_data(b));
    }

    void share_private_data(
            DuetVM& vm, std::vector<RegisterAddress>& private_reg_addr, std::vector<RegisterAddress>& share_reg_addr) {
        Instruction share_inst({"share", "pdm", "am"});
        for (size_t i = 0; i < private_reg_addr.size(); ++i) {
            std::vector<RegisterAddress> tmp;
            share_reg_addr.push_back(vm.new_data<ArithMatrix>());
            tmp.push_back(private_reg_addr[i]);
            tmp.push_back(share_reg_addr[i]);
            vm.exec_code(share_inst, tmp);
        }
    }

    void delete_data(DuetVM& vm, std::vector<RegisterAddress>& reg_addr) {
        for (std::size_t i = 0; i < reg_addr.size(); ++i) {
            vm.delete_data(reg_addr[i]);
        }
    }

    void test_set_get_data(std::size_t party_id) {
        DuetVM vm = init_vm(party_id);
        RegisterAddress addr_0 = vm.new_data<double>();
        const std::shared_ptr<double>& ptr_0 = vm.get_data<double>(addr_0);
        EXPECT_EQ(*ptr_0, 0);
        vm.delete_data(addr_0);
        EXPECT_EQ(true, vm.is_registr_empty());
        RegisterAddress addr_1 = vm.set_data<double>(std::make_shared<double>(1.1));
        const std::shared_ptr<double>& ptr_1 = vm.get_data<double>(addr_1);
        EXPECT_EQ(*ptr_1, 1.1);
        vm.delete_data(addr_1);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_airth(std::size_t party_id, std::string test_type) {
        DuetVM vm = init_vm(party_id);

        std::vector<RegisterAddress> private_reg_addr;
        prepare_private_data(party_id, vm, private_reg_addr, 5, 1);
        // now [0] is 0's private, [1] is 1's private

        std::vector<RegisterAddress> share_reg_addr;
        share_private_data(vm, private_reg_addr, share_reg_addr);
        // now [0][1] are both share

        share_reg_addr.push_back(vm.new_data<ArithMatrix>());
        // now [2] is the ret
        Instruction airth_inst({test_type, "am", "am", "am"});
        vm.exec_code(airth_inst, share_reg_addr);

        Instruction reveal_inst({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> tmp;
        private_reg_addr.push_back(vm.set_data(std::make_shared<PrivateMatrix<double>>(0)));
        tmp.push_back(share_reg_addr[2]);
        tmp.push_back(private_reg_addr[2]);
        vm.exec_code(reveal_inst, tmp);

        if (party_id == 0) {
            Matrix<double> input_0 = vm.get_data<PrivateMatrix<double>>(private_reg_addr[0])->matrix();
            Matrix<double> input_1 = vm.get_data<PrivateMatrix<double>>(private_reg_addr[1])->matrix();
            Matrix<double> aim_ret;
            if (test_type == "add") {
                aim_ret = input_0 + input_1;
            } else if (test_type == "sub") {
                aim_ret = input_0 - input_1;
            } else if (test_type == "mul") {
                aim_ret = input_0.cwiseProduct(input_1);
            }

            Matrix<double> real_ret = vm.get_data<PrivateMatrix<double>>(private_reg_addr[2])->matrix();

            EXPECT_EQ(true, is_equal_plain_matrix(aim_ret, real_ret, 0.001));
        }

        delete_data(vm, private_reg_addr);
        delete_data(vm, share_reg_addr);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_bool(std::size_t party_id, std::string test_type) {
        DuetVM vm = init_vm(party_id);

        std::vector<RegisterAddress> private_reg_addr;
        prepare_private_data(party_id, vm, private_reg_addr, 5, 1);
        // now [0] is 0's private, [1] is 1's private

        std::vector<RegisterAddress> share_reg_addr;
        share_private_data(vm, private_reg_addr, share_reg_addr);
        // now [0][1] are both share

        share_reg_addr.push_back(vm.new_data<BoolMatrix>());
        // now [2] is the ret
        Instruction bool_inst({test_type, "am", "am", "bm"});

        vm.exec_code(bool_inst, share_reg_addr);

        Instruction reveal_inst({"reveal", "bm", "pbm"});
        std::vector<RegisterAddress> tmp;
        private_reg_addr.push_back(vm.set_data(std::make_shared<PrivateMatrixBool>(0)));
        tmp.push_back(share_reg_addr[2]);
        tmp.push_back(private_reg_addr[2]);
        vm.exec_code(reveal_inst, tmp);

        if (party_id == 0) {
            Matrix<double> input_0 = vm.get_data<PrivateMatrix<double>>(private_reg_addr[0])->matrix();
            Matrix<double> input_1 = vm.get_data<PrivateMatrix<double>>(private_reg_addr[1])->matrix();
            Matrix<std::int64_t> aim_ret(5, 1);
            for (size_t i = 0; i < 5; ++i) {
                if (test_type == "gt") {
                    aim_ret(i) = input_0(i) > input_1(i) ? 1 : 0;
                } else if (test_type == "lt") {
                    aim_ret(i) = input_0(i) < input_1(i) ? 1 : 0;
                } else if (test_type == "ge") {
                    aim_ret(i) = input_0(i) >= input_1(i) ? 1 : 0;
                } else if (test_type == "eq") {
                    aim_ret(i) = (input_0(i) == input_1(i)) ? 1 : 0;
                }
            }

            Matrix<std::int64_t> real_ret = vm.get_data<PrivateMatrixBool>(private_reg_addr[2])->matrix();

            EXPECT_EQ(true, is_equal_plain_matrix(aim_ret, real_ret));
        }

        delete_data(vm, private_reg_addr);
        delete_data(vm, share_reg_addr);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    petace::network::NetParams net_params0;
    petace::network::NetParams net_params1;
};

TEST(Register, case0) {
    Register r;
    std::shared_ptr<Matrix<double>> in = std::make_shared<Matrix<double>>(5, 1);
    std::size_t matrix_size = 5;
    for (std::size_t i = 0; i < matrix_size; ++i) {
        (*in)(i) = static_cast<double>(i);
    }
    RegisterAddress address = r.set_data(in);
    const std::shared_ptr<Matrix<double>>& out0 = r.get_data<Matrix<double>>(address);
    (*in)(0) = 233;
    EXPECT_EQ(233, (*out0)(0));
    r.delete_data(address);
    EXPECT_THROW(auto out1 = r.get_data<Matrix<double>>(address), std::out_of_range);
}

TEST(Instruction, case0) {
    Instruction inst0({"123"});
    Instruction inst1({"123"});
    Instruction inst2({"1234"});
    EXPECT_EQ(inst0, inst1);
    EXPECT_EQ((inst0 == inst2), false);
    std::unordered_map<Instruction, int> map;
    map[inst0] = 0;
    EXPECT_EQ(map.at(inst0), 0);
    EXPECT_THROW(map.at(inst2), std::out_of_range);
}

TEST_F(DuetVMTest, test_set_get_data) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_set_get_data(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_add) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_airth(i, "add"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_sub) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_airth(i, "sub"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_mul) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_airth(i, "mul"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_greater) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_bool(i, "gt"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_less) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_bool(i, "lt"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_greater_eq) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_bool(i, "ge"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_equal) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_bool(i, "eq"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

}  // namespace duet
}  // namespace petace
