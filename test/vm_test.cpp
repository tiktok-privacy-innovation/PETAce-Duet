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

#include <random>
#include <thread>

#include "gtest/gtest.h"
#include "print_utils.h"
#include "test_utils.h"

#include "network/net_factory.h"
#include "network/network.h"

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
        std::shared_ptr<petace::network::Network> net =
                petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params);

        DuetVM vm(net, party_id);
        return vm;
    }

    void prepare_private_data(std::size_t party_id, DuetVM& vm, std::vector<RegisterAddress>& reg_addr,
            std::size_t rows, std::size_t cols) {
        std::shared_ptr<PrivateMatrix<double>> a = std::make_shared<PrivateMatrix<double>>(rows, cols, 0);
        std::shared_ptr<PrivateMatrix<double>> b = std::make_shared<PrivateMatrix<double>>(rows, cols, 0);
        if (party_id == 0) {
            for (std::size_t i = 0; i < rows * cols; ++i) {
                (*a)(i) = static_cast<double>(i);
            }
            for (std::size_t i = 0; i < rows * cols; ++i) {
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

    void test_mat_airth(std::size_t party_id, std::string test_type) {
        DuetVM vm = init_vm(party_id);

        std::vector<RegisterAddress> private_reg_addr;
        prepare_private_data(party_id, vm, private_reg_addr, 4, 4);
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
            if (test_type == "mat_mul") {
                aim_ret = input_0 * input_1;
            }

            Matrix<double> real_ret = vm.get_data<PrivateMatrix<double>>(private_reg_addr[2])->matrix();

            EXPECT_EQ(true, is_equal_plain_matrix(aim_ret, real_ret, 0.01));
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

        Instruction not_inst({"not", "bm", "bm"});
        std::vector<RegisterAddress> not_reg_addr;
        not_reg_addr.push_back(share_reg_addr[2]);
        not_reg_addr.push_back(vm.new_data<BoolMatrix>());
        vm.exec_code(not_inst, not_reg_addr);

        Instruction reveal_inst({"reveal", "bm", "pbm"});
        std::vector<RegisterAddress> tmp;
        private_reg_addr.push_back(vm.set_data(std::make_shared<PrivateMatrixBool>(0)));
        tmp.push_back(not_reg_addr[1]);
        tmp.push_back(private_reg_addr[2]);
        vm.exec_code(reveal_inst, tmp);

        if (party_id == 0) {
            Matrix<double> input_0 = vm.get_data<PrivateMatrix<double>>(private_reg_addr[0])->matrix();
            Matrix<double> input_1 = vm.get_data<PrivateMatrix<double>>(private_reg_addr[1])->matrix();
            Matrix<std::int64_t> aim_ret(5, 1);
            for (size_t i = 0; i < 5; ++i) {
                if (test_type == "gt") {
                    aim_ret(i) = input_0(i) > input_1(i) ? 0 : 1;
                } else if (test_type == "lt") {
                    aim_ret(i) = input_0(i) < input_1(i) ? 0 : 1;
                } else if (test_type == "ge") {
                    aim_ret(i) = input_0(i) >= input_1(i) ? 0 : 1;
                } else if (test_type == "eq") {
                    aim_ret(i) = (input_0(i) == input_1(i)) ? 0 : 1;
                }
            }

            Matrix<std::int64_t> real_ret = vm.get_data<PrivateMatrixBool>(private_reg_addr[2])->matrix();
            EXPECT_EQ(true, is_equal_plain_matrix(aim_ret, real_ret));
        }

        delete_data(vm, private_reg_addr);
        delete_data(vm, share_reg_addr);
        delete_data(vm, not_reg_addr);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_arith_bool(std::size_t party_id, std::string test_type) {
        DuetVM vm = init_vm(party_id);

        PrivateMatrixBool private_bool_matrix0(4, 3, 0);
        for (std::size_t i = 0; i < private_bool_matrix0.size(); i++) {
            private_bool_matrix0(i) = i % 2;
        }
        Instruction share_inst({"share", "pbm", "bm"});
        std::vector<RegisterAddress> share_reg_addr0;
        share_reg_addr0.push_back(vm.set_data(std::make_shared<PrivateMatrixBool>(private_bool_matrix0)));
        share_reg_addr0.push_back(vm.new_data<BoolMatrix>());
        vm.exec_code(share_inst, share_reg_addr0);

        PrivateMatrixBool private_bool_matrix1(4, 3, 0);
        for (std::size_t i = 0; i < private_bool_matrix1.size(); i++) {
            private_bool_matrix1(i) = (16 - i) ^ (16 - i + 1);
        }
        std::vector<RegisterAddress> share_reg_addr1;
        share_reg_addr1.push_back(vm.set_data(std::make_shared<PrivateMatrixBool>(private_bool_matrix1)));
        share_reg_addr1.push_back(vm.new_data<BoolMatrix>());
        vm.exec_code(share_inst, share_reg_addr1);

        std::vector<RegisterAddress> inst_reg_addr;
        inst_reg_addr.push_back(share_reg_addr0[1]);
        inst_reg_addr.push_back(share_reg_addr1[1]);
        inst_reg_addr.push_back(vm.new_data<BoolMatrix>());

        Instruction bool_inst({test_type, "bm", "bm", "bm"});
        vm.exec_code(bool_inst, inst_reg_addr);

        Instruction reveal_inst({"reveal", "bm", "pbm"});
        std::vector<RegisterAddress> tmp;
        tmp.push_back(inst_reg_addr[2]);
        tmp.push_back(vm.set_data(std::make_shared<PrivateMatrixBool>(0)));
        vm.exec_code(reveal_inst, tmp);

        if (party_id == 0) {
            Matrix<std::int64_t> aim_ret(4, 3);
            for (size_t i = 0; i < aim_ret.size(); ++i) {
                if (test_type == "and") {
                    aim_ret(i) = private_bool_matrix0(i) & private_bool_matrix1(i);
                } else if (test_type == "or") {
                    aim_ret(i) = private_bool_matrix0(i) | private_bool_matrix1(i);
                } else if (test_type == "xor") {
                    aim_ret(i) = private_bool_matrix0(i) ^ private_bool_matrix1(i);
                }
            }

            Matrix<std::int64_t> real_ret = vm.get_data<PrivateMatrixBool>(tmp[1])->matrix();
            EXPECT_EQ(true, is_equal_plain_matrix(aim_ret, real_ret));
        }

        delete_data(vm, inst_reg_addr);
        delete_data(vm, share_reg_addr0);
        delete_data(vm, share_reg_addr1);
        delete_data(vm, tmp);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_airth_float(std::size_t party_id, std::string test_type) {
        DuetVM vm = init_vm(party_id);

        std::vector<RegisterAddress> private_reg_addr;
        prepare_private_data(party_id, vm, private_reg_addr, 5, 1);
        // now [0] is 0's private, [1] is 1's private

        std::vector<RegisterAddress> share_reg_addr;
        share_private_data(vm, private_reg_addr, share_reg_addr);
        // now [0][1] are both share

        std::vector<RegisterAddress> inst_reg_addr;
        inst_reg_addr.push_back(share_reg_addr[0]);
        auto tmp_double = std::make_shared<PublicDouble>(3);
        inst_reg_addr.push_back(vm.set_data(tmp_double));
        inst_reg_addr.push_back(vm.new_data<ArithMatrix>());
        // now [2] is the ret
        Instruction airth_inst({test_type, "am", "cd", "am"});
        vm.exec_code(airth_inst, inst_reg_addr);

        Instruction reveal_inst({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> tmp;
        private_reg_addr.push_back(vm.set_data(std::make_shared<PrivateMatrix<double>>(0)));
        tmp.push_back(inst_reg_addr[2]);
        tmp.push_back(private_reg_addr[2]);
        vm.exec_code(reveal_inst, tmp);

        if (party_id == 0) {
            Matrix<double> input_0 = vm.get_data<PrivateMatrix<double>>(private_reg_addr[0])->matrix();
            PublicDouble input_1 = *vm.get_data<PublicDouble>(inst_reg_addr[1]);
            Matrix<double> aim_ret;
            if (test_type == "add") {
                aim_ret = input_0.array() + input_1;
            } else if (test_type == "sub") {
                aim_ret = input_0.array() - input_1;
            } else if (test_type == "mul") {
                aim_ret = input_0.array() * input_1;
            } else if (test_type == "div") {
                aim_ret = input_0.array() / input_1;
            }

            Matrix<double> real_ret = vm.get_data<PrivateMatrix<double>>(private_reg_addr[2])->matrix();

            EXPECT_EQ(true, is_equal_plain_matrix(aim_ret, real_ret, 0.001));
        }

        delete_data(vm, private_reg_addr);
        delete_data(vm, share_reg_addr);
        delete_data(vm, inst_reg_addr);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_airth_float_matrix(std::size_t party_id, std::string test_type) {
        DuetVM vm = init_vm(party_id);

        std::vector<RegisterAddress> private_reg_addr;
        prepare_private_data(party_id, vm, private_reg_addr, 4, 3);
        // now [0] is 0's private, [1] is 1's private

        std::vector<RegisterAddress> share_reg_addr;
        share_private_data(vm, private_reg_addr, share_reg_addr);
        // now [0][1] are both share

        std::vector<RegisterAddress> inst_reg_addr;
        inst_reg_addr.push_back(share_reg_addr[0]);
        PublicMatrix<double> double_matrix(4, 3);
        double_matrix.matrix() << -2.13039, -0.875821, -0.436358, -1.1604, -0.710526, 1.44121, -0.233295, -1.7699,
                0.390775, -0.23, 0.61128, -0.464325;
        auto tmp_double = std::make_shared<PublicMatrix<double>>(double_matrix);
        inst_reg_addr.push_back(vm.set_data(tmp_double));
        inst_reg_addr.push_back(vm.new_data<ArithMatrix>());
        // now [2] is the ret

        Instruction airth_inst({test_type, "am", "cdm", "am"});
        vm.exec_code(airth_inst, inst_reg_addr);

        Instruction reveal_inst({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> tmp;
        private_reg_addr.push_back(vm.set_data(std::make_shared<PrivateMatrix<double>>(0)));
        tmp.push_back(inst_reg_addr[2]);
        tmp.push_back(private_reg_addr[2]);
        vm.exec_code(reveal_inst, tmp);

        if (party_id == 0) {
            Matrix<double> input_0 = vm.get_data<PrivateMatrix<double>>(private_reg_addr[0])->matrix();
            Matrix<double> input_1 = vm.get_data<PublicMatrix<double>>(inst_reg_addr[1])->matrix();
            Matrix<double> aim_ret;
            if (test_type == "add") {
                aim_ret = input_0 + input_1;
            } else if (test_type == "sub") {
                aim_ret = input_0 - input_1;
            } else if (test_type == "mul") {
                aim_ret = input_0.cwiseProduct(input_1);
            } else if (test_type == "div") {
                aim_ret.resize(input_0.rows(), input_0.cols());
                for (std::size_t i = 0; i < input_0.size(); i++) {
                    aim_ret(i) = input_0(i) / input_1(i);
                }
            }

            Matrix<double> real_ret = vm.get_data<PrivateMatrix<double>>(private_reg_addr[2])->matrix();
            EXPECT_EQ(true, is_equal_plain_matrix(aim_ret, real_ret, 0.001));
        }

        delete_data(vm, private_reg_addr);
        delete_data(vm, share_reg_addr);
        delete_data(vm, inst_reg_addr);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_mat_airth_float_matrix(std::size_t party_id, std::string test_type) {
        DuetVM vm = init_vm(party_id);

        std::vector<RegisterAddress> private_reg_addr;
        prepare_private_data(party_id, vm, private_reg_addr, 4, 4);

        std::vector<RegisterAddress> share_reg_addr;
        share_private_data(vm, private_reg_addr, share_reg_addr);

        std::vector<RegisterAddress> inst_reg_addr0;
        std::vector<RegisterAddress> inst_reg_addr1;
        PublicMatrix<double> double_matrix(4, 4);
        double_matrix.matrix() << -2.13039, -0.875821, -0.436358, -1.1604, -0.710526, 1.44121, -0.233295, -1.7699,
                0.390775, -0.23, 0.61128, -0.464325, 0.609836, 0.140385, 0.919206, 0.918789;
        auto tmp_double = std::make_shared<PublicMatrix<double>>(double_matrix);

        inst_reg_addr0.push_back(vm.set_data(tmp_double));
        inst_reg_addr0.push_back(share_reg_addr[0]);
        inst_reg_addr0.push_back(vm.new_data<ArithMatrix>());

        Instruction airth_inst0({test_type, "cdm", "am", "am"});
        vm.exec_code(airth_inst0, inst_reg_addr0);

        Instruction reveal_inst0({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> tmp0;
        private_reg_addr.push_back(vm.set_data(std::make_shared<PrivateMatrix<double>>(0)));
        tmp0.push_back(inst_reg_addr0[2]);
        tmp0.push_back(private_reg_addr[2]);
        vm.exec_code(reveal_inst0, tmp0);

        inst_reg_addr1.push_back(share_reg_addr[0]);
        inst_reg_addr1.push_back(vm.set_data(tmp_double));
        inst_reg_addr1.push_back(vm.new_data<ArithMatrix>());

        Instruction airth_inst1({test_type, "am", "cdm", "am"});
        vm.exec_code(airth_inst1, inst_reg_addr1);

        Instruction reveal_inst1({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> tmp1;
        private_reg_addr.push_back(vm.set_data(std::make_shared<PrivateMatrix<double>>(0)));
        tmp1.push_back(inst_reg_addr1[2]);
        tmp1.push_back(private_reg_addr[3]);
        vm.exec_code(reveal_inst1, tmp1);

        if (party_id == 0) {
            Matrix<double> input_0 = vm.get_data<PrivateMatrix<double>>(private_reg_addr[0])->matrix();
            Matrix<double> input_1 = vm.get_data<PublicMatrix<double>>(inst_reg_addr0[0])->matrix();
            Matrix<double> aim_ret[2];
            Matrix<double> real_ret[2];
            aim_ret[0] = input_1 * input_0;
            real_ret[0] = vm.get_data<PrivateMatrix<double>>(private_reg_addr[2])->matrix();

            aim_ret[1] = input_0 * input_1;
            real_ret[1] = vm.get_data<PrivateMatrix<double>>(private_reg_addr[3])->matrix();

            EXPECT_EQ(true, is_equal_plain_matrix(aim_ret[0], real_ret[0], 0.1));
            EXPECT_EQ(true, is_equal_plain_matrix(aim_ret[1], real_ret[1], 0.1));
        }

        delete_data(vm, private_reg_addr);
        delete_data(vm, share_reg_addr);
        delete_data(vm, inst_reg_addr0);
        delete_data(vm, inst_reg_addr1);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_bool_float(std::size_t party_id, std::string test_type) {
        DuetVM vm = init_vm(party_id);

        std::vector<RegisterAddress> private_reg_addr;
        prepare_private_data(party_id, vm, private_reg_addr, 5, 1);
        // now [0] is 0's private, [1] is 1's private

        std::vector<RegisterAddress> share_reg_addr;
        share_private_data(vm, private_reg_addr, share_reg_addr);
        // now [0][1] are both share

        std::vector<RegisterAddress> inst_reg_addr;
        inst_reg_addr.push_back(share_reg_addr[0]);
        auto tmp_double = std::make_shared<PublicDouble>(3);
        inst_reg_addr.push_back(vm.set_data(tmp_double));
        inst_reg_addr.push_back(vm.new_data<BoolMatrix>());
        // now [2] is the ret
        Instruction bool_inst({test_type, "am", "cd", "bm"});
        vm.exec_code(bool_inst, inst_reg_addr);

        Instruction reveal_inst({"reveal", "bm", "pbm"});
        std::vector<RegisterAddress> tmp;
        private_reg_addr.push_back(vm.set_data(std::make_shared<PrivateMatrixBool>(0)));
        tmp.push_back(inst_reg_addr[2]);
        tmp.push_back(private_reg_addr[2]);
        vm.exec_code(reveal_inst, tmp);

        if (party_id == 0) {
            Matrix<double> input_0 = vm.get_data<PrivateMatrix<double>>(private_reg_addr[0])->matrix();
            PublicDouble input_1 = *vm.get_data<PublicDouble>(inst_reg_addr[1]);
            Matrix<std::int64_t> aim_ret(5, 1);
            for (size_t i = 0; i < 5; ++i) {
                if (test_type == "gt") {
                    aim_ret(i) = input_0(i) > input_1 ? 1 : 0;
                } else if (test_type == "lt") {
                    aim_ret(i) = input_0(i) < input_1 ? 1 : 0;
                }
            }
            Matrix<std::int64_t> real_ret = vm.get_data<PrivateMatrixBool>(private_reg_addr[2])->matrix();

            EXPECT_EQ(true, is_equal_plain_matrix(aim_ret, real_ret));
        }

        delete_data(vm, private_reg_addr);
        delete_data(vm, share_reg_addr);
        delete_data(vm, inst_reg_addr);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_bool_float_matrix(std::size_t party_id, std::string test_type) {
        DuetVM vm = init_vm(party_id);

        std::vector<RegisterAddress> private_reg_addr;
        prepare_private_data(party_id, vm, private_reg_addr, 4, 3);
        // now [0] is 0's private, [1] is 1's private

        std::vector<RegisterAddress> share_reg_addr;
        share_private_data(vm, private_reg_addr, share_reg_addr);
        // now [0][1] are both share

        std::vector<RegisterAddress> inst_reg_addr;
        inst_reg_addr.push_back(share_reg_addr[0]);
        PublicMatrix<double> double_matrix(4, 3);
        double_matrix.matrix() << -2.13039, -0.875821, -0.436358, -1.1604, -0.710526, 1.44121, -0.233295, -1.7699,
                0.390775, -0.23, 0.61128, -0.464325;
        auto tmp_double = std::make_shared<PublicMatrix<double>>(double_matrix);
        inst_reg_addr.push_back(vm.set_data(tmp_double));
        inst_reg_addr.push_back(vm.new_data<BoolMatrix>());
        // now [2] is the ret
        Instruction bool_inst({test_type, "am", "cdm", "bm"});
        vm.exec_code(bool_inst, inst_reg_addr);

        Instruction reveal_inst({"reveal", "bm", "pbm"});
        std::vector<RegisterAddress> tmp;
        private_reg_addr.push_back(vm.set_data(std::make_shared<PrivateMatrixBool>(0)));
        tmp.push_back(inst_reg_addr[2]);
        tmp.push_back(private_reg_addr[2]);
        vm.exec_code(reveal_inst, tmp);

        if (party_id == 0) {
            Matrix<double> input_0 = vm.get_data<PrivateMatrix<double>>(private_reg_addr[0])->matrix();
            Matrix<double> input_1 = vm.get_data<PublicMatrix<double>>(inst_reg_addr[1])->matrix();
            Matrix<std::int64_t> aim_ret(4, 3);
            for (size_t i = 0; i < aim_ret.size(); ++i) {
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
        delete_data(vm, inst_reg_addr);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_quick_sort(std::size_t party_id) {
        DuetVM vm = init_vm(party_id);
        std::size_t matrix_size = 20;
        std::shared_ptr<PrivateMatrix<double>> a = std::make_shared<PrivateMatrix<double>>(matrix_size, 1, 0);
        std::shared_ptr<PrivateMatrix<double>> real_ret = std::make_shared<PrivateMatrix<double>>(matrix_size, 1, 0);
        unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
        std::default_random_engine generator(seed);
        PrivateMatrix<double> aim_ret(matrix_size, 1, 0);

        std::vector<double> arr(matrix_size, 0);
        for (std::size_t i = 0; i < arr.size(); ++i) {
            arr[i] = static_cast<double>(i);
        }

        std::shuffle(arr.begin(), arr.end(), generator);

        if (party_id == 0) {
            for (std::size_t i = 0; i < arr.size(); ++i) {
                (*a)(i) = arr[i];
                aim_ret(i) = static_cast<double>(i);
            }
        }
        Instruction share_inst({"share", "pdm", "am"});
        std::vector<RegisterAddress> share_reg;
        share_reg.push_back(vm.set_data(a));
        share_reg.push_back(vm.new_data<ArithMatrix>());
        vm.exec_code(share_inst, share_reg);
        Instruction sort_inst({"quick_sort", "am"});
        std::vector<RegisterAddress> sort_reg{share_reg[1]};
        vm.exec_code(sort_inst, sort_reg);

        Instruction reveal_inst({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> reveal_reg{sort_reg[0], vm.set_data(real_ret)};
        vm.exec_code(reveal_inst, reveal_reg);
        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(aim_ret, *real_ret, 0.001));
        }
        delete_data(vm, sort_reg);
        delete_data(vm, share_reg);
        delete_data(vm, reveal_reg);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_quick_sort_by_column(std::size_t party_id) {
        DuetVM vm = init_vm(party_id);
        std::size_t matrix_size = 20;
        std::shared_ptr<PrivateMatrix<double>> a = std::make_shared<PrivateMatrix<double>>(matrix_size, 1, 0);
        std::shared_ptr<PrivateMatrix<double>> real_ret = std::make_shared<PrivateMatrix<double>>(matrix_size, 1, 0);
        std::shared_ptr<PublicIndex> index = std::make_shared<PublicIndex>(0);
        unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
        std::default_random_engine generator(seed);
        PrivateMatrix<double> aim_ret(matrix_size, 1, 0);

        std::vector<double> arr(matrix_size, 0);
        for (std::size_t i = 0; i < arr.size(); ++i) {
            arr[i] = static_cast<double>(i);
        }

        std::shuffle(arr.begin(), arr.end(), generator);

        if (party_id == 0) {
            for (std::size_t i = 0; i < arr.size(); ++i) {
                (*a)(i) = arr[i];
                aim_ret(i) = static_cast<double>(i);
            }
        }
        Instruction share_inst({"share", "pdm", "am"});
        std::vector<RegisterAddress> share_reg;
        share_reg.push_back(vm.set_data(a));
        share_reg.push_back(vm.new_data<ArithMatrix>());
        vm.exec_code(share_inst, share_reg);
        Instruction sort_inst({"quick_sort", "am", "ci", "am"});
        std::vector<RegisterAddress> sort_reg{share_reg[1], vm.set_data(index), vm.new_data<ArithMatrix>()};
        vm.exec_code(sort_inst, sort_reg);

        Instruction reveal_inst({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> reveal_reg{sort_reg[2], vm.set_data(real_ret)};
        vm.exec_code(reveal_inst, reveal_reg);
        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(aim_ret, *real_ret, 0.001));
        }
        delete_data(vm, sort_reg);
        delete_data(vm, share_reg);
        delete_data(vm, reveal_reg);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_sigmoid(std::size_t party_id) {
        DuetVM vm = init_vm(party_id);
        std::size_t matrix_size = 20;
        std::shared_ptr<PrivateMatrix<double>> a = std::make_shared<PrivateMatrix<double>>(matrix_size, 1, 0);
        std::shared_ptr<PrivateMatrix<double>> real_ret = std::make_shared<PrivateMatrix<double>>(matrix_size, 1, 0);
        unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
        std::default_random_engine generator(seed);
        PrivateMatrix<double> aim_ret(matrix_size, 1, 0);

        std::vector<double> arr(matrix_size, 0);
        for (std::size_t i = 0; i < arr.size(); ++i) {
            (*a)(i) = -1.0 + static_cast<double>(i) * 0.1;
            aim_ret(i) = 1.0 / (1.0 + std::exp(-(*a)(i)));
        }

        Instruction share_inst({"share", "pdm", "am"});
        std::vector<RegisterAddress> share_reg;
        share_reg.push_back(vm.set_data(a));
        share_reg.push_back(vm.new_data<ArithMatrix>());
        vm.exec_code(share_inst, share_reg);
        Instruction sigmoid_inst({"sigmoid", "am", "am"});
        std::vector<RegisterAddress> sigmoid_reg{share_reg[1], vm.new_data<ArithMatrix>()};
        vm.exec_code(sigmoid_inst, sigmoid_reg);

        Instruction reveal_inst({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> reveal_reg{sigmoid_reg[1], vm.set_data(real_ret)};
        vm.exec_code(reveal_inst, reveal_reg);
        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(aim_ret, *real_ret, 0.001));
        }
        delete_data(vm, sigmoid_reg);
        delete_data(vm, share_reg);
        delete_data(vm, reveal_reg);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_transpose(std::size_t party_id) {
        DuetVM vm = init_vm(party_id);
        std::size_t matrix_size = 4;
        std::shared_ptr<PrivateMatrix<double>> a = std::make_shared<PrivateMatrix<double>>(matrix_size, 1, 0);
        std::shared_ptr<PrivateMatrix<double>> real_ret = std::make_shared<PrivateMatrix<double>>();
        unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
        std::default_random_engine generator(seed);
        PrivateMatrix<double> aim_ret(1, matrix_size, 0);

        if (party_id == 0) {
            for (std::size_t i = 0; i < matrix_size; ++i) {
                (*a)(i) = static_cast<double>(i);
                aim_ret(i) = static_cast<double>(i);
            }
        }
        Instruction share_inst({"share", "pdm", "am"});
        std::vector<RegisterAddress> share_reg;
        share_reg.push_back(vm.set_data(a));
        share_reg.push_back(vm.new_data<ArithMatrix>());
        vm.exec_code(share_inst, share_reg);
        Instruction transpose_inst({"transpose", "am", "am"});
        std::vector<RegisterAddress> transpose_reg{share_reg[1], vm.new_data<ArithMatrix>()};
        vm.exec_code(transpose_inst, transpose_reg);

        Instruction reveal_inst({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> reveal_reg{transpose_reg[1], vm.set_data(real_ret)};
        vm.exec_code(reveal_inst, reveal_reg);

        if (party_id == 0) {
            ASSERT_TRUE(aim_ret.cols() == real_ret->cols());
            ASSERT_TRUE(aim_ret.rows() == real_ret->rows());
            ASSERT_TRUE(is_equal_private_matrix(aim_ret, *real_ret, 0.001));
        }
        delete_data(vm, share_reg);
        delete_data(vm, transpose_reg);
        delete_data(vm, reveal_reg);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_reshape(std::size_t party_id) {
        DuetVM vm = init_vm(party_id);
        std::size_t matrix_size = 20;
        std::size_t aim_rows = 4;
        std::size_t aim_cols = 5;
        std::shared_ptr<PrivateMatrix<double>> a = std::make_shared<PrivateMatrix<double>>(matrix_size, 1, 0);
        std::shared_ptr<PrivateMatrix<double>> real_ret = std::make_shared<PrivateMatrix<double>>();
        unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
        std::default_random_engine generator(seed);
        PrivateMatrix<double> aim_ret(aim_rows, aim_cols, 0);

        if (party_id == 0) {
            for (std::size_t i = 0; i < matrix_size; ++i) {
                (*a)(i) = static_cast<double>(i);
                aim_ret(i) = static_cast<double>(i);
            }
        }
        Instruction share_inst({"share", "pdm", "am"});
        std::vector<RegisterAddress> share_reg;
        share_reg.push_back(vm.set_data(a));
        share_reg.push_back(vm.new_data<ArithMatrix>());
        vm.exec_code(share_inst, share_reg);
        Instruction reshape_inst({"reshape", "am", "ci", "ci", "am"});
        std::vector<RegisterAddress> reshape_reg{share_reg[1], vm.set_data(std::make_shared<std::size_t>(aim_rows)),
                vm.set_data(std::make_shared<std::size_t>(aim_cols)), vm.new_data<ArithMatrix>()};
        vm.exec_code(reshape_inst, reshape_reg);

        Instruction reveal_inst({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> reveal_reg{reshape_reg[3], vm.set_data(real_ret)};
        vm.exec_code(reveal_inst, reveal_reg);

        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(aim_ret, *real_ret, 0.001));
        }
        delete_data(vm, share_reg);
        delete_data(vm, reshape_reg);
        delete_data(vm, reveal_reg);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_set_item(std::size_t party_id) {
        DuetVM vm = init_vm(party_id);
        std::size_t rows = 4;
        std::size_t cols = 5;
        std::shared_ptr<PrivateMatrix<double>> a = std::make_shared<PrivateMatrix<double>>(rows, cols, 0);
        std::shared_ptr<PrivateMatrix<double>> real_ret = std::make_shared<PrivateMatrix<double>>();
        unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
        std::default_random_engine generator(seed);
        PrivateMatrix<double> aim_ret(rows, cols, 0);
        aim_ret.matrix() << 1, 2, 3, 4, 5, 2, 1, 0, 0, 0, 3, 0, 0, 0, 0, 4, 0, 0, 0, 0;

        if (party_id == 0) {
            for (std::size_t i = 0; i < rows * cols; ++i) {
                (*a)(i) = static_cast<double>(0);
            }
        }
        Instruction share_inst({"share", "pdm", "am"});
        std::vector<RegisterAddress> share_reg;
        share_reg.push_back(vm.set_data(a));
        share_reg.push_back(vm.new_data<ArithMatrix>());
        vm.exec_code(share_inst, share_reg);

        a->resize(1, 1);
        if (party_id == 0) {
            for (std::size_t i = 0; i < a->size(); ++i) {
                (*a)(i) = static_cast<double>(i + 1);
            }
        }
        std::vector<RegisterAddress> set_reg;
        set_reg.push_back(vm.set_data(a));
        set_reg.push_back(vm.new_data<ArithMatrix>());
        vm.exec_code(share_inst, set_reg);

        Instruction set_item_inst({"set_item", "am", "am", "ci", "ci", "ci", "ci"});
        std::vector<RegisterAddress> set_item_reg{set_reg[1], share_reg[1],
                vm.set_data(std::make_shared<PublicIndex>(1)), vm.set_data(std::make_shared<PublicIndex>(1)),
                vm.set_data(std::make_shared<PublicIndex>(1)), vm.set_data(std::make_shared<PublicIndex>(1))};
        vm.exec_code(set_item_inst, set_item_reg);

        Instruction reveal_inst({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> reveal_reg{set_item_reg[1], vm.set_data(real_ret)};
        vm.exec_code(reveal_inst, reveal_reg);

        a->resize(1, 5);
        if (party_id == 0) {
            for (std::size_t i = 0; i < a->size(); ++i) {
                (*a)(i) = static_cast<double>(i + 1);
            }
        };
        std::vector<RegisterAddress> set_reg_row;
        set_reg_row.push_back(vm.set_data(a));
        set_reg_row.push_back(vm.new_data<ArithMatrix>());
        vm.exec_code(share_inst, set_reg_row);

        std::vector<RegisterAddress> set_item_reg_row{set_reg_row[1], share_reg[1],
                vm.set_data(std::make_shared<PublicIndex>(0)), vm.set_data(std::make_shared<PublicIndex>(0)),
                vm.set_data(std::make_shared<PublicIndex>(1)), vm.set_data(std::make_shared<PublicIndex>(5))};
        vm.exec_code(set_item_inst, set_item_reg_row);

        Instruction reveal_inst_row({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> reveal_reg_row{set_item_reg[1], vm.set_data(real_ret)};
        vm.exec_code(reveal_inst_row, reveal_reg_row);

        a->resize(4, 1);
        if (party_id == 0) {
            for (std::size_t i = 0; i < a->size(); ++i) {
                (*a)(i) = static_cast<double>(i + 1);
            }
        };
        std::vector<RegisterAddress> set_reg_col;
        set_reg_col.push_back(vm.set_data(a));
        set_reg_col.push_back(vm.new_data<ArithMatrix>());
        vm.exec_code(share_inst, set_reg_col);
        std::vector<RegisterAddress> set_item_reg_col{set_reg_col[1], share_reg[1],
                vm.set_data(std::make_shared<PublicIndex>(0)), vm.set_data(std::make_shared<PublicIndex>(0)),
                vm.set_data(std::make_shared<PublicIndex>(4)), vm.set_data(std::make_shared<PublicIndex>(1))};
        vm.exec_code(set_item_inst, set_item_reg_col);

        Instruction reveal_inst_col({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> reveal_reg_col{set_item_reg[1], vm.set_data(real_ret)};
        vm.exec_code(reveal_inst_col, reveal_reg_col);

        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(aim_ret, *real_ret, 0.001));
        }
        delete_data(vm, share_reg);
        delete_data(vm, set_reg);
        delete_data(vm, set_reg_row);
        delete_data(vm, set_reg_col);
        delete_data(vm, set_item_reg);
        delete_data(vm, set_item_reg_row);
        delete_data(vm, set_item_reg_col);
        delete_data(vm, reveal_reg);
        delete_data(vm, reveal_reg_row);
        delete_data(vm, reveal_reg_col);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_resize(std::size_t party_id) {
        DuetVM vm = init_vm(party_id);
        std::size_t matrix_size = 20;
        std::size_t aim_rows = 4;
        std::size_t aim_cols = 8;
        std::shared_ptr<PrivateMatrix<double>> a = std::make_shared<PrivateMatrix<double>>(matrix_size, 1, 0);
        std::shared_ptr<PrivateMatrix<double>> real_ret = std::make_shared<PrivateMatrix<double>>();
        unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
        std::default_random_engine generator(seed);
        PrivateMatrix<double> aim_ret(aim_rows, aim_cols, 0);

        if (party_id == 0) {
            for (std::size_t i = 0; i < aim_ret.size(); ++i) {
                (*a)(i % matrix_size) = static_cast<double>(i % matrix_size);
                aim_ret(i) = static_cast<double>(i % matrix_size);
            }
        }
        Instruction share_inst({"share", "pdm", "am"});
        std::vector<RegisterAddress> share_reg;
        share_reg.push_back(vm.set_data(a));
        share_reg.push_back(vm.new_data<ArithMatrix>());
        vm.exec_code(share_inst, share_reg);
        Instruction resize_inst({"resize", "am", "ci", "ci", "am"});
        std::vector<RegisterAddress> resize_reg{share_reg[1], vm.set_data(std::make_shared<std::size_t>(aim_rows)),
                vm.set_data(std::make_shared<std::size_t>(aim_cols)), vm.new_data<ArithMatrix>()};
        vm.exec_code(resize_inst, resize_reg);

        Instruction reveal_inst({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> reveal_reg{resize_reg[3], vm.set_data(real_ret)};
        vm.exec_code(reveal_inst, reveal_reg);

        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(aim_ret, *real_ret, 0.001));
        }
        delete_data(vm, share_reg);
        delete_data(vm, resize_reg);
        delete_data(vm, reveal_reg);
        EXPECT_EQ(true, vm.is_registr_empty());
    }

    void test_groupby(std::size_t party_id, std::string test_type) {
        DuetVM vm = init_vm(party_id);

        std::size_t row = 6;
        std::size_t col = 2;
        std::shared_ptr<PrivateMatrix<double>> plain0 = std::make_shared<PrivateMatrix<double>>(row, col, 0);
        std::shared_ptr<PrivateMatrix<double>> one_hot_matrix = std::make_shared<PrivateMatrix<double>>(row, col, 0);

        plain0->matrix() << -15.812, -14.7387, 7.8120, -9.7387, -2.8120, 6.7387, 1.8120, 1.7387, 5.8120, 0.7387,
                -7.8120, -9.7387;
        one_hot_matrix->matrix() << 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1;

        Instruction share_inst({"share", "pdm", "am"});
        std::vector<RegisterAddress> share_data_reg;
        share_data_reg.push_back(vm.set_data(plain0));
        share_data_reg.push_back(vm.new_data<ArithMatrix>());
        vm.exec_code(share_inst, share_data_reg);

        std::vector<RegisterAddress> share_index_reg;
        share_index_reg.push_back(vm.set_data(one_hot_matrix));
        share_index_reg.push_back(vm.new_data<ArithMatrix>());
        vm.exec_code(share_inst, share_index_reg);

        std::vector<RegisterAddress> inst_reg_addr;
        inst_reg_addr.push_back(share_data_reg[1]);
        inst_reg_addr.push_back(share_index_reg[1]);
        inst_reg_addr.push_back(vm.new_data<ArithMatrix>());
        // now [2] is the ret
        Instruction groupby_inst({test_type, "am", "am", "am"});
        vm.exec_code(groupby_inst, inst_reg_addr);

        Instruction reveal_inst({"reveal", "am", "pdm"});
        std::vector<RegisterAddress> reveal_reg;
        reveal_reg.push_back(inst_reg_addr[2]);
        reveal_reg.push_back(vm.set_data(std::make_shared<PrivateMatrix<double>>(0)));
        vm.exec_code(reveal_inst, reveal_reg);

        if (party_id == 0) {
            Matrix<double> aim_ret(2, 2);
            if (test_type == "groupby_sum") {
                aim_ret.matrix() << -14, -8.188, -32.4774, -31.7387;
            } else if (test_type == "groupby_count") {
                aim_ret.matrix() << 4, 5, 4, 5;
            } else if (test_type == "groupby_max") {
                aim_ret.matrix() << 7.8120, 7.8120, 1.73869, 1.7387;
            } else if (test_type == "groupby_min") {
                aim_ret.matrix() << -15.812, -15.812, -14.7387, -14.7387;
            }
            Matrix<double> real_ret = vm.get_data<PrivateMatrix<double>>(reveal_reg[1])->matrix();
            EXPECT_EQ(true, is_equal_plain_matrix(aim_ret, real_ret, 0.001));
        }

        delete_data(vm, share_data_reg);
        delete_data(vm, share_index_reg);
        delete_data(vm, inst_reg_addr);
        delete_data(vm, reveal_reg);
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

TEST_F(DuetVMTest, test_mat_airth) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_mat_airth(i, "mat_mul"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_add_float) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_airth_float(i, "add"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_sub_float) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_airth_float(i, "sub"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_mul_float) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_airth_float(i, "mul"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_div_float) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_airth_float(i, "div"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_add_float_matrix) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_airth_float_matrix(i, "add"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_sub_float_matrix) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_airth_float_matrix(i, "sub"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_mul_float_matrix) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_airth_float_matrix(i, "mul"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_mat_mul_float_matrix) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_mat_airth_float_matrix(i, "mat_mul"); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_div_float_matrix) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_airth_float_matrix(i, "div"); });
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

TEST_F(DuetVMTest, test_bool_float) {
    std::vector<std::string> test_types = {"gt", "lt"};
    for (std::size_t i = 0; i < test_types.size(); i++) {
        std::vector<std::thread> threads;
        std::cout << "test for: " << test_types[i] << std::endl;
        for (std::size_t j = 0; j < 2; ++j) {
            threads.emplace_back([&, j]() { test_bool_float(j, test_types[i]); });
        }
        for (auto& thread : threads) {
            thread.join();
        }
    }
}

TEST_F(DuetVMTest, test_bool_float_matrix) {
    std::vector<std::string> test_types = {"gt", "lt", "ge", "eq"};
    for (std::size_t i = 0; i < test_types.size(); i++) {
        std::vector<std::thread> threads;
        std::cout << "test for: " << test_types[i] << std::endl;
        for (std::size_t j = 0; j < 2; ++j) {
            threads.emplace_back([&, j]() { test_bool_float_matrix(j, test_types[i]); });
        }
        for (auto& thread : threads) {
            thread.join();
        }
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

TEST_F(DuetVMTest, test_arith_bool) {
    std::vector<std::string> test_types = {"and", "or", "xor"};
    for (std::size_t i = 0; i < test_types.size(); i++) {
        std::vector<std::thread> threads;
        std::cout << "test for: " << test_types[i] << std::endl;
        for (std::size_t j = 0; j < 2; ++j) {
            threads.emplace_back([&, j]() { test_arith_bool(j, test_types[i]); });
        }
        for (auto& thread : threads) {
            thread.join();
        }
    }
}

TEST_F(DuetVMTest, test_quick_sort) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_quick_sort(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_quick_sort_by_column) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_quick_sort_by_column(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_sigmoid) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_sigmoid(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_transpose) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_transpose(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_reshape) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_reshape(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_set_item) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_set_item(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_resize) {
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { test_resize(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(DuetVMTest, test_groupby) {
    std::vector<std::string> test_types = {"groupby_sum", "groupby_count", "groupby_max", "groupby_min"};
    for (std::size_t i = 0; i < test_types.size(); i++) {
        std::vector<std::thread> threads;
        std::cout << "test for: " << test_types[i] << std::endl;
        for (std::size_t j = 0; j < 2; ++j) {
            threads.emplace_back([&, j]() { test_groupby(j, test_types[i]); });
        }
        for (auto& thread : threads) {
            thread.join();
        }
    }
}

}  // namespace duet
}  // namespace petace
