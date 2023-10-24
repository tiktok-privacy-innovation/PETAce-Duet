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

#include "duet/vm/duet_vm.h"

namespace petace {
namespace duet {

void DuetVM::registr_default_op() {
    Instruction share_double_p_a({"share", "pdm", "am"});
    inst_op_map_[share_double_p_a] = std::make_shared<ShareDoublePA>();

    Instruction reveal_double_a_p({"reveal", "am", "pdm"});
    inst_op_map_[reveal_double_a_p] = std::make_shared<RevealDoubleAP>();

    Instruction reveal_double_a_c({"reveal", "am", "cdm"});
    inst_op_map_[reveal_double_a_c] = std::make_shared<RevealDoubleAC>();

    Instruction reveal_bool_b_p({"reveal", "bm", "pbm"});
    inst_op_map_[reveal_bool_b_p] = std::make_shared<RevealBoolBP>();

    Instruction add_a_a_a({"add", "am", "am", "am"});
    inst_op_map_[add_a_a_a] = std::make_shared<AddAAA>();

    Instruction sub_a_a_a({"sub", "am", "am", "am"});
    inst_op_map_[sub_a_a_a] = std::make_shared<SubAAA>();

    Instruction sub_a_cd_a({"sub", "am", "cd", "am"});
    inst_op_map_[sub_a_cd_a] = std::make_shared<SubACdA>();

    Instruction mul_b_b_b({"mul", "bm", "bm", "bm"});
    inst_op_map_[mul_b_b_b] = std::make_shared<MulBBB>();

    Instruction or_cbm_b_b({"or", "cbm", "bm", "bm"});
    inst_op_map_[or_cbm_b_b] = std::make_shared<OrCBB>();

    Instruction mul_a_a_a({"mul", "am", "am", "am"});
    inst_op_map_[mul_a_a_a] = std::make_shared<MulAAA>();

    Instruction mul_a_cd_a({"mul", "am", "cd", "am"});
    inst_op_map_[mul_a_cd_a] = std::make_shared<MulACdA>();

    Instruction div_a_cdm_a({"div", "am", "cdm", "am"});
    inst_op_map_[div_a_cdm_a] = std::make_shared<DivACdA>();

    Instruction div_a_a_a({"div", "am", "am", "am"});
    inst_op_map_[div_a_a_a] = std::make_shared<DivAAA>();

    Instruction greater_a_a_b({"gt", "am", "am", "bm"});
    inst_op_map_[greater_a_a_b] = std::make_shared<GreaterAAB>();

    Instruction less_a_a_b({"lt", "am", "am", "bm"});
    inst_op_map_[less_a_a_b] = std::make_shared<LessAAB>();

    Instruction greater_eq_a_a_b({"ge", "am", "am", "bm"});
    inst_op_map_[greater_eq_a_a_b] = std::make_shared<GreatEqualAAB>();

    Instruction equal_a_a_b({"eq", "am", "am", "bm"});
    inst_op_map_[equal_a_a_b] = std::make_shared<EqualAAB>();

    Instruction sum_a_a({"sum", "am", "am"});
    inst_op_map_[sum_a_a] = std::make_shared<SumAA>();

    Instruction less_than_zero_a_b({"less_than_zero", "am", "bm"});
    inst_op_map_[less_than_zero_a_b] = std::make_shared<LessThanZeroAB>();

    Instruction multiplexer_b_a_a({"multiplexer", "bm", "am", "am"});
    inst_op_map_[multiplexer_b_a_a] = std::make_shared<MultiplexerBAA>();

    Instruction multiplexer_b_a_a_a({"multiplexer", "bm", "am", "am", "am"});
    inst_op_map_[multiplexer_b_a_a_a] = std::make_shared<MultiplexerBAAA>();

    Instruction shuffle_p_a({"shuffle", "pdm", "am"});
    inst_op_map_[shuffle_p_a] = std::make_shared<ShufflePA>();

    Instruction shuffle_a_a({"shuffle", "am", "am"});
    inst_op_map_[shuffle_a_a] = std::make_shared<ShuffleAA>();

    Instruction row_major_argmax_and_max_a_a_a({"row_major_argmax_and_max", "am", "am", "am"});
    inst_op_map_[row_major_argmax_and_max_a_a_a] = std::make_shared<RowMajorArgmaxAndMaxAAA>();
}

DuetVM::DuetVM(const petace::network::NetParams& net_params, petace::network::NetScheme net_type, std::size_t party_id)
        : party_id_(party_id) {
    net_ = petace::network::NetFactory::get_instance().build(net_type, net_params);
    impl_ = std::make_shared<Duet>(net_, party_id_);
    reg_ = std::make_shared<Register>();
    registr_default_op();
}

void DuetVM::registr_op(const Instruction& inst, const std::shared_ptr<OperatorBase>& op) {
    inst_op_map_[inst] = op;
}

void DuetVM::delete_data(RegisterAddress address) {
    reg_->delete_data(address);
}

void DuetVM::exec_code(const Instruction& inst, const std::vector<RegisterAddress>& addr) {
    inst_op_map_.at(inst)->run(impl_, net_, reg_, addr);
}

bool DuetVM::is_registr_empty() {
    return reg_->empty();
}

std::size_t DuetVM::party_id() {
    return party_id_;
}

}  // namespace duet
}  // namespace petace
