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

#include "network/net_factory.h"

namespace petace {
namespace duet {

void DuetVM::registr_default_op() {
    Instruction share_double_p_a({"share", "pdm", "am"});
    inst_op_map_[share_double_p_a] = std::make_shared<ShareDoublePA>();

    Instruction reveal_double_a_p({"reveal", "am", "pdm"});
    inst_op_map_[reveal_double_a_p] = std::make_shared<RevealDoubleAP>();

    Instruction reveal_double_a_c({"reveal", "am", "cdm"});
    inst_op_map_[reveal_double_a_c] = std::make_shared<RevealDoubleAC>();

    Instruction share_bool_p_a({"share", "pbm", "bm"});
    inst_op_map_[share_bool_p_a] = std::make_shared<ShareBoolPA>();

    Instruction reveal_b_p({"reveal", "bm", "pbm"});
    inst_op_map_[reveal_b_p] = std::make_shared<RevealBoolBP>();

    Instruction add_a_a_a({"add", "am", "am", "am"});
    inst_op_map_[add_a_a_a] = std::make_shared<AddAAA>();

    Instruction add_a_cd_a({"add", "am", "cd", "am"});
    inst_op_map_[add_a_cd_a] = std::make_shared<AddACdA>();

    Instruction add_a_cdm_a({"add", "am", "cdm", "am"});
    inst_op_map_[add_a_cdm_a] = std::make_shared<AddACdmA>();

    Instruction sub_a_a_a({"sub", "am", "am", "am"});
    inst_op_map_[sub_a_a_a] = std::make_shared<SubAAA>();

    Instruction sub_a_cd_a({"sub", "am", "cd", "am"});
    inst_op_map_[sub_a_cd_a] = std::make_shared<SubACdA>();

    Instruction sub_a_cdm_a({"sub", "am", "cdm", "am"});
    inst_op_map_[sub_a_cdm_a] = std::make_shared<SubACdmA>();

    Instruction mul_b_b_b({"mul", "bm", "bm", "bm"});
    inst_op_map_[mul_b_b_b] = std::make_shared<MulBBB>();

    Instruction and_b_b_b({"and", "bm", "bm", "bm"});
    inst_op_map_[and_b_b_b] = std::make_shared<AndBBB>();

    Instruction and_cbm_b_b({"and", "cbm", "bm", "bm"});
    inst_op_map_[and_cbm_b_b] = std::make_shared<AndCBB>();

    Instruction or_b_b_b({"or", "bm", "bm", "bm"});
    inst_op_map_[or_b_b_b] = std::make_shared<OrBBB>();

    Instruction or_cbm_b_b({"or", "cbm", "bm", "bm"});
    inst_op_map_[or_cbm_b_b] = std::make_shared<OrCBB>();

    Instruction xor_b_b_b({"xor", "bm", "bm", "bm"});
    inst_op_map_[xor_b_b_b] = std::make_shared<XorBBB>();

    Instruction xor_cbm_b_b({"xor", "cbm", "bm", "bm"});
    inst_op_map_[xor_cbm_b_b] = std::make_shared<XorCBB>();

    Instruction not_b_b({"not", "bm", "bm"});
    inst_op_map_[not_b_b] = std::make_shared<NotBB>();

    Instruction mul_a_a_a({"mul", "am", "am", "am"});
    inst_op_map_[mul_a_a_a] = std::make_shared<MulAAA>();

    Instruction mat_mul_a_a_a({"mat_mul", "am", "am", "am"});
    inst_op_map_[mat_mul_a_a_a] = std::make_shared<MatMulAAA>();

    Instruction mul_a_cd_a({"mul", "am", "cd", "am"});
    inst_op_map_[mul_a_cd_a] = std::make_shared<MulACdA>();

    Instruction mul_a_cdm_a({"mul", "am", "cdm", "am"});
    inst_op_map_[mul_a_cdm_a] = std::make_shared<MulACdmA>();

    Instruction mat_mul_cdm_a_a({"mat_mul", "cdm", "am", "am"});
    inst_op_map_[mat_mul_cdm_a_a] = std::make_shared<MatMulCdmAA>();

    Instruction mat_mul_a_cdm_a({"mat_mul", "am", "cdm", "am"});
    inst_op_map_[mat_mul_a_cdm_a] = std::make_shared<MatMulACdmA>();

    Instruction div_a_cdm_a({"div", "am", "cdm", "am"});
    inst_op_map_[div_a_cdm_a] = std::make_shared<DivACdmA>();

    Instruction div_a_cd_a({"div", "am", "cd", "am"});
    inst_op_map_[div_a_cd_a] = std::make_shared<DivACdA>();

    Instruction div_a_a_a({"div", "am", "am", "am"});
    inst_op_map_[div_a_a_a] = std::make_shared<DivAAA>();

    Instruction greater_a_a_b({"gt", "am", "am", "bm"});
    inst_op_map_[greater_a_a_b] = std::make_shared<GreaterAAB>();

    Instruction greater_a_cd_b({"gt", "am", "cd", "bm"});
    inst_op_map_[greater_a_cd_b] = std::make_shared<GreaterACdB>();

    Instruction greater_a_cdm_b({"gt", "am", "cdm", "bm"});
    inst_op_map_[greater_a_cdm_b] = std::make_shared<GreaterACdmB>();

    Instruction less_a_a_b({"lt", "am", "am", "bm"});
    inst_op_map_[less_a_a_b] = std::make_shared<LessAAB>();

    Instruction less_a_cd_b({"lt", "am", "cd", "bm"});
    inst_op_map_[less_a_cd_b] = std::make_shared<LessACdB>();

    Instruction less_a_cdm_b({"lt", "am", "cdm", "bm"});
    inst_op_map_[less_a_cdm_b] = std::make_shared<LessACdmB>();

    Instruction greater_eq_a_a_b({"ge", "am", "am", "bm"});
    inst_op_map_[greater_eq_a_a_b] = std::make_shared<GreatEqualAAB>();

    Instruction greater_eq_a_cdm_b({"ge", "am", "cdm", "bm"});
    inst_op_map_[greater_eq_a_cdm_b] = std::make_shared<GreatEqualACdmB>();

    Instruction equal_a_a_b({"eq", "am", "am", "bm"});
    inst_op_map_[equal_a_a_b] = std::make_shared<EqualAAB>();

    Instruction equal_a_cdm_b({"eq", "am", "cdm", "bm"});
    inst_op_map_[equal_a_cdm_b] = std::make_shared<EqualACdmB>();

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

    Instruction argmax_and_max_a_a_a({"argmax_and_max", "am", "am", "am"});
    inst_op_map_[argmax_and_max_a_a_a] = std::make_shared<RowMajorArgmaxAndMaxAAA>();

    Instruction group_then_sum_by_grouped_count_cdm_a_a({"group_then_sum_by_grouped_count", "cdm", "am", "am"});
    inst_op_map_[group_then_sum_by_grouped_count_cdm_a_a] = std::make_shared<GroupThenSumByGroupedCountCdAA>();

    Instruction quick_sort_a({"quick_sort", "am"});
    inst_op_map_[quick_sort_a] = std::make_shared<QuikSortA>();

    Instruction quick_sort({"quick_sort", "am", "ci", "am"});
    inst_op_map_[quick_sort] = std::make_shared<QuikSortByColumnAA>();

    Instruction split_by_condition({"split_by_condition", "bm", "am", "am", "am"});
    inst_op_map_[split_by_condition] = std::make_shared<SplitByConditionBAAA>();

    Instruction sigmoid({"sigmoid", "am", "am"});
    inst_op_map_[sigmoid] = std::make_shared<SigmoidAA>();

    Instruction groupby_sum({"groupby_sum", "am", "am", "am"});
    inst_op_map_[groupby_sum] = std::make_shared<GroupBySumAA>();

    Instruction groupby_count({"groupby_count", "am", "am", "am"});
    inst_op_map_[groupby_count] = std::make_shared<GroupByCountAA>();

    Instruction groupby_max({"groupby_max", "am", "am", "am"});
    inst_op_map_[groupby_max] = std::make_shared<GroupByMaxAA>();

    Instruction groupby_min({"groupby_min", "am", "am", "am"});
    inst_op_map_[groupby_min] = std::make_shared<GroupByMinAA>();

    Instruction transpose({"transpose", "am", "am"});
    inst_op_map_[transpose] = std::make_shared<TransposeAA>();

    Instruction reshape({"reshape", "am", "ci", "ci", "am"});
    inst_op_map_[reshape] = std::make_shared<Reshape<ArithMatrix>>();

    Instruction reshape_bool({"reshape", "bm", "ci", "ci", "bm"});
    inst_op_map_[reshape_bool] = std::make_shared<Reshape<BoolMatrix>>();

    Instruction set_item({"set_item", "am", "am", "ci", "ci", "ci", "ci"});
    inst_op_map_[set_item] = std::make_shared<SetItemA>();

    Instruction resize({"resize", "am", "ci", "ci", "am"});
    inst_op_map_[resize] = std::make_shared<Resize<ArithMatrix>>();

    Instruction resize_bool({"resize", "bm", "ci", "ci", "bm"});
    inst_op_map_[resize_bool] = std::make_shared<Resize<BoolMatrix>>();
}

DuetVM::DuetVM(const std::shared_ptr<network::Network>& net, std::size_t party_id) : net_(net), party_id_(party_id) {
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

void DuetVM::send_buffer(const std::vector<uint8_t>& input) {
    net_->send_data(input.data(), input.size());
}

std::vector<uint8_t> DuetVM::recv_buffer(std::size_t buffer_size) {
    std::vector<uint8_t> output(buffer_size);
    net_->recv_data(output.data(), buffer_size);
    return output;
}

}  // namespace duet
}  // namespace petace
