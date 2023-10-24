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

#include <memory>
#include <unordered_map>
#include <vector>

#include "network/net_factory.h"
#include "network/network.h"

#include "duet/duet.h"
#include "duet/util/defines.h"
#include "duet/util/matrix.h"
#include "duet/vm/instruction.h"
#include "duet/vm/operator.h"
#include "duet/vm/register.h"

namespace petace {
namespace duet {

class DuetVM {
public:
    explicit DuetVM(
            const petace::network::NetParams& net_params, petace::network::NetScheme net_type, std::size_t party_id);

    void registr_op(const Instruction& inst, const std::shared_ptr<OperatorBase>& op);

    template <typename DataType>
    RegisterAddress set_data(const std::shared_ptr<DataType>& data) {
        return reg_->set_data(data);
    }

    template <typename DataType>
    RegisterAddress new_data() {
        return reg_->set_data(std::make_shared<DataType>());
    }

    template <typename DataType>
    const std::shared_ptr<DataType>& get_data(RegisterAddress address) {
        return reg_->get_data<DataType>(address);
    }

    void delete_data(RegisterAddress address);

    void exec_code(const Instruction& inst, const std::vector<RegisterAddress>& addr);

    template <typename DataType>
    std::vector<std::size_t> shape(RegisterAddress address) {
        return reg_->get_data<DataType>(address)->shape();
    }

    template <typename DataType>
    void matrix_block(RegisterAddress address_0, RegisterAddress address_1, std::size_t begin_row,
            std::size_t begin_col, std::size_t row_num, std::size_t col_num) {
        const std::shared_ptr<DataType>& x = reg_->get_data<DataType>(address_0);
        const std::shared_ptr<DataType>& y = reg_->get_data<DataType>(address_1);
        petace::duet::matrix_block<DataType>(*x, *y, begin_row, begin_col, row_num, col_num);
    }

    template <typename DataType>
    void private_matrix_block(RegisterAddress address_0, RegisterAddress address_1, std::size_t begin_row,
            std::size_t begin_col, std::size_t row_num, std::size_t col_num) {
        const std::shared_ptr<PrivateMatrix<DataType>>& x = reg_->get_data<PrivateMatrix<DataType>>(address_0);
        const std::shared_ptr<PrivateMatrix<DataType>>& y = reg_->get_data<PrivateMatrix<DataType>>(address_1);
        petace::duet::matrix_block<DataType>(*x, *y, begin_row, begin_col, row_num, col_num, party_id_);
    }

    template <typename DataType>
    void vstack(RegisterAddress address_0, RegisterAddress address_1, RegisterAddress address_2) {
        const std::shared_ptr<DataType>& x = reg_->get_data<DataType>(address_0);
        const std::shared_ptr<DataType>& y = reg_->get_data<DataType>(address_1);
        const std::shared_ptr<DataType>& z = reg_->get_data<DataType>(address_2);
        petace::duet::vstack<DataType>(*x, *y, *z);
    }

    template <typename DataType>
    void private_vstack(RegisterAddress address_0, RegisterAddress address_1, RegisterAddress address_2) {
        const std::shared_ptr<PrivateMatrix<DataType>>& x = reg_->get_data<PrivateMatrix<DataType>>(address_0);
        const std::shared_ptr<PrivateMatrix<DataType>>& y = reg_->get_data<PrivateMatrix<DataType>>(address_1);
        const std::shared_ptr<PrivateMatrix<DataType>>& z = reg_->get_data<PrivateMatrix<DataType>>(address_2);
        petace::duet::vstack<DataType>(*x, *y, *z, party_id_);
    }

    template <typename DataType>
    void hstack(RegisterAddress address_0, RegisterAddress address_1, RegisterAddress address_2) {
        const std::shared_ptr<DataType>& x = reg_->get_data<DataType>(address_0);
        const std::shared_ptr<DataType>& y = reg_->get_data<DataType>(address_1);
        const std::shared_ptr<DataType>& z = reg_->get_data<DataType>(address_2);
        petace::duet::hstack<DataType>(*x, *y, *z);
    }

    template <typename DataType>
    void private_hstack(RegisterAddress address_0, RegisterAddress address_1, RegisterAddress address_2) {
        const std::shared_ptr<PrivateMatrix<DataType>>& x = reg_->get_data<PrivateMatrix<DataType>>(address_0);
        const std::shared_ptr<PrivateMatrix<DataType>>& y = reg_->get_data<PrivateMatrix<DataType>>(address_1);
        const std::shared_ptr<PrivateMatrix<DataType>>& z = reg_->get_data<PrivateMatrix<DataType>>(address_2);
        petace::duet::hstack<DataType>(*x, *y, *z, party_id_);
    }

    bool is_registr_empty();

    std::size_t party_id();

private:
    void registr_default_op();

    std::size_t party_id_ = 0;
    std::shared_ptr<Duet> impl_ = nullptr;
    std::shared_ptr<network::Network> net_ = nullptr;
    std::shared_ptr<Register> reg_ = nullptr;
    std::unordered_map<Instruction, std::shared_ptr<OperatorBase>> inst_op_map_{};
};

}  // namespace duet
}  // namespace petace
