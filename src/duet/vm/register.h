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

#include "duet/util/defines.h"
#include "duet/vm/register_content.h"

namespace petace {
namespace duet {

class Register {
public:
    template <typename DataType>
    RegisterAddress set_data(const std::shared_ptr<DataType>& data) {
        current_count_++;
        reg_[current_count_] = std::make_shared<RegisterContent<DataType>>(data);
        return current_count_;
    }

    template <typename DataType>
    const std::shared_ptr<DataType>& get_data(RegisterAddress address) {
        return std::dynamic_pointer_cast<RegisterContent<DataType>>(reg_.at(address))->data();
    }

    void delete_data(RegisterAddress address);

    bool empty();

private:
    RegisterAddress current_count_ = 0;
    std::unordered_map<RegisterAddress, std::shared_ptr<RegisterContentBase>> reg_{};
};

}  // namespace duet
}  // namespace petace
