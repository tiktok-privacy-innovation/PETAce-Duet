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

#include "duet/util/defines.h"

namespace petace {
namespace duet {

class RegisterContentBase {
public:
    virtual ~RegisterContentBase() = default;
};

template <typename DataType>
class RegisterContent : public RegisterContentBase {
public:
    explicit RegisterContent(const std::shared_ptr<DataType>& value) : value_(value) {
    }

    const std::shared_ptr<DataType>& data() {
        return value_;
    }

private:
    std::shared_ptr<DataType> value_ = nullptr;
};

}  // namespace duet
}  // namespace petace
