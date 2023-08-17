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

#include "network/net_socket.h"

#include "duet/duet.h"

void add_test(const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet);

void mul_test(const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet);

void plain_shuffle_test(
        const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet);

void less_than_zero_test(
        const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet);

void millionaires(
        const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet);

void ppml(const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet);

void distance(const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet);
