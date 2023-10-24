# Copyright 2023 TikTok Pte. Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import unittest
from multiprocessing import Process
import numpy as np

from petace.duet.pyduet import NetParams, DuetVM, NetScheme
from petace.duet.securenumpy import Private, Public
from petace.duet.securenumpy import vstack, hstack, multiplexer
from petace.duet.securenumpy import row_major_argmax_and_max


class TestDuetVM(unittest.TestCase):
    def add_p_p(self, vm):
        data_a = np.empty((0, 0))
        data_b = np.array([[9, 9], [9, 9]], dtype=np.float64)
        if vm.party_id() == 0:
            data_a = np.array([[2, 2], [3, 4]], dtype=np.float64)

        p0 = Private(data_a, 0, vm)
        p1 = Private(data_b, 1, vm)
        a0 = p0 + p1
        p2 = a0.reveal_to(0)
        if vm.party_id() == 0:
            self.assertTrue(np.allclose(
                data_a + data_b, p2.to_numpy(), atol=0.01))

    def add_p_a(self, vm):
        data_a = np.empty((0, 0))
        data_b = np.array([[9, 9], [9, 9]], dtype=np.float64)
        if vm.party_id() == 0:
            data_a = np.array([[2, 2], [3, 4]], dtype=np.float64)

        p0 = Private(data_a, 0, vm)
        p1 = Private(data_b, 1, vm)
        a1 = p1.to_share()
        a2 = p0 + a1
        p2 = a2.reveal_to(0)
        if vm.party_id() == 0:
            self.assertTrue(np.allclose(
                data_a + data_b, p2.to_numpy(), atol=0.01))

    def mul_a_c(self, vm):
        data_a = np.empty((0, 0))
        if vm.party_id() == 0:
            data_a = np.array([[2, 2], [3, 4]], dtype=np.float64)
        p0 = Private(data_a, 0, vm)
        a0 = p0.to_share()
        a1 = a0 * 2.0
        p1 = a1.reveal_to(0)

        if vm.party_id() == 0:
            self.assertTrue(np.allclose(
                data_a * 2.0, p1.to_numpy(), atol=0.001))

    def gt_p_a(self, vm):
        data_a = np.empty((0, 0))
        data_b = np.array([[9, 9], [9, 9]], dtype=np.float64)
        if vm.party_id() == 0:
            data_a = np.array([[2, 2], [3, 4]], dtype=np.float64)
        p0 = Private(data_a, 0, vm)
        p1 = Private(data_b, 1, vm)
        b0 = p0 > p1
        p2 = b0.reveal_to(0)
        if vm.party_id() == 0:
            self.assertTrue(np.array_equal(data_a > data_b, p2.to_numpy()))

    def slice(self, vm):
        data_a = np.empty((0, 0))
        if vm.party_id() == 0:
            data_a = np.array([[2, 2], [3, 4]], dtype=np.float64)
        p0 = Private(data_a, 0, vm)
        a0 = p0.to_share()
        a1 = a0[:, :1]
        p1 = a1.reveal_to(0)
        if vm.party_id() == 0:
            self.assertTrue(np.allclose(
                data_a[:, :1], p1 .to_numpy(), atol=0.001))

    def stack(self, vm):
        data_a = np.empty((0, 0))
        data_b = np.array([[9, 9], [9, 9]], dtype=np.float64)
        if vm.party_id() == 0:
            data_a = np.array([[2, 2], [3, 4]], dtype=np.float64)
        p0 = Private(data_a, 0, vm)
        p1 = Private(data_b, 1, vm)
        a0 = p0.to_share()
        a1 = p1.to_share()
        a2 = vstack(a0, a1)
        a3 = hstack(a0, a1)
        p2 = a2.reveal_to(0)
        p3 = a3.reveal_to(0)
        if vm.party_id() == 0:
            self.assertTrue(np.allclose(
                np.vstack((data_a, data_b)), p2.to_numpy(), atol=0.001))
            self.assertTrue(np.allclose(
                np.hstack((data_a, data_b)), p3.to_numpy(), atol=0.001))

    def div_a_a(self, vm):
        data_a = np.empty((0, 0))
        data_b = np.array([[9, 9], [9, 9]], dtype=np.float64)
        if vm.party_id() == 0:
            data_a = np.array([[2, 2], [3, 4]], dtype=np.float64)
        p0 = Private(data_a, 0, vm)
        p1 = Private(data_b, 1, vm)
        a0 = p0.to_share()
        a1 = p1.to_share()
        a2 = a0 / a1
        p2 = a2.reveal_to(0)
        if vm.party_id() == 0:
            self.assertTrue(np.allclose(
                data_a / data_b, p2.to_numpy(), atol=0.001))

    def multiplexer_b_a_a_a(self, vm):
        data_a = np.empty((0, 0))
        data_b = np.array([[9, 9], [9, 9]], dtype=np.float64)
        if vm.party_id() == 0:
            data_a = np.array([[2, 2], [3, 4]], dtype=np.float64)
        p0 = Private(data_a, 0, vm)
        p1 = Private(data_b, 1, vm)
        a0 = p0.to_share()
        a1 = p1.to_share()
        a2 = a0 > a1
        a3 = multiplexer(a2, a1, a0)
        p2 = a3.reveal_to(0)
        if vm.party_id() == 0:
            self.assertTrue(np.allclose(np.maximum(
                data_a, data_b), p2.to_numpy(), atol=0.001))

    def row_major_argmax_and_max_a(self, vm):
        data_a = np.empty((0, 0))
        if vm.party_id() == 0:
            data_a = np.array([[2, 2], [3, 4]], dtype=np.float64)
        p0 = Private(data_a, 0, vm)
        a0 = p0.to_share()
        a1, a2 = row_major_argmax_and_max(a0)
        p1 = a1.reveal_to(0)
        p2 = a2.reveal_to(0)
        if vm.party_id() == 0:
            self.assertTrue(np.allclose(np.argmax(data_a, axis=0).reshape(
                (1, 2)), p1.to_numpy().astype(np.int64), atol=0.001))
            self.assertTrue(np.allclose(
                np.max(data_a, axis=0), p2.to_numpy(), atol=0.001))

    def vstack(self, vm):
        data_a = [[0, 1], [1, 1], [1, 1], [1, 1], [0, 1], [0, 5]]
        data_b = [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]]
        if vm.party_id() == 0:
            data = data_a
        else:
            data = data_b
        data_array = np.array(data, dtype=np.float64)
        max_list = data_array[:, 0].reshape((1, -1))
        max_list_priv1 = Private(max_list, 0, vm).to_share()
        max_list_priv2 = Private(max_list, 1, vm).to_share()
        max_list_priv = vstack(max_list_priv1, max_list_priv2)
        res_plain = max_list_priv.reveal_to(0).to_numpy()
        if vm.party_id() == 0:
            data_a_array = np.array(data_a, dtype=np.float64)[
                :, 0].reshape((1, -1))
            data_b_array = np.array(data_b, dtype=np.float64)[
                :, 0].reshape((1, -1))
            self.assertTrue(np.allclose(
                np.vstack([data_a_array, data_b_array]), res_plain, atol=0.001))

    def or_pub(self, vm):
        data_c = np.array([[True, True], [False, False]], dtype=np.bool_)
        data_a = np.empty((0, 0))
        data_b = np.array([[9, 9], [9, 9]], dtype=np.float64)
        if vm.party_id() == 0:
            data_a = np.array([[2, 2], [3, 4]], dtype=np.float64)
        p0 = Private(data_a, 0, vm)
        p1 = Private(data_b, 1, vm)
        c0 = Public(data_c, vm)
        a0 = p0.to_share()
        a1 = p1.to_share()
        b0 = a0 > a1
        b1 = b0 | c0
        p2 = b1.reveal_to(0)
        if vm.party_id() == 0:
            print(p2.to_numpy())
            print(np.bitwise_or((data_a > data_b), data_c))

    def run_process(self, party):
        net_params = NetParams()
        if party == 0:
            net_params.remote_addr = "127.0.0.1"
            net_params.remote_port = 8890
            net_params.local_port = 8891
        else:
            net_params.remote_addr = "127.0.0.1"
            net_params.remote_port = 8891
            net_params.local_port = 8890

        duet = DuetVM(net_params, NetScheme.SOCKET, party)

        self.add_p_p(duet)
        self.gt_p_a(duet)
        self.add_p_a(duet)
        self.mul_a_c(duet)
        self.slice(duet)
        self.stack(duet)
        self.div_a_a(duet)
        self.vstack(duet)
        self.multiplexer_b_a_a_a(duet)
        self.row_major_argmax_and_max_a(duet)
        self.or_pub(duet)
        print(party, duet.is_registr_empty())

    def test_vm_init(self):

        processes = []

        for i in range(2):
            p = Process(target=self.run_process, args=(i,))
            p.start()
            processes.append(p)

        for p in processes:
            p.join()


if __name__ == '__main__':
    unittest.main()
