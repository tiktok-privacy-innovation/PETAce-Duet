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

from typing import Union, Tuple
import operator

import numpy as np

from petace.duet.pyduet import DuetVM, Instruction
from petace.duet.exception import SecureArrayException


class Private:
    """
    The data belongs to one party privately.

    This is the data perspective of the Petace-Duet.

    Attributes:
    reg_addr (int64_t): The data address in vm impl.
    impl (VM): The VM instance. Now only support the duet vm.
    party (int): Which party have this data.
    data_type (str): The type of the underline data.

    Args:
    data (np.ndarray, None): Now only support np.ndarray.
    party (int): Which party have this data.
    impl (VM): The VM instance. Now only support the duet vm.
    data_type (str): The type of the underline data. Must set if the data is None.
    """
    reg_addr = None
    impl = None
    party = 0
    data_type = None

    def __init__(self, data: Union[np.ndarray, None], party: int, impl: "DuetVM", data_type: str = None):
        if party is None:
            raise SecureArrayException("must set party id for a private var")
        if data is None:
            if data_type is None:
                raise SecureArrayException(
                    "must set data_type unless set data")
            else:
                self.data_type = data_type
        else:
            if isinstance(data, np.ndarray):
                if len(data.shape) != 2:
                    raise SecureArrayException("only support 2 dim data")
                if data.dtype == np.float64:
                    self.data_type = "pdm"
                elif data.dtype == np.bool_:
                    self.data_type = "pbm"
                else:
                    raise SecureArrayException("unsupported data type")
            else:
                raise SecureArrayException("unsupported data type")

        if self.data_type == "pdm":
            self.reg_addr = impl.new_private_double_matrix(party)
        elif self.data_type == "pbm":
            self.reg_addr = impl.new_private_bool_matrix(party)
        else:
            raise SecureArrayException("unsupported data type")

        if data is not None:
            impl.set_private_double_matrix(data, self.reg_addr)

        self.impl = impl
        self.party = party

    def __del__(self):
        self.impl.delete_data(self.reg_addr)

    def to_numpy(self) -> np.ndarray:
        """
        Transport a Private object to a numpy array.
        """
        if self.data_type == "pdm":
            return self.impl.get_private_double_matrix(self.reg_addr)
        elif self.data_type == "pbm":
            return self.impl.get_private_bool_matrix(self.reg_addr)

    def to_share(self) -> "Share":
        """
        Transport a Private object to a Share object.
        """
        share = None
        if self.data_type == "pdm":
            share = Share(self.impl, "am")
        elif self.data_type == "pbm":
            share = Share(self.impl, "bm")
        inst = Instruction(["share", self.data_type, share.data_type])
        self.impl.exec_code(inst, [self.reg_addr, share.reg_addr])
        return share

    def _operator(self, other: Union["Private", "Share"], operation: operator) -> "Share":
        if isinstance(other, Share):
            return operation(self.to_share(), other)
        elif isinstance(other, Private):
            return operation(self.to_share(), other.to_share())
        else:
            raise SecureArrayException(f"unsupported type {type(other)}")

    def __add__(self, other) -> "Share":
        """
        ret = self + other

        Returns:
        The result Share object.
        """
        return self._operator(other, operator.add)

    def __sub__(self, other) -> "Share":
        """
        ret = self - other

        Returns:
        The result Share object.
        """
        return self._operator(other, operator.sub)

    def __mul__(self, other) -> "Share":
        """
        ret = self * other

        Returns:
        The result Share object.
        """
        return self._operator(other, operator.mul)

    def __gt__(self, other) -> "Share":
        """
        ret = self > other

        Returns:
        The result Share object.
        """
        return self._operator(other, operator.gt)

    def __lt__(self, other) -> "Share":
        """
        ret = self < other

        Returns:
        The result Share object.
        """
        return self._operator(other, operator.lt)

    def __ge__(self, other) -> "Share":
        """
        ret = self >= other

        Returns:
        The result Share object.
        """
        return self._operator(other, operator.ge)

    def __eq__(self, other) -> "Share":
        """
        ret = self == other

        Returns:
        The result Share object.
        """
        return self._operator(other, operator.eq)


class Share:
    """
    The Share data object.

    This is the data perspective of the Petace-Duet.

    Attributes:
    reg_addr (int64_t): The data address in vm impl.
    impl (VM): The VM instance. Now only support the duet vm.
    data_type (str): The type of the underline data.

    Args:
    impl (VM): The VM instance. Now only support the duet vm.
    data_type (str): The type of the underline data. Must set if the data is None.
    share (np.ndarray, None): Reserved interface. To set share data.
    """
    reg_addr = None
    impl = None
    data_type = None

    def __init__(self, impl: "DuetVM", data_type: str = None, share: np.ndarray = None):
        if share is not None:
            raise SecureArrayException("not support set share now")
        if data_type == "am":
            self.reg_addr = impl.new_airth_matrix()
        elif data_type == "bm":
            self.reg_addr = impl.new_bool_matrix()
        else:
            raise SecureArrayException("unsupported share type")
        self.data_type = data_type
        self.impl = impl

    def __del__(self):
        self.impl.delete_data(self.reg_addr)

    def reveal_to(self, party) -> "Private":
        """
        Transport a Share object to a Private object.

        Args:
        party (int): Which party have this data.
        """
        private = None
        if self.data_type == "am":
            private = Private(None, party, self.impl, "pdm")
        elif self.data_type == "bm":
            private = Private(None, party, self.impl, "pbm")
        else:
            raise SecureArrayException("unsupported share type")
        inst = Instruction(["reveal", self.data_type, private.data_type])
        self.impl.exec_code(inst, [self.reg_addr, private.reg_addr])
        return private

    def reveal(self) -> "Public":
        """
        Transport a Share object to a Public object.
        """
        public = None
        if self.data_type == "am":
            public = Public(None, self.impl, "cdm")
        else:
            raise SecureArrayException("unsupported share type")
        inst = Instruction(["reveal", self.data_type, public.data_type])
        self.impl.exec_code(inst, [self.reg_addr, public.reg_addr])
        return public

    def _operator(self, other: Union["Private", "Share", "Public", "float"], ret: "Share", operation: operator):
        if isinstance(other, Private):
            return operation(self, other.to_share())
        elif isinstance(other, Share) or isinstance(other, Public):
            inst = Instruction(
                [operation.__name__, self.data_type, other.data_type, ret.data_type])
            self.impl.exec_code(
                inst, [self.reg_addr, other.reg_addr, ret.reg_addr])
        elif isinstance(other, float):
            return operation(self, Public(other, self.impl))
        else:
            raise SecureArrayException(f"unsupported type {type(other)}")
        return ret

    def __add__(self, other: Union["Private", "Share"]) -> "Share":
        """
        ret = self + other

        Returns:
        The result Share object.
        """
        ret = Share(self.impl, "am")
        return self._operator(other, ret, operator.add)

    def __sub__(self, other: Union["Private", "Share", "float"]) -> "Share":
        """
        ret = self - other

        Returns:
        The result Share object.
        """
        ret = Share(self.impl, "am")
        return self._operator(other, ret, operator.sub)

    def __mul__(self, other: Union["Private", "Share", "Public", "float"]) -> "Share":
        """
        ret = self * other

        Returns:
        The result Share object.
        """
        ret = Share(self.impl, "am")
        return self._operator(other, ret, operator.mul)

    def __truediv__(self, other: Union["Share", "Public"]):
        """
        ret = self / other

        Returns:
        The result Share object.
        """
        if isinstance(other, Share) or isinstance(other, Public):
            ret = Share(self.impl, self.data_type)
            inst = Instruction(
                ["div", self.data_type, other.data_type, ret.data_type])
            self.impl.exec_code(
                inst, [self.reg_addr, other.reg_addr, ret.reg_addr])
            return ret
        else:
            raise SecureArrayException("unsupported share type")

    def __gt__(self, other: Union["Private", "Share"]) -> "Share":
        """
        ret = self > other

        Returns:
        The result Share object.
        """
        ret = Share(self.impl, "bm")
        return self._operator(other, ret, operator.gt)

    def __lt__(self, other: Union["Private", "Share"]) -> "Share":
        """
        ret = self < other

        Returns:
        The result Share object.
        """
        ret = Share(self.impl, "bm")
        return self._operator(other, ret, operator.lt)

    def __ge__(self, other: Union["Private", "Share"]) -> "Share":
        """
        ret = self >= other

        Returns:
        The result Share object.
        """
        ret = Share(self.impl, "bm")
        return self._operator(other, ret, operator.ge)

    def __eq__(self, other: Union["Private", "Share"]) -> "Share":
        """
        ret = self == other

        Returns:
        The result Share object.
        """
        ret = Share(self.impl, "bm")
        return self._operator(other, ret, operator.eq)

    def __ror__(self, other: "Public") -> "Share":
        """
        ret = self | other

        Returns:
        The result Share object.
        """
        ret = Share(self.impl, "bm")
        inst = Instruction(
            ["or", other.data_type, self.data_type, ret.data_type])
        self.impl.exec_code(
            inst, [other.reg_addr, self.reg_addr, ret.reg_addr])
        return ret

    def __or__(self, other: "Public") -> "Share":
        """
        ret = self | other

        Returns:
        The result Share object.
        """
        return other | self

    def __getitem__(self, index) -> "Share":
        """
        Get a submatrix.

        Args:
        index (tuple): The index of the submatrix, such as a[0: 1, 0: 1]

        Returns:
        The result submatrix.
        """
        if isinstance(index, tuple):
            if len(index) != 2:
                raise SecureArrayException(f"unsupported index type")
            else:
                if isinstance(index[0], slice) and isinstance(index[1], slice):
                    if index[0].step is not None or index[1].step is not None:
                        raise SecureArrayException(f"unsupported index type")
                    else:
                        shape = self.impl.get_airth_share_matrix_shape(
                            self.reg_addr)
                        begin_row = 0
                        begin_col = 0
                        row_num = shape[0]
                        col_num = shape[1]

                        if index[0].start is not None:
                            begin_row = index[0].start
                        if index[0].stop is not None:
                            row_num = index[0].stop - begin_row
                        if index[1].start is not None:
                            begin_col = index[1].start
                        if index[1].stop is not None:
                            col_num = index[1].stop - begin_col
                        ret = Share(self.impl, self.data_type)
                        self.impl.airth_share_matrix_block(self.reg_addr, ret.reg_addr, begin_row, begin_col, row_num,
                                                           col_num)
                        return ret

                else:
                    raise SecureArrayException(f"unsupported index type")
        else:
            raise SecureArrayException(f"unsupported index type")


class Public:
    """
    The Public data object.

    This is the data perspective of the Petace-Duet.

    Attributes:
    reg_addr (int64_t): The data address in vm impl.
    impl (VM): The VM instance. Now only support the duet vm.
    data_type (str): The type of the underline data.

    Args:
    data (np.ndarray, float, None): The data.
    impl (VM): The VM instance. Now only support the duet vm.
    data_type (str): The type of the underline data. Must set if the data is None.
    """
    reg_addr = None
    impl = None
    data_type = None

    def __init__(self, data: Union[np.ndarray, float, None], impl: "DuetVM", data_type: str = None):
        if data is None:
            if data_type is None:
                raise SecureArrayException(
                    "must set data_type unless set data")
            else:
                self.data_type = data_type
        else:
            if isinstance(data, np.ndarray):
                if len(data.shape) != 2:
                    raise SecureArrayException("only support 2 dim data")
                if data.dtype == np.float64:
                    self.data_type = "cdm"
                elif data.dtype == np.bool_:
                    self.data_type = "cbm"
                else:
                    raise SecureArrayException("unsupported data type")
            elif isinstance(data, float):
                self.data_type = "cd"
            else:
                raise SecureArrayException("unsupported data type")

        if self.data_type == "cdm":
            self.reg_addr = impl.new_public_double_matrix()
            if data is not None:
                impl.set_public_double_matrix(data, self.reg_addr)
        elif self.data_type == "cd":
            self.reg_addr = impl.new_public_double()
            if data is not None:
                impl.set_public_double(data, self.reg_addr)
        elif self.data_type == "cbm":
            self.reg_addr = impl.new_public_bool_matrix()
            if data is not None:
                impl.set_public_bool_matrix(data, self.reg_addr)
        else:
            raise SecureArrayException("unsupported data type")
        self.impl = impl

    def __del__(self):
        self.impl.delete_data(self.reg_addr)


def vstack(first: Share, second: Share) -> Share:
    """
    Stack two share object vertically.

    Returns:
    The result Share object.
    """
    ret = Share(first.impl, first.data_type)
    first.impl.airth_share_vstack(
        first.reg_addr, second.reg_addr, ret.reg_addr)
    return ret


def hstack(first: Share, second: Share) -> Share:
    """
    Stack two share object horizontally.

    Returns:
    The result Share object.
    """
    ret = Share(first.impl, first.data_type)
    first.impl.airth_share_hstack(
        first.reg_addr, second.reg_addr, ret.reg_addr)
    return ret


def multiplexer(cond: "Share", a: "Share", b: "Share") -> "Share":
    """
    if cond == 0 then retrun a else return b

    Args:
    cond Share: The boolenshare, 0/1.

    Returns:
    The result Share object.
    """
    ret = Share(cond.impl, a.data_type)
    inst = Instruction(["multiplexer", cond.data_type,
                       a.data_type, b.data_type, ret.data_type])
    cond.impl.exec_code(
        inst, [cond.reg_addr, a.reg_addr, b.reg_addr, ret.reg_addr])
    return ret


def row_major_argmax_and_max(share: "Share") -> Tuple["Share", "Share"]:
    """
    Row major argmax and max. Better performance than computing separately.

    Returns:
    max_index: The max index of this matrix(argmax).
    max_value: The max value of this matrix(max).
    """
    max_index = Share(share.impl, share.data_type)
    max_value = Share(share.impl, share.data_type)
    inst = Instruction(["row_major_argmax_and_max", share.data_type,
                       max_index.data_type, max_value.data_type])
    share.impl.exec_code(
        inst, [share.reg_addr, max_index.reg_addr, max_value.reg_addr])
    return max_index, max_value


def sum(share: "Share") -> "Share":
    """
    Sum of all values column-wise in a arithmetic secret shared matrix.

    Returns:
    The result Share object.
    """
    ret = Share(share.impl, share.data_type)
    inst = Instruction(["sum", share.data_type, ret.data_type])
    share.impl.exec_code(inst, [share.reg_addr, ret.reg_addr])
    return ret
