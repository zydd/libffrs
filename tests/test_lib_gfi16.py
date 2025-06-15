
#  test_lib_gfi16.py
#
#  Copyright 2025 Gabriel Machado
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

import itertools
import pytest
import random

import ffrs
import ffrs.reference.ntt as ref_ntt
from ffrs.reference.util import *


GF = ffrs.GFi16(3)
GF_ref = ffrs.reference.GF(65537, 1, 3)
random.seed(42)


@pytest.mark.parametrize('fn, id', [
    (GF.add, 0),
    (GF.sub, 0),
    (GF.mul, 1),
    (GF.div, 1),
])
def test_scalar_identity(fn, id):
    for a in sample_field(GF):
        assert fn(a, id) == a, a


@pytest.mark.parametrize('fn1, fn2', [
    (GF.exp, GF.log),
])
def test_scalar_inverse1(fn1, fn2):
    for a in sample_field(GF, start=1):
        assert fn2(fn1(a)) == a, a


@pytest.mark.parametrize('fn1, fn2', [
    (GF.add, GF.sub),
    (GF.mul, GF.div),
])
def test_scalar_inverse2(fn1, fn2):
    for a, b in itertools.product(sample_field(GF), sample_field(GF, start=1)):
        res = fn1(a, b)
        res_inv = fn2(res, b)
        assert res_inv == a, (a, b)


@pytest.mark.parametrize('fn', [
    GF.add,
    GF.mul,
])
def test_scalar_commutativity(fn):
    for a, b in itertools.product(sample_field(GF), sample_field(GF)):
        assert fn(a, b) == fn(b, a)


def test_scalar_inv_div():
    for a in sample_field(GF, start=1):
        assert GF.inv(a) == GF.div(1, a) != 0


def test_scalar_inv_inv():
    for a in sample_field(GF, start=1):
        assert GF.inv(a) == GF.div(1, a) != 0


def test_scalar_exp_pow():
    for a in sample_field(GF):
        assert GF.exp(a) == GF.pow(GF.primitive, a) != 0

