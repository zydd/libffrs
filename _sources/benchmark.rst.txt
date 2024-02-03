.. _benchmarks:

Benchmarks
==========


Encoding
--------

**libffrs** `v0.1-11-g5bf2eb3 <https://github.com/zydd/libffrs/tree/v0.1>`_


:py:meth:`ffrs.RS256.encode_blocks`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
    :file: benchmark_encode_blocks.csv


Decoding
--------

TBD


Methodology
-----------

The encoding speeds were measured by repeatedly encoding 100 MB buffers over a period of 5 seconds.
A minimum of 5 measurements taken, then the maximum encoding speed is used.

See `tests/benchmark.py <https://github.com/zydd/libffrs/blob/master/tests/benchmark.py>`_.

System specs
^^^^^^^^^^^^

Processor: 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz

| ``gcc``/``Clang`` results were measured on WSL2 Linux VM
| ``MSVC``/``Clang-cl`` results were measured on Windows 11
