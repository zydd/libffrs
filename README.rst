libffrs
=======

Fairly *F*\ ast & *F*\ lexible *R*\ eed-*S*\ olomon Coding

This library uses a variation Intel's Slicing-by-8 algorithm for the encoding routines,
achieving some hefty speeds depending on the code sizes. See `Benchmarks <https://zydd.github.io/libffrs/benchmark.html>`_.

The core implementation is done in C++ and exported to Python with `pybind11 <https://github.com/pybind/pybind11>`_.


Basic usage
===========

.. code-block:: python

    import ffrs

    message = bytearray(b'Hello World!')

    # Create encoder/decoder for 8 parity bytes
    rs8 = ffrs.RS256(ecc_len=8)

    # Allocate space for correction code
    data = message + bytearray(rs8.ecc_len)

    # Compute parity bytes
    rs8.encode(data)

    # Corrupt some bytes (up to ecc_len/2)
    data[3] ^= 0xaa
    data[9] ^= 0x55

    # Restore partially corrupted message
    rs8.decode(data)

    # Remove parity bytes
    recovered_message = data[:-rs8.ecc_len]

    assert message == recovered_message



Installation
============

.. code-block:: bash

    git clone --recursive https://github.com/zydd/libffrs.git
    pip install ./libffrs


Documentation
=============

Auto-generated documentation can be found at:

- https://zydd.github.io/libffrs/


Building the documentation
--------------------------

Documentation for this project is generated automatically using Sphinx based on class/method signatures and doc strings.

.. code-block:: bash

    cd ./libffrs
    pip install 'sphinx>=7.2' 'furo==2023.9.10'
    make -C ./docs/ html

Documentation files will be created in ``./build/html``

License
=======

**libffrs** is licensed under the `Apache License, Version 2.0 <https://www.apache.org/licenses/LICENSE-2.0>`_
