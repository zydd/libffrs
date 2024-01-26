libffrs
-------

**ffrs** -- Fairly Fast & Flexible Reed-Solomon coding

The core of this library is implemented in C++ and exported to Python with `pybind11 <https://github.com/pybind/pybind11>`_.

For the time being this project has only been tested on Python 3.9.2


Basic usage
-----------

.. code-block:: python

    import ffrs

    message = bytearray(b'Hello World!')

    # Create encoder/decoder for 8 parity bytes
    RS256 = ffrs.RS256(ecc_len=8)

    # Allocate space for correction code
    data = message + bytearray(RS256.ecc_len)

    # Compute parity bytes
    RS256.encode(data)

    # Corrupt some bytes (up to ecc_len/2)
    data[3] ^= 0xaa
    data[9] ^= 0x55

    # Restore partially corrupted message
    RS256.decode(data)

    # Remove parity bytes
    recovered_message = data[:-RS256.ecc_len]


Installation
------------

.. code-block:: bash

    git clone --recursive https://github.com/zydd/libffrs.git
    pip install ./libffrs


Documentation
-------------

Auto-generated documentation can be found at:

- https://zydd.github.io/libffrs/


Building the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

Documentation for this project is generated automatically using Sphinx based on class/method signatures and docstrings.

.. code-block:: bash

    cd ./libffrs
    pip install 'sphinx>=7.2' 'furo==2023.9.10'
    make -C ./docs/ html

Documentation files will be created in ``./build/html``

License
-------

**libffrs** is licensed under the `Apache License, Version 2.0 <https://www.apache.org/licenses/LICENSE-2.0>`_
