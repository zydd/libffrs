Cross-interleaving
------------------

.. code-block:: python

    CIRC16(
        inner_block_len=16,
        inner_ecc_len=1,
        outer_block_len=8,
        outer_ecc_len=2
    )

::

                inner_block_len          ╭─── inner_ecc_len
    ╭──────────────────┴─────────────────╮
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨 ⎫
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨 ⎪
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨 ⎪
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨 ⎪
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨 ⎬ outer_block_len
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨 ⎪
    🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧 🟫 ⎪ ⎫ outer_ecc_len
    🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧 🟫 ⎭ ⎭ 

    🟦 = data
    🟨 = inner ECC parity bytes - 1 per row
    🟧 = outer ECC parity bytes - 2 per column
    🟫 = parity bytes for outer ECC - 1 per ECC row


Interleaved CIRC
----------------

.. code-block:: python

    CIRC16(
        inner_block_len=4,
        inner_ecc_len=1,
        outer_block_len=6,
        outer_ecc_len=2
    )

::

    🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨
    🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨
    🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨
    🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨
    🟧🟧🟧🟫 🟧🟧🟧🟫 🟧🟧🟧🟫 🟧🟧🟧🟫
    🟧🟧🟧🟫 🟧🟧🟧🟫 🟧🟧🟧🟫 🟧🟧🟧🟫
    ╰──────────────────┬─────────────────╯
        inner_block_len * interleave


