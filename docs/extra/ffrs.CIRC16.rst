Normal encoding & interleaving
------------------------------


.. code-block:: python

    RSi16(block_len=8, ecc_len=2)

::

             block_len
         ╭───────┴───────╮
         🟦🟦🟦🟦🟦🟦🟨🟨
         ╰─────┬────╯╰─┬─╯
        message_len  ecc_len



Interleaved encoding
^^^^^^^^^^^^^^^^^^^^


.. code-block:: python

    RSi16(block_len=8, ecc_len=2, interleave=16)


::

                  interleave
    ╭─────────────────┴─────────────────╮
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 ⎫
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 ⎪
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 ⎪
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 ⎪
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 ⎬ block_len
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 ⎪
    🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 ⎪ ⎫ ecc_len
    🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 ⎭ ⎭


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
    🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩 🟫 ⎪ ⎫ outer_ecc_len
    🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩 🟫 ⎭ ⎭

    🟦 = data
    🟨 = inner ECC parity bytes - 1 per row
    🟩 = outer ECC parity bytes - 2 per column
    🟫 = parity bytes for outer ECC - 1 per ECC row

Some data corruption scenarios
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Few random errors:
""""""""""""""""""
::

    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨
    🟦🟦🟦🟦 🟦🟥🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟧 ── inner parity bytes indicate error
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨    locations in the outer block
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟥🟦🟦 🟧 ───╯
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨
    🟩🟩🟩🟩 🟩🟧🟩🟩 🟩🟩🟩🟩 🟧🟩🟩 🟫
    🟩🟩🟩🟩 🟩🟧🟩🟩 🟩🟩🟩🟩 🟧🟩🟩 🟫

    🟥 = corrupted byte
    🟧 = parity check fail

✅ Recovery possible

------------------------------------------------------------------------

Sequential:
"""""""""""
::

    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨
    🟥🟥🟥🟥 🟥🟥🟥🟥 🟥🟥🟥🟥 🟥🟥🟥 🟧 ─ error locations are known
    🟥🟥🟥🟥 🟥🟥🟥🟥 🟥🟥🟥🟥 🟥🟥🟥 🟧 ──╯
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨
    🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧 🟫
    🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧 🟫

✅ Recovery possible

-------------------------------------------------------------------------------

Periodic:
"""""""""
::

    🟥🟦🟦🟦 🟦🟦🟥🟦 🟦🟦🟦🟦 🟥🟦🟦 🟧 ⎫
    🟦🟥🟦🟦 🟦🟦🟦🟥 🟦🟦🟦🟦 🟦🟥🟦 🟧 ⎪
    🟦🟦🟥🟦 🟦🟦🟦🟦 🟥🟦🟦🟦 🟦🟦🟥 🟧 ⎬ all inner check failed
    🟦🟦🟦🟥 🟦🟦🟦🟦 🟦🟥🟦🟦 🟦🟦🟦 🟧 ⎪   cannot locate errors
    🟦🟦🟦🟦 🟥🟦🟦🟦 🟦🟦🟥🟦 🟦🟦🟦 🟧 ⎪
    🟦🟦🟦🟦 🟦🟥🟦🟦 🟦🟦🟦🟥 🟦🟦🟦 🟧 ⎭
    🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧 🟫 ⎫ 2 outer parity bytes can still locate
    🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧 🟫 ⎭ and repair up to 1 error (ecc_len/2)

✅ Recovery possible

------------------------------------------------------------------------

Aligned:
""""""""
::

    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨
    🟥🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟧 ⎫
    🟥🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟧 ⎬ error locations are known but
    🟥🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟧 ⎭ exceed number of parity bytes
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦 🟨
    🟧🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩 🟫
    🟧🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩 🟫

❌ Cannot be recovered

------------------------------------------------------------------------

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
    🟩🟩🟩🟫 🟩🟩🟩🟫 🟩🟩🟩🟫 🟩🟩🟩🟫
    🟩🟩🟩🟫 🟩🟩🟩🟫 🟩🟩🟩🟫 🟩🟩🟩🟫
    ╰─────────────────┬─────────────────╯
        inner_block_len * interleave


