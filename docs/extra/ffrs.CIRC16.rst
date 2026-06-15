Normal encoding & interleaving
------------------------------

With normal RS encoding it's possible to repair up to ``ecc_len`` errors in known locations.
or ``ecc_len/2`` errors in unknown locations.


:Example:

.. code-block:: python

    RSi16(block_len=8, ecc_len=2)

::

             block_len
         ╭───────┴───────╮
         🟦🟦🟦🟦🟦🟦🟩🟩
         ╰─────┬────╯╰─┬─╯
        message_len ecc_len

    🟦 = message
    🟩 = parity data

``RSi16`` interprets the data as 16-bit integers, therefore ``block_len=8`` corresponds to 16 byte blocks.


The encoding above is able to recover 1 random error per 16-byte block (12 message bytes + 4 parity bytes).

------------------------------------------------------------------------

:Example:

    6 random errors over 96 16-bit symbols (192 bytes) encoded with ``RSi16(block_len=8, ecc_len=2)``:

::

    🟦🟦🟥🟦🟦🟦🟧🟧 🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟥🟧🟧 🟦🟦🟦🟦🟦🟦🟩🟩
    🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟥🟧 🟦🟦🟦🟦🟦🟦🟩🟩
    🟦🟦🟦🟦🟦🟦🟧🟥 🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟥🟦🟦🟦🟧🟧 🟦🟦🟦🟦🟥🟦🟧🟧
    🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩

    🟥 = corrupted
    🟧 = parity check fail

✅ Recovery possible

------------------------------------------------------------------------


:This encoding is weak against sequential errors which is a common data loss pattern:

    6 sequential errors over 96 × 16-bit symbols (192 bytes) encoded with ``RSi16(block_len=8, ecc_len=2)``:

::

    🟦🟦🟦🟦🟦🟦🟩🟩 🟥🟥🟥🟥🟥🟥🟧🟧 🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩
    🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩
    🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩
    🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩 🟦🟦🟦🟦🟦🟦🟩🟩


❌ Cannot be recovered

------------------------------------------------------------------------


Interleaved encoding
^^^^^^^^^^^^^^^^^^^^

If we `interleave <https://en.wikipedia.org/wiki/Burst_error-correcting_code#Interleaved_codes>`_ the codec,
i.e.: encode over columns instead of sequentially, it becomes more robust against sequential errors.

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

With this setup, sequential errors appear in different blocks allowing recovery of up to ``interleave`` symbols with the same amount of parity data.

:Example:

    6 sequential errors over 96 × 16-bit symbols (192 bytes) encoded with ``RSi16(block_len=8, ecc_len=2, interleave=16)``:

::

    🟦🟦🟦🟦 🟦🟦🟥🟥 🟥🟥🟥🟥 🟦🟦🟦🟦
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦
    🟩🟩🟩🟩 🟩🟩🟧🟧 🟧🟧🟧🟧 🟩🟩🟩🟩
    🟩🟩🟩🟩 🟩🟩🟧🟧 🟧🟧🟧🟧 🟩🟩🟩🟩

✅ Recovery possible



Cross-interleaving
------------------

If we add an extra parity check on each 'row' of the interleaved code we can double the interleaved codec's ability to repair sequential errors. This *inner*/row parity can be used to locate errors in the *outer*/column blocks. This is the principle behind `Cross-interleaved Reed-Solomon coding <https://en.wikipedia.org/wiki/Cross-interleaved_Reed%E2%80%93Solomon_coding>`_ or CIRC.

.. code-block:: python

    CIRC16(
        inner_block_len=17,
        inner_ecc_len=1,
        outer_block_len=8,
        outer_ecc_len=2
    )

::

                 inner_block_len          ╭─── inner_ecc_len
    ╭───────────────────┴─────────────────╮
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨 ⎫
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨 ⎪
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨 ⎪
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨 ⎪
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨 ⎬ outer_block_len
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨 ⎪
    🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟫 ⎪ ⎫ outer_ecc_len
    🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟫 ⎭ ⎭

    🟦 = data
    🟨 = inner ECC parity data - 1 per row
    🟩 = outer ECC parity data - 2 per column
    🟫 = parity data for outer ECC - 1 per ECC row


::

    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨
    🟥🟥🟥🟥 🟥🟥🟥🟥 🟥🟥🟥🟥 🟥🟥🟥🟥 🟧 ─╮ inner parity check indicates error
    🟥🟥🟥🟥 🟥🟥🟥🟥 🟥🟥🟥🟥 🟥🟥🟥🟥 🟧 ─╯ locations in the outer block
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨
    🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟫 ─ outer blocks can repair
    🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟫   up to 2 errors (ecc_len)

    🟥 = corrupted
    🟧 = parity check fail


✅ Recovery possible

------------------------------------------------------------------------

Other data loss scenarios
^^^^^^^^^^^^^^^^^^^^^^^^^

Random errors:
""""""""""""""
::

    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨
    🟦🟦🟦🟦 🟦🟥🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟧 ── inner parity check indicates error
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨    locations in the outer block
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟥🟦🟦🟦 🟧 ───╯
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨
    🟩🟩🟩🟩 🟩🟧🟩🟩 🟩🟩🟩🟩 🟧🟩🟩🟩 🟫
    🟩🟩🟩🟩 🟩🟧🟩🟩 🟩🟩🟩🟩 🟧🟩🟩🟩 🟫

    🟥 = corrupted
    🟧 = parity check fail

✅ Recovery possible

------------------------------------------------------------------------

Periodic:
"""""""""
::

    🟥🟦🟦🟦 🟦🟦🟥🟦 🟦🟦🟦🟦 🟥🟦🟦🟦 🟧 ⎫
    🟦🟥🟦🟦 🟦🟦🟦🟥 🟦🟦🟦🟦 🟦🟥🟦🟦 🟧 ⎪
    🟦🟦🟥🟦 🟦🟦🟦🟦 🟥🟦🟦🟦 🟦🟦🟥🟦 🟧 ⎬ all inner check failed
    🟦🟦🟦🟥 🟦🟦🟦🟦 🟦🟥🟦🟦 🟦🟦🟦🟥 🟧 ⎪   cannot locate errors
    🟦🟦🟦🟦 🟥🟦🟦🟦 🟦🟦🟥🟦 🟦🟦🟦🟦 🟧 ⎪
    🟦🟦🟦🟦 🟦🟥🟦🟦 🟦🟦🟦🟥 🟦🟦🟦🟦 🟧 ⎭
    🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟫 ⎫ 2 outer parity symbols can still locate
    🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟧🟧🟧🟧 🟫 ⎭ and repair up to 1 error (ecc_len/2)

✅ Recovery possible

------------------------------------------------------------------------

Aligned:
""""""""
::

    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨
    🟥🟥🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟧 ⎫
    🟥🟥🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟧 ⎬ error locations are known but exceed
    🟥🟥🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟧 ⎭ outer code correction ability (outer_ecc_len)
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨
    🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟨
    🟧🟧🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟫
    🟧🟧🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟩🟩🟩🟩 🟫

❌ Cannot be recovered

------------------------------------------------------------------------

Interleaved CIRC
----------------

We can go one step further and interleave the CIRC blocks. This can bring several advantages:
- Allows achieving an effective interleaving greater than 65536 on the input data (maximum ``block_len`` for ``RSi16``)
- More resilience against random errors as the inner codec can also be used for error correction (for ``inner_ecc_len >= 2``)
- Faster encoding with the same amount of effective interleaving and little additional parity data

.. code-block:: python

    CIRC16(
        inner_block_len=4,
        inner_ecc_len=1,
        outer_block_len=6,
        outer_ecc_len=2,
        interleave=5
    )

::

    🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨
    🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨
    🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨
    🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨 🟦🟦🟦🟨
    🟩🟩🟩🟫 🟩🟩🟩🟫 🟩🟩🟩🟫 🟩🟩🟩🟫 🟩🟩🟩🟫
    🟩🟩🟩🟫 🟩🟩🟩🟫 🟩🟩🟩🟫 🟩🟩🟩🟫 🟩🟩🟩🟫
    ╰─────────────────────┬─────────────────────╯
          inner_block_len * interleave


------------------------------------------------------------------------

ffrs.CIRC16
-----------