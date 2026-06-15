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


:Examples:
    - ``CIRC16(2048+4, 4, 2048+128, 128)``

        Creates an encoder with 2048 interleaved blocks of 4 KiB (8 MiB in total)
        
        This codec can repair 500 KiB in the best case scenario.

    - ``CIRC16(2048+4, 4, 4096+256, 256, 128)``

        Creates an encoder with 4096×128 interleaved blocks of 4 KiB (2 GiB in total)

        This codec can repair 128 MiB in the best case scenario.

:Parameters:
    - ``inner_block_len``
        | Number of 16-bit elements per row

        Must be a multiple of ``inner_ecc_len``

    - ``inner_ecc_len``
        | Number of parity 16-bit elements per row

        Must be a power of 2

    - ``outer_block_len``
        | Number of rows per CIRC block

        Must be a multiple of ``outer_ecc_len``

    - ``outer_ecc_len``
        | Number of parity rows per CIRC block

        Must be a power of 2

    - ``interleave``
        | Number of interleaved CIRC blocks
