:Examples:
    - ``RS256(255, 223)``
        Equivalent to ``RS256(ecc_len=32)``

        Creates an encoder for 32 bytes of parity,
        capable of correcting up to 16 errors in a 255-byte block.

    - ``RS256(6, 4)``
        ::

                block_len
             ╭──────┴──────╮
             🟦🟦🟦🟦🟨🟨
             ╰────┬───╯╰─┬─╯
            message_len ecc_len

:Parameters:
    - block_len
        | :math:`n` -- default block size used by :py:meth:`encode_blocks`
        | ``block_len = message_len + ecc_len``
        | Defaults to ``255``

        .. note::
            The maximum possible ``block_len`` for :py:class:`RS256` is 255 bytes

    - message_len
        | :math:`k` -- number of actual data bytes in a block
        | Can be omitted if ``ecc_len`` is supplied

    - ecc_len
        | :math:`(n - k)` -- number of parity bytes in a block
        | Can be omitted if ``message_len`` is supplied

    - primitive
        | :math:`a` -- primitive value for :py:class:`GF256`
        | Defaults to ``2``

    - polynomial
        | :math:`P` -- irreducible polynomial for :py:class:`GF256`
        | Defaults to ``0x11b``
