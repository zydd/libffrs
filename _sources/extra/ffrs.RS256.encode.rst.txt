:Example:
    - Considering ``RS256(6, 4)`` and 32 bytes of input:

    .. code-block:: python

        import ffrs, random
        rs62 = ffrs.RS256(6, 4)
        input_buffer = random.randbytes(32)
        output_buffer = rs62.encode_blocks(input_buffer)

    ::

        Input buffer (32 bytes):
        🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦
        🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦

        Output buffer (32 data bytes + 16 parity bytes):
        🟦🟦🟦🟦🟨🟨🟦🟦🟦🟦🟨🟨🟦🟦🟦🟦🟨🟨🟦🟦🟦🟦🟨🟨
        🟦🟦🟦🟦🟨🟨🟦🟦🟦🟦🟨🟨🟦🟦🟦🟦🟨🟨🟦🟦🟦🟦🟨🟨


.. note::
    To allow concatenating the results of multiple calls to :py:meth:`encode_blocks`,
    the size of ``buffer`` should be a multiple of :py:attr:`RS256.message_len`
    (or ``block_size - RS256.ecc_len`` if the argument ``block_size`` is not ``None``)

If the size of the input buffer is not a multiple of :py:attr:`RS256.message_len`
the remaining bytes at the end of the buffer will be encoded as if the block size were smaller.

:Example:
    - Considering ``RS256(6, 4)`` and 29 bytes of input:

    ::

        Input buffer (29 bytes):
        🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦
        🟦🟦🟦🟦 🟦🟦🟦🟦 🟦🟦🟦🟦 🟦

        Output buffer (29 data bytes + 16 parity bytes):
        🟦🟦🟦🟦🟨🟨🟦🟦🟦🟦🟨🟨🟦🟦🟦🟦🟨🟨🟦🟦🟦🟦🟨🟨
        🟦🟦🟦🟦🟨🟨🟦🟦🟦🟦🟨🟨🟦🟦🟦🟦🟨🟨🟦🟨🟨
