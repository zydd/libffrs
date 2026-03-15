for src_size in [48]:
    src = list(range(src_size))
    dst_cols = 8
    dst_rows = 6
    interleave = 2
    dst = [0] * ((dst_rows + 2) * dst_cols)


    def interleave1(src, dst):
        for i in range(len(src)):
            col = i % interleave + interleave * (i // (interleave * dst_rows))
            row = (i // interleave) % dst_rows
            assert dst[row * dst_cols + col] == 0
            dst[row * dst_cols + col] = src[i]

    def interleave2(src, dst):
        count = 0
        for int_col in range(0, dst_cols, interleave):
            max_rows = min(dst_rows, len(src) // interleave + bool(len(src) % interleave))
            assert max_rows == min(dst_rows, round(len(src) / interleave + 0.5))
            for row in range(max_rows):
                for col in range(min(interleave, dst_cols)):
                    src_pos = int_col * dst_rows + row * interleave + col
                    dst_pos = int_col            + row * dst_cols   + col
                    assert dst[dst_pos] == 0
                    assert src_pos == count

                    if src_pos >= len(src):
                        return

                    dst[dst_pos] = src[src_pos]

                    count += 1

    def copy_msg_transposed(src, dst, vec_size, msg_size):
        for j in range(vec_size):
            for i in range(msg_size):
                dst[i * vec_size + j] = src[i + j * msg_size]

    def deinterleave(src, dst):
        interleave = 4
        count = 0
        for j in range(dst_cols // interleave):
            int_col = j * interleave
            for i in range(dst_rows * interleave):
                col = i % interleave
                row = i // interleave
                dst[count] = src[int_col + row * dst_cols + col]
                count += 1
                if count >= len(src):
                    return

    try:
        interleave2(src, dst)
        dst2 = list(dst)
        # copy_msg_transposed(src, dst2, dst_cols, dst_rows)
        # assert dst == dst2
        # deinterleave(list(dst), dst)
    except:
        import traceback
        traceback.print_exc()
        print("src_size:", src_size)
        raise

    for i in range(len(dst) // dst_cols):
        print(dst[i * dst_cols:(i + 1) * dst_cols])

    print()
