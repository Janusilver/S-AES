import numpy as np


class SAES:
    """
    S-AES (Simplified AES) 加密算法实现
    支持16位数据加密/解密，密钥扩展，以及多种工作模式
    """

    # S-Box 和逆 S-Box
    S_BOX = [
        [0x9, 0x4, 0xA, 0xB],
        [0xD, 0x1, 0x8, 0x5],
        [0x6, 0x2, 0x0, 0x3],
        [0xC, 0xE, 0xF, 0x7]
    ]

    INV_S_BOX = [
        [0xA, 0x5, 0x9, 0xB],
        [0x1, 0x7, 0x8, 0xF],
        [0x6, 0x0, 0x2, 0x3],
        [0xC, 0x4, 0xD, 0xE]
    ]

    # 轮常数
    RCON = [0x80, 0x30]

    def __init__(self):
        """初始化 S-AES 算法"""
        pass

    def gf_mult(self, a, b):
        """
        在 GF(2^4) 上的乘法运算，模数为 x^4 + x + 1

        Args:
            a, b: 要相乘的两个4位值

        Returns:
            GF(2^4)上的乘积结果
        """
        # 乘法表 (预计算好的)
        mult_table = [
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0xA, 0xB, 0xC, 0xD, 0xE, 0xF],
            [0, 2, 4, 6, 8, 0xA, 0xC, 0xE, 3, 1, 7, 5, 0xB, 9, 0xF, 0xD],
            [0, 3, 6, 5, 0xC, 0xF, 0xA, 9, 0xB, 8, 0xD, 0xE, 7, 4, 1, 2],
            [0, 4, 8, 0xC, 3, 7, 0xB, 0xF, 6, 2, 0xE, 0xA, 5, 1, 0xD, 9],
            [0, 5, 0xA, 0xF, 7, 2, 0xD, 8, 0xE, 0xB, 4, 1, 9, 0xC, 3, 6],
            [0, 6, 0xC, 0xA, 0xB, 0xD, 7, 1, 5, 3, 9, 0xF, 0xE, 8, 2, 4],
            [0, 7, 0xE, 9, 0xF, 8, 1, 6, 0xD, 0xA, 3, 4, 2, 5, 0xC, 0xB],
            [0, 8, 3, 0xB, 6, 0xE, 5, 0xD, 0xC, 4, 0xF, 7, 0xA, 2, 9, 1],
            [0, 9, 1, 8, 2, 0xB, 3, 0xA, 4, 0xD, 5, 0xC, 6, 0xF, 7, 0xE],
            [0, 0xA, 7, 0xD, 0xE, 4, 9, 3, 0xF, 5, 8, 2, 1, 0xB, 6, 0xC],
            [0, 0xB, 5, 0xE, 0xA, 1, 0xF, 4, 7, 0xC, 2, 9, 0xD, 6, 8, 3],
            [0, 0xC, 0xB, 7, 5, 9, 0xE, 2, 0xA, 6, 1, 0xD, 0xF, 3, 4, 8],
            [0, 0xD, 9, 4, 1, 0xC, 8, 5, 2, 0xF, 0xB, 6, 3, 0xE, 0xA, 7],
            [0, 0xE, 0xF, 1, 0xD, 3, 2, 0xC, 9, 7, 6, 8, 4, 0xA, 0xB, 5],
            [0, 0xF, 0xD, 2, 9, 6, 4, 0xB, 1, 0xE, 0xC, 3, 8, 7, 5, 0xA]
        ]
        return mult_table[a][b]

    def key_expansion(self, key):
        """
        密钥扩展函数，将16位密钥扩展为48位（3个轮密钥）

        Args:
            key: 16位初始密钥

        Returns:
            扩展后的密钥列表 [K0, K1, K2]
        """
        # 将16位密钥分成两个8位字
        w0 = (key >> 8) & 0xFF
        w1 = key & 0xFF

        # 计算 w2 = w0 ⊕ g(w1)
        g_w1 = self.g_function(w1, 1)
        w2 = w0 ^ g_w1

        # 计算 w3 = w2 ⊕ w1
        w3 = w2 ^ w1

        # 计算 w4 = w2 ⊕ g(w3)
        g_w3 = self.g_function(w3, 2)
        w4 = w2 ^ g_w3

        # 计算 w5 = w4 ⊕ w3
        w5 = w4 ^ w3

        # 构建轮密钥
        k0 = (w0 << 8) | w1
        k1 = (w2 << 8) | w3
        k2 = (w4 << 8) | w5

        return [k0, k1, k2]

    def g_function(self, word, round_num):
        """
        G函数用于密钥扩展

        Args:
            word: 8位输入字
            round_num: 轮数（1或2）

        Returns:
            变换后的8位字
        """
        # 半字节循环移位
        rotated = ((word & 0x0F) << 4) | ((word & 0xF0) >> 4)

        # 半字节代替
        nibble1 = (rotated >> 4) & 0x0F
        nibble2 = rotated & 0x0F

        sub_nibble1 = self.sub_nibble(nibble1, self.S_BOX)
        sub_nibble2 = self.sub_nibble(nibble2, self.S_BOX)

        substituted = (sub_nibble1 << 4) | sub_nibble2

        # 与轮常数异或
        result = substituted ^ self.RCON[round_num - 1]

        return result

    def sub_nibble(self, nibble, s_box):
        """
        半字节代替

        Args:
            nibble: 4位输入
            s_box: 使用的S盒

        Returns:
            代替后的4位输出
        """
        row = (nibble >> 2) & 0x03
        col = nibble & 0x03
        return s_box[row][col]

    def nibble_substitution(self, state, s_box):
        """
        对状态矩阵进行半字节代替

        Args:
            state: 2x2状态矩阵
            s_box: 使用的S盒

        Returns:
            代替后的状态矩阵
        """
        result = [[0, 0], [0, 0]]
        for i in range(2):
            for j in range(2):
                result[i][j] = self.sub_nibble(state[i][j], s_box)
        return result

    def shift_rows(self, state):
        """
        行移位变换

        Args:
            state: 2x2状态矩阵

        Returns:
            行移位后的状态矩阵
        """
        # 第二行循环左移一个半字节
        return [
            [state[0][0], state[0][1]],
            [state[1][1], state[1][0]]
        ]

    def mix_columns(self, state):
        """
        列混淆变换

        Args:
            state: 2x2状态矩阵

        Returns:
            列混淆后的状态矩阵
        """
        result = [[0, 0], [0, 0]]

        # 使用固定矩阵进行列混淆
        for j in range(2):  # 对每列操作
            result[0][j] = state[0][j] ^ self.gf_mult(4, state[1][j])
            result[1][j] = self.gf_mult(4, state[0][j]) ^ state[1][j]

        return result

    def inv_mix_columns(self, state):
        """
        逆列混淆变换

        Args:
            state: 2x2状态矩阵

        Returns:
            逆列混淆后的状态矩阵
        """
        result = [[0, 0], [0, 0]]

        # 使用逆矩阵进行列混淆
        for j in range(2):  # 对每列操作
            result[0][j] = self.gf_mult(9, state[0][j]) ^ self.gf_mult(2, state[1][j])
            result[1][j] = self.gf_mult(2, state[0][j]) ^ self.gf_mult(9, state[1][j])

        return result

    def add_round_key(self, state, round_key):
        """
        轮密钥加变换

        Args:
            state: 2x2状态矩阵
            round_key: 16位轮密钥

        Returns:
            轮密钥加后的状态矩阵
        """
        # 将轮密钥转换为2x2矩阵
        key_matrix = [
            [(round_key >> 12) & 0x0F, (round_key >> 4) & 0x0F],
            [(round_key >> 8) & 0x0F, round_key & 0x0F]
        ]

        result = [[0, 0], [0, 0]]
        for i in range(2):
            for j in range(2):
                result[i][j] = state[i][j] ^ key_matrix[i][j]

        return result

    def state_to_int(self, state):
        """
        将状态矩阵转换为16位整数

        Args:
            state: 2x2状态矩阵

        Returns:
            16位整数表示的状态
        """
        return (state[0][0] << 12) | (state[1][0] << 8) | (state[0][1] << 4) | state[1][1]

    def int_to_state(self, value):
        """
        将16位整数转换为状态矩阵

        Args:
            value: 16位整数

        Returns:
            2x2状态矩阵
        """
        return [
            [(value >> 12) & 0x0F, (value >> 4) & 0x0F],
            [(value >> 8) & 0x0F, value & 0x0F]
        ]

    def encrypt(self, plaintext, key):
        """
        S-AES加密函数

        Args:
            plaintext: 16位明文
            key: 16位密钥

        Returns:
            16位密文
        """
        # 密钥扩展
        round_keys = self.key_expansion(key)

        # 初始状态
        state = self.int_to_state(plaintext)

        # 第0轮：轮密钥加
        state = self.add_round_key(state, round_keys[0])

        # 第1轮：完整轮
        state = self.nibble_substitution(state, self.S_BOX)
        state = self.shift_rows(state)
        state = self.mix_columns(state)
        state = self.add_round_key(state, round_keys[1])

        # 第2轮：简化轮
        state = self.nibble_substitution(state, self.S_BOX)
        state = self.shift_rows(state)
        state = self.add_round_key(state, round_keys[2])

        return self.state_to_int(state)

    def decrypt(self, ciphertext, key):
        """
        S-AES解密函数

        Args:
            ciphertext: 16位密文
            key: 16位密钥

        Returns:
            16位明文
        """
        # 密钥扩展
        round_keys = self.key_expansion(key)

        # 初始状态
        state = self.int_to_state(ciphertext)

        # 第2轮逆操作
        state = self.add_round_key(state, round_keys[2])
        state = self.shift_rows(state)  # 逆行移位与行移位相同
        state = self.nibble_substitution(state, self.INV_S_BOX)

        # 第1轮逆操作
        state = self.add_round_key(state, round_keys[1])
        state = self.inv_mix_columns(state)
        state = self.shift_rows(state)
        state = self.nibble_substitution(state, self.INV_S_BOX)

        # 第0轮逆操作
        state = self.add_round_key(state, round_keys[0])

        return self.state_to_int(state)


class SAESExtended(SAES):
    """
    S-AES扩展功能类，支持多重加密和工作模式
    """

    def __init__(self):
        """初始化扩展功能"""
        super().__init__()

    def ascii_encrypt(self, text, key):
        """
        ASCII字符串加密

        Args:
            text: 要加密的ASCII字符串
            key: 16位密钥

        Returns:
            加密后的字节串
        """
        # 确保文本长度为偶数
        if len(text) % 2 != 0:
            text += ' '  # 填充空格

        result = b''
        for i in range(0, len(text), 2):
            # 将两个字符组合成16位块
            block = (ord(text[i]) << 8) | ord(text[i + 1])
            # 加密
            encrypted = self.encrypt(block, key)
            # 将加密结果拆分为两个字节
            result += bytes([(encrypted >> 8) & 0xFF, encrypted & 0xFF])

        return result

    def ascii_decrypt(self, data, key):
        """
        ASCII字符串解密

        Args:
            data: 要解密的字节串
            key: 16位密钥

        Returns:
            解密后的ASCII字符串
        """
        result = ""
        for i in range(0, len(data), 2):
            # 将两个字节组合成16位块
            block = (data[i] << 8) | data[i + 1]
            # 解密
            decrypted = self.decrypt(block, key)
            # 将解密结果拆分为两个字符
            result += chr((decrypted >> 8) & 0xFF) + chr(decrypted & 0xFF)

        return result

    def double_encrypt(self, plaintext, key):
        """
        双重加密

        Args:
            plaintext: 16位明文
            key: 32位密钥 (K1 || K2)

        Returns:
            16位密文
        """
        k1 = (key >> 16) & 0xFFFF
        k2 = key & 0xFFFF

        # 第一次加密
        intermediate = self.encrypt(plaintext, k1)
        # 第二次加密
        ciphertext = self.encrypt(intermediate, k2)

        return ciphertext

    def double_decrypt(self, ciphertext, key):
        """
        双重解密

        Args:
            ciphertext: 16位密文
            key: 32位密钥 (K1 || K2)

        Returns:
            16位明文
        """
        k1 = (key >> 16) & 0xFFFF
        k2 = key & 0xFFFF

        # 第一次解密
        intermediate = self.decrypt(ciphertext, k2)
        # 第二次解密
        plaintext = self.decrypt(intermediate, k1)

        return plaintext

    def meet_in_the_middle_attack(self, plaintext_ciphertext_pairs):
        """
        中间相遇攻击

        Args:
            plaintext_ciphertext_pairs: 明密文对列表 [(plaintext, ciphertext), ...]

        Returns:
            可能的密钥对列表 [(k1, k2), ...]
        """
        if not plaintext_ciphertext_pairs:
            return []

        # 使用第一对明密文进行攻击
        plaintext, ciphertext = plaintext_ciphertext_pairs[0]

        # 存储所有可能的中间值
        forward_table = {}

        # 遍历所有可能的K1
        for k1 in range(65536):
            intermediate = self.encrypt(plaintext, k1)
            forward_table[intermediate] = k1

        # 查找匹配的K2
        possible_keys = []
        for k2 in range(65536):
            intermediate = self.decrypt(ciphertext, k2)
            if intermediate in forward_table:
                k1 = forward_table[intermediate]
                # 验证其他明密文对
                valid = True
                for pt, ct in plaintext_ciphertext_pairs[1:]:
                    if self.double_encrypt(pt, (k1 << 16) | k2) != ct:
                        valid = False
                        break
                if valid:
                    possible_keys.append((k1, k2))

        return possible_keys

    def triple_encrypt(self, plaintext, key, mode=1):
        """
        三重加密

        Args:
            plaintext: 16位明文
            key: 密钥 (32位或48位)
            mode: 模式 (1: 32位密钥, 2: 48位密钥)

        Returns:
            16位密文
        """
        if mode == 1:
            # 32位密钥模式: E-D-E
            k1 = (key >> 16) & 0xFFFF
            k2 = key & 0xFFFF

            # 加密-解密-加密
            intermediate1 = self.encrypt(plaintext, k1)
            intermediate2 = self.decrypt(intermediate1, k2)
            ciphertext = self.encrypt(intermediate2, k1)
        else:
            # 48位密钥模式: E-E-E
            k1 = (key >> 32) & 0xFFFF
            k2 = (key >> 16) & 0xFFFF
            k3 = key & 0xFFFF

            # 加密-加密-加密
            intermediate1 = self.encrypt(plaintext, k1)
            intermediate2 = self.encrypt(intermediate1, k2)
            ciphertext = self.encrypt(intermediate2, k3)

        return ciphertext

    def triple_decrypt(self, ciphertext, key, mode=1):
        """
        三重解密

        Args:
            ciphertext: 16位密文
            key: 密钥 (32位或48位)
            mode: 模式 (1: 32位密钥, 2: 48位密钥)

        Returns:
            16位明文
        """
        if mode == 1:
            # 32位密钥模式: D-E-D
            k1 = (key >> 16) & 0xFFFF
            k2 = key & 0xFFFF

            # 解密-加密-解密
            intermediate1 = self.decrypt(ciphertext, k1)
            intermediate2 = self.encrypt(intermediate1, k2)
            plaintext = self.decrypt(intermediate2, k1)
        else:
            # 48位密钥模式: D-D-D
            k1 = (key >> 32) & 0xFFFF
            k2 = (key >> 16) & 0xFFFF
            k3 = key & 0xFFFF

            # 解密-解密-解密
            intermediate1 = self.decrypt(ciphertext, k3)
            intermediate2 = self.decrypt(intermediate1, k2)
            plaintext = self.decrypt(intermediate2, k1)

        return plaintext

    def cbc_encrypt(self, plaintext_blocks, key, iv):
        """
        CBC模式加密

        Args:
            plaintext_blocks: 明文块列表
            key: 16位密钥
            iv: 16位初始向量

        Returns:
            密文块列表
        """
        ciphertext_blocks = []
        previous_block = iv

        for block in plaintext_blocks:
            # 与前一密文块异或
            xored = block ^ previous_block
            # 加密
            encrypted = self.encrypt(xored, key)
            ciphertext_blocks.append(encrypted)
            previous_block = encrypted

        return ciphertext_blocks

    def cbc_decrypt(self, ciphertext_blocks, key, iv):
        """
        CBC模式解密

        Args:
            ciphertext_blocks: 密文块列表
            key: 16位密钥
            iv: 16位初始向量

        Returns:
            明文块列表
        """
        plaintext_blocks = []
        previous_block = iv

        for block in ciphertext_blocks:
            # 解密
            decrypted = self.decrypt(block, key)
            # 与前一密文块异或
            xored = decrypted ^ previous_block
            plaintext_blocks.append(xored)
            previous_block = block

        return plaintext_blocks


def test_basic_encryption():
    """测试基本加密功能"""
    print("=== 第1关：基本测试 ===")
    saes = SAES()

    # 测试用例
    plaintext = 0x1234
    key = 0x5678

    ciphertext = saes.encrypt(plaintext, key)
    decrypted = saes.decrypt(ciphertext, key)

    print(f"明文: 0x{plaintext:04X}")
    print(f"密钥: 0x{key:04X}")
    print(f"密文: 0x{ciphertext:04X}")
    print(f"解密: 0x{decrypted:04X}")
    print(f"加解密成功: {plaintext == decrypted}")
    print()


def test_cross_platform():
    """测试交叉平台兼容性"""
    print("=== 第2关：交叉测试 ===")
    saes = SAES()

    # 使用相同的测试向量
    test_vectors = [
        (0x1234, 0x5678),
        (0xABCD, 0xEF01),
        (0x0000, 0xFFFF),
        (0xFFFF, 0x0000)
    ]

    for plaintext, key in test_vectors:
        ciphertext = saes.encrypt(plaintext, key)
        decrypted = saes.decrypt(ciphertext, key)
        success = plaintext == decrypted

        print(f"明文: 0x{plaintext:04X}, 密钥: 0x{key:04X}, "
              f"密文: 0x{ciphertext:04X}, 解密: 0x{decrypted:04X}, "
              f"成功: {success}")

    print()


def test_ascii_encryption():
    """测试ASCII加密功能"""
    print("=== 第3关：扩展功能 - ASCII加密 ===")
    saes_ext = SAESExtended()

    text = "Hello S-AES!"
    key = 0x2D55  # 使用文档中的示例密钥

    print(f"原始文本: '{text}'")
    print(f"密钥: 0x{key:04X}")

    # 加密
    encrypted_data = saes_ext.ascii_encrypt(text, key)
    print(f"加密后 (十六进制): {encrypted_data.hex()}")

    # 解密
    decrypted_text = saes_ext.ascii_decrypt(encrypted_data, key)
    print(f"解密后: '{decrypted_text}'")

    print()


def test_multiple_encryption():
    """测试多重加密功能"""
    print("=== 第4关：多重加密 ===")
    saes_ext = SAESExtended()

    # 双重加密测试
    print("--- 双重加密 ---")
    plaintext = 0x1234
    double_key = 0x2D55ABCD  # 32位密钥

    double_encrypted = saes_ext.double_encrypt(plaintext, double_key)
    double_decrypted = saes_ext.double_decrypt(double_encrypted, double_key)

    print(f"明文: 0x{plaintext:04X}")
    print(f"双重加密密钥: 0x{double_key:08X}")
    print(f"双重加密密文: 0x{double_encrypted:04X}")
    print(f"双重解密结果: 0x{double_decrypted:04X}")
    print(f"双重加解密成功: {plaintext == double_decrypted}")

    # 三重加密测试
    print("\n--- 三重加密 (32位密钥模式) ---")
    triple_key_32 = 0x2D55ABCD  # 32位密钥
    triple_encrypted_32 = saes_ext.triple_encrypt(plaintext, triple_key_32, mode=1)
    triple_decrypted_32 = saes_ext.triple_decrypt(triple_encrypted_32, triple_key_32, mode=1)

    print(f"三重加密密文 (32位): 0x{triple_encrypted_32:04X}")
    print(f"三重解密结果 (32位): 0x{triple_decrypted_32:04X}")
    print(f"三重加解密成功: {plaintext == triple_decrypted_32}")

    print("\n--- 三重加密 (48位密钥模式) ---")
    triple_key_48 = 0x2D55ABCD1234  # 48位密钥
    triple_encrypted_48 = saes_ext.triple_encrypt(plaintext, triple_key_48, mode=2)
    triple_decrypted_48 = saes_ext.triple_decrypt(triple_encrypted_48, triple_key_48, mode=2)

    print(f"三重加密密文 (48位): 0x{triple_encrypted_48:04X}")
    print(f"三重解密结果 (48位): 0x{triple_decrypted_48:04X}")
    print(f"三重加解密成功: {plaintext == triple_decrypted_48}")

    # 中间相遇攻击演示
    print("\n--- 中间相遇攻击演示 ---")
    # 生成一个测试密钥对
    test_k1, test_k2 = 0x2D55, 0xABCD
    test_double_key = (test_k1 << 16) | test_k2

    # 创建几个明密文对
    pairs = []
    for i in range(2):
        pt = (0x1234 + i * 0x111) & 0xFFFF
        ct = saes_ext.double_encrypt(pt, test_double_key)
        pairs.append((pt, ct))

    print(f"真实密钥: K1=0x{test_k1:04X}, K2=0x{test_k2:04X}")

    # 执行中间相遇攻击
    possible_keys = saes_ext.meet_in_the_middle_attack(pairs)

    print(f"找到的可能密钥: {len(possible_keys)} 个")
    for i, (k1, k2) in enumerate(possible_keys[:3]):  # 只显示前3个
        print(f"  可能密钥 {i + 1}: K1=0x{k1:04X}, K2=0x{k2:04X}")

    print()


def test_cbc_mode():
    """测试CBC工作模式"""
    print("=== 第5关：工作模式 - CBC ===")
    saes_ext = SAESExtended()

    # 生成测试数据
    plaintext_blocks = [0x1234, 0x5678, 0x9ABC, 0xDEF0]
    key = 0x2D55
    iv = 0xFFFF  # 初始向量

    print("原始明文块:", [f"0x{block:04X}" for block in plaintext_blocks])
    print(f"密钥: 0x{key:04X}")
    print(f"初始向量: 0x{iv:04X}")

    # CBC加密
    ciphertext_blocks = saes_ext.cbc_encrypt(plaintext_blocks, key, iv)
    print("CBC加密后:", [f"0x{block:04X}" for block in ciphertext_blocks])

    # CBC解密
    decrypted_blocks = saes_ext.cbc_decrypt(ciphertext_blocks, key, iv)
    print("CBC解密后:", [f"0x{block:04X}" for block in decrypted_blocks])

    # 验证加解密成功
    success = plaintext_blocks == decrypted_blocks
    print(f"CBC加解密成功: {success}")

    # 篡改密文测试
    print("\n--- 密文篡改测试 ---")
    if len(ciphertext_blocks) > 1:
        # 篡改第二个密文块
        tampered_ciphertext = ciphertext_blocks.copy()
        tampered_ciphertext[1] ^= 0x0F0F  # 修改部分比特

        # 尝试解密被篡改的密文
        tampered_decrypted = saes_ext.cbc_decrypt(tampered_ciphertext, key, iv)
        print("篡改后解密:", [f"0x{block:04X}" for block in tampered_decrypted])

        # 显示错误传播
        errors = []
        for i, (orig, tampered) in enumerate(zip(plaintext_blocks, tampered_decrypted)):
            if orig != tampered:
                errors.append(i)

        print(f"错误传播到块: {errors}")

    print()


def main():
    """主函数，运行所有测试"""
    print("S-AES 加密系统测试")
    print("=" * 50)

    # 运行所有测试
    test_basic_encryption()
    test_cross_platform()
    test_ascii_encryption()
    test_multiple_encryption()
    test_cbc_mode()

    print("所有测试完成！")


if __name__ == "__main__":
    main()