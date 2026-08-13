#!/usr/bin/env python3
"""SAND の合成 S-box SSb = G0 || G1。

G0, G1 は 4 ビット写像で SAND-64 と SAND-128 で完全に同一なので,
結合 DDT とその CNF 符号化の厳密性検査も両変種で共通である。
"""


class SSbMixin:
    """SSb 部分。変種固有のモデルクラスに mixin する。"""

    def ssb_s_box(self, a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, w1, w2, w3):
        """
        Joint differential constraint for the synthetic S-box SSb = (G0, G1).

        Number of constraints: 58
        Input:      a0||a1||a2||a3; a0: msb   (input difference of the shared nibble)
        Output G0:  b0||b1||b2||b3; b0: msb
        Output G1:  c0||c1||c2||c3; c0: msb
        Weight:     w1 + w2 + w3, encoded in unary with w1 >= w2 >= w3.

        The nonzero entries of the joint DDT are 2, 4 and 16, so the weight
        -log2(Pr) of a valid transition is 3, 2 or 0; a weight of 1 never occurs
        and is excluded by the clauses below.  The clause set was generated from
        ssb_joint_ddt() and its models were exhaustively checked against the
        table by verify_ssb_model().

        SAND-64 and SAND-128 use the same 4-bit G0 and G1, so this clause set is
        shared verbatim; the variants differ only in the number of parallel
        S-boxes and the nibble rotation amount, which s_box_model() supplies.
        """
        return [[-a0,c0],[-a0,w2],[-a1,b1],[-a2,b2],[-a3,c3],[-b3,w2],[-c1,w1],[-c2,w1],[-w1,w2],[a0,-c0],[a1,-b1],[a2,-b2],[a3,-c3],[w1,-w2],[w1,-w3],[-c0,-c1,w3],[b2,c1,-w3],[-a0,b0,b2,b3],[-a2,-a3,c1,w3],[-a2,c1,c2,c3],[-b0,-c1,-c2,w3],[-b0,b2,b3,c0],[-b1,-b2,-c3,w3],[a1,b0,-c0,w3],[a2,-b1,c0,c1],[a3,-c0,-c2,w3],[b0,-b2,c2,w3],[b1,-b3,c0,c3],[b1,b3,-c1,w3],[b1,b3,c0,-c3],[b2,-b3,-c1,w3],[c0,c1,c3,-w3],[-a0,a2,b0,-b1,-c1],[-a1,-a3,-b3,c2,w3],[-a1,-b3,c0,c1,-w3],[-a1,b3,c0,-c1,-w3],[-a2,-b3,c0,-c1,c3],[-a3,b1,c0,-c1,-w3],[-b0,-b3,-c0,-c2,w3],[-b0,b1,b2,-b3,-c1],[-b1,b3,c1,-c2,w3],[a1,-a3,b3,c2,w3],[a1,b0,b2,b3,-w1],[a2,-b0,-b1,-b3,w3],[-a0,-b0,-b2,b3,c2,c3],[-a0,b0,b1,b3,-c1,-c2],[-a1,-a2,b3,-c0,c1,c2],[-a1,-b2,-b3,c1,-c2,-c3],[-a1,b0,-c0,-c1,-c2,c3],[-a3,b1,-b2,b3,c1,-c2],[-b0,b1,-b3,-c1,-c2,c3],[a1,-b3,-c0,c1,c2,-w3],[a1,a3,b0,-b2,-b3,c2],[b0,-b1,b3,-c0,c2,-c3],[-a0,-a1,-a2,-b0,-b3,-c1,c2],[-a0,b0,-b2,-b3,-c1,-c2,-c3],[-a1,-a2,-b0,b3,-c1,-c2,-c3],[-a2,-b0,b1,-c0,-c1,c2,-c3]]

    @staticmethod
    def g0_value(x):
        """G0 acting on a 4-bit value x (bit 3 = msb). Used only for verification."""
        x3, x2, x1, x0 = (x >> 3) & 1, (x >> 2) & 1, (x >> 1) & 1, x & 1
        y0 = (x3 & x2) ^ x0
        y1 = x1
        y2 = x2
        y3 = (y0 & x1) ^ x3
        return (y3 << 3) | (y2 << 2) | (y1 << 1) | y0

    @staticmethod
    def g1_value(x):
        """G1 acting on a 4-bit value x (bit 3 = msb). Used only for verification."""
        x3, x2, x1, x0 = (x >> 3) & 1, (x >> 2) & 1, (x >> 1) & 1, x & 1
        y3 = x3
        y2 = (x3 & x1) ^ x2
        y1 = (y2 & x0) ^ x1
        y0 = x0
        return (y3 << 3) | (y2 << 2) | (y1 << 1) | y0

    @classmethod
    def ssb_joint_ddt(cls):
        """
        Differential distribution table of the 4x8 synthetic S-box

            SSb(x) = G0(x) || G1(x),

        i.e. ddt[din][(dG0 << 4) | dG1] counts the 4-bit values x with
        G0(x) ^ G0(x ^ din) = dG0 and G1(x) ^ G1(x ^ din) = dG1 simultaneously.

        G0 and G1 consume the *same* 4-bit input value, so their output
        differences must not be drawn from two independent DDTs: the product of
        the marginal DDTs admits 157 (din, dG0, dG1) triples whereas only 97 are
        jointly realizable.
        """
        ddt = [[0] * 256 for _ in range(16)]
        for din in range(16):
            for x in range(16):
                d0 = cls.g0_value(x) ^ cls.g0_value(x ^ din)
                d1 = cls.g1_value(x) ^ cls.g1_value(x ^ din)
                ddt[din][(d0 << 4) | d1] += 1
        return ddt

    def verify_ssb_model(self, verbose=True):
        """
        Exhaustively check that the models of ssb_s_box() are exactly the valid
        (input difference, G0 output difference, G1 output difference, weight)
        tuples of the joint DDT.  All 2^15 assignments of the 15 involved
        variables are enumerated, so this is a proof, not a sample.
        """
        clauses = self.ssb_s_box(*range(1, 16))  # variables 1..15 in signature order
        ddt = self.ssb_joint_ddt()

        expected = set()
        for din in range(16):
            for dout, count in enumerate(ddt[din]):
                if count == 0:
                    continue
                weight = 4 - count.bit_length() + 1
                assert count == 1 << (count.bit_length() - 1), "non power-of-two DDT entry"
                bits = []
                for value in (din, (dout >> 4) & 0xF, dout & 0xF):
                    bits += [(value >> 3) & 1, (value >> 2) & 1, (value >> 1) & 1, value & 1]
                bits += [1 if k < weight else 0 for k in range(3)]
                expected.add(tuple(bits))

        models = set()
        for code in range(1 << 15):
            bits = [(code >> (14 - k)) & 1 for k in range(15)]
            assign = {k + 1: bits[k] for k in range(15)}
            if all(any(assign[abs(l)] == (1 if l > 0 else 0) for l in clause) for clause in clauses):
                models.add(tuple(bits))

        if models != expected:
            raise AssertionError(
                f"SSb CNF is not exact: {len(models - expected)} spurious and "
                f"{len(expected - models)} missing assignments"
            )
        if verbose:
            print(f"SSb joint-DDT model verified: {len(expected)} valid transitions, "
                  f"{len(clauses)} clauses, all 2^15 assignments checked")
        return True
