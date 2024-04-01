from Bio.Seq import Seq

CLIP = 2**128 - 1


class RollingHash:
    def __init__(self, k: int = 501, hbase: int = 239):
        self.k = k
        self.hbase = hbase
        self.d = {"A": 0, "C": 1, "G": 2, "T": 3}

    def hash_pair(self, sequence: str, pos: int) -> tuple[int, int]:
        h = 0
        for i in range(pos, pos + self.k):
            h = h * self.hbase + self.d[sequence[i]]
            h = h & CLIP

        h_rc = 0
        rc_sequence = str(Seq(sequence).reverse_complement())
        L = len(rc_sequence)
        for i in range(L - pos - self.k, L):
            h_rc = h_rc * self.hbase + self.d[rc_sequence[i]]
            h_rc = h_rc & CLIP

        return h, h_rc

    def hash(self, sequence: str, pos: int = 0) -> str:
        h, h_rc = self.hash_pair(sequence, pos)
        if h <= h_rc:
            return str(h) + "1" + sequence[pos + self.k : pos + self.k + 1]

        return str(h_rc) + "0" + sequence[pos + self.k : pos + self.k + 1]
