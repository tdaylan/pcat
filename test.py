class cntrstrt():
    def incr(self):
        self.cntr += 1
        return self.cntr
    def __init__(self):
        self.cntr = 0
cntr = cntrstrt()
