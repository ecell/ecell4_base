from ecell4.core import *

import  unittest


class ContextText(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        sp1 = Species("A")
        self.assertEqual(sp1.count(Species("A")), 1)
        self.assertEqual(sp1.count(Species("A.A")), 2)

    def test2(self):
        sp1 = Species("A.B")
        self.assertEqual(sp1.count(Species("A.B")), 1)
        self.assertEqual(sp1.count(Species("B.A")), 1)

    def test3(self):
        sp1 = Species("A(p=u^_)")
        self.assertEqual(sp1.count(Species("A(p=u^1).B(b^1)")), 1)


if __name__ == "__main__":
    unittest.main()
