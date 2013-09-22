import unittest

from decorator2 import create_species


class TestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test_species(self):
        sp1 = create_species("A")

        self.assertTrue(str(sp1) == "A()")
        self.assertTrue(sp1.match(sp1))


if __name__ == '__main__':
    unittest.main()
