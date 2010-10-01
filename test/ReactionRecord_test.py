#!/usr/bin/python

import _gfrd
import unittest

class ReactionRecordTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_instantiation(self):
        _gfrd.ReactionRecord()
        try:
            _gfrd.ReactionRecord(1, 1, 1)
            self.fail()
        except:
            pass
        try:
            _gfrd.ReactionRecord(1, 1, 1, 1)
            self.fail()
        except:
            pass
        try:
            _gfrd.ReactionRecord(1, (1,), 1)
            self.fail()
        except:
            pass
        try:
            _gfrd.ReactionRecord(1, (1,), 1, 1)
            self.fail()
        except:
            pass
        r = _gfrd.ReactionRecord(1, (_gfrd.ParticleID(0, 0),), _gfrd.ParticleID(0, 1))
        self.assertEqual(r.reaction_rule_id, 1)
        self.assertEqual(r.products, (_gfrd.ParticleID(0, 0),))
        self.assertEqual(r.reactants, (_gfrd.ParticleID(0, 1),))

        r = _gfrd.ReactionRecord(1, (_gfrd.ParticleID(0, 0),), _gfrd.ParticleID(0, 1), _gfrd.ParticleID(0, 2))
        self.assertEqual(r.reaction_rule_id, 1)
        self.assertEqual(r.products, (_gfrd.ParticleID(0, 0),))
        self.assertEqual(r.reactants, (_gfrd.ParticleID(0, 1),
                                       _gfrd.ParticleID(0, 2)))

        r = _gfrd.ReactionRecord(1, (_gfrd.ParticleID(0, 0), _gfrd.ParticleID(0, 1)), _gfrd.ParticleID(0, 1), _gfrd.ParticleID(0, 2))
        self.assertEqual(r.reaction_rule_id, 1)
        self.assertEqual(r.products, (_gfrd.ParticleID(0, 0),
                                      _gfrd.ParticleID(0, 1)))
        self.assertEqual(r.reactants, (_gfrd.ParticleID(0, 1),
                                       _gfrd.ParticleID(0, 2)))

if __name__ == "__main__":
    unittest.main()
