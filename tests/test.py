from selectiontest import *
from iotest import *
from analysistest import *
from frametest import *
import unittest


if __name__ == '__main__':
    test_io = unittest.TestLoader().loadTestsFromTestCase(TestIOFunctions)
    unittest.TextTestRunner(verbosity=2).run(test_io)

    test_selection = unittest.TestLoader().loadTestsFromTestCase(SelectionTest)
    unittest.TextTestRunner(verbosity=2).run(test_selection)

    test_analysis = unittest.TestLoader().loadTestsFromTestCase(TestFrameAnalysis)
    unittest.TextTestRunner(verbosity=2).run(test_analysis)

    test_frame = unittest.TestLoader().loadTestsFromTestCase(TestFullAnalysis)
    unittest.TextTestRunner(verbosity=2).run(test_frame)
