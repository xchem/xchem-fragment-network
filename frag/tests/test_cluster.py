import unittest
from frag.cluster import dpmeans,kmeans



class ParserTest(unittest.TestCase):

    def test_dp_cluster(self):
        """
        Basic test to show clustering works.
        :return:
        """
        dp = dpmeans([(1.0,1.0,1.0),(3.0,3.0,3.0)], 1.0, 0, False)
        dp.run()
        self.assertEqual(len(dp.clusters),2)

    def test_kmean_cluster(self):
        """
        Simple test to show clustering works.
        :return:
        """
        kmean = kmeans([(1.0,1.0,1.0),(3.0,3.0,3.0)], 2, 0, False)
        kmean.run()
        self.assertEqual(len(kmean.clusters),2)