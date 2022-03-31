import unittest
from metadensity.metadensity import *

class Test(unittest.TestCase):

    def test_order_multifeature(self):

        exons = defaultdict(set, {'exon': set([
            (100,110), (160,190), (130,150)
        ])
        
        }
        )

        # PLUS STRAND
        m = Metagene(esnt='FAKE_GENE', chro='chr1', start=100, end=200, strand='+', transcript_type='protein_coding'
        , features = exons)

        sorted_exon = m.order_multi_feature(feature='exon')

        self.assertEqual(sorted_exon, [(100,110), (130,150), (160,190)])

        # MINUS STRAND
        m = Metagene(esnt='FAKE_GENE', chro='chr1', start=100, end=200, strand='-', transcript_type='protein_coding'
        , features = exons)

        sorted_exon = m.order_multi_feature(feature='exon')

        self.assertEqual(sorted_exon, [(100,110), (130,150), (160,190)][::-1])
    
    def test_create_feature(self):

        exons = defaultdict(set, {'exon': set([
            (100,110), (160,190), (130,150)
        ])})

        # PLUS STRAND
        m = Metagene(esnt='FAKE_GENE', chro='chr1', start=100, end=200, strand='+', transcript_type='protein_coding'
        , features = exons)

        # test point feature
        m.create_feature((150,155), 'random_window')

        self.assertIn('random_window', m.featnames)
        self.assertEqual(m.features['random_window'], set([(150,155)]))
    
    def test_create_downstream(self):
        
        exons = defaultdict(set, {'exon': set([
            (100,110), (160,190), (130,150)
        ])})

        # PLUS STRAND
        m = Metagene(esnt='FAKE_GENE', chro='chr1', start=100, end=200, strand='+', transcript_type='protein_coding'
        , features = exons)

        m.create_downstream_feature((150,155), feature_type = 'exon', feature_name = 'random')

        self.assertIn('random_downstream_exon', m.featnames)
        self.assertEqual(m.features['random_downstream_exon'], set([(160,190)]))

        # when no downstream exsit
        m.create_downstream_feature((190,200), feature_type = 'exon', feature_name = 'random2')
        self.assertEqual(m.features['random2_downstream_exon'], set())

        # MINUS STRAND
        m = Metagene(esnt='FAKE_GENE', chro='chr1', start=100, end=200, strand='-', transcript_type='protein_coding'
        , features = exons)

        m.create_downstream_feature((150,155), feature_type = 'exon', feature_name = 'random3')
        self.assertEqual(m.features['random3_downstream_exon'], set([(130,150)]))

        m.create_downstream_feature((100,105), feature_type = 'exon', feature_name = 'random4')
        self.assertEqual(m.features['random4_downstream_exon'], set())
    

    def test_create_downstream(self):
        
        exons = defaultdict(set, {'exon': set([
            (100,110), (160,190), (130,150)
        ])})

        # PLUS STRAND
        m = Metagene(esnt='FAKE_GENE', chro='chr1', start=100, end=200, strand='+', transcript_type='protein_coding'
        , features = exons)

        m.create_upstream_feature((150,155), feature_type = 'exon', feature_name = 'random')

        self.assertIn('random_upstream_exon', m.featnames)
        self.assertEqual(m.features['random_upstream_exon'], set([(130,150)]))

        # when no downstream exsit
        m.create_upstream_feature((100,105), feature_type = 'exon', feature_name = 'random2')
        self.assertEqual(m.features['random2_upstream_exon'], set())

        # MINUS STRAND
        m = Metagene(esnt='FAKE_GENE', chro='chr1', start=100, end=200, strand='-', transcript_type='protein_coding'
        , features = exons)

        m.create_upstream_feature((150,155), feature_type = 'exon', feature_name = 'random3')
        self.assertEqual(m.features['random3_upstream_exon'], set([(160,190)]))

        m.create_upstream_feature((190,195), feature_type = 'exon', feature_name = 'random4')
        self.assertEqual(m.features['random4_upstream_exon'], set())


    def test_relative_feature(self):

        exons = defaultdict(set, {'exon': set([
            (100,110), (160,190), (130,150)
        ])})

        # PLUS STRAND
        m = Metagene(esnt='FAKE_GENE', chro='chr1', start=100, end=200, strand='+', transcript_type='protein_coding'
        , features = exons)

        rela_pos = m.to_relative_axis((100,150))
        self.assertEqual(rela_pos[0], 0)
        self.assertEqual(rela_pos[1], 50)

        # MINUS STRAND
        m = Metagene(esnt='FAKE_GENE', chro='chr1', start=100, end=200, strand='-', transcript_type='protein_coding'
        , features = exons)

        rela_pos = m.to_relative_axis((100,150))
        self.assertEqual(rela_pos[0], 50)
        self.assertEqual(rela_pos[1], 100)

if __name__ == '__main__':
    unittest.main()

        