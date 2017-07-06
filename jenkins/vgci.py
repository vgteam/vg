#!/usr/bin/env python2.7
import logging
import shutil
import subprocess
import tempfile
import textwrap
import filecmp
import pytest
from unittest import TestCase, skip, TestLoader, TextTestRunner
from urlparse import urlparse
from uuid import uuid4
import os, sys
import argparse
import timeout_decorator
from boto.s3.connection import S3Connection

log = logging.getLogger(__name__)

class VGCITest(TestCase):
    """
    Continuous Integration VG tests.  All depend on toil-vg being installed.  Along with 
    toil[aws,mesos].  They are somewhat derived from the toil-vg unittests, but are
    much slower.  
    """
    def setUp(self):
        self.workdir = tempfile.mkdtemp()
        self.tempdir = tempfile.mkdtemp()
        
        self.f1_threshold = 0.005
        self.auc_threshold = 0.02
        self.input_store = 's3://cgl-pipeline-inputs/vg_cgl/bakeoff'
        self.vg_docker = None
        self.container = None # Use default in toil-vg, which is Docker
        self.verify = True
        self.do_teardown = True
        self.baseline = 's3://cgl-pipeline-inputs/vg_cgl/vg_ci/jenkins_regression_baseline'
        self.cores = 8

        self.loadCFG()
                
    def tearDown(self):
        shutil.rmtree(self.tempdir)        
        if self.do_teardown:
            shutil.rmtree(self.workdir)

    def loadCFG(self):
        """ It's a hassle passing parameters through pytest.  Hack
        around for now by loading from a file of key/value pairs. """
        if os.path.isfile('vgci_cfg.tsv'):
            with open('vgci_cfg.tsv') as f:
                for line in f:
                    toks = line.split()
                    if len(toks) == 2 and toks[0][0] != '#':
                        # override vg docker (which defaults to value from vg_config.py)
                        if toks[0] == 'vg-docker-version':
                            self.vg_docker = toks[1]
                        # can use "Docker", "Singularity" or "None" (the string) as a container system
                        if toks[0] == 'container':
                            self.container = toks[1]
                        # dont verify output.  tests will pass if they dont crash or timeout
                        elif toks[0] == 'verify' and toks[1].lower() == 'false':
                            self.verify = False
                        # dont delete the working directory
                        elif toks[0] == 'teardown' and toks[1].lower() == 'false':
                            self.do_teardown = False
                        # override the working directory (defaults to temp)
                        elif toks[0] == 'workdir':
                            self.workdir = toks[1]
                        elif toks[0] == 'baseline':
                            self.baseline = toks[1]
                        elif toks[0] == 'cores':
                            self.cores = int(toks[1])

    def _jobstore(self, tag = ''):
        return os.path.join(self.workdir, 'jobstore{}'.format(tag))

    def _outstore(self, tag = ''):
        return os.path.join(self.workdir, 'outstore-{}'.format(tag))

    def _input(self, filename):
        return os.path.join(self.input_store, filename)

    def _bakeoff_coords(self, region):
        if region == 'BRCA1':
            return 17, 43044293
        elif region == 'BRCA2':
            return 13, 32314860
        elif region == 'SMA':
            return 5, 69216818
        elif region == 'LRC-KIR':
            return 19, 54025633
        elif region == 'MHC':
            return 6, 28510119

    def _read_baseline_file(self, tag, path):
        """ read a (small) text file from the baseline store """
        if self.baseline.startswith('s3://'):
            toks = self.baseline[5:].split('/')
            bname = toks[0]
            keyname = '/{}/outstore-{}/{}'.format('/'.join(toks[1:]), tag, path)
            bucket = S3Connection().get_bucket(bname)
            key = bucket.get_key(keyname)
            return key.get_contents_as_string()
        else:
            with open(os.path.join(self.baseline, 'outstore-{}'.format(tag), path)) as f:
                return f.read()

    def _get_remote_file(self, src, tgt):
        """ get a file from a store """
        if not os.path.exists(os.path.dirname(tgt)):
            os.makedirs(os.path.dirname(tgt))
        if src.startswith('s3://'):
            toks = src[5:].split('/')
            bname = toks[0]
            keyname = '/' + '/'.join(toks[1:])
            bucket = S3Connection().get_bucket(bname)
            key = bucket.get_key(keyname)
            with open(tgt, 'w') as f:
                return key.get_contents_to_file(f)
        else:
            shutil.copy2(src, tgt_file)

    def _toil_vg_index(self, chrom, graph_path, xg_path, gcsa_path, misc_opts, dir_tag, file_tag):
        """ Wrap toil-vg index.  Files passed are copied from store instead of computed """
        job_store = self._jobstore(dir_tag)
        out_store = self._outstore(dir_tag)
        opts = '--realTimeLogging --logInfo --config jenkins/toil_vg_config.yaml '
        if self.vg_docker:
            opts += '--vg_docker {} '.format(self.vg_docker)
        if self.container:
            opts += '--container {} '.format(self.container)
        if chrom:
            opts += '--chroms {} '.format(chrom)
        if graph_path:
            opts += '--graphs {} '.format(graph_path)
        if xg_path:
            opts += '--skip_xg '
            self._get_remote_file(xg_path, os.path.join(out_store, os.path.basename(xg_path)))
        if gcsa_path:
            opts += '--skip_gcsa '
            self._get_remote_file(gcsa_path, os.path.join(out_store, os.path.basename(gcsa_path)))
            self._get_remote_file(gcsa_path + '.lcp', os.path.join(out_store, os.path.basename(gcsa_path) + '.lcp'))
        opts += '--index_name {}'.format(file_tag)
        if misc_opts:
            opts += misc_opts + ' '
        
        cmd = 'toil-vg index {} {} {}'.format(job_store, out_store, opts)
        
        subprocess.check_call(cmd, shell=True)        
        
        
    def _toil_vg_run(self, sample_name, chrom, graph_path, xg_path, gcsa_path, fq_path,
                     true_vcf_path, fasta_path, interleaved, misc_opts, tag):
        """ Wrap toil-vg run as a shell command.  Expects reads to be in single fastq
        inputs can be None if toil-vg supports not having them (ie don't need to 
        include gcsa_path if want to reindex)
        """

        job_store = self._jobstore(tag)
        out_store = self._outstore(tag)
        opts = '--realTimeLogging --logInfo --config jenkins/toil_vg_config.yaml '
        if self.vg_docker:
            opts += '--vg_docker {} '.format(self.vg_docker)
        if self.container:
            opts += '--container {} '.format(self.container)
        if chrom:
            opts += '--chroms {} '.format(chrom)
        if graph_path:
            opts += '--graphs {} '.format(graph_path)
        if xg_path:
            opts += '--xg_index {} '.format(xg_path)
        if gcsa_path:
            opts += '--gcsa_index {} '.format(gcsa_path)
        if fq_path:
            opts += '--fastq {} '.format(fq_path)
        if true_vcf_path:
            opts += '--vcfeval_baseline {} '.format(true_vcf_path)
            opts += '--vcfeval_fasta {} '.format(fasta_path)
            opts += '--vcfeval_opts \" --ref-overlap\" '
        if interleaved:
            opts += '--interleaved '
        if misc_opts:
            opts += misc_opts + ' '
        opts += '--gcsa_index_cores {} --kmers_cores {} \
        --alignment_cores {} --calling_cores {} --vcfeval_cores {} '.format(
            self.cores, self.cores, self.cores, self.cores, self.cores)
        
        cmd = 'toil-vg run {} {} {} {}'.format(job_store, sample_name, out_store, opts)
        
        subprocess.check_call(cmd, shell=True)        

    def _verify_f1(self, sample, tag='', threshold=None):
        # grab the f1.txt file from the output store
        if sample:
            f1_name = '{}_vcfeval_output_f1.txt'.format(sample)
        else:
            f1_name = 'vcfeval_output_f1.txt'
        f1_path = os.path.join(self._outstore(tag), f1_name)
        with open(f1_path) as f1_file:
            f1_score = float(f1_file.readline().strip())
        baseline_f1 = float(self._read_baseline_file(tag, f1_name).strip())
        
        # compare with threshold
        if not threshold:
            threshold = self.f1_threshold
            
        print 'F1: {}  Baseline: {}  Threshold: {}'.format(
            f1_score, baseline_f1, threshold)
        self.assertTrue(f1_score >= baseline_f1 - threshold)

    def _test_bakeoff(self, region, graph, skip_indexing):
        """ Run bakeoff F1 test for NA12878 """
        tag = '{}-{}'.format(region, graph)
        chrom, offset = self._bakeoff_coords(region)        
        if skip_indexing:
            xg_path = None
            gcsa_path = self._input('{}-{}.gcsa'.format(graph, region))
        else:
            xg_path = None
            gcsa_path = None            
        self._toil_vg_run('NA12878', chrom,
                          self._input('{}-{}.vg'.format(graph, region)),
                          xg_path, gcsa_path,
                          self._input('platinum_NA12878_{}.fq.gz'.format(region)),
                          self._input('platinum_NA12878_{}.vcf.gz'.format(region)),
                          self._input('chr{}.fa.gz'.format(chrom)), True,
                          '--vcf_offsets {}'.format(offset), tag)

        if self.verify:
            self._verify_f1('NA12878', tag)

    def _mapeval_vg_run(self, reads, base_xg_path, fasta_path, test_index_bases,
                        test_names, score_baseline_name, tag):
        """ Wrap toil-vg mapeval as a shell command. 
        
        Evaluates realignments (to the linear reference and to a set of graphs)
        of reads simulated from a single "base" graph. Realignments are
        evaluated based on how close the realignments are to the original
        simulated source position. Simulations are done inside this function.
        
        Simulates the given number of reads (reads), from the given XG file
        (base_xg_path). Uses the given FASTA (fasta_path) as a BWA reference for
        comparing vg and BWA alignments within mapeval. (Basically, BWA against
        the linear reference functions as a negative control "graph" to compare
        against the real test graphs.)
        
        test_index_bases specifies a list of basenames (without extension) for a
        .xg, .gcsa, and .gcsa.lcp file set, one per of graph that is to be
        compared.
        
        test_names has one entry per graph to be compared, and specifies where
        the realigned read GAM files should be saved.
        
        score_baseline_name, if not None, is a name from test_names to be used
        as a score baseline for comparing all the realignment scores against.
        
        tag is a unique slug for this test/run, which determines the Toil job
        store name to use, and the location where the output files should be
        saved.
        
        """

        job_store = self._jobstore(tag)
        out_store = self._outstore(tag)

        # start by simulating some reads
        # TODO: why are we using strings here when we could use much safer lists???
        opts = '--realTimeLogging --logInfo --config jenkins/toil_vg_config.yaml '
        if self.vg_docker:
            opts += '--vg_docker {} '.format(self.vg_docker)
        if self.container:
            opts += '--container {} '.format(self.container)
        # note, using the same seed only means something if using same
        # number of chunks.  we make that explicit here
        opts += '--maxCores {} --sim_chunks {} --seed {} '.format(self.cores, self.cores, self.cores)
        opts += '--sim_opts \'-l 150 -p 500 -v 50 -e 0.05 -i 0.01\' '
        cmd = 'toil-vg sim {} {} {} {} --gam {}'.format(
            job_store, base_xg_path, reads / 2, out_store, opts)
        subprocess.check_call(cmd, shell=True)

        # then run mapeval
        opts = '--realTimeLogging --logInfo '        
        if self.vg_docker:
            opts += '--vg_docker {} '.format(self.vg_docker)
        if self.container:
            opts += '--container {} '.format(self.container)
        opts += '--maxCores {} '.format(self.cores)
        opts += '--bwa --bwa-paired --vg-paired '
        opts += '--fasta {} '.format(fasta_path)
        opts += '--index-bases {} '.format(' '.join(test_index_bases))
        opts += '--gam-names {} '.format(' '.join(test_names))
        opts += '--gam_input_reads {} '.format(os.path.join(out_store, 'sim.gam'))
        opts += '--alignment_cores {} '.format(self.cores)
        if score_baseline_name is not None:
            opts += '--compare-gam-scores {} '.format(score_baseline_name)
        
        cmd = 'toil-vg mapeval {} {} {} {}'.format(
            job_store, out_store, os.path.join(out_store, 'true.pos'), opts)
        subprocess.check_call(cmd, shell=True)

    def _tsv_to_dict(self, stats, row_1 = 1):
        """ convert tsv string into dictionary """
        stats_dict = dict()
        for line in stats.split('\n')[row_1:]:
            toks = line.split()
            if len(toks) > 1:
                stats_dict[toks[0]] = [float(x) for x in toks[1:]]
        return stats_dict

    def _verify_mapeval(self, reads, tag):
        """ Check the simulated mapping evaluation results """

        stats_path = os.path.join(self._outstore(tag), 'stats.tsv')
        with open(stats_path) as stats:
            stats_tsv = stats.read()
        baseline_tsv = self._read_baseline_file(tag, 'stats.tsv')

        stats_dict = self._tsv_to_dict(stats_tsv)
        baseline_dict = self._tsv_to_dict(baseline_tsv)

        for key, val in baseline_dict.items():
            print '{}  Acc: {} Baseline: {}  Auc: {} Baseline: {}  Threshold: {}'.format(
                key, stats_dict[key][1], val[1], stats_dict[key][2], val[2], self.auc_threshold)
            self.assertTrue(stats_dict[key][0] == reads)
            self.assertTrue(stats_dict[key][1] >= val[1] - self.auc_threshold)
            # disable roc test for now
            #self.assertTrue(stats_dict[key][2] >= val[2] - self.auc_threshold)
        
    def _test_mapeval(self, reads, region, baseline_graph, test_graphs, score_baseline_graph=None):
        """ Run simulation on a bakeoff graph
        
        Simulate the given number of reads from the given baseline_graph
        (snp1kg, primary, etc.) and realign them against all the graphs in the
        test_graph list.
        
        Needs to know the bekeoff region that is being run, in order to look up
        the actual graphs files for each graph type.
        
        Verifies that the realignments are sufficiently good.
        
        If score_baseline is set to a graph name from test_graphs, computes
        score differences for reach read against that baseline.
        
        """
        tag = 'sim-{}-{}'.format(region, baseline_graph)
        
        # compute the xg indexes from scratch
        index_bases = []
        for graph in set([baseline_graph] + test_graphs):
            chrom, offset = self._bakeoff_coords(region)        
            vg_path = self._input('{}-{}.vg'.format(graph, region))
            self._toil_vg_index(chrom, vg_path, None, self._input('{}-{}.gcsa'.format(graph, region)),
                                None, tag, '{}-{}'.format(graph, region))
            
        fasta_path = self._input('{}.fa'.format(region))
        xg_path = os.path.join(self._outstore(tag), '{}-{}'.format(baseline_graph, region) + '.xg')
        test_index_bases = []
        for test_graph in test_graphs:
            test_tag = '{}-{}'.format(test_graph, region)
            test_index_bases.append(os.path.join(self._outstore(tag), test_tag))
        test_xg_paths = os.path.join(self._outstore(tag), tag + '.xg')
        self._mapeval_vg_run(reads, xg_path, fasta_path, test_index_bases,
                             test_graphs, score_baseline_graph, tag)

        if self.verify:
            self._verify_mapeval(reads, tag)

    @timeout_decorator.timeout(3600)
    def test_sim_brca2_snp1kg(self):
        """ Mapping and calling bakeoff F1 test for BRCA1 primary graph """
        # Using 50k simulated reads from snp1kg BRCA1, realign against all these
        # other BRCA1 graphs and make sure the realignments are sufficiently
        # good.
        # Compare all realignment scores agaisnt the scores for the primary
        # graph.
        self._test_mapeval(50000, 'BRCA1', 'snp1kg',
                           ['primary', 'snp1kg', 'cactus'],
                           score_baseline_graph='primary')

    @timeout_decorator.timeout(3600)
    def test_sim_mhc_snp1kg(self):
        """ Mapping and calling bakeoff F1 test for BRCA1 primary graph """        
        self._test_mapeval(50000, 'MHC', 'snp1kg',
                           ['primary', 'snp1kg', 'cactus'])    

    @timeout_decorator.timeout(200)
    def test_map_brca1_primary(self):
        """ Mapping and calling bakeoff F1 test for BRCA1 primary graph """
        self._test_bakeoff('BRCA1', 'primary', True)

    @timeout_decorator.timeout(200)        
    def test_map_brca1_snp1kg(self):
        """ Mapping and calling bakeoff F1 test for BRCA1 snp1kg graph """
        self._test_bakeoff('BRCA1', 'snp1kg', True)

    @timeout_decorator.timeout(200)        
    def test_map_brca1_cactus(self):
        """ Mapping and calling bakeoff F1 test for BRCA1 cactus graph """
        self._test_bakeoff('BRCA1', 'cactus', True)

    @timeout_decorator.timeout(900)        
    def test_full_brca2_primary(self):
        """ Indexing, mapping and calling bakeoff F1 test for BRCA2 primary graph """
        self._test_bakeoff('BRCA2', 'primary', False)

    @timeout_decorator.timeout(900)        
    def test_full_brca2_snp1kg(self):
        """ Indexing, mapping and calling bakeoff F1 test for BRCA2 snp1kg graph """
        self._test_bakeoff('BRCA2', 'snp1kg', False)

    @timeout_decorator.timeout(900)        
    def test_full_brca2_cactus(self):
        """ Indexing, mapping and calling bakeoff F1 test for BRCA2 cactus graph """
        self._test_bakeoff('BRCA2', 'cactus', False)

    @skip("skipping test to keep runtime down")
    @timeout_decorator.timeout(2000)        
    def test_map_sma_primary(self):
        """ Indexing, mapping and calling bakeoff F1 test for SMA primary graph """
        self._test_bakeoff('SMA', 'primary', True)

    @skip("skipping test to keep runtime down")        
    @timeout_decorator.timeout(2000)        
    def test_map_sma_snp1kg(self):
        """ Indexing, mapping and calling bakeoff F1 test for SMA snp1kg graph """
        self._test_bakeoff('SMA', 'snp1kg', True)

    @skip("skipping test to keep runtime down")        
    @timeout_decorator.timeout(2000)        
    def test_map_sma_cactus(self):
        """ Indexing, mapping and calling bakeoff F1 test for SMA cactus graph """
        self._test_bakeoff('SMA', 'cactus', True)

    @skip("skipping test to keep runtime down")         
    @timeout_decorator.timeout(2000)        
    def test_map_lrc_kir_primary(self):
        """ Indexing, mapping and calling bakeoff F1 test for LRC-KIR primary graph """
        self._test_bakeoff('LRC-KIR', 'primary', True)

    @skip("skipping test to keep runtime down")         
    @timeout_decorator.timeout(2000)        
    def test_map_lrc_kir_snp1kg(self):
        """ Indexing, mapping and calling bakeoff F1 test for LRC-KIR snp1kg graph """
        self._test_bakeoff('LRC-KIR', 'snp1kg', True)
        
    @skip("skipping test to keep runtime down")         
    @timeout_decorator.timeout(2000)        
    def test_map_lrc_kir_cactus(self):
        """ Indexing, mapping and calling bakeoff F1 test for LRC-KIR cactus graph """
        self._test_bakeoff('LRC-KIR', 'cactus', True)

    @timeout_decorator.timeout(10000)        
    def test_map_mhc_primary(self):
        """ Indexing, mapping and calling bakeoff F1 test for MHC primary graph """
        self._test_bakeoff('MHC', 'primary', True)

    @timeout_decorator.timeout(10000)        
    def test_map_mhc_snp1kg(self):
        """ Indexing, mapping and calling bakeoff F1 test for MHC snp1kg graph """
        self._test_bakeoff('MHC', 'snp1kg', True)

    @skip("skipping test to keep runtime down (baseline missing as well)")          
    @timeout_decorator.timeout(10000)        
    def test_map_mhc_cactus(self):
        """ Indexing, mapping and calling bakeoff F1 test for MHC cactus graph """
        self._test_bakeoff('MHC', 'cactus', True)
