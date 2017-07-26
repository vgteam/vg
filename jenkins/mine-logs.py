#!/usr/bin/env python2.7
"""
Mine the Jenkins test log XML and the test output files and generate a report
of how good VG is at various tasks.
"""
import subprocess
import tempfile
import os, sys
import argparse
import xml.etree.ElementTree as ET
import textwrap

def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__, 
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # General options
    parser.add_argument('xml_in',
                        help='XML result from PyTest in JUnit format')
    parser.add_argument('work_dir',
                        help='vgci work dir')
    parser.add_argument('html_out',
                        help='output html')
    parser.add_argument('md_out',
                        help='output markdown')

    args = args[1:]
        
    return parser.parse_args(args)

def testname_to_outstore(test_name):
    """
    convert something like test_full_brca2_cactus to oustore-full-BRCA2-cactus
    """
    toks = test_name.replace('lrc_kir', 'lrc-kir').split('_')
    assert toks[0] == 'test' and len(toks) == 4
    return 'outstore-{}-{}-{}'.format(
        toks[1], toks[2].replace('lrc-kir', 'lrc_kir').upper(), toks[3])
    
def parse_junit_xml(xml_path, work_dir, suite_html_fn, case_html_fn, out_html,
                    suite_md_fn, case_md_fn, out_md):
    """
    Scan through the JUnit Report.  Run the callbacks on the testsuitres as
    well as each testcase for markdown and html
    """
    xml_tree = ET.parse(xml_path)
    xml_root = xml_tree.getroot()
    
    for testsuite in xml_root.iter('testsuite'):
        ts_name = testsuite.attrib['name']
        ts_tests = int(testsuite.attrib['tests'])
        ts_fails = int(testsuite.attrib['failures'])
        ts_skips = int(testsuite.attrib['skips'])
        ts_time = int(float(testsuite.attrib['time']))

        if suite_html_fn:
            suite_html_fn(out_html, ts_name, ts_tests, ts_fails, ts_skips, ts_time)
        if suite_md_fn:
            suite_md_fn(out_md, ts_name, ts_tests, ts_fails, ts_skips, ts_time)
    
        for testcase in xml_root.iter('testcase'):
            tc_name = testcase.attrib['name']
            tc_out = testcase.find('system-out').text
            tc_err = testcase.find('system-err').text
            tc_time = int(float(testcase.attrib['time']))
            failure = testcase.find('failure')
            if failure is not None:
                tc_failure_txt = failure.text
                tc_failure_msg = failure.attrib['message']
            else:
                tc_failure_txt = None
                tc_failure_msg = None

            if case_html_fn:
                case_html_fn(out_html, tc_name, os.path.join(work_dir, testname_to_outstore(tc_name)),
                        tc_out, tc_err, tc_time, tc_failure_txt, tc_failure_msg)
            if case_md_fn:
                case_md_fn(out_md, tc_name, os.path.join(work_dir, testname_to_outstore(tc_name)),
                        tc_out, tc_err, tc_time, tc_failure_txt, tc_failure_msg)
            

def write_md_header(md_file):
    """
    Write the top of the MarkDown file
    """
    md_file.write('## vg Test Report Summary\n')
    md_file.write('The full report')

    if os.getenv('ghprbPullTitle'):
        md_file.write(' for [PR \'{}\']({})'.format(os.getenv('ghprbPullTitle'), os.getenv('ghprbPullLink')))
    elif os.getenv('ghprbActualCommit'):
        md_file.write(' for commit {}'.format('ghprbActualCommit'))

    md_file.write(' is available here [here]({{REPORT_URL}})\n\n')


def write_md_testsuite(md_file, name, num_tests, num_fails, num_skips, total_time):
    """
    Write some statistics for the whole test suite
    """
    md_file.write('{} tests passed, {} tests failed and {} tests skipped in {} seconds\n\n'.format(
        num_tests-num_fails-num_skips, num_fails, num_skips, total_time))
    
def write_md_testcase(html_file, name, outstore, stdout, stderr, seconds, fail_txt, fail_msg):
    """
    Given some information extracted from the XML, write a MarkDown summary
    """
    md = '### {}\n\n'.format(name)
    if fail_txt is None:
        md += 'Passed in {} seconds\n\n'.format(seconds)
    else:
        md += 'Failed in {} seconds\n\n'.format(seconds)
        md += 'Failure Message: `{}`\n\n'.format(fail_msg)
    md += 'Standard Output:\n\n'
    for line in textwrap.wrap(stdout, 80):
        md += '     ' + line + '\n'
    md += '\n'
    html_file.write(md)


def write_html_header(html_file):
    """
    Write the top of the HTML file
    """
    html_file.write('''\
<!DOCTYPE html>
<html>
    <head>
        <title>vg Test Report</title>
    </head>
    <body>
''')
    html_file.write('<h2>vg Test Report</h2>\n')
    if os.getenv('ghprbPullTitle'):
        html_file.write('Pull Request \"{}\" <a href="{}"> (Link) </a>\n'.format(
            os.getenv('ghprbPullTitle'), os.getenv('ghprbPullLink')))
    elif os.getenv('ghprbActualCommit'):
        html_file.write('<p>Commit ID {}</p>'.format('ghprbActualCommit'))
    
def write_html_testsuite(html_file, name, num_tests, num_fails, num_skips, total_time):
    """
    Write some statistics for the whole test suite
    """
    html_file.write('<p>{} tests passed, {} tests failed and {} tests skipped in {} seconds</p>\n'.format(
        num_tests-num_fails-num_skips, num_fails, num_skips, total_time))

def write_html_testcase(html_file, name, outstore, stdout, stderr, seconds, fail_txt, fail_msg):
    """
    Given some information extracted from the XML, write a HTML view
    """
    html_file.write('<h3>{}</h3>\n'.format(name))
    if fail_txt is None:
        html_file.write('<p>Passed in {} seconds</p>\n'.format(seconds))
    else:
        html_file.write('<p>Failed in {} seconds</p>\n'.format(seconds))
        html_file.write('<p>Failure Message: `{}`</p>\n'.format(fail_msg))
    html_file.write('<p>Standard Output:</p>\n<p>')
    for line in textwrap.wrap(stdout, 80):
        html_file.write(line + '&#10')
    html_file.write('</p>')
    

def main(args):
    """
    Scrape out some information from the PyTest logs.  Present as MarkDown summary
    and HTML details
    """    
    options = parse_args(args)

    # just scrape stdout and a few stats for for now 
    with open(options.md_out, 'w') as md_file, \
         open(options.html_out, 'w') as html_file:

        write_md_header(md_file)
        write_html_header(html_file)
                                  
        parse_junit_xml(options.xml_in, options.work_dir,
                        write_html_testsuite, write_html_testcase, html_file,
                        write_md_testsuite, write_md_testcase, md_file)

        html_file.write('</body>\n</html>\n')
    
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        

