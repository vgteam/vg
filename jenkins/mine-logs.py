#!/usr/bin/env python2.7
"""
Mine the Jenkins test log XML and the test output files and generate a report
of how good VG is at various tasks.
"""
import logging
import subprocess
import tempfile
import os, sys
import argparse
import xml.etree.ElementTree as ET
import textwrap
import shutil
import datetime

def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__, 
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # General options
    parser.add_argument('xml_in',
                        help='XML result from PyTest in JUnit format')
    parser.add_argument('work_dir',
                        help='vgci work dir')
    parser.add_argument('html_out_dir',
                        help='output html directory')
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


def parse_testsuite_xml(testsuite):
    """
    Flatten fields of interset from a TestSuite XML element into a dict
    """
    ts = dict()
    ts['name'] = testsuite.get('name')
    ts['tests'] = int(testsuite.get('tests', 0))
    ts['fails'] = int(testsuite.get('failures', 0))
    ts['skips'] = int(testsuite.get('skips', 0))
    ts['errors'] = int(testsuite.get('errors', 0))
    ts['passes'] = ts['tests'] - ts['fails'] - ts['skips']
    ts['time'] = int(float(testsuite.get('time', 0)))
    return ts
        
    
def parse_testcase_xml(testcase):
    """
    Flatten fields of interest from a TestCase XML element into a dict,
    where anything not found is None
    """
    tc = dict()
    tc['name'] = testcase.get('name')
    tc['time'] = int(float(testcase.get('time', 0)))

    if testcase.find('system-out') is not None:
        tc['stdout'] = testcase.find('system-out').text
    else:
        tc['stdout'] = None

    if testcase.find('system-err') is not None:
        tc['stderr'] = testcase.find('system-err').text
    else:
        tc['stderr'] = None

    failure = testcase.find('failure')
    if failure is not None:
        tc['fail-txt'] = failure.text
        tc['fail-msg'] = failure.get('message')
        tc['failed'] = True
    else:
        tc['fail-txt'] = None
        tc['fail-msg'] = None
        tc['failed'] = False

    return tc

def md_summary(xml_root):
    """
    Make a brief summary in Markdown of a testsuite
    """
    build_number = os.getenv('BUILD_NUMBER')
    if build_number:
        md = '[Jenkins vg tests](http://jenkins.cgcloud.info/job/vg/{}/) complete'.format(
            build_number)
    else:
        md = 'Jenkins vg tests complete'
    try:
        if os.getenv('ghprbPullId'):
            md += ' for [PR {}]({})'.format(os.getenv('ghprbPullId'), os.getenv('ghprbPullLink'))
        elif build_number:
            md += ' for merge to master'
        md += '.  View the [full report here]({{REPORT_URL}}).\n\n'

        # will have to modify this if we ever add another test suite
        testsuite = xml_root
        
        ts = parse_testsuite_xml(testsuite)

        md += '{} tests passed, {} tests failed and {} tests skipped in {} seconds\n\n'.format(
            ts['passes'], ts['fails'], ts['skips'], ts['time'])

        if ts['fails'] is not None and int(ts['fails']) > 1:
            md += 'Failed tests:\n\n'

        for testcase in testsuite.iter('testcase'):
            try:
                tc = parse_testcase_xml(testcase)
                if tc['failed']:
                    md += '* {} ({} seconds)\n'.format(tc['name'], tc['time'])
            except:
                md += '**Error parsing Test Case XML**\n'                
    except:
        md += ' **Error parsing Test Suite XML**\n'
    return md
        
def html_header(xml_root):
    """
    Make an HTML header for the test suite XML
    """
    try:
        report = '<!DOCTYPE html><html><head>\n'
        report += '<title>vg Test Report</title>\n'
        report += '</head><body>\n'

        # will have to modify this if we ever add another test suite
        testsuite = xml_root
        
        ts = parse_testsuite_xml(testsuite)

        build_number = os.getenv('BUILD_NUMBER')

        report += '<h2>vg Test Report'
        if os.getenv('ghprbPullId'):
            report += ' for <a href={}>PR {}</a>'.format(os.getenv('ghprbPullLink'),
                                                         os.getenv('ghprbPullId'))
        elif build_number:
            report += ' for merge to master'
        report += '</h2>\n'

        report += '<p> {} tests passed, {} tests failed and {} tests skipped in {} seconds'.format(
            ts['passes'], ts['fails'], ts['skips'], ts['time'])
        report += '.</p>\n'

        if build_number:
            report += '<p> <a href=http://jenkins.cgcloud.info/job/vg/{}/>'.format(build_number)
            report += 'Jenkins test page </a></ap>\n'

        report += '<p> Report generated on {} </p>\n'.format(str(datetime.datetime.now()))
        
    except:
        report += ' Error parsing Test Suite XML\n'
    return report

def html_testcase(testcase, work_dir, report_dir):
    """
    Make an HTML report for a single test case
    """
    report = ''
    try:
        tc = parse_testcase_xml(testcase)

        report += '<h3>{}</h3>\n'.format(tc['name'])
        if tc['failed']:
            report += '<p>Failed in {} seconds</p>\n'.format(tc['time'])
            report += '<p>Failure Message: `{}`</p>\n'.format(tc['fail-msg'])
            report += '<p>Failure Text: `{}`</p>\n'.format(tc['fail-txt'])
        else:
            report += '<p>Passed in {} seconds</p>\n'.format(tc['time'])

        outstore = os.path.join(work_dir, testname_to_outstore(tc['name']))

        for plot_name in 'roc', 'qq':
            plot_path = os.path.join(outstore, '{}.pdf'.format(plot_name))
            if os.path.isfile(plot_path):
                new_name = '{}-{}.pdf'.format(tc['name'], plot_name)
                shutil.copy2(plot_path, os.path.join(report_dir, new_name))
                report += '<p><a href={}>{} Plot</a></p>\n'.format(new_name, plot_name.upper())

        if tc['stdout']:
            report += '<p>Standard Output:</p>\n<p>'
            for line in textwrap.wrap(tc['stdout'], 80):
                report += line + '&#10'
            report += '</p>\n'

        if tc['stderr']:
            err_name = '{}-stderr.txt'.format(tc['name'])
            with open(os.path.join(report_dir, err_name), 'w') as err_file:
                err_file.write(tc['stderr'])
            report += '<p><a href={}>Standard Error</a></p>\n'.format(err_name)

    except:
        report += 'Error parseing Test Case XML\n'

    return report

def write_html_report(xml_root, work_dir, html_dir, html_name = 'index.html'):
    """ Write the HTML report in a given directory
    """
    if not os.path.isdir(html_dir):
        os.makedirs(html_dir)
    
    html_path = os.path.join(html_dir, html_name)
    with open(html_path, 'w') as html_file:
        header = html_header(xml_root)
        html_file.write(header)
        for testcase in xml_root.iter('testcase'):
            tc_body = html_testcase(testcase, work_dir, html_dir)
            html_file.write(tc_body)
        html_file.write('</body>\n</html>\n')
        

def main(args):
    """
    Scrape out some information from the PyTest logs.  Present as MarkDown summary
    and HTML details
    """    
    options = parse_args(args)

    # Load up the XML
    xml_tree = ET.parse(options.xml_in)
    xml_root = xml_tree.getroot()

    # Write our Markdown summary
    markdown = md_summary(xml_root)
    with open(options.md_out, 'w') as md_file:
        md_file.write(markdown)

    # Write our HTML report
    write_html_report(xml_root, options.work_dir, options.html_out_dir)
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        

