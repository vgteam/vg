#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Mine the Jenkins test log XML and the test output files and generate a report
of how good VG is at various tasks.
"""
from __future__ import unicode_literals
import logging
import subprocess
import tempfile
import os
import sys
import re
import argparse
# We need to use cElementTree because normal ElementTree thinks perfectly
# respectable UTF-8 characters aren't allowed in UTF-8 XML.
import xml.etree.cElementTree as ET
import textwrap
import shutil
import datetime
import cgi
import io
import traceback
from collections import defaultdict



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
    assert toks[0] == 'test' and len(toks) >= 4
    return 'outstore-{}-{}-{}'.format(
        toks[1], toks[2].replace('lrc-kir', 'lrc_kir').upper(), '-'.join(toks[3:]))


def scrape_mapeval_runtimes(text):
    """ toil-vg mapeval uses the RealtimeLogger to print the running times of 
    each call to vg map (or vg mpmap).  We piece that information together
    here 
    
    I would rather this be in vgci.py, but don't have access to the log there. 
    So scrape it out here then add it to the table from vgci.py with join_runtimes below

    Todo: save the log into the workdir and mine from vgci.py
    """
    
    runtimes = defaultdict(int)
    for line in text.split('\n'):
        # we want to parse something like
        #host 2017-08-22 11:09:45,362 MainThread INFO toil-rt: Aligned /tmp/toil-55082/aligned-snp1kg_HG00096_0.gam. Process took 4.37375712395 seconds with single-end vg-map
        search_tok = 'Aligned'
        i = line.find(search_tok)
        if i >= 0 and line.find('Process took') > 0:
            try:
                toks = line[i + len(search_tok) + 1: ].rstrip().split()
                outputfile = toks[0]
                seconds = round(float(toks[3]), 1)
                read_type = toks[-2]
                method = toks[-1]
                # name will be something like tmpdir/aligned-primary_0.gam.
                # so we cut down to just primary_0
                name = os.path.basename(outputfile)[8:-5]
                # then strip the _0
                name = name[0:name.rfind('_')]
                # and any tags for paired end
                name = name.replace('-pe','').replace('-se','').replace('-mp','')
                key = method if 'bwa' in method else name
                if 'mpmap' in method:
                    key += '-mp'
                if 'single-end' in read_type:
                    key += '-se'
                elif 'paired-end' in read_type:
                    key += '-pe'
                runtimes[key] += seconds
            except:
                pass

    return runtimes
        
def join_runtimes(mapeval_table, runtime_dict):
    """ Tack on some mapping runtimes that we mined above to the mapeval table that was printed 
    by the jenkins tests """
    mapeval_table[0].append("Map Time (s)")
    for i in range(1, len(mapeval_table)):
        row = mapeval_table[i]
        if row[0].rstrip('*') in runtime_dict:
            row.append(runtime_dict[row[0].rstrip('*')])
        else:
            row.append(-1)

    return mapeval_table
    
def parse_begin_message(line):
    """ Return True, name, is_tsv if tag found otherwise False, None, None """
    try:
        if 'VGCI' in line:
            elem = ET.fromstring(line.encode('utf8').rstrip() + '</VGCI>')
            if elem.tag == 'VGCI':
                return True, elem.get('name'), elem.get('tsv', 'false').lower() == 'true'
    except:
        pass
    return False, None, None

def parse_end_message(line):
    """ Return True if line marks end of a report """
    try:
        if 'VGCI' in line:
            elem = ET.fromstring('<VGCI>' + line.encode('utf8').rstrip())
            if elem.tag == 'VGCI':
                return True
    except:
        pass
    return False

def scrape_messages(text):
    """ Use the tags above to extract blocks and tables from a buffer.
    A message is a tuple of (name, table) or (name, text) where name 
    can be none """

    messages = []
    msg = None
    name = None
    
    for line in text.split('\n'):
        is_msg, msg_name, is_tsv = parse_begin_message(line)
        # start message
        if is_msg:
            if msg is not None:
                messages.append((name, msg))
                msg, name = None, None
            name = msg_name            
            if is_tsv:
                msg = []
            else:
                msg = ''
        # end message
        elif parse_end_message(line):
            if msg is not None:
                messages.append((name, msg))
                msg, name = None, None
        # continue message
        elif msg is not None:
            if isinstance(msg, list):
                row = line.rstrip().split('\t')
                if len(row):
                    msg.append(row)
            else:
                msg += line
    if msg:
        messages.append((name, msg))
        msg, name = None, None

    return messages
    

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
    
    # We need to emit only Unicode things, but we may get str or Unicode
    # depending on if the parser thinks we have UTF-8 or ASCII data in a field.
    
    tc = dict()
    tc['name'] = unicode(testcase.get('name'))
    tc['time'] = int(float(testcase.get('time', 0)))

    if testcase.find('system-out') is not None:
        tc['stdout'] = unicode(testcase.find('system-out').text)
    else:
        tc['stdout'] = None

    if testcase.find('system-err') is not None:
        tc['stderr'] = unicode(testcase.find('system-err').text)
    else:
        tc['stderr'] = None

    if testcase.find('skipped') is not None:
        tc['skipped'] = True
        tc['skip-msg'] = unicode(testcase.find('skipped').text)
    else:
        tc['skipped'] = False
        tc['skip-msg'] = None

    failure = testcase.find('failure')
    if failure is not None:
        tc['fail-txt'] = unicode(failure.text)
        tc['fail-msg'] = unicode(failure.get('message'))
        tc['failed'] = True
    else:
        tc['fail-txt'] = None
        tc['fail-msg'] = None
        tc['failed'] = False

    return tc

def get_vgci_warnings(tc):
    """
    extract warnings from the parsed test case stderr
    """
    warnings = []
    if tc['stderr']:
        for line in tc['stderr'].split('\n'):
            # For each line
            if "WARNING" in line and "vgci" in line:
                # If it is a CI warning, keep it
                warnings.append(line.rstrip())
    return warnings

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
            md += ' for [PR {}]({})'.format(os.getenv('ghprbPullId'), os.getenv('ghprbPullLink').decode('utf8'))
        elif build_number:
            md += ' for merge to master'
        md += '.  View the [full report here]({{REPORT_URL}}).\n\n'

        # will have to modify this if we ever add another test suite
        testsuite = xml_root
        
        ts = parse_testsuite_xml(testsuite)

        md += '{} tests passed, {} tests failed and {} tests skipped in {} seconds\n\n'.format(
            ts['passes'], ts['fails'], ts['skips'], ts['time'])

        if ts['fails'] is not None and int(ts['fails']) >= 1:
            md += 'Failed tests:\n\n'

        warnings = []
        
        for testcase in testsuite.iter('testcase'):
            try:
                tc = parse_testcase_xml(testcase)
                if tc['failed']:
                    md += '* {} ({} seconds)\n'.format(tc['name'], tc['time'])
                warnings += get_vgci_warnings(tc)
            except:
                md += '**Error parsing Test Case XML**\n'
            
        if len(warnings) > 0:
            if ts['fails'] is not None and int(ts['fails']) >= 1:
                md += '\n'
            md += 'Tests produced {} warnings. {} were for lower-than-expected alignment scores\n\n'.format(
                len(warnings), len([w for w in warnings if 'negative score' in w]))
            
    except:
        md += ' **Error parsing Test Suite XML**\n'
    return md
    
def escape(string):
    """
    Return an HTML-escaped version of the given string. Should be called on
    every user-supplied string before being included in HTML output. Escapes
    quotes, so the string can be used as an HTML attribute value.
    """
    return cgi.escape(string, quote=True)
        
def html_header(xml_root):
    """
    Make an HTML header for the test suite XML
    """
    
    report = ''
    
    try:
        report += '''
<!DOCTYPE html><html><head>
<meta charset="utf-8"/>
<title>''' + escape(os.getenv('ghprbPullTitle', 'master').decode('utf8')) + ''': vg Test Report</title>
<style> table {
font-family: arial, sans-serif; border-collapse: collapse; width: 100%;}
td, dh { border: 1px solid #dddddd; text-align: left; padding: 8px;}
tr:nth-child(even) { background-color: #dddddd;}
html { background-image: url("https://source.unsplash.com/random"); background-size: cover; }
body { background-color: rgba(255, 255, 255, 0.8); width: 80%; padding: 1em; border-radius: 1em; margin: auto; margin-top: 1em; margin-bottom: 1em; }
h1, h2, h3, h4, h5, h6 { font-family: sans-serif; } 
</style></head><body>
'''

        # will have to modify this if we ever add another test suite
        testsuite = xml_root
        
        ts = parse_testsuite_xml(testsuite)

        build_number = os.getenv('BUILD_NUMBER')

        report += '<h2>vg Test Report'
        if os.getenv('ghprbPullId'):
            report += ' for <a href={}>PR {}: {}</a>'.format(os.getenv('ghprbPullLink'),
                                                             os.getenv('ghprbPullId'),
                                                             escape(os.getenv('ghprbPullTitle').decode('utf8')))
        elif build_number:
            report += ' for merge to master'
            
        if ts['fails'] is not None and int(ts['fails']) >= 1:
            # Some tests failed
            report += ' 👿' + '🔥' * int(ts['fails'])
        else:
            # All tests passed
            report += ' 😄'
            
        report += '</h2>\n'

        report += '<p> {} tests passed, {} tests failed and {} tests skipped in {} seconds'.format(
            ts['passes'], ts['fails'], ts['skips'], ts['time'])
        report += '.</p>\n'

        if build_number:
            report += '<p> <a href=http://jenkins.cgcloud.info/job/vg/{}/>'.format(build_number)
            report += 'Jenkins test page </a></p>\n'

        report += '<p> Report generated on {} </p>\n'.format(str(datetime.datetime.now()))
        
    except:
        report += ' Error parsing Test Suite XML\n'
    return report

def html_table(table, caption = None):
    """
    take a table (list of lists) as scraped above and write a HTML table.
    """
    t = '<table>\n'
    if caption:
        t += '<caption align="bottom"><i>{}</i></caption>\n'.format(caption)                
    for i, row in enumerate(table):
        t += '<tr>\n'
        tag = 'th' if i == 0 else 'td'
        for col in row:
            t += '<{}>{}</{}>\n'.format(tag, col, tag)
    t += '</table>\n'
    return t

def html_testcase(tc, work_dir, report_dir, max_warnings = 10):
    """
    Make an HTML report for a single test case
    """
    report = '\n<hr>\n'
    try:
        if tc['failed']:
            report += '<h3><font color="FF0000">{}</font></h3>\n'.format(tc['name'])
        else:
            report += '<h3>{}</h3>\n'.format(tc['name'])            
        if tc['failed']:
            report += '<p><b><font color="FF0000">Failed</font></b> in {} seconds</p>\n'.format(tc['time'])
            report += '<p>Failure Text</p>\n<pre>'
            for line in tc['fail-txt'].split('\n'):
                for subline in textwrap.wrap(line, 80):
                    report += escape(subline) + '\n'
            report += '</pre>\n'

            report += '<p>Failure Message</p>\n<pre>'
            for line in tc['fail-msg'].split('\n'):
                for subline in textwrap.wrap(line, 80):
                    report += escape(subline) + '\n'
            report += '</pre>\n'

        else:
            report += '<p>Passed in {} seconds</p>\n'.format(tc['time'])

        outstore = os.path.join(work_dir, testname_to_outstore(tc['name']))

        images = []
        captions = []
        baseline_images = []
        for plot_name in 'pr', 'pr.control', 'pr.primary.filter', 'pr.control.primary.filter', 'roc', 'qq', 'qq.control':
            plot_path = os.path.join(outstore, '{}.svg'.format(plot_name))
            if os.path.isfile(plot_path):
                new_name = '{}-{}.svg'.format(tc['name'], plot_name)
                shutil.copy2(plot_path, os.path.join(report_dir, new_name))
                images.append(new_name)
                captions.append(plot_name.upper())
                baseline_images.append(os.path.join(
                    'https://cgl-pipeline-inputs.s3.amazonaws.com/vg_cgl/vg_ci/jenkins_regression_baseline',
                    os.path.basename(outstore),
                    os.path.basename(plot_path)))

        if len(images) > 0:
            report += '<table class="image">\n'
            i = 0
            row_size = 2
            while i < len(images):
                report += '<tr>'
                for j in range(i, i + row_size):
                    if j < len(images):
                        report += '<td><img width=\"300px\" src=\"{}\"></td>'.format(images[j])
                report += '</tr>\n<tr>'
                for j in range(i, i + row_size):
                    if j < len(images):
                        report += '<td><a href=\"{}\">{} Plot</a>'.format(images[j], captions[j])
                        report += ' <a href=\"{}\">(baseline)</a></td>'.format(baseline_images[j])
                report += '</tr>\n'
                i += row_size
            report += '</table>\n'

        # extract running times from stderr (only for mapeval)
        map_time_table = None
        if tc['stderr'] and 'sim' in tc['name']:
            map_time_table = scrape_mapeval_runtimes(tc['stderr'])
            
        if tc['stdout']:
            # extract only things in <VCGI> tags from stdout
            messages = scrape_messages(tc['stdout'])
            for message in messages:                
                name, body = message[0], message[1]                    
                if isinstance(body, list):
                    if name and name.startswith('map eval results') and map_time_table:
                        body = join_runtimes(body, map_time_table)
                    report += html_table(body, name)
                else:
                    if name:
                        report += '<h5>{}</h5>\n'.format(name)                    
                    report += '<pre>'
                    for line in body.split('\n'):
                        for subline in textwrap.wrap(line, 80):
                            report += escape(subline) + '\n'
                    report += '</pre>\n'

        if tc['stderr']:
            # Look for warning lines in stderr
            warnings = get_vgci_warnings(tc)

            # report top warnings in HTML
            if len(warnings) > 0:
                report += '<p>{} CI warnings found'.format(len(warnings))
                if len(warnings) > max_warnings:
                    report += '; showing first {}'.format(max_warnings)
                report += ':'
                report += '</p>\n'
                for warning in warnings[:max_warnings]:            
                    report += '<pre>{}</pre>\n'.format(escape('\n'.join(textwrap.wrap(warning, 80))))

                # warnings file
                if len(warnings) > max_warnings:
                    warn_name = '{}-warnings.txt'.format(tc['name'])
                    with io.open(os.path.join(report_dir, warn_name), 'w', encoding='utf8') as warn_file:
                        for warning in warnings:
                            warn_file.write(warning + '\n')

            # entire stderr output (which also includes warnings, should we filter?)
            err_name = '{}-stderr.txt'.format(tc['name'])
            with io.open(os.path.join(report_dir, err_name), 'w', encoding='utf8') as err_file:
                err_file.write(tc['stderr'])

            # link to warnings and stderr
            report += '<p>'
            if len(warnings) > max_warnings:
                report += '<a href={}>All CI Warnings</a>, '.format(warn_name)
            report += '<a href={}>Standard Error</a></p>\n'.format(err_name)

    except int as e:
        report += 'Error parsing Test Case XML\n'

    return report

def write_html_report(xml_root, work_dir, html_dir, html_name = 'index.html'):
    """ Write the HTML report in a given directory
    """
    if not os.path.isdir(html_dir):
        os.makedirs(html_dir)
        
    html_path = os.path.join(html_dir, html_name)
    with io.open(html_path, 'w', encoding='utf8') as html_file:
        header = html_header(xml_root)
        html_file.write(header)
        
        parsed_testcases = []
        for testcase in xml_root.iter('testcase'):
            try:
                tc = parse_testcase_xml(testcase)
                parsed_testcases.append(tc)
            except Exception as err:
                logging.warning('Unexpected error parsing testcase XML {}'.format(testcase))
                print(traceback.format_exc())

        # sort test cases so failed first, then sim, the by name
        def sort_key(tc):
            key = ('0' if 'sim' in tc['name'] else '1') 
            key = ('0' if tc['failed'] else '1') + key
            return key + tc['name']
                
        for tc in sorted(parsed_testcases, key=sort_key):
            if not tc['skipped']:
                tc_body = html_testcase(tc, work_dir, html_dir)
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
    with io.open(options.md_out, 'w', encoding='utf8') as md_file:
        md_file.write(markdown)

    # Write our HTML report
    write_html_report(xml_root, options.work_dir, options.html_out_dir)
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        

