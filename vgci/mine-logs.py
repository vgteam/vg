#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Mine the CI test log XML and the test output files and generate a report
of how good VG is at various tasks.
"""

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
import html
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

def load_mapeval_runtimes(map_times_path):
    """
    read the runtimes out of map_times.txt, assuming first line is a header
    """
    try:
        map_times_dict = {}
        with open(map_times_path) as map_times_file:
            lines = [line for line in map_times_file]
            for line in lines[1:]:
                toks = line.split('\t')
                map_times_dict[toks[0]] = toks[1]
        return map_times_dict
    except:
        pass
    return None        
        
def join_runtimes(mapeval_table, runtime_dict):
    """ Tack on some mapping runtimes that we mined above to the mapeval table that was printed 
    by the vgci tests """
    mapeval_table[0].append("Map Time (s)")
    for i in range(1, len(mapeval_table)):
        row = mapeval_table[i]
        row_name = row[0].rstrip('*')
        if row_name.endswith('-se') and row_name not in runtime_dict:
            row_name = row_name[:-3]
        if row_name in runtime_dict:
            row.append(runtime_dict[row_name])
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
    

def parse_all_testsuite_xml(xml_root):
    """
    Flatten fields of interest from all <testsuite> elements in the document into a dict
    """
   
    # We want to aggregate test counts and times over all suites.
    ts = {
        'tests': 0,
        'fails': 0,
        'skips': 0,
        'errors': 0,
        'time': 0.0
    }
   
    for testsuite in xml_root.iter('testsuite'):
        # Add in each test suite
        ts['tests'] += int(testsuite.get('tests', 0))
        ts['fails'] += int(testsuite.get('failures', 0))
        ts['skips'] += int(testsuite.get('skips', 0))
        ts['errors'] += int(testsuite.get('errors', 0))
        ts['time'] += float(testsuite.get('time', 0))
       
    # Compute passing test count
    ts['passes'] = ts['tests'] - ts['fails'] - ts['skips']
    
    # Round off the time to integer seconds
    ts['time'] = int(ts['time'])
    
    return ts
        
    
def parse_testcase_xml(testcase):
    """
    Flatten fields of interest from a TestCase XML element into a dict,
    where anything not found is None
    """
    
    # We need to emit only Unicode things, but we may get str or Unicode
    # depending on if the parser thinks we have UTF-8 or ASCII data in a field.
    
    tc = dict()
    tc['name'] = str(testcase.get('name'))
    tc['time'] = int(float(testcase.get('time', 0)))

    if testcase.find('system-out') is not None:
        tc['stdout'] = str(testcase.find('system-out').text)
    else:
        tc['stdout'] = None

    if testcase.find('system-err') is not None:
        tc['stderr'] = str(testcase.find('system-err').text)
    else:
        tc['stderr'] = None

    if testcase.find('skipped') is not None:
        tc['skipped'] = True
        tc['skip-msg'] = str(testcase.find('skipped').text)
    else:
        tc['skipped'] = False
        tc['skip-msg'] = None

    failure = testcase.find('failure')
    if failure is not None:
        tc['fail-txt'] = str(failure.text)
        tc['fail-msg'] = str(failure.get('message'))
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
    
    # Get the URL of the pipeline being run, or None
    pipeline_url = os.getenv('CI_PIPELINE_URL')
    # And the branch being built
    branch = os.getenv('CI_COMMIT_REF_NAME')
    # Check if on CI
    in_ci = bool(os.getenv('CI'))
    
    if pipeline_url:
        md = '[vg CI tests]({})'.format(pipeline_url)
    else:
        md = 'vg CI tests'
    if xml_root:
        md += ' complete'
    else:
        md += ' never ran'
        if pipeline_url:
            md += ' ([Check the logs for setup or build errors]'
            md += '({}))'.format(pipeline_url)
    try:
        if branch == 'master':
            md += ' for merge to master'
        elif branch is not None:
            md += ' for branch {}'.format(escape(branch))
        elif in_ci:
            md += ' for no branch'
        md += '.  View the [full report here]({{REPORT_URL}}).\n\n'

        ts = parse_all_testsuite_xml(xml_root)

        md += '{} tests passed, {} tests failed and {} tests skipped in {} seconds\n\n'.format(
            ts['passes'], ts['fails'], ts['skips'], ts['time'])

        if ts['fails'] is not None and int(ts['fails']) >= 1:
            md += 'Failed tests:\n\n'

        warnings = []
        
        for testcase in xml_root.iter('testcase'):
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
        if xml_root:
            md += ' **Error parsing Test Suite XML**\n'
    return md
    
def escape(string):
    """
    Return an HTML-escaped version of the given string. Should be called on
    every user-supplied string before being included in HTML output. Escapes
    quotes, so the string can be used as an HTML attribute value.
    """
    return html.escape(string, quote=True)
        
def html_header(xml_root):
    """
    Make an HTML header for the test suite XML
    """
    
    report = ''
    
    report += '''
<!DOCTYPE html><html><head>
<meta charset="utf-8"/>
<title>''' + escape(os.getenv('CI_COMMIT_REF_NAME', 'master')) + ''': vg Test Report</title>
<style> table {
font-family: arial, sans-serif; border-collapse: collapse; width: 100%;}
td, dh { border: 1px solid #dddddd; text-align: left; padding: 8px;}
tr:nth-child(even) { background-color: #dddddd;}
html { background-image: url("https://source.unsplash.com/random"); background-size: cover; }
body { background-color: rgba(255, 255, 255, 0.8); width: 80%; padding: 1em; border-radius: 1em; margin: auto; margin-top: 1em; margin-bottom: 1em; }
h1, h2, h3, h4, h5, h6 { font-family: sans-serif; } 
</style></head><body>
'''

    # Get the URL of the pipeline being run, or None
    pipeline_url = os.getenv('CI_PIPELINE_URL')
    # And the branch being built
    branch = os.getenv('CI_COMMIT_REF_NAME')
    # Check if on CI
    in_ci = bool(os.getenv('CI'))

    report += '<h2>vg Test Report'
    if branch == 'master':
        report += ' for merge to master'
    elif branch is not None:
        report += ' for branch {}'.format(escape(branch))
    elif in_ci:
        report += ' for no branch'

    if xml_root:
        try:
            ts = parse_all_testsuite_xml(xml_root)
        except:
            report += ' Error parsing Test Suite XML\n'
            ts = defaultdict(lambda: None)
            
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
    else:
        # Tests never ran
        report += ' 👿'
        report += '</h2>\n'
        report += '<p> tests never ran. </p>'
        if pipeline_url:
            report += '<p><a href="{}">'.format(pipeline_url)
            report += 'Check the logs for setup or build errors. </a></p>\n'

    if pipeline_url:
        report += '<p> <a href="{}">'.format(pipeline_url)
        report += 'CI test page </a></p>\n'

    report += '<p> Report generated on {} </p>\n'.format(str(datetime.datetime.now()))
        
    return report

def html_table(table, caption = None, rounding = None):
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
            try:
                rcol = round(float(col), rounding) if rounding else col
            except:
                rcol = col
            t += '<{}>{}</{}>\n'.format(tag, rcol, tag)
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
        possible_plot_names = ['pr', 'pr.control', 'pr.primary.filter', 'pr.control.primary.filter', 'roc', 'qq',
            'qq.control', 'roc-snp', 'roc-non_snp', 'roc-weighted']
        possible_plot_directories = ['', 'plots/']
        for plot_name in possible_plot_names:
            for plot_directory in possible_plot_directories:
                plot_path = os.path.join(outstore, '{}{}.svg'.format(plot_directory, plot_name))
                if os.path.isfile(plot_path):
                    new_name = '{}-{}.svg'.format(tc['name'], plot_name)
                    shutil.copy2(plot_path, os.path.join(report_dir, new_name))
                    images.append(new_name)
                    captions.append(plot_name.upper())
                    baseline_images.append(os.path.join(
                        's3://vg-data/vg_ci/vgci_regression_baseline',
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

        # extract running times from map_times.tsv (only for mapeval)
        map_time_table = None
        if tc['stderr'] and 'sim' in tc['name']:
            table_path = os.path.join(outstore, 'map_times.tsv')
            if os.path.isfile(table_path):
                map_time_table = load_mapeval_runtimes(table_path)

        stdout_name, err_name, warn_name = None, None, None
                
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

            # entire stdout output (since for some reason mapping logs only there sometimes)
            stdout_name = '{}-stdout.txt'.format(tc['name'])
            with io.open(os.path.join(report_dir, stdout_name), 'w', encoding='utf8') as out_file:
                out_file.write(tc['stdout'])

        # throw the calling times into their own table
        table_path = os.path.join(outstore, 'call_times.tsv')
        if os.path.isfile(table_path):
            with open(table_path) as table_file:
                call_time_table = [line.strip().split('\t') for line in table_file]
                report += '<p>\n{}\n</p>'.format(html_table(call_time_table, 'calling times', rounding=3))
                
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

        # link to warnings and stderr and stdout
        if stdout_name or warn_name or err_name:
            report += '<p>'
            if stdout_name:
                report += '<a href={}>Standard Output</a></p>\n'.format(stdout_name)
            if warn_name and len(warnings) > max_warnings:
                report += '<a href={}>All CI Warnings</a>, '.format(warn_name)
            if err_name:
                report += '<a href={}>Standard Error</a></p>\n'.format(err_name)
            report += '</p>'

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

        if xml_root:
            parsed_testcases = []
            for testcase in xml_root.iter('testcase'):
                try:
                    tc = parse_testcase_xml(testcase)
                    parsed_testcases.append(tc)
                except Exception as err:
                    logging.warning('Unexpected error parsing testcase XML {}'.format(testcase))
                    print((traceback.format_exc()))

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
    try:
        xml_tree = ET.parse(options.xml_in)
        xml_root = xml_tree.getroot()
    except Exception as e:
        xml_root = None
    
    # Write our Markdown summary
    markdown = md_summary(xml_root)
    with io.open(options.md_out, 'w', encoding='utf8') as md_file:
        md_file.write(markdown)

    # Write our HTML report
    write_html_report(xml_root, options.work_dir, options.html_out_dir)
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        

