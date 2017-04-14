#!/usr/bin/env python

"""
Print the latest tag of a Quay repository.  If no tag exists that's new enough, print "None"
"""
import os
import sys
import requests
from dateutil import parser, relativedelta, tz
import datetime
import json
import argparse

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('quay_repo', help='name of repository on quay.io')
    parser.add_argument('--max-age', help='maximum age in minutes of repository',
                        type=int, default=60)
    parser.add_argument('--tag-suffix', help='only consider tags that end with suffix',
                        default='-run')
    args = args[1:]        
    return parser.parse_args(args)

def latest_tag_info(repository, tag_suffix):
    """ Returns the latest tag and its age in minutes as a tuple
    """
    
    # pull in all the tags that end with suffix
    tags = []
    get_tags = None
    page = 1
    while not get_tags or response['has_additional'] is True:
        get_tags = 'https://quay.io/api/v1/repository/{}/tag/?page={}'.format(repository, page)
        response = requests.get(get_tags).json()
        for tag in response['tags']:
            if tag['name'].endswith(tag_suffix):
                tags.append(tag['name'])
        page += 1

    # find the most recent one by looking at the metadata for each
    # (because I'm not sure if API guarantees any kind of chronological ordering)
    newest_tag = None
    newest_age = sys.maxint
    for tag in tags:
        get_image = 'https://quay.io/api/v1/repository/{}/tag/{}/images'.format(repository, tag)
        response = requests.get(get_image).json()
    
        # use dateutil for timezone-aware comparison
        creation_time = parser.parse(response['images'][0]['created'])
        now = datetime.datetime.now(tz.tzlocal())
        delta = relativedelta.relativedelta(now, creation_time)
        age_minutes = delta.years * 525600 + delta.days * 1440 + delta.hours * 60 + delta.minutes

        if age_minutes < newest_age:
            newest_tag = tag
            newest_age = age_minutes
    
    return newest_tag, newest_age


def main(args):

    options = parse_args(args)

    try:
        tag, age = latest_tag_info(options.quay_repo, options.tag_suffix)
        if age <= options.max_age:
            print 'quay.io/{}:{}'.format(options.quay_repo, tag)
            return 0
    except int as e:
        pass
    print "None"
    return 1
    

if __name__ == "__main__" :
    sys.exit(main(sys.argv))

