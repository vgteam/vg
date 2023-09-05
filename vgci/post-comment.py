#!/usr/bin/python3
# post-comment.py: post standard input as a Github comment

"""
Post a Github comment.
"""

import argparse
import os
import sys

import github

def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    
    parser.add_argument("--in_file", type=argparse.FileType("r"),
        default=sys.stdin,
        help="File to read comment Markdown from (default: standard input)")
    parser.add_argument("--gh_token", default=os.environ.get("GH_TOKEN"),
        help="Github token with public_repo permission to use to log in (default: GH_TOKEN env var)")
    parser.add_argument("repo",
        help="repo (user/repoName) to comment on things in")
    # Then we need either a PR number or a commmit hash to comment on
    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument("--commit",
        help="commit hash to comment on")
    target_group.add_argument("--pr", type=int,
        help="PR issue number to comment on")
    
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
   
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    if options.gh_token is None:
        raise RuntimeError("No Github token found; cannot post comments")
    
    # Connect to Github
    client = github.Github(options.gh_token)
    
    # Get the repo
    repo = client.get_repo(options.repo)
    
    if options.commit is not None:
        # Comment on this commit
        target = repo.get_commit(options.commit)
        comment = target.create_comment(options.in_file.read())
    elif options.pr is not None:
        # Comment on this PR
        target = repo.get_pull(options.pr)
        # Need to create_issue_comment because create_comment is for fancy line
        # by line commenting
        comment = target.create_issue_comment(options.in_file.read())
    else:
        raise RuntimeError("No comment target specified")
        
    
    # Report that it worked.
    # Need to use html_url for user reporting because url is just JSON.
    print(("Created comment: {}".format(comment.html_url)))
    
    return 0
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

