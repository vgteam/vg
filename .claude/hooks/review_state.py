#!/usr/bin/env python3
"""
Remember which states of the code have had an automated review, and remind
Claude Code to review code that has not been reviewed before it is sent off.

A state is a git tree: the tree that committing every change to a tracked or
staged file would produce. Committing reviewed changes therefore keeps them
reviewed, while any further edit makes the state new and unreviewed. Untracked
files that were never staged are not part of the state, and neither is a
submodule's checked-out commit until it is staged, because inspecting every
submodule checkout takes seconds.

Reviewed trees are listed, one per line, in the file `claude-reviewed-trees`
in the git directory of the current worktree, so the record never shows up in
`git status` and is never pushed.

Usage:
    review_state.py mark          Record the current state as reviewed.
    review_state.py status        Print "reviewed" or "unreviewed" and exit 0 or 1.
    review_state.py prompt-hook   UserPromptSubmit hook; reads hook JSON on stdin.
    review_state.py send-hook     PreToolUse hook for Bash; reads hook JSON on stdin.

The hooks never fail the action they watch: on any git error they stay silent.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

RECORD_NAME = "claude-reviewed-trees"

REVIEW_SKILL = "exhaustive-code-review"

# Phrases by which a user says code is done, should go out, or is about to get
# their own review. A false match only costs one injected reminder.
READY_PATTERN = re.compile(
    r"pull request|\bPRs?\b|\bready (for|to) (review|merge|land|ship|go|push)"
    r"|\b(code|it|this|branch|change|changes) (is|are|looks?) (ready|done|good to go)"
    r"|\bgood to go\b|\bship it\b|\bsend (it|this|them) (off|out|upstream)"
    r"|\b(I'll|I will|let me|I'm going to|I want to) (review|look over|read over|check over)"
    r"|\bpush (it|this|the branch)\b",
    re.IGNORECASE,
)

# Shell commands that send code to other people.
SEND_PATTERN = re.compile(r"\bgh\s+pr\s+(create|ready)\b|\bgit\s+push\b")


def git(*args, env=None):
    """Run git and return its stripped standard output, or raise CalledProcessError."""
    return subprocess.run(
        ["git", *args], check=True, capture_output=True, text=True, env=env
    ).stdout.strip()


def record_path():
    return git("rev-parse", "--path-format=absolute", "--git-path", RECORD_NAME)


def working_tree_hash():
    """
    Return the tree that committing all tracked and staged changes would make.

    Stage into a copy of the index so the user's own index is left alone.
    """
    index = git("rev-parse", "--path-format=absolute", "--git-path", "index")
    fd, scratch_index = tempfile.mkstemp(
        prefix="claude-review-index-", dir=os.path.dirname(index)
    )
    os.close(fd)
    try:
        shutil.copyfile(index, scratch_index)
        env = dict(os.environ, GIT_INDEX_FILE=scratch_index)
        git("-c", "diff.ignoreSubmodules=all", "add", "--update", "--", ":/", env=env)
        return git("write-tree", env=env)
    finally:
        os.unlink(scratch_index)


def head_tree_hash():
    return git("rev-parse", "HEAD^{tree}")


def reviewed_trees():
    try:
        with open(record_path()) as record:
            return {line.strip() for line in record if line.strip()}
    except FileNotFoundError:
        return set()


def mark():
    tree = working_tree_hash()
    if tree not in reviewed_trees():
        with open(record_path(), "a") as record:
            record.write(tree + "\n")
    print(f"Recorded tree {tree} as reviewed.")


def status():
    reviewed = working_tree_hash() in reviewed_trees()
    print("reviewed" if reviewed else "unreviewed")
    return 0 if reviewed else 1


def prompt_hook(hook_input):
    if not READY_PATTERN.search(hook_input.get("prompt", "")):
        return
    if working_tree_hash() in reviewed_trees():
        return
    context = (
        "The code in this checkout has not had an automated review in its "
        "current state. If the user is asking to prepare a pull request, to "
        "push, or says the code is ready or that they will review it now, "
        f"run the `{REVIEW_SKILL}` skill on the branch's changes FIRST, before "
        "preparing the PR or handing the code over, and tell the user you are "
        "doing so. Skip it only if the user has explicitly said not to review. "
        "This note is added whenever the prompt matches loose keywords. If "
        "the prompt is not about sending or handing over code, act as if this "
        "note were absent, and do not mention it or this hook to the user."
    )
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "UserPromptSubmit",
            "additionalContext": context,
        }
    }))


def send_hook(hook_input):
    command = hook_input.get("tool_input", {}).get("command", "")
    if not SEND_PATTERN.search(command):
        return
    # What leaves the machine is the committed HEAD, not the working tree.
    if head_tree_hash() in reviewed_trees():
        return
    reason = (
        "The commit being sent has not had an automated review. Strongly "
        f"recommended: deny this, and have Claude run the `{REVIEW_SKILL}` "
        "skill first. Approve to send it unreviewed."
    )
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "ask",
            "permissionDecisionReason": reason,
            "additionalContext": (
                f"If the user denies this, run the `{REVIEW_SKILL}` skill, "
                "then try again."
            ),
        }
    }))


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    mode = sys.argv[1]
    if mode == "mark":
        mark()
    elif mode == "status":
        sys.exit(status())
    elif mode in ("prompt-hook", "send-hook"):
        try:
            hook_input = json.load(sys.stdin)
            if mode == "prompt-hook":
                prompt_hook(hook_input)
            else:
                send_hook(hook_input)
        except (subprocess.CalledProcessError, OSError, ValueError):
            # Outside a git checkout, or with no commits yet, there is nothing to
            # check, and a hook error must not block the user's work.
            pass
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
