#!/usr/bin/env bash
# publish-docs.sh: Entrypoint for Travis to build and publish Doxygen documentation for vg
# Based on the method of https://gist.github.com/domenic/ec8b0fc8ab45f39403dd

set -e

# Configuration
# What branch do the docs go on?
DEST_BRANCH="gh-pages"
# What repo do they go in? Must be an SSH repo specifier for committing to.
DEST_REPO="git@github.com:vgteam/vg.git"
# What directory, relative to the source repo's root, do the built docs come from?
# Probably needs to have a trailing slash
SOURCE_DIR="doc/doxygen/html/"
# What directory, relative to the dest repo's root, do the built docs go to?
# Also probably needs a trailing slash
DEST_DIR="./"
# Who should be seen as making the commits?
COMMIT_AUTHOR_NAME="VG Doc Bot"
COMMIT_AUTHOR_EMAIL="anovak+vgdocbot@soe.ucsc.edu"

# We expect GITLAB_SECRET_FILE_DOCS_SSH_KEY to come in from the environment,
# specifying the private deploy key we will use to get at the docs repo.

# Make sure no submodules have untracked files from caching
# See <https://gist.github.com/nicktoumpelis/11214362#file-repo-rinse-sh-L2>
git submodule foreach --recursive git clean -xfd

# Find all the submodules that Doxygen wants to look at and make sure we have
# them.
#
# Gitlab does some exciting reconfiguration of our Git repository (see
# <https://docs.gitlab.com/ci/runners/git_submodules/#check-out-nested-submodules>
# and the part about "externalizing" the Git configuration). It's probably
# something with `url.<local-path>.insteadOf`, but I can't find a source to
# back that up.
#
# This (along with Git refusing to let you work with local path source
# repositories because it can't figure out how to do it securely) is probably
# the reason that attempts to update to some commits in submodules fail with:
# 
#     fatal: transport 'file' not allowed
#
# (See <https://gitlab.com/gitlab-org/gitlab-runner/-/work_items/38908>.) It's
# not clear exactly why Git thinks we need to fetch the commits, or why this
# happens for some commits and not others.
#
# The security concern here is CVE-2022-39253 in Git, where, when you clone
# from a local path, Git makes a copy of all the `.git/objects` files, and so
# if those files are symlinks to your secret files, then you can end up making
# copies of those files with more permissive permissions in your clone, which
# other users on the same system could then read. That's not a big concern in
# our CI environment, because there aren't any lower-permission users inside
# the CI container and so anyone who could read those files could read the
# secrets anyway, but it's still nice to not flag off security things.
#
# We could just set the Gitlab job to recursively clone, and not touch the
# submodules here, which works, but that's slow.
#
# We could use GIT_SUBMODULE_PATHS to tell the job which submodules to populate
# (see
# <https://docs.gitlab.com/ci/runners/configure_runners/#sync-or-exclude-specific-submodules-from-ci-jobs>),
# but then that would need to always match the Doxygen config.
#
# So instead we do some rocket surgery on the Git repo. THIS WILL DESTROY LOCAL
# SUBMODULE WORKING COPIES AND THEIR GIT HISTORY, so DO NOT run this script if
# you have any commits that aren't pushed elsewhere!
#
# This ritual was devised by Anthropic Claude, and it works, but we're not
# actually fully able to explain why.
GITDIR=$(git rev-parse --git-dir)
DOXYGEN_DEPS=$(cat Doxyfile | grep "^INPUT *=" | cut -f2 -d'=' | tr ' ' '\n' | grep "^ *deps" | sed 's_ *\(deps/[^/]*\).*_\1_' | sort | uniq)
for dep in ${DOXYGEN_DEPS}; do
    # Tell Git to stop maintaining the work tree for the submodule, if it is.
    git submodule deinit -f -- "${dep}" || true
    # Remove the work tree and the Git repository information for the
    # submodule. This is where we think the file:// origins might be hiding.
    rm -rf "${dep}" "${GITDIR}/modules/${dep}"
done
# Now we re-clone those submodules, which should use the URLs they usually use
# instead of whatever Gitlab did.
#
# Explicitly refusing to use file protocols here (instead of the default "user"
# mode of using them when asked by the user, as documented at
# <https://git-scm.com/docs/git-config#Documentation/git-config.txt-protocolallow>)
# *might* protect us from the case that a submodule in vg directly references a
# specially prepared Git repo that has been smuggled into CI, alongside some
# kind of cache misconfiguration that would let the cloned repo smuggle a
# secret into the cache, where it could be read by other malicious code
# smuggled in that can only read the cache and not secrets, all orchestrated by
# someone who can make exactly this limited set of malicious commits. This is
# obviously absurd, and it's not clear that the "user" mode allows "file" for
# root submodules anyway, but I'm not about to let Anthropic Claude do more
# security theater than me, so it stays.
echo "${DOXYGEN_DEPS}" | xargs -n 1 git -c protocol.file.allow=never submodule update --init --recursive

# Build the documentation.
# Assumes we are running in the repo root.
make docs

# Get ready to deploy the docs

# Make a scratch directory *outside* our normal git repo
SCRATCH_DIR="$(mktemp -d)"
# And clean it up when we stop
function cleanup {
  rm -Rf ${SCRATCH_DIR}
}
trap cleanup EXIT



# Set up our SSH key
touch "${SCRATCH_DIR}/deploy_key"

# Protect it so the agent is happy
chmod 600 "${SCRATCH_DIR}/deploy_key"

# Fill it in with NO COMMAND ECHO
set +x
echo "${GITLAB_SECRET_FILE_DOCS_SSH_KEY}" > ${SCRATCH_DIR}/deploy_key

# Turn on echo so we can see what we're doing.
# This MUST happen only AFTER we are done touching the encryption stuff.
set -x

# Make sure we have an known_hosts
mkdir -p ~/.ssh
touch ~/.ssh/known_hosts
cat ~/.ssh/known_hosts

# Clone the dest repo, now that we can authenticate.
# Don't check it out, so we can get just the branch we want or start a new branch with a clean working copy.
git -c "core.sshCommand=ssh -i ${SCRATCH_DIR}/deploy_key -o 'UserKnownHostsFile=/dev/null' -o 'StrictHostKeyChecking=no'" clone --no-checkout "${DEST_REPO}" "${SCRATCH_DIR}/dest"

# Go in and get/make the destination branch
pushd "${SCRATCH_DIR}/dest"
git checkout "${DEST_BRANCH}" || git checkout --orphan "${DEST_BRANCH}"
popd

# Drop the files in
# See https://explainshell.com/explain?cmd=rsync+-aqr+--delete+--exclude
# We need to not clobber any .git in the destination.
rsync -avr "${SOURCE_DIR}" "${SCRATCH_DIR}/dest/${DEST_DIR}" --delete --exclude .git

# Go back in to make the commit
pushd "${SCRATCH_DIR}/dest"

# Disable Jeckyll processing for Github Pages since we did it already
touch .nojekyll
git add .nojekyll

# Add all the files here (except hidden ones) and add deletions
git add -A

# Become the user we want to be
git config user.name "${COMMIT_AUTHOR_NAME}"
git config user.email "${COMMIT_AUTHOR_EMAIL}"

# Make the commit. Tolerate failure because this fails when there is nothing to commit.
git commit -m "Commit new auto-generated docs" || true

if [[ -z "${CI_COMMIT_BRANCH}" || "${CI_COMMIT_BRANCH}" != "${CI_DEFAULT_BRANCH}" ]]; then
    # If we're not a real mainline commit, we just make sure the docs build.
    echo "Documentation should not be deployed because this is not a mainline build"
    exit 0
fi

# If we are on the right branch, actually push the commit.
# Push the commit. This does not fail if there is no commit.
git -c "core.sshCommand=ssh -i ${SCRATCH_DIR}/deploy_key -o 'UserKnownHostsFile=/dev/null' -o 'StrictHostKeyChecking=no'" push origin "${DEST_BRANCH}"




