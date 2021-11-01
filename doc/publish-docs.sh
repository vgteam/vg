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

# Find all the submodules that Doxygen wants to look at and make sure we have
# those.
cat Doxyfile  | grep "^INPUT *=" | cut -f2 -d'=' | tr ' ' '\n' | grep "^ *deps" | sed 's_ *\(deps/[^/]*\).*_\1_' | sort | uniq | xargs -n 1 git submodule update --init --recursive

# Build the documentation.
# Assumes we are running in the repo root.
make docs

# Get ready to deploy the docs

# Make a scratch directory *outside* our normal git repo
SCRATCH_DIR="$(pwd)/../tmp"
mkdir -p "${SCRATCH_DIR}"

# Despite the _SECRET_FILE ending we used to use the contents of
# GITLAB_SECRET_FILE_DOCS_SSH_KEY as the key. But now we want to support it as
# a file instead.

# Turn off echo to protect key
set +x

if [[ -e "${GITLAB_SECRET_FILE_DOCS_SSH_KEY}" ]] ; then
    # It's a file, so copy it into place
    cp "${GITLAB_SECRET_FILE_DOCS_SSH_KEY}" "${SCRATCH_DIR}/deploy_key"
else
    # Old style CI: it's a value.

    # Set up our SSH key
    touch "${SCRATCH_DIR}/deploy_key"

    # Protect it so the agent is happy
    chmod 600 "${SCRATCH_DIR}/deploy_key"

    # Fill it in
    echo "${GITLAB_SECRET_FILE_DOCS_SSH_KEY}" > ${SCRATCH_DIR}/deploy_key
fi

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




