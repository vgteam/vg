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
COMMIT_AUTHOR_NAME="Travis Doc Bot"
COMMIT_AUTHOR_EMAIL="anovak+travisdocbot@soe.ucsc.edu"
# What SSH key, relative to this repo's root, should we decrypt and use for doc deployment?
ENCRYPTED_SSH_KEY_FILE="doc/deploy_key.enc"

# We expect DOCS_KEY_ENCRYPTION_LABEL to come in from the environment, specifying the ID
# of the encrypted deploy key we will use to get at the docs repo.

# Build the documentation.
# Assumes we are running in the repo root.
make docs

if [[ ! -z "${TRAVIS_PULL_REQUEST_SLUG}" && "${TRAVIS_PULL_REQUEST_SLUG}" != "${TRAVIS_REPO_SLUG}" ]]; then
    # This is an external PR. We have no access to the encryption keys for the encrypted deploy SSH key.
    # We want to check out the dest repo with that key because it's much simpler than hacking the remote from https to ssh.
    # So we won't even test copying the docs over to the destination repo.
    echo "Not testing deploy; no encryption keys available for external PRs."
    exit 0
fi

# Get ready to deploy the docs

# Make a scratch directory *outside* our normal git repo
SCRATCH_DIR="../tmp"
mkdir -p "${SCRATCH_DIR}"

# Get our encryption key and IV variable names
ENCRYPTION_KEY_VAR="encrypted_${DOCS_KEY_ENCRYPTION_LABEL}_key"
ENCRYPTION_IV_VAR="encrypted_${DOCS_KEY_ENCRYPTION_LABEL}_iv"

echo "Want to decrypt ${ENCRYPTED_SSH_KEY_FILE} using key from variable ${ENCRYPTION_KEY_VAR} and IV from variable ${ENCRYPTION_IV_VAR}"

if [[ -z "${!ENCRYPTION_KEY_VAR}" ]]; then
    echo "Encryption key not found!"
    exit 1
fi

if [[ -z "${!ENCRYPTION_IV_VAR}" ]]; then
    echo "Encryption IV not found!"
    exit 1
fi

# Decrypt the encrypted deploy SSH key
# Get the key and IV from the variables we have the names of.
openssl aes-256-cbc -K "${!ENCRYPTION_KEY_VAR}" -iv "${!ENCRYPTION_IV_VAR}" -in "${ENCRYPTED_SSH_KEY_FILE}" -out "${SCRATCH_DIR}/deploy_key" -d
# Protect it so the agent is happy
chmod 600 "${SCRATCH_DIR}/deploy_key"

# Start an agent and add the key
eval "$(ssh-agent -s)"
ssh-add "${SCRATCH_DIR}/deploy_key"

# Turn on echo so we can see what we're doing.
# This MUST happen only AFTER we are done toucking the encryption stuff.
set -x

# Clone the dest repo, now that we can authenticate.
# Don't check it out, so we can get just the branch we want or start a new branch with a clean working copy.
git clone --no-checkout "${DEST_REPO}" "${SCRATCH_DIR}/dest"

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

# Add all the files here (except hidden ones) and add deletions
git add -A

# Become the user we want to be
git config user.name "${COMMIT_AUTHOR_NAME}"
git config user.email "${COMMIT_AUTHOR_EMAIL}"

# Make the commit. Tolerate failure because this fails when there is nothing to commit.
git commit -m "Commit new auto-generated docs" || true

if [[ "${TRAVIS_PULL_REQUEST}" != "false" || "${TRAVIS_BRANCH}" != "master" ]]; then
    # If we're not a real master commit, we just make sure the docs build.
    # Also, unless we're a branch in the main vgteam/vg repo, we don't have access to the encryption keys anyway.
    # So we can't even try to deploy.
    echo "Documentation should not be deployed because this is not a mainline master build"
    exit 0
fi

# If we are on the right branch, actually push the commit.
# Push the commit. This does not fail if there is no commit.
git push origin "${DEST_BRANCH}"




