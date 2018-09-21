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
SOURCE_DIR="doc/doxygen/"
# What directory, relative to the dest repo's root, do the built docs go to?
# Also probably needs a trailing slash
DEST_DIR="./"
# Who should be seen as making the commits?
COMMIT_AUTHOR_NAME="Travis Doc Bot"
COMMIT_AUTHOR_EMAIL="anovak+travisdocbot@soe.ucsc.edu"
# What SSH key, relative to this repo's root, should we decrypt and use for doc deployment?
ENCRYPTED_SSH_KEY_FILE="doc/deploy_key.enc"

# We expect ENCRYPTION_LABEL to come in from the environment, specifying the ID
# of the encrypted deploy key we will use to get at the docs repo.

# Build the documentation.
# Assumes we are running in the repo root.
make docs

# Get ready to deploy the docs

# Make a scratch directory
mkdir -p ./tmp

# Get our encryption key and IV variable names
ENCRYPTION_KEY_VAR="encrypted_${ENCRYPTION_LABEL}_key"
ENCRYPTION_IV_VAR="encrypted_${ENCRYPTION_LABEL}_iv"

# Decrypt the encrypted deploy SSH key
# Get the key and IV from the variables we have the names of.
openssl aes-256-cbc -K "${!ENCRYPTION_KEY_VAR}" -iv "${!ENCRYPTION_IV_VAR}" -in "${ENCRYPTED_SSH_KEY_FILE}" -out ./tmp/deploy_key -d
# Protect it so the agent is happy
chmod 600 ./tmp/deploy_key

# Start an agent and add the key
eval "$(ssh-agent -s)"
ssh-add ./tmp/deploy_key

# Check out the dest repo, now that we can authenticate, shallow-ly to avoid getting all history
git clone "${DEST_REPO}" ./tmp/dest

# Go in and get/make the destination branch
cd ./tmp/dest
git checkout "${DEST_BRANCH}" || git checkout --orphan "${DEST_BRANCH}"
cd ../..

# Drop the files in
# See https://explainshell.com/explain?cmd=rsync+-aqr+--delete+--filter
# We need to not clobber any .git in the destination.
rsync -aqr "${SOURCE_DIR}" "tmp/dest/${DEST_DIR}" --delete --filter='protect .git/**/*'

cd ./tmp/dest

if git diff --quiet; then
    # Don't commit nothing
    echo "No changes to commit"
    exit 0
fi

# Become the user we want to be
git config user.name "${COMMIT_AUTHOR_NAME}"
git config user.email "${COMMIT_AUTHOR_EMAIL}"

# Make the commit
git commit -am "Commit new auto-generated docs"

if [[ "${TRAVIS_PULL_REQUEST}" != "false" || "${TRAVIS_BRANCH}" != "master" ]]; then
    # If we're not a real master commit, we just did all this for testing purposes.
    echo "Documentation should not be deployed"
    exit 0
fi

# If we are on the right branch, actually push the commit.
# Push the commit
git push origin "${DEST_BRANCH}"




