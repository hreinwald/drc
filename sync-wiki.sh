#!/bin/bash
# Script to sync wiki content from the main repository to the GitHub Wiki repository
# This script clones the wiki repo, copies the content, and pushes the changes

set -e

echo "Syncing wiki content to GitHub Wiki..."

# Create temporary directory
TEMP_DIR=$(mktemp -d)
cd "$TEMP_DIR"

# Clone the wiki repository
echo "Cloning wiki repository..."
git clone https://github.com/hreinwald/drc.wiki.git
cd drc.wiki

# Copy wiki files from main repo
echo "Copying wiki files..."
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cp "$SCRIPT_DIR/wiki/Home.md" .
cp "$SCRIPT_DIR/wiki/_Sidebar.md" .

# Check if there are changes
if git diff --quiet && git diff --cached --quiet; then
    echo "No changes to sync"
    cd /
    rm -rf "$TEMP_DIR"
    exit 0
fi

# Commit and push changes
echo "Committing changes..."
git add Home.md _Sidebar.md
git commit -m "Update wiki with documentation links and sidebar

- Add comprehensive Home.md with drc package introduction
- Add links to documentation site at hreinwald.github.io/drc
- Include quick links table for documentation, workflows, repo, issues, and discussions
- Add _Sidebar.md for persistent navigation
- Include installation instructions"

echo "Pushing to GitHub Wiki..."
git push origin master

# Cleanup
cd /
rm -rf "$TEMP_DIR"

echo "✓ Wiki sync completed successfully!"
echo "Visit https://github.com/hreinwald/drc/wiki to see the updated wiki"
